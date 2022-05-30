"""@package codegen.codegen_util.py
Util functions used during code generation.

@author: Johannes Seiferth
"""

from multiprocessing import cpu_count, Pool
from re import sub
from shlex import split as shlex_split
from shutil import which
from subprocess import run, PIPE, CalledProcessError
from sys import maxsize as sys_maxsize
from typing import Dict, List, Optional, Tuple, Union

from pathlib2 import Path

import offsite.config
from offsite.config import Config
from offsite.util.math_utils import eval_math_expr

# Indention used for code formatting.
INDENTION = '  '


def indent(lvl: int) -> str:
    """Write indention string.

    Parameters:
    -----------
    indention: int
        Used indention level.

    Returns:
    --------
    str
        Written indention string.
    """
    return INDENTION * lvl


def eval_loop_boundary(boundary: str, constants: Optional[List[Tuple[str, Union[str, float, int]]]]) -> str:
    """Evaluate loop boundary expression.

    Parameters:
    -----------
    method: ODEMethod
        Used ODE method.
    boundary: str
        Boundary expression.

    Returns:
    --------
    str
        Evaluated boundary expression.
    """
    boundary = eval_math_expr(str(boundary), constants)
    try:
        boundary = int(boundary)
    except TypeError:
        pass
    return str(boundary)


def replace_var_with_factor(var: str, computation: str, factor: str) -> str:
    """Replace loop variable term with an integer term.

    Parameters:
    -----------
    var: str
        variable as string.
    computation: str
        computation string
    factor: str
        factor to replace var

    Returns:
    --------
    str
        Computation string with replaced var.
    """
    return computation.replace('[' + var + ']', '[' + str(factor) + ']')


def replace_incr_with_assign_op(computation: str) -> str:
    """Replace increment expression statements with an assignment statement.

    Parameters:
    -----------
    computation: str
        Computation string.

    Returns:
    --------
    str
        Computation string with replaced operator.
    """
    return sub(r'(\+=|-)', '=', computation)


def substitute_rhs_call(computation: str, component: str, input_vector: str,
                        constants: Optional[List[Tuple[str, Union[str, float, int]]]]) -> str:
    """Substitute all RHS calls of the given computation code with the passed IVP component.

    Parameters:
    -----------
    computation: str
        Code that contains one or more RHS calls.
    component: str
        IVP component used to replace RHS calls.
    input_vector: str
        Input vector used in RHS call.
    constants: list of Tuple (str, str)
        Constants of the IVP.

    Returns:
    --------
    str
        Code with substituted RHS calls.
    """
    assert '%RHS' in computation
    config: Config = offsite.config.offsiteConfig
    component = component.strip()
    # Substitute keywords %in and %last_idx.
    component = component.replace('%in', input_vector)
    component = component.replace('%last_idx', 'last')
    # Replace index keyword %idx and %last_idx with the variable names defined in the given config.
    component = component.replace('%idx', config.var_idx)
    # Replace constants in component string.
    for constant_name, constant_value in constants:
        regex = r'(?![a-zA-Z0-9]-_)' + constant_name + r'(?![a-zA-Z0-9-_])'
        if constant_name in ['g', 'n']:
            component = sub(regex, str(eval_math_expr(constant_value, constants, cast_to=int)), component)
        else:
            component = sub(regex, str(eval_math_expr(constant_value, constants)), component)
    # Substitute RHS call.
    return computation.replace('%RHS', component)


def substitute_rhs_func_call(computation: str, rhs_func: str, input_vector: str, butcher_node: str) -> str:
    """Substitute all RHS calls of the given computation code with the passed RHS function call.

    Parameters:
    -----------
    computation: str
        Code that contains one or more RHS calls.
    rhs_func: str
        Name of the called function type.
    input_vector: str
        Input vector passed to the RHS function call.
    butcher_node: str
        Used butcher node coefficient.

    Returns:
    --------
    str
        Code with substituted RHS calls.
    """
    assert '%RHS' in computation
    config: Config = offsite.config.offsiteConfig
    # Build RHS function call.
    if rhs_func == 'eval_range':
        # Remove loop over system dimension. Eval_range runs over the system range internally.
        lhs = shlex_split(computation)[0].replace('[{}]'.format(config.var_idx), '')
        # Construct RHS call.
        computation = 'eval_range(first, last, t + {} * h, {}, {})'.format(butcher_node, input_vector, lhs)
    elif rhs_func == 'eval_component':
        rhs_func = 'eval_component({}, t + {} * h, {})'.format(config.var_idx, butcher_node, input_vector)
        computation = computation.replace('%RHS', rhs_func)
    else:
        assert False
    return computation


def substitute_stencil_call(computation: str, constants: Optional[List[Tuple[str, Union[str, float, int]]]],
                            input_vector: str, butcher_node: str, stencil_path: str) -> str:
    """Substitute all RHS calls of the given computation code with a Yasksite-style stencil call.

    Parameters:
    -----------
    computation: str
        Code that contains one or more RHS calls.
    constants: list of Tuple (str, str)
        Constants of the IVP.
    input_vector: str
        Input vector passed to the RHS function call.
    butcher_node: str
        Used butcher node coefficient.
    stencil_path: str
        Path to stencil code.

    Returns:
    --------
    str
        Code with substituted RHS calls.
    """
    assert '%RHS' in computation
    # Write std::map string that contains all constants.
    constant_map = '{'
    for constant_name, constant_value in constants:
        constant_map += '{\"' + constant_name + '\", ' + constant_value + '}'
    constant_map += '}'
    # Build Yasksite stencil.
    stencil = 'RHS({}, {}, {}, {})'.format(input_vector, constant_map, Path(stencil_path).stem, butcher_node)
    # Substitute RHS call.
    return computation.replace('%RHS', stencil)


def write_closing_bracket(indent_lvl: int) -> str:
    """Write closing bracket.

    Parameters:
    -----------
    indent_lvl: int
        Current indention level of the surrounding code.

    Returns:
    --------
    str
        Written closing bracket code line.
    """
    return indent(indent_lvl) + '}' + '\n'


def write_tiling_loop(block_varname: str, indent_lvl: int) -> str:
    """Write tiling loop.

    Parameters:
    -----------
    block_varname: str
        Name suffix of the block size variable.
    indent_lvl: int
        Current indention level of the surrounding code.

    Returns:
    --------
    str
        Written tiling loop code line.
    """
    bs_var = 'bs_{}'.format(block_varname)
    limit_var = 'limit_{}'.format(block_varname)
    block_var = 'B'.format(block_varname)
    #
    loop = indent(indent_lvl) + 'int {}, {};\n'.format(bs_var, limit_var)
    loop += indent(indent_lvl) + 'for (int jj=first, {0}=imin({1}, last-first+1), '.format(bs_var, block_var)
    loop += '{0}=imax(first, last+1-{1}); jj<=last; {2}=last+1-jj, {0}=last)\n'.format(limit_var, block_var, bs_var)
    loop += indent(indent_lvl) + '{\n'
    indent_lvl += 1
    loop += indent(indent_lvl) + 'for (; jj <= {}; jj += {})'.format(limit_var, bs_var)
    loop += indent(indent_lvl) + '{\n'
    return loop


def write_instrument_kernel_start(indent_lvl: int) -> str:
    """Write instrumentation code that would be run before kernel execution.

    Parameters:
    -----------
    indent_lvl: int
        Current indention level of the surrounding code.

    Returns:
    --------
    str
        Written instrumentation code.
    """
    # Write code to instrument the kernel code.
    instr_start = '#ifdef INSTRUMENT\n'
    instr_start += indent(indent_lvl) + '{\n'
    indent_lvl += 1
    instr_start += indent(indent_lvl) + '#pragma omp barrier\n'
    instr_start += indent(indent_lvl) + 'time_snap_t time;\n'
    instr_start += indent(indent_lvl) + 'if (me == 0)\n'
    instr_start += indent(indent_lvl) + '{\n'
    instr_start += indent(indent_lvl + 1) + 'time_snap_start(&time);\n'
    indent_lvl -= 1
    instr_start += indent(indent_lvl) + '}\n'
    instr_start += '#endif\n'
    return instr_start


def write_instrument_kernel_end(indent_lvl: int, kernel_id: int) -> str:
    """Write instrumentation code that would be run after kernel execution.

    Parameters:
    -----------
    indent_lvl: int
        Current indention level of the surrounding code.
    kernel_id: int
        Database id of the kernel.

    Returns:
    --------
    str
        Written instrumentation code.
    """
    # Write code to instrument the kernel code.
    instr_end = '#ifdef INSTRUMENT\n'
    instr_end += indent(indent_lvl) + '#pragma omp barrier\n'
    instr_end += indent(indent_lvl) + 'if (me == 0)\n'
    instr_end += indent(indent_lvl) + '{\n'
    indent_lvl += 1
    instr_end += indent(indent_lvl) + 'double T = time_snap_stop(&time);\n'
    # Has to be raw to keep the newline character in printf.
    instr_end += '#ifdef _OPENMP\n'
    instr_end += indent(indent_lvl) + \
                 r'printf("#Kernel={}\§t#Threads=%d\§t%.20e\§n", omp_get_num_threads(), T / 1e9 / n);'.format(kernel_id)
    instr_end += '\n#else\n'
    instr_end += indent(indent_lvl) + \
                 r'printf("#Kernel={}\§t#Threads=1\§t%.20e\§n", T / 1e9 / n);'.format(kernel_id)
    instr_end += '\n#endif\n'
    indent_lvl -= 1
    instr_end += indent(indent_lvl) + '}\n'
    indent_lvl -= 1
    instr_end += indent(indent_lvl) + '}\n'
    instr_end += '#endif\n'
    return instr_end


def create_variant_name(variant: List['Kernel'], skeleton: str) -> str:
    """Create filename of an implementation variant.

    Parameters:
    -----------
    variant: list of Kernel
        Kernels associated with this implementation variant.
    skeleton: str
        Name of the associated ImplSkeleton.

    Returns:
    --------
    str
        Created filename.
    """
    return '{}_{}'.format(skeleton, '_'.join((kernel.name for kernel in variant)))


def format_codefile(path: Path):
    """Format code file with command line tool indent.

    Parameters:
    -----------
    path: Path
        Path to code file.

    Returns:
    --------
    -
    """
    if which('indent'):
        try:
            cmd = ['indent', '{}'.format(path)]
            run(cmd, check=True, stdout=PIPE).stdout
        except CalledProcessError as error:
            raise RuntimeWarning('Formatting file {} with indent failed: {}'.format(path, error))


def write_codes_to_file(codes: Dict[str, str], folder: Path = None, suffix: str = '.c', num_workers: int = sys_maxsize,
                        format_code: bool = True) -> List[Path]:
    """Write codes to file.

    Parameters:
    -----------
    codes: dict (key=str, val=str)
        Code strings (val) and their associated file names (key).
    folder: Path
        Relative path to folder that will contain the written code files.
    suffix: str
        File suffix used for all generated codes files.
    num_workers: int
        Number of worker threads used to write code files.
    format_code: bool
        Use indent to format the code.

    Returns:
    --------
    list of Path
        Written code files.
    """
    written_files: List[Path] = list()
    if num_workers > 1:
        # Create worker threads ...
        num_workers = min(cpu_count() - 1, len(codes), num_workers)
        num_workers = max(num_workers, 1)
        # ... and init the worker pool.
        pool = Pool(num_workers)
        errors = list()
        for file_name, code in codes.items():
            pool.apply_async(write_code_to_file, args=(code, file_name, folder, suffix, format_code),
                             callback=written_files.append, error_callback=errors.append)
        # Wait for threads.
        pool.close()
        pool.join()
        # Raise error if failed.
        if errors:
            raise RuntimeError('Failed to write code to files: Error in worker threads.')
    else:
        for file_name, code in codes.items():
            written_files.append(write_code_to_file(code, file_name, folder, suffix, format_code))
    return written_files


def write_code_to_file(
        code_str: str, file_name: str, folder: Path = None, suffix: str = '.c', format_code: bool = True) -> Path:
    """Write code to file.

    Parameters:
    -----------
    code_str: str
        Source code string.
    file_name: str
        Requested file name.
    folder: Path
        Relative path to folder that will contain the written code file.
    file_suffix: str
        File suffix used for the generated codes file.
    format_code: bool
        Use indent to format the code.

    Returns:
    --------
    Path
        Path to written code file.
    """
    # Create folder if it does not yet exist.
    if folder is not None:
        folder = Path(folder)
        if folder and not folder.exists():
            folder.mkdir(parents=True)
    # Write code to file.
    file_path = Path('{}{}'.format(file_name, suffix))
    if folder is not None:
        file_path = folder / file_path
    with file_path.open('w') as file_handle:
        file_handle.write(code_str)
    # Replace newline characters in printf commands.
    try:
        cmd = ['sed', '-i', 's/§n/n/g', '{}'.format(file_path)]
        run(cmd, check=True, stdout=PIPE).stdout
    except CalledProcessError as error:
        raise RuntimeError('Unable to substitute newline character stubs in {}: {}'.format(file_path, error))
    # Replace tabulator characters in printf commands.
    try:
        cmd = ['sed', '-i', 's/§t/t/g', '{}'.format(file_path)]
        run(cmd, check=True, stdout=PIPE).stdout
    except CalledProcessError as error:
        raise RuntimeError('Unable to substitute tabulator character stubs in {}: {}'.format(file_path, error))
    # Format code file.
    if format_code:
        format_codefile(file_path)
    return file_path
