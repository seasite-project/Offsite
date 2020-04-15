"""@package impl_generator
Definition of class ImplCodeGenerator.
"""

from pathlib import Path
from subprocess import run, PIPE, CalledProcessError
from typing import Dict, List

import attr
from regex import sub

from offsite.codegen.code_generator import CodeGenerator
from offsite.codegen.code_object import CodeStr
from offsite.codegen.codegen_util import format_codefile, gen_loop_string, indent
from offsite.codegen.kernel_generator import KernelCodeGenerator
from offsite.descriptions.impl_skeleton import ImplSkeleton
from offsite.descriptions.ivp import IVP
from offsite.descriptions.kernel_template import KernelTemplate, Kernel
from offsite.descriptions.ode_method import ODEMethod
from offsite.evaluation.math_utils import eval_math_expr


@attr.s(kw_only=True)
class ImplCodeGenerator(CodeGenerator):
    """Representation of the code generator for implementation variant code.

    Attributes:
    -----------
    db_session : sqlalchemy.orm.session.Session
        Used database session.
    folder : Path
        Relative path to the location of the code files generated.
    filepaths : xxx
        xxx
    skeleton : ImplSkeleton
        ImplSkeleton to which the impl variants belong for which code is generated. Set by generate() function.
    loaded_templates : Dict of KernelTemplate
        Dict of KernelTemplate for which code is generated.
    """
    db_session = attr.ib(type='sqlalchemy.orm.session.Session')
    folder = attr.ib(type=Path, default=Path('tmp/variants/'))
    filepaths = attr.ib(type=Dict, default=dict())
    skeleton = attr.ib(type='ImplSkeleton', default=None, init=False)
    loaded_templates = attr.ib(type=Dict, default=dict())
    kernel_stubs = attr.ib(type=Dict, default=dict())
    stub_idx = attr.ib(type=int, default=0, init=False)

    def __attrs_post_init__(self):
        """Create this ImplCodeGenerator object.

        Parameters:
        -----------
        -

        Returns:
        --------
        -
        """
        self.parser = {
            '%LOOP_START': self.resolve_loop_start,
            '%LOOP_END': self.resolve_loop_end,
            '%COM': self.resolve_comm,
            '%KERNEL': self.resolve_kernel_template,
            '%CMD': self.resolve_cmd
        }

    def generate(self, skeleton: ImplSkeleton, impl_variants: List[List[str]], ivp: IVP,
                 method: ODEMethod) -> Dict[Path, str]:
        # Set some members to match current code generation.
        self.skeleton = skeleton
        self.ivp = ivp
        self.method = method
        # Set default indention level.
        ind_lvl = 1
        # Reset kernel stub dictionary.
        self.kernel_stubs = dict()
        # Parse impl skeleton object and create impl code.
        impl = self.parse_skeleton(ind_lvl)
        # Apply required loop splits.
        if 'loop splits' in skeleton.codegen:
            for loop_split in skeleton.codegen['loop splits']:
                loop_split = loop_split.split(' ')
                impl.split_loop(int(loop_split[1]) + 1, loop_split[0])
        # Write generated implementation variants to files.
        generated_files = {}
        for vid, variant in impl_variants:
            kernels = self.derive_impl_variant_kernels(variant)
            impl_variant_str = self.generate_impl_variant(impl, kernels, vid)
            path = Path('{}/{}.h'.format(self.folder, self.create_filename(kernels)))
            generated_files[path] = impl_variant_str
        return generated_files

    def generate_impl_variant(self, impl: CodeStr, kernels: List[Kernel], variant_id: int) -> List[str]:
        # Insert kernel codes into implementation skeleton.
        impl_str = impl.get()
        for stub_name, stub in self.kernel_stubs.items():
            # Retrieve Kernel object.
            template_name = stub[0]
            kernel = next(filter(lambda x: x.template.name == template_name, kernels))
            # Generate kernel code.
            ode_solution_vector = stub[1]
            code = KernelCodeGenerator(db_session=self.db_session, config=self.config).generate(
                kernel, self.method, self.ivp, ode_solution_vector)
            # Substitute kernel code stub with actual kernel code.
            impl_str = self.substitute_kernel_stub(impl_str, stub_name, code.code_string, template_name, kernel.db_id)
        # Write implementation variant.
        return self.write_skeleton_includes(variant_id) + self.write_function_description() + \
               self.write_instrument_impl(1) + impl_str + '}'

    def parse_skeleton(self, lvl_ind: int) -> CodeStr:
        """Parse ImplSkeleton object and create its corresponding impl code object.

        Parameters:
        -----------
        lvl_ind : int
            Indention level.

        Returns:
        --------
        CodeStr
            Generated impl code.
        """
        # Initialize CodeStr object.
        code = CodeStr('', lvl_ind)
        # Parse impl skeleton.
        for line in self.skeleton.code.split('\n'):
            func_tuple = self.parse_line(line)
            if func_tuple[0]:
                func_tuple[0](code, *func_tuple[1])
        return code

    def resolve_loop_start(self, code_str: CodeStr, *args):
        """Resolve loop_start code line.

        Parameters:
        -----------
        code_str : CodeStr
            Container to store the created code.
        args
            Required function parameters.

        Returns:
        --------
        -
        """
        loop_var = args[0]
        loop_iterations = args[1]
        # Write loop header code.
        self.write_loop(loop_var, loop_iterations, code_str)

    def resolve_loop_end(self, code_str: CodeStr, *args):
        """Resolve loop_end code line.

        Parameters:
        -----------
        code_str : CodeStr
            Container to store the created code.
        args
            Required function parameters.

        Returns:
        --------
        -
        """
        self.write_closing_brackets(1, code_str)

    def resolve_comm(self, code_str: CodeStr, *args):
        """Resolve communication code line.

        Parameters:
        -----------
        code_str : CodeStr
            Container to store the created code.
        args
            Required function parameters.

        Returns:
        --------
        -
        """
        command_split = args[0].split('_')
        lib = command_split[0]
        lib_command = command_split[1]
        if lib == 'omp':
            self.write_omp_operation(lib_command, code_str)
        elif lib == 'pthreads':
            self.write_pthreads_operation(lib_command, code_str)
        else:
            assert False

    def resolve_kernel_template(self, code_str: CodeStr, *args):
        """Resolve kernel template code line.

        Parameters:
        -----------
        code_str : CodeStr
            Container to store the created code.
        args
            Required function parameters.

        Returns:
        --------
        -
        """
        template_name = args[0]
        # Load kernel template from database.
        if template_name not in self.loaded_templates.keys():
            template = KernelTemplate.from_database(self.db_session, template_name)
            self.loaded_templates[template_name] = template
        # Write kernel template code.
        self.write_kernel_template_stub(template_name, args[1:], code_str)

    def resolve_cmd(self, code_str: CodeStr, *args):
        """Resolve command code line.

        Parameters:
        -----------
        code_str : CodeStr
            Container to store the created code.
        args
            Required function parameters.

        Returns:
        --------
        -
        """
        cmd_name = args[0]
        if cmd_name == 'SWAP':
            assert len(args) == 4
            self.write_double_pointer_swap(args[1], args[2], args[3], code_str)
        else:
            assert False

    def substitute_kernel_stub(self, impl_str: str, stub: str, kernel_code: CodeStr, template: str, kernel_id: int):
        # Write code to instrument the kernel code.
        indent_lvl = kernel_code.get_indent()
        instrument_start = self.write_instrument_kernel_start(indent_lvl)
        instrument_end = self.write_instrument_kernel_end(indent_lvl, kernel_id)
        # Substitute kernel stub.
        impl_str = sub(r'(?![\w\d\s]){}(?![\w\d\S])'.format(stub), r'//{} %{}\n{}{}{}'.format(
            template, kernel_id, instrument_start, kernel_code.string, instrument_end), impl_str)
        return impl_str

    @staticmethod
    def write_instrument_impl(impl_variant_id: int):
        indent_lvl = 1
        instr = '#ifdef INSTRUMENT\n'
        instr += indent(indent_lvl) + 'if (me == 0)\n'
        instr += indent(indent_lvl) + '{\n'
        # instr += indent(indent_lvl + 1) + r'printf("\n#ImplVariant-{}\n");'.format('')
        instr += indent(indent_lvl) + '}\n'
        instr += '#endif\n'
        return instr

    @staticmethod
    def write_instrument_kernel_start(indent_lvl: int):
        # Write code to instrument the kernel code.
        instr_start = '#ifdef INSTRUMENT\n'
        instr_start += indent(indent_lvl) + '{\n'
        indent_lvl += 1
        instr_start += indent(indent_lvl) + '#pragma omp barrier\n'
        instr_start += indent(indent_lvl) + 'time_snap_t time;\n'
        instr_start += indent(indent_lvl) + 'if (me == 0)\n'
        instr_start += indent(indent_lvl) + '{\n'
        instr_start += indent(indent_lvl + 1) + 'time_snap_start(&time);\n'
        instr_start += indent(indent_lvl) + '}\n'
        instr_start += '#endif\n'
        return instr_start

    @staticmethod
    def write_instrument_kernel_end(indent_lvl: int, kernel_id: int):
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
                     r'printf("#Kernel={}\§t#Threads=%d\§t%.20e\§n", omp_get_num_threads(), T / 1e9 / n);\n'.format(
                         kernel_id)
        instr_end += '#else\n'
        instr_end += indent(indent_lvl) + \
                     r'printf("#Kernel={}\§t#Threads=1\§t%.20e\§n", T / 1e9 / n);\n'.format(kernel_id)
        instr_end += '#endif\n'
        indent_lvl -= 1
        instr_end += indent(indent_lvl) + '}\n'
        indent_lvl -= 1
        instr_end += indent(indent_lvl) + '}\n'
        instr_end += '#endif\n'
        return instr_end

    def write_loop(self, loop_var: str, loop_iterations: str, code_str: CodeStr):
        """Write for loop with increment 1.

        Parameters:
        -----------
        loop_var : str
            Name of the loop variable.
        loop_iterations : str
            Number of iterations executed.
        code_str : CodeStr
            Container to store the created code.

        Returns:
        --------
        str
            Code string passed with attached for loop code.
        """
        code_str.add(indent(code_str.get_indent()))
        code_str.add(gen_loop_string(self.method, loop_var, 0, loop_iterations))
        code_str.inc_indent()

    @staticmethod
    def write_omp_operation(op_name: str, code_str: CodeStr):
        """Write OpenMP operation.

        Parameters:
        -----------
        op_name : str
            Name of the operation executed.
        code_str : CodeStr
            Container to store the created code.

        Returns:
        --------
        -
        """
        if op_name == 'barrier':
            code_str.add('#pragma omp barrier\n')
        else:
            assert False

    @staticmethod
    def write_pthreads_operation(op_name: str, code_str: CodeStr):
        """Write pthreads operation.

        Parameters:
        -----------
        op_name : str
            Name of the operation executed.
        code_str : CodeStr
            Container to store the created code.

        Returns:
        --------
        -
        """
        assert False

    def write_kernel_template_stub(self, template_name: str, args, code_str: CodeStr):
        """Write kernel template stub.

        Parameters:
        -----------
        template_name : str
            Name of the called kernel template.
        code_str : CodeStr
            Container to store the created code.

        Returns:
        --------
        -
        """
        # Create unique name for this kernel template stub.
        stub_name = '%{}_{}'.format(template_name, str(self.stub_idx))
        self.stub_idx += 1
        # Determine input vector variable name.
        if args:
            solution_vector = str(args[0])
        elif 'RHS_input' in self.loaded_templates[template_name].codegen:
            solution_vector = self.loaded_templates[template_name].codegen['RHS_input']
        else:
            solution_vector = self.config.ode_solution_vector
        # Store kernel stub.
        self.kernel_stubs[stub_name] = (template_name, solution_vector)
        # Write kernel stub code.
        code_str.add(indent(code_str.get_indent()) + '{}\n'.format(stub_name))

    @staticmethod
    def write_double_pointer_swap(ptr_one: str, ptr_two: str, ptr_type: str, code_str: CodeStr):
        """Write pointer swap code.

        Parameters:
        -----------
        ptr_first : str
            Name of the first pointer pointer swapped.
        ptr_second : str
            Name of the first pointer pointer swapped.
        ptr_type : str
            Type of the pointers swapped.
        code_str : CodeStr
            Container to store the created code.

        Returns:
        --------
        -
        """
        code_str.add('\n')
        code_str.add(indent(code_str.get_indent()) + '{} tmp = {};\n'.format(ptr_type, ptr_one))
        code_str.add(indent(code_str.get_indent()) + '{} = {};\n'.format(ptr_one, ptr_two))
        code_str.add(indent(code_str.get_indent()) + '{} = tmp;\n'.format(ptr_two))
        code_str.add('\n')

    def derive_impl_variant_kernels(self, impl_variant: List[str]) -> List[Kernel]:
        return [self.map_kernel_id_to_object(kid) for kid in impl_variant]

    def map_kernel_id_to_object(self, kernel_id: str) -> 'Kernel':
        # Query kernel object from database.
        kernel = Kernel.select(self.db_session, kernel_id)
        # Select corresponding kernel template.
        template = self.loaded_templates[kernel.template.name]
        try:
            return next(filter(lambda x: x.name == kernel.name, template.variants))
        except StopIteration:
            raise RuntimeError('')

    def write_skeleton_includes(self, variant_id: int):
        """Write header includes.

        Parameters:
        -----------
        variant_id : int
            Database ID associated to the generated implementation variant.

        Returns:
        --------
        str
            Generated code.
        """
        codegen = self.skeleton.codegen
        includes = ''
        # Add inclusion guard first.
        includes += '#pragma once\n\n'
        # Add variant database ID as define.
        includes += '#define VARIANT_ID {}\n'.format(variant_id)
        includes += '\n'
        # IVP.
        if self.skeleton.isIVPdependent:
            includes += '#include <math.h>\n'
            includes += '#include "RHS_{}.h"\n'.format(self.ivp.name)
            self.write_rhs_function()
        # ODE method.
        includes += '#include "ODE_{}.h"\n'.format(self.method.name)
        self.write_ode_method()
        # Data structures.
        if 'datastructs' in codegen:
            includes += '#include "{}"\n'.format(codegen['datastructs'])
        # Instrumentation.
        includes += '#ifdef INSTRUMENT\n#include "timesnap.h"\n#endif\n'
        return includes + '\n'

    @staticmethod
    def write_function_description() -> str:
        """Write function description.

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Generated code.
        """
        func_args = ('const int me', 'const int first', 'const int last', 'double t', 'double h')
        header = 'void timestep('
        header += ', '.join(func_args)
        return header + ') {\n'

    @staticmethod
    def write_2d_array(name: str, data: List[List[str]]) -> str:
        """Write 2D double array.

        Parameters:
        -----------
        name : str
            Name of the array.
        data : List of List of str
            Array data.

        Returns:
        --------
        str
            Generated code.
        """
        array_str = 'double ' + name + '[' + str(len(data)) + ']' + '[' + str(len(data[0])) + ']' + ' = {'
        for row in data:
            array_str += '{'
            for elem in row:
                array_str += str(eval_math_expr(elem)) + ', '
            array_str = array_str[:-2]
            array_str += '}, '
        array_str = array_str[:-2]
        array_str += '};\n'
        return array_str

    @staticmethod
    def write_1d_array(name: str, data: List[str]) -> str:
        """Write 1D double array.

        Parameters:
        -----------
        name : str
            Name of the array.
        data : List of str
            Array data.

        Returns:
        --------
        str
            Generated code.
        """
        array_str = 'double ' + name + '[' + str(len(data)) + ']' + ' = {'
        for elem in data:
            array_str += str(eval_math_expr(elem)) + ', '
        array_str = array_str[:-2]
        array_str += '};\n'
        return array_str

    def write_rhs_function(self):
        """Write RHS Function to file.

        Parameters:
        -----------
        -

        Returns:
        --------
        -
        """
        # Create folder if it does not yet exist.
        if self.folder and not self.folder.exists():
            self.folder.mkdir(parents=True)
        # Write inclusion guard first.
        rhs_func_str = '#pragma once\n\n'
        # Write IVP name.
        rhs_func_str += '#define PROBLEM_NAME "{}"\n'.format(self.ivp.name)
        rhs_func_str += '#define PROBLEM_ID {}\n'.format(self.ivp.db_id)
        rhs_func_str += '\n'
        # Write eval_range function.
        rhs_func_str += 'static inline void eval_range('
        rhs_func_str += 'int first, int last, double t, const double *%in, double *f) {\n'
        rhs_func_str += self.ivp.code_eval_range
        rhs_func_str += '}\n'
        # Write eval_component function.
        rhs_func_str += '\nstatic inline double eval_component('
        rhs_func_str += 'int j, double t, const double *%in) {\n'
        rhs_func_str += self.ivp.code_eval_component
        rhs_func_str += '}\n'
        # Write initial solution function.
        rhs_func_str += '\nstatic inline void initial_values('
        rhs_func_str += 'const double t, double *%in) {\n'
        rhs_func_str += self.ivp.code_initial_values
        rhs_func_str += '}\n'
        # Replace solution vector stubs.
        rhs_func_str = sub('%in', self.config.ode_solution_vector, rhs_func_str)
        # Replace IVP constants.
        constants = self.ivp.constants.as_tuple()
        for name, desc in self.ivp.constants.items():
            name_regex = r'(?![a-zA-Z0-9]-_)' + name + r'(?![a-zA-Z0-9-_])'
            rhs_func_str = sub(name_regex, str(eval_math_expr(desc.value, constants)), rhs_func_str)
        # Write RHS function to file.
        path = Path('{}/RHS_{}.h'.format(self.folder, self.ivp.name))
        with path.open('w') as file_handle:
            file_handle.write(rhs_func_str)
        # Format code with indent tool if available.
        format_codefile(path)

    def write_butcher_table(self) -> str:
        """Write butcher table arrays.

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Generated code.
        """
        params_string = self.write_2d_array('A', self.method.coefficientsA)
        params_string += self.write_1d_array('b', self.method.coefficientsB)
        params_string += self.write_1d_array('c', self.method.coefficientsC)
        return params_string + '\n'

    def write_ode_method(self):
        """Write ODE method to file.

        Parameters:
        -----------
        -

        Returns:
        --------
        -
        """
        # Create folder if it does not yet exist.
        if self.folder and not self.folder.exists():
            self.folder.mkdir(parents=True)
        # Write inclusion guard first.
        method_str = '#pragma once\n\n'
        # Write method name and database id.
        method_str += '#define METHOD_NAME "{}"\n'.format(self.method.name)
        method_str += '#define METHOD_ID {}\n'.format(self.method.db_id)
        method_str += '\n'
        # Write method characteristics.
        method_str += '#define s {}\n'.format(self.method.stages)
        method_str += '#define m {}\n'.format(self.method.correctorSteps)
        method_str += '#define ORDER {}\n'.format(self.method.order_)
        method_str += '\n'
        # Write butcher coefficients.
        method_str += self.write_butcher_table()
        # Write ODE method to file.
        path = Path('{}/ODE_{}.h'.format(self.folder, self.method.name))
        with path.open('w') as file_handle:
            file_handle.write(method_str)
        # Format code with indent tool if available.
        format_codefile(path)

    def create_filename(self, variant: List[Kernel]) -> str:
        """Create filename of an implementation variant.

        Parameters:
        -----------
        variant : list of Kernel
            Kernels associated with this implementation variant.

        Returns:
        --------
        str
            Created filename.
        """
        return '{}_{}'.format(self.skeleton.name, '_'.join((kernel.name for kernel in variant)))

    @staticmethod
    def write_to_file(generated_files: Dict[Path, str]):
        """Write generated implementation variant codes to files.

        Parameters:
        -----------
        generated_files : dict (key=Path, value=str)
            Generated implementation variant codes and their associated file paths.

        Returns:
        --------
        -
        """
        for path, code in generated_files.items():
            with path.open('w') as file_handle:
                file_handle.write(code)
            # Replace newline characters in printf commands.
            try:
                cmd = ['sed', '-i', 's/§n/n/g', '{}'.format(path)]
                run(cmd, check=True, stdout=PIPE).stdout
            except CalledProcessError as error:
                raise RuntimeError('Unable to substitute newline character stubs in: {}'.format(path, error))
            # Replace tabulator characters in printf commands.
            try:
                cmd = ['sed', '-i', 's/§t/t/g', '{}'.format(path)]
                run(cmd, check=True, stdout=PIPE).stdout
            except CalledProcessError as error:
                raise RuntimeError('Unable to substitute tabulator character stubs in: {}'.format(path, error))
            # Format code with indent tool if available.
            format_codefile(path)
