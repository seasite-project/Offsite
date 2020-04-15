"""@package codegen_utils.py
Util functions, classes used during code generation.
"""

from enum import Enum
from pathlib import Path
from shutil import which
from subprocess import run, PIPE, CalledProcessError

from pcre import sub

from offsite.descriptions.ode_method import ODEMethod
from offsite.evaluation.math_utils import eval_math_expr, corrector_steps, stages


class GeneratorType(Enum):
    """Defines what type of code is generated.

    - PMODEL
        Code of a kernel used during performance prediction.
    - IMPL
        Code of an implementation variant.
    - METHOD
        Code of a kernel.
    """
    PMODEL = 'pmodel'
    IMPL = 'impl'
    METHOD = 'method'


def escape_string(string):
    escaped_string = ''
    for char in string:
        if char in '[]+-)(*{}':
            escaped_string += '\\'
        escaped_string += char
    return escaped_string


def replace_var_with_factor(var: str, computation: str, factor: str):
    """

    Parameters:
    -----------
    var : str
        variable as string.
    computation : str
        computation string
    factor : str
        factor to replace var

    Returns:
    --------
    str
        Computation string with replaced var.
    """
    new_computation = sub(' ' + var + ';', ' ' + str(factor) + ';', computation)
    return sub(r'\[' + var + r'\]', '[' + str(factor) + ']', new_computation)


def indent(level: int):
    """Create indentation string for formatting of the generated code.

    Parameters:
    -----------
    level : int
        Level of the indentation.

    Returns:
    --------
    str
        Indentation.
    """
    indention = '  '
    return indention * level if level > 0 else ''


def gen_loop_string(method, loop_var, start_idx, end_idx, comparator: str = '<'):
    loop_str = 'for(int ' + loop_var + ' = ' + str(start_idx) + '; ' + loop_var + ' ' + comparator + ' ' + \
               str(eval_iterations(method, end_idx)) + '; ++' + loop_var + ') {\n'
    return loop_str


def eval_iterations(method: ODEMethod, iterations: str):
    constants = [corrector_steps(method), stages(method)]
    iterations = eval_math_expr(str(iterations), constants)
    try:
        iterations = int(iterations)
    except TypeError:
        pass
    return str(iterations)


def cut_string_edges(pre_char, post_char, string):
    pre_splitted = pre_char.join(string.split(pre_char)[1:])
    post_splitted = post_char.join(pre_splitted.split(post_char)[:-1])
    return post_splitted


def format_codefile(path: Path):
    """Format given code file with command line tool indent.

    Parameters:
    -----------
    path : Path
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
            raise RuntimeError('Formatting file {} with indent failed: {}'.format(path, error))
