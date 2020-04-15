"""@package pmodel_code
Definitions of classes Code, CodeStr, PModelCode, BenchCode, KernelCode, YasksiteCode.
"""

from abc import abstractmethod
from pathlib import Path
from typing import List

import attr
from pcre import compile as pcre_compile, findall, M
from regex import sub
from sortedcontainers import SortedDict

from offsite.codegen.codegen_util import escape_string, replace_var_with_factor, cut_string_edges, indent, format_codefile
from offsite.descriptions.ode_method import ODEMethod
from offsite.evaluation.math_utils import eval_math_expr


@attr.s
class CodeStr:
    """Representation of a code string object.

    Attributes:
    -----------
    content : list of str
        Code lines of this object.
    indent : int
        Indent level of the latest code string line.
    loop_stack : List[str]
        Stack of all loops of the code string that should be unrolled.
    """
    # code_lines = attr.ib(type=List[str])
    string = attr.ib(type=str)
    indent = attr.ib(type=int)
    loop_stack = attr.ib(type=List[str], default=list())

    def unroll(self, unroll_stack: SortedDict, method: ODEMethod):
        #
        for loops_on_level in unroll_stack.values():
            # Unroll loops.
            for loop in loops_on_level:
                idx = loop[0]
                iterations = loop[1]
                assign = loop[2]
                #
                loop_nests = self.find_loop_nest(idx)
                for loop_nest in loop_nests:
                    loop = loop_nest.group(0)
                    # Read loop boundary.
                    boundary = method.stages if iterations == 's' else 0
                    # Determine body of the outermost loop (without curly brackets).
                    loop_split = loop.split('\n')
                    comp_str = ''
                    for elem in range(1, len(loop_split) - 1):
                        comp_str += loop_split[elem] + '\n'
                    # Determine start index of outermost loop.
                    start_idx = int(findall(r'\d+', findall(r'for\([\w\s\=<>]*', str(loop_split))[0])[0])
                    # Unroll body of the outermost loop.
                    unrolled_comp = self.construct_unrolled_loop(start_idx, idx, boundary, comp_str, assign)
                    #
                    comp_str = escape_string(comp_str)
                    #
                    self.string = sub(r'[\t]*for\(int ' + idx + r'[\+\-\(\)<>;\s\=\w]*\{[\s]*' + comp_str + r'[\s]*}',
                                      unrolled_comp, self.string)

    def find_loop_nest(self, idx: str):
        pattern = r'[\t]*for\(int ' + r'{}'.format(idx) + r'[\s]+' + \
                  r'[\w\s;=<>\+\-\(\)]* (?<rec>\{(?:[^\{\}]+|(?&rec))*\})'
        return pcre_compile(pattern, M).finditer(self.string)

    def split_loop(self, iteration_count, idx):
        splitted_loop = ''
        for loop_nest in self.find_loop_nest(idx):
            loop_nest = loop_nest.group(0)
            escaped_loop = escape_string(loop_nest)
            loop_content = cut_string_edges('{', '}', loop_nest)
            loop_header = loop_nest.split('{')[0]
            for i in range(0, iteration_count):
                replaced_loop_content = replace_var_with_factor(idx, loop_content, i)
                if i == 0:
                    replaced_loop_content = self.initialize_loop_content(replaced_loop_content)
                splitted_loop += replaced_loop_content
            splitted_loop += replace_var_with_factor('0', loop_header, iteration_count) + '{\n' + loop_content + '}'
            self.string = sub(escaped_loop, splitted_loop, self.string)

    def initialize_loop_content(self, loop_content):
        loop_content = sub(r'(\+=|\-=)', '=', loop_content)
        return loop_content

    def construct_unrolled_loop(self, start_idx, loop_var, steps, computation, assign):
        computation = computation.strip()
        new_computation = ''
        for iteration in range(start_idx, steps):
            replaced = replace_var_with_factor(loop_var, computation, iteration)
            if iteration == assign:
                replaced = self.initialize_loop_content(replaced)
            if iteration != start_idx:
                replaced = indent(self.indent + 1) + replaced
            if iteration != (steps - 1):
                replaced += '\n'
            new_computation += replaced
        return new_computation

    def inc_indent(self):
        self.indent += 1

    def dec_indent(self):
        self.indent -= 1

    def add(self, string):
        self.string += string

    def prepend(self, string):
        self.string = string + self.string

    def get_indent(self):
        return self.indent

    def get(self):
        return self.string

    def set(self, string):
        self.string = string

    def __add__(self, another_variant_str):
        self.string += another_variant_str.string

    def reset(self):
        self.string = ''
        self.indent = 0

    def butcher(self, method: ODEMethod):
        """Replace butcher table entries in the code string with its actual values as given by the ODE method passed.

         Parameters:
         -----------
         method : ODEMethod
             ODE method used.

         Returns:
         --------
         -
         """
        # Replace A coefficients.
        for idx_row, row_a in enumerate(method.coefficientsA):
            for idx_col, a in enumerate(row_a):
                self.string = self.string.replace('A[{}][{}]'.format(idx_row, idx_col), eval_math_expr(a, cast_to=str))
        # Replace b coefficients.
        for idx, b in enumerate(method.coefficientsB):
            self.string = self.string.replace('b[{}]'.format(idx), eval_math_expr(b, cast_to=str))
        # Replace c coefficients.
        for idx, c in enumerate(method.coefficientsC):
            self.string = self.string.replace('c[{}]'.format(idx), eval_math_expr(c, cast_to=str))


@attr.s
class Code:
    """Representation of a code object.

    Attributes:
    -----------
    name : str
        Name of this object.
    ivp_name : str
        Name of the included IVP.
    frame_code : str
        Frame code that is prepended to each single code string.
    """
    name = attr.ib(type=str)
    ivp_name = attr.ib(type=str)
    frame_code = attr.ib(type=str)

    @abstractmethod
    def add_code_str(self, code_str: CodeStr):
        """Add an additional code string to this code object.

        Parameters:
        -----------
        code_str : CodeString
            Code string representation.

        Returns:
        --------
        -
        """

    @abstractmethod
    def write_to_file(self, folder: Path = None):
        """Write code to file(s) and store.

        Parameters:
        -----------
        folder : Path
            Relative path to the location of the code files generated.

        Returns:
        --------
        -
        """

    @abstractmethod
    def create_filename(self) -> str:
        """Create unique filename of a particular code file.

        Parameters:
        -----------
        -

        Returns:
        --------
        -
        """


class PModelCode(Code):
    """Representation of a PModelCode object.

    Attributes:
    -----------
    code_strings : CodeStr
        String representations of the code of this object.
    paths : list of Paths
        Relative paths of the code files generated from the code string.
    """
    code_strings = attr.ib(type=List[CodeStr], default=[])
    paths = attr.ib(type=List[Path], default=[])

    def add_code_str(self, code_str: CodeStr):
        """Add an additional code string to this code object.

        Parameters:
        -----------
        code_str : CodeString
            Code string representation.

        Returns:
        --------
        -
        """
        code_str.prepend(self.frame_code)
        self.code_strings.append(code_str)

    def write_to_file(self, folder: Path = None):
        """Write pmodel codes to file.

        Parameters:
        -----------
        folder : Path
            Relative path to folder that will contain the code files written.

        Returns:
        --------
        -
        """
        # Create folder if it does not yet exist.
        if folder and not folder.exists():
            folder.mkdir(parents=True)
        # Write code.
        self.paths = list()
        for pid, code_str in enumerate(self.code_strings):
            if isinstance(code_str, list):
                for rid, pmodel_rhs_str in enumerate(code_str):
                    filename = self.create_filename(pid, rid)
                    path = Path('{}.c'.format(filename))
                    if folder:
                        path = folder / path
                    with path.open('w') as file_handle:
                        file_handle.write(pmodel_rhs_str.get())
                    self.paths.append(path)
            else:
                filename = self.create_filename(pid)
                path = Path('{}.c'.format(filename))
                if folder:
                    path = folder / path
                with path.open('w') as file_handle:
                    file_handle.write(code_str.get())
                self.paths.append(path)
        # Format code files.
        for path in self.paths:
            format_codefile(path)

    def create_filename(self, pid: int, rid: int = -1) -> str:
        """Create unique filename of a particular code file.

        Parameters:
        -----------
        pid : int
            ID of the pmodel kernel used.
        rid : int
            ID of the IVP component used.

        Returns:
        --------
        -
        """
        name_postfix = '_pmodel'
        if len(self.code_strings) > 1:
            name_postfix += '_' + str(pid + 1)
        if rid > -1:
            name_postfix += '_' + self.ivp_name
            if len(self.code_strings[pid]) > 1:
                name_postfix += '_' + str(rid + 1)
        return self.name + name_postfix


class BenchCode(PModelCode):
    """Representation of a BenchCode object.

    Attributes:
    -----------
    -
    """

    def write_to_file(self, folder: Path = None):
        """Write pmodel codes to file.

        Parameters:
        -----------
        folder : pathlib.Path
            Relative path to folder that will contain the code files written.

        Returns:
        --------
        -
        """
        self.paths = list()
        for pid, code_string in enumerate(self.code_strings):
            if isinstance(code_string, list):
                for rid, rhs_bench_str in enumerate(code_string):
                    filename = self.create_filename(pid, rid)
                    path = Path('{}.c'.format(filename))
                    if folder:
                        path = folder / path
                    with path.open('w') as file_handle:
                        file_handle.write(rhs_bench_str.get())
                    self.paths.append(path)
            else:
                filename = self.create_filename(pid)
                path = Path('{}.c'.format(filename))
                if folder:
                    path = folder / path
                with path.open('w') as file_handle:
                    file_handle.write(code_string.get())
                self.paths.append(path)
        # Format code files.
        for path in self.paths:
            format_codefile(path)


@attr.s
class KernelCode(Code):
    """Representation of a KernelCode object.

    Attributes:
    -----------
    code_string : CodeStr
        String representation of the code of this object.
    path : Paths
        Relative path of the code file generated from the code string.
    """
    code_string = attr.ib(type=CodeStr, default=None)
    path = attr.ib(type=Path, default=None)

    def add_code_str(self, code_str: CodeStr):
        """Add a code string to this code object.

        Parameters:
        -----------
        code_str : CodeString
            Code string representation.

        Returns:
        --------
        -
        """
        code_str.prepend(self.frame_code)
        self.code_string = code_str

    def write_to_file(self, folder: Path = None):
        """Write kernel code to file.

        Parameters:
        -----------
        folder : Path
            Relative path to folder that will contain the code files written.

        Returns:
        --------
        -
        """
        filename = self.create_filename()
        path = Path('{}.h'.format(filename))
        if folder:
            path = folder / path
        with path.open('w') as file_handle:
            file_handle.write(self.code_string.get())
        self.path = path
        # Format code file.
        format_codefile(self.path)

    def create_filename(self, rid: int = -1) -> str:
        """Create unique filename of a particular code file.

        Parameters:
        -----------
        rid : int
            ID of the IVP component used.

        Returns:
        --------
        -
        """
        name_postfix = ''
        if rid > -1:
            name_postfix += '_' + self.ivp_name
        return self.name + name_postfix


@attr.s
class YasksiteCode(Code):
    """Representation of a YasksiteCode object.

    Attributes:
    -----------
    code_string : CodeStr
        String representation of the code of this object.
    path : Paths
        Relative path of the code file generated from the code string.
    """
    code_string = attr.ib(type=CodeStr, default=CodeStr('', 0))
    path = attr.ib(type=Path, default=Path())

    def add_code_str(self, code_str: CodeStr):
        """Add a code string to this code object.

        Parameters:
        -----------
        code_str : CodeString
            Code string representation.

        Returns:
        --------
        -
        """
        code_str.prepend(self.frame_code)
        self.code_string = code_str

    def write_to_file(self, folder: Path = None):
        """Write yasksite code to file.

        Parameters:
        -----------
        folder : Path
            Relative path to folder that will contain the code files written.

        Returns:
        --------
        -
        """
        # Create folder if it does not yet exist.
        if folder and not folder.exists():
            folder.mkdir(parents=True)
        # Write to file.
        filename = self.create_filename()
        path = Path('{}.c'.format(filename))
        if folder:
            path = folder / path
        with path.open('w') as file_handle:
            file_handle.write(self.code_string.get())
        self.path = path

    def create_filename(self) -> str:
        """Create unique filename of a particular code file.

        Parameters:
        -----------
        -

        Returns:
        --------
        -
        """
        name_postfix = '_yasksite'
        if self.ivp_name:
            name_postfix += '_' + self.ivp_name
        return self.name + name_postfix

    def contains(self, substr: str) -> bool:
        """Check if code string contains a particular substring.

        Parameters:
        -----------
        substr : str
            Checked substring.

        Returns:
        --------
        Bool
            True if code contains the substring else False.
        """
        return substr in self.code_string.get()
