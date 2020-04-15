"""@package kernel_generator
Definition of class KernelCodeGenerator.
"""

from shlex import split as shlex_split
from typing import List, Tuple

import attr
from pcre import sub
from sortedcontainers import SortedDict

from offsite.codegen.code_generator import CodeGenerator
from offsite.codegen.code_object import CodeStr, KernelCode
from offsite.codegen.codegen_util import gen_loop_string, indent
from offsite.descriptions.ivp import IVP
from offsite.descriptions.kernel_template import Kernel
from offsite.descriptions.ode_method import ODEMethod


@attr.s(kw_only=True)
class KernelCodeGenerator(CodeGenerator):
    """Representation of the code generator for kernel code.

    Attributes:
    -----------
    loop_stack : List of tuple of str
        Stack of opened loops.
    unroll_stack : SortedDict
        Stack of loops to be unrolled.
    kernel : Kernel
        Kernel for which code is generated. Set by generate() function.
    template : KernelTemplate
        KernelTemplate to which the kernel belongs for which code is generated. Set by generate() function.
    db_session : sqlalchemy.orm.session.Session
        Used database session.
    """
    loop_stack = attr.ib(type=List[Tuple[str, str]], default=list())
    unroll_stack = attr.ib(type=SortedDict, default=SortedDict())
    kernel = attr.ib(type='Kernel', default=None, init=False)
    template = attr.ib(type='KernelTemplate', default=None, init=False)
    db_session = attr.ib(type='sqlalchemy.orm.session.Session')

    def __attrs_post_init__(self):
        """Create this KernelCodeGenerator object.

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
            '%PRAGMA': self.resolve_pragma,
            '%COMP': self.resolve_comp
        }

    def generate(self, kernel: Kernel, method: ODEMethod, ivp: IVP, input_vector: str) -> KernelCode:
        """Generate kernel code for a particular Kernel object.

        Parameters:
        -----------
        kernel : Kernel
            Representation of a Kernel object.
        method : ODEMethod
            ODE method for which code is generated.
        ivp : IVP
            IVP possibly evaluated by the kernel for which code is generated.
        input_vector : str
            Name of the used input vector.

        Returns:
        --------
        KernelCode
            Generated Kernel code.
        """
        # Set some members to match current code generation.
        self.kernel = kernel
        self.template = kernel.template
        self.ivp = ivp
        self.method = method
        # Reset loop unroll stack.
        self.unroll_stack.clear()

        # Initialize kernel code object.
        kernel_code = KernelCode(kernel.name, ivp.name if ivp else '', '')
        # Parse kernel object and create kernel code.
        self.parse_kernel(kernel_code)
        # Insert RHS evaluations.
        if self.template.isIVPdependent:
            if 'RHS' in kernel_code.code_string.get():
                kernel_code.code_string = self.replace_rhs_tags(input_vector, kernel_code.code_string)
        # Apply required loop splits.
        if 'loop splits' in self.template.codegen:
            for loop_split in self.template.codegen['loop splits']:
                splitted_loop = loop_split.split(' ')
                kernel_code.code_string.split_loop(int(splitted_loop[1]) + 1, splitted_loop[0])
        # Further optimize the generated code ...
        # ... unroll loops.
        kernel_code.code_string.unroll(self.unroll_stack, self.method)
        # ... insert butcher table coefficients.
        kernel_code.code_string.butcher(self.method)
        return kernel_code

    def parse_kernel(self, kernel_code: KernelCode):
        """Parse Kernel object and create its corresponding kernel code object.

        Parameters:
        -----------
        pmodel_code : PModelCode
            Container to store the created kernel code.

        Returns:
        --------
        -
        """
        # Initialize CodeStr object for kernel code.
        code = CodeStr('', 1)
        # Parse kernel.
        for line in self.kernel.code.split('\n'):
            func_tuple = self.parse_line(line)
            if func_tuple[0]:
                func_tuple[0](code, *func_tuple[1])
        code.set(self.replace_indices(code.string))
        kernel_code.code_string = code

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
        # Add loop to stack.
        self.loop_stack.append((loop_var, loop_iterations))
        # Write loop header code.
        self.write_loop(loop_var, loop_iterations, code_str)
        if len(args) > 2:
            # Unroll loop body.
            if args[2] == 'unroll':
                lvl = 0
                try:
                    lvl = int(args[3])
                except IndexError:
                    pass
                # Replace loop iteration by assignment statement.
                assign = None
                if len(args) > 4 and args[4] == 'assign':
                    try:
                        assign = int(args[5])
                    except IndexError:
                        pass
                self.push_unroll_stack(lvl, loop_var, loop_iterations, assign)

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
        self.loop_stack.pop()

    def resolve_pragma(self, code_str: CodeStr, *args):
        """Resolve pragma code line.

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
        code_str.add('#pragma {}\n'.format(args[0]))

    def resolve_comp(self, code_str: CodeStr, *args):
        """Resolve computation code line.

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
        self.write_computation(args[0], code_str)

    def push_unroll_stack(self, lvl: int, loop_name: str, loop_iterations: str, assign: int = None):
        """Push loop to unroll stack.

        Parameters:
        -----------
        lvl : int
            Used to define the order in which loops are unrolled.
        loop_var : str
            Name of the loop variable.
        loop_iterations : str
            Number of iterations executed.
        assign : int
            Replace this loop iteration with an assignment statement.

        Returns:
        --------
        -
        """
        if lvl not in self.unroll_stack:
            self.unroll_stack[lvl] = set()
        self.unroll_stack[lvl].add((loop_name, loop_iterations, assign))

    def write_computation(self, cid: str, code_str: CodeStr):
        """Write computation.

        Parameters:
        -----------
        cid : str
            Identifier of the computation used.
        code_str : CodeStr
            Container to store the created code.

        Returns:
        --------
        -
        """
        comp_str = indent(code_str.get_indent()) + self.template.computations.get_computation(cid).computation + ';\n'
        code_str.add(comp_str)

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
        # Adjust loop borders and loop comparator for loops that run over the ODE system.
        if loop_var == self.config.var_idx:
            first_iter = self.config.var_first_idx
            last_iter = self.config.var_last_idx
            comparator = '<='
        else:
            first_iter = 0
            last_iter = loop_iterations
            comparator = '<'
        # Write loop code.
        code_str.add(indent((code_str.get_indent())))
        code_str.add(gen_loop_string(self.method, loop_var, first_iter, last_iter, comparator))
        # Increment code indent.
        code_str.inc_indent()

    def replace_indices(self, string: str) -> str:
        """Replace the template specific indices in string.

        Parameters:
        -----------
        string : str
            String with not replaced indices.

        Returns:
        --------
        str
            String with replaced indices.
        """
        string = sub('%idx', self.config.var_idx, string)
        string = sub('%last_idx', self.config.var_last_idx, string)
        return string

    def replace_rhs_tags(self, input_vector: str, code_str: CodeStr) -> CodeStr:
        """Replace %RHS tags in code string with the given IVP evaluation function.

        Parameters:
        -----------
        input_vector : str
            Name of the used input vector.
        code_str : CodeStr
            Container storing the kernel code string with RHS tags.

        Returns:
        --------
        CodeStr
            Container that stores the kernel code string with replaced RHS tags.
        """
        #
        func_type = ''.join(self.template.codegen['RHS'])
        skip_next = False
        #
        replaced_lines = list()
        for line in code_str.get().split('\n'):
            if skip_next is True:
                skip_next = False
                continue
            elif '%RHS' in line:
                # Write RHS term.
                rhs_term = indent(code_str.get_indent() + 1)
                if func_type == 'eval_range':
                    replaced_lines = replaced_lines[:-1]
                    #
                    elem_lhs = shlex_split(line)[0]
                    elem_lhs = sub(r'\[j\]', '', elem_lhs)
                    #
                    rhs_term += 'eval_range(first, last, t + {} * h, {}, {})'.format(
                        self.template.codegen['RHS_butcher_nodes'], input_vector, elem_lhs)
                    #
                    line = line.split('=')[1].strip()
                    #
                    skip_next = True
                elif func_type == 'eval_component':
                    rhs_term += 'eval_component({}, t + {} * h, {})'.format(
                        self.config.var_idx, self.template.codegen['RHS_butcher_nodes'], input_vector)
                else:
                    assert False
                # Replace RHS tags.
                line = sub('%RHS', rhs_term, line)
            replaced_lines.append(line + '\n')
        code_str.set(''.join(replaced_lines))
        return code_str
