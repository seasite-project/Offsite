"""@package pmodel_generator
Definition of class TemplateCodeGenerator.
"""

from abc import abstractmethod
from typing import List, Tuple

import attr
from sortedcontainers import SortedDict

from offsite.codegen.code_generator import CodeGenerator
from offsite.codegen.code_object import CodeStr
from offsite.codegen.codegen_util import indent


@attr.s(kw_only=True)
class TemplateCodeGenerator(CodeGenerator):
    """Representation of the code generator for kernel template code.

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
    """
    loop_stack = attr.ib(type=List[Tuple[str, str]], default=list())
    unroll_stack = attr.ib(type=SortedDict, default=SortedDict())
    kernel = attr.ib(type='Kernel', default=None, init=False)
    template = attr.ib(type='KernelTemplate', default=None, init=False)

    @abstractmethod
    def construct_variable_defs(self) -> str:
        """Construct variable definitions of a particular kernel.

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Variable definitions string.
        """

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
