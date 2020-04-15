"""@package code_generator
Definition of class CodeGenerator.
"""

from abc import abstractmethod
from typing import Dict, Tuple

import attr

from offsite.codegen.code_object import CodeStr
from offsite.codegen.codegen_util import indent

# Indention used for code formatting.
INDENTION = '  '


@attr.s(kw_only=True)
class CodeGenerator:
    """Representation of the code generator base class.

    Attributes:
    -----------
    config : Config
        Program configuration.
    parser : Dict of functions
        Parser functions used to resolve code tags.
    ivp : IVP
        IVP for which code is generated. Set by generate() function.
    method : ODEMethod
        ODEMethod for which code is generated. Set by generate() function.
    """
    config = attr.ib(type='Config')
    parser = attr.ib(type=Dict, init=False)
    ivp = attr.ib(type='IVP', default=None, init=False)
    method = attr.ib(type='ODEMethod', default=None, init=False)

    def parse_line(self, line: str) -> Tuple[str, str]:
        """Parse code line and return its information as tuple.

        Parameters:
        -----------
        line : str
            Parsed code line.

        Returns:
        --------
        tuple of str
            Function executed in this code line and its parameters.
        """
        line = line.strip()
        line_split = line.split(' ')
        command = line_split[0]
        func_params = line_split[1:]
        func_to_token = self.parser.get(command)
        try:
            return func_to_token, func_params
        except:
            return None

    @abstractmethod
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

    @abstractmethod
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

    @abstractmethod
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

    def write_closing_brackets(self, open_brackets: int, code_str: CodeStr):
        """Write closing brackets.

        Parameters:
        -----------
        open_brackets: int
            Number of brackets to be closed.
        code_str : CodeStr
            Container to store the created code.

        Returns:
        --------
        -
        """
        code_str.dec_indent()
        bracket_string = '{}}}\n'.format(indent(code_str.get_indent())) * open_brackets
        code_str.add(bracket_string)
