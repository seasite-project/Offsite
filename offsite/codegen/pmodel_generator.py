"""@package pmodel_generator
Definition of class PModelCodeGenerator.
"""

from copy import deepcopy
from typing import List

import attr
from pcre import sub

from offsite.codegen.code_object import CodeStr, PModelCode
from offsite.codegen.codegen_util import gen_loop_string, indent
from offsite.codegen.template_generator import TemplateCodeGenerator
from offsite.descriptions.ivp import IVP
from offsite.descriptions.ode_method import ODEMethod
from offsite.evaluation.math_utils import eval_math_expr


@attr.s(kw_only=True)
class PModelCodeGenerator(TemplateCodeGenerator):
    """Representation of the code generator for pmodel kernel code.

    Attributes:
    -----------
    -
    """

    def __attrs_post_init__(self):
        """Create this PModelCodeGenerator object.

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
            '%COMP': self.resolve_comp,
            '%PMODEL': self.resolve_pmodel
        }

    def generate(self, kernel: 'Kernel', method: ODEMethod, ivp: IVP) -> PModelCode:
        """Generate pmodel code for a particular Kernel object.

        Parameters:
        -----------
        kernel : Kernel
            Kernel for which code is generated.
        method : ODEMethod
            ODE method for which code is generated.
        ivp : IVP
            IVP possibly evaluated by the kernel for which code is generated.

        Returns:
        --------
        PModelCode
            Generated PModel code.
        """
        # Set some members to match current code generation.
        self.kernel = kernel
        self.template = kernel.template
        self.ivp = ivp
        self.method = method
        # Reset loop unroll stack.
        self.unroll_stack.clear()

        # Generate variable definitions.
        var_defs = self.construct_variable_defs()
        # Initialize pmodel code object.
        pmodel_code = PModelCode(kernel.name, ivp.name if ivp else '', var_defs)
        pmodel_code.code_strings = list()
        # Parse kernel object and create pmodel code.
        self.parse_kernel(pmodel_code)
        # Insert RHS evaluations.
        for idx, pmodel_str in enumerate(pmodel_code.code_strings):
            if '%RHS' in pmodel_str.get():
                assert self.ivp
                pmodel_code.code_strings[idx] = self.replace_rhs_tags(pmodel_str)
        # Further optimize the generated code ...
        for pmodel_str in pmodel_code.code_strings:
            if isinstance(pmodel_str, list):  # Multiple RHS components.
                for pmodel_rhs_str in pmodel_str:
                    # ... unroll loops.
                    pmodel_rhs_str.unroll(self.unroll_stack, self.method)
                    # ... insert butcher table coefficients.
                    pmodel_rhs_str.butcher(self.method)
            else:  # None or a single RHS component.
                # ... unroll loops.
                pmodel_str.unroll(self.unroll_stack, self.method)
                # ... insert butcher table coefficients.
                pmodel_str.butcher(self.method)
        return pmodel_code

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
        var_defs = ''
        for name, desc in self.template.datastructs.items():
            var_defs += '{} {}{};\n'.format(desc.datatype, name, desc.dimensions_to_string())
        return var_defs + '\n'

    def parse_kernel(self, pmodel_code: 'PModelCode'):
        """Parse Kernel object and create its corresponding pmodel code objects.

        Parameters:
        -----------
        pmodel_code : PModelCode
            Container to store the created pmodel codes.

        Returns:
        --------
        -
        """
        # Initialize CodeStr object for next pmodel code.
        code = CodeStr('', 0)
        # Parse kernel.
        for line in self.kernel.code.split('\n'):
            func_tuple = self.parse_line(line)
            if func_tuple[0]:
                evaluated_func = func_tuple[0](code, *func_tuple[1])
                # Kernel requires multiple pmodel codes.
                if isinstance(evaluated_func, CodeStr):
                    pmodel_code.add_code_str(evaluated_func)
                    # Switch to next pmodel code.
                    code.reset()
                    for loop in self.loop_stack:
                        self.write_loop(loop[0], loop[1], code)
        pmodel_code.add_code_str(code)

    def resolve_pmodel(self, code_str: 'CodeStr', *args) -> 'CodeStr':
        """Resolve pmodel code line.

        Parameters:
        -----------
        code_str : CodeStr
            Container to store the created code.
        args
            Required function parameters.

        Returns:
        --------
        CodeStr:
            Copy of the passed code string.
        """
        self.write_closing_brackets(len(self.loop_stack), code_str)
        return deepcopy(code_str)

    def write_loop(self, loop_var: str, loop_iterations: str, code_str: 'CodeStr'):
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

    def replace_rhs_tags(self, code_str: 'CodeStr') -> List['CodeStr']:
        """Replace %RHS tags in code string with the given IVP component(s).

        If the IVP contains multiple components an individual code string is created for each component.

        Parameters:
        -----------
        code_str : CodeStr
            Container storing the pmodel code string with RHS tags.

        Returns:
        --------
        List of CodeStr
            Container that stores the pmodel code string(s) with replaced RHS tags.
        """
        rhs_code_strs = list()
        for name, component in self.ivp.components.items():
            replaced_code_str = ''
            # Write RHS term.
            rhs_term = self.write_rhs_term(component)
            # Replace RHS tags.
            for line in code_str.get().split('\n'):
                if '%RHS' in line:
                    line = sub('%RHS', rhs_term, line)
                replaced_code_str += line + '\n'
            pmodel_cpy = deepcopy(code_str)
            pmodel_cpy.set(replaced_code_str)
            rhs_code_strs.append(pmodel_cpy)
        return rhs_code_strs

    def write_rhs_term(self, component: 'IVPComponent') -> str:
        """Write RHS term.

        Parameters:
        -----------
        component : IVPComponent
            Evaluated IVP component.

        Returns:
        --------
        str
            Code string generated.
        """
        # Collect IVP constants and its values.
        constants = self.ivp.constants.as_tuple()
        # Generate RHS term.
        rhs_term = ''
        for line in component.code.split('\n'):
            # Replace constants.
            for name, desc in self.ivp.constants.items():
                name_regex = r'(?![a-zA-Z0-9]-_)' + name + r'(?![a-zA-Z0-9-_])'
                line = sub(name_regex, str(eval_math_expr(desc.value, constants)), line)
            rhs_term += line + '\n'
        rhs_term = rhs_term.rstrip()
        # Substitute %in with ODE solution vector variable.
        rhs_term = sub("%in", self.config.ode_solution_vector, rhs_term)
        # Replaces indices.
        rhs_term = sub('%idx', self.config.var_idx, rhs_term)
        rhs_term = sub('%last_idx', self.config.var_last_idx, rhs_term)
        return rhs_term
