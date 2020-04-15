"""@package yasksite_generator
Definition of class YasksiteCodeGenerator.
"""

from pathlib import Path

import attr
from pcre import sub

from offsite.codegen.code_object import CodeStr, YasksiteCode
from offsite.codegen.codegen_util import gen_loop_string, indent
from offsite.codegen.template_generator import TemplateCodeGenerator
from offsite.descriptions.ivp import IVP
from offsite.descriptions.ode_method import ODEMethod
from offsite.descriptions.parser_utils import DatastructType
from offsite.evaluation.math_utils import corrector_steps, stages


@attr.s(kw_only=True)
class YasksiteCodeGenerator(TemplateCodeGenerator):
    """Representation of the code generator for pmodel kernel code.

    Attributes:
    -----------
    -
    """

    def __attrs_post_init__(self):
        """Create this YasksiteCodeGenerator object.

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
        }

    def generate(self, kernel: 'Kernel', method: ODEMethod, ivp: IVP) -> YasksiteCode:
        """Generate yasksite code for a particular Kernel object.

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
        YasksiteCode
            Generated yasksite code.
        """
        # Set some members to match current code generation.
        self.kernel = kernel
        self.template = kernel.template
        self.ivp = ivp
        self.method = method
        # Reset loop unroll stack.
        self.unroll_stack.clear()

        # var_defs = self.construct_variable_defs()
        #
        # Initialize yasksite code object.
        code = YasksiteCode(kernel.name, ivp.name if ivp else '', '')
        code.code_string.set(self.construct_variable_defs())
        # Parse kernel object and create pmodel code.
        self.parse_kernel(code)
        # Insert RHS evaluations.
        if code.contains('%RHS'):
            assert self.ivp
            assert self.ivp.characteristic.isStencil
            self.replace_rhs_tags(code)
        # Further optimize the generated code ...
        # ... unroll loops.
        code.code_string.unroll(self.unroll_stack, self.method)
        # ... insert butcher table coefficients.
        code.code_string.butcher(self.method)
        return code

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
            if desc.isYasksiteParam:
                var_type = 'PARAM'
            elif desc.struct_type == DatastructType.scalar:
                var_type = 'GRID_POINT'
            else:
                # Is spatial grid.
                if 'n' in desc.size:
                    if len(desc.size) == 1:
                        var_type = 'GRID'
                    else:
                        var_type = 'GRID_ARRAY'
                else:
                    var_type = 'GRID_POINT_ARRAY'
            # Write variable definition.
            var_defs += '{} {}{};\n'.format(var_type, name, desc.dimensions_to_string_yasksite())
        # Replace ODE method constants.
        corrector_steps_tup = corrector_steps(self.method)
        if '[{}]'.format(corrector_steps_tup[0]) in var_defs:
            var_defs = sub(r'\[{}\]'.format(corrector_steps_tup[0]), '[{}]'.format(corrector_steps_tup[1]), var_defs)
        stages_tup = stages(self.method)
        if '[{}]'.format(stages_tup[0]) in var_defs:
            var_defs = sub(r'\[{}\]'.format(stages_tup[0]), '[{}]'.format(stages_tup[1]), var_defs)
        return var_defs + '\n'

    def parse_kernel(self, yasksite_code: YasksiteCode):
        """Parse Kernel object and create its corresponding yasksite code objects.

        Parameters:
        -----------
        yasksite_code : YasksiteCode
            Container to store the created yasksite codes.

        Returns:
        --------
        -
        """
        code = yasksite_code.code_string
        # Parse kernel.
        for line in self.kernel.code.split('\n'):
            func_tuple = self.parse_line(line)
            if func_tuple[0]:
                func_tuple[0](code, *func_tuple[1])

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

    def replace_rhs_tags(self, code: YasksiteCode):
        """Replace %RHS tags in code string with the given IVP component.

        Parameters:
        -----------
        code : YasksiteCode
            Container storing the yasksite code string with RHS tags.

        Returns:
        --------
        """
        rhs_code_str = ''
        for line in code.code_string.get().split('\n'):
            if '%RHS' in line:
                input_vec = self.template.codegen['RHS_input']
                # Write IVP constants to std::map.
                constant_map = '{'
                for name, desc in self.ivp.constants.items():
                    constant_map += '{\"' + name + '\", ' + desc.value + '}'
                constant_map += '}'
                repl = 'RHS({}, {}, {}, {})'.format(
                    input_vec, constant_map, Path(self.ivp.code_stencil_path).stem,
                    self.template.codegen['RHS_butcher_nodes'])
                # TODO FIX bug in constant_map creation --> single }
                line = sub(r'%RHS', repl, line)
            rhs_code_str += line + '\n'
        code.code_string.set(rhs_code_str)
