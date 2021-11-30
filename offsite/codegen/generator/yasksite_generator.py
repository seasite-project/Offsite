"""@package codegen.generator.yasksite_generator
Definition of class YasksiteCodeGenerator.
"""

from copy import deepcopy
from typing import Dict, List, Optional, Tuple, Union

import attr
from sortedcontainers import SortedDict

from offsite.codegen.code_dsl.code_node import CodeNode, CodeNodeType, ComputationNode, LoopNode
from offsite.codegen.code_dsl.code_tree import CodeTree
from offsite.codegen.codegen_util import eval_loop_boundary, substitute_stencil_call
from offsite.descriptions.ode import IVP, ODEMethod, corrector_steps, stages
from offsite.descriptions.parser import DatastructDict, DatastructType


@attr.s
class YasksiteCodeGenerator:
    """Representation of the code generator for Yasksite kernel code.

    Attributes:
    -----------
    unroll_stack: SortedDict
        Stack of loops to be unrolled.
    """
    unroll_stack = attr.ib(type=SortedDict, default=SortedDict())

    def collect_loop_meta_data(self, node: CodeNode):
        """Traverse tree depth first and collect all unroll tasks.

        Parameters:
        -----------
        node: CodeNode
            Root node of the code tree.

        Returns:
        --------
        -
        """
        if node.type == CodeNodeType.LOOP:
            node: LoopNode
            if node.optimize_unroll is not None:
                self.unroll_stack[node.optimize_unroll] = node
        # Traverse tree depth first.
        if node.child:
            self.collect_loop_meta_data(node.child)
        if node.next:
            self.collect_loop_meta_data(node.next)

    def unroll_tree(self, tree: CodeNode):
        """Apply all applicable loop unroll optimizations to the passed code tree.

        Parameters:
        -----------
        tree: CodeNode
            Root node of the code tree.

        Returns:
        --------
        -
        """
        while True:
            self.unroll_stack.clear()
            self.collect_loop_meta_data(tree)
            if len(self.unroll_stack.keys()) == 0:
                break
            CodeTree.unroll_loop(self.unroll_stack.values()[0])

    def generate(self, kernel: 'Kernel', method: ODEMethod, ivp: IVP) -> Dict[str, str]:
        """Generate Yasksite code for a particular Kernel object.

        Parameters:
        -----------
        kernel: Kernel
            Kernel for which code is generated.
        method: ODEMethod
            ODE method for which code is generated.
        ivp: IVP
            IVP possibly evaluated by the kernel for which code is generated.

        Returns:
        --------
        dict of str
            Generated Yasksite codes.
        """
        # Set some members to match current code generation.
        code_tree: CodeNode = deepcopy(kernel.code_tree).root
        code_tree.name = kernel.name

        self.unroll_stack.clear()

        # Generate Yasksite code tree.
        self.generate_yasksite_tree(code_tree, kernel.template, method, ivp)
        # Collect loop data.
        self.collect_loop_meta_data(code_tree)
        # Optimize tree: unroll
        self.unroll_tree(code_tree)
        # Substitute butcher table coefficients.
        CodeTree.substitute_butcher_coefficients(
            code_tree, method.coefficientsA, method.coefficientsB, method.coefficientsC)
        # Generate code from tree and write to string.
        codes: Dict[str, str] = dict()
        if kernel.template.isIVPdependent:
            assert ivp is not None
            assert ivp.characteristic.isStencil
            code_name: str = code_tree.name + '_{}'.format(ivp.name)
            codes[code_name] = self.generate_yasksite_code(code_tree)
        else:
            codes[code_tree.name] = self.generate_yasksite_code(code_tree)
        # Prepend variable definitions to written code.
        var_defs: str = self.construct_variable_defs(kernel.template.datastructs, method)
        codes: Dict[str, str] = {name: var_defs + code for name, code in codes.items()}
        return codes

    @staticmethod
    def construct_variable_defs(datastructs: DatastructDict, method: ODEMethod) -> str:
        """Construct variable definitions.

        Parameters:
        -----------
        datastructs: DatastructDict
            Used data structures.
        method: ODEMethod
            Used ODE method.

        Returns:
        --------
        str
            Variable definitions string.
        """
        var_defs = ''
        for name, desc in datastructs.items():
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
        if method is not None:
            corrector_steps_tup: Tuple[str, str] = corrector_steps(method)
            if '[{}]'.format(corrector_steps_tup[0]) in var_defs:
                var_defs: str = var_defs.replace('[{}]'.format(corrector_steps_tup[0]),
                                                 '[{}]'.format(corrector_steps_tup[1]))
            stages_tup = stages(method)
            if '[{}]'.format(stages_tup[0]) in var_defs:
                var_defs: str = var_defs.replace('[{}]'.format(stages_tup[0]), '[{}]'.format(stages_tup[1]))
        return var_defs + '\n'

    @staticmethod
    def generate_yasksite_tree(node: CodeNode, template: 'KernelTemplate', method: ODEMethod, ivp: IVP):
        """Generate Yasksite code tree.

        Parameters:
        -----------
        node: CodeNode
            Root node of the code tree.
        template: KernelTemplate
            Used kernel template.
        method: ODEMethod
            Used ODE method.
        ivp: IVP
            Used IVP.

        Returns:
        --------
        -
        """
        if node.type == CodeNodeType.LOOP:
            node: LoopNode
            # Evaluate loop boundary expression.
            constants: Optional[List[Tuple[str, Union[str, float, int]]]] = list()
            if method is not None:
                constants.extend([corrector_steps(method), stages(method)])
            node.boundary = eval_loop_boundary(node.boundary, constants)
        elif node.type == CodeNodeType.COMPUTATION:
            node: ComputationNode
            # Substitute computation.
            node.computation = template.computations[node.computation].computation
            # Check if computation includes RHS calls.
            if ivp is not None and '%RHS' in node.computation:
                assert 'RHS_input' in template.codegen
                assert 'RHS_butcher_nodes' in template.codegen
                node.computation = substitute_stencil_call(
                    node.computation, ivp.constants.as_tuple(), template.codegen['RHS_input'],
                    template.codegen['RHS_butcher_nodes'], ivp.code_stencil_path)
        # Traverse tree depth first.
        if node.child:
            YasksiteCodeGenerator.generate_yasksite_tree(node.child, template, method, ivp)
        if node.next:
            YasksiteCodeGenerator.generate_yasksite_tree(node.next, template, method, ivp)

    @staticmethod
    def generate_yasksite_code(node: CodeNode, code: str = '') -> str:
        """Write Yasksite code.

        Parameters:
        -----------
        node: CodeNode
            Root node of the code tree.
        code: str
            Current code string.

        Returns:
        --------
        str
            Generated Yasksite code.
        """
        if node.type != CodeNodeType.ROOT and node.type != CodeNodeType.LOOP:
            code += node.to_yasksite_codeline()
        # Traverse tree depth first.
        if node.child:
            code += YasksiteCodeGenerator.generate_yasksite_code(node.child)
        if node.next:
            code += YasksiteCodeGenerator.generate_yasksite_code(node.next)
        return code
