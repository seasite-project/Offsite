"""@package codegen.generator.kerncraft_generator
Definition of class KerncraftCodeGenerator.
"""

from copy import deepcopy
from re import sub
from typing import Dict, List, Optional, Tuple, Union

import attr
from sortedcontainers import SortedDict

import offsite.config
from offsite.codegen.code_dsl.code_node import CodeNode, CodeNodeType, ComputationNode, LoopNode, RootNode
from offsite.codegen.code_dsl.code_tree import CodeTree
from offsite.codegen.codegen_util import eval_loop_boundary, write_closing_bracket
from offsite.config import Config
from offsite.descriptions.ode import IVP, ODEMethod, corrector_steps, ivp_grid_size, ivp_system_size, stages
from offsite.descriptions.parser import ComponentDesc, DatastructDict
from offsite.util.math_utils import eval_math_expr


@attr.s
class KerncraftCodeGenerator:
    """Representation of the code generator for kerncraft kernel code.

    Attributes:
    -----------
    loop_stack: List of Tuple(str, LoopNode)
        Stack of opened loops.
    unroll_stack: SortedDict
        Stack of loops to be unrolled.
    unroll_counter: Dict
        Index counter for each of the unrolled loops.
    pmodel_trees: list of CodeNode
        Pmodel trees associated with the current kernel object.
    switch_to_nxt_pmodel: bool
        If true switch to next pmodel tree.
    """
    loop_stack = attr.ib(type=List[Tuple[str, LoopNode]], default=list())
    unroll_stack = attr.ib(type=SortedDict, default=SortedDict())
    unroll_counter = attr.ib(type=Dict, default=dict())
    pmodel_trees = attr.ib(type=List[CodeNode], default=list())
    switch_to_nxt_pmodel = attr.ib(type=bool, default=False, init=False)

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
        self.unroll_counter = dict()
        while True:
            self.unroll_stack.clear()
            self.collect_loop_meta_data(tree)
            if len(self.unroll_stack.keys()) == 0:
                break
            loop = self.unroll_stack.values()[0]
            CodeTree.unroll_loop(loop)
            self.unroll_counter[loop.var] = loop.start

    def generate(self, kernel: 'Kernel', method: Optional[ODEMethod] = None, ivp: Optional[IVP] = None,
                 system_size: Optional[int] = None) -> Dict[str, str]:
        """Generate pmodel code for a particular Kernel object.

        Parameters:
        -----------
        kernel: Kernel
            Kernel for which code is generated.
        method: ODEMethod
            ODE method for which code is generated.
        ivp: IVP
            IVP possibly evaluated by the kernel for which code is generated.
        system_size: int
            Used ODE system size.

        Returns:
        --------
        dict of str
            Generated pmodel codes.
        """
        # Set some members to match current code generation.
        code_tree: CodeNode = deepcopy(kernel.code_tree).root
        self.pmodel_trees = [code_tree]

        self.unroll_stack.clear()
        self.unroll_counter = dict()
        self.loop_stack = list()

        # Generate pmodel code trees.
        self.generate_pmodel_trees(code_tree, kernel.template, method, ivp, system_size)
        # Collect loop data.
        for tree in self.pmodel_trees:
            self.collect_loop_meta_data(tree)
        # Optimize tree: unroll
        for tree in self.pmodel_trees:
            self.unroll_tree(tree)
        # Substitute butcher table coefficients.
        if method is not None:
            for tree in self.pmodel_trees:
                CodeTree.substitute_butcher_coefficients(
                    tree, method.coefficientsA, method.coefficientsB, method.coefficientsC, True)
        # Update pmodel kernel names.
        idx: int = 1
        for tree in self.pmodel_trees:
            tree.name = kernel.name + ('_{}'.format(idx) if len(self.pmodel_trees) > 1 else '')
            idx += 1
        # Generate code from tree and write to string.
        codes: Dict[str, str] = dict()
        for tree in self.pmodel_trees:
            if kernel.template.isIVPdependent:
                assert ivp is not None
                config: Config = offsite.config.offsiteConfig
                # Determine RHS input variable.
                rhs_input: str = kernel.template.codegen[
                    'RHS_input'] if 'RHS_input' in kernel.template.codegen else config.ode_solution_vector
                # Determine IVP constants.
                ivp_constants: Optional[List[Tuple[str, Union[str, float, int]]]] = ivp.constants.as_tuple()
                if system_size is not None:
                    ivp_constants.extend([ivp_grid_size(eval_math_expr(ivp.gridSize, [ivp_system_size(system_size)],
                                                                       cast_to=int)), ivp_system_size(system_size)])
                # Generate code for each RHS component.
                idx: int = 1
                for rhs in ivp.components.values():
                    code_name: str = tree.name + '_{}'.format(ivp.name) + (
                        '_{}'.format(idx) if len(ivp.components.keys()) > 1 else '')
                    # Deepcopy required to assure that unroll counter is reset for each IVP component.
                    codes[code_name] = self.generate_kerncraft_code(
                        tree, rhs, rhs_input, ivp_constants, deepcopy(self.unroll_counter))
                    idx += 1
            else:
                codes[tree.name] = self.generate_kerncraft_code(tree)
        # Prepend variable definitions to written code.
        var_defs: str = self.construct_variable_defs(kernel.template.datastructs)
        codes = {name: var_defs + code for name, code in codes.items()}
        return codes

    @staticmethod
    def construct_variable_defs(datastructs: DatastructDict) -> str:
        """Construct variable definitions of a particular kernel.

        Parameters:
        -----------
        datastruct: DatastructDict
            Used data structures.

        Returns:
        --------
        str
            Variable definitions string.
        """
        var_defs: str = ''
        for name, desc in datastructs.items():
            var_defs += '{} {}{};\n'.format(desc.datatype, name, desc.dimensions_to_string())
        return var_defs + '\n'

    def generate_pmodel_trees(self, node: CodeNode, template: 'KernelTemplate', method: Optional[ODEMethod] = None,
                              ivp: Optional[IVP] = None, system_size: Optional[int] = None):
        """Generate pmodel code trees.

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
        system_size: int
            Used ODE system size.

        Returns:
        --------
        -
        """
        if self.switch_to_nxt_pmodel:
            assert node.parent is None
            # Cut off preceding pmodel node.
            node.prev.next = None
            node.prev = None
            # Prepend current loop stack.
            root: CodeNode = node
            if self.loop_stack:
                root = node
                for loop_node in reversed(self.loop_stack):
                    root.parent = deepcopy(loop_node)
                    assert loop_node.parent is None
                    root.parent.set_relatives(None, None, None, root)
                    root = root.parent
            # Add new root node.
            new_root = RootNode(0, CodeNodeType.ROOT, None, None, None, root, '')
            root.parent = new_root
            #
            self.pmodel_trees.append(new_root)
            self.switch_to_nxt_pmodel = False
        if node.type == CodeNodeType.LOOP:
            # Evaluate loop boundary expression.
            constants: Optional[List[Tuple[str, Union[str, float, int]]]] = list()
            if method is not None:
                constants.extend([corrector_steps(method), stages(method)])
            if ivp is not None and system_size is not None:
                constants.extend([ivp_grid_size(ivp.gridSize), ivp_system_size(system_size)])
            node.boundary = eval_loop_boundary(node.boundary, constants)
            #
            self.loop_stack.append(deepcopy(node))
        elif node.type == CodeNodeType.COMPUTATION:
            node: ComputationNode
            # Substitute computation.
            node.computation = template.computations[node.computation].computation
            # Replace constants.
            constants: Optional[List[Tuple[str, Union[str, float, int]]]] = list()
            if method is not None:
                constants.extend([corrector_steps(method), stages(method)])
            for constant in constants:
                node.computation = sub(r'\b{}\b'.format(constant[0]), str(constant[1]), node.computation)
            # Check if computation includes RHS calls.
            if ivp is not None and '%RHS' in node.computation:
                node.isIVPdependent = True
        elif node.type in [CodeNodeType.COMMUNICATION_CLUST, CodeNodeType.COMMUNICATION_NODE]:
            pass
        elif node.type == CodeNodeType.PMODEL:
            assert node.parent is None
            # Cut off tree before PMODEL node.
            node.prev.next = None
            #
            self.switch_to_nxt_pmodel = True
        # Traverse tree depth first.
        if node.child:
            self.generate_pmodel_trees(node.child, template, method, ivp, system_size)
            self.loop_stack = self.loop_stack[:-1]
        if node.next:
            self.generate_pmodel_trees(node.next, template, method, ivp, system_size)

    @staticmethod
    def generate_kerncraft_code(node: CodeNode, rhs: ComponentDesc = None, rhs_input: str = None,
                                constants: Optional[List[Tuple[str, Union[str, float, int]]]] = None,
                                unroll_counter: Optional[Dict[str, str]] = None, code: str = ''):
        """Write kerncraft code.

        Parameters:
        -----------
        node: CodeNode
            Root node of the code tree.
        system_size: int
            Used ODE system size.
        rhs: ComponentDesc
            IVP component to replace RHS calls.
        rhs_input: str
            Input vector used in RHS calls.
        constants: list of Tuple (str, str)
            Constants available.
        unroll_counter: dict (key=str, val=str)
            Index counter for each of the unrolled loops.
        code: str
            Current code string.

        Returns:
        --------
        str
            Generated kerncraft code.
        """
        config: Config = offsite.config.offsiteConfig
        if node.type == CodeNodeType.COMPUTATION and node.isIVPdependent:
            node: ComputationNode
            assert rhs is not None
            assert rhs_input is not None
            assert constants is not None
            assert unroll_counter is not None
            # Update RHS input to match unrolled loop indices.
            cur_rhs_input = rhs_input
            for var, counter in unroll_counter.items():
                # Update loop's unroll counter.
                if '[{}]'.format(var) in rhs_input:
                    cur_rhs_input = rhs_input.replace('[{}]'.format(var), '[{}]'.format(counter))
                    unroll_counter[var] = str(eval_math_expr(counter + '+ 1', cast_to=int))
            # Generate code.
            code += node.to_kerncraft_codeline(rhs.code, cur_rhs_input, constants)
        elif node.type == CodeNodeType.LOOP and node.var == config.var_idx and rhs is not None:
            node: LoopNode
            # Adjust node start and end to fit the IVP component considered.
            node.start = rhs.first
            node.boundary = '{} + {}'.format(node.start, rhs.equations)
            #
            code += node.to_kerncraft_codeline()
        elif node.type != CodeNodeType.ROOT:
            code += node.to_kerncraft_codeline()
        # Traverse tree depth first.
        if node.child:
            code += KerncraftCodeGenerator.generate_kerncraft_code(
                node.child, rhs, rhs_input, constants, unroll_counter)
        if node.type == CodeNodeType.LOOP:
            code += write_closing_bracket(node.indent)
        if node.next:
            code += KerncraftCodeGenerator.generate_kerncraft_code(node.next, rhs, rhs_input, constants, unroll_counter)
        return code
