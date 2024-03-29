"""@package codegen.generator.kernel_generator
Definition of class KernelCodeGenerator.

@author: Johannes Seiferth
"""

from copy import deepcopy
from typing import Dict, List, Optional, Tuple, Union

import attr
from sortedcontainers import SortedDict

import offsite.config
from offsite.codegen.code_dsl.code_node import CodeNode, CodeNodeType, ComputationNode, LoopNode, RootNode
from offsite.codegen.code_dsl.code_tree import CodeTree
from offsite.codegen.codegen_util import eval_loop_boundary, substitute_rhs_func_call, write_closing_bracket, \
    write_tiling_loop
from offsite.config import Config
from offsite.descriptions.impl.kernel_template import Kernel, KernelTemplate
from offsite.descriptions.ode import ODEMethod, corrector_steps, stages


@attr.s
class KernelCodeGenerator:
    """Representation of the code generator for kernel code.

    Attributes:
    -----------
    unroll_stack: SortedDict
        Stack of loops to be unrolled.
    split_stack: dict (key=str, val=list of CodeNode)
        Stack of loops to be split.
    """
    unroll_stack = attr.ib(type=SortedDict, default=SortedDict())
    split_stack = attr.ib(type=dict, default=dict())

    def collect_loop_meta_data(self, node: CodeNode, loop_splits: Optional[Dict[str, int]] = None):
        """Traverse tree depth first and collect all unroll and split tasks.

        Parameters:
        -----------
        node: CodeNode
            Root node of the code tree.
        loop_splits: dict (key=str, val=str)
            Requested loop splits (key=loop var; val=split at loop index)

        Returns:
        --------
        -
        """
        if node.type == CodeNodeType.LOOP:
            node: LoopNode
            # Collect unroll tasks.
            if node.optimize_unroll is not None:
                if node.optimize_unroll in self.unroll_stack:
                    self.unroll_stack[node.optimize_unroll].append(node)
                else:
                    self.unroll_stack[node.optimize_unroll] = [node]
            # Collect split tasks.
            if loop_splits:
                if node.var in loop_splits:
                    node.split_at = loop_splits[node.var]
                    if node.var in self.split_stack:
                        self.split_stack[node.var].append(node)
                    else:
                        self.split_stack[node.var] = [node]
        # Traverse tree depth first.
        if node.child:
            self.collect_loop_meta_data(node.child, loop_splits)
        if node.next:
            self.collect_loop_meta_data(node.next, loop_splits)

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
            CodeTree.unroll_loop(self.unroll_stack.values()[0][0])

    def generate(self, kernel: Kernel, method: ODEMethod, input_vector: str, gen_tiled_code: bool) -> str:
        """Generate kernel code for a particular Kernel object.

        Parameters:
        -----------
        kernel: Kernel
            Kernel for which code is generated.
        method: ODEMethod
            ODE method for which code is generated.
        input_vector: str
            Name of the used input vector name.
        gen_tiled_code: bool
            Generated tiled code version.

        Returns:
        --------
        str
            Generated kernel code.
        """
        # Set some members to match current code generation.
        code_tree: CodeNode = deepcopy(kernel.code_tree).root

        self.unroll_stack.clear()
        self.split_stack = dict()

        # Generate kernel code tree.
        self.generate_kernel_tree(code_tree, kernel, method, input_vector)
        # Read loop split information from template.
        loop_splits: Optional[Dict[str, int]] = dict()
        if 'loop splits' in kernel.template.codegen:
            for ls in kernel.template.codegen['loop splits']:
                var, split_at = ls.split(' ')
                loop_splits[var] = int(split_at)
        # Collect loop data.
        self.collect_loop_meta_data(code_tree, loop_splits)
        # Optimize tree: split
        if self.split_stack:
            for loops in self.split_stack.values():
                for loop in loops:
                    CodeTree.split_loop(loop)
        # Collect loop data again to include possibly split loops.
        self.collect_loop_meta_data(code_tree, loop_splits)
        # Optimize tree: unroll
        self.unroll_tree(code_tree)
        # Substitute butcher table coefficients.
        if method is not None:
            CodeTree.substitute_butcher_coefficients(
                code_tree, method.coefficientsA, method.coefficientsB, method.coefficientsC)
        #
        rhs_func: Optional[str]
        if kernel.template.isIVPdependent:
            rhs_func = kernel.template.codegen['RHS'] if 'RHS' not in kernel.codegen else kernel.codegen['RHS']
        else:
            rhs_func = None
        # Generate code from tree and write to string.
        return self.generate_kernel_code(code_tree, gen_tiled_code, rhs_func)

    @staticmethod
    def generate_kernel_tree(node: CodeNode, kernel: Kernel, method: ODEMethod, input_vector: str):
        """Generate yasksite code tree.

        Parameters:
        -----------
        node: CodeNode
            Root node of the code tree.
        kernel: Kernel
            Used kernel.
        method: ODEMethod
            Used ODE method.
        input_vector: str
            Name of the used input vector.

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
            template: KernelTemplate = kernel.template
            # Substitute computation.
            node.computation = template.computations[node.computation].computation
            # Check if computation includes RHS calls.
            if '%RHS' in node.computation:
                assert ('RHS' in template.codegen or 'RHS' in kernel.codegen)
                assert 'RHS_butcher_nodes' in template.codegen
                rhs_func: str = template.codegen['RHS'] if 'RHS' not in kernel.codegen else kernel.codegen['RHS']
                butcher_node: str = template.codegen['RHS_butcher_nodes']
                node.computation = substitute_rhs_func_call(node.computation, rhs_func, input_vector, butcher_node)
                if rhs_func == 'eval_range':
                    node.parent.flag = True
        elif node.type in [CodeNodeType.COMMUNICATION_CLUST, CodeNodeType.COMMUNICATION_NODE]:
            pass
        elif node.type == CodeNodeType.ROOT:
            node.name = kernel.name
        # Traverse tree depth first.
        if node.child:
            KernelCodeGenerator.generate_kernel_tree(node.child, kernel, method, input_vector)
        if node.next:
            KernelCodeGenerator.generate_kernel_tree(node.next, kernel, method, input_vector)

    @staticmethod
    def generate_kernel_code(node: CodeNode, gen_tiled_code: bool, rhs_func: Optional[str] = None,
                             kernel_name: str = '', code: str = ''):
        """Write kernel code.

        Parameters:
        -----------
        node: CodeNode
            Root node of the code tree.
        rhs_func: str
            Name of the used RHS function.
        gen_tiled_code: bool
            Generated tiled code version.
        kernel_name: str
            Name of the generated kernel.
        code: str
            Current code string.

        Returns:
        --------
        str
            Generated kernel code.
        """
        config: Config = offsite.config.offsiteConfig
        loop_skipped = False
        if node.type == CodeNodeType.LOOP and node.var == config.var_idx:
            node: LoopNode
            if node.flag:
                # Loop is replaced by RHS eval_range call and can thus be skipped.
                loop_skipped = True
            else:
                if gen_tiled_code:
                    code += node.to_tiling_kernel_codeline(kernel_name)
                else:
                    code += node.to_kernel_codeline()
        elif node.type == CodeNodeType.ROOT and gen_tiled_code:
            node: RootNode
            kernel_name = node.name
            code += write_tiling_loop(kernel_name, node.indent)
        elif node.type != CodeNodeType.ROOT and node.type != CodeNodeType.PMODEL:
            code += node.to_kernel_codeline()
        # Traverse tree depth first.
        if node.child:
            code += KernelCodeGenerator.generate_kernel_code(node.child, gen_tiled_code, rhs_func, kernel_name)
        if node.type == CodeNodeType.LOOP:
            if not loop_skipped:
                code += write_closing_bracket(node.indent)
        if node.next:
            code += KernelCodeGenerator.generate_kernel_code(node.next, gen_tiled_code, rhs_func, kernel_name)
        if node.type == CodeNodeType.ROOT and gen_tiled_code:
            # Close previously opened tiling loops.
            code += write_closing_bracket(node.indent)
            code += write_closing_bracket(node.indent)
        return code
