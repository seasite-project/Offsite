"""@package kernel_generator
Definition of class KernelCodeGenerator.
"""

from copy import deepcopy
from typing import Dict

import attr
from sortedcontainers import SortedDict

import offsite.config
from offsite.codegen.code_tree import CodeTree, CodeNode, CodeNodeType
from offsite.codegen.codegen_util import eval_loop_boundary, substitute_rhs_func, write_closing_bracket
from offsite.descriptions.kernel_template import Kernel
from offsite.descriptions.ode_method import ODEMethod


@attr.s
class KernelCodeGenerator:
    """Representation of the code generator for kernel code.

    Attributes:
    -----------
    unroll_stack : SortedDict
        Stack of loops to be unrolled.
    split_stack : dict (key=str, val=list of CodeNode)
        Stack of loops to be split.
    """
    unroll_stack = attr.ib(type=SortedDict, default=SortedDict())
    split_stack = attr.ib(type=dict, default=dict())

    def collect_loop_meta_data(self, node: CodeNode, loop_splits: Dict[str, str] = None):
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

    def generate(self, kernel: 'Kernel', method: 'ODEMethod', input_vector: str) -> str:
        """Generate kernel code for a particular Kernel object.

        Parameters:
        -----------
        kernel: Kernel
            Kernel for which code is generated.
        method: ODEMethod
            ODE method for which code is generated.
        input_vector: str
            Name of the used input vector name.

        Returns:
        --------
        str
            Generated kernel code.
        """
        # Set some members to match current code generation.
        code_tree = deepcopy(kernel.code_tree).root

        self.unroll_stack.clear()
        self.split_stack = dict()

        # Generate kernel code tree.
        self.generate_kernel_tree(code_tree, kernel, method, input_vector)
        # Read loop split information from template.
        loop_splits = dict()
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
        # Collect loop data again to include possibly splitted loops.
        self.collect_loop_meta_data(code_tree, loop_splits)
        # Optimize tree: unroll
        self.unroll_tree(code_tree)
        # Substitute butcher table coefficients.
        CodeTree.substitute_butcher_coefficients(code_tree, method)
        #
        if kernel.template.isIVPdependent:
            rhs_func = kernel.template.codegen['RHS'] if 'RHS' not in kernel.codegen else kernel.codegen['RHS']
        else:
            rhs_func = None
        # Generate code from tree and write to string.
        code = self.generate_kernel_code(code_tree, rhs_func)
        return code

    @staticmethod
    def generate_kernel_tree(node: CodeNode, kernel: 'Kernel', method: 'ODEMethod', input_vector: str):
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
            # Evaluate loop boundary expression.
            node.boundary = eval_loop_boundary(method, node.boundary)
        elif node.type == CodeNodeType.COMPUTATION:
            template = kernel.template
            # Substitute computation.
            node.computation = template.computations[node.computation].computation
            # Check if computation includes RHS calls.
            if '%RHS' in node.computation:
                assert 'RHS' in template.codegen
                assert 'RHS_butcher_nodes' in template.codegen
                rhs_func = template.codegen['RHS'] if 'RHS' not in kernel.codegen else kernel.codegen['RHS']
                butcher_node = template.codegen['RHS_butcher_nodes']
                node.computation = substitute_rhs_func(node.computation, rhs_func, input_vector, butcher_node)
                if rhs_func == 'eval_range':
                    node.parent.flag = True
        elif node.type == CodeNodeType.COMMUNICATION:
            # TODO add support later
            pass
        # Traverse tree depth first.
        if node.child:
            KernelCodeGenerator.generate_kernel_tree(node.child, kernel, method, input_vector)
        if node.next:
            KernelCodeGenerator.generate_kernel_tree(node.next, kernel, method, input_vector)

    @staticmethod
    def generate_kernel_code(node: CodeNode, rhs_func: str = None, code: str = ''):
        """Write kernel code.

        Parameters:
        -----------
        node: CodeNode
            Root node of the code tree.
        rhs_func: str
            Name of the used RHS function.
        code: str
            Current code string.

        Returns:
        --------
        str
            Generated kernel code.
        """
        config = offsite.config.offsiteConfig
        loop_skipped = False
        if node.type == CodeNodeType.LOOP and node.var == config.var_idx:
            if node.flag:
                # Loop is replaced by RHS eval_range call and can thus be skipped.
                loop_skipped = True
            else:
                code += node.to_kernel_codeline()
        elif node.type != CodeNodeType.ROOT and node.type != CodeNodeType.PMODEL:
            code += node.to_kernel_codeline()
        # Traverse tree depth first.
        if node.child:
            code += KernelCodeGenerator.generate_kernel_code(node.child, rhs_func)
        if node.type == CodeNodeType.LOOP:
            if not loop_skipped:
                code += write_closing_bracket(node.indent)
        if node.next:
            code += KernelCodeGenerator.generate_kernel_code(node.next, rhs_func)
        return code
