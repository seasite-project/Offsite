"""@package codegen.code_dsl.code_tree
Definitions of classes CodeTree, CodeTreeGenerator.
"""

from copy import deepcopy
from re import sub
from typing import Dict, List, Optional, Tuple, Union

import attr
from lark import Visitor, tree as _tree, lexer as _lexer

from offsite.codegen.code_dsl.code_node import CodeNodeType, CodeNode, CommunicationClustLvlNode, \
    CommunicationNodeLvlNode, ComputationNode, KernelNode, LoopNode, PModelNode, RootNode, SwapNode
from offsite.codegen.codegen_util import replace_var_with_factor, replace_incr_with_assign_op
from offsite.codegen.computation_dsl.computation_tree import ComputationTreeVisitor
from offsite.descriptions.parser import ComputationDict
from offsite.util.math_utils import eval_math_expr


@attr.s
class CodeTree:
    root = attr.ib(type=CodeNode, default=RootNode(0, CodeNodeType.ROOT, None, None, None, None, ''))

    @staticmethod
    def substitute_butcher_coefficients(node: CodeNode, coeffs_a: List[List[str]], coeffs_b: List[str],
                                        coeffs_c: List[str], use_dummy_values=False):
        """
        Recursively substitute all occurrences of butcher table entries in the passed code_dsl tree with the coefficient
        values provided by the given ODE method.

        Parameters:
        -----------
        node: CodeNode
            Root node of the code_dsl tree.
        coeffs_a: list of list of str
            Coefficient matrix A of the ODE method.
        coeffs_b: list of str
            Coefficient vector b of the ODE method.
        coeffs_c: list of str
            Coefficient vector c of the ODE method.

        Returns:
        --------
        -
        """
        if node.type == CodeNodeType.COMPUTATION:
            node: ComputationNode
            # Replace A coefficients.
            for idx_row, A_row in enumerate(coeffs_a):
                for idx_col, elem in enumerate(A_row):
                    evaluated: str = eval_math_expr(elem, cast_to=str)
                    node.computation = node.computation.replace('A[{}][{}]'.format(idx_row, idx_col), evaluated)
            # Replace b coefficients.
            for idx, elem in enumerate(coeffs_b):
                evaluated: str = eval_math_expr(elem, cast_to=str)
                node.computation = node.computation.replace('b[{}]'.format(idx), evaluated)
            # Replace c coefficients.
            for idx, elem in enumerate(coeffs_c):
                evaluated: str = eval_math_expr(elem, cast_to=str)
                node.computation = node.computation.replace('c[{}]'.format(idx), evaluated)
            # Replace not already substituted coefficients with dummy values. Helps enabling running some kernels
            # in kerncraft LC mode.
            if use_dummy_values:
                node.computation = sub(r'A\[[\w|\s|\d][\w|\s|\d|+|-|*/]?\]\[[\w|\s|\d][\w|\s|\d|+-|*/]?\]',
                                       '0.123', node.computation)
                node.computation = sub(r'b\[[\w|\s|\d][\w|\s|\d|+|-|*|/]?\]', '0.456', node.computation)
                node.computation = sub(r'c\[[\w|\s|\d][\w|\s|\d|+|-|*|/]?\]', '0.789', node.computation)
        # Traverse tree depth first.
        if node.child:
            CodeTree.substitute_butcher_coefficients(
                node.child, coeffs_a, coeffs_b, coeffs_c, use_dummy_values)
        if node.next:
            CodeTree.substitute_butcher_coefficients(
                node.next, coeffs_a, coeffs_b, coeffs_c, use_dummy_values)

    @staticmethod
    def unroll_loop(loop: LoopNode):
        """Unroll the passed loop code_dsl tree.

        Parameters:
        -----------
        loop: LoopNode
            Root node of the code_dsl tree.
        skip_assign: bool
            If true do not apply assign optimization transformations.

        Returns:
        --------
        -
        """
        assert (loop.parent is not None) ^ (loop.prev is not None)
        assert loop.child is not None
        if loop.parent:
            parent: CodeNode = loop.parent
            parent.child = None
            # Cut loop body from original loop.
            loop_bdy: CodeNode = loop.child
            loop_bdy.parent = None
            # Update indention.
            CodeTree.update_indent(loop_bdy, -1)
            # Copy loop body and concatenate all unrolled loop bodies.
            parent.child = deepcopy(loop_bdy)
            parent.child.parent = parent
            # ..... unroll iteration 'start'.
            cur_loop_bdy: CodeNode = parent.child
            # Replace loop variable with current loop iteration index.
            cur: CodeNode = cur_loop_bdy
            while cur is not None:
                cur.computation = replace_var_with_factor(loop.var, cur.computation, loop.start)
                cur.computation = replace_incr_with_assign_op(cur.computation) if (
                        loop.optimize_assign is not None and loop.optimize_assign == int(
                    loop.start)) else cur.computation
                cur = cur.next
            # ..... unroll all remaining iterations.
            for idx in range(int(loop.start) + 1, int(loop.boundary)):
                loop_bdy_no_idx: CodeNode = deepcopy(loop_bdy)
                # Replace loop variable with current loop iteration index.
                cur: CodeNode = loop_bdy_no_idx
                while cur is not None:
                    cur.computation = replace_var_with_factor(loop.var, cur.computation, str(idx))
                    cur.computation = replace_incr_with_assign_op(cur.computation) if (
                            loop.optimize_assign is not None and loop.optimize_assign == idx) else cur.computation
                    cur = cur.next
                # Attach to end of loop body block.
                while cur_loop_bdy.next is not None:
                    cur_loop_bdy = cur_loop_bdy.next
                cur_loop_bdy.next = loop_bdy_no_idx
                loop_bdy_no_idx.prev = cur_loop_bdy
                cur_loop_bdy = loop_bdy_no_idx
            # Attach sibling nodes of old loop node (now unrolled).
            while cur_loop_bdy.next is not None:
                cur_loop_bdy = cur_loop_bdy.next
            cur_loop_bdy.next = loop.next
            if loop.next:
                loop.next.prev = cur_loop_bdy
        elif loop.prev:
            sibling: CodeNode = loop.prev
            sibling.next = None
            # Cut loop body from original loop.
            loop_bdy: CodeNode = loop.child
            loop_bdy.parent = None
            # Update indention.
            CodeTree.update_indent(loop_bdy, -1)
            # Copy loop body and concatenate all unrolled loop bodies.
            sibling.next = deepcopy(loop_bdy)
            sibling.next.prev = sibling
            # ..... unroll iteration 'start'.
            cur_loop_bdy: CodeNode = sibling.next
            # Replace loop variable with current loop iteration index.
            cur: CodeNode = cur_loop_bdy
            while cur is not None:
                cur.computation = replace_var_with_factor(loop.var, cur.computation, loop.start)
                cur.computation = replace_incr_with_assign_op(cur.computation) if (
                        loop.optimize_assign is not None and loop.optimize_assign == int(
                    loop.start)) else cur.computation
                cur = cur.next
            # ..... unroll iteration 'start'.
            for idx in range(int(loop.start) + 1, int(loop.boundary)):
                loop_bdy_no_idx: CodeNode = deepcopy(loop_bdy)
                # Replace loop variable with current loop iteration index.
                cur: CodeNode = loop_bdy_no_idx
                while cur is not None:
                    cur.computation = replace_var_with_factor(loop.var, cur.computation, str(idx))
                    cur.computation = replace_incr_with_assign_op(cur.computation) if (
                            loop.optimize_assign is not None and loop.optimize_assign == idx) else cur.computation
                    cur = cur.next
                while cur_loop_bdy.next is not None:
                    cur_loop_bdy = cur_loop_bdy.next
                cur_loop_bdy.next = loop_bdy_no_idx
                loop_bdy_no_idx.prev = cur_loop_bdy
                cur_loop_bdy = loop_bdy_no_idx
            # Attach sibling nodes of old loop node (now unrolled).
            while cur_loop_bdy.next is not None:
                cur_loop_bdy = cur_loop_bdy.next
            cur_loop_bdy.next = loop.next
            if loop.next:
                loop.next.prev = cur_loop_bdy

    @staticmethod
    def split_loop(loop: LoopNode):
        """Split the passed loop code_dsl tree.

        Parameters:
        -----------
        loop: LoopNode
            Root node of the code_dsl tree.

        Returns:
        --------
        -
        """
        assert (loop.parent is not None) ^ (loop.prev is not None)
        assert loop.child is not None
        assert loop.split_at == 0
        # At the moment only supports splitting the first iteration.
        if loop.parent:
            parent: CodeNode = loop.parent
            parent.child = None
            loop.parent = None
            # Copy loop body of the original loop.
            loop_bdy: CodeNode = deepcopy(loop.child)
            loop_bdy.parent = None
            # Update indention.
            CodeTree.update_indent(loop_bdy, 1)
            # Update loop variable with splitted loop iteration index.
            CodeTree.replace_loop_var_with_iteration_idx(loop_bdy, loop)
            # Update loop start index.
            loop.start = loop.split_at + 1
            # Attach loop body as first child to parent loop.
            parent.child = loop_bdy
            loop_bdy.parent = parent
            # Prepend split loop body to original loop.
            cur_loop_bdy = loop_bdy
            while cur_loop_bdy.next is not None:
                cur_loop_bdy = cur_loop_bdy.next
            cur_loop_bdy.next = loop
            loop.prev = cur_loop_bdy
        elif loop.prev:
            sibling: CodeNode = loop.prev
            sibling.next = None
            loop.prev = None
            # Copy loop body of the original loop.
            loop_bdy: CodeNode = deepcopy(loop.child)
            loop_bdy.parent = None
            # Update indention.
            CodeTree.update_indent(loop_bdy, 1)
            # Update loop variable with splitted loop iteration index.
            CodeTree.replace_loop_var_with_iteration_idx(loop_bdy, loop)
            # Update loop start index.
            loop.start = loop.split_at + 1
            # Prepend loop body to original loop.
            sibling.next = loop_bdy
            loop_bdy.prev = sibling
            cur_loop_bdy = loop_bdy
            while cur_loop_bdy.next is not None:
                cur_loop_bdy = cur_loop_bdy.next
            cur_loop_bdy.next = loop
            loop.prev = cur_loop_bdy

    @staticmethod
    def visit(node: CodeNode):
        """Traverse the passed loop code_dsl tree and print information all visited nodes.

        Parameters:
        -----------
        node: CodeNode
            Root node of the code_dsl tree.

        Returns:
        --------
        -
        """
        if node.child:
            CodeTree.visit(node.child)
        if node.next:
            CodeTree.visit(node.next)

    @staticmethod
    def replace_loop_var_with_iteration_idx(node: CodeNode, loop: LoopNode):
        """
        Replace within the passed code_dsl tree, all occurrences of a given loop variable with a given loop iteration
        index.

        Parameters:
        -----------
        node: CodeNode
            Root node of the code_dsl tree.
        loop: LoopNode
            Used loop node that determines the loop variable name and loop iteration index.

        Returns:
        --------
        -
        """
        if node is None:
            return
        if node.type == CodeNodeType.COMPUTATION:
            node: ComputationNode
            node.computation = replace_var_with_factor(loop.var, node.computation, str(loop.split_at))
            node.computation = replace_incr_with_assign_op(node.computation) if loop.split_at == 0 else node.computation
        if node.child:
            CodeTree.replace_loop_var_with_iteration_idx(node.child, loop)
        if node.next:
            CodeTree.replace_loop_var_with_iteration_idx(node.next, loop)

    @staticmethod
    def iteration_count(node: CodeNode, iter_count: str = '', ignore_unroll: bool = False) -> str:
        """Count the total number of loop iteration of the current node's subtree.

        Parameters:
        -----------
        node: CodeNode
            Start node.
        iter_count: str
            Current loop iteration count.

        Returns:
        --------
        str
            Total number of loop iterations of the start node's subtree.
        """
        if node is None:
            return iter_count
        if node.type == CodeNodeType.LOOP and ('unroll' not in node.optional_args or ignore_unroll):
            node: LoopNode
            if iter_count:
                iter_count = '{} * ({} - {})'.format(iter_count, node.boundary, node.start)
            else:
                iter_count = '({} - {})'.format(node.boundary, node.start)
        if node.child:
            iter_count = CodeTree.iteration_count(node.child, iter_count, ignore_unroll)
        if node.next and node.type is not CodeNodeType.PMODEL:
            iter_count = CodeTree.iteration_count(node.next, iter_count, ignore_unroll)
        return iter_count

    @staticmethod
    def iteration_count_up_to_root(node: CodeNode, iter_count: str = '', ignore_unroll: bool = False) -> str:
        """Count the total number of loop iteration from the current node up to the root node.

        Parameters:
        -----------
        node: CodeNode
            Start node.
        iter_count: str
            Current loop iteration count.

        Returns:
        --------
        str
            Total number of loop iterations from the start node to the root node.
        """
        if node is None:
            return iter_count
        if node.type == CodeNodeType.LOOP and ('unroll' not in node.optional_args or ignore_unroll):
            node: LoopNode
            if iter_count:
                iter_count = '{} * ({} - {})'.format(iter_count, node.boundary, node.start)
            else:
                iter_count = '({} - {})'.format(node.boundary, node.start)
        if node.parent:
            iter_count = CodeTree.iteration_count(node.parent, iter_count, ignore_unroll)
        return iter_count

    @staticmethod
    def find_pmodel_node(node: CodeNode, requested_idx: int, nxt_visited_idx=0,
                         found_node: Optional[PModelNode] = None) -> Tuple[Optional[PModelNode], int]:
        """Traverse node's subtree and return the 'requested_idx' occurrence of a PModelNode object.

        Parameters:
        -----------
        node: CodeNode
            Root node of the code_dsl tree.
        requested_idx: int
            xxx
        nxt_visited: int
            xyy
        found_node: CodeNode
            Used loop node that determines the loop variable name and loop iteration index.

        Returns:
        --------
        PModelNode
            Found node.
            Returns None if no proper CodeNode object was found.
        """
        if node.type == CodeNodeType.ROOT:
            if requested_idx == 0:
                return node, nxt_visited_idx + 1
            nxt_visited_idx += 1
        elif node.type == CodeNodeType.PMODEL:
            if nxt_visited_idx == requested_idx:
                return node, nxt_visited_idx + 1
            nxt_visited_idx += 1
        if node.child:
            found_node, nxt_visited_idx = CodeTree.find_pmodel_node(
                node.child, requested_idx, nxt_visited_idx, found_node)
        if node.next:
            found_node, nxt_visited_idx = CodeTree.find_pmodel_node(
                node.next, requested_idx, nxt_visited_idx, found_node)
        return found_node, nxt_visited_idx

    @staticmethod
    def count_clust_lvl_communication(node: CodeNode, comm_count: Dict[str, str] = dict()) -> Dict[str, str]:
        """Count the total number of cluster-level communication operation executions in the current node's subtree.

        Parameters:
        -----------
        node: CodeNode
            Start node.
        comm_count: Dict (key=str, val=str)
            Current communication operation execution count.

        Returns:
        --------
        Dict (key=str, val=str)
            Total communication operation executions count of the start node's subtree.
        """
        if node is None:
            return comm_count
        if node.type == CodeNodeType.COMMUNICATION_CLUST:
            node: CommunicationClustLvlNode
            # Count number of executions of this operation.
            execs = CodeTree.iteration_count_up_to_root(node.my_parent(), ignore_unroll=True)
            execs = execs if execs else '1'
            #
            op: str = node.operation
            op_count: str = comm_count[op] if op in comm_count else ''
            comm_count[op] = eval_math_expr('{}+{}'.format(op_count, execs), cast_to=str)
        if node.child:
            comm_count = CodeTree.count_clust_lvl_communication(node.child, comm_count)
        if node.next:
            comm_count = CodeTree.count_clust_lvl_communication(node.next, comm_count)
        return comm_count

    @staticmethod
    def count_node_lvl_communication(node: CodeNode, comm_count: Dict[str, str] = dict()) -> Dict[str, str]:
        """Count the total number of node-level communication operation executions in the current node's subtree.

        Parameters:
        -----------
        node: CodeNode
            Start node.
        comm_count: Dict (key=str, val=str)
            Current communication operation execution count.

        Returns:
        --------
        Dict (key=str, val=str)
            Total communication operation executions count of the start node's subtree.
        """
        if node is None:
            return comm_count
        if node.type == CodeNodeType.COMMUNICATION_NODE:
            node: CommunicationNodeLvlNode
            # Count number of executions of this operation.
            execs = CodeTree.iteration_count_up_to_root(node.my_parent(), ignore_unroll=True)
            execs = execs if execs else '1'
            #
            op: str = node.operation
            op_count: str = comm_count[op] if op in comm_count else ''
            comm_count[op] = eval_math_expr('{}+{}'.format(op_count, execs), cast_to=str)
        if node.child:
            comm_count = CodeTree.count_node_lvl_communication(node.child, comm_count)
        if node.next:
            comm_count = CodeTree.count_node_lvl_communication(node.next, comm_count)
        return comm_count

    @staticmethod
    def count_kernel(node: CodeNode, kernel_count: Dict[str, str] = dict()) -> Dict[str, str]:
        """Count the total number of kernel executions in the current node's subtree.

        Parameters:
        -----------
        node: CodeNode
            Start node.
        kernel_count: Dict (key=str, val=str)
            Current kernel execution count.

        Returns:
        --------
        Dict (key=str, val=str)
            Total kernel executions count of the start node's subtree.
        """
        if node is None:
            return kernel_count
        if node.type == CodeNodeType.KERNEL:
            node: KernelNode
            # Count number of executions of this kernel.
            execs = CodeTree.iteration_count_up_to_root(node.my_parent(), ignore_unroll=True)
            execs = execs if execs else '1'
            #
            name: str = node.template_name
            count: str = kernel_count[name] if name in kernel_count else ''
            kernel_count[name] = eval_math_expr('{}+{}'.format(count, execs), cast_to=str)
        if node.child:
            kernel_count = CodeTree.count_kernel(node.child, kernel_count)
        if node.next:
            kernel_count = CodeTree.count_kernel(node.next, kernel_count)
        return kernel_count

    @staticmethod
    def update_indent(node: CodeNode, change_by_val: int):
        """Update indention of the passed code_dsl tree by a given value.

        Parameters:
        -----------
        loop: CodeNode
            Root node of the code_dsl tree.
        change_by_val: int
            Increment the indention level by 'change_by_val'.

        Returns:
        --------
        -
        """
        node.indent = max(node.indent + change_by_val, 0)
        if node.child:
            CodeTree.update_indent(node.child, change_by_val)
        if node.next:
            CodeTree.update_indent(node.next, change_by_val)


class CodeTreeGenerator(Visitor):
    def __init__(self):
        super().__init__()

        self.code_tree = CodeTree()

        self.last_visited_node: CodeNode = self.code_tree.root
        self.switch_to_nxt_lvl: bool = False
        self.current_loop_lvl: Optional[LoopNode] = None
        self.current_indent_lvl: int = 0
        self.current_pragmas: List[str] = list()

        self.frame_lvl: int = 2

        self.available_computations: List[ComputationDict] = list()

    def get_relatives(self) -> Tuple[Optional[CodeNode], Optional[CodeNode], Optional[CodeNode], Optional[CodeNode]]:
        if self.switch_to_nxt_lvl:
            prev = None
            parent = self.current_loop_lvl
        else:
            prev = self.last_visited_node
            parent = None
        nxt = None
        child = None
        return prev, nxt, parent, child

    def update_loop_lvl(self, tree: _tree.Tree):
        lvl_gap = tree.lvl - self.current_indent_lvl
        while lvl_gap != 1 and self.current_loop_lvl is not None:
            # Update global information: last visited node.
            self.last_visited_node = self.current_loop_lvl
            # Move up a level in the loop hierarchy.
            self.current_loop_lvl = self.move_up_loop_lvl(self.current_loop_lvl)
            # Update level gap.
            lvl_gap += 1
            # Can't move above root level.
            assert self.current_indent_lvl > 0

    def move_up_loop_lvl(self, node: CodeNode) -> Optional[LoopNode]:
        self.current_indent_lvl -= 1
        #
        parent: LoopNode = node.parent
        if parent:
            return parent
        while True:
            node = node.prev
            if not node:
                return None
            if node.parent:
                parent = node.parent
                return parent

    @staticmethod
    def read_loop(tree: _tree.Tree):
        prev_is_pragma = False
        base_idx = 0
        # Read pragmas.
        for entry in tree.children:
            if isinstance(entry, _tree.Tree) and entry.data == 'pragma':
                base_idx += 1
                prev_is_pragma = True
            else:
                if prev_is_pragma and isinstance(entry, _lexer.Token) and entry.value == '\n':
                    base_idx += 1
                prev_is_pragma = False
        # Read variable name.
        var = str(tree.children[base_idx])
        # Read loop boundary expression.
        bound: Optional[Union[_tree.Tree, str]] = tree.children[base_idx + 1]
        if isinstance(bound, _tree.Tree):
            bound = ''.join(bound.children)
        else:
            bound = str(tree.children[base_idx + 1])
        # Read further optional arguments.
        opt_args = [str(entry) for entry in tree.children[(base_idx + 2):]]
        return var, bound, opt_args

    def loop(self, tree: _tree.Tree):
        self.update_loop_lvl(tree)
        # Read tree information.
        loop_var, loop_boundary, optional_args = self.read_loop(tree)
        # Create loop node.
        prev, nxt, parent, child = self.get_relatives()
        node = LoopNode(tree.lvl - self.frame_lvl, CodeNodeType.LOOP, prev, nxt, parent, child, loop_var, str(0),
                        loop_boundary, optional_args, pragmas=self.current_pragmas)
        self.current_pragmas = list()
        # Update information.
        if not self.switch_to_nxt_lvl:
            # Update local information: previous sibling node
            self.last_visited_node.next = node
        else:
            # Update global information: double nested loop -> add current node as innermost loop level
            if self.current_loop_lvl:
                self.current_loop_lvl.child = node
        # Update global information: last visited node
        self.last_visited_node = node
        # Move down a level in the loop hierarchy.
        self.switch_to_nxt_lvl = True
        self.current_loop_lvl = node
        self.current_indent_lvl = tree.lvl

    def pragma(self, tree: _tree.Tree):
        # Read tree information.
        pragma_id = str(tree.children[0])
        # Create pragma string.
        self.current_pragmas.append('#pragma ' + pragma_id)

    def comp(self, tree: _tree.Tree):
        self.update_loop_lvl(tree)
        # Read tree information.
        comp_id = str(tree.children[0])
        # Get all references.
        references = ComputationTreeVisitor().get_references(self.available_computations[comp_id].computation)
        # Create computation node.
        prev, nxt, parent, child = self.get_relatives()
        node = ComputationNode(tree.lvl - self.frame_lvl, CodeNodeType.COMPUTATION, prev, nxt, parent, child, comp_id,
                               references)
        # Update information.
        if self.switch_to_nxt_lvl:
            self.switch_to_nxt_lvl = False
            # Update global information: double nested loop -> add current node as innermost loop level
            if self.current_loop_lvl:
                self.current_loop_lvl.child = node
        else:
            # Update local information: previous sibling node
            self.last_visited_node.next = node
        # Update global information: last visited node
        self.last_visited_node = node

    def pmodel(self, tree: _tree.Tree):
        self.update_loop_lvl(tree)
        # Create pmodel node.
        prev, nxt, parent, child = self.get_relatives()
        node = PModelNode(self.current_indent_lvl, CodeNodeType.PMODEL, prev, nxt, parent, child)
        # Update information.
        if self.switch_to_nxt_lvl:
            self.switch_to_nxt_lvl = False
            # Update global information: double nested loop -> add current node as innermost loop level
            if self.current_loop_lvl:
                self.current_loop_lvl.child = node
        else:
            # Update local information: previous sibling node
            self.last_visited_node.next = node
        # Update global information: last visited node
        self.last_visited_node = node

    def comm_clust(self, tree: _tree.Tree):
        self.update_loop_lvl(tree)
        # Read tree information.
        comm_id = str(tree.children[0])
        # Read optional input vector name argument.
        input_vec: Optional[Union[_tree.Tree, str]] = tree.children[1] if len(tree.children) > 1 else None
        if input_vec is not None:
            if isinstance(input_vec, _tree.Tree):
                input_vec = '{}[{}]'.format(input_vec.children[0], input_vec.children[1])
            else:
                input_vec = str(input_vec)
        # Create communication node.
        prev, nxt, parent, child = self.get_relatives()
        node = CommunicationClustLvlNode(
            tree.lvl - self.frame_lvl, CodeNodeType.COMMUNICATION_CLUST, prev, nxt, parent, child, comm_id, input_vec)
        # Update information.
        if self.switch_to_nxt_lvl:
            self.switch_to_nxt_lvl = False
            # Update global information: double nested loop -> add current node as innermost loop level
            if self.current_loop_lvl:
                self.current_loop_lvl.child = node
        else:
            # Update local information: previous sibling node
            self.last_visited_node.next = node
        # Update global information: last visited node
        self.last_visited_node = node

    def comm_node(self, tree: _tree.Tree):
        self.update_loop_lvl(tree)
        # Read tree information.
        comm_id = str(tree.children[0])
        # Read optional input vector name argument.
        input_vec: Optional[Union[_tree.Tree, str]] = tree.children[1] if len(tree.children) > 1 else None
        if input_vec is not None:
            if isinstance(input_vec, _tree.Tree):
                input_vec = '{}[{}]'.format(input_vec.children[0], input_vec.children[1])
            else:
                input_vec = str(input_vec)
        # Create communication node.
        prev, nxt, parent, child = self.get_relatives()
        node = CommunicationNodeLvlNode(
            tree.lvl - self.frame_lvl, CodeNodeType.COMMUNICATION_NODE, prev, nxt, parent, child, comm_id, input_vec)
        # Update information.
        if self.switch_to_nxt_lvl:
            self.switch_to_nxt_lvl = False
            # Update global information: double nested loop -> add current node as innermost loop level
            if self.current_loop_lvl:
                self.current_loop_lvl.child = node
        else:
            # Update local information: previous sibling node
            self.last_visited_node.next = node
        # Update global information: last visited node
        self.last_visited_node = node

    def kernel(self, tree: _tree.Tree):
        self.update_loop_lvl(tree)
        # Read tree information.
        kernel_id = str(tree.children[0])
        input_vec: Optional[Union[_tree.Tree, str]] = tree.children[1] if len(tree.children) > 1 else None
        if input_vec is not None:
            if isinstance(input_vec, _tree.Tree):
                input_vec = '{}[{}]'.format(input_vec.children[0], input_vec.children[1])
            else:
                input_vec = str(input_vec)
        # Create kernel node.
        prev, nxt, parent, child = self.get_relatives()
        node = KernelNode(
            tree.lvl - self.frame_lvl, CodeNodeType.KERNEL, prev, nxt, parent, child, kernel_id, input_vec)
        # Update information.
        if self.switch_to_nxt_lvl:
            self.switch_to_nxt_lvl = False
            # Update global information: double nested loop -> add current node as innermost loop level
            if self.current_loop_lvl:
                self.current_loop_lvl.child = node
        else:
            # Update local information: previous sibling node
            self.last_visited_node.next = node
        # Update global information: last visited node
        self.last_visited_node = node

    def swap(self, tree: _tree.Tree):
        self.update_loop_lvl(tree)
        # Read tree information.
        arg1 = str(tree.children[0])
        arg2 = str(tree.children[1])
        datatype = str(tree.children[2])
        # Create swap node.
        prev, nxt, parent, child = self.get_relatives()
        node = SwapNode(tree.lvl - self.frame_lvl, CodeNodeType.SWAP, prev, nxt, parent, child, arg1, arg2, datatype)
        # Update information.
        if self.switch_to_nxt_lvl:
            self.switch_to_nxt_lvl = False
            # Update global information: double nested loop -> add current node as innermost loop level
            if self.current_loop_lvl:
                self.current_loop_lvl.child = node
        else:
            # Update local information: previous sibling node
            self.last_visited_node.next = node
        # Update global information: last visited node
        self.last_visited_node = node

    def generate(self, lark_tree: _tree.Tree, computations: ComputationDict) -> CodeTree:
        self.available_computations = computations
        #
        lark_tree.lvl = 0
        LevelCounter().visit_topdown(lark_tree)
        #
        self.visit_topdown(lark_tree)
        return self.code_tree


class LevelCounter(Visitor):
    def __default__(self, tree):
        for subtree in tree.children:
            if isinstance(subtree, _tree.Tree):
                assert not hasattr(subtree, 'parent')
                assert not hasattr(subtree, 'lvl')
                try:
                    subtree.lvl = tree.lvl + 1
                except AttributeError:
                    subtree.lvl = 1
