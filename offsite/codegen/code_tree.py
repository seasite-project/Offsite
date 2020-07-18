"""@package code_tree
Definitions of classes CodeNodeType, CodeNode, RootNode, LoopNode, CommunicationNode, ComputationNode, SwapNode,
KernelNode, PModelNode, CodeTree, CodeTreeGenerator.
"""

from abc import abstractmethod
from copy import deepcopy
from enum import Enum
from re import sub
from typing import List, Tuple

import attr
from lark import Visitor, tree as _tree, lexer as _lexer

import offsite.config
from offsite.codegen.codegen_util import indent, replace_var_with_factor, replace_incr_with_assign_op, \
    substitute_rhs_call, write_instrument_kernel_start, write_instrument_kernel_end
from offsite.descriptions.ode_method import ODEMethod
from offsite.evaluation.math_utils import eval_math_expr


class CodeNodeType(Enum):
    LOOP = 'LOOP'
    COMPUTATION = 'COMPUTATION'
    COMMUNICATION = 'COMMUNICATION'
    KERNEL = 'KERNEL'
    PMODEL = 'PMODEL'
    SWAP = 'SWAP'
    ROOT = 'ROOT'


@attr.s
class CodeNode:
    """Representation of a code node object.

    Attributes:
    -----------
    indent : int
        Indent level of this object.
    type: CodeNodeType
        Node type of this object
    prev: CodeNode
        Left sibling node of this object.
    next: CodeNode
        Right sibling node of this object.
    parent: CodeNode
        Parent node of this object.
    child: CodeNode
        Child node of this object.
    """
    indent = attr.ib(type=int)
    type = attr.ib(type=CodeNodeType)
    prev = attr.ib(type='CodeNode')
    next = attr.ib(type='CodeNode')
    parent = attr.ib(type='CodeNode')
    child = attr.ib(type='CodeNode')

    def set_relatives(self, prev: 'CodeNode', nxt: 'CodeNode', parent: 'CodeNode', child: 'CodeNode'):
        """Set relations of this object.

        Parameters:
        -----------
        prev: CodeNode
            New left sibling node.
        nxt: CodeNode
            New right sibling node.
        parent: CodeNode
            New parent sibling node.
        child: CodeNode
            New child sibling node.

        Returns:
        --------
        -
        """
        self.prev = prev
        self.next = nxt
        self.parent = parent
        self.child = child

    @abstractmethod
    def to_codeline(self) -> str:
        pass

    @abstractmethod
    def to_implementation_codeline(self) -> str:
        pass

    @abstractmethod
    def to_kerncraft_codeline(self) -> str:
        pass

    @abstractmethod
    def to_yasksite_codeline(self) -> str:
        pass


@attr.s
class RootNode(CodeNode):
    """Representation of a root node object.

    Attributes:
    -----------
    indent : int
        Indent level of this object.
    type: CodeNodeType
        Node type of this object
    prev: CodeNode
        Left sibling node of this object.
    next: CodeNode
        Right sibling node of this object.
    parent: CodeNode
        Parent node of this object.
    child: CodeNode
        Child node of this object.
    name: str
        Name identifier of the code tree.
    """
    name = attr.ib(str)

    def to_codeline(self) -> str:
        """Write node content to code line.

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Written code line.
        """
        return ''

    def to_implementation_codeline(self) -> str:
        """Write node content to implementation code line.

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Written code line.
        """
        return self.to_codeline()

    def to_kerncraft_codeline(self) -> str:
        """Write node content to kerncraft code line.

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Written code line.
        """
        return self.to_codeline()

    def to_kernel_codeline(self) -> str:
        """Write node content to kernel code line.

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Written code line.
        """
        return self.to_codeline()

    def to_yasksite_codeline(self) -> str:
        """Write node content to yasksite code line.

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Written code line.
        """
        return self.to_codeline()


@attr.s
class LoopNode(CodeNode):
    """Representation of a loop node object.

    Attributes:
    -----------
    indent : int
        Indent level of this object.
    type: CodeNodeType
        Node type of this object
    prev: CodeNode
        Left sibling node of this object.
    next: CodeNode
        Right sibling node of this object.
    parent: CodeNode
        Parent node of this object.
    child: CodeNode
        Child node of this object.
    var: str
        Loop variable.
    start: str
        Start value of the loop variable.
    boundary: str
        Break condition of the the loop. (assuming < condition).
    optional_args: list of str
        Optional loop arguments (unroll, assign, ...)
    optimize_unroll: int
        Global unroll id used to determine in which order loop nests are unrolled.
    optimize_assign: int
        Turn all statements in iteration 'optimize_assign' into assignment statements.
    pragmas: list of str
        Pragmas attached to this loop object.
    flag: bool
        Helper flag used during code generation (e.g. to skip loops)
    """
    var = attr.ib(type=str)
    start = attr.ib(type=str)
    boundary = attr.ib(type=str)

    optional_args = attr.ib(type=List['str'], default=list())
    optimize_unroll = attr.ib(type=int, default=None)
    optimize_assign = attr.ib(type=int, default=None)

    pragmas = attr.ib(type=List['str'], default=list())

    flag = attr.ib(type=bool, default=False, init=False)

    def __attrs_post_init__(self):
        """Parse optional arguments during object creation.

        Parameters:
        -----------
        -

        Returns:
        --------
        -
        """
        num_args = len(self.optional_args)
        # Parse optional arguments.
        if self.optional_args:
            if self.optional_args[0] == 'unroll':
                if num_args == 1:
                    self.optimize_unroll = 0
                else:
                    self.optimize_unroll = int(self.optional_args[1])
            else:
                raise RuntimeError('Invalid argument for loop node: ' + self.optional_args[0])
            if num_args > 2:
                if self.optional_args[2] == 'assign':
                    if num_args > 3:
                        self.optimize_assign = int(self.optional_args[3])
                    else:
                        raise RuntimeError('Missing argument! \'assign\' requires iteration number (e.g. assign 0): ')
            if num_args > 4:
                raise RuntimeError('Too many arguments (max==4) for loop node: ' + self.optional_args)
        # Strip empty space from members.
        self.var = self.var.strip()
        if isinstance(self.boundary, str):
            self.boundary = self.boundary.strip()

    def to_codeline(self) -> str:
        """Write node content to code line.

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Written code line.
        """
        pragmas = '\n'.join(self.pragmas) + ('\n' if self.pragmas else '')
        return pragmas + indent(self.indent) + 'for (int ' + self.var + '=' + str(self.start) + '; ' + self.var + '<' + \
               self.boundary + '; ++' + self.var + ') {' + '\n'

    def to_implementation_codeline(self) -> str:
        """Write node content to implementation code line.

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Written code line.
        """
        return self.to_codeline()

    def to_kerncraft_codeline(self) -> str:
        """Write node content to kerncraft code line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code line.
         """
        return self.to_codeline()

    def to_kernel_codeline(self) -> str:
        """Write node content to kernel code line.

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Written code line.
        """
        config = offsite.config.offsiteConfig
        # Adjust loop boundaries if the loop runs over the system dimension.
        if self.var == config.var_idx:
            first = config.var_first_idx
            last = config.var_last_idx
            comparator = '<='
        else:
            first = self.start
            last = self.boundary
            comparator = '<'
        # Write pragmas
        pragmas = '\n'.join(self.pragmas) + ('\n' if self.pragmas else '')
        # Write code
        return pragmas + indent(self.indent) + 'for (int ' + self.var + '=' + str(
            first) + '; ' + self.var + comparator + \
               last + '; ++' + self.var + ') {' + '\n'

    def to_tiling_kernel_codeline(self, block_varname: str) -> str:
        """Write node content to tiling kernel code line.

        Parameters:
        -----------
        block_varname: str
        Name suffix of the block size variable.

        Returns:
        --------
        str
            Written code line.
        """
        config = offsite.config.offsiteConfig
        # We only support tiling the system dimension loop right now.
        assert self.var == config.var_idx
        # Write pragmas.
        pragmas = '\n'.join(self.pragmas) + ('\n' if self.pragmas else '')
        # Write code.
        bs_var = 'bs_{}'.format(block_varname)
        return pragmas + indent(
            self.indent) + 'for (int ' + self.var + '= jj; ' + self.var + ' < ' + bs_var + ' + jj; ++' + self.var + \
               ') {' + '\n'

    def to_yasksite_codeline(self) -> str:
        """Write node content to yasksite code line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code line.
         """
        raise NotImplementedError('Yasksite code generation does not support loop nodes!')


@attr.s
class CommunicationNode(CodeNode):
    """Representation of a communication node object.

    Attributes:
    -----------
    indent : int
        Indent level of this object.
    type: CodeNodeType
        Node type of this object
    prev: CodeNode
        Left sibling node of this object.
    next: CodeNode
        Right sibling node of this object.
    parent: CodeNode
        Parent node of this object.
    child: CodeNode
        Child node of this object.
    operation: str
        Name of the communication operation.
    input_var_name: str
        Name of the used input vector variable.
    """
    operation = attr.ib(type=str)
    input_var_name = attr.ib(type=str, default=None)

    def to_codeline(self, ivp: 'IVP' = None) -> str:
        """Write node content to code line.

         Parameters:
         -----------
         ivp: IVP
            Used IVP.

         Returns:
         --------
         str
             Written code line.
         """
        if self.operation == 'omp_barrier':
            return '#pragma omp barrier \n'
        elif self.operation == 'mpi_allgather':
            assert self.input_var_name is not None
            code = '#pragma omp master\n'
            code += '{\n'
            code += 'MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, {}, NP, MPI_DOUBLE, MPI_COMM_WORLD);'.format(
                self.input_var_name)
            code += '}\n'
            return code
        elif self.operation == 'mpi_communicate':
            assert self.input_var_name is not None
            code = '#pragma omp master\n'
            code += '{\n'
            code += 'communicate({});'.format(self.input_var_name)
            code += '}\n'
            return code
        raise RuntimeError('Communication node contains unknown/unsupported communication operation!')

    def to_implementation_codeline(self, ivp: 'IVP' = None) -> str:
        """Write node content to implementation code line.

         Parameters:
         -----------
         ivp: IVP
            Used IVP.

         Returns:
         --------
         str
             Written code line.
         """
        return self.to_codeline(ivp)

    def to_kerncraft_codeline(self) -> str:
        """Write node content to kerncraft code line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code line.
         """
        raise NotImplementedError('Kerncraft code generation does not support communication nodes!')

    def to_kernel_codeline(self) -> str:
        """Write node content to kernel code line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code line.
         """
        raise NotImplementedError('Kernel code generation does not support communication nodes!')

    def to_yasksite_codeline(self) -> str:
        """Write node content to yasksite code line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code line.
         """
        raise NotImplementedError('Yasksite code generation does not support communication nodes!')


@attr.s
class ComputationNode(CodeNode):
    """Representation of a computation node object.

    Attributes:
    -----------
    indent : int
        Indent level of this object.
    type: CodeNodeType
        Node type of this object
    prev: CodeNode
        Left sibling node of this object.
    next: CodeNode
        Right sibling node of this object.
    parent: CodeNode
        Parent node of this object.
    child: CodeNode
        Child node of this object.
    computation: str
        Computation (statement, function call, ...).
    isIVPdepent: bool
        Computation contains RHS evaluations.
    """
    computation = attr.ib(type=str)
    isIVPdependent = attr.ib(type=bool, default=False)  # TODO deleteable?!

    def __attrs_post_init__(self):
        """Strip empty spaces from computation during object creation.

        Parameters:
        -----------
        -

        Returns:
        --------
        -
        """
        # Strip empty space from attribute computation.
        self.computation = self.computation.strip()

    def to_codeline(self) -> str:
        """Write node content to code line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code line.
         """
        # for debugging purposes only (--> used in CodeTree.visit())
        return indent(self.indent) + self.computation + ';' + '\n'

    def to_implementation_codeline(self) -> str:
        """Write node content to implementation code line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code line.
         """
        raise NotImplementedError('Implementation code generation does not support computation nodes!')

    def to_kerncraft_codeline(self, rhs: str = None, ivp_constants: List[Tuple[str, str]] = None) -> str:
        """Write node content to kerncraft code line.

         Parameters:
         -----------
         rhs: str
            Executed IVP component.
         ivp_constants: list of tuple(str, str)
            Available IVP constants.

         Returns:
         --------
         str
             Written code line.
         """
        if rhs and ivp_constants:
            computation = substitute_rhs_call(self.computation, rhs, ivp_constants)
        else:
            computation = self.computation
        return indent(self.indent) + computation + ';' + '\n'

    def to_kernel_codeline(self) -> str:
        """Write node content to kernel code line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code line.
         """
        return self.to_codeline()

    def to_yasksite_codeline(self) -> str:
        """Write node content to yasksite code line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code line.
         """
        return self.to_codeline()


@attr.s
class SwapNode(CodeNode):
    """Representation of a pointer swap node object.

    Attributes:
    -----------
    indent : int
        Indent level of this object.
    type: CodeNodeType
        Node type of this object
    prev: CodeNode
        Left sibling node of this object.
    next: CodeNode
        Right sibling node of this object.
    parent: CodeNode
        Parent node of this object.
    child: CodeNode
        Child node of this object.
    arg1: str
        First argument of the command.
    arg2: str
        Second argument of the command.
    datatype: str
        Used datatype.
    """
    arg1 = attr.ib(str)
    arg2 = attr.ib(str)
    datatype = attr.ib(str)

    def to_codeline(self) -> str:
        """Write node content to code line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code line.
         """
        string = '\n'
        string += '#pragma omp master\n'
        string += '{\n'
        string += indent(self.indent) + self.datatype + '** tmp = ' + self.arg1 + ' ;\n'
        string += indent(self.indent) + self.arg1 + ' = ' + self.arg2 + ' ;\n'
        string += indent(self.indent) + self.arg2 + ' = tmp;\n'
        string += '}\n'
        string += '\n'
        return string

    def to_implementation_codeline(self) -> str:
        """Write node content to implementation code line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code line.
         """
        return self.to_codeline()

    def to_kerncraft_codeline(self) -> str:
        """Write node content to kerncraft code line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code line.
         """
        raise NotImplementedError('Kerncraft code generation does not support swap nodes!')

    def to_kernel_codeline(self) -> str:
        """Write node content to kernel code line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code line.
         """
        raise NotImplementedError('Kernel code generation does not support swap nodes!')

    def to_yasksite_codeline(self) -> str:
        """Write node content to yasksite code line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code line.
         """
        raise NotImplementedError('Yasksite code generation does not support swap nodes!')


@attr.s
class KernelNode(CodeNode):
    """Representation of a kernel node object.

    Attributes:
    -----------
    indent : int
        Indent level of this object.
    type: CodeNodeType
        Node type of this object
    prev: CodeNode
        Left sibling node of this object.
    next: CodeNode
        Right sibling node of this object.
    parent: CodeNode
        Parent node of this object.
    child: CodeNode
        Child node of this object.
    template_name: str
        Name of the associated kernel template.
    input_var_name: str
        Name of the used input vector variable.
    """
    template_name = attr.ib(str)
    input_var_name = attr.ib(str)

    def to_codeline(self) -> str:
        """Write node content to code line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code line.
         """
        return indent(self.indent) + 'KERNEL ' + self.template_name + '\n'

    def to_implementation_codeline(self, kernel_code: str, kernel_id: int, incl_instrumentation: bool = True) -> str:
        """Write node content to implementation code line.

         Parameters:
         -----------
         kernel_code: str
            Written kernel code.
        kernel_id: int
            Database id of the written kernel.

         Returns:
         --------
         str
             Written code line.
         """
        # Write code.
        code = '//{} %{}\n'.format(self.template_name, kernel_id)
        code += write_instrument_kernel_start(self.indent) if incl_instrumentation else ''
        code += kernel_code
        code += write_instrument_kernel_end(self.indent, kernel_id) if incl_instrumentation else ''
        return code

    def to_kerncraft_codeline(self) -> str:
        """Write node content to kerncraft code line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code line.
         """
        raise NotImplementedError('Kerncraft code generation does not support kernel nodes!')

    def to_kernel_codeline(self) -> str:
        """Write node content to kernel code line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code line.
         """
        raise NotImplementedError('Kernel code generation does not support kernel nodes!')

    def to_yasksite_codeline(self) -> str:
        """Write node content to yasksite code line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code line.
         """
        raise NotImplementedError('Yasksite code generation does not support kernel nodes!')


@attr.s
class PModelNode(CodeNode):
    """Representation of a pmodel node object.

    Attributes:
    -----------
    indent : int
        Indent level of this object.
    type: CodeNodeType
        Node type of this object
    prev: CodeNode
        Left sibling node of this object.
    next: CodeNode
        Right sibling node of this object.
    parent: CodeNode
        Parent node of this object.
    child: CodeNode
        Child node of this object.
    """

    def to_codeline(self) -> str:
        """Write node content to code line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code line.
         """
        return ''

    def to_implementation_codeline(self) -> str:
        """Write node content to implementation code line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code line.
         """
        self.to_codeline()

    def to_kerncraft_codeline(self) -> str:
        """Write node content to kerncraft code line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code line.
         """
        self.to_codeline()

    def to_kernel_codeline(self) -> str:
        """Write node content to kernel code line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code line.
         """
        self.to_codeline()

    def to_yasksite_codeline(self) -> str:
        """Write node content to yasksite code line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code line.
         """
        raise NotImplementedError('Yasksite code generation does not support pmodel nodes!')


@attr.s
class CodeTree:
    root = attr.ib(type=CodeNode, default=RootNode(0, CodeNodeType.ROOT, None, None, None, None))

    @staticmethod
    def substitute_butcher_coefficients(node: CodeNode, method: ODEMethod, use_dummy_values=False):
        """
        Recursively substitute all occurrences of butcher table entries in the passed code tree with the coefficient
        values provided by the given ODE method.

        Parameters:
        -----------
        node: CodeNode
            Root node of the code tree.
        method: ODE method
            Used ODE method.

        Returns:
        --------
        -
        """
        if node.type == CodeNodeType.COMPUTATION:
            # Replace A coefficients.
            for idx_row, A_row in enumerate(method.coefficientsA):
                for idx_col, elem in enumerate(A_row):
                    evaluated = eval_math_expr(elem, cast_to=str)
                    node.computation = node.computation.replace('A[{}][{}]'.format(idx_row, idx_col), evaluated)
            # Replace b coefficients.
            for idx, elem in enumerate(method.coefficientsB):
                evaluated = eval_math_expr(elem, cast_to=str)
                node.computation = node.computation.replace('b[{}]'.format(idx), evaluated)
            # Replace c coefficients.
            for idx, elem in enumerate(method.coefficientsC):
                evaluated = eval_math_expr(elem, cast_to=str)
                node.computation = node.computation.replace('c[{}]'.format(idx), evaluated)
            # Replace not already substituted coefficients with dummy values. Helps enabling running some kernels
            # in kerncraft LC mode.
            if use_dummy_values:
                node.computation = sub(r'A\[[\w|\s|\d][\w|\s|\d|+|-|\*/]?\]\[[\w|\s|\d][\w|\s|\d|\+-|\*/]?\]',
                                       '0.123', node.computation)
                node.computation = sub(r'b\[[\w|\s|\d][\w|\s|\d|+|-|*|/]?\]', '0.456', node.computation)
                node.computation = sub(r'c\[[\w|\s|\d][\w|\s|\d|+|-|*|/]?\]', '0.789', node.computation)
        # Traverse tree depth first.
        if node.child:
            CodeTree.substitute_butcher_coefficients(node.child, method, use_dummy_values)
        if node.next:
            CodeTree.substitute_butcher_coefficients(node.next, method, use_dummy_values)

    @staticmethod
    def unroll_loop(loop: CodeNode):
        """Unroll the passed loop code tree.

        Parameters:
        -----------
        loop: CodeNode
            Root node of the code tree.
        skip_assign: bool
            If true do not apply assign optimization transformations.

        Returns:
        --------
        -
        """
        assert (loop.parent is not None) ^ (loop.prev is not None)
        assert loop.child is not None
        if loop.parent:
            parent = loop.parent
            parent.child = None
            # Cut loop body from original loop.
            loop_bdy = loop.child
            loop_bdy.parent = None
            # Update indention.
            CodeTree.update_indent(loop_bdy, -1)
            # Copy loop body and concatenate all unrolled loop bodies.
            parent.child = deepcopy(loop_bdy)
            parent.child.parent = parent
            # ..... unroll iteration 'start'.
            cur_loop_bdy = parent.child
            # Replace loop variable with current loop iteration index.
            cur = cur_loop_bdy
            while cur is not None:
                cur.computation = replace_var_with_factor(loop.var, cur.computation, loop.start)
                cur.computation = replace_incr_with_assign_op(cur.computation) if (
                        loop.optimize_assign is not None and loop.optimize_assign == int(
                    loop.start)) else cur.computation
                cur = cur.next
            # ..... unroll all remaining iterations.
            for idx in range(int(loop.start) + 1, int(loop.boundary)):
                loop_bdy_no_idx = deepcopy(loop_bdy)
                # Replace loop variable with current loop iteration index.
                cur = loop_bdy_no_idx
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
            sibling = loop.prev
            sibling.next = None
            # Cut loop body from original loop.
            loop_bdy = loop.child
            loop_bdy.parent = None
            # Update indention.
            CodeTree.update_indent(loop_bdy, -1)
            # Copy loop body and concatenate all unrolled loop bodies.
            sibling.next = deepcopy(loop_bdy)
            sibling.next.prev = sibling
            # ..... unroll iteration 'start'.
            cur_loop_bdy = sibling.next
            # Replace loop variable with current loop iteration index.
            cur = cur_loop_bdy
            while cur is not None:
                cur.computation = replace_var_with_factor(loop.var, cur.computation, loop.start)
                cur.computation = replace_incr_with_assign_op(cur.computation) if (
                        loop.optimize_assign is not None and loop.optimize_assign == int(
                    loop.start)) else cur.computation
                cur = cur.next
            # ..... unroll iteration 'start'.
            for idx in range(int(loop.start) + 1, int(loop.boundary)):
                loop_bdy_no_idx = deepcopy(loop_bdy)
                # Replace loop variable with current loop iteration index.
                cur = loop_bdy_no_idx
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
    def split_loop(loop: CodeNode):
        """Split the passed loop code tree.

        Parameters:
        -----------
        loop: CodeNode
            Root node of the code tree.

        Returns:
        --------
        -
        """
        assert (loop.parent is not None) ^ (loop.prev is not None)
        assert loop.child is not None
        assert loop.split_at == 0
        # TODO: at the moment only supports splitting the first iteration
        if loop.parent:
            parent = loop.parent
            parent.child = None
            loop.parent = None
            # Copy loop body of the original loop.
            loop_bdy = deepcopy(loop.child)
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
            sibling = loop.prev
            sibling.next = None
            loop.prev = None
            # Copy loop body of the original loop.
            loop_bdy = deepcopy(loop.child)
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
        """Traverse the passed loop code tree and print information all visited nodes.

        Parameters:
        -----------
        node: CodeNode
            Root node of the code tree.

        Returns:
        --------
        -
        """
        print(node.to_codeline())
        if node.child:
            CodeTree.visit(node.child)
        if node.next:
            CodeTree.visit(node.next)

    @staticmethod
    def replace_loop_var_with_iteration_idx(node: CodeNode, loop: LoopNode):
        """
        Replace within a passed code tree, all occurrences of a given loop variable with a given loop iteration index.

        Parameters:
        -----------
        node: CodeNode
            Root node of the code tree.
        loop: CodeNode
            Used loop node that determines the loop variable name and loop iteration index.

        Returns:
        --------
        -
        """
        if node.type == CodeNodeType.COMPUTATION:
            node.computation = replace_var_with_factor(loop.var, node.computation, str(loop.split_at))
            node.computation = replace_incr_with_assign_op(node.computation) if loop.split_at == 0 else node.computation
        if node.child:
            CodeTree.replace_loop_var_with_iteration_idx(node.child, loop)
        if node.next:
            CodeTree.replace_loop_var_with_iteration_idx(node.next, loop)

    @staticmethod
    def iteration_count(node: CodeNode, iter_count: str = '') -> str:
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
        if node.type == CodeNodeType.LOOP and 'unroll' not in node.optional_args:
            if iter_count:
                iter_count = '{} * ({} - {})'.format(iter_count, node.boundary, node.start)
            else:
                iter_count = '({} - {})'.format(node.boundary, node.start)
        if node.child:
            iter_count = CodeTree.iteration_count(node.child, iter_count)
        if node.next and node.type is not CodeNodeType.PMODEL:
            iter_count = CodeTree.iteration_count(node.next, iter_count)
        return iter_count

    @staticmethod
    def iteration_count_up_to_root(node: CodeNode, iter_count: str = '') -> str:
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
        if node.type == CodeNodeType.LOOP and 'unroll' not in node.optional_args:
            if iter_count:
                iter_count = '{} * ({} - {})'.format(iter_count, node.boundary, node.start)
            else:
                iter_count = '({} - {})'.format(node.boundary, node.start)
        if node.parent:
            iter_count = CodeTree.iteration_count(node.parent, iter_count)
        return iter_count

    @staticmethod
    def find_pmodel_node(node: CodeNode, requested_idx: int, nxt_visited_idx=0,
                         found_node: CodeNode = None) -> PModelNode:
        """Traverse node's subtree and return the 'requested_idx' occurrence of a PModelNode object.

        Parameters:
        -----------
        node: CodeNode
            Root node of the code tree.
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
            found_node, nxt_visited_idx = CodeTree.find_pmodel_node(node.child, requested_idx, nxt_visited_idx,
                                                                    found_node)
        if node.next:
            found_node, nxt_visited_idx = CodeTree.find_pmodel_node(node.next, requested_idx, nxt_visited_idx,
                                                                    found_node)
        return found_node, nxt_visited_idx

    @staticmethod
    def update_indent(node: CodeNode, change_by_val: int):
        """Update indention of the passed code tree by a given value.

        Parameters:
        -----------
        loop: CodeNode
            Root node of the code tree.
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

        self.last_visited_node = self.code_tree.root
        self.switch_to_nxt_lvl = False
        self.current_loop_lvl = None
        self.current_indent_lvl = 0
        self.current_pragmas = list()

    def get_relatives(self) -> Tuple[CodeNode, CodeNode, CodeNode, CodeNode]:
        if self.switch_to_nxt_lvl:
            prev = None
            parent = self.current_loop_lvl
        else:
            prev = self.last_visited_node
            parent = None
        nxt = None
        child = None
        return prev, nxt, parent, child

    def move_up_loop_lvl(self, node: CodeNode) -> CodeNode:
        self.current_indent_lvl -= 1
        #
        parent = node.parent
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
    def read_loop_s(tree: 'lark.tree.Tree'):
        prev_is_pragma = False
        base_idx = 1
        # Read pragmas.
        for entry in tree.children:
            if isinstance(entry, _tree.Tree) and entry.data == 'pragma':
                base_idx += 1
                prev_is_pragma = True
            else:
                if prev_is_pragma and isinstance(entry, _lexer.Token) and entry.value == '\n':
                    base_idx += 1
                prev_is_pragma = False
        # Read variable and boundary.
        var = str(tree.children[base_idx])
        bound = ''.join(str(x) for x in tree.children[base_idx + 1].children)
        # Read further optional arguments.
        opt_args = [str(entry) for entry in tree.children[(base_idx + 2):]]
        return var, bound, opt_args

    def loop_s(self, tree: 'lark.tree.Tree'):
        # Read tree information.
        loop_var, loop_boundary, optional_args = self.read_loop_s(tree)
        # Create loop node.
        prev, nxt, parent, child = self.get_relatives()
        node = LoopNode(self.current_indent_lvl, CodeNodeType.LOOP, prev, nxt, parent, child, loop_var, str(0),
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
        self.current_indent_lvl += 1

    def loop_e(self, tree: 'lark.tree.Tree'):
        # Read tree information.
        loop_var = str(tree.children[1])
        # Check if loop nesting is valid.
        assert loop_var == self.current_loop_lvl.var
        # Update global information: last visited node.
        self.last_visited_node = self.current_loop_lvl
        # Move up a level in the loop hierarchy.
        self.current_loop_lvl = self.move_up_loop_lvl(self.current_loop_lvl)

    def pragma(self, tree: 'lark.tree.Tree'):
        # Read tree information.
        pragma_id = str(tree.children[1])
        # Create pragma string.
        self.current_pragmas.append('#pragma ' + pragma_id)

    def comp(self, tree: 'lark.tree.Tree'):
        # Read tree information.
        comp_id = str(tree.children[1])
        # Create computation node.
        prev, nxt, parent, child = self.get_relatives()
        node = ComputationNode(self.current_indent_lvl, CodeNodeType.COMPUTATION, prev, nxt, parent, child, comp_id)
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

    def pmodel(self, tree: 'lark.tree.Tree'):
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

    def comm(self, tree: 'lark.tree.Tree'):
        # Read tree information.
        comm_id = str(tree.children[1])
        # Read optional input vector name argument.
        input_vec = tree.children[2] if len(tree.children) > 2 else None
        if input_vec is not None:
            if isinstance(input_vec, _tree.Tree):
                input_vec = '{}[{}]'.format(input_vec.children[0], input_vec.children[1])
            else:
                input_vec = str(input_vec)
        # Create computation node.
        prev, nxt, parent, child = self.get_relatives()
        node = CommunicationNode(
            self.current_indent_lvl, CodeNodeType.COMMUNICATION, prev, nxt, parent, child, comm_id, input_vec)
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

    def kernel(self, tree: 'lark.tree.Tree'):
        # Read tree information.
        kernel_id = str(tree.children[1])
        input_vec = tree.children[2] if len(tree.children) > 2 else None
        if input_vec is not None:
            if isinstance(input_vec, _tree.Tree):
                input_vec = '{}[{}]'.format(input_vec.children[0], input_vec.children[1])
            else:
                input_vec = str(input_vec)
        # Create computation node.
        prev, nxt, parent, child = self.get_relatives()
        node = KernelNode(self.current_indent_lvl, CodeNodeType.KERNEL, prev, nxt, parent, child, kernel_id, input_vec)
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

    def swap(self, tree: 'lark.tree.Tree'):
        # Read tree information.
        arg1 = str(tree.children[1])
        arg2 = str(tree.children[2])
        datatype = str(tree.children[3])
        # Create computation node.
        prev, nxt, parent, child = self.get_relatives()
        node = SwapNode(self.current_indent_lvl, CodeNodeType.SWAP, prev, nxt, parent, child, arg1, arg2, datatype)
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

    def generate(self, lark_tree: 'lark.tree.Tree') -> CodeTree:
        self.visit_topdown(lark_tree)
        return self.code_tree
