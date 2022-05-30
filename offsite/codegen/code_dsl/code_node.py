"""@package codegen.code_dsl.code_node
Definitions of classes CodeNodeType, CodeNode, RootNode, LoopNode, CommunicationClustLvlNode, CommunicationNodeLvlNode,
ComputationNode, SwapNode, KernelNode, PModelNode.

@author: Johannes Seiferth
"""

from abc import abstractmethod
from enum import Enum
from typing import Dict, List, Optional, Tuple, Union

import attr

import offsite.config
from offsite.codegen.codegen_util import indent, substitute_rhs_call, write_instrument_kernel_start, \
    write_instrument_kernel_end
from offsite.config import Config


class CodeNodeType(Enum):
    LOOP = 'LOOP'
    COMPUTATION = 'COMPUTATION'
    COMMUNICATION_NODE = 'COMMUNICATION_NODE'
    COMMUNICATION_CLUST = 'COMMUNICATION_CLUST'
    KERNEL = 'KERNEL'
    PMODEL = 'PMODEL'
    SWAP = 'SWAP'
    ROOT = 'ROOT'

    def __eq__(self, other):
        return self.value == other.value


@attr.s
class CodeNode:
    """Representation of a code_dsl node object.

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

    def set_relatives(self, prev: Optional['CodeNode'], nxt: Optional['CodeNode'], parent: Optional['CodeNode'],
                      child: Optional['CodeNode']):
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

    def my_parent(self) -> Optional['CodeNode']:
        if self.type == CodeNodeType.ROOT:
            return None
        if self.parent is not None:
            return self.parent
        assert self.prev is not None
        return self.prev.my_parent()

    @abstractmethod
    def to_codeline(self) -> str:
        pass

    @abstractmethod
    def to_implementation_codeline(self) -> str:
        pass

    @abstractmethod
    def to_kernel_codeline(self) -> str:
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
        Name identifier of the code_dsl tree.
    """
    name = attr.ib(str)

    def to_codeline(self) -> str:
        """Write node content to code_dsl line.

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Written code_dsl line.
        """
        return ''

    def to_implementation_codeline(self) -> str:
        """Write node content to implementation code_dsl line.

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Written code_dsl line.
        """
        return self.to_codeline()

    def to_kerncraft_codeline(self) -> str:
        """Write node content to kerncraft code_dsl line.

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Written code_dsl line.
        """
        return self.to_codeline()

    def to_kernel_codeline(self) -> str:
        """Write node content to kernel code_dsl line.

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Written code_dsl line.
        """
        return self.to_codeline()

    def to_yasksite_codeline(self) -> str:
        """Write node content to yasksite code_dsl line.

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Written code_dsl line.
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
        Break condition of the loop. (assuming < condition).
    optional_args: list of str
        Optional loop arguments (unroll, assign, ...)
    optimize_unroll: int
        Global unroll id used to determine in which order loop nests are unrolled.
    optimize_assign: int
        Turn all statements in iteration 'optimize_assign' into assignment statements.
    pragmas: list of str
        Pragmas attached to this loop object.
    flag: bool
        Helper flag used during code_dsl generation (e.g. to skip loops)
    """
    var = attr.ib(type=str)
    start = attr.ib(type=str)
    boundary = attr.ib(type=str)

    optional_args = attr.ib(type=List['str'], default=list())
    optimize_unroll = attr.ib(type=int, default=None)
    optimize_assign = attr.ib(type=int, default=None)

    pragmas = attr.ib(type=List['str'], default=list())

    flag = attr.ib(type=bool, default=False, init=False)

    split_at = attr.ib(type=int, default=None, init=False)

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
                raise RuntimeError('Too many arguments (max==4) for loop node: ' + ';'.join(self.optional_args))
        # Strip empty space from members.
        self.var = self.var.strip()
        if isinstance(self.boundary, str):
            self.boundary = self.boundary.strip()

    def to_codeline(self) -> str:
        """Write node content to code_dsl line.

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Written code_dsl line.
        """
        pragmas = '\n'.join(self.pragmas) + ('\n' if self.pragmas else '')
        return pragmas + indent(self.indent) + 'for (int ' + self.var + '=' + str(
            self.start) + '; ' + self.var + '<' + self.boundary + '; ++' + self.var + ') {' + '\n'

    def to_implementation_codeline(self) -> str:
        """Write node content to implementation code_dsl line.

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Written code_dsl line.
        """
        return self.to_codeline()

    def to_kerncraft_codeline(self) -> str:
        """Write node content to kerncraft code_dsl line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code_dsl line.
         """
        return self.to_codeline()

    def to_kernel_codeline(self) -> str:
        """Write node content to kernel code_dsl line.

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Written code_dsl line.
        """
        config: Config = offsite.config.offsiteConfig
        # Adjust loop boundaries if the loop runs over the system dimension.
        if self.var == config.var_idx:
            first: str = 'first'
            last: str = 'last'
            comparator = '<='
        else:
            first: str = self.start
            last: str = self.boundary
            comparator = '<'
        # Write pragmas
        pragmas = '\n'.join(self.pragmas) + ('\n' if self.pragmas else '')
        # Write code_dsl
        code = pragmas + indent(self.indent)
        code += 'for (int ' + self.var + '=' + str(
            first) + '; ' + self.var + comparator + last + '; ++' + self.var + ') {' + '\n'
        return code

    def to_tiling_kernel_codeline(self, block_varname: str) -> str:
        """Write node content to tiling kernel code_dsl line.

        Parameters:
        -----------
        block_varname: str
        Name suffix of the block size variable.

        Returns:
        --------
        str
            Written code_dsl line.
        """
        config: Config = offsite.config.offsiteConfig
        # We only support tiling the system dimension loop right now.
        assert self.var == config.var_idx
        # Write pragmas.
        pragmas = '\n'.join(self.pragmas) + ('\n' if self.pragmas else '')
        # Write code_dsl.
        bs_var = 'bs_{}'.format(block_varname)
        code = pragmas + indent(self.indent)
        code += 'for (int ' + self.var + '= jj; ' + self.var + ' < ' + bs_var + ' + jj; ++' + self.var + ') {' + '\n'
        return code

    def to_yasksite_codeline(self) -> str:
        """Write node content to yasksite code_dsl line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code_dsl line.
         """
        raise NotImplementedError('Yasksite code_dsl generation does not support loop nodes!')


@attr.s
class CommunicationClustLvlNode(CodeNode):
    """Representation of a cluster-level communication node object.

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

    def to_codeline(self) -> str:
        """Write node content to code line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code_dsl line.
         """
        if self.operation == 'mpi_allgather':
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
        raise RuntimeError('Cluster-level communication node contains unknown/unsupported communication operation!')

    def to_implementation_codeline(self) -> str:
        """Write node content to implementation code_dsl line.

         Parameters:
         -----------
         -

         Returns:
         --------
         str
             Written code_dsl line.
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
        # Skipped when generating kerncraft code.
        return ''

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
        raise NotImplementedError('Yasksite code generation does not support cluster-level communication nodes!')


@attr.s
class CommunicationNodeLvlNode(CodeNode):
    """Representation of a node-level communication node object.

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
        if self.operation == 'omp_barrier':
            return '#pragma omp barrier \n'
        raise RuntimeError('Node-level communication node contains unknown/unsupported communication operation!')

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
        # Skipped when generating kerncraft code.
        return ''

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
        raise NotImplementedError('Yasksite code generation does not support node-level communication nodes!')


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
    references: List[dict]
        Array references included in this computation.
    isIVPdependent: bool
        Computation contains RHS evaluations.
    """
    computation = attr.ib(type=str)
    references = attr.ib(type=List[Dict])
    isIVPdependent = attr.ib(type=bool, default=False)

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
        # For debugging purposes only (used by CodeTree.visit()).
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

    def to_kerncraft_codeline(self, rhs: Optional[str] = None, rhs_input: Optional[str] = None,
                              ivp_constants: Optional[List[Tuple[str, Union[str, float, int]]]] = None) -> str:
        """Write node content to kerncraft code line.

         Parameters:
         -----------
         rhs: str
            Executed IVP component.
        rhs_input: str
            Used IVP input vector.
         ivp_constants: list of tuple(str, str)
            Available IVP constants.

         Returns:
         --------
         str
             Written code line.
         """
        if rhs and ivp_constants:
            computation: str = substitute_rhs_call(self.computation, rhs, rhs_input, ivp_constants)
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
        string += indent(self.indent) + self.datatype + '** swp_tmp = ' + self.arg1 + ' ;\n'
        string += indent(self.indent) + self.arg1 + ' = ' + self.arg2 + ' ;\n'
        string += indent(self.indent) + self.arg2 + ' = swp_tmp;\n'
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

    # noinspection PyMethodOverriding
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
        raise NotImplementedError('Yasksite code generation does not support pmodel nodes!')
