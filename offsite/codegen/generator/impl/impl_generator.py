"""@package codegen.generator.impl.impl_generator
Definition of abstract class ImplCodeGenerator.

@author: Johannes Seiferth
"""

from abc import ABC, abstractmethod
from copy import deepcopy
from typing import Dict, List, Optional, Tuple, Union

import attr
from pathlib2 import Path
from sortedcontainers import SortedDict

import offsite.config
from offsite.codegen.code_dsl.code_node import CodeNode, CodeNodeType, KernelNode, LoopNode
from offsite.codegen.code_dsl.code_tree import CodeTree
from offsite.codegen.codegen_util import eval_loop_boundary, write_closing_bracket, create_variant_name
from offsite.codegen.generator.kernel_generator import KernelCodeGenerator
from offsite.config import Config
from offsite.descriptions.impl.impl_skeleton import ImplSkeleton
from offsite.descriptions.impl.kernel_template import KernelTemplate, Kernel
from offsite.descriptions.ode import IVP, ODEMethod, corrector_steps, stages
from offsite.descriptions.parser import DatastructDict
from offsite.util.math_utils import eval_math_expr


@attr.s
class ImplCodeGenerator(ABC):
    """Representation of the base class for implementation variant code generators in different languages.

    Attributes:
    -----------
    db_session: sqlalchemy.orm.session.Session
        Used database session.
    folder_impl: Path
        Folder all created implementation variant code files are stored in.
    folder_ivp: Path
        Folder all created IVP problem code files are stored in.
    folder_method: Path
        Folder all created ODE method code files are stored in.
    unroll_stack: SortedDict
        Stack of loops to be unrolled.
    loaded_templates: dict (key=str, val=KernelTemplate)
        Required and already loaded kernel templates.
    required_datastructs: DatastructDict
        Required and used data structures.
    """
    db_session = attr.ib(type='sqlalchemy.orm.session.Session')
    folder_impl = attr.ib(type=Path)
    folder_ivp = attr.ib(type=Path)
    folder_method = attr.ib(type=Path)
    unroll_stack = attr.ib(type=SortedDict, init=False)
    loaded_templates = attr.ib(type=Dict, init=False)
    required_datastructs = attr.ib(type=DatastructDict, init=False)

    def generate(self, skeleton: ImplSkeleton, impl_variants: List[Tuple[int]], ivp: Optional[IVP] = None,
                 method: Optional[ODEMethod] = None) -> Dict[Path, str]:
        """Generate implementation variant codes.

        Parameters:
        -----------
        skeleton: ImplSkeleton
            ImplSkeleton for which code is generated.
        impl_variants: list of int tuples
            Implementation variants for which code is generated.
        ivp: IVP
            IVP for which code is generated.
        method: ODEMethod
            ODE method for which code is generated.

        Returns:
        --------
        dict of str
            Generated implementation variant codes.
        """
        config: Config = offsite.config.offsiteConfig
        # Set some members to match current code generation.
        code_tree: CodeNode = deepcopy(skeleton.code_tree).root
        code_tree.name = skeleton.name

        self.unroll_stack = SortedDict()
        self.unroll_stack.clear()
        self.loaded_templates = dict()
        self.loaded_templates.clear()
        self.required_datastructs = DatastructDict()
        self.required_datastructs.clear()

        # Create folder if it does not yet exist.
        if self.folder_impl and not self.folder_impl.exists():
            self.folder_impl.mkdir(parents=True)

        # Generate implementation code trees.
        self._generate_implementation_trees(code_tree, skeleton, method)
        # Optimize tree: unroll
        self._unroll_tree(code_tree)
        # Generate code from tree and write to string.
        tiled_str: str = 'tiled ' if config.args.tile else ''
        codes: Dict[Path, str] = dict()
        for vid, variant in impl_variants:
            # Generate required kernel codes.
            kernels: List[Kernel] = self._derive_impl_variant_kernels(variant)
            # Create implementation variant name.
            variant_name: str = create_variant_name(kernels, skeleton.name)
            if config.args.verbose:
                print('Generating {}implementation variant {} (id = {})'.format(tiled_str, variant_name, vid))
            # Generate implementation variant.
            name = Path('{}/{}.h'.format(self.folder_impl, create_variant_name(kernels, skeleton.name)))
            codes[name] = self._generate_impl_variant(code_tree, vid, kernels, skeleton, method, ivp, config.args.tile)
        return codes

    def _collect_loop_meta_data(self, node: CodeNode):
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
        if node.child:
            self._collect_loop_meta_data(node.child)
        if node.next:
            self._collect_loop_meta_data(node.next)

    def _unroll_tree(self, tree: CodeNode):
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
            self._collect_loop_meta_data(tree)
            if len(self.unroll_stack.keys()) == 0:
                break
            CodeTree.unroll_loop(self.unroll_stack.values()[0])

    def _generate_implementation_trees(self, node: CodeNode, skeleton: ImplSkeleton, method: ODEMethod):
        """Generate implementation variant code tree.

        Parameters:
        -----------
        node: CodeNode
            Root node of the code tree.
        skeleton: ImplSkeleton
            Used implementation skeleton.
        method: ODEMethod
            Used ODE method.

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
        elif node.type == CodeNodeType.KERNEL:
            node: KernelNode
            # Load kernel template from database.
            if node.template_name not in self.loaded_templates.keys():
                template: KernelTemplate = KernelTemplate.from_database(self.db_session, node.template_name)
                self.loaded_templates[node.template_name] = template
                self.required_datastructs = {**template.datastructs, **self.required_datastructs}
        # Traverse tree depth first.
        if node.child:
            self._generate_implementation_trees(node.child, skeleton, method)
        if node.next:
            self._generate_implementation_trees(node.next, skeleton, method)

    def _derive_impl_variant_kernels(self, impl_variant: List[int]) -> List[Kernel]:
        return [self._map_kernel_id_to_object(kid) for kid in impl_variant]

    def _map_kernel_id_to_object(self, kernel_id: int) -> 'Kernel':
        # Query kernel object from database.
        kernel: Kernel = Kernel.select(self.db_session, kernel_id)
        # Select corresponding kernel template.
        template: KernelTemplate = self.loaded_templates[kernel.template.name]
        try:
            return next(filter(lambda x: x.name == kernel.name, template.variants))
        except StopIteration:
            raise RuntimeError('Failed to find kernel \'{}\'!'.format(kernel.name))

    @abstractmethod
    def _generate_impl_variant(self, impl: CodeNode, variant_id: int, kernels: List[Kernel], skeleton: ImplSkeleton,
                               method: Optional[ODEMethod] = None, ivp: Optional[IVP] = None,
                               gen_tiled_code: bool = False) -> str:
        """Write implementation variant code.

        Parameters:
        -----------
        impl: CodeNode
            Root node of the code tree.
        variant_id: int
            Database id of the generated variant.
        kernels: list of Kernel
            Used kernels.
        skeleton: ImplSkeleton
            Used impl skeleton.
        method: ODEMethod
            Used ODE method.
        ivp: IVP
            Used IVP.
        gen_tiled_code: bool
            Generated tiled code version.

        Returns:
        --------
        str
            Generated impl variant code.
        """
        pass

    @staticmethod
    def _generate_implementation_code(node: CodeNode, kernels: List[Kernel], method: Optional[ODEMethod] = None,
                                      gen_tiled_code: bool = False, code: str = '') -> str:
        """Write impl variant code.

        Parameters:
        -----------
        node: CodeNode
            Root node of the code tree.
        kernels: list of Kernel
            Used kernels.
        method: ODEMethod
            Used ODE method.
        gen_tiled_code: bool
            Generated tiled code version.
        code: str
            Current code string.

        Returns:
        --------
        str
            Generated implementation variant code.
        """
        if node.type == CodeNodeType.KERNEL:
            node: KernelNode
            try:
                kernel: Kernel = next(filter(lambda x: x.template.name == node.template_name, kernels))
            except StopIteration:
                raise RuntimeError('Failed to find kernel template \'{}\'!'.format(node.template_name))
            # Determine the used input vector.
            input_vector = ''
            if kernel.template.isIVPdependent:
                if node.input_var_name:
                    input_vector = node.input_var_name
                elif 'RHS_input' in kernel.template.codegen:
                    input_vector = kernel.template.codegen['RHS_input']
                else:
                    raise RuntimeError(
                        'Kernel template {} is missing required entry \'RHS_input\'.'.format(kernel.template.name))
            # Create kernel code
            kernel_code = KernelCodeGenerator().generate(kernel, method, input_vector, gen_tiled_code)
            code += node.to_implementation_codeline(kernel_code, kernel.db_id)
        elif node.type != CodeNodeType.ROOT:
            code += node.to_implementation_codeline()
        # Traverse tree depth first.
        if node.child:
            code += ImplCodeGenerator._generate_implementation_code(node.child, kernels, method, gen_tiled_code)
        # Close opened brackets.
        if node.type == CodeNodeType.LOOP:
            code += write_closing_bracket(node.indent)
        if node.next:
            code += ImplCodeGenerator._generate_implementation_code(node.next, kernels, method, gen_tiled_code)
        return code

    @abstractmethod
    def _write_skeleton_includes(self, variant_id: int, skeleton: ImplSkeleton, method: Optional[ODEMethod] = None,
                                 ivp: Optional[IVP] = None) -> str:
        """Write header includes.

        Parameters:
        -----------
        variant_id: int
            Database ID associated to the generated impl variant.
        skeleton: ImplSkeleton
            Used impl skeleton.
        method: ODEMethod
            Used ODE method.
        ivp: IVP
            Used IVP.

        Returns:
        --------
        str
            Generated code.
        """
        pass

    @abstractmethod
    def _write_datastructures(self, skeleton: str):
        """Write data structures to file.

         Parameters:
         -----------
         skeleton: str
             Name of the used impl skeleton.

         Returns:
         --------
         -
         """
        pass

    @abstractmethod
    def _write_rhs_function(self, ivp: IVP):
        """Write RHS functions to file.

        Parameters:
        -----------
        ivp: IVP
            Used IVP.

        Returns:
        --------
        -
        """
        pass

    @abstractmethod
    def _write_ode_method(self, method: ODEMethod):
        """Write ODE method to file.

        Parameters:
        -----------
        method: ODEMethod
            Used ODE method.

        Returns:
        --------
        -
        """
        pass

    @staticmethod
    def _write_butcher_table(method: ODEMethod) -> str:
        """Write butcher table arrays.

        Parameters:
        -----------
        method: ODEMethod
            Used ODE method.

        Returns:
        --------
        str
            Generated code.
        """
        butcher = ImplCodeGenerator._write_2d_array('A', method.coefficientsA)
        butcher += ImplCodeGenerator._write_1d_array('b', method.coefficientsB)
        butcher += ImplCodeGenerator._write_1d_array('c', method.coefficientsC)
        return butcher + '\n'

    @staticmethod
    def _write_2d_array(name: str, data: List[List[str]]) -> str:
        """Write 2D double array.

        Parameters:
        -----------
        name: str
            Name of the array.
        data: List of list of str
            Array data.

        Returns:
        --------
        str
            Generated code.
        """
        array_str = 'double ' + name + '[' + str(len(data)) + ']' + '[' + str(len(data[0])) + ']' + ' = {'
        for row in data:
            array_str += '{'
            for elem in row:
                array_str += str(eval_math_expr(elem)) + ', '
            array_str = array_str[:-2]
            array_str += '}, '
        array_str = array_str[:-2]
        array_str += '};\n'
        return array_str

    @staticmethod
    def _write_1d_array(name: str, data: List[str]) -> str:
        """Write 1D double array.

        Parameters:
        -----------
        name: str
            Name of the array.
        data: List of str
            Array data.

        Returns:
        --------
        str
            Generated code.
        """
        array_str = 'double ' + name + '[' + str(len(data)) + ']' + ' = {'
        for elem in data:
            array_str += str(eval_math_expr(elem)) + ', '
        array_str = array_str[:-2]
        array_str += '};\n'
        return array_str
