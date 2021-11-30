"""@package codegen.generator.impl_generator_cpp
Definition of class ImplCodeGeneratorCPP.
"""

from copy import deepcopy
from pathlib import Path
from re import sub
from typing import Dict, List, Optional, Tuple, Union

import attr
from sortedcontainers import SortedDict

import offsite.config
from offsite.codegen.code_dsl.code_node import CodeNode, CodeNodeType, KernelNode, LoopNode
from offsite.codegen.code_dsl.code_tree import CodeTree
from offsite.codegen.codegen_util import eval_loop_boundary, format_codefile, indent, write_closing_bracket, \
    create_variant_name
from offsite.codegen.generator.kernel_generator import KernelCodeGenerator
from offsite.config import Config
from offsite.descriptions.impl.impl_skeleton import ImplSkeleton
from offsite.descriptions.impl.kernel_template import KernelTemplate, Kernel
from offsite.descriptions.ode import IVP, ODEMethod, corrector_steps, stages
from offsite.descriptions.parser import DatastructDict, DatastructType
from offsite.util.math_utils import eval_math_expr, solve_equation


@attr.s
class ImplCodeGeneratorCPP:
    """Representation of the code generator for implementation variant CPP code.

    Attributes:
    -----------
    db_session: sqlalchemy.orm.session.Session
        Used database session.
    folder_ds: Path
        Folder all created datastructure code files are stored in.
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
    folder_ds = attr.ib(type=Path)
    folder_impl = attr.ib(type=Path)
    folder_ivp = attr.ib(type=Path)
    folder_method = attr.ib(type=Path)
    unroll_stack = attr.ib(type=SortedDict, default=SortedDict())
    loaded_templates = attr.ib(type=Dict, default=dict())
    required_datastructs = attr.ib(type=DatastructDict, default=DatastructDict())

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

    def generate(self, skeleton: ImplSkeleton, impl_variants: List[Tuple[int]], ivp: IVP,
                 method: ODEMethod) -> Dict[Path, str]:
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

        self.unroll_stack.clear()
        self.loaded_templates.clear()
        self.required_datastructs.clear()

        # Create folder if it does not yet exist.
        if self.folder_impl and not self.folder_impl.exists():
            self.folder_impl.mkdir(parents=True)

        # Generate implementation code trees.
        self.generate_implementation_trees(code_tree, skeleton, method)
        # Optimize tree: unroll
        self.unroll_tree(code_tree)
        # Generate code from tree and write to string.
        codes: Dict[Path, str] = dict()
        for vid, variant in impl_variants:
            if config.args.verbose:
                print('Generating implementation variant ', end='')
            # Generate required kernel codes.
            kernels: List[Kernel] = self.derive_impl_variant_kernels(variant)
            # Create implementation variant name.
            variant_name: str = create_variant_name(kernels, skeleton.name)
            if config.args.verbose:
                print('{} (id = {}) ... '.format(variant_name, vid), end='')
            # Generate implementation variant.
            name = Path('{}/{}.c'.format(self.folder_impl, variant_name))
            codes[name] = self.generate_impl_variant(code_tree, vid, kernels, skeleton, method, ivp, False)
            # .. generate tiled version too?
            if config.args.tile:
                if config.args.verbose:
                    print('adding tiled version ... '.format(variant_name, vid), end='')
                name = Path('{}/{}_tiled.c'.format(self.folder_impl, variant_name))
                vid = 123456
                codes[name] = self.generate_impl_variant(code_tree, vid, kernels, skeleton, method, ivp, True)
            if config.args.verbose:
                print('done.')
        return codes

    def generate_implementation_trees(self, node: CodeNode, skeleton: ImplSkeleton, method: ODEMethod):
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
            self.generate_implementation_trees(node.child, skeleton, method)
        if node.next:
            self.generate_implementation_trees(node.next, skeleton, method)

    def derive_impl_variant_kernels(self, impl_variant: List[int]) -> List[Kernel]:
        return [self.map_kernel_id_to_object(kid) for kid in impl_variant]

    def map_kernel_id_to_object(self, kernel_id: int) -> 'Kernel':
        # Query kernel object from database.
        kernel: Kernel = Kernel.select(self.db_session, kernel_id)
        # Select corresponding kernel template.
        template: KernelTemplate = self.loaded_templates[kernel.template.name]
        try:
            return next(filter(lambda x: x.name == kernel.name, template.variants))
        except StopIteration:
            raise RuntimeError('Failed to find kernel \'{}\'!'.format(kernel.name))

    def generate_impl_variant(self, impl: CodeNode, variant_id: int, kernels: List[Kernel], skeleton: ImplSkeleton,
                              method: ODEMethod, ivp: IVP, gen_tiled_code: bool = False) -> str:
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
        # Generate impl variant code.
        code = ImplCodeGeneratorCPP.generate_implementation_code(impl, kernels, method, gen_tiled_code)
        #
        local_var_defs = ''
        for name, desc in self.required_datastructs.items():
            # Skip names reserved for Butcher table entries.
            if name in ('A', 'b', 'c'):
                continue
            # Skip names reserved for pre-defined variables.
            if name in ('h', 't', 'g', 'y'):
                continue
            if desc.struct_type != DatastructType.scalar:
                code = sub(r'\b{}\b'.format(name), 'ds->{}'.format(name), code)
            else:
                local_var_defs += 'double {};\n'.format(name)
        # Write frame code.
        frame = self.write_skeleton_includes(variant_id, skeleton, method, ivp)
        frame += ImplCodeGeneratorCPP.write_function_description(skeleton.name)
        if local_var_defs != '':
            frame += local_var_defs + '\n'
        frame += ImplCodeGeneratorCPP.write_instrument_impl(variant_id)
        frame += code
        frame += 'return STEP_SUCCESS;\n}\n'  # Concludes code function block.
        return frame

    @staticmethod
    def generate_implementation_code(
            node: CodeNode, kernels: List[Kernel], method: ODEMethod, gen_tiled_code: bool, code: str = '') -> str:
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
            Generated yasksite code.
        """
        config: Config = offsite.config.offsiteConfig
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
                    input_vector = config.ode_solution_vector
            # Create kernel code
            kernel_code: str = KernelCodeGenerator().generate(kernel, method, input_vector, gen_tiled_code)
            code += node.to_implementation_codeline(kernel_code, kernel.db_id)
        elif node.type != CodeNodeType.ROOT:
            code += node.to_implementation_codeline()
        # Traverse tree depth first.
        if node.child:
            code += ImplCodeGeneratorCPP.generate_implementation_code(node.child, kernels, method, gen_tiled_code)
        # Close opened brackets.
        if node.type == CodeNodeType.LOOP:
            code += write_closing_bracket(node.indent)
        if node.next:
            code += ImplCodeGeneratorCPP.generate_implementation_code(node.next, kernels, method, gen_tiled_code)
        return code

    def write_skeleton_includes(self, variant_id: int, skeleton: ImplSkeleton, method: ODEMethod, ivp: IVP) -> str:
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
        includes = ''
        # Add variant database ID as define.
        includes += '#define VARIANT_ID {}\n'.format(variant_id)
        includes += '\n'
        # IVP.
        if skeleton.isIVPdependent:
            assert IVP is not None
            includes += '#include <math.h>\n'
            includes += '#include <{}/RHS_{}.h>\n'.format(self.folder_ivp, ivp.name)
            self.write_rhs_function(ivp)
        # ODE method.
        if method is not None:
            includes += '#include <{}/ODE_{}.h>\n'.format(self.folder_method, method.name)
            self.write_ode_method(method)
        # Data structures.
        self.write_datastructures(skeleton.name)
        includes += '#include <{}/DS_{}.h>\n'.format(self.folder_ds, skeleton.name)
        # Instrumentation.
        includes += '#ifdef INSTRUMENT\n#include "timesnap.h"\n#endif\n'
        # Step cancellation.
        includes += '#ifdef STEP_CANCELLATION\n'
        includes += '#include <step_cancellation.hpp>\n'
        includes += '#endif\n'
        return includes + '\n'

    def write_datastructures(self, skeleton: str):
        """Write data structures to file.

         Parameters:
         -----------
         skeleton: str
             Name of the used impl skeleton.

         Returns:
         --------
         -
         """
        # Create folder if it does not yet exist.
        if self.folder_ds and not self.folder_ds.exists():
            self.folder_ds.mkdir(parents=True)
        # Write inclusion guard first.
        dstruct_str = '#pragma once\n\n'
        # Write data structures.
        ds_alloc_str = ''
        ds_free_str = ''
        ds_typedef_str = ''
        for name, desc in self.required_datastructs.items():
            # Skip names reserved for Butcher table entries.
            if name in ('A', 'b', 'c', 'y'):
                continue
            # Skip names reserved for pre-defined variables.
            if name in ('h', 't', 'g', 'y'):
                continue
            # Add to alloc, free and struct.
            if desc.struct_type == DatastructType.scalar:
                ds_typedef_str += 'double {};\n'.format(name)
            elif desc.struct_type == DatastructType.array1D:
                ds_alloc_str += 'ds->{} = alloc1d({});\n'.format(name, desc.size[0])
                ds_free_str += 'free1d(ds->{});\n'.format(name)
                ds_typedef_str += 'double *{};\n'.format(name)
            elif desc.struct_type == DatastructType.array2D:
                ds_alloc_str += 'ds->{} = alloc2d({}, {});\n'.format(name, desc.size[0], desc.size[1])
                ds_free_str += 'free2d(ds->{});\n'.format(name)
                ds_typedef_str += 'double **{};\n'.format(name)
            elif desc.struct_type == DatastructType.array3D:
                ds_alloc_str += 'ds->{} = alloc3d({}, {}, {});\n'.format(name, desc.size[0], desc.size[1], desc.size[2])
                ds_free_str += 'free3d(ds->{});\n'.format(name)
                ds_typedef_str += 'double **{};\n'.format(name)
            else:
                assert False
        # Write typedef struct.
        dstruct_str += 'typedef struct {\n'
        dstruct_str += ds_typedef_str
        dstruct_str += '} '
        dstruct_str += 'DS_{}_t;\n\n'.format(skeleton)
        # Write C++ guard and namespace.
        dstruct_str += '#ifdef __cplusplus\n'
        dstruct_str += 'namespace DS_' + skeleton + ' {\n'
        # Write allocate function.
        dstruct_str += 'inline DS_{}_t *allocate_datastructures(const int n, const int s) '.format(skeleton)
        dstruct_str += '{\n'
        dstruct_str += 'DS_{0}_t *ds = new DS_{0}_t;\n'.format(skeleton)
        dstruct_str += ds_alloc_str
        dstruct_str += 'return ds;\n'
        dstruct_str += '}\n\n'
        # Write free function.
        dstruct_str += 'inline void free_datastructures(DS_{}_t *ds) '.format(skeleton)
        dstruct_str += '{\n'
        dstruct_str += ds_free_str
        dstruct_str += '}\n\n'
        # Close namespace and C++ guard.
        dstruct_str += '}\n'
        dstruct_str += '#endif\n'
        # Write to file.
        path = Path('{}/DS_{}.h'.format(self.folder_ds, skeleton))
        with path.open('w') as file_handle:
            file_handle.write(dstruct_str)
        # Format code with indent tool if available.
        format_codefile(path)

    def write_rhs_function(self, ivp: IVP):
        """Write RHS functions to file.

        Parameters:
        -----------
        ivp: IVP
            Used IVP.

        Returns:
        --------
        -
        """
        config: Config = offsite.config.offsiteConfig
        # Create folder if it does not yet exist.
        if self.folder_ivp and not self.folder_ivp.exists():
            self.folder_ivp.mkdir(parents=True)
        # Write inclusion guard first.
        rhs_func_str = '#pragma once\n\n'
        # Write includes.
        rhs_func_str += '#include <math.h>\n\n'
        # Write macro definitions.
        rhs_func_str += '#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))\n'
        rhs_func_str += '#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))\n'
        rhs_func_str += '\n'
        # Write IVP name.
        rhs_func_str += '#define PROBLEM_NAME "{}"\n'.format(ivp.name)
        rhs_func_str += '#define PROBLEM_ID {}\n'.format(ivp.db_id)
        rhs_func_str += '\n'
        # Write ODE size.
        rhs_func_str += '#define n {}\n\n'.format(str(solve_equation('g', ivp.gridSize, 'n')[0]).replace('g', 'N'))
        # Write eval_range function.
        rhs_func_str += 'static inline void eval_range('
        rhs_func_str += 'int first, int last, double t, const double *%in, double *f) {\n'
        rhs_func_str += ivp.code_eval_range
        rhs_func_str += '}\n'
        # Write eval_component function.
        rhs_func_str += '\nstatic inline double eval_component('
        rhs_func_str += 'int j, double t, const double *%in) {\n'
        rhs_func_str += ivp.code_eval_component
        rhs_func_str += '}\n'
        # Write initial solution function.
        rhs_func_str += '\nstatic inline void initial_values('
        rhs_func_str += 'const double t, double *%in) {\n'
        rhs_func_str += ivp.code_initial_values
        rhs_func_str += '}\n'
        # Replace solution vector stubs.
        rhs_func_str = rhs_func_str.replace('%in', config.ode_solution_vector)
        # Replace IVP constants.
        constants: Optional[List[Tuple[str, Union[str, float, int]]]] = ivp.constants.as_tuple()
        for name, desc in ivp.constants.items():
            regex = r'(?![a-zA-Z0-9]-_)' + name + r'(?![a-zA-Z0-9-_])'
            rhs_func_str = sub(regex, str(eval_math_expr(desc.value, constants)), rhs_func_str)
        # Replace MAX, MIN.
        rhs_func_str = sub(r'[di]max', 'MAX', rhs_func_str)
        rhs_func_str = sub(r'[di]min', 'MIN', rhs_func_str)
        # Write RHS function to file.
        path = Path('{}/RHS_{}.h'.format(self.folder_ivp, ivp.name))
        with path.open('w') as file_handle:
            file_handle.write(rhs_func_str)
        # Format code with indent tool if available.
        format_codefile(path)

    def write_ode_method(self, method: ODEMethod):
        """Write ODE method to file.

        Parameters:
        -----------
        method: ODEMethod
            Used ODE method.

        Returns:
        --------
        -
        """
        # Create folder if it does not yet exist.
        if self.folder_method and not self.folder_method.exists():
            self.folder_method.mkdir(parents=True)
        # Write inclusion guard first.
        method_str = '#pragma once\n\n'
        # Write method name and database id.
        method_str += '#define METHOD_NAME "{}"\n'.format(method.name)
        method_str += '#define METHOD_ID {}\n'.format(method.db_id)
        method_str += '\n'
        # Write method characteristics.
        method_str += '#define s {}\n'.format(method.stages)
        method_str += '#define m {}\n'.format(method.correctorSteps)
        method_str += '#define ORDER {}\n'.format(method.order_)
        method_str += '\n'
        # Write butcher coefficients.
        method_str += self.write_butcher_table(method)
        # Write ODE method to file.
        path = Path('{}/ODE_{}.h'.format(self.folder_method, method.name))
        with path.open('w') as file_handle:
            file_handle.write(method_str)
        # Format code with indent tool if available.
        format_codefile(path)

    @staticmethod
    def write_butcher_table(method: ODEMethod) -> str:
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
        params_string = ImplCodeGeneratorCPP.write_2d_array('A', method.coefficientsA)
        params_string += ImplCodeGeneratorCPP.write_1d_array('b', method.coefficientsB)
        params_string += ImplCodeGeneratorCPP.write_1d_array('c', method.coefficientsC)
        return params_string + '\n'

    @staticmethod
    def write_2d_array(name: str, data: List[List[str]]) -> str:
        """Write 2D double array.

        Parameters:
        -----------
        name: str
            Name of the array.
        data: List of List of str
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
    def write_1d_array(name: str, data: List[str]) -> str:
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

    @staticmethod
    def write_function_description(skeleton: str) -> str:
        """Write function description.

        Parameters:
        -----------
        template: str
            Name of the associated impl skeleton.

        Returns:
        --------
        str
            Generated code.
        """
        func_args = ('const int me', 'const int first', 'const int last', 'double t', 'double h', 'double *y',
                     'DS_{}_t *ds'.format(skeleton))
        header = 'int step('
        header += ', '.join(func_args)
        header += '\n'
        header += '#ifdef STEP_CANCELLATION\n'
        header += ', double T_best_impl\n'
        header += '#endif\n'
        header += ') {\n'
        return header

    @staticmethod
    def write_instrument_impl(impl_variant_id: int):
        indent_lvl = 1
        instr = '#ifdef INSTRUMENT\n'
        instr += indent(indent_lvl) + 'if (me == 0)\n'
        instr += indent(indent_lvl) + '{\n'
        instr += indent(indent_lvl + 1) + r'printf("\n#ImplVariant-{}\n");'.format(impl_variant_id)
        instr += indent(indent_lvl) + '}\n'
        instr += '#endif\n'
        return instr
