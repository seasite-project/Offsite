"""@package impl_generator_cpp
Definition of class ImplCodeGeneratorCpp.
"""

from copy import deepcopy
from pathlib import Path
from re import sub
from subprocess import run, PIPE, CalledProcessError
from typing import Dict, List, Tuple

import attr
from sortedcontainers import SortedDict

import offsite.config
from offsite.codegen.code_tree import CodeTree, CodeNode, CodeNodeType, KernelNode, LoopNode
from offsite.codegen.codegen_util import eval_loop_boundary, format_codefile, write_closing_bracket, create_variant_name
from offsite.codegen.kernel_generator import KernelCodeGenerator
from offsite.descriptions.impl_skeleton import ImplSkeleton
from offsite.descriptions.ivp import IVP
from offsite.descriptions.kernel_template import KernelTemplate, Kernel
from offsite.descriptions.ode_method import ODEMethod
from offsite.descriptions.parser_utils import DatastructDict
from offsite.evaluation.math_utils import eval_math_expr, solve_equation, corrector_steps, stages


@attr.s
class ImplCodeGeneratorCpp:
    """Representation of the code generator for kerncraft kernel code.

    Attributes:
    -----------
    db_session : sqlalchemy.orm.session.Session
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
        TODO
    """
    db_session = attr.ib(type='sqlalchemy.orm.session.Session')
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
                 method: ODEMethod) -> Dict[str, str]:
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
        config = offsite.config.offsiteConfig

        # Set some members to match current code generation.
        code_tree = deepcopy(skeleton.code_tree).root
        code_tree.name = skeleton.name

        self.unroll_stack.clear()
        self.loaded_templates.clear()
        self.required_datastructs.clear()

        # Create folder if it does not yet exist.
        if self.folder_impl and not self.folder_impl.exists():
            self.folder_impl.mkdir(parents=True)

        # Generate implementation code tree.
        self.generate_implementation_tree(code_tree, skeleton, method)
        # Optimize tree: unroll
        self.unroll_tree(code_tree)
        # Generate code from tree and write to string.
        codes = dict()
        for vid, variant in impl_variants:
            # Generate required kernel codes.
            kernels = self.derive_impl_variant_kernels(variant)
            # Generate implementation variant.
            # .. generate code.
            name = Path('{}/{}.hpp'.format(self.folder_impl, create_variant_name(kernels, skeleton.name)))
            codes[name] = self.generate_impl_variant(code_tree, vid, kernels, skeleton, method, ivp)
            # .. generate tiled version too?
            if config.args.tile:
                name = Path('{}/{}_tiled.hpp'.format(self.folder_impl, create_variant_name(kernels, skeleton.name)))
                vid = 123456
                # TODO fix variant number
                codes[name] = self.generate_impl_variant(code_tree, vid, kernels, skeleton, method, ivp, True)
        return codes

    def generate_implementation_tree(self, node: CodeNode, skeleton: 'ImplSkeleton', method: 'ODEMethod'):
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
            constants = [corrector_steps(method), stages(method)]
            node.boundary = eval_loop_boundary(node.boundary, constants)
        elif node.type == CodeNodeType.KERNEL:
            node: KernelNode
            # Load kernel template from database.
            if node.template_name not in self.loaded_templates.keys():
                template = KernelTemplate.from_database(self.db_session, node.template_name)
                self.loaded_templates[node.template_name] = template
                self.required_datastructs = {**template.datastructs, **self.required_datastructs}
        # Traverse tree depth first.
        if node.child:
            self.generate_implementation_tree(node.child, skeleton, method)
        if node.next:
            self.generate_implementation_tree(node.next, skeleton, method)

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
            Used implementation skeleton.
        method: ODEMethod
            Used ODE method.
        ivp: IVP
            Used IVP.
        gen_tiled_code: bool
            Generated tiled code version.

        Returns:
        --------
        str
            Generated implementation variant code.
        """
        variant_name: str = create_variant_name(kernels, skeleton.name)
        # Generate implementation variant code.
        code: str = ImplCodeGeneratorCpp.generate_implementation_code(impl, kernels, method, gen_tiled_code)
        # Write frame code.
        return self.write_skeleton_includes(skeleton, method, ivp) + self.write_class(variant_name, variant_id,
                                                                                      method.name, code)

    @staticmethod
    def generate_implementation_code(
            node: CodeNode, kernels: List[Kernel], method: ODEMethod, gen_tiled_code: bool, code: str = '') -> str:
        """Write implementation variant code.

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
        config = offsite.config.offsiteConfig
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
            code += node.to_implementation_codeline(kernel_code, kernel.db_id, False)
        elif node.type != CodeNodeType.ROOT:
            code += node.to_implementation_codeline()
        # Traverse tree depth first.
        if node.child:
            code += ImplCodeGeneratorCpp.generate_implementation_code(node.child, kernels, method, gen_tiled_code)
        # Close opened brackets.
        if node.type == CodeNodeType.LOOP:
            code += write_closing_bracket(node.indent)
        if node.next:
            code += ImplCodeGeneratorCpp.generate_implementation_code(node.next, kernels, method, gen_tiled_code)
        return code

    def write_skeleton_includes(self, skeleton: ImplSkeleton, method: ODEMethod, ivp: IVP) -> str:
        """Write header includes.

        Parameters:
        -----------
        variant_id: int
            Database ID associated to the generated implementation variant.
        skeleton: ImplSkeleton
            Used implementation skeleton.
        method: ODEMethod
            Used ODE method.
        ivp: IVP
            Used IVP.

        Returns:
        --------
        str
            Generated code.
        """
        # Add inclusion guard first.
        includes = '#pragma once\n\n'
        # Define used datastructure.
        includes += '#define DATASTRUCTURE DS_Super\n\n'
        # Include base class.
        includes += '#include "../IRKSolver.hpp"\n'
        # Include datastructure, method and RHS classes.
        includes += '#include "../../ds/DS_Super.hpp"\n'
        includes += '#include "../../method/{}.hpp"\n'.format(method.name)
        includes += '#include "../../ODE/ODE.hpp"\n'
        # IVP.
        if skeleton.isIVPdependent:
            self.write_rhs_function(ivp)
        # ODE method.
        self.write_ode_method(method)
        # Data structures.
        # self.write_datastructures(skeleton.name)
        return includes + '\n'

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
        config = offsite.config.offsiteConfig
        # Create folder if it does not yet exist.
        if self.folder_ivp and not self.folder_ivp.exists():
            self.folder_ivp.mkdir(parents=True)
        # Determine IVP size from IVP grid size.
        n = str(solve_equation('g', ivp.gridSize, 'n')[0])
        n = n.replace('g**2', 'g * g')  # TODO fix
        # Write inclusion guard first.
        rhs_func_str = '#pragma once\n\n'
        # Write defines.
        if True:  # TODO condition
            rhs_func_str += '#define ODE_MPI_SPARSE\n\n'
        # Write includes.
        rhs_func_str += '#include <algorithm>\n'
        rhs_func_str += '#include <cmath>\n'
        rhs_func_str += '#include <iostream>\n'
        rhs_func_str += '#include "ODE.hpp"\n'
        rhs_func_str += '\n'
        # Write IVP name.
        rhs_func_str += '#define PROBLEM_NAME "{}"\n'.format(ivp.name)
        rhs_func_str += '#define PROBLEM_ID {}\n'.format(ivp.db_id)
        rhs_func_str += '\n'
        # Open class.
        rhs_func_str += 'class {} : public IVP\n'.format(ivp.name)
        rhs_func_str += '{\n'
        # Write members.
        rhs_func_str += 'public:\n\n'
        rhs_func_str += 'int N;\n'
        # Write constructor.
        rhs_func_str += '{}('.format(ivp.name)
        rhs_func_str += 'const int N_,'
        rhs_func_str += 'const double h_,'
        rhs_func_str += 'const double te_,'
        rhs_func_str += 'const double t_,'
        rhs_func_str += 'const double TOL_'
        rhs_func_str += ')\n'
        rhs_func_str += '{\n'
        rhs_func_str += 'N  = N_;\n'
        rhs_func_str += 'n = {};\n'.format(n)
        rhs_func_str += 'acc_dist = {};\n'.format(ivp.characteristic.access_distance)
        rhs_func_str += 'stencil_radius = {};\n'.format(
            -1 if not ivp.characteristic.isStencil else ivp.characteristic.stencil_radius)
        rhs_func_str += 'y = new double[n];\n'
        rhs_func_str += 'initial_value();\n'
        rhs_func_str += 'h = h_;\n'
        rhs_func_str += 'te = te_;\n'
        rhs_func_str += 't = t_;\n'
        rhs_func_str += 'TOL = TOL_;\n'
        # End constructor.
        rhs_func_str += '}\n\n'
        # Write destructor.
        rhs_func_str += '~{}()\n'.format(ivp.name)
        rhs_func_str += '{\n'
        rhs_func_str += 'delete[] y;\n'
        rhs_func_str += '}\n\n'
        # Write initial solution function.
        rhs_func_str += 'inline void initial_value()\n'
        rhs_func_str += '{\n'
        rhs_func_str += ivp.code_initial_values
        rhs_func_str += '}\n\n'
        # Write eval_range function.
        rhs_func_str += 'inline void eval_range('
        rhs_func_str += 'const int first, const int last, const double t, const double *%in, double *f) {\n'
        rhs_func_str += ivp.code_eval_range
        rhs_func_str += '}\n\n'
        # Write eval_component function.
        rhs_func_str += 'inline double eval_component('
        rhs_func_str += 'const int j, const double t, const double *%in) {\n'
        rhs_func_str += ivp.code_eval_component
        rhs_func_str += '}\n\n'
        # Write required_indices function.
        rhs_func_str += 'inline std::set<int> required_indices('
        rhs_func_str += 'const int first, const int last) {\n'
        rhs_func_str += ivp.code_required_indices
        rhs_func_str += '}\n\n'
        # End class
        rhs_func_str += '};\n'
        # Replace solution vector stubs.
        rhs_func_str = rhs_func_str.replace('%in', config.ode_solution_vector)
        # Replace IVP constants.
        constants = ivp.constants.as_tuple()
        for name, desc in ivp.constants.items():
            regex = r'(?![a-zA-Z0-9]-_)' + name + r'(?![a-zA-Z0-9-_])'
            rhs_func_str = sub(regex, str(eval_math_expr(desc.value, constants)), rhs_func_str)
        # Substitute some code parts to get the desired code format.
        # Replace 'g' with 'N'.
        # TODO Hack
        rhs_func_str = rhs_func_str.replace(' g', ' N')
        rhs_func_str = rhs_func_str.replace('g ', 'N ')
        rhs_func_str = rhs_func_str.replace('[g', '[N')
        rhs_func_str = rhs_func_str.replace('g]', 'N]')
        rhs_func_str = rhs_func_str.replace(' g ', ' N ')
        rhs_func_str = rhs_func_str.replace('/g', '/N')
        rhs_func_str = rhs_func_str.replace('*g', '*N')
        rhs_func_str = rhs_func_str.replace('imin', 'std::min')
        rhs_func_str = rhs_func_str.replace('imax', 'std::max')
        rhs_func_str = rhs_func_str.replace('dmin', 'std::min')
        rhs_func_str = rhs_func_str.replace('dmax', 'std::max')
        # Write RHS function to file.
        path = Path('{}/RHS_{}.hpp'.format(self.folder_ivp, ivp.name))
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
        # Write include of base class.
        method_str += '#include "method.hpp"\n'
        method_str += '\n'
        # Open struct.
        method_str += 'struct {} : Method\n'.format(method.name)
        method_str += '{\n'
        # Write constructor.
        method_str += '{}()\n'.format(method.name)
        method_str += '{\n'
        # Write method characteristics
        method_str += 's = {};\n'.format(method.stages)
        method_str += 'm = {};\n'.format(method.correctorSteps)
        method_str += 'ORD = {};\n'.format(method.order_)
        method_str += '\n'
        # Write butcher coefficients.
        method_str += self.write_butcher_table(method)
        # End constructor.
        method_str += '}\n'
        method_str += '\n'
        # Write destructor.
        method_str += '~{}()\n'.format(method.name)
        method_str += '{\n'
        method_str += 'for (int i=0; i<s; ++i)\n'
        method_str += 'delete[] A[i];\n'
        method_str += 'delete A;\n'
        method_str += 'delete b;\n'
        method_str += 'delete c;\n'
        method_str += '}\n'
        # End struct.
        method_str += '};'
        # Write to file.
        path = Path('{}/{}.hpp'.format(self.folder_method, method.name))
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
        butcher = ImplCodeGeneratorCpp.write_2d_array('A', method.coefficientsA)
        butcher += ImplCodeGeneratorCpp.write_1d_array('b', method.coefficientsB)
        butcher += ImplCodeGeneratorCpp.write_1d_array('c', method.coefficientsC)
        return butcher + '\n'

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
        array_str = name + ' = new double*[' + str(len(data)) + '];\n'
        for idx, row in enumerate(data):
            array_str += name + '[' + str(idx) + '] = new double[' + str(len(data)) + '] {'
            for elem in row:
                array_str += str(eval_math_expr(elem)) + ', '
            array_str = array_str[:-2]
            array_str += '};\n'
        array_str += '\n'
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
        array_str = name + ' = new double[' + str(len(data)) + '] {'
        for elem in data:
            array_str += str(eval_math_expr(elem)) + ', '
        array_str = array_str[:-2]
        array_str += '};\n'
        array_str += '\n'
        return array_str

    def write_class(self, variant_name: str, variant_id: int, method_name: str, step_code: str) -> str:
        """Write class.

        Parameters:
        -----------
        variant_id: int
            Database ID associated to the generated implementation variant.
        step_code: str
            Generated step function code.

        Returns:
        --------
        str
            Generated code.
        """
        class_str = 'class {} : public IRKSolver\n'.format(variant_name)
        class_str += '{\n'
        class_str += 'std::shared_ptr<DS_Super> ds;\n\n'
        class_str += 'public:\n\n'
        # Write constructor.
        class_str += '{}('.format(variant_name)
        class_str += 'std::shared_ptr<IVP> ode,'
        class_str += 'std::shared_ptr<DS_Super> ds,'
        class_str += 'std::shared_ptr<{}> method,'.format(method_name)
        class_str += 'std::shared_ptr<HBD_1D> part'
        class_str += ')\n'
        class_str += '{\n'
        class_str += 'this->ode = ode;\n'
        class_str += 'this->ds = ds;\n'
        class_str += 'this->method = method;\n'
        class_str += 'this->part = part;\n'
        class_str += '}\n\n'
        # Write getter for variant name.
        class_str += 'const std::string variantName() const\n'
        class_str += '{\n'
        class_str += 'return \"{}\";\n'.format(variant_name)
        class_str += '}\n\n'
        # Write getter for variant id.
        class_str += 'const int variantID() const\n'
        class_str += '{\n'
        class_str += 'return {};\n'.format(variant_id)
        class_str += '}\n\n'
        # Write step function.
        func_args = ('const int me', 'const int first', 'const int last', 'double t', 'double h')
        class_str += 'void step('
        class_str += ', '.join(func_args)
        class_str += ')\n'
        class_str += '{\n'
        class_str += 'boost::mpi::communicator world;\n\n'
        # class_str += 'int NP = ode->n / world.size();\n\n'
        class_str += step_code
        class_str += '}\n'
        # End class.
        class_str += '};\n'
        # Substitute some code parts to get the desired code format.
        class_str = sub(r'\bA\[\b', 'method->A[', class_str)
        class_str = sub(r'\bb\[\b', 'method->b[', class_str)
        class_str = sub(r'\bc\[\b', 'method->c[', class_str)
        class_str = sub(r'\beval_', 'ode->eval_', class_str)
        for name, desc in self.required_datastructs.items():
            # Skip names reserved for Butcher table entries.
            if name in ('A', 'b', 'c'):
                continue
            # Skip names reserved for pre-defined variables.
            if name in ('h', 't', 'g', 'y'):
                continue
            class_str = sub(r'\b{}\b'.format(name), 'ds->{}'.format(name), class_str)
        class_str = sub(r'\by\b', 'ode->y', class_str)
        return class_str + '\n'

    @staticmethod
    def write_codes_to_file(generated_files: Dict[Path, str]):
        """Write generated implementation variant codes to files.

        Parameters:
        -----------
        generated_files: dict (key=Path, value=str)
            Generated implementation variant codes and their associated file paths.

        Returns:
        --------
        -
        """
        for path, code in generated_files.items():
            with path.open('w') as file_handle:
                file_handle.write(code)
            # Replace newline characters in printf commands.
            # TODO: Is there a Python solution for this?
            try:
                cmd = ['sed', '-i', 's/§n/n/g', '{}'.format(path)]
                run(cmd, check=True, stdout=PIPE).stdout
            except CalledProcessError as error:
                raise RuntimeError('Unable to substitute newline character stubs in {}: {}'.format(path, error))
            # Replace tabulator characters in printf commands.
            try:
                cmd = ['sed', '-i', 's/§t/t/g', '{}'.format(path)]
                run(cmd, check=True, stdout=PIPE).stdout
            except CalledProcessError as error:
                raise RuntimeError('Unable to substitute tabulator character stubs in {}: {}'.format(path, error))
            # Format code with indent tool if available.
            format_codefile(path)
