"""@package impl_generator
Definition of class ImplCodeGenerator.
"""

from copy import deepcopy
from pathlib import Path
from re import sub
from subprocess import run, PIPE, CalledProcessError
from typing import Dict, List, Tuple

import attr
from sortedcontainers import SortedDict

import offsite.config
from offsite.codegen.code_tree import CodeTree, CodeNode, CodeNodeType
from offsite.codegen.codegen_util import eval_loop_boundary, format_codefile, indent, write_closing_bracket
from offsite.codegen.kernel_generator import KernelCodeGenerator
from offsite.descriptions.impl_skeleton import ImplSkeleton
from offsite.descriptions.ivp import IVP
from offsite.descriptions.kernel_template import KernelTemplate, Kernel
from offsite.descriptions.ode_method import ODEMethod
from offsite.evaluation.math_utils import eval_math_expr


@attr.s
class ImplCodeGenerator:
    """Representation of the code generator for implementation variant C code.

    Attributes:
    -----------
    db_session:
        TODO
    folder:
        TODO
    unroll_stack: SortedDict
        Stack of loops to be unrolled.
    loaded_templates: dict (key=str, val=KernelTemplate)
        TODO

    """
    db_session = attr.ib(type='sqlalchemy.orm.session.Session')
    folder = attr.ib(type=Path, default=Path('tmp/variants/'))
    unroll_stack = attr.ib(type=SortedDict, default=SortedDict())
    loaded_templates = attr.ib(type=Dict, default=dict())

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
        # Set some members to match current code generation.
        code_tree = deepcopy(skeleton.code_tree).root
        code_tree.name = skeleton.name

        self.unroll_stack.clear()

        # Generate implementation code trees.
        self.generate_implementation_trees(code_tree, skeleton, method)
        # Optimize tree: unroll
        self.unroll_tree(code_tree)
        # Generate code from tree and write to string.
        codes = dict()
        for vid, variant in impl_variants:
            # Generate required kernel codes.
            kernels = self.derive_impl_variant_kernels(variant)
            # Generate implementation variant.
            name = Path('{}/{}.h'.format(self.folder, self.create_variantname(kernels, skeleton)))
            codes[name] = self.generate_impl_variant(code_tree, vid, kernels, skeleton, method, ivp)
        return codes

    @staticmethod
    def create_variantname(variant: List[Kernel], skeleton) -> str:
        """Create filename of an implementation variant.

        Parameters:
        -----------
        variant : list of Kernel
            Kernels associated with this implementation variant.

        Returns:
        --------
        str
            Created filename.
        """
        return '{}_{}'.format(skeleton.name, '_'.join((kernel.name for kernel in variant)))

    def generate_implementation_trees(self, node: CodeNode, skeleton: 'ImplSkeleton', method: 'ODEMethod'):
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
            # Evaluate loop boundary expression.
            node.boundary = eval_loop_boundary(method, node.boundary)
        elif node.type == CodeNodeType.KERNEL:
            # Load kernel template from database.
            if node.template_name not in self.loaded_templates.keys():
                template = KernelTemplate.from_database(self.db_session, node.template_name)
                self.loaded_templates[node.template_name] = template
        # Traverse tree depth first.
        if node.child:
            self.generate_implementation_trees(node.child, skeleton, method)
        if node.next:
            self.generate_implementation_trees(node.next, skeleton, method)

    def derive_impl_variant_kernels(self, impl_variant: List[str]) -> List[Kernel]:
        return [self.map_kernel_id_to_object(kid) for kid in impl_variant]

    def map_kernel_id_to_object(self, kernel_id: str) -> 'Kernel':
        # Query kernel object from database.
        kernel = Kernel.select(self.db_session, kernel_id)
        # Select corresponding kernel template.
        template = self.loaded_templates[kernel.template.name]
        try:
            return next(filter(lambda x: x.name == kernel.name, template.variants))
        except StopIteration:
            raise RuntimeError('')

    def generate_impl_variant(self, impl: CodeNode, variant_id: int, kernels: List[Kernel], skeleton: ImplSkeleton, method: ODEMethod, ivp: IVP) -> str:
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

        Returns:
        --------
        str
            Generated implementation variant code.
        """
        # Generate implementation variant code.
        code = ImplCodeGenerator.generate_implementation_code(impl, kernels, method)
        # Write frame code.
        code = self.write_skeleton_includes(variant_id, skeleton, method, ivp) + \
               ImplCodeGenerator.write_function_description() + ImplCodeGenerator.write_instrument_impl(1) + code + '}'
        return code

    @staticmethod
    def generate_implementation_code(
            node: CodeNode, kernels: Dict[str, 'Kernel'], method: ODEMethod, code: str = '') -> str:
        """Write implementation variant code.

        Parameters:
        -----------
        node: CodeNode
            Root node of the code tree.
        kernels: dict (key=str, val=Kernel)
            Used kernels.
        method: ODEMethod
            Used ODE method.
        code: str
            Current code string.

        Returns:
        --------
        str
            Generated yasksite code.
        """
        config = offsite.config.offsiteConfig
        if node.type == CodeNodeType.KERNEL:
            try:
                kernel = next(filter(lambda x: x.template.name == node.template_name, kernels))
            except StopIteration:
                print('')
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
            kernel_code = KernelCodeGenerator().generate(kernel, method, input_vector)
            code += node.to_implementation_codeline(kernel_code, kernel.db_id)
        elif node.type != CodeNodeType.ROOT:
            code += node.to_implementation_codeline()
        # Traverse tree depth first.
        if node.child:
            code += ImplCodeGenerator.generate_implementation_code(node.child, kernels, method)
        # Close opened brackets.
        if node.type == CodeNodeType.LOOP:
            code += write_closing_bracket(node.indent)
        if node.next:
            code += ImplCodeGenerator.generate_implementation_code(node.next, kernels, method)
        return code

    def write_skeleton_includes(self, variant_id: int, skeleton: ImplSkeleton, method: ODEMethod, ivp: IVP) -> str:
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
        codegen = skeleton.codegen
        includes = ''
        # Add inclusion guard first.
        includes += '#pragma once\n\n'
        # Add variant database ID as define.
        includes += '#define VARIANT_ID {}\n'.format(variant_id)
        includes += '\n'
        # IVP.
        if skeleton.isIVPdependent:
            includes += '#include <math.h>\n'
            includes += '#include "RHS_{}.h"\n'.format(ivp.name)
            self.write_rhs_function(ivp)
        # ODE method.
        includes += '#include "ODE_{}.h"\n'.format(method.name)
        self.write_ode_method(method)
        # Data structures.
        if 'datastructs' in codegen:
            includes += '#include "{}"\n'.format(codegen['datastructs'])
        # Instrumentation.
        includes += '#ifdef INSTRUMENT\n#include "timesnap.h"\n#endif\n'
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
        if self.folder and not self.folder.exists():
            self.folder.mkdir(parents=True)
        # Write inclusion guard first.
        rhs_func_str = '#pragma once\n\n'
        # Write IVP name.
        rhs_func_str += '#define PROBLEM_NAME "{}"\n'.format(ivp.name)
        rhs_func_str += '#define PROBLEM_ID {}\n'.format(ivp.db_id)
        rhs_func_str += '\n'
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
        constants = ivp.constants.as_tuple()
        for name, desc in ivp.constants.items():
            regex = r'(?![a-zA-Z0-9]-_)' + name + r'(?![a-zA-Z0-9-_])'
            rhs_func_str = sub(regex, str(eval_math_expr(desc.value, constants)), rhs_func_str)
        # Write RHS function to file.
        path = Path('{}/RHS_{}.h'.format(self.folder, ivp.name))
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
        if self.folder and not self.folder.exists():
            self.folder.mkdir(parents=True)
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
        path = Path('{}/ODE_{}.h'.format(self.folder, method.name))
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
        params_string = ImplCodeGenerator.write_2d_array('A', method.coefficientsA)
        params_string += ImplCodeGenerator.write_1d_array('b', method.coefficientsB)
        params_string += ImplCodeGenerator.write_1d_array('c', method.coefficientsC)
        return params_string + '\n'

    @staticmethod
    def write_2d_array(name: str, data: List[List[str]]) -> str:
        """Write 2D double array.

        Parameters:
        -----------
        name : str
            Name of the array.
        data : List of List of str
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
        name : str
            Name of the array.
        data : List of str
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
    def write_function_description() -> str:
        """Write function description.

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Generated code.
        """
        func_args = ('const int me', 'const int first', 'const int last', 'double t', 'double h')
        header = 'void timestep('
        header += ', '.join(func_args)
        return header + ') {\n'

    @staticmethod
    def write_instrument_impl(impl_variant_id: int):
        indent_lvl = 1
        instr = '#ifdef INSTRUMENT\n'
        instr += indent(indent_lvl) + 'if (me == 0)\n'
        instr += indent(indent_lvl) + '{\n'
        # TODO
        # instr += indent(indent_lvl + 1) + r'printf("\n#ImplVariant-{}\n");'.format('')
        instr += indent(indent_lvl) + '}\n'
        instr += '#endif\n'
        return instr

    @staticmethod
    def write_codes_to_file(generated_files: Dict[Path, str]):
        """Write generated implementation variant codes to files.

        Parameters:
        -----------
        generated_files : dict (key=Path, value=str)
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
