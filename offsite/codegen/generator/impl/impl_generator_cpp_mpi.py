"""@package codegen.generator.impl.impl_generator_cpp_mpi
Definition of class ImplCodeGeneratorCppMPI.

@author: Johannes Seiferth
"""
from re import sub
from typing import List, Optional

import attr
from pathlib2 import Path

from offsite.codegen.code_dsl.code_node import CodeNode
from offsite.codegen.codegen_util import format_codefile, create_variant_name
from offsite.codegen.generator.impl.impl_generator import ImplCodeGenerator
from offsite.descriptions.impl.impl_skeleton import ImplSkeleton
from offsite.descriptions.impl.kernel_template import Kernel
from offsite.descriptions.ode import IVP, ODEMethod
from offsite.util.math_utils import eval_math_expr, solve_equation


@attr.s
class ImplCodeGeneratorCppMPI(ImplCodeGenerator):
    """Representation of the code generator for impl variant CPP code used by the MPI driver framework.

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
    folder_ds: Path
        Folder all created datastructure code files are stored in.
    unroll_stack: SortedDict
        Stack of loops to be unrolled.
    loaded_templates: dict (key=str, val=KernelTemplate)
        Required and already loaded kernel templates.
    required_datastructs: DatastructDict
        Required and used data structures.
    """
    folder_ds = attr.ib(type=Path)

    def _generate_impl_variant(self, impl: CodeNode, variant_id: int, kernels: List[Kernel], skeleton: ImplSkeleton,
                               method: Optional[ODEMethod] = None, ivp: Optional[IVP] = None,
                               gen_tiled_code: bool = False) -> str:
        """Write impl variant code.

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
        variant_name: str = create_variant_name(kernels, skeleton.name)
        # Generate impl variant code.
        code: str = self._generate_implementation_code(impl, kernels, method, gen_tiled_code)
        # Write frame code.
        return self._write_skeleton_includes(variant_id, skeleton, method, ivp) + self._write_class(
            variant_name, variant_id, method, ivp, code)

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
        # Add inclusion guard first.
        includes = '#pragma once\n\n'
        # Define used datastructure.
        includes += '#define DATASTRUCTURE DS_Super\n\n'
        # Include base class.
        includes += '#include "../IRKSolver.hpp"\n'
        # Include datastructure class.
        includes += '#include "../../ds/DS_Super.hpp"\n'
        # IVP.
        if skeleton.isIVPdependent:
            assert IVP is not None
            includes += '#include "../../ODE/ODE.hpp"\n'
            self._write_rhs_function(ivp)
        # ODE method.
        if method is not None:
            includes += '#include "../../method/{}.hpp"\n'.format(method.name)
            self._write_ode_method(method)
        return includes + '\n'

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
        # TODO

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
        # Create folder if it does not yet exist.
        if self.folder_ivp and not self.folder_ivp.exists():
            self.folder_ivp.mkdir(parents=True)
        # Determine IVP size from IVP grid size.
        n = str(solve_equation('g', ivp.gridSize, 'n')[0])
        n = n.replace('g**2', 'g * g')  # TODO fix
        # Write inclusion guard first.
        rhs_func_str = '#pragma once\n\n'
        # Write defines.
        if ivp.characteristic.isSparse:
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
        rhs_func_str += 'g  = N_;\n'
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
        input_vec = 'y'
        rhs_func_str = rhs_func_str.replace('%in', input_vec)
        # Replace IVP constants.
        constants = ivp.constants.as_tuple()
        for name, desc in ivp.constants.items():
            regex = r'(?![a-zA-Z0-9]-_)' + name + r'(?![a-zA-Z0-9-_])'
            rhs_func_str = sub(regex, str(eval_math_expr(desc.value, constants)), rhs_func_str)
        # Replace MAX, MIN.
        rhs_func_str = sub(r'[di]max', 'MAX', rhs_func_str)
        rhs_func_str = sub(r'[di]min', 'MIN', rhs_func_str)
        # Write RHS function to file.
        path = Path('{}/RHS_{}.hpp'.format(self.folder_ivp, ivp.name))
        with path.open('w') as file_handle:
            file_handle.write(rhs_func_str)
        # Format code with indent tool if available.
        format_codefile(path)

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
        method_str += self._write_butcher_table(method)
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
        array_str = name + ' = new double[' + str(len(data)) + '] {'
        for elem in data:
            array_str += str(eval_math_expr(elem)) + ', '
        array_str = array_str[:-2]
        array_str += '};\n'
        array_str += '\n'
        return array_str

    def _write_class(self, variant_name: str, variant_id: int, method: ODEMethod, ivp: IVP, step_code: str) -> str:
        """Write solver class.

        Parameters:
        -----------
        variant_id: int
            Database ID associated to the generated impl variant.
        step_code: str
            Generated step function code.

        Returns:
        --------
        str
            Generated code.
        """
        class_str = 'class {} : public Solver\n'.format(variant_name)
        class_str += '{\n'
        class_str += 'std::shared_ptr<DS_Super> ds;\n\n'
        class_str += 'public:\n\n'
        # Write constructor.
        class_str += '{}('.format(variant_name)
        if ivp is not None:
            class_str += 'std::shared_ptr<IVP> ode,'
        class_str += 'std::shared_ptr<DS_Super> ds,'
        if method is not None:
            class_str += 'std::shared_ptr<{}> method,'.format(method.name)
        class_str += 'std::shared_ptr<HBD_1D> part'
        class_str += ')\n'
        class_str += '{\n'
        if ivp is not None:
            class_str += 'this->ode = ode;\n'
        class_str += 'this->ds = ds;\n'
        if method is not None:
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
        class_str += step_code
        class_str += '}\n'
        # End class.
        class_str += '};\n'
        # Substitute some code parts to get the desired code format.
        if method is not None:
            class_str = sub(r'\bA\[\b', 'method->A[', class_str)
            class_str = sub(r'\bb\[\b', 'method->b[', class_str)
            class_str = sub(r'\bc\[\b', 'method->c[', class_str)
        if ivp is not None:
            class_str = sub(r'\beval_', 'ode->eval_', class_str)
        for name, desc in self.required_datastructs.items():
            # Skip names reserved for Butcher table entries.
            if name in ('A', 'b', 'c') and method is not None:
                continue
            # Skip names reserved for pre-defined variables.
            if name in ('h', 't', 'g', 'y') and ivp is not None:
                continue
            class_str = sub(r'\b{}\b'.format(name), 'ds->{}'.format(name), class_str)
        class_str = sub(r'\by\b', 'ode->y', class_str)
        return class_str + '\n'
