"""@package codegen.generator.impl.impl_generator_c
Definition of class ImplCodeGeneratorC.

@author: Johannes Seiferth
"""

from re import sub
from typing import List, Optional, Tuple, Union

import attr
from pathlib2 import Path

from offsite.codegen.code_dsl.code_node import CodeNode
from offsite.codegen.codegen_util import format_codefile, indent
from offsite.codegen.generator.impl.impl_generator import ImplCodeGenerator
from offsite.descriptions.impl.impl_skeleton import ImplSkeleton
from offsite.descriptions.impl.kernel_template import Kernel
from offsite.descriptions.ode import IVP, ODEMethod
from offsite.descriptions.parser import DatastructType
from offsite.util.math_utils import eval_math_expr


@attr.s
class ImplCodeGeneratorC(ImplCodeGenerator):
    """Representation of the code generator for implementation variant C code.

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
        # Generate impl variant code.
        code = self._generate_implementation_code(impl, kernels, method, gen_tiled_code)
        # Write frame code.
        return self._write_skeleton_includes(variant_id, skeleton, method, ivp) + \
               self._write_function_description() + self._write_instrument_impl(variant_id) + \
               code + '}'

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
        includes = ''
        # Add inclusion guard first.
        includes += '#pragma once\n\n'
        # Add variant database ID as define.
        includes += '#define VARIANT_ID {}\n'.format(variant_id)
        includes += '\n'
        # IVP.
        if skeleton.isIVPdependent:
            assert ivp is not None
            includes += '#include <math.h>\n'
            includes += '#include "RHS_{}.h"\n'.format(ivp.name)
            self._write_rhs_function(ivp)
        # ODE method.
        if method is not None:
            includes += '#include "ODE_{}.h"\n'.format(method.name)
            self._write_ode_method(method)
        # Data structures.
        self._write_datastructures(skeleton.name)
        includes += '#include "DS_{}.h"\n'.format(skeleton.name)
        # Instrumentation.
        includes += '#ifdef INSTRUMENT\n#include "timesnap.h"\n#endif\n'
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
        # Create folder if it does not yet exist.
        if self.folder_impl and not self.folder_impl.exists():
            self.folder_impl.mkdir(parents=True)
        # Write inclusion guard first.
        dstruct_str = '#pragma once\n\n'
        # Write data structures.
        ds_variables = ''
        ds_alloc_str = ''
        ds_free_str = ''
        for name, desc in self.required_datastructs.items():
            # Skip names reserved for Butcher table entries.
            if name in ('A', 'b', 'c'):
                continue
            # Skip names reserved for pre-defined variables.
            if name in ('h', 't', 'g'):
                continue
            # Add to variables.
            ds_variables += '{} {}{};\n'.format(desc.datatype, '*' * desc.struct_type, name)
            # Add to alloc and free.
            if desc.struct_type == DatastructType.scalar:
                continue
            if desc.struct_type == DatastructType.array1D:
                ds_alloc_str += '{} = alloc1d({});\n'.format(name, desc.size[0])
                ds_free_str += 'free1d({});\n'.format(name)
            elif desc.struct_type == DatastructType.array2D:
                ds_alloc_str += '{} = alloc2d({}, {});\n'.format(name, desc.size[0], desc.size[1])
                ds_free_str += 'free2d({});\n'.format(name)
            elif desc.struct_type == DatastructType.array3D:
                ds_alloc_str += '{} = alloc3d({}, {}, {});\n'.format(name, desc.size[0], desc.size[1], desc.size[2])
                ds_free_str += 'free3d({});\n'.format(name)
            else:
                assert False
        # Write variables.
        dstruct_str += ds_variables + '\n'
        # Write allocate function.
        dstruct_str += 'static void allocate_data_structures() {\n'
        dstruct_str += ds_alloc_str
        dstruct_str += '}\n\n'
        # Write free function.
        dstruct_str += 'static void free_data_structures() {\n'
        dstruct_str += ds_free_str
        dstruct_str += '}\n\n'
        # Write to file.
        path = Path('{}/DS_{}.h'.format(self.folder_impl, skeleton))
        with path.open('w') as file_handle:
            file_handle.write(dstruct_str)
        # Format code with indent tool if available.
        format_codefile(path)

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
        # Write inclusion guard first.
        rhs_func_str = '#pragma once\n\n'
        # Write IVP name.
        rhs_func_str += '#define PROBLEM_NAME "{}"\n'.format(ivp.name)
        rhs_func_str += '#define PROBLEM_ID {}\n'.format(ivp.db_id)
        rhs_func_str += '\n'
        # Write eval_range function.
        rhs_func_str += 'static inline void eval_range('
        rhs_func_str += 'int first, int last, double t, const double *{}, double *f) {\n'
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
        rhs_func_str = rhs_func_str.replace('%in', 'y')
        # Replace IVP constants.
        constants: Optional[List[Tuple[str, Union[str, float, int]]]] = ivp.constants.as_tuple()
        for name, desc in ivp.constants.items():
            regex = r'(?![a-zA-Z0-9]-_)' + name + r'(?![a-zA-Z0-9-_])'
            rhs_func_str = sub(regex, str(eval_math_expr(desc.value, constants)), rhs_func_str)
        # Write RHS function to file.
        path = Path('{}/RHS_{}.h'.format(self.folder_ivp, ivp.name))
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
        # Write method characteristics.
        method_str += '#define s {}\n'.format(method.stages)
        method_str += '#define m {}\n'.format(method.correctorSteps)
        method_str += '#define ORDER {}\n'.format(method.order_)
        method_str += '\n'
        # Write butcher coefficients.
        method_str += self._write_butcher_table(method)
        # Write ODE method to file.
        path = Path('{}/ODE_{}.h'.format(self.folder_method, method.name))
        with path.open('w') as file_handle:
            file_handle.write(method_str)
        # Format code with indent tool if available.
        format_codefile(path)

    @staticmethod
    def _write_function_description() -> str:
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
    def _write_instrument_impl(impl_variant_id: int):
        indent_lvl = 1
        instr = '#ifdef INSTRUMENT\n'
        instr += indent(indent_lvl) + 'if (me == 0)\n'
        instr += indent(indent_lvl) + '{\n'
        instr += indent(indent_lvl + 1) + r'printf("\n#ImplVariant-{}\n");'.format(impl_variant_id)
        instr += indent(indent_lvl) + '}\n'
        instr += '#endif\n'
        return instr
