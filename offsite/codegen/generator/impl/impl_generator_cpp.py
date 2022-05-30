"""@package codegen.generator.impl.impl_generator_cpp
Definition of class ImplCodeGeneratorCPP.

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
from offsite.util.math_utils import eval_math_expr, solve_equation


@attr.s
class ImplCodeGeneratorCPP(ImplCodeGenerator):
    """Representation of the code generator for implementation variant CPP code.

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
        #
        code = sub(r'\bimin\b', 'MIN', code)
        code = sub(r'\bimax\b', 'MAX', code)
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
        frame = self._write_skeleton_includes(variant_id, skeleton, method, ivp)
        frame += ImplCodeGeneratorCPP._write_function_description(skeleton.name, gen_tiled_code)
        if local_var_defs != '':
            frame += local_var_defs + '\n'
        frame += ImplCodeGeneratorCPP._write_instrument_impl(variant_id)
        frame += code
        frame += 'return STEP_SUCCESS;\n}\n'  # Concludes code function block.
        return frame

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
        # Add variant database ID as define.
        includes += '#define VARIANT_ID {}\n'.format(variant_id)
        includes += '\n'
        # IVP.
        if skeleton.isIVPdependent:
            assert IVP is not None
            includes += '#include <math.h>\n'
            includes += '#include <{}/RHS_{}.h>\n'.format(self.folder_ivp, ivp.name)
            self._write_rhs_function(ivp)
        # ODE method.
        if method is not None:
            includes += '#include <{}/ODE_{}.h>\n'.format(self.folder_method, method.name)
            self._write_ode_method(method)
        # Data structures.
        self._write_datastructures(skeleton.name)
        includes += '#include <{}/DS_{}.h>\n'.format(self.folder_ds, skeleton.name)
        # Instrumentation.
        includes += '#ifdef INSTRUMENT\n#include "timesnap.h"\n#endif\n'
        # Step cancellation.
        includes += '#ifdef STEP_CANCELLATION\n'
        includes += '#include <step_cancellation.hpp>\n'
        includes += '#endif\n'
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
        rhs_func_str = rhs_func_str.replace('%in', 'y')
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
    def _write_function_description(skeleton: str, gen_tiled_code: bool) -> str:
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
        if gen_tiled_code:
            func_args = ('const int me', 'const int first', 'const int last', 'double t', 'double h', 'int B',
                         'double *y', 'DS_{}_t *ds'.format(skeleton))
        else:
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
    def _write_instrument_impl(impl_variant_id: int):
        indent_lvl = 1
        instr = '#ifdef INSTRUMENT\n'
        instr += indent(indent_lvl) + 'if (me == 0)\n'
        instr += indent(indent_lvl) + '{\n'
        instr += indent(indent_lvl + 1) + r'printf("\n#ImplVariant-{}\n");'.format(impl_variant_id)
        instr += indent(indent_lvl) + '}\n'
        instr += '#endif\n'
        return instr
