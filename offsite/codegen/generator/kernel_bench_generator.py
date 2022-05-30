"""@package codegen.generator.kernel_bench_generator
Definition of class KernelBenchCodeGenerator.

@author: Johannes Seiferth
"""

from copy import deepcopy
from re import sub
from typing import Dict, List, Optional, Tuple, Union

import attr
from pathlib2 import Path
from sortedcontainers import SortedDict

import offsite.config
from offsite.codegen.code_dsl.code_node import CodeNode, CodeNodeType, ComputationNode, LoopNode
from offsite.codegen.code_dsl.code_tree import CodeTree
from offsite.codegen.codegen_util import eval_loop_boundary, substitute_rhs_func_call, write_closing_bracket, \
    write_code_to_file
from offsite.config import Config
from offsite.descriptions.ode import IVP, ODEMethod, corrector_steps, stages
from offsite.descriptions.parser import DatastructType
from offsite.util.math_utils import eval_math_expr


@attr.s
class KernelBenchFiles:
    main_src = attr.ib(type=Path)
    kernel_src = attr.ib(type=Path)
    dummy_src = attr.ib(type=Path)


@attr.s
class KernelBenchCodeGenerator:
    """Representation of the code generator for kernel code.

    Attributes:
    -----------
    unroll_stack: SortedDict
        Stack of loops to be unrolled.
    split_stack: dict (key=str, val=list of CodeNode)
        Stack of loops to be split.
    """
    unroll_stack = attr.ib(type=SortedDict, default=SortedDict())
    split_stack = attr.ib(type=dict, default=dict())

    def collect_loop_meta_data(self, node: CodeNode, loop_splits: Optional[Dict[str, int]] = None):
        """Traverse tree depth first and collect all unroll and split tasks.

        Parameters:
        -----------
        node: CodeNode
            Root node of the code tree.
        loop_splits: dict (key=str, val=str)
            Requested loop splits (key=loop var; val=split at loop index)

        Returns:
        --------
        -
        """
        if node.type == CodeNodeType.LOOP:
            node: LoopNode
            # Collect unroll tasks.
            if node.optimize_unroll is not None:
                if node.optimize_unroll in self.unroll_stack:
                    self.unroll_stack[node.optimize_unroll].append(node)
                else:
                    self.unroll_stack[node.optimize_unroll] = [node]
            # Collect split tasks.
            if loop_splits:
                if node.var in loop_splits:
                    node.split_at = loop_splits[node.var]
                    if node.var in self.split_stack:
                        self.split_stack[node.var].append(node)
                    else:
                        self.split_stack[node.var] = [node]
        # Traverse tree depth first.
        if node.child:
            self.collect_loop_meta_data(node.child, loop_splits)
        if node.next:
            self.collect_loop_meta_data(node.next, loop_splits)

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
            CodeTree.unroll_loop(self.unroll_stack.values()[0][0])

    def generate(self, kernel: 'Kernel', method: ODEMethod, ivp: IVP) -> KernelBenchFiles:
        """Generate kernel bench code for name particular Kernel object.

        Parameters:
        -----------
        kernel: Kernel
            Kernel for which code is generated.
        method: ODEMethod
            ODE method for which code is generated.
        ivp: IVP
            IVP for which code is generated.

        Returns:
        --------
        KernelBenchFiles
            Paths to generated code files.
        """
        # Set some members to match current code generation.
        code_tree: CodeNode = deepcopy(kernel.code_tree).root

        self.unroll_stack.clear()
        self.split_stack = dict()

        # Determine the used input vector.
        try:
            input_vector = kernel.template.codegen['RHS_input']
        except KeyError:
            raise RuntimeError(
                'Kernel template {} is missing required entry \'RHS_input\'.'.format(kernel.template.name))

        # Generate kernel code tree.
        self.generate_kernel_tree(code_tree, kernel, method, input_vector)
        # Read loop split information from template.
        loop_splits: Optional[Dict[str, int]] = dict()
        if 'loop splits' in kernel.template.codegen:
            for ls in kernel.template.codegen['loop splits']:
                var, split_at = ls.split(' ')
                loop_splits[var] = int(split_at)
        # Collect loop data.
        self.collect_loop_meta_data(code_tree, loop_splits)
        # Optimize tree: split
        if self.split_stack:
            for loops in self.split_stack.values():
                for loop in loops:
                    CodeTree.split_loop(loop)
        # Collect loop data again to include possibly splitted loops.
        self.collect_loop_meta_data(code_tree, loop_splits)
        # Optimize tree: unroll
        self.unroll_tree(code_tree)
        # Substitute butcher table coefficients.
        CodeTree.substitute_butcher_coefficients(
            code_tree, method.coefficientsA, method.coefficientsB, method.coefficientsC)
        # Write header files.
        self.write_util_header()
        self.write_kernel_header(kernel)
        # Write source code files.
        path_main = self.write_main_src(kernel, ivp)
        path_kernel = self.write_kernel_src(code_tree, ivp, kernel)
        path_dummy = self.write_dummy_src()
        return KernelBenchFiles(path_main, path_kernel, path_dummy)

    @staticmethod
    def generate_kernel_tree(node: CodeNode, kernel: 'Kernel', method: ODEMethod, input_vector: str):
        """Generate yasksite code tree.

        Parameters:
        -----------
        node: CodeNode
            Root node of the code tree.
        kernel: Kernel
            Used kernel.
        method: ODEMethod
            Used ODE method.
        input_vector: str
            Name of the used input vector.

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
        elif node.type == CodeNodeType.COMPUTATION:
            node: ComputationNode
            template = kernel.template
            # Substitute computation.
            node.computation = template.computations[node.computation].computation
            # Check if computation includes RHS calls.
            if '%RHS' in node.computation:
                assert ('RHS' in template.codegen or 'RHS' in kernel.codegen)
                assert 'RHS_butcher_nodes' in template.codegen
                rhs_func: str = template.codegen['RHS'] if 'RHS' not in kernel.codegen else kernel.codegen['RHS']
                butcher_node: str = template.codegen['RHS_butcher_nodes']
                node.computation = substitute_rhs_func_call(node.computation, rhs_func, input_vector, butcher_node)
                if rhs_func == 'eval_range':
                    node.parent.flag = True
        elif node.type in [CodeNodeType.COMMUNICATION_CLUST, CodeNodeType.COMMUNICATION_NODE]:
            pass
        elif node.type == CodeNodeType.ROOT:
            node.name = kernel.name
        # Traverse tree depth first.
        if node.child:
            KernelBenchCodeGenerator.generate_kernel_tree(node.child, kernel, method, input_vector)
        if node.next:
            KernelBenchCodeGenerator.generate_kernel_tree(node.next, kernel, method, input_vector)

    @staticmethod
    def generate_kernel_code(node: CodeNode, code: str = ''):
        """Write kernel code.

        Parameters:
        -----------
        node: CodeNode
            Root node of the code tree.
        rhs_func: str
            Name of the used RHS function.
        kernel_name: str
            Name of the generated kernel.
        code: str
            Current code string.

        Returns:
        --------
        str
            Generated kernel code.
        """
        config: Config = offsite.config.offsiteConfig
        loop_skipped = False
        if node.type == CodeNodeType.LOOP and node.var == config.var_idx:
            node: LoopNode
            if node.flag:
                # Loop is replaced by RHS eval_range call and can thus be skipped.
                loop_skipped = True
            else:
                code += node.to_kernel_codeline()
        elif node.type != CodeNodeType.ROOT and node.type != CodeNodeType.PMODEL:
            code += node.to_kernel_codeline()
        # Traverse tree depth first.
        if node.child:
            code += KernelBenchCodeGenerator.generate_kernel_code(node.child)
        if node.type == CodeNodeType.LOOP:
            if not loop_skipped:
                code += write_closing_bracket(node.indent)
        if node.next:
            code += KernelBenchCodeGenerator.generate_kernel_code(node.next)
        return code

    @staticmethod
    def write_rhs_function(ivp: IVP):
        """Write RHS functions to str.

        Parameters:
        -----------
        ivp: IVP
            Used IVP.

        Returns:
        --------
        -
        """
        # Write eval_range function.
        rhs_func_str = 'inline void eval_range('
        rhs_func_str += 'int first, int last, double t, const double *%in, double *f) {\n'
        rhs_func_str += ivp.code_eval_range
        rhs_func_str += '}\n'
        # Write eval_component function.
        rhs_func_str += '\ninline double eval_component('
        rhs_func_str += 'int j, double t, const double *%in) {\n'
        rhs_func_str += ivp.code_eval_component
        rhs_func_str += '}\n'
        # Replace solution vector stubs.
        input_vec = 'y'
        rhs_func_str = rhs_func_str.replace('%in', input_vec)
        # Replace IVP constants.
        constants: Optional[List[Tuple[str, Union[str, float, int]]]] = ivp.constants.as_tuple()
        for name, desc in ivp.constants.items():
            regex = r'(?![a-zA-Z0-9]-_)' + name + r'(?![a-zA-Z0-9-_])'
            rhs_func_str = sub(regex, str(eval_math_expr(desc.value, constants)), rhs_func_str)
        return rhs_func_str

    @staticmethod
    def write_kernel_src(code_tree, ivp: IVP, kernel: 'Kernel') -> Path:
        code = '#include "kernel_{}.h"\n'.format(kernel.name)
        code += '\n'
        code += '#include <math.h>\n'
        code += '\n'
        # Write RHS functions.
        code += KernelBenchCodeGenerator.write_rhs_function(ivp)
        code += '\n'
        # Write kernel function.
        code_args = 'int first, int last'
        for name, desc in kernel.template.datastructs.items():
            if desc.struct_type == DatastructType.array1D:
                code_args += ', {} *{}'.format(desc.datatype, name)
            elif desc.struct_type == DatastructType.array2D:
                code_args += ', {} **{}'.format(desc.datatype, name)
            elif desc.struct_type == DatastructType.array3D:
                code_args += ', {} ***{}'.format(desc.datatype, name)
        kernel_code = 'inline void kernel({})\n'.format(code_args)
        kernel_code += '{\n'
        kernel_code += KernelBenchCodeGenerator.generate_kernel_code(code_tree)
        kernel_code += '}\n'
        code += kernel_code
        # Write to file.
        return write_code_to_file(code, 'tmp/kernel_{}'.format(kernel.name))

    @staticmethod
    def write_kernel_header(kernel: 'Kernel'):
        """Write kernel header to file.

         Parameters:
         -----------
         kernel: str
             Name of the used kernel.

         Returns:
         --------
         -
         """
        args = 'int first, int last'
        for name, desc in kernel.template.datastructs.items():
            if desc.struct_type == DatastructType.array1D:
                args += ', {} *{}'.format(desc.datatype, name)
            elif desc.struct_type == DatastructType.array2D:
                args += ', {} **{}'.format(desc.datatype, name)
            elif desc.struct_type == DatastructType.array3D:
                args += ', {} **{}'.format(desc.datatype, name)
        code = '#ifndef KERNEL_H_\n'
        code += '#define KERNEL_H_\n'
        code += '\n'
        code += 'extern int n;\n'
        code += 'extern int g;\n'
        code += '\n'
        code += '#include <math.h>\n'
        code += '\n'
        code += 'inline double dmin(double a, double b)\n'
        code += '{\n'
        code += 'if (a < b) return a;\n'
        code += 'else return b;\n'
        code += '}\n'
        code += '\n'
        code += 'inline double dmax(double a, double b)\n'
        code += '{\n'
        code += 'if (a < b) return b;\n'
        code += 'else return a;\n'
        code += '}\n'
        code += '\n'
        code += 'inline double imin(int a, int b)\n'
        code += '{\n'
        code += 'if (a < b) return a;\n'
        code += 'else return b;\n'
        code += '}\n'
        code += '\n'
        code += 'inline int imax(int a, int b)\n'
        code += '{\n'
        code += 'if (a < b) return b;\n'
        code += 'else return a;\n'
        code += '}\n'
        code += '\n'
        code += 'void kernel({});\n'.format(args)
        code += '\n'
        for name, desc in kernel.template.datastructs.items():
            if desc.struct_type == DatastructType.scalar:
                code += '{} {};\n'.format(desc.datatype, name)
        code += '\n'
        code += '#endif\n'
        # Write to file.
        write_code_to_file(code, 'tmp/kernel_{}'.format(kernel.name), suffix='.h')

    @staticmethod
    def write_main_src(kernel: 'Kernel', ivp: IVP) -> Path:
        # Write header includes.
        code = '#include <omp.h>\n'
        code += '#include <stdio.h>\n'
        code += '#include "kernel_{}.h"\n'.format(kernel.name)
        code += '#include "offsite_util.h"\n'
        code += '\n'
        # Write extern definitions.
        for name, desc in kernel.template.datastructs.items():
            if desc.struct_type == DatastructType.scalar:
                code += 'extern double {};\n'.format(name)
        code += 'void dummy(void *);\n'
        code += '\n'
        code += 'extern int var_false;\n'
        code += 'int n = 0;\n'
        code += 'int g = 0;\n'
        code += '\n'

        # Add main function code.
        code += 'int main(int argc, char **argv) {\n'
        code += 'const int threads = atoi(argv[2]);\n'
        code += 'n = atoi(argv[3]);\n'
        code += 'g = {};\n'.format(ivp.gridSize)
        code += 'const int s = atoi(argv[4]);\n'
        code += '\n'
        code += 'double times[atoi(argv[1])];\n'
        code += '\n'
        # Allocate and init (scalar) datastructures.
        for name, desc in kernel.template.datastructs.items():
            if desc.struct_type == DatastructType.array1D:
                code += '{} *{} = aligned_malloc1d({}, 32);\n'.format(desc.datatype, name, desc.size[0])
            elif desc.struct_type == DatastructType.array2D:
                code += '{} **{} = aligned_malloc2d({}, {}, 32);\n'.format(
                    desc.datatype, name, desc.size[0], desc.size[1])
            elif desc.struct_type == DatastructType.array3D:
                code += '{} ***{} = aligned_malloc3d({}, {}, {}, 32);\n'.format(
                    desc.datatype, name, desc.size[0], desc.size[1], desc.size[2])
        code += '\n'
        for name, desc in kernel.template.datastructs.items():
            if desc.struct_type == DatastructType.scalar:
                code += '{} = 0.048910108990;\n'.format(name)
        code += '\n'
        # Init time measurement.
        code += 'time_snap_t ts;\n'
        code += 'init_time_snap();\n'
        code += '\n'
        # Start parallel region.
        code += '#pragma omp parallel\n'
        code += '{\n'
        for name, desc in kernel.template.datastructs.items():
            if desc.struct_type == DatastructType.array1D:
                code += '#pragma omp barrier\n'
                code += '#pragma omp for\n'
                code += '#pragma nounroll_and_jam\n'
                code += '#pragma ivdep\n'
                code += 'for (int j = 0; j < {}; ++j)\n'.format(desc.size[0])
                code += '{}[j] = 0.048910108990;\n'.format(name)
            elif desc.struct_type == DatastructType.array2D:
                code += 'for (int i = 0; i < {}; ++i)\n'.format(desc.size[0])
                code += '{\n'
                code += '#pragma omp barrier\n'
                code += '#pragma omp for\n'
                code += '#pragma nounroll_and_jam\n'
                code += '#pragma ivdep\n'
                code += 'for (int j = 0; j < {}; ++j)\n'.format(desc.size[1])
                code += '{}[i][j] = 0.048910108990;\n'.format(name)
                code += '}\n'
            elif desc.struct_type == DatastructType.array3D:
                code += 'for (int i = 0; i < {}; ++i)\n'.format(desc.size[0])
                code += '{\n'
                code += 'for (int j = 0; j < {}; ++j)\n'.format(desc.size[1])
                code += '{\n'
                code += '#pragma omp barrier\n'
                code += '#pragma omp for\n'
                code += '#pragma nounroll_and_jam\n'
                code += '#pragma ivdep\n'
                code += 'for (int k = 0; k < {}; ++k)\n'.format(desc.size[2])
                code += '{}[i][j][k] = 0.048910108990;\n'.format(name)
                code += '}\n'
                assert False
        code += '\n'
        code_dummy = 'if (var_false)\n'
        code_dummy += '{\n'
        for name, desc in kernel.template.datastructs.items():
            if desc.struct_type == DatastructType.scalar:
                code_dummy += 'dummy(&{});\n'.format(name)
            elif desc.struct_type == DatastructType.array1D:
                code_dummy += 'dummy({});\n'.format(name)
            elif desc.struct_type == DatastructType.array2D:
                code_dummy += 'dummy({}[0]);\n'.format(name)
            elif desc.struct_type == DatastructType.array3D:
                code_dummy += 'dummy({}[0][0]);\n'.format(name)
        code_dummy += '}\n'
        code_dummy += '\n'
        code += code_dummy
        # Data distribution.
        code += 'const int me = omp_get_thread_num();\n'
        code += 'const int first = me * n / threads;\n'
        code += 'const int last = (me + 1) * n / threads - 1;\n'
        code += '\n'
        code += '#pragma omp barrier\n'
        code += '\n'
        #
        code += 'for (int warmup = 1; warmup >= 0; --warmup)\n'
        code += '{\n'
        code += 'int repeat = 2;\n'
        code += 'if (warmup == 0)\n'
        code += '{\n'
        code += 'repeat = atoi(argv[1]);\n'
        # code += 'likwid_markerStartRegion("loop");\n'
        code += '}\n'
        #
        code += '\n'
        code += 'for (; repeat > 0; --repeat)\n'
        code += '{\n'
        code += '#pragma omp barrier\n'
        code += 'if (me == 0)\n'
        code += 'time_snap_start(&ts);\n'
        code_args = 'first, last'
        for name, desc in kernel.template.datastructs.items():
            if desc.struct_type in [DatastructType.array1D, DatastructType.array2D, DatastructType.array3D]:
                code_args += ', {}'.format(name)
        code += 'kernel({});\n'.format(code_args)
        code += '\n'
        code += '#pragma omp barrier\n'
        code += 'if (me == 0)\n'
        code += '{\n'
        code += 'times[repeat-1] = time_snap_stop(&ts) / 1e9;\n'
        code += '}\n'
        code += code_dummy
        code += '}\n'
        code += '}\n'
        code += '\n'
        code += '}\n'
        code += '\n'
        # Print results.
        code += '\n'
        code += 'double total = 0.0;\n'
        code += 'for (int i = 0; i < atoi(argv[1]); ++i)\n'
        code += '{\n'
        # time per component.
        code += 'printf("%.15f\§n", times[i]);\n'
        code += 'total += times[i] / (double) n;\n'
        code += '}\n'
        code += '\n'
        # Deallocate datastructures.
        for name, desc in kernel.template.datastructs.items():
            if desc.struct_type == DatastructType.array1D:
                code += 'free1d({});\n'.format(name)
            elif desc.struct_type == DatastructType.array2D:
                code += 'free2d({});\n'.format(name)
            elif desc.struct_type == DatastructType.array3D:
                code += 'free3d({});\n'.format(name)
        code += '\n'
        code += 'return 0;\n'
        code += '}\n'
        code += '\n'
        # Write to file.
        return write_code_to_file(code, 'tmp/bench_{}'.format(kernel.name))

    @staticmethod
    def write_util_header():
        code = '#ifndef OFFSITE_UTIL_H_\n'
        code += '#define OFFSITE_UTIL_H_\n'
        code += '\n'
        code += '#include <assert.h>\n'
        code += '#include <stddef.h>\n'
        code += '#include <stdint.h>\n'
        code += '#include <stdlib.h>\n'
        code += '#include <sys/time.h>\n'
        code += '\n'
        code += '#define ALIGN(X, Y) ((unsigned long) ((((X) + (Y) - 1)/(Y)) * (Y)) % 256 < (Y) ? ((((X) + (Y) - 1)/(Y)) * (Y)) + (Y) : ((((X) + (Y) - 1)/(Y)) * (Y)))\n'
        code += '\n'
        code += 'inline void *aligned_malloc(size_t size, size_t align)\n'
        code += '{\n'
        code += 'void *result;\n'
        code += '#if defined(__INTEL_COMPILER)\n'
        code += 'result = _mm_malloc(size, align);\n'
        code += '#else\n'
        code += 'if (posix_memalign(&result, align, size)) result = 0;\n'
        code += '#endif\n'
        code += 'return result;\n'
        code += '}\n'
        code += '\n'
        code += 'double *aligned_malloc1d(size_t a, const size_t alignment)\n'
        code += '{\n'
        code += 'return (double *) aligned_alloc(alignment, a * sizeof(double));\n'
        code += '}\n'
        code += '\n'
        code += 'double **aligned_malloc2d(size_t a, size_t b, const size_t alignment)\n'
        code += '{\n'
        code += 'size_t i, row_size, row_count;\n'
        code += 'double ** x;\n'
        code += 'row_size = ALIGN(b * sizeof(double), alignment);\n'
        code += 'row_count = row_size / sizeof(double);\n'
        code += 'x = (double **)\n'
        code += 'aligned_alloc(alignment, a * sizeof(double *)); \n'
        code += 'x[0] = (double *)\n'
        code += 'aligned_alloc(alignment, a * row_size);\n'
        code += 'for (i = 1; i < a; i++)\n'
        code += 'x[i] = x[0] + i * row_count;\n'
        code += 'return x;\n'
        code += '}\n'
        code += '\n'
        code += 'double ***aligned_malloc3d(size_t a, size_t b, size_t c, const size_t alignment)\n'
        code += '{\n'
        code += 'size_t i, row_size, row_count;\n'
        code += 'double ***x;\n'
        code += 'assert ((alignment % sizeof(double)) == 0);\n'
        code += 'row_size = ALIGN(c * sizeof(double), alignment);\n'
        code += 'row_count = row_size / sizeof(double);\n'
        code += 'x = (double ** *)\n'
        code += 'aligned_alloc(alignment, a * sizeof(double **));\n'
        code += 'x[0] = (double **)\n'
        code += 'aligned_alloc(alignment, a * b * sizeof(double *));\n'
        code += 'x[0][0] = (double *)\n'
        code += 'aligned_alloc(alignment, a * b * row_size);\n'
        code += 'for (i = 1; i < a; i++)\n'
        code += 'x[i] = x[0] + i * b;\n'
        code += 'for (i = 1; i < a * b; i++)\n'
        code += 'x[0][i] = x[0][0] + i * row_count;\n'
        code += 'return x;\n'
        code += '}\n'
        code += 'static void free1d(double * p)\n'
        code += '{\n'
        code += 'free((void *) p);\n'
        code += '}\n'
        code += '\n'
        code += 'static void free2d(double ** p)\n'
        code += '{\n'
        code += 'free((void *) p[0]);\n'
        code += 'free((void *) p);\n'
        code += '}\n'
        code += '\n'
        code += 'static void free3d(double ** * p)\n'
        code += '{\n'
        code += 'free((void *) p[0][0]);\n'
        code += 'free((void *) p[0]);\n'
        code += 'free((void *) p);\n'
        code += '}\n'
        code += '\n'
        code += 'typedef struct timeval time_snap_t;\n'
        code += '\n'
        code += '#define init_time_snap()\n'
        code += '\n'
        code += 'inline void time_snap_start(time_snap_t * ts)\n'
        code += '{\n'
        code += 'gettimeofday(ts, NULL);\n'
        code += '}\n'
        code += '\n'
        code += 'inline uint64_t time_snap_stop(const time_snap_t * ts1)\n'
        code += '{\n'
        code += 'time_snap_t ts2;\n'
        code += 'gettimeofday( & ts2, NULL);\n'
        code += '\n'
        code += 'return ((uint64_t)(ts2.tv_sec - ts1->tv_sec) * 1000000 + (uint64_t) ts2.tv_usec - (uint64_t) ts1->tv_usec) * 1000;\n'
        code += '}\n'
        code += '\n'
        code += '#endif\n'
        # Write to file.
        write_code_to_file(code, 'tmp/offsite_util', suffix='.h')

    @staticmethod
    def write_dummy_src() -> Path:
        code = 'void dummy(void* a) {}\n'
        code += 'int var_false = 0;\n'
        # Write to file.
        return write_code_to_file(code, 'tmp/dummy')
