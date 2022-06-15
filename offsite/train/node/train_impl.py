"""@package train.node.train_impl
Functions to train the tuning database with implementation variant predictions.

@author: Johannes Seiferth
"""

from copy import deepcopy
from datetime import datetime
from getpass import getuser
from multiprocessing import cpu_count, Pool
from statistics import mean
from traceback import print_exc
from typing import Dict, List, Optional, Tuple, Union

from pandas import DataFrame, concat
from pathlib2 import Path
from sqlalchemy.orm import Session

import offsite.config
from offsite import __version__
from offsite.codegen.codegen_util import write_code_to_file, write_codes_to_file
from offsite.codegen.generator.impl.impl_generator_c import ImplCodeGeneratorC
from offsite.config import Config, ModelToolType
from offsite.database import commit
from offsite.descriptions.impl.impl_skeleton import ImplSkeleton
from offsite.descriptions.impl.kernel_template import Kernel
from offsite.descriptions.machine import MachineState
from offsite.descriptions.ode import IVP, ODEMethod, ivp_system_size, ivp_grid_size, corrector_steps, stages
from offsite.train.communication.openmp.omp_barrier import OmpBarrierBenchmark
from offsite.train.impl_variant import ImplVariant, ImplVariantRecord, fuse_equal_impl_records
from offsite.train.node.util.communication_costs import compute_node_lvl_communication_costs
from offsite.train.node.util.performance_model import compute_impl_variant_pred
from offsite.train.train_utils import deduce_available_impl_variants, deduce_impl_variant_sample_intervals, \
    fetch_kernel_runtime_prediction_data
from offsite.util.math_utils import eval_math_expr, remove_outliers
from offsite.util.process_utils import run_process
from offsite.util.time import start_timer, stop_timer


def train_impl_variant_predictions(
        db_session: Session, machine: MachineState, skeletons: List[ImplSkeleton], predicted_kernels: List[Kernel],
        methods: Optional[List[ODEMethod]] = None, ivps: Optional[List[IVP]] = None):
    """Train database with implementation variant runtime predictions.

    Parameters:
    -----------
    db_session: sqlalchemy.orm.session.Session
        Used database session.
    machine: MachineState
        Used machine.
    skeletons: List of ImplSkeleton
        Trained ImplSkeleton objects.
    methods: List of ODEMethod.
        Used ODe methods.
    ivps: List of IVP
        Used IVPs.
    predicted_kernels: List of Kernel
        List of kernels for which we obtained new predictions in this run.

    Returns:
    --------
    -
    """
    config: Config = offsite.config.offsiteConfig
    ts = -1
    if config.args.verbose:
        print('  * Implementation variant predictions...', end='', flush=True)
        ts = start_timer()
    # Initialize worker thread pool.
    pool: Pool = Pool(max(cpu_count() - 1, 1))
    # Compute the implementation variant predictions for all implementation skeletons:
    # * Main thread fetches all the data required
    # * while the worker thread pool handles the actual computation of the predictions.
    records: List[ImplVariantRecord] = list()
    errors: List[str] = list()
    for skeleton in skeletons:
        # ... fetch the benchmark data required to compute the communication costs of the impl skeleton.
        benchmark_data: DataFrame = DataFrame()
        for _ in skeleton.communicationOperationsNodeLvl.keys():
            data = OmpBarrierBenchmark().select(db_session, machine)
            if benchmark_data.empty:
                benchmark_data = data
            else:
                benchmark_data = concat((benchmark_data, data), axis=1)
        # ... determine all available impl variants and store in database.
        available_variants, kernel_executions = deduce_available_impl_variants(
            db_session, skeleton, True, predicted_kernels, config.args.filter_yasksite_opt)
        impl_variants: List[ImplVariant] = [impl.to_database(db_session) for impl in available_variants]
        impl_variant_ids: List[int] = [impl.db_id for impl in impl_variants]
        if not impl_variant_ids:
            continue
        for ivp in ivps:
            if not ivp.characteristic.isStencil and skeleton.modelTool == ModelToolType.YASKSITE:
                continue
            if skeleton.modelTool is not ivp.modelTool:
                continue
            for method in methods:
                __push_train_impl_tasks(db_session, pool, machine, skeleton, kernel_executions, impl_variants,
                                        impl_variant_ids, benchmark_data, records, errors, method, ivp)
        commit(db_session)
    # Wait for all threads and collect results.
    pool.close()
    pool.join()
    # Raise error if impl variant prediction failed.
    if errors:
        db_session.rollback()
        raise RuntimeError('Failed to train impl variant predictions: Error in worker threads.')
    # Before training the database, reduce the total number of intervals by combining adjacent intervals giving the
    # same prediction.
    records_df: DataFrame = DataFrame(records)
    records_df = fuse_equal_impl_records(records_df)
    # Train database with impl variant runtime predictions.
    ImplVariantRecord.update(db_session, records_df)
    # Commit to database.
    commit(db_session)
    if config.args.verbose:
        print('done: Duration {}.'.format(stop_timer(ts)))


def __push_train_impl_tasks(db_session: Session, pool: Pool, machine: MachineState, skeleton: ImplSkeleton,
                            kernel_executions,
                            impl_variants: List[ImplVariant], impl_variant_ids: List[int], benchmark_data, records,
                            errors,
                            method: Optional[ODEMethod] = None, ivp: Optional[IVP] = None):
    config: Config = offsite.config.offsiteConfig
    # Max number of the machine's CPU cores used.
    max_cores: int = machine.coresPerSocket + 1
    core_count: List[int] = [*range(1, max_cores)]
    # Method and IVP database ID if available.
    meth_id = -100 if method is None else method.db_id
    ivp_id = -100 if ivp is None else ivp.db_id
    # Compute implementation variant prediction
    # ... remove old, (possibly obsolete) impl variant records from the database.
    ImplVariantRecord.remove_records(db_session, impl_variant_ids, machine.db_id, machine.compiler.db_id, meth_id,
                                     ivp_id, machine.clock, core_count, 'MODEL')
    # ... fetch the kernel runtime prediction data from the database and
    pred_data = {cores: fetch_kernel_runtime_prediction_data(
        db_session, machine.db_id, machine.compiler.db_id, meth_id, ivp_id, cores, machine.clock, mode=None,
        ode_size=config.args.ode_size) for cores in core_count}
    # ... use the worker threads to compute impl variant predictions.
    pool.apply_async(
        compute_impl_variant_predictions, callback=records.extend, error_callback=errors.append,
        args=(impl_variants, skeleton, method, ivp, machine, benchmark_data, kernel_executions, pred_data))


def compute_impl_variant_predictions(
        impl_variants: List[ImplVariant], skeleton: ImplSkeleton, method: ODEMethod, ivp: IVP, machine: MachineState,
        bench_data: DataFrame, kernel_executions: Dict[int, Union[str, int, float]],
        prediction_data: Dict[int, Dict[int, DataFrame]]) -> List[Tuple]:
    """
    Compute impl variant runtime predictions for a set of impl variants and a given configuration of machine state, ODE
    method and IVP.

    Parameters:
    -----------
    impl_variants: list of ImplVariant
        List of available impl variants derived from the used impl skeleton.
    skeleton: ImplSkeleton
        Used impl skeleton.
    method: ODE method
        Used ODE method.
    ivp: IPV
        Used IVP.
    machine: MachineState
        Used machine.
    bench_data: pandas.DataFrame
        Benchmark data obtained for the communication operations required by the impl skeleton used.
    kernel_executions: dict (key=int, value=sympy.Basic)
        Number of times the single kernels associated with this impl skeleton are executed in a single
        iteration step.
    prediction_data: dict of dict (value=pandas.DataFrame)
        Kernel prediction data required to compute impl variant runtime predictions for the used
        impl skeleton.

    Returns:
    --------
    List of ImplVariant
        Computed impl variant runtime predictions.
    """
    # MachineState frequency used.
    frequency = machine.clock
    # Max number of the machine's CPU cores used.
    max_cores = machine.coresPerSocket + 1
    try:
        # Compute the communication costs of this impl skeleton.
        comm_costs_skeleton = compute_node_lvl_communication_costs(
            bench_data, skeleton.communicationOperationsNodeLvl, max_cores, method, ivp)
        # Compute impl variant run predictions.
        records = list()
        for impl in impl_variants:
            comm_costs_variant: Dict[int, float] = deepcopy(comm_costs_skeleton)
            # If available ...
            if impl.kernel_communication_node_lvl:
                # ... add the kernel's communication costs ...
                comm_costs_kernel = compute_node_lvl_communication_costs(
                    bench_data, impl.kernel_communication_node_lvl, max_cores, method, ivp)
                # ... to the skeleton's costs to get the total communication costs of the impl variant.
                assert (len(comm_costs_kernel) == len(comm_costs_skeleton))
                for cores, costs in comm_costs_kernel.items():
                    comm_costs_variant[cores] = eval_math_expr('{} + {}'.format(comm_costs_skeleton[cores], costs))
            # Pre-compute some constants if available and required.
            constants = list()
            if method is not None:
                constants = [corrector_steps(method), stages(method)]
            if ivp is not None:
                constants.append(ivp_grid_size(ivp.gridSize))
            #
            for cores in range(1, max_cores):
                # Determine number of times each kernel has to be executed for the given configuration.
                executions = {
                    i: eval_math_expr(e, constants) for i, e in kernel_executions.items() if i in impl.kernels}
                # Select the kernel runtime prediction data of this variant's kernels.
                impl_pred_data = prediction_data[cores].loc[prediction_data[cores]['kernel'].isin(impl.kernels)]
                # ... for all sample intervals.
                for interval, data in deduce_impl_variant_sample_intervals(impl_pred_data).items():
                    # Compute the impl variant runtime prediction (in seconds) of this variant using the
                    # kernel runtime predictions of its kernels as well as its communication costs.
                    pred: str = compute_impl_variant_pred(data, executions, comm_costs_variant[cores])
                    #
                    record = {'impl': impl.db_id, 'machine': machine.db_id, 'compiler': machine.compiler.db_id,
                              'method': method.db_id, 'ivp': ivp.db_id, 'cores': cores, 'frequency': frequency,
                              'first': interval.first, 'last': interval.last, 'prediction': str(pred),
                              'pred_mode': 'MODEL', 'updatedIn': __version__,
                              'updatedOn': datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f'), 'updatedBy': getuser()}
                    records.append(record)
        return records
    except Exception as e:
        # Print stack trace of the executing worker thread.
        print_exc()
        print('')
        raise e


def train_impl_variant_runtimes(db_session: Session, machine: MachineState, skeletons: List[ImplSkeleton],
                                methods: Optional[List[ODEMethod]] = None, ivps: Optional[List[IVP]] = None):
    """Train database with benchmarked implementation variant runtimes.

    Parameters:
    -----------
    db_session: sqlalchemy.orm.session.Session
        Used database session.
    machine: MachineState
        Used machine.
    skeletons: List of ImplSkeleton
        Trained ImplSkeleton objects.
    methods: List of ODEMethod.
        Used ODe methods.
    ivps: List of IVP
        Used IVPs.

    Returns:
    --------
    -
    """
    config: Config = offsite.config.offsiteConfig
    n: int = config.args.ode_size
    cores: int = config.args.cores
    # TODO kmp_affinity
    # TODO check if ODE size is applicable.

    # Write header file and main benchmark source code file.
    folder: Path = Path('tmp/bench_variants')
    _write_bench_utils_header(folder)
    bench_src_stub: str = _write_impl_bench_main_code()

    # Initialize implementation variant code generator.
    codegen: ImplCodeGeneratorC = ImplCodeGeneratorC(db_session, folder, folder, folder)

    # Pin CPU frequency for benchmarking.
    if config.args.verbose:
        print('\nPinning CPU frequency to {} Hz for implementation variant benchmark.'.format(machine.clock))
    # MachineState.pin_cpu_frequency(machine.clock)

    # Benchmark the implementation variant runtimes for all implementation skeletons:
    records: List[ImplVariantRecord] = list()
    for skeleton in skeletons:
        # ... determine all available impl variants and store in database.
        available_variants, _ = deduce_available_impl_variants(db_session, skeleton, False)
        impl_variants: List[ImplVariant] = [impl.to_database(db_session) for impl in available_variants]
        if not impl_variants:
            continue
        impl_variants_kernels = [(impl.db_id, impl.kernels) for impl in impl_variants]
        print(impl_variants_kernels)
        for ivp in ivps:
            if not ivp.characteristic.isStencil and skeleton.modelTool == ModelToolType.YASKSITE:
                continue
            if skeleton.modelTool is not ivp.modelTool:
                continue
            for method in methods:
                # Generate instrumented implementation variant code.
                codes = codegen.generate(skeleton, impl_variants_kernels, False, ivp, method)
                write_codes_to_file(codes, suffix='')

                # Benchmark implementation variants.
                records: List[Dict] = list()
                for impl in impl_variants:
                    # Write specific benchmark file for this implementation variant.
                    variant_name: str = ImplVariant.fetch_impl_variant_name(db_session, impl.db_id)
                    bench_src: str = bench_src_stub.replace('VARIANT_HEADER', '"{}.h"'.format(variant_name))
                    bench_file: Path = write_code_to_file(bench_src, 'bench_impl', folder)

                    # Compile code.
                    binary: Path = folder / Path('bench.out')
                    # Construct compiler call.
                    cmd: List[str] = [machine.compiler.name]
                    cmd.extend((flag for flag in machine.compiler.flags.split(' ')))
                    cmd.append(str(bench_file))
                    cmd.append('-Dn={}'.format(n))
                    cmd.append('-Dg={}'.format(eval_math_expr(ivp.gridSize, [ivp_system_size(n)], cast_to=int)))
                    cmd.append('-lm')
                    cmd.append('-o{}'.format(binary))
                    # Compile.
                    run_process(cmd)

                    # Run benchmark for current configuration.
                    cmd = ['./{}'.format(binary), str(config.repetitions_implementation_variants), str(cores)]
                    output: str = run_process(cmd)

                    # Process benchmark results.
                    times: List[float] = remove_outliers([float(x) for x in output.split('\n') if x != ''])
                    runtime: float = mean(times)
                    runtime_per_component: str = eval_math_expr('{} * n'.format(runtime / n), cast_to=str)
                    # ... and store them.
                    record = {'impl': impl.db_id, 'machine': machine.db_id, 'compiler': machine.compiler.db_id,
                              'method': method.db_id, 'ivp': ivp.db_id, 'cores': cores, 'frequency': machine.clock,
                              'first': n, 'last': n, 'prediction': runtime_per_component, 'pred_mode': 'BENCH',
                              'updatedIn': __version__, 'updatedOn': datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f'),
                              'updatedBy': getuser()}
                    records.append(record)
        commit(db_session)

    # Reset CPU frequency limits after benchmarking.
    if config.args.verbose:
        print('\nResetting CPU frequency limits after implement variant benchmark')
    # MachineState.reset_cpu_frequency_limits()

    # Train database with implementation variant runtimes.
    records_df: DataFrame = DataFrame(records)
    ImplVariantRecord.update(db_session, records_df)
    # Commit to database.
    commit(db_session)


def _write_impl_bench_main_code() -> str:
    # Write header includes.
    code = '#include <omp.h>\n'
    code += '#include <stdio.h>\n'
    code += '#include "offsite_util.h"\n'
    code += '\n'
    code += '#include VARIANT_HEADER\n'
    code += '\n'
    # Add main function code.
    code += 'int main(int argc, char **argv) {\n'
    code += 'int threads = atoi(argv[2]);\n'
    code += '\n'
    code += 'double times[atoi(argv[1])];\n'
    code += '\n'
    code += 'time_snap_t ts;\n'
    code += '\n'
    code += '#ifdef _OPENMP\n'
    code += '#pragma omp parallel num_threads(threads)\n'
    code += '{\n'
    code += '#pragma omp single\n'
    code += 'threads = omp_get_num_threads();\n'
    code += '}\n'
    code += '#else\n'
    code += '{\n'
    code += 'cpu_set_t mask;\n'
    code += 'CPU_ZERO( & mask);\n'
    code += 'CPU_SET(0, & mask);\n'
    code += 'sched_setaffinity(0, sizeof(cpu_set_t), & mask);\n'
    code += '}\n'
    code += '#endif\n'
    code += '\n'
    code += 'allocate_data_structures();\n'
    code += '\n'
    code += 'initial_values(0.123, y);\n'
    code += '\n'
    code += 'init_time_snap();\n'
    code += '\n'
    # Open parallel region.
    code += '#pragma omp parallel num_threads(threads)\n'
    code += '{\n'
    # Data distribution.
    code += '#ifdef _OPENMP\n'
    code += 'const int me = omp_get_thread_num();\n'
    code += '# else\n'
    code += 'const int me = 0;\n'
    code += '#endif\n'
    code += '\n'
    code += 'const int first = me * n / threads;\n'
    code += 'const int last = (me + 1) * n / threads - 1;\n'
    code += '\n'
    # Run benchmark.
    code += '#pragma omp barrier\n'
    code += 'for (int warmup = 1; warmup >= 0; --warmup)\n'
    code += '{\n'
    code += 'timestep(me, first, last, 1.234, 0.123);\n'
    code += '}\n'
    code += '\n'
    code += 'int repeat = atoi(argv[1]);\n'
    code += 'for (; repeat > 0; --repeat)\n'
    code += '{\n'
    code += '#pragma omp barrier\n'
    code += 'if (me == 0)\n'
    code += 'time_snap_start(&ts);\n'
    code += '\n'
    code += 'timestep(me, first, last, 1.234, 0.123);\n'
    code += '#pragma omp barrier\n'
    code += '\n'
    code += 'if (me == 0)\n'
    code += '{\n'
    code += 'times[repeat-1] = time_snap_stop(&ts) / 1e9;\n'
    code += '}\n'
    code += '}\n'
    code += '}\n'
    code += '\n'
    # Print results.
    code += 'double total = 0.0;\n'
    code += 'for (int i = 0; i < atoi(argv[1]); ++i)\n'
    code += '{\n'
    code += 'printf("%.15f\§n", times[i]);\n'
    code += 'total += times[i] / (double) n;\n'
    code += '}\n'
    code += '\n'
    code += 'free_data_structures();\n'
    code += '\n'
    code += 'return 0;\n'
    code += '}\n'
    code += '\n'
    return code


def _write_bench_utils_header(folder: Path):
    code = '#ifndef OFFSITE_UTIL_H_\n'
    code += '#define OFFSITE_UTIL_H_\n'
    code += '\n'
    code += '#include <assert.h>\n'
    code += '#include <stddef.h>\n'
    code += '#include <stdint.h>\n'
    code += '#include <stdlib.h>\n'
    code += '#include <sys/time.h>\n'
    code += '\n'
    code += '#define ALIGNMENT 32\n'
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
    code += 'double *alloc1d(size_t a)\n'
    code += '{\n'
    code += 'return (double *) aligned_alloc(ALIGNMENT, a * sizeof(double));\n'
    code += '}\n'
    code += '\n'
    code += 'double **alloc2d(size_t a, size_t b)\n'
    code += '{\n'
    code += 'size_t i, row_size, row_count;\n'
    code += 'double ** x;\n'
    code += 'row_size = ALIGN(b * sizeof(double), ALIGNMENT);\n'
    code += 'row_count = row_size / sizeof(double);\n'
    code += 'x = (double **)\n'
    code += 'aligned_alloc(ALIGNMENT, a * sizeof(double *)); \n'
    code += 'x[0] = (double *)\n'
    code += 'aligned_alloc(ALIGNMENT, a * row_size);\n'
    code += 'for (i = 1; i < a; i++)\n'
    code += 'x[i] = x[0] + i * row_count;\n'
    code += 'return x;\n'
    code += '}\n'
    code += '\n'
    code += 'double ***alloc3d(size_t a, size_t b, size_t c)\n'
    code += '{\n'
    code += 'size_t i, row_size, row_count;\n'
    code += 'double ***x;\n'
    code += 'assert ((ALIGNMENT % sizeof(double)) == 0);\n'
    code += 'row_size = ALIGN(c * sizeof(double), ALIGNMENT);\n'
    code += 'row_count = row_size / sizeof(double);\n'
    code += 'x = (double ** *)\n'
    code += 'aligned_alloc(ALIGNMENT, a * sizeof(double **));\n'
    code += 'x[0] = (double **)\n'
    code += 'aligned_alloc(ALIGNMENT, a * b * sizeof(double *));\n'
    code += 'x[0][0] = (double *)\n'
    code += 'aligned_alloc(ALIGNMENT, a * b * row_size);\n'
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
    code += '#endif\n'
    # Write to file.
    write_code_to_file(code, 'offsite_util', folder, suffix='.h')
