"""@package train.node.train_impl
Functions to train the tuning database with implementation variant predictions.
"""

from copy import deepcopy
from datetime import timedelta
from multiprocessing import cpu_count, Pool
from time import time
from traceback import print_exc
from typing import Dict, List, Optional, Union

from pandas import DataFrame, concat
from sqlalchemy.orm import Session

import offsite.config
from offsite.config import Config, ModelToolType
from offsite.database import commit
from offsite.descriptions.impl.impl_skeleton import ImplSkeleton
from offsite.descriptions.impl.kernel_template import Kernel
from offsite.descriptions.machine import MachineState
from offsite.descriptions.ode import IVP, ivp_grid_size, ODEMethod, corrector_steps, stages
from offsite.solver import SolverType
from offsite.train.communication.openmp.omp_barrier import OmpBarrierBenchmark
from offsite.train.impl_variant import ImplVariant, ImplVariantRecord
from offsite.train.node.util.communication_costs import compute_node_lvl_communication_costs
from offsite.train.node.util.performance_model import compute_impl_variant_pred
from offsite.train.train_utils import deduce_available_impl_variants, deduce_impl_variant_sample_intervals, \
    fetch_kernel_runtime_prediction_data
from offsite.util.math_utils import eval_math_expr


def train_impl_variant(
        db_session: Session, machine: MachineState, skeletons: List[ImplSkeleton], predicted_kernels: List[Kernel],
        ivps: Optional[List[IVP]] = None, methods: List[ODEMethod] = None):
    """Train database with implementation variant runtime predictions.

    Parameters:
    -----------
    db_session: sqlalchemy.orm.session.Session
        Used database session.
    machine: MachineState
        Used machine.
    skeletons: List of ImplSkeleton
        Trained ImplSkeleton objects.
    predicted_kernels: List of Kernel
        Kernels for which new predictions were obtained in this run.
    methods: List of ODEMethod.
        Used ODe methods.
    ivps: List of IVP
        Used IVPs.

    Returns:
    --------
    -
    """
    config: Config = offsite.config.offsiteConfig
    if config.args.verbose:
        print('  * Implementation variant predictions...', end='', flush=True)
        start_time_impl = time()
    if config.scenario.solver.type == SolverType.ODE:
        train_impl_variant_predictions_ode(db_session, machine, skeletons, methods, ivps, predicted_kernels)
    elif config.scenario.solver.type == SolverType.GENERIC:
        train_impl_variant_predictions(db_session, machine, skeletons, predicted_kernels)
    # Commit to database.
    commit(db_session)
    if config.args.verbose:
        stop_time_impl = time()
        print('done: Duration {}.'.format(timedelta(seconds=round(stop_time_impl - start_time_impl, 0))))


def train_impl_variant_predictions_ode(db_session: Session, machine: MachineState, skeletons: List[ImplSkeleton],
                                       methods: List[ODEMethod], ivps: List[IVP], predicted_kernels: List[Kernel]):
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
    # Max number of the machine's CPU cores used.
    max_cores: int = machine.coresPerSocket + 1
    core_count: List[int] = [*range(1, max_cores)]
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
                # ... remove old, (possibly obsolete) impl variant records from the database.
                ImplVariantRecord.remove_records(db_session, impl_variant_ids, machine.db_id, machine.compiler.db_id,
                                                 method.db_id, ivp.db_id, machine.clock, core_count)
                # ... fetch the kernel runtime prediction data from the database.
                pred_data = {cores: fetch_kernel_runtime_prediction_data(
                    db_session, machine.db_id, machine.compiler.db_id, method.db_id, ivp.db_id, cores, machine.clock,
                    mode=None, ode_size=config.args.ode_size) for cores in core_count}
                # ... use worker threads to compute impl variant predictions.
                pool.apply_async(
                    compute_impl_variant_predictions_ode, callback=records.extend, error_callback=errors.append,
                    args=(impl_variants, skeleton, method, ivp, machine, benchmark_data, kernel_executions, pred_data))
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
    records = ImplVariantRecord.fuse_equal_records(records)
    # Train database with impl variant runtime predictions.
    ImplVariantRecord.update(db_session, records)


def compute_impl_variant_predictions_ode(
        impl_variants: List[ImplVariant], skeleton: ImplSkeleton, method: ODEMethod, ivp: IVP, machine: MachineState,
        bench_data: DataFrame, kernel_executions: Dict[int, Union[str, int, float]],
        prediction_data: Dict[int, Dict[int, DataFrame]]) -> List[ImplVariantRecord]:
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
            #
            for cores in range(1, max_cores):
                # Determine number of times each kernel has to be executed for the given configuration.
                executions = {
                    i: eval_math_expr(e, [corrector_steps(method), stages(method), ivp_grid_size(ivp.gridSize)]) for
                    i, e in kernel_executions.items() if i in impl.kernels}
                # Select the kernel runtime prediction data of this variant's kernels.
                impl_pred_data = prediction_data[cores].loc[prediction_data[cores]['kernel'].isin(impl.kernels)]
                # ... for all sample intervals.
                for interval, data in deduce_impl_variant_sample_intervals(impl_pred_data).items():
                    # Compute the impl variant runtime prediction (in seconds) of this variant using the
                    # kernel runtime predictions of its kernels as well as its communication costs.
                    pred = compute_impl_variant_pred(data, executions, comm_costs_variant[cores])
                    records.append(ImplVariantRecord(impl.db_id, machine.db_id, machine.compiler.db_id, method.db_id,
                                                     ivp.db_id, cores, frequency, interval, str(pred), 'MODEL'))
        return records
    except Exception as e:
        # Print stack trace of the executing worker thread.
        print_exc()
        print('')
        raise e


def train_impl_variant_predictions(
        db_session: Session, machine: MachineState, skeletons: List[ImplSkeleton], predicted_kernels: List[Kernel]):
    """Train database with impl variant runtime predictions.

    Parameters:
    -----------
    db_session: sqlalchemy.orm.session.Session
        Used database session.
    machine: MachineState
        Used machine.
    skeletons: List of ImplSkeleton
        Trained ImplSkeleton objects.
    predicted_kernels: List of Kernel
        List of kernels for which we obtained new predictions in this run.

    Returns:
    --------
    -
    """
    config: Config = offsite.config.offsiteConfig
    # Max number of the machine's CPU cores used.
    max_cores: int = machine.coresPerSocket + 1
    core_count: List[int] = [*range(1, max_cores)]
    # Initialize worker thread pool.
    pool: Pool = Pool(max(cpu_count() - 1, 1))
    # Compute the impl variant predictions for all impl skeletons:
    # * Main thread fetches all the data required
    # * while the worker thread pool handles the actual computation of the predictions.
    records: List[ImplVariantRecord] = list()
    errors: List[str] = list()
    for skeleton in skeletons:
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
        # ... remove old, (possibly obsolete) ranking records from the database.
        ImplVariantRecord.remove_records(db_session, impl_variant_ids, machine.db_id, machine.compiler.db_id, -100,
                                         -100, machine.clock, [*range(1, max_cores)])
        # ... fetch the kernel runtime prediction data from the database.
        pred_data = {cores: fetch_kernel_runtime_prediction_data(
            db_session, machine.db_id, machine.compiler.db_id, -100, -100, cores, machine.clock,
            mode=None, ode_size=config.args.ode_size) for cores in core_count}
        # ... use worker threads to compute impl variant predictions.
        pool.apply_async(
            compute_impl_variant_predictions, callback=records.extend, error_callback=errors.append,
            args=(impl_variants, skeleton, machine, benchmark_data, kernel_executions, pred_data))
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
    records = ImplVariantRecord.fuse_equal_records(records)
    # Train database with impl variant runtime predictions.
    ImplVariantRecord.update(db_session, records)


def compute_impl_variant_predictions(impl_variants: List[ImplVariant], skeleton: ImplSkeleton, machine: MachineState,
                                     bench_data: DataFrame, kernel_executions: Dict[int, Union[str, int, float]],
                                     prediction_data: Dict[int, DataFrame]) -> List[ImplVariantRecord]:
    """
    Compute impl variant runtime predictions for a set of impl variants and a given machine.

    Parameters:
    -----------
    impl_variants: list of ImplVariant
        List of available impl variants derived from the used impl skeleton.
    skeleton: ImplSkeleton
        Used impl skeleton.
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
    config: Config = offsite.config.offsiteConfig
    # MachineState frequency used.
    frequency = machine.clock
    # Max number of the machine's CPU cores used.
    max_cores = machine.coresPerSocket + 1
    try:
        # Compute the communication costs of this impl skeleton.
        comm_costs_skeleton = compute_node_lvl_communication_costs(
            bench_data, skeleton.communicationOperationsNodeLvl, max_cores, None, None)
        # Compute impl variant run predictions
        records = list()
        for impl in impl_variants:
            comm_costs_variant: Dict[int, float] = deepcopy(comm_costs_skeleton)
            # If available ...
            if impl.kernel_communication_node_lvl:
                # ... add the kernel's communication costs ...
                comm_costs_kernel = compute_node_lvl_communication_costs(
                    bench_data, impl.kernel_communication_node_lvl, max_cores, None, None)
                # ... to the skeleton's costs to get the total communication costs of the impl variant.
                assert (len(comm_costs_kernel) == len(comm_costs_skeleton))
                for cores, costs in comm_costs_kernel.items():
                    comm_costs_variant[cores] = eval_math_expr('{} + {}'.format(comm_costs_skeleton[cores], costs))
            #
            for cores in range(1, max_cores):
                # Determine number of times each kernel has to be executed for the given configuration.
                executions = {i: eval_math_expr(e, []) for i, e in kernel_executions.items() if
                              i in impl.kernels}
                # Select the kernel runtime prediction data of this variant's kernels.
                impl_pred_data = prediction_data[cores].loc[prediction_data[cores]['kernel'].isin(impl.kernels)]
                # ... for all sample intervals.
                for interval, data in deduce_impl_variant_sample_intervals(impl_pred_data).items():
                    # Compute the impl variant runtime prediction (in seconds) of this variant using
                    # the kernel runtime predictions of its kernels as well as its communication costs.
                    pred = compute_impl_variant_pred(data, executions, comm_costs_variant[cores])
                    records.append(ImplVariantRecord(impl.db_id, machine.db_id, machine.compiler.db_id, -1010, -100,
                                                     cores, frequency, interval, str(pred), config.args.mode.value))
        return records
    except Exception as e:
        # Print stack trace of the executing worker thread.
        print_exc()
        print('')
        raise e
