"""@package train_impl
Functions to train the tuning database with implementation variant predictions.
"""

from multiprocessing import cpu_count, Pool
from traceback import print_exc
from typing import Dict, List

from sqlalchemy.orm import Session

import offsite.config
from offsite.config import ModelToolType
from offsite.db.db import commit
from offsite.descriptions.impl_skeleton import ImplSkeleton
from offsite.descriptions.ivp import IVP
from offsite.descriptions.machine import Machine
from offsite.descriptions.ode_method import ODEMethod
from offsite.evaluation.benchmark import BenchmarkRecord
from offsite.evaluation.math_utils import eval_math_expr, corrector_steps, ivp_grid_size, stages
from offsite.evaluation.performance_model import compute_impl_variant_runtime_pred, ImplVariantRecord
from offsite.train.train_utils import deduce_available_impl_variants, deduce_impl_variant_sample_intervals, \
    fetch_and_sort_kernel_runtime_prediction_data


def train_impl_variant_predictions(db_session: Session, machine: Machine, skeletons: List[ImplSkeleton],
                                   methods: List[ODEMethod], ivps: List[IVP]):
    """Train database with implementation variant runtime predictions.

    Parameters:
    -----------
    db_session : sqlalchemy.orm.session.Session
        Used database session.
    machine : Machine
        Used machine.
    skeletons : List of ImplSkeleton
        Trained ImplSkeleton objects.
    methods : List of ODEMethod.
        Used ode methods.
    ivps : List of IVP
        Used IVPs.

    Returns:
    --------
    -
    """
    config = offsite.config.offsiteConfig
    # Machine frequency used.
    frequency = machine.clock
    # Max number of the machine's CPU cores used.
    max_cores = machine.coresPerSocket + 1
    # Fetch required data ...
    benchmark_data = dict()
    impl_variants = dict()
    kernel_executions = dict()
    prediction_data = dict()
    for skeleton in skeletons:
        # ... fetch the benchmark data required to compute the communication costs of the implementation skeleton.
        benchmark_data[skeleton.db_id] = BenchmarkRecord.select(
            db_session, machine, [bench_name for bench_name in skeleton.communicationOperations.keys()])
        # ... determine all available implementation variants and store in database.
        available_variants, kernel_executions[skeleton.db_id] = deduce_available_impl_variants(
            db_session, skeleton, config.args.filter_yasksite_opt)
        impl_variants[skeleton.db_id] = [impl.to_database(db_session) for impl in available_variants]
        # ... fetch and sort kernel runtime prediction data from the database.
        prediction_data[skeleton.db_id] = dict()
        for ivp in ivps:
            if not ivp.characteristic.isStencil and skeleton.modelTool == ModelToolType.YASKSITE:
                continue
            if skeleton.modelTool is not ivp.modelTool:
                continue
            prediction_data[skeleton.db_id][ivp.db_id] = dict()
            for method in methods:
                # Remove old, (possibly obsolete) ranking records from the database.
                ImplVariantRecord.remove_records(
                    db_session, [impl.db_id for impl in available_variants if impl], machine.db_id,
                    machine.compiler.db_id, method.db_id,
                    ivp.db_id, frequency, [*range(1, max_cores)])
                #
                prediction_data[skeleton.db_id][ivp.db_id][method.db_id] = dict()
                for cores in range(1, max_cores):
                    prediction_data[skeleton.db_id][ivp.db_id][method.db_id][cores] = \
                        fetch_and_sort_kernel_runtime_prediction_data(
                            db_session, machine.db_id, machine.compiler.db_id, method.db_id, ivp.db_id, cores,
                            frequency, config.args.mode.value, config.args.ode_size)
        commit(db_session)
    # Train records.
    records = list()
    errors = list()
    # Initialize worker thread pools.
    pool = Pool(cpu_count())
    # Compute implementation variant runtime predictions for ...
    # ... all implementation skeletons
    for skeleton in skeletons:
        # ... all IVPs
        for ivp in ivps:
            if not ivp.characteristic.isStencil and skeleton.modelTool == ModelToolType.YASKSITE:
                continue
            if skeleton.modelTool is not ivp.modelTool:
                continue
            # ... all ODE methods
            for method in methods:
                pool.apply_async(compute_impl_variant_predictions,
                                 args=(impl_variants[skeleton.db_id], skeleton, method, ivp, machine,
                                       benchmark_data[skeleton.db_id], kernel_executions[skeleton.db_id],
                                       prediction_data[skeleton.db_id][ivp.db_id][method.db_id]),
                                 callback=records.extend, error_callback=errors.append)
    # Wait for all threads and collect results.
    pool.close()
    pool.join()
    # Raise error if implementation variant prediction failed.
    if errors:
        db_session.rollback()
        raise RuntimeError('Failed to train implementation variant predictions: Error in worker threads.')
    # Train database with implementation variant runtime predictions.
    ImplVariantRecord.update(db_session, records)


def compute_impl_variant_predictions(
        impl_variants: List['ImplVariant'], skeleton: ImplSkeleton, method: ODEMethod, ivp: IVP, machine: Machine,
        benchmark_data: 'pandas.DataFrame', kernel_executions: Dict[int, 'sympy.Basic'],
        prediction_data: Dict[int, Dict[int, List['KernelRecord']]]) -> List[ImplVariantRecord]:
    """
    Compute implementation variant runtime predictions for a set of implementation variants and a given configuration
    of machine, ODE method and IVP.

    Parameters:
    -----------
    impl_variants
        List of available implementation variants derived from the used implementation skeleton.
    skeleton : ImplSkeleton
        Used implementation skeleton.
    method : ODE method
        Used ODE method.
    ivp : IPV
        Used IVP.
    machine : Machine
        Used machine.
    benchmark_data : pandas.DataFrame
        Benchmark data obtained for the communication operations required by the implementation skeleton used.
    kernel_executions : dict (key=int, value=sympy.Basic)
        Number of times the single kernels associated with this implementation skeleton are executed in a single
        iteration step.
    prediction_data : dict of dict (value=list of KernelRecord)
        Kernel prediction data required to compute implementation variant runtime predictions for the used
        implementation skeleton.

    Returns:
    --------
    List of ImplVariant
        Computed implementation variant runtime predictions.
    """
    config = offsite.config.offsiteConfig
    # Machine frequency used.
    frequency = machine.clock
    # Max number of the machine's CPU cores used.
    max_cores = machine.coresPerSocket + 1
    try:
        # Compute the communication costs of a single iteration step of this implementation skeleton for the given
        # configuration of ODE method and IVP.
        communication_costs = skeleton.compute_communication_costs(benchmark_data, method, ivp)
        # Compute implementation variant run predictions
        records = list()
        # .. all number of cores.
        for cores in range(1, max_cores):
            pred_data = prediction_data[cores]
            for impl in impl_variants:
                # Determine number of times each kernel has to be executed for the given configuration.
                executions = {
                    i: eval_math_expr(e, [corrector_steps(method), stages(method), ivp_grid_size(ivp.gridSize)]) for
                    i, e in kernel_executions.items() if i in impl.kernels}
                # Select the kernel runtime prediction data of this variant's kernels.
                impl_pred_data = {key: pred_data[key] for key in pred_data if key in impl.kernels}
                # ... for all sample intervals.
                for interval, data in deduce_impl_variant_sample_intervals(impl_pred_data).items():
                    # Compute the implementation variant runtime prediction (in seconds) of this variant using
                    # the kernel runtime predictions of its kernels as well as its communication costs.
                    pred = compute_impl_variant_runtime_pred(impl.kernels, data, executions, communication_costs[cores])
                    records.append(
                        ImplVariantRecord(impl.db_id, machine.db_id, machine.compiler.db_id, method.db_id, ivp.db_id,
                                          cores, frequency, interval, str(pred), config.args.mode.value))
        return records
    except Exception as e:
        # Print stack trace of the executing worker thread.
        print_exc()
        print('')
        raise e
