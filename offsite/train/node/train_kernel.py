"""@package train.node.train_kernel
Functions to train the tuning database with kernel predictions.

@author: Johannes Seiferth
"""

from multiprocessing import cpu_count, Pool
from os import putenv
from statistics import mean
from traceback import print_exc
from typing import Dict, List, Optional, Set, Tuple

from pathlib2 import Path
from sqlalchemy.orm import Session

import offsite.config
from offsite.codegen.generator.kernel_bench_generator import KernelBenchFiles
from offsite.config import Config, ModelToolType, ProgramModeType
from offsite.database import commit
from offsite.descriptions.impl.kernel_template import Kernel, KernelTemplate, PModelKernel, KernelRecord, PModelRecord
from offsite.descriptions.machine import MachineState
from offsite.descriptions.ode import IVP, ODEMethod
from offsite.solver import SolverType
from offsite.train.node.util.kerncraft_utils import execute_kerncraft_bench_mode, execute_kerncraft_ecm_mode, \
    parse_kerncraft_output_bench_mode, parse_kerncraft_output_ecm_mode
from offsite.train.node.util.performance_model import compute_pmodel_kernel_pred, compute_kernel_pred
from offsite.train.node.util.yasksite_utils import execute_yasksite_bench_mode, execute_yasksite_ecm_mode, \
    parse_yasksite_output
from offsite.train.train_utils import fuse_equal_records
from offsite.util.math_utils import remove_outliers
from offsite.util.process_utils import run_process
from offsite.util.sample_interval import SampleInterval, SampleType
from offsite.util.time import start_timer, stop_timer

IntervalRecordList = List[Tuple[SampleInterval, str]]
PModelPredictionDict = Dict[int, List[str]]


def train_kernel(db_session: Session, machine: MachineState, templates: List[KernelTemplate],
                 ivps: Optional[List[IVP]] = None, methods: Optional[List[ODEMethod]] = None):
    """Train database with kernel runtime predictions.

    Parameters:
    -----------
    db_session: sqlalchemy.orm.session.Session
        Used database session.
    machine: MachineState
        Used MachineState.
    templates: List of KernelTemplate
        Kernel templates to be trained.
    methods: List of ODEMethod.
        Used ODE methods.
    ivps: List of IVP
        Used IVPs.

    Returns:
    --------
    predicted_kernels: List of Kernel
        List of kernels for which we obtained new predictions in this run.
    """
    config: Config = offsite.config.offsiteConfig
    ts = -1
    if config.args.verbose:
        print('  * Kernel predictions...', end='', flush=True)
        ts = start_timer()
    predicted_kernels: Set[Kernel] = set()
    if config.args.mode == ProgramModeType.MODEL:
        if config.scenario.solver.type == SolverType.ODE:
            predicted_kernels = train_kernel_predictions_ode(db_session, machine, templates, methods, ivps)
        elif config.scenario.solver.type == SolverType.GENERIC:
            predicted_kernels = train_kernel_predictions_generic(db_session, machine, templates)
    elif config.args.mode == ProgramModeType.RUN:
        if config.scenario.solver.type == SolverType.ODE:
            train_kernel_runtimes_ode(db_session, machine, templates, methods, ivps)
        elif config.scenario.solver.type == SolverType.GENERIC:
            train_kernel_runtimes_generic(db_session, machine, templates)
    # Commit to database.
    commit(db_session)
    if config.args.verbose:
        print(' done: Duration {}.'.format(stop_timer(ts)))
    return predicted_kernels


def train_kernel_predictions_ode(db_session: Session, machine: MachineState, templates: List[KernelTemplate],
                                 methods: List[ODEMethod], ivps: List[IVP]) -> Set[Kernel]:
    """Train database with kernel runtime predictions.

    Compute the kernel runtime predictions for all configuration of kernels, IVP, ODE method and machine state.

    Parameters:
    -----------
    db_session: sqlalchemy.orm.session.Session
        Used database session.
    machine: MachineState
        Used MachineState.
    templates: List of KernelTemplate
        Kernel templates to be trained.
    methods: List of ODEMethod.
        Used ODE methods.
    ivps: List of IVP
        Used IVPs.

    Returns:
    --------
    predicted_kernels: Set of Kernel
        Set of kernels for which we obtained new predictions in this run.
    """
    config: Config = offsite.config.offsiteConfig
    data_kc = list()
    data_ys = list()
    errors = list()
    # Initialize worker thread pools.
    kernel_count = sum((len(t.variants) for t in templates))
    pool_kc = None
    pool_ys = None
    if not config.pred_model_tool:
        # Determine number of workers.
        # One worker handles YASKSITE tasks ...
        num_workers_ys = 1
        # ... while remaining workers handle KERNCRAFT task.
        num_workers_kc = min(cpu_count() - 1, kernel_count)
        num_workers_kc = max(num_workers_kc - 1, 1)
        # Spawn worker pools.
        pool_kc = Pool(num_workers_kc)
        pool_ys = Pool(num_workers_ys)
    elif config.pred_model_tool == ModelToolType.KERNCRAFT:
        # Determine number of workers.
        num_workers_kc = min(cpu_count() - 1, kernel_count)
        # Spawn worker pool.
        pool_kc = Pool(num_workers_kc)
    elif config.pred_model_tool == ModelToolType.YASKSITE:
        num_workers_ys = 1
        # Spawn worker pool.
        pool_ys = Pool(num_workers_ys)
    assert (pool_kc is not None or pool_ys is not None)
    # Compute kernel runtime predictions for ...
    benched_rhs_kernels = list()
    predicted_kernels = set()
    # ... all ODE methods.
    for method in methods:
        # ... all kernel templates.
        for template in templates:
            # ... all kernels.
            for kernel in template.variants:
                tool = template.modelTool
                # ... all IVPs.
                if template.isIVPdependent:
                    # If the template contains IVP calls, a separate kernel prediction has to be computed for each
                    # evaluated IVP.
                    for ivp in ivps:
                        if tool == ModelToolType.YASKSITE and not ivp.characteristic.isStencil:
                            print('\nSkipping kernel template \'{}\'! Unable to predict '.format(template.name), end='')
                            print('using YASKSITE: {} is not a stencil!'.format(ivp.name))
                            continue
                        if tool is not ivp.modelTool:
                            print('\nSkipping kernel template \'{}\'! Unable to predict '.format(template.name), end='')
                            print('using model tool  \'{}\': IVP \'{}\' requires model tool \'{}\'!'.format(
                                tool.value, ivp.name, ivp.modelTool.value))
                            continue
                        # If 'no_update_mode' is set, check if ECM predictions for this kernel are already in the
                        # database for the given run configuration.
                        if config.args.update is True:
                            if KernelRecord.contains(
                                    db_session, kernel.db_id, machine.db_id, machine.compiler.db_id, method.db_id,
                                    ivp.db_id, machine.clock, machine.coresPerSocket, config.args.mode.value,
                                    config.args.ode_size, [SampleType.MODEL_INNER, SampleType.MODEL_BORDER]):
                                continue
                        predicted_kernels.add(kernel.db_id)
                        # Push kernel prediction computation task to thread pool.
                        if tool == ModelToolType.KERNCRAFT:
                            if config.args.bench_rhs:
                                benched_rhs_kernels.append((kernel, machine, method, ivp))
                            else:
                                pool_kc.apply_async(compute_kernel_prediction, args=(kernel, machine, method, ivp),
                                                    callback=data_kc.append, error_callback=errors.append)
                        elif tool == ModelToolType.YASKSITE:
                            pool_ys.apply_async(compute_kernel_prediction, args=(kernel, machine, method, ivp),
                                                callback=data_ys.append, error_callback=errors.append)
                        else:
                            raise RuntimeError('Unsupported model tool type \'{}\' detected!'.format(tool))
                else:
                    # If 'no_update_mode' is set, check if ECM predictions for this kernel are already in the database
                    # for the given run configuration.
                    if config.args.update is True:
                        if KernelRecord.contains(
                                db_session, kernel.db_id, machine.db_id, machine.compiler.db_id, method.db_id, -1,
                                machine.clock, machine.coresPerSocket, config.args.mode.value, config.args.ode_size,
                                [SampleType.MODEL_INNER, SampleType.MODEL_BORDER]):
                            continue
                    predicted_kernels.add(kernel.db_id)
                    # Push kernel prediction computation task to thread pool.
                    if tool == ModelToolType.KERNCRAFT:
                        pool_kc.apply_async(compute_kernel_prediction, args=(kernel, machine, method, None),
                                            callback=data_kc.append, error_callback=errors.append)
                    elif tool == ModelToolType.YASKSITE:
                        pool_ys.apply_async(compute_kernel_prediction, args=(kernel, machine, method, None),
                                            callback=data_ys.append, error_callback=errors.append)
                    else:
                        raise RuntimeError('Unsupported model tool type \'{}\' detected!'.format(tool))
    # Wait for all threads and collect results.
    accumulated_data = list()
    if not config.pred_model_tool:
        # Wait for threads.
        pool_kc.close()
        pool_ys.close()
        pool_kc.join()
        pool_ys.join()
        # Collect results.
        accumulated_data = data_kc + data_ys
    elif config.pred_model_tool == ModelToolType.KERNCRAFT:
        # Wait for threads.
        pool_kc.close()
        pool_kc.join()
        # Collect results.
        accumulated_data = data_kc
    elif config.pred_model_tool == ModelToolType.YASKSITE:
        # Wait for threads.
        pool_ys.close()
        pool_ys.join()
        # Collect results.
        accumulated_data = data_ys
    # Raise error if kernel prediction failed.
    if errors:
        db_session.rollback()
        raise RuntimeError('Failed to train kernel predictions: Error in worker threads.')
    if config.args.bench_rhs:
        for bench_task in benched_rhs_kernels:
            data = bench_rhs_kernel_runtime_ode(bench_task[0], bench_task[1], bench_task[2], bench_task[3])
            accumulated_data.append(data)
    # Before training the database, reduce the total number of intervals by:
    # * combining adjacent intervals giving the same prediction,
    # * interpolating border intervals.
    for kernel_db_id, run_config_db, kernel_predictions, pmodel_predictions in accumulated_data:
        for pmodel, records in pmodel_predictions.items():
            for record in records:
                PModelRecord.update(db_session, pmodel, run_config_db[0], run_config_db[1], run_config_db[2],
                                    run_config_db[3], record[0], record[1], config.pred_incore_tool.value)
        for cores, records in kernel_predictions.items():
            records = fuse_equal_records(records)
            KernelRecord.update(
                db_session, kernel_db_id, run_config_db[0], run_config_db[1], run_config_db[2], run_config_db[3],
                cores, machine.clock, records, config.args.mode.value)
    return predicted_kernels


def train_kernel_predictions_generic(
        db_session: Session, machine: MachineState, templates: List[KernelTemplate]) -> Set[Kernel]:
    """Train database with kernel runtime predictions.

    Compute the kernel runtime predictions for all configuration of kernels and machine state.

    Parameters:
    -----------
    db_session: sqlalchemy.orm.session.Session
        Used database session.
    machine: MachineState
        Used MachineState.
    templates: List of KernelTemplate
        Kernel templates to be trained.

    Returns:
    --------
    predicted_kernels: Set of Kernel
        Set of kernels for which new predictions were obtained in this run.
    """
    config: Config = offsite.config.offsiteConfig
    data_kc = list()
    data_ys = list()
    errors = list()
    # Initialize worker thread pools.
    kernel_count = sum((len(t.variants) for t in templates))
    pool_kc = None
    pool_ys = None
    if not config.pred_model_tool:
        # Determine number of workers.
        # One worker handles YASKSITE tasks ...
        num_workers_ys = 1
        # ... while remaining workers handle KERNCRAFT task.
        num_workers_kc = min(cpu_count(), kernel_count)
        num_workers_kc = max(num_workers_kc - 1, 1)
        # Spawn worker pools.
        pool_kc = Pool(num_workers_kc)
        pool_ys = Pool(num_workers_ys)
    elif config.pred_model_tool == ModelToolType.KERNCRAFT:
        # Determine number of workers.
        num_workers_kc = min(cpu_count(), kernel_count)
        # Spawn worker pool.
        pool_kc = Pool(num_workers_kc)
    elif config.pred_model_tool == ModelToolType.YASKSITE:
        num_workers_ys = 1
        # Spawn worker pool.
        pool_ys = Pool(num_workers_ys)
    assert (pool_kc is not None or pool_ys is not None)
    # Compute kernel runtime predictions for ...
    predicted_kernels = set()
    # ... all kernel templates.
    for template in templates:
        # ... all kernels.
        for kernel in template.variants:
            tool = template.modelTool
            # If 'no_update_mode' is set, check if ECM predictions for this kernel are already in the database for the
            # given run configuration.
            if config.args.update is True:
                if KernelRecord.contains(
                        db_session, kernel.db_id, machine.db_id, machine.compiler.db_id, -100, -100, machine.clock,
                        machine.coresPerSocket, config.args.mode.value, config.args.ode_size,
                        [SampleType.MODEL_INNER, SampleType.MODEL_BORDER]):
                    continue
            predicted_kernels.add(kernel.db_id)
            # Push kernel prediction computation task to thread pool.
            if tool == ModelToolType.KERNCRAFT:
                pool_kc.apply_async(compute_kernel_prediction, args=(kernel, machine), callback=data_kc.append,
                                    error_callback=errors.append)
            elif tool == ModelToolType.YASKSITE:
                pool_ys.apply_async(compute_kernel_prediction, args=(kernel, machine), callback=data_ys.append,
                                    error_callback=errors.append)
            else:
                raise RuntimeError('Unsupported model tool type \'{}\' detected!'.format(tool))
    # Wait for all threads and collect results.
    accumulated_data = list()
    if not config.pred_model_tool:
        # Wait for threads.
        pool_kc.close()
        pool_ys.close()
        pool_kc.join()
        pool_ys.join()
        # Collect results.
        accumulated_data = data_kc + data_ys
    elif config.pred_model_tool == ModelToolType.KERNCRAFT:
        # Wait for threads.
        pool_kc.close()
        pool_kc.join()
        # Collect results.
        accumulated_data = data_kc
    elif config.pred_model_tool == ModelToolType.YASKSITE:
        # Wait for threads.
        pool_ys.close()
        pool_ys.join()
        # Collect results.
        accumulated_data = data_ys
    # Raise error if kernel prediction failed.
    if errors:
        db_session.rollback()
        raise RuntimeError('Failed to train kernel predictions: Error in worker threads.')
    # Before training the database, reduce the total number of intervals by:
    # * combining adjacent intervals giving the same prediction,
    # * interpolating border intervals.
    for kernel_db_id, run_config_db, kernel_predictions, pmodel_predictions in accumulated_data:
        for pmodel, records in pmodel_predictions.items():
            for record in records:
                PModelRecord.update(db_session, pmodel, run_config_db[0], run_config_db[1], -100, -100, record[0],
                                    record[1], config.pred_incore_tool.value)
        for cores, records in kernel_predictions.items():
            records = fuse_equal_records(records)
            KernelRecord.update(db_session, kernel_db_id, run_config_db[0], run_config_db[1], -100, -100, cores,
                                machine.clock, records, config.args.mode.value)
    return predicted_kernels


def compute_kernel_prediction(
        kernel: Kernel, machine: MachineState, method: Optional[ODEMethod] = None, ivp: Optional[IVP] = None):
    """
    Compute the kernel runtime prediction for a given configuration of kernel, machine state and (if required) IVP and
    ODE method.

    Parameters:
    -----------
    kernel: Kernel
        Trained kernel.
    machine: MachineState
        Trained machine.
    method: ODE Method
        Trained ODE method.
    ivp: IVP
        Trained IVP.

    Returns:
    --------
    -
    """
    config: Config = offsite.config.offsiteConfig
    # MachineState frequency used.
    frequency = machine.clock
    try:
        # If no fixed ODE system size is given, deduce the set of significant sample intervals from the pmodel kernels'
        # working sets.
        if config.args.ode_size:
            ode_size = config.args.ode_size
            intervals = [SampleInterval(ode_size, ode_size, SampleType.MODEL_INNER, ode_size)]
        elif config.args.nrange:
            intervals = config.args.nrange
        else:
            intervals = kernel.deduce_relevant_samples(machine, method, ivp)
        # Initialize lists to store the kernel predictions for different numbers of cores.
        kernel_predictions = dict()
        ecm_records = dict()
        for interval in intervals:
            kernel.generate_pmodel_code(method, ivp, interval.sample)
            # Determine the ECM result for each pmodel kernel using the kerncraft tool.
            pmodel_predictions = dict()
            for pmodel in kernel.pmodel_kernels:
                if pmodel.db_id not in ecm_records:
                    ecm_records[pmodel.db_id] = list()
                pmodel_predictions, ecm_predictions = train_pmodel_prediction(
                    pmodel, machine, interval.sample, pmodel_predictions, method, ivp)
                ecm_records[pmodel.db_id].append((interval, ecm_predictions))
            # Compute the kernel runtime prediction of this variant by combining its pmodel kernel predictions obtained.
            for cores, predictions in pmodel_predictions.items():
                pred = compute_kernel_pred(predictions, frequency)
                if cores not in kernel_predictions:
                    kernel_predictions[cores] = list()
                kernel_predictions[cores].append((interval, pred))
        return kernel.db_id, (machine.db_id, machine.compiler.db_id, -1 if not method else method.db_id,
                              - 1 if not ivp else ivp.db_id), kernel_predictions, ecm_records
    except Exception as e:
        # Print stack trace of the executing worker thread.
        print_exc()
        print('')
        raise e


def train_pmodel_prediction(
        pmodel: PModelKernel, machine: MachineState, trained_size: int, pmodel_kernel_predictions: PModelPredictionDict,
        method: Optional[ODEMethod], ivp: Optional[IVP]) -> PModelPredictionDict:
    """Train database with pmodel kernel runtime prediction.

    Compute the pmodel kernel runtime prediction for a given configuration of kernel, machine state and (if required)
    IVP and ODE method.

    Parameters:
    -----------
    pmodel: PModel
        Trained pmodel kernel.
    machine: MachineState
        Trained machine.
    trained_size: int
        Trained system size.
    pmodel_kernel_predictions: dict (key=int, value=str)
        Trained pmodel kernel runtime predictions sorted by number of CPU cores used.
    method: ODE method
        Trained ODE method.
    ivp: IVP
        Trained IVP.

    Returns:
    --------
    dict (key=int, value=str)
        Trained pmodel kernel runtime predictions.
    """
    # Compute the number of iterations.
    iteration_count = pmodel.iteration_count(method, ivp)
    # Compute the ECM predictions.
    ecm_predictions: Dict[int, float]
    if isinstance(pmodel.code_path, dict):
        ecm_predictions = dict()
        for code_path in pmodel.code_path.values():
            # Execute kerncraft and parse its output.
            out = execute_kerncraft_ecm_mode(code_path, machine, method, ivp, trained_size)
            preds = parse_kerncraft_output_ecm_mode(out)
            ecm_predictions = {k: ecm_predictions.get(k, 0) + preds.get(k, 0)
                               for k in set(ecm_predictions) | set(preds)}
    else:
        if pmodel.kernel.template.modelTool == ModelToolType.KERNCRAFT:
            # Execute kerncraft and parse its output.
            out = execute_kerncraft_ecm_mode(pmodel.code_path, machine, method, ivp, trained_size)
            ecm_predictions = parse_kerncraft_output_ecm_mode(out)
        elif pmodel.kernel.template.modelTool == ModelToolType.YASKSITE:
            out = execute_yasksite_ecm_mode(
                pmodel.code_path, machine, method, ivp, trained_size, pmodel.kernel.optimization_parameters)
            ecm_predictions = parse_yasksite_output(out)
        else:
            assert False
    # Compute pmodel kernel predictions using the ECM results and iteration counts.
    for cores, ecm in ecm_predictions.items():
        pred = compute_pmodel_kernel_pred(iteration_count, ecm, machine.elements_per_cacheline)
        # Store pmodel predictions for the computation of the variant's kernel runtime prediction.
        if cores not in pmodel_kernel_predictions:
            pmodel_kernel_predictions[cores] = list()
        pmodel_kernel_predictions[cores].append(pred)
    return pmodel_kernel_predictions, ecm_predictions


def train_kernel_runtimes_ode(db_session: Session, machine: MachineState, templates: List[KernelTemplate],
                              methods: List[ODEMethod], ivps: List[IVP]):
    """Train database with kernel runtimes.

    Compute the kernel runtimes for all configuration of kernels, IVP, ODE method and machine state.

    Parameters:
    -----------
    db_session: sqlalchemy.orm.session.Session
        Used database session.
    machine: MachineState
        Used MachineState.
    templates: List of KernelTemplate
        Kernel templates to be trained.
    methods: List of ODEMethod.
        Used ODE methods.
    ivps: List of IVP
        Used IVPs.

    Returns:
    --------
    -
    """
    config: Config = offsite.config.offsiteConfig
    # Benchmark kernel runtimes for ...
    accumulated_data = list()
    # ... all ODE methods.
    for method in methods:
        # ... all kernel templates.
        for template in templates:
            # ... all kernels.
            for kernel in template.variants:
                tool = template.modelTool
                if template.isIVPdependent:
                    # ... all IVPs.
                    # If the template contains IVP calls, a separate kernel prediction has to be computed for each
                    # evaluated IVP.
                    for ivp in ivps:
                        if tool == ModelToolType.YASKSITE and not ivp.characteristic.isStencil:
                            print('\nSkipping kernel template \'{}\'! Unable to '.format(template.name), end='')
                            print('benchmark using YASKSITE: {} is not a stencil!'.format(ivp.name))
                            continue
                        if tool is not ivp.modelTool:
                            print('\nSkipping kernel template \'{}\'! Unable to predict '.format(template.name), end='')
                            print('using model tool  \'{}\': IVP \'{}\' requires model tool \'{}\'!'.format(
                                tool.value, ivp.name, ivp.modelTool.value))
                            continue
                        # If 'no_update_mode' is set, check if ECM predictions for this kernel are already in the
                        # database for the given run configuration.
                        if config.args.update is True:
                            if KernelRecord.contains(
                                    db_session, kernel.db_id, machine.db_id, machine.compiler.db_id, method.db_id,
                                    ivp.db_id, machine.clock, machine.coresPerSocket, config.args.mode.value,
                                    config.args.ode_size, [SampleType.BENCH_INNER, SampleType.BENCH_BORDER]):
                                continue
                        # Benchmark kernel.
                        data = bench_kernel_runtime(kernel, machine, method, ivp)
                        accumulated_data.append(data)
                else:
                    # If 'no_update_mode' is set, check if ECM predictions for this kernel are already in the database
                    # for the given run configuration.
                    if config.args.update is True:
                        if KernelRecord.contains(
                                db_session, kernel.db_id, machine.db_id, machine.compiler.db_id, method.db_id, -1,
                                machine.clock, machine.coresPerSocket, config.args.mode.value, config.args.ode_size,
                                [SampleType.BENCH_INNER, SampleType.BENCH_BORDER]):
                            continue
                    # Benchmark kernel.
                    data = bench_kernel_runtime(kernel, machine, method, None)
                    accumulated_data.append(data)
    # Before training the database, reduce the total number of intervals by:
    # * combining adjacent intervals giving the same prediction,
    # * interpolating border intervals.
    for kernel_db_id, run_config_db, kernel_predictions in accumulated_data:
        for cores, records in kernel_predictions.items():
            records = fuse_equal_records(records)
            KernelRecord.update(
                db_session, kernel_db_id, run_config_db[0], run_config_db[1], run_config_db[2], run_config_db[3],
                cores, machine.clock, records, config.args.mode.value)


def train_kernel_runtimes_generic(db_session: Session, machine: MachineState, templates: List[KernelTemplate]):
    """Train database with kernel runtimes.

    Compute the kernel runtimes for all configuration of kernels, IVP, ODE method and machine state.

    Parameters:
    -----------
    db_session: sqlalchemy.orm.session.Session
        Used database session.
    machine: MachineState
        Used MachineState.
    templates: List of KernelTemplate
        Kernel templates to be trained.

    Returns:
    --------
    -
    """
    config: Config = offsite.config.offsiteConfig
    # Benchmark kernel runtimes for ...
    accumulated_data = list()
    # ... all kernel templates.
    for template in templates:
        # ... all kernels.
        for kernel in template.variants:
            # If 'no_update_mode' is set, check if ECM predictions for this kernel are already in the database for the
            # given run configuration.
            if config.args.update is True:
                if KernelRecord.contains(
                        db_session, kernel.db_id, machine.db_id, machine.compiler.db_id, -100, -100, machine.clock,
                        machine.coresPerSocket, config.args.mode.value, config.args.ode_size,
                        [SampleType.BENCH_INNER, SampleType.BENCH_BORDER]):
                    continue
            # Benchmark kernel.
            data = bench_kernel_runtime(kernel, machine)
            accumulated_data.append(data)
    # Before training the database, reduce the total number of intervals by:
    # * combining adjacent intervals giving the same prediction,
    # * interpolating border intervals.
    for kernel_db_id, run_config_db, kernel_predictions in accumulated_data:
        for cores, records in kernel_predictions.items():
            records = fuse_equal_records(records)
            KernelRecord.update(db_session, kernel_db_id, run_config_db[0], run_config_db[1], -100, -100, cores,
                                machine.clock, records, config.args.mode.value)


def bench_kernel_runtime(
        kernel: Kernel, machine: MachineState, method: Optional[ODEMethod] = None, ivp: Optional[IVP] = None):
    """
    Benchmark the kernel runtime prediction for a given configuration of kernel, IVP, ODE method and machine state.

    Parameters:
    -----------
    kernel: Kernel
        Trained kernel.
    machine: MachineState
        Trained machine.
    method: ODE Method
        Trained ODE method.
    ivp: IVP
        Trained IVP.

    Returns:
    --------
    -
    """
    config: Config = offsite.config.offsiteConfig
    # MachineState frequency used.
    frequency = machine.clock
    # Max number of CPU cores used.
    max_cores = machine.coresPerSocket + 1
    # Generate pmodel code for all pmodel kernels of this kernel.
    kernel.generate_pmodel_code(method, ivp)
    # Initialize lists to store the benchmarked kernel runtimes for different numbers of cores.
    kernel_runtimes = {cores: list() for cores in range(1, max_cores)}
    # If no fixed ODE system size is given, deduce the set of significant sample intervals from the pmodel kernels'
    # working sets.
    if config.args.ode_size:
        ode_size = config.args.ode_size
        intervals = [SampleInterval(ode_size, ode_size, SampleType.BENCH_INNER, ode_size)]
    elif config.args.nrange:
        intervals = config.args.nrange
    else:
        intervals = kernel.deduce_relevant_samples(machine, method, ivp)
    for interval in intervals:
        # Determine the benchmarked runtimes for each pmodel using the designated tool.
        pmodel_runtimes = dict()
        # Benchmark pmodel kernels.
        for pmodel in kernel.pmodel_kernels:
            pmodel_runtimes = train_pmodel_runtime(pmodel, machine, interval, pmodel_runtimes, method, ivp)
            # Compute the kernel runtime prediction of this variant by combining its benchmarked pmodel kernel
            # predictions.
            for cores, predictions in pmodel_runtimes.items():
                pred = compute_kernel_pred(predictions, frequency)
                if cores not in kernel_runtimes:
                    kernel_runtimes[cores] = list()
                kernel_runtimes[cores].append((interval, pred))
    return kernel.db_id, (
        machine.db_id, machine.compiler.db_id, -1 if not method else method.db_id,
        - 1 if not ivp else ivp.db_id), kernel_runtimes


def train_pmodel_runtime(
        pmodel: PModelKernel, machine: MachineState, interval: SampleInterval, pmodel_runtimes: PModelPredictionDict,
        method: Optional[ODEMethod] = None, ivp: Optional[IVP] = None):
    """Train database with benchmarked pmodel kernel runtime prediction.

    Benchmark the pmodel kernel runtime prediction for a given configuration of kernel, IVP, ODE method and machine
    state.

    Parameters:
    -----------
    pmodel: PModel
        Trained pmodel kernel.
    machine: MachineState
        Trained machine.
    method: ODE method
        Trained ODE method.
    ivp: IVP
        Trained IVP.
    interval: SampleInterval
        Trained system size interval.
    pmodel_predictions: dict (key=int, value=str)
        Trained benchmarked pmodel kernel runtime predictions sorted by number of CPU cores used.

    Returns:
    --------
    dict (key=int, value=str)
        Trained pmodel kernel runtime predictions.
    """
    # Compute the number of iterations.
    iteration_count = pmodel.iteration_count(method, ivp)
    # Max number of cores used.
    max_cores = machine.coresPerSocket + 1
    # Compute ECM predictions.
    ecm_predictions = dict()
    for cores in range(1, max_cores):
        if isinstance(pmodel.code_path, dict):
            pred = 0.0
            for code_path in pmodel.code_path.values():
                # Execute kerncraft and parse its output.
                out = execute_kerncraft_bench_mode(code_path, machine, method, ivp, interval.sample, cores)
                pred += parse_kerncraft_output_bench_mode(out)
            ecm_predictions[cores] = pred
        else:
            if pmodel.kernel.template.modelTool == ModelToolType.KERNCRAFT:
                # Execute kerncraft and parse its output.
                out = execute_kerncraft_bench_mode(pmodel.code_path, machine, method, ivp, interval.sample, cores)
                ecm_predictions[cores] = parse_kerncraft_output_bench_mode(out)
            elif pmodel.kernel.template.modelTool == ModelToolType.YASKSITE:
                out = execute_yasksite_bench_mode(
                    pmodel.code_path, machine, method, ivp, interval.sample, pmodel.kernel.optimization_parameters)
                ecm_predictions = parse_yasksite_output(out)
            else:
                assert False
    # Compute pmodel kernel predictions using the ECM results and iteration counts.
    predictions = dict()
    for cores, ecm in ecm_predictions.items():
        pred = compute_pmodel_kernel_pred(iteration_count, ecm, machine.elements_per_cacheline)
        predictions[cores] = pred
        # Store pmodel predictions for the computation of the variant's kernel runtime prediction.
        if cores not in pmodel_runtimes:
            pmodel_runtimes[cores] = list()
        pmodel_runtimes[cores].append(pred)
    return pmodel_runtimes


def bench_rhs_kernel_runtime_ode(kernel: Kernel, machine: MachineState, method: ODEMethod, ivp: IVP):
    """
    Benchmark the kernel runtime prediction for a given configuration of kernel, IVP, ODE method and machine state.

    Parameters:
    -----------
    kernel: Kernel
        Trained kernel.
    machine: MachineState
        Trained machine.
    method: ODE Method
        Trained ODE method.
    ivp: IVP
        Trained IVP.

    Returns:
    --------
    -
    """
    config: Config = offsite.config.offsiteConfig
    # Max number of CPU cores used.
    max_cores = machine.coresPerSocket + 1
    # Generate instrumented kernel code.
    code: KernelBenchFiles = kernel.generate_rhs_bench_code(method, ivp)
    # Initialize lists to store the benchmarked kernel runtimes for different numbers of cores.
    kernel_runtimes = {cores: list() for cores in range(1, max_cores)}
    # If no fixed ODE system size is given, deduce the set of significant sample intervals from the pmodel kernels'
    # working sets.
    if config.args.ode_size:
        ode_size = config.args.ode_size
        intervals = [SampleInterval(ode_size, ode_size, SampleType.BENCH_INNER, ode_size)]
    elif config.args.nrange:
        intervals = config.args.nrange
    else:
        intervals = kernel.deduce_relevant_samples(machine, method, ivp)
    # Compile code.
    binary: Path = Path('tmp/a.out')
    # Construct compiler call.
    cmd = [machine.compiler.name]
    cmd.extend((flag for flag in machine.compiler.flags.split(' ')))
    cmd.extend((str(code.kernel_src), str(code.dummy_src), str(code.main_src)))
    cmd.append('-lm')
    cmd.append('-o{}'.format(binary))
    # Compile.
    run_process(cmd)
    # Run benchmark.
    for cores in range(1, max_cores):
        putenv('OMP_NUM_THREADS', str(cores))
        for interval in intervals:
            cmd = ['./{}'.format(binary), str(config.repetitions_communication_operations), str(cores),
                   str(interval.sample), str(method.stages)]
            output = run_process(cmd)
            # Process benchmark results.
            data = remove_outliers([float(x) for x in output.split('\n') if x != ''])
            kernel_runtimes[cores].append((interval, mean(data)))
    # We pass an empty dict here since we don't collect pmodel ECM records but use the same data formats.
    empty_dict = dict()
    return kernel.db_id, (machine.db_id, machine.compiler.db_id, method.db_id, - 1 if not ivp else ivp.db_id), \
           kernel_runtimes, empty_dict
