"""@package train_kernel
Functions to train the tuning database with kernel predictions.
"""

from multiprocessing import cpu_count, Pool
from traceback import print_exc
from typing import Dict, List, Tuple

from sqlalchemy.orm import Session

import offsite.config
from offsite.config import ModelToolType
from offsite.descriptions.ivp import IVP
from offsite.descriptions.kernel_template import Kernel, KernelTemplate, PModelKernel
from offsite.descriptions.machine import Machine
from offsite.descriptions.ode_method import ODEMethod
from offsite.evaluation.kerncraft_utils import execute_kerncraft_bench_mode, execute_kerncraft_ecm_mode, \
    parse_kerncraft_output_bench_mode, parse_kerncraft_output_ecm_mode
from offsite.evaluation.performance_model import compute_pmodel_kernel_pred, compute_kernel_runtime_pred, \
    KernelRecord, SampleInterval
from offsite.evaluation.yasksite_utils import execute_yasksite_bench_mode, execute_yasksite_ecm_mode, \
    parse_yasksite_output
from offsite.train.train_utils import reduce_records

IntervalRecordList = List[Tuple[SampleInterval, str]]


def train_kernel_predictions(db_session: Session, machine: Machine, templates: List[KernelTemplate],
                             methods: List[ODEMethod], ivps: List[IVP]):
    """Train database with kernel runtime predictions.

    Compute the kernel runtime predictions for all configuration of kernels, IVP, ODE method and available machine.

    Parameters:
    -----------
    db_session : sqlalchemy.orm.session.Session
        Used database session.
    machine : Machine
        Used Machine.
    templates : List of KernelTemplate
        Kernel templates to be trained.
    methods : List of ODEMethod.
        Used ODE methods.
    ivps : List of IVP
        Used IVPs.

    Returns:
    --------
    -
    """
    config = offsite.config.offsiteConfig
    data_kc = list()
    data_ys = list()
    errors = list()
    # Initialize worker thread pools.
    kernel_count = sum((len(t.variants) for t in templates))
    if not config.args.tool:
        # Determine number of workers.
        # One worker handles YASKSITE tasks ...
        num_workers_ys = 1
        # ... while remaining workers handle KERNCRAFT task.
        num_workers_kc = min(cpu_count(), kernel_count)
        num_workers_kc = max(num_workers_kc - 1, 1)
        # Spawn worker pools.
        pool_kc = Pool(num_workers_kc)
        pool_ys = Pool(num_workers_ys)
    elif config.args.tool == ModelToolType.KERNCRAFT:
        # Determine number of workers.
        num_workers_kc = min(cpu_count(), kernel_count)
        # Spawn worker pool.
        pool_kc = Pool(num_workers_kc)
    elif config.args.tool == ModelToolType.YASKSITE:
        num_workers_ys = 1
        # Spawn worker pool.
        pool_ys = Pool(num_workers_ys)
    else:
        assert False
    # Compute kernel runtime predictions for ...
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
                                    config.args.ode_size):
                                continue
                        # Push kernel prediction computation task to thread pool.
                        if tool == ModelToolType.KERNCRAFT:
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
                                machine.clock, machine.coresPerSocket, config.args.mode.value, config.args.ode_size):
                            continue
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
    if not config.args.tool:
        # Wait for threads.
        pool_kc.close()
        pool_ys.close()
        pool_kc.join()
        pool_ys.join()
        # Collect results.
        accumulated_data = data_kc + data_ys
    elif config.args.tool == ModelToolType.KERNCRAFT:
        # Wait for threads.
        pool_kc.close()
        pool_kc.join()
        # Collect results.
        accumulated_data = data_kc
    elif config.args.tool == ModelToolType.YASKSITE:
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
    for kernel_db_id, run_config_db, kernel_predictions in accumulated_data:
        for cores, records in kernel_predictions.items():
            records = reduce_records(records)
            KernelRecord.update(
                db_session, kernel_db_id, run_config_db[0], run_config_db[1], run_config_db[2], run_config_db[3],
                cores, machine.clock, records, config.args.mode.value)


def compute_kernel_prediction(kernel: Kernel, machine: Machine, method: ODEMethod, ivp: IVP):
    """Compute the kernel runtime prediction for a given configuration of kernel, IVP, ODE method and available machine.

    Parameters:
    -----------
    kernel : Kernel
        Trained kernel.
    machine : Machine
        Trained machine.
    method : ODE Method
        Trained ODE method.
    ivp : IVP
        Trained IVP.

    Returns:
    --------
    -
    """
    config = offsite.config.offsiteConfig
    # Machine frequency used.
    frequency = machine.clock
    try:
        # Generate pmodel code for all pmodel kernels of this kernel.
        kernel.generate_pmodel_code(method, ivp)
        # Initialize lists to store the kernel predictions for different numbers of cores.
        kernel_predictions = dict()
        # If no fixed ODE system size is given, deduce the set of significant sample intervals from the pmodel kernels'
        # working sets.
        if config.args.ode_size:
            ode_size = config.args.ode_size
            intervals = [SampleInterval(ode_size, ode_size, ode_size)]
        else:
            intervals = kernel.deduce_relevant_samples(machine, method, ivp)
        # TODO check if ODE size is applicable
        for interval in intervals:
            # Determine the ECM result for each pmodel kernel using the kerncraft tool.
            pmodel_predictions = dict()
            for pmodel in kernel.pmodel_kernels:
                pmodel_predictions = train_pmodel_prediction(pmodel, machine, method, ivp, interval, pmodel_predictions)
            # Compute the kernel runtime prediction of this variant by combining its pmodel kernel predictions obtained.
            for cores, predictions in pmodel_predictions.items():
                pred = compute_kernel_runtime_pred(predictions, frequency)
                if cores not in kernel_predictions:
                    kernel_predictions[cores] = list()
                kernel_predictions[cores].append((interval, pred))
        return kernel.db_id, (machine.db_id, machine.compiler.db_id, method.db_id, - 1 if not ivp else ivp.db_id), \
               kernel_predictions
    except Exception as e:
        # Print stack trace of the executing worker thread.
        print_exc()
        print('')
        raise e


def train_pmodel_prediction(pmodel: PModelKernel, machine: Machine, method: ODEMethod, ivp: IVP,
                            interval: SampleInterval, pmodel_predictions: Dict[int, str]):
    """Train database with pmodel kernel runtime prediction.

    Compute the pmodel kernel runtime prediction for a given configuration of kernel, IVP, ODE method and available
    machine.

    Parameters:
    -----------
    pmodel : PModel
        Trained pmodel kernel.
    machine: Machine
        Trained machine.
    method: ODE method
        Trained ODE method.
    ivp: IVP
        Trained IVP.
    interval: SampleInterval
        Trained system size interval.
    pmodel_predictions : dict (key=int, value=str)
        Trained pmodel kernel runtime predictions sorted by number of CPU cores used.

    Returns:
    --------
    dict (key=int, value=str)
        Trained pmodel kernel runtime predictions.
    """
    # Compute the number of iterations.
    iteration_count = pmodel.iteration_count(method, ivp)
    if isinstance(pmodel.code_path, dict):
        ecm_predictions = dict()
        for code_path in pmodel.code_path.values():
            # Execute kerncraft and parse its output.
            out = execute_kerncraft_ecm_mode(code_path, machine, method, ivp, interval.sample)
            preds = parse_kerncraft_output_ecm_mode(out)
            ecm_predictions = {k: ecm_predictions.get(k, 0) + preds.get(k, 0)
                               for k in set(ecm_predictions) | set(preds)}
    else:
        if pmodel.kernel.template.modelTool == ModelToolType.KERNCRAFT:
            # Execute kerncraft and parse its output.
            out = execute_kerncraft_ecm_mode(
                pmodel.code_path, machine, method, ivp, interval.sample)
            ecm_predictions = parse_kerncraft_output_ecm_mode(out)
        elif pmodel.kernel.template.modelTool == ModelToolType.YASKSITE:
            out = execute_yasksite_ecm_mode(
                pmodel.code_path, machine, method, ivp, interval.sample, pmodel.kernel.optimization_parameters)
            ecm_predictions = parse_yasksite_output(out)
        else:
            assert False
    # Compute pmodel kernel predictions using the ECM results and iteration counts.
    predictions = dict()
    for cores, ecm in ecm_predictions.items():
        pred = compute_pmodel_kernel_pred(iteration_count, machine, ecm)
        predictions[cores] = pred
        # Store pmodel predictions for the computation of the variant's kernel runtime prediction.
        if cores not in pmodel_predictions:
            pmodel_predictions[cores] = list()
        pmodel_predictions[cores].append(pred)
    return pmodel_predictions


def train_kernel_runtimes(db_session: Session, machine: Machine, templates: List[KernelTemplate],
                          methods: List[ODEMethod], ivps: List[IVP]):
    """Train database with kernel runtimes.

    Compute the kernel runtimes for all configuration of kernels, IVP, ODE method and available machine.

    Parameters:
    -----------
    db_session : sqlalchemy.orm.session.Session
        Used database session.
    machine : Machine
        Used Machine.
    templates : List of KernelTemplate
        Kernel templates to be trained.
    methods : List of ODEMethod.
        Used ODE methods.
    ivps : List of IVP
        Used IVPs.

    Returns:
    --------
    -
    """
    config = offsite.config.offsiteConfig
    # TODO: Check if CPU frequency is fixed to the value specified in the machine description.
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
                                    config.args.ode_size):
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
                                machine.clock, machine.coresPerSocket, config.args.mode.value, config.args.ode_size):
                            continue
                    # Benchmark kernel.
                    data = bench_kernel_runtime(kernel, machine, method, None)
                    accumulated_data.append(data)
    # Before training the database, reduce the total number of intervals by:
    # * combining adjacent intervals giving the same prediction,
    # * interpolating border intervals.
    for kernel_db_id, run_config_db, kernel_predictions in accumulated_data:
        for cores, records in kernel_predictions.items():
            records = reduce_records(records)
            KernelRecord.update(
                db_session, kernel_db_id, run_config_db[0], run_config_db[1], run_config_db[2], run_config_db[3],
                cores, machine.clock, records, config.args.mode.value)


def bench_kernel_runtime(kernel: Kernel, machine: Machine, method: ODEMethod, ivp: IVP):
    """
    Benchmark the kernel runtime prediction for a given configuration of kernel, IVP, ODE method and available machine.

    Parameters:
    -----------
    kernel : Kernel
        Trained kernel.
    machine : Machine
        Trained machine.
    method : ODE Method
        Trained ODE method.
    ivp : IVP
        Trained IVP.

    Returns:
    --------
    -
    """
    config = offsite.config.offsiteConfig
    # Machine frequency used.
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
        intervals = [SampleInterval(ode_size, ode_size, ode_size)]
    else:
        intervals = kernel.deduce_relevant_samples(machine, method, ivp)
    # TODO check if ODE size is applicable
    for interval in intervals:
        # Determine the benchmarked runtimes for each pmodel using the designated tool.
        pmodel_runtimes = dict()
        # Benchmark pmodel kernels.
        for pmodel in kernel.pmodel_kernels:
            pmodel_runtimes = train_pmodel_runtime(pmodel, machine, method, ivp, interval, pmodel_runtimes)
            # Compute the kernel runtime prediction of this variant by combining its benchmarked pmodel kernel
            # predictions.
            for cores, predictions in pmodel_runtimes.items():
                pred = compute_kernel_runtime_pred(predictions, frequency)
                if cores not in kernel_runtimes:
                    kernel_runtimes[cores] = list()
                kernel_runtimes[cores].append((interval, pred))
    return kernel.db_id, (
        machine.db_id, machine.compiler.db_id, method.db_id, - 1 if not ivp else ivp.db_id), kernel_runtimes


def train_pmodel_runtime(pmodel: PModelKernel, machine: Machine, method: ODEMethod, ivp: IVP, interval: SampleInterval,
                         pmodel_runtimes: Dict[int, str]):
    """Train database with benchmarked pmodel kernel runtime prediction.

    Benchmark the pmodel kernel runtime prediction for a given configuration of kernel, IVP, ODE method and available
    machine.

    Parameters:
    -----------
    pmodel : PModel
        Trained pmodel kernel.
    machine: Machine
        Trained machine.
    method: ODE method
        Trained ODE method.
    ivp: IVP
        Trained IVP.
    interval: SampleInterval
        Trained system size interval.
    pmodel_predictions : dict (key=int, value=str)
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
        pred = compute_pmodel_kernel_pred(iteration_count, machine, ecm)
        predictions[cores] = pred
        # Store pmodel predictions for the computation of the variant's kernel runtime prediction.
        if cores not in pmodel_runtimes:
            pmodel_runtimes[cores] = list()
        pmodel_runtimes[cores].append(pred)
    return pmodel_runtimes
