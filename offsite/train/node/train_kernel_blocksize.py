"""@package train.node.train_kernel_blocksize
Functions to train the tuning database with kernel block size predictions.
"""

from datetime import timedelta
from multiprocessing import cpu_count, Pool
from pathlib import Path
from time import time
from traceback import print_exc
from typing import List, Optional

from sqlalchemy.orm import Session

import offsite.config
from offsite.config import Config, ModelToolType
from offsite.descriptions.impl.kernel_template import Kernel, KernelTemplate, PModelKernel
from offsite.descriptions.machine import MachineState
from offsite.descriptions.ode import IVP, ODEMethod, ivp_grid_size, corrector_steps, stages
from offsite.solver import SolverType
from offsite.train.node.util.kerncraft_utils import execute_kerncraft_lc_mode, parse_kerncraft_output_lc_mode
from offsite.util.math_utils import eval_math_expr, solve_equation


def train_kernel_blocksizes(
        db_session: Session, machine: MachineState, templates: List[KernelTemplate],
        methods: Optional[List[ODEMethod]] = None, ivps: Optional[List[IVP]] = None) -> List[KernelTemplate]:
    """Train database with kernel blocksize predictions.

    Parameters:
    -----------
    db_session: sqlalchemy.orm.session.Session
        Used database session.
    machine: MachineState
        Used machine.
    templates: List of KernelTemplate
        Kernel templates to be trained.
    methods: List of ODEMethod.
        Used ODE methods.
    ivps: List of IVP
        Used IVPs.

    Returns:
    --------
    List of KernelTemplate
        Trained Kernel Templates.
    """
    config: Config = offsite.config.offsiteConfig
    if config.args.verbose:
        print('  * Kernel block sizes...', end='', flush=True)
    start_time_block = time()
    if config.scenario.solver.type == SolverType.ODE:
        templates = train_kernel_blocksizes_ode(db_session, machine, templates, methods, ivps)
    elif config.scenario.solver.type == SolverType.GENERIC:
        templates = train_kernel_blocksizes_generic(db_session, machine, templates)
    templates = [template.to_database(db_session) for template in templates]
    stop_time_block = time()
    if config.args.verbose:
        print(' done. Duration {}.'.format(timedelta(seconds=round(stop_time_block - start_time_block, 0))))
    return templates


def train_kernel_blocksizes_ode(db_session: Session, machine: MachineState, templates: List[KernelTemplate],
                                methods: List[ODEMethod], ivps: List[IVP]):
    """Train database with kernel blocksize predictions.

    Compute the kernel block size predictions for all configuration of kernels, IVP, ODE method and machine state.

    Parameters:
    -----------
    db_session: sqlalchemy.orm.session.Session
        Used database session.
    machine: MachineState
        Used machine.
    templates: List of KernelTemplate
        Kernel templates to be trained.
    methods: List of ODEMethod.
        Used ODE methods.
    ivps: List of IVP
        Used IVPs.

    Returns:
    --------
    List of KernelTemplate
        Trained Kernel Templates.
    """
    data = list()
    errors = list()
    # Initialize worker thread pools.
    kernel_count = sum((len(t.variants) for t in templates))
    # Determine number of workers.
    num_workers = min(cpu_count() - 1, kernel_count)
    # Spawn worker pool.
    pool = Pool(num_workers)
    # Compute kernel block size predictions for ...
    # ... all ODE methods.
    for method in methods:
        # ... all kernel templates.
        for template in templates:
            # ... all kernels.
            for kernel in template.variants:
                tool = template.modelTool
                if tool != ModelToolType.KERNCRAFT:
                    continue
                # ... all IVPs.
                if template.isIVPdependent:
                    # If the template contains IVP calls, a separate kernel prediction has to be computed for each
                    # evaluated IVP.
                    for ivp in ivps:
                        if tool is not ivp.modelTool:
                            print('\nSkipping kernel template \'{}\'! Unable to predict '.format(template.name), end='')
                            print('block sizes using model tool \'{}\': IVP \'{}\' requires model tool \'{}\'!'.format(
                                tool.value, ivp.name, ivp.modelTool.value))
                            continue
                        pool.apply_async(train_kernel_blocksize_prediction, args=(kernel, machine, method, ivp),
                                         callback=data.append, error_callback=errors.append)
                else:
                    pool.apply_async(train_kernel_blocksize_prediction, args=(kernel, machine, method, None),
                                     callback=data.append, error_callback=errors.append)
    # Wait for threads.
    pool.close()
    pool.join()
    # Raise error if kernel prediction failed.
    if errors:
        db_session.rollback()
        raise RuntimeError('Failed to train kernel block size predictions: Error in worker threads.')
    # Convert to dictionary.
    data = dict(data)
    for template in templates:
        tool = template.modelTool
        for kernel in template.variants:
            if tool != ModelToolType.KERNCRAFT:
                continue
            kernel.optimization_parameters = dict()
            kernel.optimization_parameters['blocking_sizes'] = data[kernel.db_id]
    return templates


def train_kernel_blocksizes_generic(db_session: Session, machine: MachineState, templates: List[KernelTemplate]):
    """Train database with kernel blocksize predictions.

    Compute the kernel block size predictions for all configuration of kernels and machine state.

    Parameters:
    -----------
    db_session: sqlalchemy.orm.session.Session
        Used database session.
    machine: MachineState
        Used machine.
    templates: List of KernelTemplate
        Kernel templates to be trained.

    Returns:
    --------
    List of KernelTemplate
        Trained Kernel Templates.
    """
    data = list()
    errors = list()
    # Initialize worker thread pools.
    kernel_count = sum((len(t.variants) for t in templates))
    # Determine number of workers.
    num_workers = min(cpu_count(), kernel_count)
    # Spawn worker pool.
    pool = Pool(num_workers)
    # Compute kernel block size predictions for ...
    # ... all kernel templates.
    for template in templates:
        # ... all kernels.
        for kernel in template.variants:
            tool = template.modelTool
            if tool != ModelToolType.KERNCRAFT:
                continue
            pool.apply_async(train_kernel_blocksize_prediction, args=(kernel, machine), callback=data.append,
                             error_callback=errors.append)
    # Wait for threads.
    pool.close()
    pool.join()
    # Raise error if kernel prediction failed.
    if errors:
        db_session.rollback()
        raise RuntimeError('Failed to train kernel block size predictions: Error in worker threads.')
    # Convert to dictionary.
    data = dict(data)
    for template in templates:
        tool = template.modelTool
        for kernel in template.variants:
            if tool != ModelToolType.KERNCRAFT:
                continue
            kernel.optimization_parameters = dict()
            kernel.optimization_parameters['blocking_sizes'] = data[kernel.db_id]
    return templates


def train_kernel_blocksize_prediction(
        kernel: Kernel, machine: MachineState, method: Optional[ODEMethod] = None, ivp: Optional[IVP] = None):
    """
    Compute the kernel block size prediction for a given configuration of kernel, IVP, ODE method and machine state.

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
    try:
        # Generate pmodel code for all pmodel kernels of this kernel.
        kernel.generate_pmodel_code(method, ivp)
        # Apply layer condition analysis to determine suitable L1/L2/L3 cache blocking sizes.
        bsizes_per_cache_lvl = None
        for pmodel in kernel.pmodel_kernels:
            if isinstance(pmodel.code_path, dict):
                for pmodel_path in pmodel.code_path.values():
                    bsizes_per_cache_lvl = compute_kernel_blocking_sizes(pmodel_path, pmodel, machine, method, ivp)
            else:
                bsizes_per_cache_lvl = compute_kernel_blocking_sizes(pmodel.code_path, pmodel, machine, method, ivp)
        assert bsizes_per_cache_lvl is not None
        return kernel.db_id, bsizes_per_cache_lvl
    except Exception as e:
        # Print stack trace of the executing worker thread.
        print_exc()
        print('')
        raise e


def compute_kernel_blocking_sizes(path: Path, pmodel_kernel: PModelKernel, machine: MachineState,
                                  method: Optional[ODEMethod] = None, ivp: Optional[IVP] = None):
    """
    Compute blocking size for a given configuration of kernel, IVP, ODE method and machine state using kerncraft's
    layer condition analysis.

    Parameters:
    -----------
    path: pathlib.Path
        Relative path to the kernel file executed.
    kernel: Kernel
        Used kernel.
    machine: MachineState
        Used machine.
    method: ODE Method
        Used ODE method.
    ivp: IVP
        Used IVP.

    Returns:
    --------
    -
    """
    lc_safety_margin = offsite.config.offsiteConfig.layer_condition_safety_margin
    elem_per_cl = machine.elements_per_cacheline
    # Execute kerncraft in LC mode to try to obtain blocking size suggestions.
    out = execute_kerncraft_lc_mode(path, machine, method, ivp)
    # Parse kerncraft output.
    layer_conditions_per_cache_lvl = dict()
    if out == 'LC failed':
        # LC mode failed on this kernel. We instead use the working set to determine blocking sizes.
        # .. for L1 cache
        bsizes = pmodel_kernel.determine_max_n_of_working_sets_for_cache_lvl(method, machine.l1CacheElements)
        for idx, bs in enumerate(bsizes):
            bsizes[idx] = adjust_blocksize(bs, lc_safety_margin, elem_per_cl)
        layer_conditions_per_cache_lvl['L1'] = bsizes
        # .. for L2 cache
        bsizes = pmodel_kernel.determine_max_n_of_working_sets_for_cache_lvl(method, machine.l2CacheElements)
        for idx, bs in enumerate(bsizes):
            bsizes[idx] = adjust_blocksize(bs, lc_safety_margin, elem_per_cl)
        layer_conditions_per_cache_lvl['L2'] = bsizes
        # .. for L3 cache
        bsizes = pmodel_kernel.determine_max_n_of_working_sets_for_cache_lvl(method, machine.l3CacheElements)
        for idx, bs in enumerate(bsizes):
            bsizes[idx] = adjust_blocksize(bs, lc_safety_margin, elem_per_cl)
        layer_conditions_per_cache_lvl['L3'] = bsizes
    else:
        constants = list()
        if method is not None:
            constants.extend([corrector_steps(method), stages(method)])
        if ivp is not None:
            constants.append((ivp_grid_size(ivp.gridSize)))
        elem_per_cl = int(machine.elements_per_cacheline)
        # Parse kerncraft output.
        layer_conditions_per_cache_lvl = parse_kerncraft_output_lc_mode(out)
        for cache, lcs in layer_conditions_per_cache_lvl.items():
            for idx, lc in enumerate(lcs):
                bs = eval_math_expr(solve_equation(lc[0], lc[1], 'n', constants)[0], cast_to=int)
                layer_conditions_per_cache_lvl[cache][idx] = adjust_blocksize(bs, lc_safety_margin, elem_per_cl)
    return layer_conditions_per_cache_lvl


def adjust_blocksize(bs: int, safety_margin: float, elem_per_cl: int):
    """
    Adjust block size by applying a safety margin factor and ensuring that the block size is a multiple of the
    cacheline size.

    Parameters:
    -----------
    bs: int
        Block size.
    safety_margin: float
        Block size will be reduced by this factor.
    elem_per_cl: int
        Number of elements fitting into a single cacheline.

    Returns:
    --------
    float
        Adjusted block size.
    """
    # Apply safety margin.
    bs = int(bs / safety_margin)
    # Adjust block size to be a multiple of the number of elements per cacheline.
    return int(bs / elem_per_cl) * elem_per_cl
