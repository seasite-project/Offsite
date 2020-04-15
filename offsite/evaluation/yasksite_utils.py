"""@package yasksite_util
Utility functions to use the yasksite tool.
"""

from os import environ
from pathlib import Path
from subprocess import run, PIPE, CalledProcessError
from typing import Dict

import offsite.config
from offsite.descriptions.ivp import IVP
from offsite.descriptions.machine import Machine
from offsite.descriptions.ode_method import ODEMethod
from offsite.evaluation.math_utils import eval_math_expr, corrector_steps, stages, ivp_system_size


def execute_yasksite(kernel: Path, machine: Machine, model: str, size_str: str, timesteps: str,
                     radius: str, optimization_parameters: Dict[str, str], stencil_file_path: str) -> str:
    """Run the yasksite tool with a given kernel code and return the output.

    Parameters:
    -----------
    kernel: pathlib.Path
        Relative path to the kernel file executed.
    machine: Machine.
        Machine the kernel code is run on.
    model : str
        Name of the yasksite mode used.
    size_str : str
        TODO DOC
    timesteps : str
        DOC
    radius : str
        DOC
    optimization_parameters : dict (key=str, value=str)
        DOC
    stencil_file_path : str
        DOC

    Returns:
    --------
    str
        Output of the yasksite tool.
    """
    # Yasksite options.
    ys_exec = 'ys_offtune'
    # Check passed optimization parameters.
    folding = '' if 'ys_folding' not in optimization_parameters else optimization_parameters['ys_folding']
    ys_folding = '' if not folding else '-f {}'.format(folding)
    ys_blocking = 'plain' if 'ys_blocking' not in optimization_parameters else optimization_parameters['ys_blocking']
    # Set required environment variables.
    env_var = {'OMP_NUM_THREADS': str(machine.coresPerSocket), 'OMP_PLACES': 'cores', 'OMP_PROC_BIND': 'close'}
    environ.update(env_var)
    taskset_str = '0-' + str(machine.coresPerSocket - 1)
    # Execute yasksite.
    cmd = ['taskset', '-c', taskset_str, ys_exec, str(kernel), '-i', timesteps, '-s', size_str, '-o', ys_blocking,
           '-r', radius, ' -M', model, '-m', str(machine.path), '-b', str(stencil_file_path), ys_folding]
    try:
        return run(cmd, check=True, encoding='utf-8', stdout=PIPE).stdout
    except CalledProcessError as error:
        raise RuntimeError('yasksite failed: {}'.format(error))


def execute_yasksite_ecm_mode(kernel: Path, machine: Machine, method: ODEMethod, ivp: IVP, iterations: int,
                              optimization_parameters: Dict[str, str]) -> str:
    """Run the yasksite tool in ECM mode with a given kernel code and return the output.

    Parameters:
    -----------
    kernel : pathlib.Path
        Relative path to the kernel file executed.
    machine: Machine
        Trained machine.
    method: ODE method
        Trained ODE method.
    ivp: IVP
        Trained IVP.
    iterations : integer
        Number of iterations executed by the kernel.
    optimization_parameters : dict (key=str, value=str)
        Applied YaskSite optimization parameters.

    Returns:
    --------
    str
        Output of the yasksite tool.
    """
    config = offsite.config.offsiteConfig
    size = ""
    timesteps = 1
    radius = str(-1)
    if method:
        size = str(iterations) + ':1'
        timesteps = str(method.correctorSteps)
    if ivp:
        grid_size = int(eval_math_expr(ivp.gridSize,
                                       [corrector_steps(method), stages(method), ivp_system_size(iterations)]))
        if ivp.characteristic.stencil_dim == 2:
            size = str(grid_size) + ':' + str(grid_size)
        elif ivp.characteristic.stencil_dim == 3:
            size = str(grid_size) + ':' + str(grid_size) + ':' + str(grid_size)
        else:
            assert False
        radius = str(ivp.characteristic.stencil_radius)
        # Yasksite stencil file path.
        # ... substitute stencil directory variable.
        if 'YASKSITE_STENCIL_DIR' in ivp.code_stencil_path:
            resolved_path = Path(ivp.code_stencil_path.replace('YASKSITE_STENCIL_DIR', config.yasksite_stencil_dir))
        # ... absolute path.
        stencil_abs_path = resolved_path.resolve()
        # ... check if path exists.
        if not stencil_abs_path.exists():
            raise RuntimeError('Stencil file of IVP \'{}\' not found: \'{}\''.format(ivp.name, stencil_abs_path))
    return execute_yasksite(kernel, machine, 'ECM', size, timesteps, radius, optimization_parameters, stencil_abs_path)


def parse_yasksite_output(output: str) -> Dict[int, float]:
    """Parse yasksite's output and return the ECM results.

    Parameters:
    -----------
    output: str
        Yasksite output.

    Returns:
    --------
    dict(int, float)
        ECM predictions obtained for different core counts (key: number of cores: value: prediction).
    """
    # Parse yasksite output.
    cores = None
    predictions = None
    for line in output.splitlines(True):
        line = line.strip()
        if line.startswith('cores'):
            line = line.split('||')[1]
            cores = [int(x.replace(' ', '')) for x in line.split('|')]
        elif line.startswith('ECM perf. (cy/CL)') or line.startswith('BENCHMARK perf. (cy/CL)'):
            line = line.split('||')[1]
            predictions = [float(x.replace(' ', '')) for x in line.split('|')]
    if not cores:
        raise RuntimeError('Unable to parse yasksite output: {}'.format(output))
    if not predictions:
        raise RuntimeError('Unable to parse yasksite output: {}'.format(output))
    if len(cores) != len(predictions):
        raise RuntimeError('Unable to parse yasksite output: {}'.format(output))
    results = dict(zip(cores, predictions))
    return results
