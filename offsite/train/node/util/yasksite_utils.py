"""@package train.node.util.yasksite_util
Utility functions to use the yasksite tool.

@author: Johannes Seiferth
"""

from os import environ
from subprocess import run, PIPE, CalledProcessError
from typing import Dict, List, Optional

from pathlib2 import Path

import offsite.config
from offsite.config import Config
from offsite.descriptions.machine import MachineState
from offsite.descriptions.ode import IVP, ODEMethod, ivp_system_size, corrector_steps, stages
from offsite.util.math_utils import eval_math_expr

StringDict = Dict[str, str]
FloatDict = Dict[int, float]


def execute_yasksite(kernel: Path, machine: MachineState, model: str, size_str: str, timesteps: str, radius: str,
                     optimization_parameters: StringDict, stencil_file_path: str) -> str:
    """Run the yasksite tool with a given kernel code and return the output.

    Parameters:
    -----------
    kernel: pathlib.Path
        Relative path to the kernel file executed.
    machine: MachineState.
        Machine the kernel code is run on.
    model: str
        Name of the yasksite mode used.

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


def execute_yasksite_ecm_mode(kernel: Path, machine: MachineState, method: Optional[ODEMethod], ivp: Optional[IVP],
                              iterations: int, optimization_parameters: StringDict) -> str:
    """Run the yasksite tool in ECM mode with a given kernel code and return the output.

    Parameters:
    -----------
    kernel: pathlib.Path
        Relative path to the kernel file executed.
    machine: MachineState
        Trained machine.
    method: ODE method
        Trained ODE method.
    ivp: IVP
        Trained IVP.
    iterations: integer
        Number of iterations executed by the kernel.
    optimization_parameters: dict (key=str, value=str)
        Applied YaskSite optimization parameters.

    Returns:
    --------
    str
        Output of the yasksite tool.
    """
    config: Config = offsite.config.offsiteConfig
    size = ""
    timesteps = 1
    radius = str(-1)
    stencil_path = ''
    if method is not None:
        size = str(iterations) + ':1'
        timesteps = str(method.correctorSteps)
    if ivp is not None:
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
        if 'YASKSITE_STENCIL_DIR' not in ivp.code_stencil_path:
            raise RuntimeError('Stencil file of IVP \'{}\' missing.'.format(ivp.name))
        resolved_path = Path(ivp.code_stencil_path.replace('YASKSITE_STENCIL_DIR', config.yasksite_stencil_dir))
        # ... absolute path.
        stencil_path = resolved_path.resolve()
        # ... check if path exists.
        if not stencil_path.exists():
            raise RuntimeError('Stencil file of IVP \'{}\' not found: \'{}\''.format(ivp.name, stencil_path))
    return execute_yasksite(kernel, machine, 'ECM', size, timesteps, radius, optimization_parameters, stencil_path)


def execute_yasksite_bench_mode(kernel: Path, machine: MachineState, method: Optional[ODEMethod], ivp: Optional[IVP],
                                iterations: int, optimization_parameters: StringDict) -> str:
    """Run the yasksite tool in benchmark mode with a given kernel code and return the output.

    Parameters:
    -----------
    kernel: pathlib.Path
        Relative path to the kernel file executed.
    machine: MachineState
        Trained machine.
    method: ODE method
        Trained ODE method.
    ivp: IVP
        Trained IVP.
    iterations: integer
        Number of iterations executed by the kernel.
    optimization_parameters: dict (key=str, value=str)
        Applied YaskSite optimization parameters.

    Returns:
    --------
    str
        Output of the yasksite tool.
    """
    config: Config = offsite.config.offsiteConfig
    size = ""
    timesteps = 1
    radius = str(-1)
    stencil_path = ''
    if method is not None:
        size = str(iterations) + ':1'
        timesteps = str(method.correctorSteps)
    if ivp is not None:
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
        if 'YASKSITE_STENCIL_DIR' not in ivp.code_stencil_path:
            raise RuntimeError('Stencil file of IVP \'{}\' missing.'.format(ivp.name))
        resolved_path = Path(ivp.code_stencil_path.replace('YASKSITE_STENCIL_DIR', config.yasksite_stencil_dir))
        # ... absolute path.
        stencil_path = resolved_path.resolve()
        # ... check if path exists.
        if not stencil_path.exists():
            raise RuntimeError('Stencil file of IVP \'{}\' not found: \'{}\''.format(ivp.name, stencil_path))
    return execute_yasksite(kernel, machine, 'BENCH', size, timesteps, radius, optimization_parameters, stencil_path)


def parse_yasksite_output(output: str) -> FloatDict:
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
    cores: List[int] = list()
    predictions: List[float] = list()
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
    results: FloatDict = dict(zip(cores, predictions))
    return results
