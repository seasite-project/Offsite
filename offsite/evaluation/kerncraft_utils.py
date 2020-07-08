"""@package kerncraft_util
Utility functions to use the kerncraft tool (https://github.com/RRZE-HPC/kerncraft).
"""

from pathlib import Path
from subprocess import run, PIPE, STDOUT, CalledProcessError
from sys import version_info
from typing import Dict, List

import offsite.config
from offsite.descriptions.ivp import IVP
from offsite.descriptions.machine import Machine
from offsite.descriptions.ode_method import ODEMethod
from offsite.evaluation.math_utils import eval_math_expr, corrector_steps, stages, ivp_system_size


def execute_kerncraft(kernel: Path, machine: Machine, model: str, defines: List[str]) -> str:
    """Run the kerncraft tool with a given kernel code and return the output.

    Parameters:
    -----------
    kernel: pathlib.Path
        Relative path to the kernel file executed.
    machine: Machine.
        Machine the kernel code is run on.
    model : str
        Name of the kerncraft mode used.
    defines: list of str
        Variable defines passed to kerncraft.

    Returns:
    --------
    str
        Output of the kerncraft tool.
    """
    config = offsite.config.offsiteConfig
    # Kerncraft options.
    cmd = ['kerncraft', '-p', model, '-m', str(machine.path), '-i', config.args.incore.value]
    # Defines.
    cmd.extend(defines)
    # Kernel name
    cmd.append(str(kernel))
    try:
        if version_info[1] > 5:
            return run(cmd, check=True, encoding='utf-8', stdout=PIPE, stderr=STDOUT).stdout
        return run(cmd, check=True, stdout=PIPE).stdout
    except CalledProcessError as error:
        if model == 'LC':
            # TODO
            return 'LC failed'
        else:
            raise RuntimeError('kerncraft failed: {}'.format(error))


def execute_kerncraft_ecm_mode(kernel: Path, machine: Machine, method: ODEMethod, ivp: IVP, iterations: int) -> str:
    """Run the kerncraft tool in ECM mode with a given kernel code and return the output.

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
    iterations : int
        Number of iterations executed by the kernel.

    Returns:
    --------
    str
        Output of the kerncraft tool.
    """
    # Defines.
    defines = ['-D', 'n', str(iterations)]
    if method:
        defines.extend(['-D', 's', str(method.stages), '-D', 'm', str(method.correctorSteps)])
    if ivp:
        grid_size = int(eval_math_expr(
            ivp.gridSize, [corrector_steps(method), stages(method), ivp_system_size(iterations)]))
        defines.extend(['-D', 'g', str(grid_size)])
    # Run kerncraft.
    return execute_kerncraft(kernel, machine, 'ECM', defines)


def execute_kerncraft_bench_mode(kernel: Path, machine: Machine, method: ODEMethod, ivp: IVP, iterations: int,
                                 cores: int) -> str:
    """Run the kerncraft tool in benchmark mode with a given kernel code and return the output.

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
    iterations : int
        Number of iterations executed by the kernel.
    cores : int
        Number of cores used.

    Returns:
    --------
    str
        Output of the kerncraft tool.
    """
    # Defines.
    defines = ['-D', 'n', str(iterations), '--ignore-warnings', '--cores', str(cores), '--no-phenoecm']
    if method:
        defines.extend(['-D', 's', str(method.stages), '-D', 'm', str(method.correctorSteps)])
    if ivp:
        grid_size = int(eval_math_expr(
            ivp.gridSize, [corrector_steps(method), stages(method), ivp_system_size(iterations)]))
        defines.extend(['-D', 'g', str(grid_size)])
    # Run kerncraft.
    return execute_kerncraft(kernel, machine, 'Benchmark', defines)


def execute_kerncraft_lc_mode(kernel: Path, machine: Machine, method: ODEMethod, ivp: IVP) -> str:
    """Run the kerncraft tool in LC (layer condition) mode with a given kernel code and return the output.

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

    Returns:
    --------
    str
        Output of the kerncraft tool.
    """

    dummy_n = 1000
    # Defines.
    defines = ['-D', 'n', str(dummy_n)]
    if method:
        defines.extend(['-D', 's', str(method.stages), '-D', 'm', str(method.correctorSteps)])
    if ivp:
        grid_size = int(eval_math_expr(
            ivp.gridSize, [corrector_steps(method), stages(method), ivp_system_size(dummy_n)]))
        defines.extend(['-D', 'g', str(grid_size)])
    # Run kerncraft.
    return execute_kerncraft(kernel, machine, 'LC', defines)


def check_kerncraft_version():
    """Check kerncraft version.

    Parameters:
    -----------
    -

    Returns:
    --------
    -
    """
    # Check kerncraft version.
    if version_info[1] > 5:
        kerncraft_version = run(['kerncraft', '--version'], check=True, encoding='utf-8', stdout=PIPE).stdout
    else:
        kerncraft_version = run(['kerncraft', '--version'], check=True, stderr=PIPE).stderr
        kerncraft_version = kerncraft_version.decode("utf-8")
    # Split version number from kerncraft's version string.
    kerncraft_version_number = kerncraft_version.split(' ')[1].strip().split('.')
    # Split version number.
    kerncraft_version_major = int(kerncraft_version_number[0])
    kerncraft_version_minor = int(kerncraft_version_number[1])
    kerncraft_version_micro = int(kerncraft_version_number[2])
    # Check kerncraft version.
    required_kerncraft_version_major = 0
    required_kerncraft_version_minor = 8
    required_kerncraft_version_micro = 5
    if required_kerncraft_version_major > kerncraft_version_major \
            or required_kerncraft_version_minor > kerncraft_version_minor \
            or required_kerncraft_version_micro > kerncraft_version_micro:
        raise RuntimeWarning(
            'Warning: offsite_tune might be unable to parse output of this kerncraft version. '
            'Supported versions are {}.{}.{} or higher'.format(required_kerncraft_version_major,
                                                               required_kerncraft_version_minor,
                                                               required_kerncraft_version_micro))


def parse_kerncraft_output_ecm_mode(output: str) -> Dict[int, float]:
    """Parse output of kerncraft's ECM mode and return the ECM results.

    Parameters:
    -----------
    output: str
        Kerncraft output.

    Returns:
    --------
    dict(int, float)
        ECM predictions obtained for different core counts (key: number of cores: value: prediction).
    """
    # Check kerncraft version.
    check_kerncraft_version()
    # Parse kerncraft output.
    cores = None
    predictions = None
    for line in output.splitlines(True):
        line = line.strip()
        if line.startswith('cores'):
            line = line.split('||')[1]
            cores = [int(x) for x in line.split('|')]
        elif line.startswith('perf. (cy/CL)'):
            line = line.split('||')[1]
            predictions = [float(x) for x in line.split('|')]
    if not cores:
        raise RuntimeError('Unable to parse kerncraft output: {}'.format(output))
    if not predictions:
        raise RuntimeError('Unable to parse kerncraft output: {}'.format(output))
    if len(cores) != len(predictions):
        raise RuntimeError('Unable to parse kerncraft output: {}'.format(output))
    results = dict(zip(cores, predictions))
    return results


def parse_kerncraft_output_bench_mode(output: str) -> Dict[int, float]:
    """Parse output of kerncraft's benchmark mode and return the ECM results.

    Parameters:
    -----------
    output: str
        Kerncraft output.

    Returns:
    --------
    dict(int, float)
        ECM predictions obtained for different core counts (key: number of cores: value: prediction).
    """
    # Check kerncraft version.
    check_kerncraft_version()
    # Parse kerncraft output.
    ecm = None
    for line in output.splitlines(True):
        line = line.strip()
        if line.startswith('Runtime (per cacheline update'):
            line = line.split(':')[1].strip()
            ecm = float(line.split('cy')[0].strip())
    if ecm is None:
        raise RuntimeError('Unable to parse kerncraft output: {}'.format(output))
    return ecm


def parse_kerncraft_output_lc_mode(output: str) -> Dict[int, float]:
    """Parse output of kerncraft's lc (layer condition) mode and return the suggested values.

    Parameters:
    -----------
    output: str
        Kerncraft output.

    Returns:
    --------
    TODO
        TODO
    """
    # Check kerncraft version.
    check_kerncraft_version()
    # Parse kerncraft output.
    lcs = dict()
    cur_key = None
    for line in output.splitlines(True):
        line = line.strip()
        if line.startswith('Layer conditions for L1'):
            cur_key = 'L1'
            lcs[cur_key] = list()
        elif line.startswith('Layer conditions for L2'):
            cur_key = 'L2'
            lcs[cur_key] = list()
        elif line.startswith('Layer conditions for L3'):
            cur_key = 'L3'
            lcs[cur_key] = list()
        elif cur_key is None or not line or line.startswith('condition') or line.startswith('True'):
            continue
        else:
            if '<=' in line:
                split = line.split('<=')
                lhs = split[0].strip()
                rhs = eval_math_expr(split[1].strip().split(' ')[0], cast_to=int)
                lcs[cur_key].append((lhs, rhs))
            else:
                split = line.split('<')
                lhs = split[0].strip()
                rhs = eval_math_expr(split[1].strip().split(' ')[0], cast_to=int)
                lcs[cur_key].append((lhs, rhs))
    return lcs
