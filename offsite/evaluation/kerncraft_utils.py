"""@package kerncraft_util
Utility functions to use the kerncraft tool (https://github.com/RRZE-HPC/kerncraft).
"""

from pathlib import Path
from subprocess import run, PIPE, STDOUT, CalledProcessError
from sys import version_info
from typing import Dict, List

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
    # Kerncraft options.
    cmd = ['kerncraft', '-p', model, '-m', str(machine.path)]
    # Defines.
    cmd.extend(defines)
    # Kernel name
    cmd.append(str(kernel))
    try:
        if version_info[1] > 5:
            return run(cmd, check=True, encoding='utf-8', stdout=PIPE, stderr=STDOUT).stdout
        return run(cmd, check=True, stdout=PIPE).stdout
    except CalledProcessError as error:
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
    required_kerncraft_version_micro = 3
    if required_kerncraft_version_major > kerncraft_version_major \
            or required_kerncraft_version_minor > kerncraft_version_minor \
            or required_kerncraft_version_micro > kerncraft_version_micro:
        raise RuntimeWarning('Warning: offsite_tune might be unable to parse output of this kerncraft version. '
                             'Supported versions are {}.{}.{} or higher'.format(
            required_kerncraft_version_major, required_kerncraft_version_minor, required_kerncraft_version_micro))


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
