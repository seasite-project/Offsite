"""@package apps.offsite_bench_openmp
Main script of the offsite_bench_openmp application.

@author: Johannes Seiferth
"""

from argparse import ArgumentParser, Namespace
from pathlib2 import Path

import offsite.config
from offsite import __version__
from offsite.config import init_config, Config
from offsite.descriptions.machine import MachineState
from offsite.train.communication.benchmark import BenchmarkType
from offsite.train.communication.openmp.omp_barrier import OmpBarrierBenchmark
from offsite.util.file_extensions import __config_ext__


def parse_program_args_app_bench_openmp() -> Namespace:
    """
    Parse the available program arguments of the offsite_bench_openmp application.

    Parameters:
    -----------
    -

    Returns:
    --------
    argparse.Namespace
        Parsed program arguments.
    """
    # Create argument parser object.
    parser = ArgumentParser(description='Run OpenMP communication cost benchmarks on this system.')
    # Available options.
    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(__version__),
                        help='Print program version and exit.')
    parser.add_argument('--verbose', action='store_true', default=False, help='Print further information on this run.')
    parser.add_argument('--benchmark', action='store', required=True, type=BenchmarkType,
                        help='Name of the benchmark executed.\nAvailable benchmarks:\n * omp_barrier')
    parser.add_argument('--machine', action='store', required=True, type=Path,
                        help='Path to YAML machine state description (.yaml) file.')
    parser.add_argument('--compiler', action='store', required=True, help='Name of the compiler.')
    parser.add_argument('--config', action='store', type=Path,
                        help='Tweak code generation process by passing a custom configuration ({}) file'.format(
                            __config_ext__))
    # Parse program arguments.
    return parser.parse_args()


def run_benchmark():
    """Run the OpenMP communication cost benchmarks.

    Parameters:
    -----------
    -

    Returns:
    -
    """
    config: Config = offsite.config.offsiteConfig
    # Parse passed machine state file.
    machine = MachineState.from_yaml(config.args.machine, config.args.compiler)
    # Run benchmark.
    if config.args.benchmark == BenchmarkType.OMP_BARRIER:
        # Run benchmark and write results to file.
        params = {'reps': config.repetitions_communication_operations}
        OmpBarrierBenchmark().run(machine, params, save_in_db=False, write_to_file=True)
    else:
        raise RuntimeWarning('Unknown benchmark {}!'.format(config.args.benchmark))


def run():
    """
    Run command line interface.

    Parameters:
    -----------
    -

    Returns:
    -
    """
    # Create parser and parse arguments.
    args = parse_program_args_app_bench_openmp()
    # Create custom or default configuration.
    init_config(args)
    # Run bench app.
    run_benchmark()
