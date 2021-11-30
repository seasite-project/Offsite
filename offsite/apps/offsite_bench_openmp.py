"""@package apps.offsite_bench_openmp
Main script of the offsite_bench_openmp application.
"""

from argparse import ArgumentParser, Namespace
from pathlib import Path

import offsite.config
from offsite import __version__
from offsite.config import BenchType, __config_ext__, init_config, Config
from offsite.descriptions.machine import MachineState
from offsite.train.communication.openmp.omp_barrier import OmpBarrierBenchmark


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
    parser.add_argument('--benchmark', action='store', required=True, type=BenchType,
                        help='Name of the benchmark executed.\nAvailable benchmarks:\n * omp_barrier')
    parser.add_argument('--machine', action='store', required=True, type=Path,
                        help='Path to YAML machine state description (.yaml) file.')
    parser.add_argument('--compiler', action='store', required=True, help='Name of the compiler.')
    parser.add_argument('--config', action='store', type=Path,
                        help='Tweak code generation process by passing a custom configuration ({}) file'.format(
                            __config_ext__))
    # Parse program arguments.
    return parser.parse_args()


def write_benchmark_file(bench: BenchType, data, machine: MachineState):
    """Write benchmark results to file.

    Parameters:
    -----------
    bench: BenchType
        Executed benchmark.
    data:
        Obtained benchmark results.
    machine: MachineState
        MachineState used to run the benchmark.

    Returns:
    -
    """
    bname = bench.value
    mname = machine.name
    cname = machine.compiler.name
    cversion = machine.compiler.version
    cflags = machine.compiler.flags
    path = Path('{}_{}_{}{}.bench'.format(bname, mname, cname, cversion))
    with path.open('w') as fhandle:
        # Write header.
        fhandle.write('benchmark: {}\n'.format(bname))
        fhandle.write('machine: {}\n'.format(mname))
        fhandle.write('compiler: {} {}\n'.format(cname, cversion))
        fhandle.write('flags: {}\n'.format(cflags))
        fhandle.write('frequency: {}\n'.format(machine.clock))
        # Write benchmark data.
        fhandle.write('data:\n')
        for rec in data:
            fhandle.write('- - {}\n'.format(rec.cores))
            fhandle.write('  - {}\n'.format(rec.data))


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
    if config.args.benchmark == BenchType.OMP_BARRIER:
        # Run benchmark.
        params = {'reps': config.repetitions_communication_operations}
        data = OmpBarrierBenchmark().run(machine, params, False)
        # Write benchmark data to file.
        write_benchmark_file(config.args.benchmark, data, machine)
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
