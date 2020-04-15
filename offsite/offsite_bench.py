"""@package offsite_bench
Main script of the offsite_bench application.
"""

from argparse import ArgumentParser
from pathlib import Path

from offsite import __version__
from offsite.config import BenchType, Config, __config_ext__
from offsite.descriptions.machine import Machine
from offsite.evaluation.benchmark import OmpBarrierBenchmark


def parse_program_args_app_bench() -> 'argparse.Namespace':
    """
    Parse the available program arguments of the offsite_bench application.

    Parameters:
    -----------
    -

    Returns:
    --------
    argparse.Namespace
        Parsed program arguments.
    """
    # Create argument parser object.
    parser = ArgumentParser(description='Run communication cost benchmark on this system.')
    # Available options.
    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(__version__),
                        help='Print program version and exit.')
    parser.add_argument('--verbose', action='store_true', default=False, help='Print further information on this run.')
    parser.add_argument('--benchmark', action='store', required=True, type=BenchType,
                        help='Name of the benchmark executed.')
    parser.add_argument('--machine', action='store', required=True, type=Path,
                        help='Path to YAML machine description (.yaml) file.')
    parser.add_argument('--compiler', action='store', required=True, help='Name of the compiler.')
    parser.add_argument('--config', action='store', type=Path,
                        help='Tweak code generation process by passing a custom configuration ({}) file'.format(
                            __config_ext__))
    # Parse program arguments.
    return parser.parse_args()


def write_benchmark_file(bench: BenchType, data, machine: 'Machine'):
    """Write benchmark results to file.

    Parameters:
    -----------
    bench : BenchType
        Executed benchmark.
    data:
        Obtained benchmark results.
    machine:
        Machine used to run the benchmark.

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


def run_benchmark(config: 'Config'):
    """Run the communication cost benchmark.

    Parameters:
    -----------
    config: Config
        Program configuration.

    Returns:
    -
    """
    # Parse passed machine file.
    machine = Machine.from_yaml(config.args.machine, config.args.compiler)
    # TODO Pin CPU frequency / check if pinned
    if config.args.benchmark == BenchType.OMP_BARRIER:
        # Run benchmark.
        bench = OmpBarrierBenchmark()
        data = bench.run(machine, config.repetitions_communication_operations, False)
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
    args = parse_program_args_app_bench()
    # Parse custom configuration file if present else use default configuration.
    config = Config.from_file(args) if args.config else Config(args)
    # Run bench app.
    run_benchmark(config)
