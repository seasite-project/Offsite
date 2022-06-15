"""@package apps.offsite_bench_impl
Main script of the offsite_bench_impl application.

@author: Johannes Seiferth
"""

from argparse import ArgumentParser, Namespace
from pathlib2 import Path
from sqlalchemy.orm import Session

import offsite.config
from offsite import __version__
from offsite.config import Config, init_config
from offsite.database import close, open_db
from offsite.database.db_mapping import mapping
from offsite.descriptions.parser_util import parse_verify_yaml_desc, print_yaml_desc
from offsite.solver import SolverType
from offsite.train.node.train_impl import train_impl_variant_runtimes
from offsite.tuning_scenario import TuningScenario
from offsite.util.file_extensions import __config_ext__, __tuning_scenario_ext__
from offsite.util.time import start_timer, stop_timer


def parse_program_args_app_bench_impl() -> Namespace:
    """Parse the passed program arguments.

    Parameters:
    -----------
    -

    Returns:
    --------
    argparse.Namespace
        Parsed program arguments.
    """
    # Create argument parser object.
    parser = ArgumentParser(description='Benchmark implementation variants.')
    # General options.
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(__version__),
                        help='Print program version and exit.')
    parser.add_argument('--verbose', action='store_true', default=False, help='Print further information on this run.')
    parser.add_argument('--db', action='store', required=True, help='Path to used database.')
    # Tuning scenario options.
    parser.add_argument('--scenario', action='store', required=True, type=Path,
                        help='Path to the used tuning scenario ({}) file which contains all required information on '
                             'the used machine configuration, compiler, implementation skeleton, kernel templates, ODE '
                             'methods and IVPs.'.format(__tuning_scenario_ext__))
    parser.add_argument('-n', '--ode-size', action='store', type=int, required=True,
                        help='Tune for a fixed ODE system size.')
    parser.add_argument('--cores', action='store', type=int, required=True, help='Tune for a fixed number of cores.')
    # Configuration options.
    parser.add_argument('--config', action='store', type=Path,
                        help='Customize autotuning process by passing a configuration ({}) file'.format(__config_ext__))
    # Parse program arguments.
    args = parser.parse_args()
    return args


def bench():
    """Run the offsite_bench_impl application.

    Parameters:
    -----------
    -

    Returns:
    --------
    -
    """
    config: Config = offsite.config.offsiteConfig

    # Start timer.
    ts = start_timer()

    # Open database.
    db_session: Session = open_db(config.args.db)

    # Parser phase.
    print('#' * 80 + '\n\nParser phase...\n')
    if config.args.verbose:
        print('#' * 80)
    # ... read tuning scenario.
    config.scenario = TuningScenario.from_file(config.args.scenario)
    config.scenario.solver = config.scenario.solver.to_database(db_session)
    assert config.scenario.solver.type in [SolverType.GENERIC, SolverType.ODE]
    # ... parse YAML descriptions and create corresponding Python objects ...
    machine, skeletons, templates, methods, ivps = parse_verify_yaml_desc(db_session, config.scenario)
    # ... print information on the parsed descriptions.
    if config.args.verbose:
        print_yaml_desc(machine, skeletons, templates, methods, ivps)
    # ... check if passed cores argument is valid for the given machine state.
    if config.args.cores > machine.coresPerSocket:
        raise RuntimeError('Error: Passed number of cores surpasses the used machine\'s maxime core count of {}'.format(
            machine.coresPerSocket))

    # Benchmarking phase.
    print('#' * 80 + '\n\nBenchmarking implementation variants...\n')
    # ... train database with benchmarked implementation variant runtimes.
    train_impl_variant_runtimes(db_session, machine, skeletons, methods, ivps)

    # Close database.
    close(db_session)

    # Stop timer.
    print('\n' + '#' * 80 + '\n')
    print('Tuning run took: {}.'.format(stop_timer(ts)))


def run():
    """Run the off_tune tuning application.

    Parameters:
    -----------
    -

    Returns:
    --------
    -
    """
    # Map database classes.
    mapping()
    # Parse the program arguments.
    args: Namespace = parse_program_args_app_bench_impl()
    # Create custom or default configuration.
    init_config(args)
    # Tune application.
    bench()


# For debugging purposes.
if __name__ == '__main__':
    run()
