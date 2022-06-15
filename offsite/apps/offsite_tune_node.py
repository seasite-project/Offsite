"""@package apps.offsite_tune_node
Main script of the offsite_tune_node autotuning application.

@author: Johannes Seiferth
"""

from typing import List

from argparse import ArgumentParser, Namespace
from pathlib2 import Path
from sqlalchemy.orm import Session

import offsite.config
from offsite import __version__
from offsite.config import Config, ModelToolType, ProgramModeType, init_config
from offsite.database import close, open_db
from offsite.database.db_mapping import mapping
from offsite.descriptions.impl.kernel_template import Kernel
from offsite.descriptions.parser_util import parse_verify_yaml_desc, print_yaml_desc
from offsite.solver import SolverType
from offsite.train.communication.train_communication import train_node_communication
from offsite.train.node.train_impl import train_impl_variant_predictions
from offsite.train.node.train_kernel import train_kernel
from offsite.train.node.train_kernel_blocksize import train_kernel_blocksizes
from offsite.tuning_scenario import TuningScenario
from offsite.util.file_extensions import __config_ext__, __tuning_scenario_ext__
from offsite.util.sample_interval import derive_samples_from_range_expr
from offsite.util.time import start_timer, stop_timer


def parse_program_args_app_tune_node() -> Namespace:
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
    parser = ArgumentParser(description='Offline autotuner for time-step based simulations of ODE solvers.')
    # General options.
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(__version__),
                        help='Print program version and exit.')
    parser.add_argument('--verbose', action='store_true', default=False, help='Print further information on this run.')
    parser.add_argument('--mode', action='store', type=ProgramModeType, default=ProgramModeType.MODEL,
                        help='Available program modes: MODEL (default), RUN')
    parser.add_argument('--db', action='store', required=True, help='Path to used database.')
    # Tuning scenario options.
    parser.add_argument('--scenario', action='store', required=True, type=Path,
                        help='Path to the used tuning scenario ({}) file which contains all required information on '
                             'the used machine configuration, compiler, implementation skeleton, kernel templates, ODE '
                             'methods and IVPs.'.format(__tuning_scenario_ext__))
    parser.add_argument('-n', '--ode-size', action='store', type=int, default=None,
                        help='Tune for a fixed ODE system size only. This argument is prioritized over argument '
                             '--nrange if both are set! If neither argument is passed the working set model is applied '
                             'to determine the set of tested ODE system sizes.')
    parser.add_argument('--nrange', action='store', type=str, default=None,
                        help='Tune for a fixed range of ODE system sizes only. Range ist passed in format: '
                             '\"first:last:increment\". Argument -n is prioritized over this argument if both are set! '
                             'If neither argument is passed the working set model is applied to determine the set of '
                             'tested ODE system sizes.')
    # Configuration options.
    parser.add_argument('--config', action='store', type=Path,
                        help='Customize autotuning process by passing a configuration ({}) file'.format(__config_ext__))
    parser.add_argument('--update', action='store_true', default=False,
                        help='Only execute autotuning steps (benchmark, performance modeling, runtime tests, ...) for '
                             'which the database is still missing fitting data for the considered autotuning '
                             'configuration. E.g., useful when adding an additional implementation skeleton, IVP, ...')
    parser.add_argument('--filter-yasksite-opt', action='store_true', default=False,
                        help='Filter out all implementation variants whose kernels do not all use the same '
                             'YaskSite optimization parameters.')
    parser.add_argument('--include-all-kernels', action='store_true', default=False,
                        help='If activated predictions will be obtained for all available kernel templates. By default '
                             'predictions are only obtained for those kernel templates actually used by any of the '
                             'given implementation skeletons.')
    parser.add_argument('--bench-rhs', action='store_true', default=False,
                        help='Benchmark the runtime of all IVP-dependent kernels instead of analytically predicting '
                             'them. This mode considers all boundary cases and computations.')
    parser.add_argument('--no-blocksizes', action='store_true', default=False,
                        help='Skip kernel block size prediction.')
    # Parse program arguments.
    args = parser.parse_args()
    # ...
    args.nrange = derive_samples_from_range_expr(args.nrange, args.mode) if args.nrange else None
    return args


def tune():
    """Run the offsite_tune_node application.

    Parameters:
    -----------
    -

    Returns:
    --------
    -
    """
    config: Config = offsite.config.offsiteConfig
    assert config.args.mode in [ProgramModeType.MODEL, ProgramModeType.RUN]

    # Start timer.
    ts = start_timer()

    # Open database.
    db_session: Session = open_db(config.args.db)

    # Parser phase.
    print('#' * 80 + '\n\nParser phase...\n')
    if config.args.verbose:
        print('#' * 80)
    # ... store solver type in database.
    config.scenario = TuningScenario.from_file(config.args.scenario)
    config.scenario.solver = config.scenario.solver.to_database(db_session)
    assert config.scenario.solver.type in [SolverType.GENERIC, SolverType.ODE]
    # ... parse YAML descriptions and create corresponding Python objects ...
    machine, skeletons, templates, methods, ivps = parse_verify_yaml_desc(db_session, config.scenario)
    # ... print information on the parsed descriptions.
    if config.args.verbose:
        print_yaml_desc(machine, skeletons, templates, methods, ivps)

    # Training phase.
    print('#' * 80 + '\n\nTraining phase...\n')
    # ... train database with communication cost benchmarks.
    if config.args.verbose:
        print('  * Communication costs...', end='', flush=True)
    train_node_communication(db_session, config.scenario, machine, skeletons)
    if config.args.verbose:
        print(' done.')
    # ... train database with kernel block size predictions.
    if not config.args.no_blocksizes and (
            config.pred_model_tool == ModelToolType.KERNCRAFT or config.pred_model_tool is None):
        templates = train_kernel_blocksizes(db_session, machine, templates, methods, ivps)
    # ... train database with kernel runtime predictions.
    predicted_kernels: List[Kernel] = train_kernel(db_session, machine, templates, ivps, methods)
    # ... train database with implementation variant runtime predictions.
    train_impl_variant_predictions(db_session, machine, skeletons, predicted_kernels, methods, ivps)

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
    args: Namespace = parse_program_args_app_tune_node()
    # Create custom or default configuration.
    init_config(args)
    # Tune application.
    tune()


# For debugging purposes.
if __name__ == '__main__':
    run()
