"""@package offsite_tune
Main script of the offsite_tune autotuning application.
"""

from argparse import ArgumentParser
from pathlib import Path
from time import time

import offsite.config
from offsite import __version__
from offsite.config import ModelToolType, ProgramModeType, IncoreToolType, __bench_ext__, __config_ext__, \
    __impl_skeleton_ext__, __ivp_ext__, __kernel_template_ext__, __ode_method_ext__, init_config
from offsite.db.db import commit, close, open_db
from offsite.db.db_mapping import mapping
from offsite.descriptions.parser import parse_verify_yaml_desc, print_yaml_desc
from offsite.evaluation.ranking import rank
from offsite.train.train_communication import train_communication_costs
from offsite.train.train_impl import train_impl_variant_predictions
from offsite.train.train_kernel import train_kernel_predictions, train_kernel_runtimes
from offsite.train.train_kernel_blocksize import train_kernel_blocksizes


def parse_program_args_app_tune() -> 'argparse.Namespace':
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
    parser.add_argument('--tool', action='store', type=ModelToolType, default=None,
                        help='Can be used to limit the performance model phase to a particular modelling tool. If set '
                             'only those implementations variants and kernels predictable by the modeling tool '
                             'selected are considered in this autotuning run. By default this option is not set and '
                             'all tools are used. Available model tools: KERNCRAFT, YASKSITE')
    parser.add_argument('--incore', action='store', type=IncoreToolType, default=IncoreToolType.OSACA,
                        help='Select which in-core model tool Kerncraft uses internally. Available in-core models: '
                             'OSACA(default), IACA, LLVM-MCA')
    parser.add_argument('--db', action='store', default='tune.db', help='Path to database. Default: tune.db.')
    parser.add_argument('--update', action='store_true', default=False,
                        help='Only execute autotuning steps (benchmark, performance modeling, runtime tests, ...) for '
                             'which the database is still missing fitting data for the considered autotuning '
                             'configuration. E.g., useful when adding an additional implementation skeleton, IVP, ...')
    parser.add_argument('-n', '--ode-size', action='store', type=int, default=False,
                        help='Tune for a fixed ODE system size only. If argument is not  passed the working set model '
                             'is applied to determine the set of tested ODE system sizes.')
    # Mandatory options.
    parser.add_argument('--machine', action='store', required=True, type=Path,
                        help='Path to YAML machine description (.yaml) file.')
    parser.add_argument('--compiler', action='store', required=True, help='Name of the compiler.')
    parser.add_argument('--kernel', action='store', required=True, type=Path,
                        help='Path to folder containing YAML kernel template description ({}) files.'.format(
                            __kernel_template_ext__))
    parser.add_argument('--impl', action='store', required=True, type=Path,
                        help='Path to folder containing YAML implementation skeleton description ({}) files.'.format(
                            __impl_skeleton_ext__))
    parser.add_argument('--method', action='store', required=True, type=Path,
                        help='Path to either a folder containing YAML ODE method description files or to a particular '
                             'ODE method description ({}) file'.format(__ode_method_ext__))
    parser.add_argument('--ivp', action='store', required=True, type=Path,
                        help='Path to either a folder containing YAML IVP description files or to a particular IVP '
                             'description ({}) file.'.format(__ivp_ext__))
    # Configuration options.
    parser.add_argument('--config', action='store', type=Path,
                        help='Customize autotuning process by passing a configuration ({}) file'.format(__config_ext__))
    parser.add_argument('--tolerance', action='store', default=5.0, type=float,
                        help='Set tolerance that controls how far implementation variants may deviate from the best '
                             'found implementation variant in order to be incorporated during the ranking phase. '
                             'Default: 5.0')
    parser.add_argument('--ws-interpolate', action='store_true', default=False,
                        help='Interpolate prediction/runtime of the sample interval created in particular border '
                             'regions.')
    parser.add_argument('--filter-yasksite-opt', action='store_true', default=False,
                        help='Filter out all implementation variants whose kernels do not all use the same'
                             'YaskSite optimization parameters.')
    # Benchmark options.
    parser.add_argument('--bench', action='store', type=Path,
                        help='Instead of running a benchmark pass results as ({}) file'.format(__bench_ext__))
    # Parse program arguments.
    return parser.parse_args()


def tune():
    """Run the offsite_tune application for ODEs.

    Parameters:
    -----------
    -

    Returns:
    --------
    -
    """
    args = offsite.config.offsiteConfig.args
    verbose = args.verbose
    # Start timer.
    if verbose:
        start_time = time()
    # Open database.
    db_session = open_db(args.db)
    # Parser phase.
    print('#' * 80 + '\n\nParser phase...\n')
    if verbose:
        print('#' * 80)
    # Parse YAML descriptions and create corresponding Python objects.
    machine, skeletons, templates, methods, ivps = parse_verify_yaml_desc(db_session)
    # Print information on the parsed descriptions.
    if verbose:
        print_yaml_desc(machine, skeletons, templates, methods, ivps)
    # Training phase.
    print('#' * 80 + '\n\nTraining phase...\n')
    # Train database with communication cost benchmarks.
    if verbose:
        print('  * Communication costs...', end='', flush=True)
    train_communication_costs(db_session, machine, skeletons)
    if verbose:
        print(' done.')
    # Train database with kernel block size predictions.
    if args.tool == ModelToolType.KERNCRAFT or args.tool is None:
        if verbose:
            print('  * Kernel block sizes...', end='', flush=True)
            start_time_block = time()
        templates = train_kernel_blocksizes(db_session, machine, templates, methods, ivps)
        templates = [template.to_database(db_session) for template in templates]
        if verbose:
            print(' done.')
            stop_time_block = time()
            print('#' * 80 + '\n')
            print('Kernel block size phase took {} seconds.'.format(round(stop_time_block - start_time_block, 3)))
    # Train database with kernel runtime predictions.
    if verbose:
        print('  * Kernel predictions...', end='', flush=True)
        start_time_kernel = time()
    if args.mode == ProgramModeType.MODEL:
        train_kernel_predictions(db_session, machine, templates, methods, ivps)
    elif args.mode == ProgramModeType.RUN:
        train_kernel_runtimes(db_session, machine, templates, methods, ivps)
    else:
        raise RuntimeError('Unsupported program mode!')
    # Commit to database.
    commit(db_session)
    if verbose:
        print(' done.')
        stop_time_kernel = time()
        print('#' * 80 + '\n')
        print('Kernel prediction phase took {} seconds.'.format(round(stop_time_kernel - start_time_kernel, 3)))
    # Train database with implementation variant runtime predictions.
    if verbose:
        print('  * Implementation variant predictions...', end='', flush=True)
        start_time_impl = time()
    train_impl_variant_predictions(db_session, machine, skeletons, methods, ivps)
    # Commit to database.
    commit(db_session)
    if verbose:
        print(' done.\n')
        stop_time_impl = time()
        print('#' * 80 + '\n')
        print('Impl prediction phase took {} seconds.'.format(round(stop_time_impl - start_time_impl, 3)))
    # Ranking phase.
    print('#' * 80 + '\n\nRanking phase...\n')
    if verbose:
        start_time_ranking = time()
    rank(db_session, machine, methods, ivps)
    if verbose:
        print(' done.\n')
        stop_time_ranking = time()
        print('#' * 80 + '\n')
        print('Ranking phase took {} seconds.'.format(round(stop_time_ranking - start_time_ranking, 3)))
    # Commit to database.
    commit(db_session)
    # Close database.
    close(db_session)
    # Stop timer.
    if verbose:
        stop_time = time()
        print('#' * 80 + '\n')
        print('Tuning run took {} seconds.'.format(round(stop_time - start_time, 3)))


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
    args = parse_program_args_app_tune()
    # Create custom or default configuration.
    init_config(args)
    # Tune application.
    tune()


# For debugging purposes.
if __name__ == '__main__':
    run()
