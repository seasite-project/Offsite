"""@package offsite_rank
Main script of the offsite_rank application.
"""

from argparse import ArgumentParser, Namespace
from pathlib import Path
from time import time
from typing import List

from sqlalchemy.orm import Session

import offsite.config
from offsite import __version__
from offsite.config import __config_ext__, __ivp_ext__, __ode_method_ext__, init_config
from offsite.db.db import commit, close, open_db
from offsite.db.db_mapping import mapping
from offsite.descriptions.parser import print_yaml_desc, parse_ivps, parse_machine, parse_methods, parse_ranking_tasks
from offsite.descriptions.ranking_task import RankTask
from offsite.evaluation.ranking import create_rankings


def parse_program_args_app_rank() -> Namespace:
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
    parser.add_argument('--db', action='store', default='tune.db', help='Path to database. Default: tune.db.')
    parser.add_argument('--update', action='store_true', default=False,
                        help='Only execute autotuning steps (benchmark, performance modeling, runtime tests, ...) for '
                             'which the database is still missing fitting data for the considered autotuning '
                             'configuration. E.g., useful when adding an additional implementation skeleton, IVP, ...')
    parser.add_argument('-n', '--ode-size', action='store', type=int, default=None,
                        help='Tune for a fixed ODE system size only. If argument is not passed the working set model '
                             'is applied to determine the set of tested ODE system sizes.')
    # Mandatory options.
    parser.add_argument('--machine', action='store', required=True, type=Path,
                        help='Path to YAML machine description (.yaml) file.')
    parser.add_argument('--compiler', action='store', required=True, help='Name of the compiler.')
    parser.add_argument('--method', action='store', required=True, type=Path,
                        help='Path to either a folder containing YAML ODE method description files or to a particular '
                             'ODE method description ({}) file'.format(__ode_method_ext__))
    parser.add_argument('--ivp', action='store', required=True, type=Path,
                        help='Path to either a folder containing YAML IVP description files or to a particular IVP '
                             'description ({}) file.'.format(__ivp_ext__))
    parser.add_argument('--tasks', action='store', required=True, type=str,
                        help='TODO')
    # Configuration options.
    parser.add_argument('--config', action='store', type=Path,
                        help='Customize autotuning process by passing a configuration ({}) file'.format(__config_ext__))
    # Parse program arguments.
    return parser.parse_args()


def rank():
    """Run the offsite_rank application for ODEs.

    Parameters:
    -----------
    -

    Returns:
    --------
    -
    """
    args: Namespace = offsite.config.offsiteConfig.args
    verbose: bool = args.verbose
    # Start timer.
    if verbose:
        start_time = time()
    # Open database.
    db_session: Session = open_db(args.db)
    # Parser phase.
    print('#' * 80 + '\n\nParser phase...\n')
    if verbose:
        print('#' * 80)
    # Parse YAML descriptions and create corresponding Python objects.
    # .. machine
    machine = parse_machine(args.machine, args.compiler)
    machine = machine.to_database(db_session)
    # .. IVP
    ivps = parse_ivps(args.ivp, None)
    ivps = [ivp.to_database(db_session) for ivp in ivps]
    # .. ODE method
    methods = parse_methods(args.method)
    methods = [method.to_database(db_session) for method in methods]
    # Print information on the parsed descriptions.
    if verbose:
        print_yaml_desc(machine, None, None, methods, ivps)
    # Ranking phase.
    print('#' * 80 + '\n\nRanking phase...\n')
    if verbose:
        start_time_ranking = time()
    # .. create ranking tasks.
    rank_tasks: List[RankTask] = parse_ranking_tasks(args.tasks)
    if len(rank_tasks) == 0:
        raise RuntimeError('Failed to create any ranking tasks!')
    # .. create rankings.
    create_rankings(db_session, machine, methods, ivps, rank_tasks, args.ode_size)
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
    args: Namespace = parse_program_args_app_rank()
    # Create custom or default configuration.
    init_config(args)
    # Rank application.
    rank()


# For debugging purposes.
if __name__ == '__main__':
    run()
