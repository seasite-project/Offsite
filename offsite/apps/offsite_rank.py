"""@package apps.offsite_rank
Main script of the offsite_rank application.
"""

from argparse import ArgumentParser, Namespace
from pathlib import Path
from time import time
from typing import Dict, List

from sqlalchemy.orm import Session

import offsite.config
from offsite import __version__
from offsite.config import __config_ext__, __ivp_ext__, __ode_method_ext__, __tuning_scenario_ext__, init_config
from offsite.database import close, commit, open_db
from offsite.database.db_mapping import mapping
from offsite.descriptions.machine.machine import parse_machine_state
from offsite.descriptions.ode import IVP, ODEMethod, parse_ivp, parse_ivps, parse_methods
from offsite.descriptions.parser_util import print_yaml_desc
from offsite.ranking.ranking import create_rankings, create_rankings_ode
from offsite.ranking.ranking_task import RankTask, parse_ranking_tasks
from offsite.solver import Solver, SolverType, SolverSpecificTableType
from offsite.tuning_scenario import TuningScenario


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
    parser.add_argument('--db', action='store', required=True, help='Path to used database.')
    parser.add_argument('--update', action='store_true', default=False,
                        help='Only execute autotuning steps (benchmark, performance modeling, runtime tests, ...) for '
                             'which the database is still missing fitting data for the considered autotuning '
                             'configuration. E.g., useful when adding an additional impl skeleton, IVP, ...')
    # Tuning scenario options.
    parser.add_argument('--scenario', action='store', required=True, type=Path,
                        help='Path to the used tuning scenario ({}) file which contains all required information on '
                             'the used machine state, compiler, impl skeleton, kernel templates, ODE methods '
                             'and IVPs.'.format(__tuning_scenario_ext__))
    parser.add_argument('-n', '--ode-size', action='store', type=int, default=None,
                        help='Tune for a fixed ODE system size only. If argument is not passed the working set model '
                             'is applied to determine the set of tested ODE system sizes.')
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
    # Read tuning scenario.
    scenario: TuningScenario = TuningScenario.from_file(args.scenario)
    solver: Solver = scenario.solver
    # Parse YAML descriptions and create corresponding Python objects.
    # .. machine state
    machine = parse_machine_state(scenario.machine, scenario.compiler)
    machine = machine.to_database(db_session)
    # ... ODE method descriptions if required by solver.
    methods = None
    if SolverSpecificTableType.ODE_METHOD in scenario.solver.specific_tables:
        methods_dict: Dict[str, ODEMethod] = dict()
        if scenario.method_path is not None:
            for p in (p for p in scenario.method_path if p.suffix == __ode_method_ext__):
                m = ODEMethod.from_yaml(p)
                methods_dict[m.name] = m
        if scenario.method_dir is not None:
            for d in scenario.method_dir:
                for m in parse_methods(d):
                    methods_dict[m.name] = m
        methods: List[ODEMethod] = [method.to_database(db_session) for method in methods_dict.values()]
        if not methods:
            raise RuntimeError('No valid ODE methods found: \'{}\''.format(args.method))
    # ... IVP descriptions if required by solver.
    ivps = None
    if SolverSpecificTableType.IVP in scenario.solver.specific_tables:
        ivps_dict: Dict[str, IVP] = dict()
        if scenario.ivp_path is not None:
            for p in (p for p in scenario.ivp_path if p.suffix == __ivp_ext__):
                i = parse_ivp(p, None)
                ivps_dict[i.name] = i
        if scenario.ivp_dir is not None:
            for d in scenario.ivp_dir:
                for i in parse_ivps(d):
                    if i.name not in methods:
                        ivps_dict[i.name] = i
        ivps: List[IVP] = [ivp.to_database(db_session) for ivp in ivps_dict.values()]
        if not ivps:
            raise RuntimeError('No valid IVPs found: \'{}\''.format(args.ivp))
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
    if solver.type == SolverType.ODE:
        create_rankings_ode(db_session, machine, methods, ivps, rank_tasks, args.ode_size)
    else:
        create_rankings(db_session, machine, rank_tasks, args.ode_size)
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
