"""@package apps.offsite_tune_node
Main script of the offsite_tune_node autotuning application.
"""

from argparse import ArgumentParser, Namespace
from datetime import timedelta
from pathlib import Path
from time import time
from typing import List

from sqlalchemy.orm import Session

import offsite.config
from offsite import __version__
from offsite.config import Config, ModelToolType, ProgramModeType, SolverType, __config_ext__, \
    __tuning_scenario_ext__, init_config
from offsite.database import close, open_db
from offsite.database.db_mapping import mapping
from offsite.descriptions.impl.kernel_template import Kernel
from offsite.descriptions.parser_util import parse_verify_yaml_desc, print_yaml_desc
from offsite.train.communication.train_communication import train_node_communication
from offsite.train.node.train_impl import train_impl_variant
from offsite.train.node.train_kernel import train_kernel
from offsite.train.node.train_kernel_blocksize import train_kernel_blocksizes
from offsite.tuning_scenario import TuningScenario
from offsite.util.sample_interval import SampleInterval, SampleType


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
    return parser.parse_args()


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
    start_time = time()
    # Open database.
    db_session: Session = open_db(config.args.db)
    # Parser phase.
    print('#' * 80 + '\n\nParser phase...\n')
    if config.args.verbose:
        print('#' * 80)
    # Read tuning scenario.
    config.scenario = TuningScenario.from_file(config.args.scenario)
    assert config.scenario.solver.type in [SolverType.GENERIC, SolverType.ODE]
    # Parse YAML descriptions and create corresponding Python objects ...
    machine, skeletons, templates, methods, ivps = parse_verify_yaml_desc(db_session, config.scenario)
    # Print information on the parsed descriptions.
    if config.args.verbose:
        print_yaml_desc(machine, skeletons, templates, methods, ivps)
    # Training phase.
    print('#' * 80 + '\n\nTraining phase...\n')
    # Train database with communication cost benchmarks.
    if config.args.verbose:
        print('  * Communication costs...', end='', flush=True)
    train_node_communication(db_session, config.scenario, machine, skeletons)
    if config.args.verbose:
        print(' done.')
    # Train database with kernel block size predictions.
    if not config.args.no_blocksizes and (
            config.pred_model_tool == ModelToolType.KERNCRAFT or config.pred_model_tool is None):
        templates = train_kernel_blocksizes(db_session, machine, templates, methods, ivps)
    # Train database with kernel runtime predictions.
    predicted_kernels: List[Kernel] = train_kernel(db_session, machine, templates, ivps, methods)
    # Train database with implementation variant runtime predictions.
    train_impl_variant(db_session, machine, skeletons, predicted_kernels, ivps, methods)
    # Close database.
    close(db_session)
    # Stop timer.
    stop_time = time()
    print('\n' + '#' * 80 + '\n')
    print('Tuning run took: {}.'.format(timedelta(seconds=round(stop_time - start_time, 0))))


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
    args.nrange = derive_samples(args.nrange, args.mode) if args.nrange else None
    # Create custom or default configuration.
    init_config(args)
    # Tune application.
    tune()


def derive_samples(intv_range: str, mode: ProgramModeType):
    if mode == ProgramModeType.MODEL:
        sample_type = SampleType.MODEL_INNER
    elif mode == ProgramModeType.RUN:
        sample_type = SampleType.BENCH_INNER
    else:
        assert False

    intervals = list()
    # Split interval range string.
    intv_range = intv_range.split(':')
    assert (len(intv_range) == 3)
    cur_sample = int(intv_range[0])
    last_sample = int(intv_range[1])
    incr = int(intv_range[2])
    # First interval starts with first sample.
    e = cur_sample + int(incr / 2)
    intervals.append(SampleInterval(cur_sample, e, sample_type, cur_sample))
    # Add inner intervals.
    cur_sample += incr
    s = e + 1
    e += incr
    while cur_sample < last_sample:
        intervals.append(SampleInterval(s, e, sample_type, cur_sample))
        cur_sample += incr
        s = e + 1
        e += incr
    # Add last interval.
    intervals.append(SampleInterval(s, last_sample, sample_type, last_sample))
    return intervals


# For debugging purposes.
if __name__ == '__main__':
    run()
