"""@package train.node.util.node_communication
Functions to train the tuning database with communication costs.

@author: Johannes Seiferth
"""

from typing import List, Tuple

from sqlalchemy.orm import Session

import offsite.config
from offsite.config import Config
from offsite.descriptions.impl.impl_skeleton import ImplSkeleton
from offsite.descriptions.machine import MachineState
from offsite.descriptions.parser import load_yaml
from offsite.train.communication import AVAIL_BENCHMARKS
from offsite.train.communication.openmp.omp_barrier import OmpBarrierRecord
from offsite.tuning_scenario import TuningScenario
from offsite.util.sample_interval import SampleInterval

IntervalRecordList = List[Tuple[SampleInterval, str]]


def train_node_communication(
        db_session: Session, scenario: TuningScenario, machine: MachineState, skeletons: List[ImplSkeleton]):
    """Train database with benchmark data for the node-level communication operations.

    Runs the communication benchmarks required to estimate the communication costs of impl variants on a
    node-level.

    Parameters:
    -----------
    db_session: sqlalchemy.orm.session.Session
        Used database session.
    scenario: Tuning scenario.
        Used tuning scenario.
    machine: MachineState
        Used machine.
    skeletons: list of ImplSkeleton
        Used ImplSkeleton objects.

    Returns:
    --------
    -
    """
    config: Config = offsite.config.offsiteConfig
    # Select all benchmarks required by at least one ImplSkeleton object.
    required_benchmarks = set()
    for skeleton in skeletons:
        for operation in skeleton.communicationOperationsNodeLvl:
            try:
                required_benchmarks.add(AVAIL_BENCHMARKS[operation])
            except KeyError:
                raise RuntimeError('Missing benchmark for {}!'.format(operation))
    if not required_benchmarks:
        return
    # Parse benchmark data files.
    # TODO support multiple benchmark data files
    bench_data = {}
    if scenario.bench_omp_barrier is not None:
        yaml = load_yaml(scenario.bench_omp_barrier)
        if OmpBarrierRecord.check_bench_file(yaml, machine):
            name, data = OmpBarrierRecord.read_bench_file(yaml, machine)
            bench_data[name] = data
    # Run benchmarks and store the benchmark data obtained in the database.
    for benchmark in required_benchmarks:
        # If 'update_mode' is set check if benchmark was already run for the given configuration of machine and compiler
        # by checking whether the database contains fitting data.
        if config.args.update is True:
            # Database contains fitting data. Don't run benchmark again.
            if benchmark.contains(db_session, machine):
                continue
        # Use passed benchmark data if available.
        bench_name = benchmark.type.value
        if bench_name in bench_data:
            data = bench_data.get(bench_name)
        # else run benchmark.
        else:
            data = benchmark.run(machine, {'reps': config.repetitions_communication_operations})
        # Train database with benchmark results.
        benchmark.update(db_session, machine, data)
