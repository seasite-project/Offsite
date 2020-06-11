"""@package train_communication
Functions to train the tuning database with communication costs.
"""

from typing import List, Tuple

from sqlalchemy.orm import Session

import offsite.config
from offsite.descriptions.impl_skeleton import ImplSkeleton
from offsite.descriptions.machine import Machine
from offsite.descriptions.parser_utils import load_yaml
from offsite.evaluation.benchmark import AVAIL_BENCHMARKS, BenchmarkRecord
from offsite.evaluation.performance_model import SampleInterval

IntervalRecordList = List[Tuple[SampleInterval, str]]


def train_communication_costs(db_session: Session, machine: Machine, skeletons: List[ImplSkeleton]):
    """Train database with communication benchmark data.

    Runs communication benchmarks required to estimate the communication costs of the variants derived from the given
    implementation skeletons.

    Parameters:
    -----------
    db_session : sqlalchemy.orm.session.Session
        Used database session.
    machine : Machine
        Used machine.
    skeletons : list of ImplSkeleton
        Used ImplSkeleton objects.

    Returns:
    --------
    -
    """
    config = offsite.config.offsiteConfig
    # Select all benchmarks required by at least one ImplSkeleton object.
    required_benchmarks = set()
    for skeleton in skeletons:
        for operation in skeleton.communicationOperations:
            if operation not in AVAIL_BENCHMARKS:
                raise RuntimeError('Missing benchmark for {}!'.format(operation))
            required_benchmarks.add(AVAIL_BENCHMARKS[operation])
    if not required_benchmarks:
        return
    # Parse benchmark data files.
    bench_data = {}
    if config.args.bench:
        data = load_yaml(config.args.bench)
        # Check machine name.
        if data['machine'] != machine.name:
            raise RuntimeError('Passed benchmark data \'{}\' were raised for a different machine \'{}\'! Data for '
                               'machine \'{}\' required!'.format(config.args.bench, data['machine'], machine.name))
        # Check compiler version.
        name, version = data['compiler'].split(' ')
        if name != machine.compiler.name:
            raise RuntimeError(
                'Passed benchmark data \'{}\' were raised for a different compiler \'{}\'! Data for compiler \'{}\' '
                'required!'.format(config.args.bench, name, machine.compiler.name))
        if version != machine.compiler.version:
            raise RuntimeError(
                'Passed benchmark data \'{}\' were raised for a different compiler version \'{}\'! Data for compiler '
                'version \'{}\' required!'.format(config.args.bench, version, machine.compiler.version))
        # Check compiler flags.
        flags = data['flags']
        if flags != machine.compiler.flags:
            raise RuntimeError(
                'Passed benchmark data \'{}\' were raised for a different set of compiler flags \'{}\'! Data for '
                'compiler flags \'{}\' required!'.format(config.args.bench, flags, machine.compiler.flags))
        # Check CPU frequency.
        if data['frequency'] != frequency:
            raise RuntimeError(
                'Passed benchmark data \'{}\' were raised for different CPU frequency \'{} Hz\'! Data for frequency '
                '\'{} Hz\'! required!'.format(config.args.bench, data['frequency'], machine.clock))
        # Convert data to Benchmark records.
        bench_name = data['benchmark']
        bench_data[bench_name] = [BenchmarkRecord(
            bench_name, machine.db_id, machine.compiler.db_id, row[1], machine.clock, row[0]) for row in data['data']]
    # Run benchmarks and store the benchmark data obtained in the database.
    for benchmark in required_benchmarks:
        # If 'update_mode' is set check if benchmark was already run for the given configuration of machine and compiler
        # by checking whether the database contains fitting data.
        if config.args.update is True:
            # Database contains fitting data. Don't run benchmark again.
            if BenchmarkRecord.contains(db_session, machine, benchmark):
                continue
        # Use passed benchmark data if available
        bench_name = benchmark.name
        if bench_name in bench_data:
            data = bench_data.get(bench_name)
        # else run benchmark.
        else:
            data = benchmark.run(machine, config.repetitions_communication_operations)
        # Train database with benchmark results.
        BenchmarkRecord.update(db_session, machine, benchmark, data)
