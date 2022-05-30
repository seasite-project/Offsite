"""@package ranking.ranking
Definition of class RankingRecord and of different ranking functions.

@author: Johannes Seiferth
"""

from copy import deepcopy
from datetime import datetime
from getpass import getuser
from multiprocessing import cpu_count, Pool
from sys import maxsize as sys_maxsize
from traceback import print_exc
from typing import List, Optional, Set, Tuple, Union

import attr
from pandas import read_sql_query, DataFrame
from sqlalchemy import Column, DateTime, Enum, Float, ForeignKey, Integer, String, Table, UniqueConstraint
from sqlalchemy.orm import Query, Session

from offsite import __version__
from offsite.database import METADATA
from offsite.database.db import insert, bulk_insert
from offsite.descriptions.machine import MachineState
from offsite.descriptions.ode import IVP, ODEMethod
from offsite.ranking.ranking_task import RankDeviationTask, RankOrderTask, RankTask, RankingCriteriaType, \
    RankingCutoffType
from offsite.train.impl_variant import ImplVariant, ImplVariantRecord
from offsite.util.sample_interval import SampleInterval

RankingData = List[Tuple[SampleInterval, List[int]]]


@attr.s
class RankingRecord:
    """Representation of an implementation variant ranking database table record.

    Attributes:
    -----------
    machine: int
        Used machine.
    compiler: int
        Used compiler.
    method: int
        Used ODE method.
    ivp: int
        Used IVP.
    cores: int
        Used number of cores.
    frequency: float
        Used CPU frequency.
    sample: SampleInterval
        Used sample interval.
    variants: Set of int
        ID's of implementation variants.
    db_id: int
        ID of associated ranking database table record.
    """
    machine = attr.ib(type=int)
    compiler = attr.ib(type=int)
    method = attr.ib(type=int)
    ivp = attr.ib(type=int)
    cores = attr.ib(type=int)
    frequency = attr.ib(type=float)
    sample = attr.ib(type=SampleInterval)
    variants = attr.ib(type=Set[int])
    rankCriteria = attr.ib(type=RankingCriteriaType)
    cutoffCriteria = attr.ib(type=RankingCutoffType)
    cutoffValue = attr.ib(type=float)
    variants_serial = attr.ib(type=str, init=False)
    first = attr.ib(type=int, init=False)
    last = attr.ib(type=int, init=False)
    db_id = attr.ib(type=int, init=False)

    # Database information.
    db_table = Table('ranking', METADATA,
                     Column('db_id', Integer, primary_key=True),
                     Column('machine', Integer, ForeignKey('machine.db_id')),
                     Column('compiler', Integer, ForeignKey('compiler.db_id')),
                     Column('method', Integer, ForeignKey('ode_method.db_id')),
                     Column('ivp', Integer, ForeignKey('ivp.db_id')),
                     Column('frequency', Float),
                     Column('cores', Integer),
                     Column('first', Integer),
                     Column('last', Integer),
                     Column('variants_serial', String),
                     Column('rankCriteria', Enum(RankingCriteriaType)),
                     Column('cutoffCriteria', Enum(RankingCutoffType)),
                     Column('cutoffValue', Float),
                     Column('updatedIn', String, default=__version__),
                     Column('updatedOn', DateTime, default=datetime.now, onupdate=datetime.now),
                     Column('updatedBy', String, default=getuser(), onupdate=getuser()),
                     UniqueConstraint('machine', 'compiler', 'method', 'ivp', 'frequency', 'cores', 'first', 'last',
                                      'rankCriteria', 'cutoffCriteria', 'cutoffValue'),
                     sqlite_autoincrement=True)

    def __attrs_post_init__(self):
        """Create this ranking record object.

        Parameters:
        -----------
        -

        Returns:
        --------
        -
        """
        self.first = self.sample.first
        self.last = self.sample.last
        self.variants_serial = ','.join(map(str, self.variants))

    def to_database(self, db_session: Session):
        """Push this ranking record object to the database.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        -
        """
        # Attribute variants_serial.
        self.variants_serial = ','.join(map(str, self.variants))
        # Insert RankingRecord object.
        insert(db_session, self)

    @staticmethod
    def update(db_session: Session, records: List['RankingRecord']):
        """ Update data records in Ranking table with new ranking data.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.
        records : list
            Ranking data.

        Returns:
        --------
        -
        """
        bulk_insert(db_session, records)

    @staticmethod
    def remove_record(db_session: Session, machine: int, compiler: int, method: int, ivp: int, frequency: float,
                      core_counts: List[int], first: int, last: int, rank_criteria: RankingCriteriaType,
                      cutoff_criteria: RankingCutoffType, cutoff_value: float):
        """Remove a single ranking data record from the database.

        Remove the record that matches the provided configuration of machine configuration, compiler, ODE method, IVP,
        CPU frequency, number of CPU cores and ODE system size.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        machine: int
            Used machine.
        compiler: int
            Used compiler.
        method: int
            Used ODE method.
        ivp: int
            Used IVP.
        frequency: float
            Used CPU frequency.
        core_counts: List of int
            Used core counts.
        rank_criteria: RankingCriteriaType
            Used ranking criteria.
        cutoff_criteria: CutoffCriteriaType
            Used cutoff criteria.
        cutoff_value: float
            Used cutoff value.

        Returns:
        --------
        -
        """
        db_session.query(RankingRecord).filter(
            RankingRecord.machine.is_(machine), RankingRecord.compiler.is_(compiler), RankingRecord.method.is_(method),
            RankingRecord.ivp.is_(ivp), RankingRecord.frequency.is_(frequency), RankingRecord.cores.in_(core_counts),
            RankingRecord.first.is_(first), RankingRecord.last.is_(last), RankingRecord.rankCriteria.is_(rank_criteria),
            RankingRecord.cutoffCriteria.is_(cutoff_criteria), RankingRecord.cutoffValue.is_(cutoff_value)).delete(
            synchronize_session=False)

    @staticmethod
    def remove_records(db_session: Session, machine: int, compiler: int, method: int, ivp: int, frequency: float,
                       core_counts: List[int], rank_criteria: RankingCriteriaType, cutoff_criteria: RankingCutoffType,
                       cutoff_value: float):
        """Remove ranking data record(s) from the database.

        Remove all records that match the provided configuration of machine state, compiler, ODE method, IVP, CPU
        frequency and number of CPU cores.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        machine: int
            Used machine.
        compiler: int
            Used compiler.
        method: int
            Used ODE method.
        ivp: int
            Used IVP.
        frequency: float
            Used CPU frequency.
        core_counts : List of int
            Used core counts.
        rank_criteria: RankingCriteriaType
            Used ranking criteria.
        cutoff_criteria: CutoffCriteriaType
            Used cutoff criteria.
        cutoff_value: float
            Used cutoff value.

        Returns:
        --------
        -
        """
        db_session.query(RankingRecord).filter(
            RankingRecord.machine.is_(machine), RankingRecord.compiler.is_(compiler), RankingRecord.method.is_(method),
            RankingRecord.ivp.is_(ivp), RankingRecord.frequency.is_(frequency), RankingRecord.cores.in_(core_counts),
            RankingRecord.rankCriteria.is_(rank_criteria), RankingRecord.cutoffCriteria.is_(cutoff_criteria),
            RankingRecord.cutoffValue.is_(cutoff_value)).delete(synchronize_session=False)

    @staticmethod
    def select(db_session: Session, machine: MachineState, method: int, ivp: int, task: RankTask,
               ode_size: int) -> DataFrame:
        """Retrieve RankingRecord table data record(s) from the database.

        Return all records that match the provided configuration of machine, compiler and benchmark name(s).

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        machine: MachineState
            Used machine.
        method: int
            Database ID used ODE method.
        ivp: int
            Database ID used IVP.
        task: RankTask
            Used ranking task.
        ode_size: int
            Used ODE system size.

        Returns:
        --------
        pandas.DataFrame
            Retrieved list of data records.
        """
        query: Query
        if ode_size:
            query = db_session.query(RankingRecord.variants_serial).filter(
                RankingRecord.machine.is_(machine.db_id), RankingRecord.compiler.is_(machine.compiler.db_id),
                RankingRecord.method.is_(method), RankingRecord.ivp.is_(ivp),
                RankingRecord.frequency.is_(machine.clock), RankingRecord.first.is_(ode_size),
                RankingRecord.last.is_(ode_size), RankingRecord.rankCriteria.is_(task.type),
                RankingRecord.cutoffCriteria.is_(task.cutoff_criteria),
                RankingRecord.cutoffValue.is_(task.cutoff_value)).distinct()
        else:
            query = db_session.query(RankingRecord.variants_serial).filter(
                RankingRecord.machine.is_(machine.db_id), RankingRecord.compiler.is_(machine.compiler.db_id),
                RankingRecord.method.is_(method), RankingRecord.ivp.is_(ivp),
                RankingRecord.frequency.is_(machine.clock), RankingRecord.rankCriteria.is_(task.type),
                RankingRecord.cutoffCriteria.is_(task.cutoff_criteria),
                RankingRecord.cutoffValue.is_(task.cutoff_value)).distinct()
        return read_sql_query(query.statement, db_session.bind)


def fuse_equal_rankings(rankings: RankingData) -> RankingData:
    """
    Reduce total number of rankings by combining neighbouring rankings. Adjacent rankings can be combined if they
    contain the same impl variants.

    Parameters:
    -----------
    rankings: list of tuple(SampleInterval, set of impl variant ids)
        Rankings of impl variants for ascending intervals.

    Returns:
    --------
    list of tuple(SampleInterval, set of impl variant ids)
        Rankings of impl variants for ascending intervals after fusing.
    """
    if len(rankings) <= 1:
        return rankings
    rankings_cpy = deepcopy(rankings)
    # Iterate over the rankings and combine those that contain the same impl variants.
    combined = list()
    r_prev = set()
    intv = None
    for ranking in rankings_cpy:
        i_cur = ranking[0]
        r_cur = ranking[1]
        if r_prev == r_cur:
            # Increase upper bound of current interval.
            intv.last = i_cur.last
        else:
            # Save interval.
            if intv:
                combined.append((intv, r_prev))
            # Start new interval.
            r_prev = r_cur
            intv = SampleInterval(i_cur.first, i_cur.last, None)
    # Save final interval.
    combined.append((intv, r_prev))
    return combined


def create_rankings_ode(db_session: Session, machine: MachineState, methods: List[ODEMethod], ivps: List[IVP],
                        rank_tasks: List[Union[RankDeviationTask, RankOrderTask]], ode_size: Optional[int] = None):
    """Rank impl variants according to the given ranking tasks and store all rankings obtained in database.

    Parameters:
    -----------
    db_session: sqlalchemy.orm.session.Session
        Used database session.
    machine: MachineState
        Used machine.
    methods: List of ODEMethod.
        Used ODE methods.
    ivps: List of IVP
        Used IVPs.
    rank_tasks: List of RankTask
        Used ranking tasks for which ranking will be created and stored in database.
    ode_size: int
        Fixed ODE system size if given.

    Returns:
    --------
    -
    """
    # CPU frequency used.
    frequency: float = machine.clock
    # Max number of CPU cores used.
    max_cores: int = machine.coresPerSocket + 1
    # Retrieve impl variant IDs from database.
    impl_variants: List[int] = ImplVariant.select_ids(db_session)
    # Initialize worker thread pools.
    pool = Pool(cpu_count())
    # Rank variants for ...
    records: List[RankingRecord] = list()
    errors = list()
    # ... all ODE methods
    for method in methods:
        # ... all IVPS:
        for ivp in ivps:
            for task in rank_tasks:
                # Remove old, (possibly obsolete) ranking records from the database.
                if ode_size is not None:
                    # Fixed ODE system size.
                    RankingRecord.remove_record(
                        db_session, machine.db_id, machine.compiler.db_id, method.db_id, ivp.db_id, frequency,
                        [*range(1, max_cores)], ode_size, ode_size, task.type, task.cutoff_criteria, task.cutoff_value)
                else:
                    RankingRecord.remove_records(
                        db_session, machine.db_id, machine.compiler.db_id, method.db_id, ivp.db_id, frequency,
                        [*range(1, max_cores)], task.type, task.cutoff_criteria, task.cutoff_value)
                # ... all core counts.
                for cores in range(1, max_cores):
                    prediction_data = ImplVariantRecord.select(
                        db_session, impl_variants, machine.db_id, machine.compiler.db_id, method.db_id, ivp.db_id,
                        frequency, cores, ode_size)
                    if prediction_data.empty:
                        pool.close()
                        pool.join()
                        raise RuntimeError('No fitting impl variant records found in database.')
                    # Push rank task to worker pool.
                    pool.apply_async(rank_variants_ode,
                                     args=(machine, method, ivp, cores, prediction_data, task, ode_size),
                                     callback=records.extend, error_callback=errors.append)
    # Wait for all threads and collect results.
    pool.close()
    pool.join()
    # Raise error if ranking failed.
    if errors:
        db_session.rollback()
        raise RuntimeError('Failed to rank impl variants: Error in worker threads.')
    # Train database with ranking records.
    RankingRecord.update(db_session, records)


def rank_variants_ode(machine: MachineState, method: ODEMethod, ivp: IVP, cores: int, data: DataFrame,
                      rank_task: Union[RankDeviationTask, RankOrderTask], ode_size: Optional[int] = None):
    records = list()
    try:
        # Fixed ODE system size.
        if ode_size is not None:
            intv = SampleInterval(ode_size, ode_size, None)
            ranking = (intv, rank_task(data, ode_size))
            # Create ranking records.
            records.append(RankingRecord(machine.db_id, machine.compiler.db_id, method.db_id, ivp.db_id, cores,
                                         machine.clock, ranking[0], ranking[1], rank_task.type,
                                         rank_task.cutoff_criteria, rank_task.cutoff_value))
        else:
            # ... split prediction data by impl variant ID.
            # Pass list [idx] here to ensure DataFrame is used instead of Series.
            prediction_data = {idx: data.loc[[idx]] for idx in data.index.unique().values}
            #
            ranking_data = list()
            start = 1
            end = 0
            while any(len(impl.index) > 1 for impl in prediction_data.values()):
                end = sys_maxsize
                # Determine lowest end value of current samples.
                for impl in prediction_data.values():
                    entry = impl.iloc[0]
                    start = entry['first']
                    end = min(entry['last'], end)
                    if start == end:
                        break
                # Add sample interval with this end value to 'rank_data'.
                sample = SampleInterval(start, end, None)
                #
                current_predictions = {idx: impl.iloc[0]['prediction'] for idx, impl in prediction_data.items()}
                current_predictions = DataFrame.from_dict(current_predictions, 'index', columns=['prediction'])
                # Rank variants.
                ranking_data.append((sample, rank_task(current_predictions, end)))
                # Update end of all current samples.
                for impl in prediction_data.values():
                    impl.iat[0, 0] = end + 1
                # Remove all current samples already represented in the list of significant samples.
                for impl_id, impl in prediction_data.items():
                    entry = impl.iloc[0]
                    if len(impl.index) > 1 and entry['first'] > entry['last']:
                        prediction_data[impl_id] = impl[1:]
            # Add interval with 'inf' end.
            sample = SampleInterval(end + 1, sys_maxsize, None)
            #
            current_predictions = {idx: impl.iloc[0]['prediction'] for idx, impl in prediction_data.items()}
            current_predictions = DataFrame.from_dict(current_predictions, 'index', columns=['prediction'])

            # Rank variants.
            ranking_data.append((sample, rank_task(current_predictions, end)))
            # Reduce number of rankings by fusing adjacent rankings if possible.
            ranking_data = fuse_equal_rankings(ranking_data)
            # Create ranking records.
            for ranking in ranking_data:
                records.append(RankingRecord(machine.db_id, machine.compiler.db_id, method.db_id, ivp.db_id, cores,
                                             machine.clock, ranking[0], ranking[1], rank_task.type,
                                             rank_task.cutoff_criteria, rank_task.cutoff_value))
        return records
    except Exception as e:
        # Print stack trace of the executing worker thread.
        print_exc()
        print('')
        raise e


def create_rankings(db_session: Session, machine: MachineState,
                    rank_tasks: List[Union[RankDeviationTask, RankOrderTask]], ode_size: Optional[int] = None):
    """Rank impl variants according to the given ranking tasks and store all rankings obtained in database.

    Parameters:
    -----------
    db_session: sqlalchemy.orm.session.Session
        Used database session.
    machine: MachineState
        Used machine.
    rank_tasks: List of RankTask
        Used ranking tasks for which ranking will be created and stored in database.
    ode_size: int
        Fixed ODE system size if given.

    Returns:
    --------
    -
    """
    # CPU frequency used.
    frequency: float = machine.clock
    # Max number of CPU cores used.
    max_cores: int = machine.coresPerSocket + 1
    # Retrieve impl variant IDs from database.
    impl_variants: List[int] = ImplVariant.select_ids(db_session)
    # Initialize worker thread pools.
    pool = Pool(cpu_count())
    # Rank variants for ...
    records: List[RankingRecord] = list()
    errors = list()
    for task in rank_tasks:
        # Remove old, (possibly obsolete) ranking records from the database.
        if ode_size is not None:
            # Fixed ODE system size.
            RankingRecord.remove_record(
                db_session, machine.db_id, machine.compiler.db_id, -100, -100, frequency, [*range(1, max_cores)],
                ode_size, ode_size, task.type, task.cutoff_criteria, task.cutoff_value)
        else:
            RankingRecord.remove_records(
                db_session, machine.db_id, machine.compiler.db_id, -100, -100, frequency, [*range(1, max_cores)],
                task.type, task.cutoff_criteria, task.cutoff_value)
        # ... all core counts.
        for cores in range(1, max_cores):
            prediction_data = ImplVariantRecord.select(db_session, impl_variants, machine.db_id,
                                                       machine.compiler.db_id, -100, -100, frequency, cores, ode_size)
            if prediction_data.empty:
                pool.close()
                pool.join()
                raise RuntimeError('No fitting impl variant records found in database.')
            # Push rank task to worker pool.
            pool.apply_async(rank_variants, args=(machine, cores, prediction_data, task, ode_size),
                             callback=records.extend, error_callback=errors.append)
    # Wait for all threads and collect results.
    pool.close()
    pool.join()
    # Raise error if ranking failed.
    if errors:
        db_session.rollback()
        raise RuntimeError('Failed to rank impl variants: Error in worker threads.')
    # Train database with ranking records.
    RankingRecord.update(db_session, records)


def rank_variants(machine: MachineState, cores: int, data: DataFrame,
                  rank_task: Union[RankDeviationTask, RankOrderTask], ode_size: Optional[int] = None):
    records = list()
    try:
        # Fixed ODE system size.
        if ode_size is not None:
            intv = SampleInterval(ode_size, ode_size, None)
            ranking = (intv, rank_task(data, ode_size))
            # Create ranking records.
            records.append(RankingRecord(
                machine.db_id, machine.compiler.db_id, -100, -100, cores, machine.clock, ranking[0], ranking[1],
                rank_task.type, rank_task.cutoff_criteria, rank_task.cutoff_value))
        else:
            # ... split prediction data by impl variant ID.
            # Pass list [idx] here to ensure DataFrame is used instead of Series.
            prediction_data = {idx: data.loc[[idx]] for idx in data.index.unique().values}
            #
            ranking_data = list()
            start = 1
            end = 0
            while any(len(impl.index) > 1 for impl in prediction_data.values()):
                end = sys_maxsize
                # Determine lowest end value of current samples.
                for impl in prediction_data.values():
                    entry = impl.iloc[0]
                    start = entry['first']
                    end = min(entry['last'], end)
                    if start == end:
                        break
                # Add sample interval with this end value to 'rank_data'.
                sample = SampleInterval(start, end, None)
                #
                current_predictions = {idx: impl.iloc[0]['prediction'] for idx, impl in prediction_data.items()}
                current_predictions = DataFrame.from_dict(current_predictions, 'index', columns=['prediction'])
                # Rank variants.
                ranking_data.append((sample, rank_task(current_predictions, end)))
                # Update end of all current samples.
                for impl in prediction_data.values():
                    impl.iat[0, 0] = end + 1
                # Remove all current samples already represented in the list of significant samples.
                for impl_id, impl in prediction_data.items():
                    entry = impl.iloc[0]
                    if len(impl.index) > 1 and entry['first'] > entry['last']:
                        prediction_data[impl_id] = impl[1:]
            # Add interval with 'inf' end.
            sample = SampleInterval(end + 1, sys_maxsize, None)
            #
            current_predictions = {idx: impl.iloc[0]['prediction'] for idx, impl in prediction_data.items()}
            current_predictions = DataFrame.from_dict(current_predictions, 'index', columns=['prediction'])
            # Rank variants.
            ranking_data.append((sample, rank_task(current_predictions, end)))
            # Reduce number of rankings by fusing adjacent rankings if possible.
            ranking_data = fuse_equal_rankings(ranking_data)
            # Create ranking records.
            for ranking in ranking_data:
                records.append(RankingRecord(
                    machine.db_id, machine.compiler.db_id, -100, -100, cores, machine.clock, ranking[0], ranking[1],
                    rank_task.type, rank_task.cutoff_criteria, rank_task.cutoff_value))
        return records
    except Exception as e:
        # Print stack trace of the executing worker thread.
        print_exc()
        print('')
        raise e
