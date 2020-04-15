"""@package ranking
Definitions of ranking functions.
"""

from datetime import datetime
from getpass import getuser
from multiprocessing import cpu_count, Pool
from sys import maxsize as sys_maxsize
from typing import List, Tuple, Set

import attr
from pandas import DataFrame
from sqlalchemy import Column, DateTime, Float, ForeignKey, Integer, String, Table, UniqueConstraint
from sqlalchemy.orm import Session

import offsite.config
from offsite import __version__
from offsite.config import ProgramModeType
from offsite.db import METADATA
from offsite.db.db import insert, bulk_insert
from offsite.descriptions.impl_skeleton import ImplVariant
from offsite.descriptions.ivp import IVP
from offsite.descriptions.machine import Machine
from offsite.descriptions.ode_method import ODEMethod
from offsite.evaluation.math_utils import eval_math_expr, ivp_system_size, percent_deviation
from offsite.evaluation.performance_model import ImplVariantRecord, SampleInterval

RankingData = List[Tuple[SampleInterval, Set[int]]]


def deduce_best_variants(prediction_data: 'pandas.DataFrame', system_size: float, tolerance: float = 5.0) -> Set[int]:
    """
    Rank implementation variants by ascending runtime prediction and return the set of ids of those implementation
    variants that are within a given tolerance from the best predicted implementation variant.

    Parameters:
    -----------
    data : dict (key: impl variant id, value: ImplVariantRecord)
        Implementation variant prediction data.
    system_size : float
        Rank variants for this system size.
    tolerance : float
        Tolerance.

    Returns:
    -------
    Set of int
        Set of implementation variant ids ranked by ascending runtime prediction.
    """
    constants = [ivp_system_size(system_size), ('x', system_size)]
    # Evaluate predictions for the given ODE system size.
    for idx, row in prediction_data.iterrows():
        prediction_data.at[idx, 'prediction'] = eval_math_expr(row[2], constants, cast_to=float)
    # Rank implementation variants by ascending runtime prediction.
    ranking = prediction_data.sort_values(by=['prediction'])
    # Return all implementation variants within the tolerance.
    index_best_impl = ranking.head(1).index.values[0]
    return set([idx for idx, impl in ranking.iterrows() \
                if (percent_deviation(impl['prediction'], ranking.at[index_best_impl, 'prediction'])) <= tolerance])


def fuse_equal_rankings(rankings: RankingData) -> RankingData:
    """
    Reduce total number of rankings by combining neighbouring rankings. Adjacent rankings can be combined if they
    contain the same implementation variants.

    Parameters:
    -----------
    rankings : list of tuple(SampleInterval, set of impl variant ids)
        Rankings of implementation variants for ascending intervals.

    Returns:
    --------
    list of tuple(SampleInterval, set of impl variant ids)
        Rankings of implementation variants for ascending intervals after fusing.
    """
    if len(rankings) <= 1:
        return rankings
    # Iterate over the rankings and combine those that contain the same implementation variants.
    combined = list()
    r_prev = set()
    intv = None
    for ranking in rankings:
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
            intv = SampleInterval(i_cur.first, i_cur.last)
    # Save final interval.
    combined.append((intv, r_prev))
    return combined


def callback_error_message(error):
    print('Error: {}'.format(error))


def rank(db_session: Session, machine: Machine, methods: List[ODEMethod], ivps: List[IVP]):
    """
    Rank implementation variants by ascending runtime prediction and return the set of ids of those implementation
    variants that are within a given tolerance from the best predicted implementation variant.

    Parameters:
    -----------
    config : Config
        Program configuration.
    db_session : sqlalchemy.orm.session.Session
        Used database session.
    machine : Machine
        Used Machine.
    methods : List of ODEMethod.
        Used ODE methods.
    ivps : List of IVP
        Used IVPs.

    Returns:
    --------
    -
    """
    config = offsite.config.offsiteConfig
    args = config.args
    if args.mode == ProgramModeType.RUN:
        mode = 'RUN'
    elif args.mode == ProgramModeType.MODEL:
        mode = 'MODEL'
    else:
        assert False
    # Machine frequency used.
    frequency = machine.clock
    # Max number of CPU cores used.
    max_cores = machine.coresPerSocket + 1
    # Retrieve implementation variants from database.
    impl_variants = ImplVariant.select_all(db_session)
    # Initialize worker thread pools.
    pool = Pool(cpu_count())
    # Rank variants for ...
    records = list()
    # ... all ODE methods
    for method in methods:
        # ... all IVPS:
        for ivp in ivps:
            if args.tool is not ivp.modelTool:
                continue
            # Remove old, (possibly obsolete) ranking records from the database.
            # Fixed ODE system size.
            if config.args.ode_size:
                RankingRecord.remove_record(db_session, machine.db_id, machine.compiler.db_id, method.db_id,
                                            ivp.db_id, frequency, [*range(1, max_cores)], config.args.ode_size,
                                            config.args.ode_size)
            else:
                RankingRecord.remove_records(db_session, machine.db_id, machine.compiler.db_id, method.db_id,
                                             ivp.db_id, frequency, [*range(1, max_cores)])
            # ... all core counts.
            for cores in range(1, max_cores):
                prediction_data = ImplVariantRecord.select(db_session, [impl.db_id for impl in impl_variants],
                                                           machine.db_id, machine.compiler.db_id,
                                                           method.db_id, ivp.db_id, frequency, cores, mode)
                # Push rank task to worker pool.
                pool.apply_async(rank_variants, args=(machine, method, ivp, cores, prediction_data),
                                 callback=records.extend, error_callback=callback_error_message)
    # Wait for all threads and collect results.
    pool.close()
    pool.join()
    # Train database with ranking records.
    RankingRecord.update(db_session, records)


def rank_variants(machine: Machine, method: ODEMethod, ivp: IVP, cores: int, data: 'pandas.DataFrame'):
    config = offsite.config.offsiteConfig
    records = list()
    # Fixed ODE system size.
    if config.args.ode_size:
        ranking = (SampleInterval(config.args.ode_size, config.args.ode_size),
                   deduce_best_variants(data, config.args.ode_size, config.args.tolerance))
        # Create ranking records.
        records.append(RankingRecord(machine.db_id, machine.compiler.db_id, method.db_id, ivp.db_id, cores,
                                     machine.clock, ranking[0], ranking[1]))
    else:
        # ... split prediction data by implementation variant ID.
        prediction_data = {idx: data.loc[idx] for idx in data.index.unique().values}
        #
        ranking_data = list()
        start = 1
        end = 0
        while any(len(impl) > 1 for impl in prediction_data.values()):
            end = sys_maxsize
            # Determine lowest end value of current samples.
            low_sample = None
            for impl in prediction_data.values():
                first_row = [x for x in impl.head(1).iterrows()][0][1]
                start = first_row['first']
                try:
                    low_sample = first_row if first_row['last'] < end else low_sample
                    end = min(first_row['last'], end)
                except TypeError:
                    low_sample = first_row
                    end = first_row['last'] if end == 'inf' else end
                if start == end:
                    break
            # Add sample interval with this end value to 'rank_data'.
            sample = SampleInterval(start, end)
            #
            current_predictions = dict()
            for idx, impl in prediction_data.items():
                current_predictions[idx] = [x for x in impl.head(1).iterrows()][0][1]
            current_predictions = DataFrame.from_dict(current_predictions, 'index')
            # Rank variants.
            ranking_data.append((sample, deduce_best_variants(current_predictions, end, config.args.tolerance)))
            # Update end of all current samples.
            for impl in prediction_data.values():
                impl.iat[0, 0] = end + 1
            # Remove all current samples already represented in the list of significant samples.
            for impl_id, impl in prediction_data.items():
                first_row = [x for x in impl.head(1).iterrows()][0][1]
                if len(impl) > 1 and (first_row['last'] == 'inf' or first_row['first'] > first_row['last']):
                    prediction_data[impl_id] = impl[1:]
        # Add interval with 'inf' end.
        sample = SampleInterval(end + 1, 'inf')
        #
        current_predictions = dict()
        for idx, impl in prediction_data.items():
            current_predictions[idx] = [x for x in impl.head(1).iterrows()][0][1]
        current_predictions = DataFrame.from_dict(current_predictions, 'index')
        # Rank variants.
        ranking_data.append((sample, deduce_best_variants(current_predictions, end, config.args.tolerance)))
        # Reduce number of rankings by fusing adjacent rankings if possible.
        ranking_data = fuse_equal_rankings(ranking_data)
        # Create ranking records.
        for ranking in ranking_data:
            records.append(RankingRecord(machine.db_id, machine.compiler.db_id, method.db_id, ivp.db_id, cores,
                                         machine.clock, ranking[0], ranking[1]))
    return records


@attr.s
class RankingRecord:
    """
    Representation of an implementation variant prediction table database record.

    Attributes:
    -----------
    machine : int
        Used machine.
    compiler : int
        Used compiler.
    method : int
        Used ODE method.
    ivp : int
        Used IVP.
    cores: int
        Used number of cores.
    frequency : float
        Used CPU frequency.
    sample : SampleInterval
        Used sample interval.
    variants : Set of int
        ID's of implementation variants.
    db_id : int
        ID of associated ranking database table record.
    """
    machine = attr.ib(type=int)
    compiler = attr.ib(type=int)
    method = attr.ib(type=int)
    ivp = attr.ib(type=int)
    cores = attr.ib(type=int)
    frequency = attr.ib(type=float)
    sample = attr.ib(type='SampleInterval')
    variants = attr.ib(type=Set[int])
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
                     Column('createdBy', String, default=getuser()),
                     Column('createdIn', String, default=__version__),
                     Column('createdOn', DateTime, default=datetime.now),
                     Column('updatedOn', DateTime, default=datetime.now, onupdate=datetime.now),
                     Column('updatedBy', String, default=getuser(), onupdate=getuser()),
                     UniqueConstraint('machine', 'compiler', 'method', 'ivp', 'frequency', 'cores', 'first', 'last'),
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
        db_session : sqlalchemy.orm.session.Session
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
                      core_counts: List[int], first: int, last: int):
        """Remove a single ranking data record from the database.

        Remove the record that matches the provided configuration of machine, compiler, ODE method, IVP, CPU frequency,
        number of CPU cores and ODE system size.

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
        frequency:
            Used CPU frequency.
        core_counts : List of int
            Used core counts.

        Returns:
        --------
        -
        """
        db_session.query(RankingRecord).filter(
            RankingRecord.machine.is_(machine), RankingRecord.compiler.is_(compiler), RankingRecord.method.is_(method),
            RankingRecord.ivp.is_(ivp), RankingRecord.frequency.is_(frequency), RankingRecord.cores.in_(core_counts),
            RankingRecord.first.is_(first), RankingRecord.last.is_(last)).delete(synchronize_session=False)

    @staticmethod
    def remove_records(db_session: Session, machine: int, compiler: int, method: int, ivp: int, frequency: float,
                       core_counts: List[int]):
        """Remove ranking data record(s) from the database.

        Remove all records that match the provided configuration of machine, compiler, ODE method, IVP, CPU frequency
        and number of CPU cores.

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
        frequency:
            Used CPU frequency.
        core_counts : List of int
            Used core counts.

        Returns:
        --------
        -
        """
        db_session.query(RankingRecord).filter(
            RankingRecord.machine.is_(machine), RankingRecord.compiler.is_(compiler), RankingRecord.method.is_(method),
            RankingRecord.ivp.is_(ivp), RankingRecord.frequency.is_(frequency),
            RankingRecord.cores.in_(core_counts)).delete(synchronize_session=False)
