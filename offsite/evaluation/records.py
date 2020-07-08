"""@package records
Performance modeling functions.
"""

from datetime import datetime
from getpass import getuser
from sys import maxsize as sys_maxsize
from typing import List, Tuple, Set

import attr
from pandas import read_sql_query
from sqlalchemy import Column, DateTime, Float, ForeignKey, Integer, String, Table, UniqueConstraint
from sqlalchemy.orm import Session
from sqlalchemy.orm.exc import MultipleResultsFound, NoResultFound

from offsite import __version__
from offsite.db import METADATA
from offsite.db.db import insert, bulk_insert
from offsite.descriptions.impl_variant import ImplVariant
from offsite.evaluation.performance_model import SampleInterval


@attr.s
class BenchmarkRecord:
    """Representation of a benchmark table database record.

    Attributes:
    -----------
    name : str
        Name of this object.
    machine : int
        Used machine.
    compiler : int
        Used compiler.
    data : float
        Benchmark result.
    frequency : float
        Used CPU frequency.
    cores : int
        Used number of CPU cores.
    db_id : int
        ID of associated benchmark result database table record.
    """
    name = attr.ib(type=str)
    machine = attr.ib(type=int)
    compiler = attr.ib(type=int)
    data = attr.ib(type=float)
    frequency = attr.ib(type=float)
    cores = attr.ib(type=int)
    db_id = attr.ib(type=int, init=False)

    # Database information.
    db_table = Table('benchmark_result', METADATA,
                     Column('db_id', Integer, primary_key=True),
                     Column('name', String),
                     Column('machine', Integer, ForeignKey('machine.db_id')),
                     Column('compiler', Integer, ForeignKey('compiler.db_id')),
                     Column('data', Float),
                     Column('frequency', Float),
                     Column('cores', Integer),
                     Column('upadtedIn', String, default=__version__),
                     Column('updatedOn', DateTime, default=datetime.now, onupdate=datetime.now),
                     Column('updatedBy', String, default=getuser(), onupdate=getuser()),
                     sqlite_autoincrement=True)

    def to_database(self, db_session: Session):
        """Push this benchmark record object to the database.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        -
        """
        insert(db_session, self)

    @staticmethod
    def contains(db_session: Session, machine: 'Machine', benchmark: 'OmpBarrierBenchmark') -> bool:
        """
        Check if the Benchmark table already contains data of a particular benchmark for a given configuration of
        machine and compiler.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.
        machine : Machine
            Used Machine.
        benchmark : Benchmark
            Used Benchmark.

        Returns:
        --------
        boolean
            True if table contains fitting data False else.
        """
        data = db_session.query(BenchmarkRecord).filter(
            BenchmarkRecord.name.like(benchmark.name), BenchmarkRecord.machine.is_(machine.db_id),
            BenchmarkRecord.compiler.is_(machine.compiler.db_id)).first()
        return bool(data)

    @staticmethod
    def update(db_session: Session, machine: 'Machine', benchmark: 'OmpBarrierBenchmark',
               benchmark_records: List['BenchmarkRecord']):
        """Update data records in Benchmark table with new benchmark data.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.
        machine : Machine
            Used Machine.
        benchmark : Benchmark
            Used Benchmark.
        benchmark_records : list of BenchmarkRecord
            Results of the executed benchmark as list of mathematical expression strings.

        Returns:
        --------
        -
        """
        for record in benchmark_records:
            # Select and remove all already included data records.
            try:
                queried_record = db_session.query(BenchmarkRecord).filter(
                    BenchmarkRecord.name.like(benchmark.name),
                    BenchmarkRecord.machine.is_(machine.db_id),
                    BenchmarkRecord.compiler.is_(machine.compiler.db_id),
                    BenchmarkRecord.frequency.is_(record.frequency),
                    BenchmarkRecord.cores.is_(record.cores)).one()
                # Update data record.
                queried_record.data = (record.data + queried_record.data) / 2
            except NoResultFound:
                # Insert new data record.
                record.to_database(db_session)
            except MultipleResultsFound:
                raise RuntimeError('Unable to update benchmark record!')

    @staticmethod
    def select(db_session: Session, machine: 'Machine', benchmarks: List[str]) -> List['KernelRecord']:
        """Retrieve BenchmarkRecord table data record(s) from the database.

        Return all records that match the provided configuration of machine, compiler and benchmark name(s).

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        machine: Machine
            Used Machine.
        benchmarks : Benchmark
            Names of used benchmarks.

        Returns:
        --------
        pandas.DataFrame
            Retrieved list of data records.
        """
        query = db_session.query(BenchmarkRecord.cores, BenchmarkRecord.name, BenchmarkRecord.data,
                                 BenchmarkRecord.frequency).filter(BenchmarkRecord.machine.is_(machine.db_id),
                                                                   BenchmarkRecord.compiler.is_(machine.compiler.db_id),
                                                                   BenchmarkRecord.name.in_(benchmarks)).order_by(
            BenchmarkRecord.cores, BenchmarkRecord.name)
        data = read_sql_query(query.statement, db_session.bind, index_col='cores')
        return data


@attr.s
class ImplVariantRecord:
    """Representation of an implementation variant prediction table database record.

    Attributes:
    -----------
    impl : int
        Used implementation variant.
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
    first : int
        Start of sample interval.
    last : int
        End of sample interval.
    prediction : str
        Implementation variant's runtime prediction.
    mode : str
        Prediction obtained in RUN or MODEL mode.
    db_id : int
        ID of associated implementation variant prediction database table
        record.
    """
    impl = attr.ib(type=int)
    machine = attr.ib(type=int)
    compiler = attr.ib(type=int)
    method = attr.ib(type=int)
    ivp = attr.ib(type=int)
    cores = attr.ib(type=int)
    frequency = attr.ib(type=float)
    sample = attr.ib(type='SampleInterval')
    first = attr.ib(type=int, init=False)
    last = attr.ib(type=int, init=False)
    prediction = attr.ib(type=str)
    mode = attr.ib(type=str)
    db_id = attr.ib(type=int, init=False)

    # Database information.
    db_table = Table('impl_variant_prediction', METADATA,
                     Column('db_id', Integer, primary_key=True),
                     Column('impl', Integer, ForeignKey('impl_variant.db_id')),
                     Column('machine', Integer, ForeignKey('machine.db_id')),
                     Column('compiler', Integer, ForeignKey('compiler.db_id')),
                     Column('method', Integer, ForeignKey('ode_method.db_id')),
                     Column('ivp', Integer, ForeignKey('ivp.db_id')),
                     Column('frequency', Float),
                     Column('cores', Integer),
                     Column('first', Integer),
                     Column('last', Integer),
                     Column('prediction', String),
                     Column('mode', String),
                     Column('updatedIn', String, default=__version__),
                     Column('updatedOn', DateTime, default=datetime.now, onupdate=datetime.now),
                     Column('updatedBy', String, default=getuser(), onupdate=getuser()),
                     sqlite_autoincrement=True)

    def __attrs_post_init__(self):
        """Create this implementation variant record object.

        Parameters:
        -----------
        -

        Returns:
        --------
        -
        """
        self.first = self.sample.first
        self.last = self.sample.last

    def to_database(self, db_session: Session):
        """Push this implementation variant record object to the database.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        -
        """
        # Attribute prediction.
        self.prediction = str(self.prediction)
        # Insert ImplVariantRecord object.
        insert(db_session, self)

    @staticmethod
    def contains(db_session: Session, impl: 'ImplVariant', machine: int, compiler: int, method: int, ivp: int,
                 cores: int, frequency: float,
                 mode: str) -> bool:
        """
        Check if the ImplVariantRecord table already contains data of a particular implementation variant for a given
        configuration of machine, compiler, ODE method, IVP CPU frequency, and number of CPU cores.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.
        impl : ImplVariant
            Used implementation variant.
        machine : int
            Used machine.
        compiler : int
            Used compiler.
        method : int
            Used ODE method.
        ivp : int
            Used IVP.
        cores : int
            Number of cores.
        frequency : float
            Used CPU frequency.
        mode : str
            Application mode used to obtain data.

        Returns:
        --------
        boolean
            True if table contains fitting data False else.
        """
        data = db_session.query(ImplVariantRecord).filter(
            ImplVariantRecord.impl.is_(impl.db_id), ImplVariantRecord.machine.is_(machine),
            ImplVariantRecord.compiler.is_(compiler), ImplVariantRecord.method.is_(method),
            ImplVariantRecord.ivp.is_(ivp), ImplVariantRecord.cores.is_(cores),
            ImplVariantRecord.frequency.is_(frequency), ImplVariantRecord.mode.is_(mode)).all()
        return bool(data)

    @staticmethod
    def update(db_session: Session, records: List['ImplVariantRecord']):
        """Update data records in ImplVariantPrediction table with new implementation variant prediction data.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.
        impl_variant_records : list
            Impl variant prediction records.
        mode : str
            Application mode used to obtain data.

        Returns:
        --------
        -
        """
        bulk_insert(db_session, records)

    @staticmethod
    def remove_records(db_session: Session, impls: List[int], machine: int, compiler: int, method: int, ivp: int,
                       frequency: float, core_counts: List[int]):
        """Remove implementation variant data record(s) from the database.

         Remove all records that match the provided configuration of machine, compiler, ODE method, IVP, CPU frequency
         and number of CPU cores.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        impl : List of int
            Used implementation variant IDs.
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
        db_session.query(ImplVariantRecord).filter(
            # ImplVariantRecord.impl.in_(impls),
            ImplVariantRecord.machine.is_(machine),
            ImplVariantRecord.compiler.is_(compiler), ImplVariantRecord.method.is_(method),
            ImplVariantRecord.ivp.is_(ivp), ImplVariantRecord.frequency.is_(frequency),
            ImplVariantRecord.cores.in_(core_counts)).delete(synchronize_session=False)

    @staticmethod
    def select(db_session: Session, impls: List[int], machine: int, compiler: int, method: int, ivp: int,
               frequency: float, cores: int) -> 'pandas.DataFrame':
        """Retrieve implementation variant data record(s) from the database.

        Return all records that match the provided configuration of machine, compiler, ODE method, IVP, CPU frequency,
        and number of CPU cores.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        impls : List of int
            Used implementation variant IDs.
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
        cores : int
            Used number of cores.

        Returns:
        --------
        pandas.DataFrame
            Retrieved implementation variant data records.
        """
        query = db_session.query(
            ImplVariantRecord.impl, ImplVariantRecord.first, ImplVariantRecord.last, ImplVariantRecord.prediction). \
            filter(ImplVariantRecord.impl.in_(impls), ImplVariantRecord.machine.is_(machine),
                   ImplVariantRecord.compiler.is_(compiler), ImplVariantRecord.method.is_(method),
                   ImplVariantRecord.ivp.is_(ivp), ImplVariantRecord.frequency.is_(frequency),
                   ImplVariantRecord.cores.is_(cores)). \
            order_by(ImplVariantRecord.impl)
        data = read_sql_query(query.statement, db_session.bind, index_col='impl')
        return data


@attr.s
class KernelRecord:
    """Representation of a kernel runtime prediction table database record.

    Attributes:
    -----------
    kernel : int
        Used kernel.
    machine : int
        Used machine.
    compiler : int
        Used compiler.
    method : int
        Used ODE method.
    ivp : int
        Used IVP.
    frequency : float
        Used CPU frequency.
    cores : int
        Used number of CPU cores.
    sample : SampleInterval
        Used sample interval.
    prediction : float
        Kernel runtime prediction.
    mode : str
        Prediction obtained in RUN or MODEL mode.
    weight : float
        TODO
    db_id : int
        ID of associated kernel prediction database table record.
    """
    kernel = attr.ib(type=int)
    machine = attr.ib(type=int)
    compiler = attr.ib(type=int)
    method = attr.ib(type=int)
    ivp = attr.ib(type=int)
    frequency = attr.ib(type=float)
    cores = attr.ib(type=int)
    sample = attr.ib(type='SampleInterval')
    prediction = attr.ib(type=str)
    mode = attr.ib(type=str)
    weight = attr.ib(type=float, init=False, default=1.0)
    first = attr.ib(type=int, init=False)
    last = attr.ib(type=int, init=False)
    db_id = attr.ib(type=int, init=False)

    # Database information.
    db_table = Table('kernel_prediction', METADATA,
                     Column('db_id', Integer, primary_key=True),
                     Column('kernel', Integer, ForeignKey('kernel.db_id')),
                     Column('machine', Integer, ForeignKey('machine.db_id')),
                     Column('compiler', Integer, ForeignKey('compiler.db_id')),
                     Column('method', Integer, ForeignKey('ode_method.db_id')),
                     Column('ivp', Integer, ForeignKey('ivp.db_id')),
                     Column('frequency', Integer),
                     Column('cores', Integer),
                     Column('first', Integer),
                     Column('last', Integer),
                     Column('prediction', String),
                     Column('mode', String),
                     Column('weight', Float),
                     Column('updatedIn', String, default=__version__),
                     Column('updatedOn', DateTime, default=datetime.now, onupdate=datetime.now),
                     Column('updatedBy', String, default=getuser(), onupdate=getuser()),
                     sqlite_autoincrement=True)

    def __attrs_post_init__(self):
        """Create this ECM Record object.

        Parameters:
        -----------
        -

        Returns:
        --------
        -
        """
        self.first = self.sample.first
        self.last = self.sample.last

    def to_database(self, db_session: Session):
        """Push this ECM record object to the database.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        -
        """
        # Attribute prediction.
        self.prediction = str(self.prediction)
        # Insert KernelRecord object.
        insert(db_session, self)

    @staticmethod
    def update(db_session: Session, kernel: int, machine: int, compiler: int, method: int, ivp: int, cores: int,
               frequency: float, records: List[Tuple[SampleInterval, str]], mode: str):
        """Update data records in PModelRecord table with new data.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.
        machine : int
            Used machine.
        compiler : int
            Used compiler.
        method : int
            Used ODE method.
        ivp : int
            Used IVP.
        cores : int
            Used CPU cores.
        frequency : float
            Used CPU frequency.
        interval : SampleInterval
            Used sample interval.
        records : list of tuple (SampleInterval, kernel prediction)
            Kernel prediction data.
        mode : str
            Application mode used to obtain data.

        Returns:
        --------
        -
        """
        for record in records:
            interval = record[0]
            prediction = record[1]
            # Select and update all already included data records.
            try:
                queried_record = db_session.query(KernelRecord).filter(
                    KernelRecord.kernel.is_(kernel), KernelRecord.machine.is_(machine),
                    KernelRecord.compiler.is_(compiler), KernelRecord.method.is_(method), KernelRecord.ivp.is_(ivp),
                    KernelRecord.cores.is_(cores), KernelRecord.frequency.is_(frequency),
                    KernelRecord.first.is_(interval.first), KernelRecord.last.is_(interval.last),
                    KernelRecord.mode.is_(mode)).one()
                # Update data record.
                queried_record.prediction = str(prediction)
                # queried_record.weight = weight
                # TODO How to handle weight properly?
                # (a) keep value of queriedRecord
                # (b) reset to 1.0
            except NoResultFound:
                # Insert new data record.
                record = KernelRecord(
                    kernel, machine, compiler, method, ivp, frequency, cores, interval, prediction, mode)
                record.to_database(db_session)
            except MultipleResultsFound:
                raise RuntimeError('Unable to update Kernel record!')

    @staticmethod
    def select(db_session: Session, machine: int, compiler: int, method: int, ivp: int, cores: int, frequency: float,
               mode: str, ode_size: int) -> List['KernelRecord']:
        """Retrieve KernelRecord table data record(s) from the database.

        Return all records that match the provided configuration of machine, compiler, ODE method, IVP, CPU frequency
        and number of cores.

        If a ODE system size 'ode_size' is passed only records for that fixed 'ode_size' are returned.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        machine : int
            Used machine.
        compiler : int
            Used compiler.
        method : int
            Used ODE method.
        ivp : int
            Used IVP.
        cores : int
            Used number of cores.
        frequency: float
            Used CPU frequency.
        mode : str
            Application mode used to obtain data.
        ode_size : int
            Used fixed ODE system size.

        Returns:
        --------
        list of KernelRecord
            Retrieved list of data records.
        """
        # Required ivp.db_id = -1 to include all none IVP-dependent kernels.
        if not ode_size:
            data = db_session.query(KernelRecord).filter(
                KernelRecord.machine.is_(machine), KernelRecord.compiler.is_(compiler), KernelRecord.method.is_(method),
                KernelRecord.ivp.in_([ivp, -1]), KernelRecord.frequency.is_(frequency), KernelRecord.cores.is_(cores),
                KernelRecord.mode.is_(mode)).all()
        else:
            data = db_session.query(KernelRecord).filter(
                KernelRecord.machine.is_(machine), KernelRecord.compiler.is_(compiler), KernelRecord.method.is_(method),
                KernelRecord.ivp.in_([ivp, -1]), KernelRecord.frequency.is_(frequency), KernelRecord.cores.is_(cores),
                KernelRecord.mode.is_(mode), KernelRecord.first >= ode_size, KernelRecord.last <= ode_size).all()
        return data

    @staticmethod
    def contains(db_session: Session, kernel: int, machine: int, compiler: int, method: int, ivp: int, frequency: float,
                 max_cores: int, mode: str, ode_size: int) -> bool:
        """
        Check if the KernelRecord table already contains data of a particular kernel for a given configuration of
        machine, compiler, ODE method, IVP, CPU frequency, and number of CPU cores.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.
        kernel : int
            Used kernel.
        machine : int
            Used machine.
        compiler : int
            Used compiler.
        method : int
            Used ODE method.
        ivp : int
            Used IVP.
        frequency : float
            Used CPU frequency.
        max_cores : int
            Maximum number of cores. Check for core counts up to this value.
        mode : str
            Application mode used to obtain data.
        ode_size : int
            Used ODE size or None.

        Returns:
        --------
        boolean
            True if table contains fitting data False else.
        """
        # ODE size is fixed.
        if ode_size:
            data = db_session.query(KernelRecord).filter(
                KernelRecord.kernel.is_(kernel), KernelRecord.machine.is_(machine), KernelRecord.compiler.is_(compiler),
                KernelRecord.method.is_(method), KernelRecord.ivp.is_(ivp), KernelRecord.frequency.is_(frequency),
                KernelRecord.cores.in_((c for c in range(1, max_cores + 1))), KernelRecord.first.is_(ode_size),
                KernelRecord.last.is_(ode_size), KernelRecord.mode.is_(mode)).all()
            return bool(len(data) == max_cores)
        # ODE size is not fixed.
        data = db_session.query(KernelRecord).filter(KernelRecord.kernel.is_(kernel), KernelRecord.machine.is_(machine),
                                                     KernelRecord.compiler.is_(compiler),
                                                     KernelRecord.method.is_(method), KernelRecord.ivp.is_(ivp),
                                                     KernelRecord.frequency.is_(frequency),
                                                     KernelRecord.cores.in_((c for c in range(1, max_cores + 1))),
                                                     KernelRecord.mode.is_(mode)).order_by(KernelRecord.cores.asc(),
                                                                                           KernelRecord.first.asc()).all()
        if not data:
            return False
        # Check if data for all possible ODE sizes is available.
        cur_core = 0
        cur_last = sys_maxsize
        for record in data:
            if record.cores != cur_core:
                if cur_last != sys_maxsize:
                    return False
                cur_core = record.cores
                if record.first != 1:
                    return False
                cur_last = record.last
            else:
                if (cur_last + 1) != record.first:
                    return False
                cur_last = record.last
        assert cur_core == max_cores
        return bool(cur_last == sys_maxsize)


@attr.s
class RankingRecord:
    """
    Representation of an implementation variant ranking database table record.

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
                     Column('updatedIn', String, default=__version__),
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
