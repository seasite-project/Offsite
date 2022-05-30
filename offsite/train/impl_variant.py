"""@package train.node.records.impl_variant
Definition of classes ImplVariant and ImplVariantRecord.

@author: Johannes Seiferth
"""

from datetime import datetime
from getpass import getuser
from multiprocessing import cpu_count, Pool
from traceback import print_exc
from typing import Dict, List, Optional, Tuple

import attr
from pandas import read_sql_query, DataFrame, Series
from sqlalchemy import Column, DateTime, Float, ForeignKey, Integer, String, Table, UniqueConstraint
from sqlalchemy.exc import NoResultFound, MultipleResultsFound
from sqlalchemy.orm import Query, Session
from sympy import simplify

from offsite import __version__
from offsite.codegen.codegen_util import create_variant_name
from offsite.database import METADATA, insert
from offsite.descriptions.impl.impl_skeleton import ImplSkeleton
from offsite.descriptions.impl.kernel_template import Kernel
from offsite.util.math_utils import eval_math_expr
from offsite.util.sample_interval import SampleInterval

KernelDict = Dict[int, Kernel]
StringDict = Dict[str, str]


@attr.s
class ImplVariant:
    """Representation of an ImplVariant table database record.

    Attributes:
    -----------
    skeleton: int
        Used implementation skeleton.
    kernels: list of int
        Used machine.
    db_id: int
        ID of associated implementation variant database table record.
    """
    skeleton = attr.ib(type=int)
    kernels = attr.ib(type=List[int])
    kernel_communication_cluster_lvl = attr.ib(type=Dict[str, str])
    kernel_communication_node_lvl = attr.ib(type=Dict[str, str])
    kernels_serial = attr.ib(type=int, init=False)
    db_id = attr.ib(type=int, init=False)

    # Database information.
    db_table = Table('impl_variant', METADATA,
                     Column('db_id', Integer, primary_key=True),
                     Column('skeleton', Integer, ForeignKey('impl_skeleton.db_id')),
                     Column('kernels_serial', Integer),
                     Column('updatedIn', String, default=__version__),
                     Column('updatedOn', DateTime, default=datetime.now, onupdate=datetime.now),
                     Column('updatedBy', String, default=getuser(), onupdate=getuser()),
                     UniqueConstraint('skeleton', 'kernels_serial'),
                     sqlite_autoincrement=True)

    def to_database(self, db_session: Session) -> 'ImplVariant':
        """Push this impl variant object to the database.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        ImplVariant
            Instance of this object connected to database session.
        """
        # Attribute kernels_serial.
        self.kernels_serial = ','.join(map(str, self.kernels))
        # Check if database already contains this ImplVariant object.
        variant: ImplVariant = db_session.query(ImplVariant).filter(
            ImplVariant.skeleton.is_(self.skeleton), ImplVariant.kernels_serial.is_(self.kernels_serial)).one_or_none()
        if variant:
            # Supplement attributes not saved in database.
            variant.kernels = self.kernels
            # Update serialized members.
            variant.kernels_serial = ','.join(map(str, variant.kernels))
            # Attribute kernel_communication.
            clust_lvl_comm, node_lvl_comm, _ = ImplVariant.count_kernel_communication(db_session, variant.kernels, {})
            variant.kernel_communication_cluster_lvl = clust_lvl_comm
            variant.kernel_communication_node_lvl = node_lvl_comm
            return variant
        # Add new object to database.
        # Attribute kernels_serial.
        self.kernels_serial = ','.join(map(str, self.kernels))
        # Insert ImplVariantRecord object.
        insert(db_session, self)
        return self

    @staticmethod
    def select(db_session: Session, variant_ids: List[int]) -> List['ImplVariant']:
        """
        Retrieve the ImplVariant table data record(s) from the database that match the given implementation variant IDs.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        variant_ids: list of int
            IDs of the implementation variants requested.

        Returns:
        --------
        list of ImplVariant
            Retrieved list of data records.
        """
        variants: List[ImplVariant] = db_session.query(ImplVariant).filter(ImplVariant.db_id.in_(variant_ids)).all()
        if len(variants) != len(variant_ids):
            raise RuntimeError('Unable to select all requested variants!')
        # Attributes not stored in database....
        loaded_kernels: KernelDict = dict()
        for variant in variants:
            # Attribute kernels.
            variant.kernels = list(map(int, variant.kernels_serial.split(',')))
            # Attribute clust/node_lvl_comm.
            clust_lvl_comm, node_lvl_comm, loaded_kernels = ImplVariant.count_kernel_communication(
                db_session, variant.kernels, loaded_kernels)
            variant.kernel_communication_cluster_lvl = clust_lvl_comm
            variant.kernel_communication_node_lvl = node_lvl_comm
        return variants

    @staticmethod
    def select_ids(db_session: Session) -> List[int]:
        """Retrieve all available ImplVariant IDs from the database.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        list of int
            Retrieved list of ImplVariant IDs.
        """
        data = db_session.query(ImplVariant.db_id).all()
        ids: List[int] = [d[0] for d in data]
        return ids

    @staticmethod
    def count_kernel_communication(db_session: Session, kernels: List[int],
                                   loaded_kernels: KernelDict) -> Tuple[StringDict, StringDict, KernelDict]:
        clust_lvl_comm: StringDict = dict()
        node_lvl_comm: StringDict = dict()
        for kid in kernels:
            if kid not in loaded_kernels:
                loaded_kernels[kid] = Kernel.select(db_session, kid)
            kernel: Kernel = loaded_kernels[kid]
            #
            for op, excs in kernel.communicationOperationsClusterLvl.items():
                if op not in clust_lvl_comm:
                    clust_lvl_comm[op] = excs
                else:
                    clust_lvl_comm[op] = eval_math_expr('{} + {}'.format(clust_lvl_comm[op], excs))
            for op, excs in kernel.communicationOperationsNodeLvl.items():
                if op not in node_lvl_comm:
                    node_lvl_comm[op] = excs
                else:
                    node_lvl_comm[op] = eval_math_expr('{} + {}'.format(node_lvl_comm[op], excs))
        return clust_lvl_comm, node_lvl_comm, loaded_kernels

    @staticmethod
    def fetch_impl_variant_name(db_session: Session, db_id: int) -> str:
        # Get implementation variant record first.
        variant: ImplVariant = ImplVariant.select(db_session, [db_id])[0]
        # Next, get records of all kernels used in the variant as well as the name of the variant's skeleton.
        try:
            skeleton: str = db_session.query(ImplSkeleton.name).filter(ImplSkeleton.db_id.is_(variant.skeleton)).one()[
                0]
        except NoResultFound:
            raise RuntimeError('Unable to load ImplSkeleton object from database!')
        except MultipleResultsFound:
            raise RuntimeError('Unable to load unique ImplSkeleton object from database!')
        # ... kernel records.
        kernels: List[Kernel] = [Kernel.select(db_session, kid) for kid in variant.kernels]
        # Derive implementation name from given information.
        name: str = create_variant_name(kernels, skeleton)
        return name


@attr.s
class ImplVariantRecord:
    """Representation of an implementation variant prediction table database record.

    Attributes:
    -----------
    impl: int
        Used implementation variant.
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
    first: int
        Start of sample interval.
    last: int
        End of sample interval.
    prediction: str
        Implementation variant's runtime prediction.
    mode: str
        Prediction obtained in RUN or MODEL mode.
    db_id: int
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
    sample = attr.ib(type=SampleInterval)
    first = attr.ib(type=int, init=False)
    last = attr.ib(type=int, init=False)
    prediction = attr.ib(type=str)
    pred_mode = attr.ib(type=str)
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
        db_session: sqlalchemy.orm.session.Session
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
    def contains(db_session: Session, impl: ImplVariant, machine: int, compiler: int, method: int, ivp: int, cores: int,
                 frequency: float, mode: str) -> bool:
        """
        Check if the ImplVariantRecord table already contains data of a particular implementation variant for a given
        configuration of machine, compiler, ODE method, IVP CPU frequency, and number of CPU cores.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        impl: ImplVariant
            Used impl variant.
        machine: int
            Used machine.
        compiler: int
            Used compiler.
        method: int
            Used ODE method.
        ivp: int
            Used IVP.
        cores: int
            Number of cores.
        frequency: float
            Used CPU frequency.
        mode: str
            Application mode used to obtain data.

        Returns:
        --------
        boolean
            True if table contains fitting data False else.
        """
        return bool(db_session.query(ImplVariantRecord).filter(
            ImplVariantRecord.impl.is_(impl.db_id), ImplVariantRecord.machine.is_(machine),
            ImplVariantRecord.compiler.is_(compiler), ImplVariantRecord.method.is_(method),
            ImplVariantRecord.ivp.is_(ivp), ImplVariantRecord.cores.is_(cores),
            ImplVariantRecord.frequency.is_(frequency), ImplVariantRecord.pred_mode.is_(mode)).all())

    @staticmethod
    def update(db_session: Session, records: DataFrame):
        """Update data records in ImplVariantPrediction table with new impl variant prediction data.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        records: DataFrame
            Impl variant prediction records.

        Returns:
        --------
        -
        """
        records.to_sql(name='impl_variant_prediction', con=db_session.connection(), if_exists='append', index=False)

    @staticmethod
    def remove_records(db_session: Session, impls: List[int], machine: int, compiler: int, method: int, ivp: int,
                       frequency: float, core_counts: List[int]):
        """Remove impl variant data record(s) from the database.

         Remove all records that match the provided configuration of machine, compiler, ODE method, IVP, CPU frequency
         and number of CPU cores.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        impl: List of int
            Used impl variant IDs.
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

        Returns:
        --------
        -
        """
        db_session.query(ImplVariantRecord).filter(
            ImplVariantRecord.impl.in_(impls), ImplVariantRecord.machine.is_(machine),
            ImplVariantRecord.compiler.is_(compiler), ImplVariantRecord.method.is_(method),
            ImplVariantRecord.ivp.is_(ivp), ImplVariantRecord.frequency.is_(frequency),
            ImplVariantRecord.cores.in_(core_counts)).delete(synchronize_session=False)

    @staticmethod
    def select(db_session: Session, impls: List[int], machine: int, compiler: int, method: int, ivp: int,
               frequency: float, cores: int, ode_size: int) -> DataFrame:
        """Retrieve impl variant data record(s) from the database.

        Return all records that match the provided configuration of machine, compiler, ODE method, IVP, CPU frequency,
        and number of CPU cores.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        impls: List of int
            Used impl variant IDs.
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
        cores: int
            Used number of cores.
        ode_size: int
            Used ODE system size.

        Returns:
        --------
        pandas.DataFrame
            Retrieved impl variant data records.
        """
        query: Query
        if ode_size:
            query: Query = db_session.query(
                ImplVariantRecord.impl, ImplVariantRecord.first, ImplVariantRecord.last,
                ImplVariantRecord.prediction).filter(
                ImplVariantRecord.impl.in_(impls), ImplVariantRecord.machine.is_(machine),
                ImplVariantRecord.compiler.is_(compiler), ImplVariantRecord.method.is_(method),
                ImplVariantRecord.ivp.is_(ivp), ImplVariantRecord.frequency.is_(frequency),
                ImplVariantRecord.cores.is_(cores), ImplVariantRecord.first.is_(ode_size),
                ImplVariantRecord.last.is_(ode_size)).order_by(ImplVariantRecord.impl)
        else:
            query: Query = db_session.query(
                ImplVariantRecord.impl, ImplVariantRecord.first, ImplVariantRecord.last,
                ImplVariantRecord.prediction).filter(
                ImplVariantRecord.impl.in_(impls), ImplVariantRecord.machine.is_(machine),
                ImplVariantRecord.compiler.is_(compiler), ImplVariantRecord.method.is_(method),
                ImplVariantRecord.ivp.is_(ivp), ImplVariantRecord.frequency.is_(frequency),
                ImplVariantRecord.cores.is_(cores)).order_by(ImplVariantRecord.impl)
        return read_sql_query(query.statement, db_session.bind, index_col='impl')


def fuse_equal_impl_records(records: DataFrame) -> DataFrame:
    """
    Reduce total number of records by combining adjacent intervals. Adjacent intervals can be combined if they give the
    same prediction.

    Parameters:
    -----------
    records: DataFrame
        Impl variant records.

    Returns:
    --------
    DataFrame
        Reduced set of impl variant records.
    """
    if records.empty:
        return records
    fused_records: List[Dict] = list()
    errors: List[str] = list()
    # Initialize worker thread pool.
    pool: Pool = Pool(max(cpu_count() - 1, 1))
    # Fuse neighbouring intervals that give the same prediction.
    groups_impl = records.groupby(records.impl)
    for group_impl in groups_impl:
        recs = groups_impl.get_group(group_impl[0])
        pool.apply_async(  # comma after args needed here!
            __fuse_equal_impl_records, callback=fused_records.extend, error_callback=errors.append, args=(recs,))
    # Wait for all threads and collect results.
    pool.close()
    pool.join()
    # Raise error if failed.
    if errors:
        raise RuntimeError('Failed to fuse records: Error in worker threads.')
    fused_records_df = DataFrame(fused_records)
    return fused_records_df


def __fuse_equal_impl_records(records: DataFrame) -> List[Dict]:
    def __series_to_record(ser: Series):
        return {'impl': int(ser.get('impl')), 'machine': int(ser.get('machine')), 'compiler': int(ser.get('compiler')),
                'method': int(ser.get('method')), 'ivp': int(ser.get('ivp')), 'cores': int(ser.get('cores')),
                'frequency': float(ser.get('frequency')), 'first': int(ser.get('first')), 'last': int(ser.get('last')),
                'prediction': str(ser.get('prediction')), 'mode': str(ser.get('mode'))}

    fused_records: List[Dict] = list()
    try:
        for group_cores in records.groupby(['cores']):
            records = group_cores[1].sort_values(by=['first'])
            #
            preds = records.prediction.to_numpy()
            if (preds[0] != preds).any():
                cur: Optional[Dict] = None
                for _, rec in records.iterrows():
                    rec: Series
                    # Test if both intervals have the same prediction.
                    if cur is not None and simplify(cur['prediction']) == simplify(rec['prediction']):
                        # Increase upper bound of current interval.
                        cur['last'] = rec.get('last')
                    else:
                        # Save interval.
                        if cur is not None:
                            fused_records.append(cur)
                        # Update current record, result and impl.
                        cur = __series_to_record(rec)
                # Save final interval.
                assert cur is not None
                fused_records.append(cur)
            else:
                rec = __series_to_record(records.iloc[0])
                rec['last'] = int(records.iloc[-1].get('last'))
                fused_records.append(rec)
        return fused_records
    except Exception as e:
        # Print stack trace of the executing worker thread.
        print_exc()
        print('')
        raise e
