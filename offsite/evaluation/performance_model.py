"""@package performance_model
Performance modeling functions.
"""

from datetime import datetime
from enum import Enum
from getpass import getuser
from sys import maxsize as sys_maxsize
from typing import Dict, List, Tuple

import attr
from pandas import read_sql_query
from sqlalchemy import Column, DateTime, Float, ForeignKey, Integer, String, Table
from sqlalchemy.orm import Session
from sqlalchemy.orm.exc import MultipleResultsFound, NoResultFound

from offsite import __version__
from offsite.db.db import METADATA
from offsite.db.db import insert, bulk_insert
from offsite.descriptions.ivp import IVP
from offsite.descriptions.machine import Machine
from offsite.evaluation.math_utils import eval_math_expr, ivp_grid_size, ivp_system_size


def compute_pmodel_kernel_pred(num_iterations: int, machine: Machine, ecm: float) -> float:
    """Compute the pmodel kernel prediction in cycles of a PModelKernel.

    The kernel prediction states the total number of cycles executed when executing a given number of iterations of
    that pmodel kernel.

    Parameters:
    -----------
    num_iterations : int
        Number of iterations executed by the PModelKernel considered.
    machine : Machine
        Machine the kernel is run on.
    ecm : float
        ECM result (in cycles per cache line) obtained by the kerncraft tool for the PModelKernel considered when
        executing 'num_iterations' iterations.

    Returns:
    --------
    float
        Kernel prediction in cycles of a PModelKernel.
    """
    constants = [('ecm', ecm), ('iterations', num_iterations), ('elements_cacheline', machine.elements_per_cacheline)]
    return eval_math_expr('ecm * iterations / elements_cacheline', constants)


def compute_kernel_runtime_pred(pmodel_predictions: List[str], frequency: float) -> float:
    """Compute the kernel runtime prediction in seconds of a Kernel.

    The kernel runtime prediction states the runtime in seconds when executing that kernel with a particular CPU
    frequency.

    Parameters:
    -----------
    pmodel_predictions : list of str
        PModel kernel predictions obtained for the pmodel kernels associated with this Kernel when running on a
        particular number of CPU cores.
    frequency : float
        CPU frequency the kernel is executed with.

    Returns:
    --------
    float
        Kernel runtime prediction in seconds.
    """
    # Construct equation.
    known_variables = [('f', frequency)]
    equation_str = ''
    for i, prediction in enumerate(pmodel_predictions):
        equation_str += '+ p{0}'.format(i)
        known_variables.extend([('p{}'.format(i), prediction)])
    equation_str = '({})/f'.format(equation_str)
    return eval_math_expr(equation_str, known_variables)


def compute_impl_variant_runtime_pred(associated_kernels: Tuple[int], kernel_runtime_preds: Dict[int, float],
                                      executions: Dict[int, float], communication_costs: float) -> float:
    """Compute the node-level runtime prediction in seconds of an implementation variant.

    The node-level runtime prediction states the runtime in seconds when executing a single time step of that
    implementation variant. The prediction is computed using the kernel runtime predictions of its associated kernels
    as well as the communication costs of this implementation variant.

    Parameters:
    -----------
    associated_kernels : tuple of int
        Ids of the Kernel objects associated with this implementation variant.
    kernel_runtime_preds : dict (key: Kernel object id)
        Kernel runtime predictions obtained for the kernels associated with this implementation variant.
    executions : dict (key: Kernel object id)
        Number of times the kernels associated with this implementation variant are run when executing a single time
        step.
    communication_costs : float
        Communication costs per step of this implementation variant.

    Returns:
    --------
    float
        Node-level runtime prediction in seconds.
    """
    # Compute the node-level runtime prediction.
    node_level_prediction = eval_math_expr(0.0)
    # Add up the single kernel runtime predictions.
    for kernel in associated_kernels:
        try:
            expr = '{} * {} * {}'.format(
                kernel_runtime_preds[kernel][0], kernel_runtime_preds[kernel][1], executions[kernel])
            node_level_prediction += eval_math_expr(expr)
        except KeyError:
            raise RuntimeError('Failed to compute impl variant runtime prediction!')
    # Add communication costs.
    node_level_prediction += communication_costs
    return node_level_prediction


class SamplePosition(Enum):
    """
    Defines if the sample lies in the border region of two adjacent intervals or not.

    - INNER
        Sample interval in the inner region of an interval.
    - BORDER
        Sample interval in the border region of two adjacent intervals.
    """
    INNER = 'INNER'
    BORDER = 'BORDER'


@attr.s(hash=True)
class SampleInterval:
    """Representation of a SampleInterval object.

    Attributes:
    -----------
    first : int
        First value included in the sample interval.
    last : int
        First value included in the sample interval.
    sample : int
        Actual value used to sample the interval.
    position : SamplePosition
        Type of the interval. Can be either 'INNER' or 'BORDER'.
    """
    first = attr.ib(type=int)
    last = attr.ib(type=int)
    sample = attr.ib(type=int, default=None)
    region = attr.ib(type=SamplePosition, default=SamplePosition.INNER)

    def __hash__(self):
        return hash((self.first, self.last))

    def __eq__(self, other):
        return (self.first, self.last) == (other.first, other.last)

    def __ne__(self, other):
        return not self == other

    def median(self, ivp: IVP = None) -> int:
        """Return median value of this interval.

        If an IVP is passed not the exact median is returned but the nearest value that satisfies the constraints of
        the IVP regarding possible system sizes(e.g. square system size for Heat2D, ...).

        Parameters:
        -----------
        ivp: IVP
            Used IVP.

        Returns:
        --------
        int
            Median value of this interval.
        """
        median = int(round(self.first + (self.last - self.first) / 2))
        # Consider restraints of the grid size of the IVP constraining possible but also sampleable, median values.
        if ivp:
            # Round up to previous ...
            grid_size = int(eval_math_expr(ivp.gridSize, [ivp_system_size(median)]))
            median_low = int(eval_math_expr(ivp.ivp_size(), [ivp_grid_size(grid_size)]))
            if median_low < self.first or median_low > self.last:
                # Round up to next ...
                grid_size = int(round(eval_math_expr(ivp.gridSize, [ivp_system_size(median)])))
                median_high = int(round(eval_math_expr(ivp.ivp_size(), [ivp_grid_size(grid_size)])))
                if median_high < self.first or median_high > self.last:
                    median = None
                else:
                    median = median_high
            else:
                median = median_low
        return median


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
                     Column('createdBy', String, default=getuser()),
                     Column('createdIn', String, default=__version__),
                     Column('createdOn', DateTime, default=datetime.now),
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
                KernelRecord.cores.in_([c for c in range(1, max_cores + 1)]), KernelRecord.first.is_(ode_size),
                KernelRecord.last.is_(ode_size), KernelRecord.mode.is_(mode)).all()
            return bool(len(data) == max_cores)
        # ODE size is not fixed.
        data = db_session.query(KernelRecord).filter(
            KernelRecord.kernel.is_(kernel), KernelRecord.machine.is_(machine), KernelRecord.compiler.is_(compiler),
            KernelRecord.method.is_(method), KernelRecord.ivp.is_(ivp), KernelRecord.frequency.is_(frequency),
            KernelRecord.cores.in_([c for c in range(1, max_cores + 1)]), KernelRecord.mode.is_(mode)).order_by(
            KernelRecord.cores.asc(), KernelRecord.first.asc()).all()
        if not data:
            return False
        # Check if data for all possible ODE sizes is available.
        cur_core = 0
        cur_last = sys_maxsize
        for record in data:
            if record.cores != cur_core:
                if cur_last != sys_maxsize and cur_last != 'inf':
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
        return bool(cur_last == sys_maxsize) or (cur_last == 'inf')


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
                     Column('createdBy', String, default=getuser()),
                     Column('createdIn', String, default=__version__),
                     Column('createdOn', DateTime, default=datetime.now),
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
               frequency: float, cores: int, mode: 'ModelToolType') -> 'pandas.DataFrame':
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
        mode : str
            Application mode used to obtain data.

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
                   ImplVariantRecord.cores.is_(cores), ImplVariantRecord.mode.is_(mode)). \
            order_by(ImplVariantRecord.impl)
        data = read_sql_query(query.statement, db_session.bind, index_col='impl')
        return data
