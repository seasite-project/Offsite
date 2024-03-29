"""@package descriptions.ode_method
Definition of class ODEMethod.

@author: Johannes Seiferth
"""

from datetime import datetime
from getpass import getuser
from typing import List, Tuple, Union

import attr
from pathlib2 import Path
from sqlalchemy import Column, DateTime, Integer, String, Table
from sqlalchemy.exc import NoResultFound, MultipleResultsFound
from sqlalchemy.orm import Session

from offsite import __version__
from offsite.database import METADATA, insert
from offsite.descriptions.parser import load_yaml
from offsite.util.file_extensions import __ode_method_ext__


@attr.s
class ODEMethod:
    """Representation of a ODE method.

    A ODEMethod object describes the characteristics of a ODE method.

    Attributes:
    -----------
    name: str
        Name of this object.
    stages: int
        Number of stages of the ODE method.
    order_: int
        Order of the ODE method.
    correctorSteps: int
        Number of corrector steps of the ODE method.
    coefficientsA: list of list of str
        Coefficient matrix A of the ODE method.
    coefficientsB: list of str
        Coefficient vector b of the ODE method.
    coefficientsC: list of str
        Coefficient vector c of the ODE method.
    db_id: int
        ID of associated ODEMethod database table record.
    """
    name = attr.ib(type=str)
    stages = attr.ib(type=int)
    order_ = attr.ib(type=int)
    correctorSteps = attr.ib(type=int)
    coefficientsA = attr.ib(type=List[List[str]])
    coefficientsB = attr.ib(type=List[str])
    coefficientsC = attr.ib(type=List[str])
    coefficientsA_serial = attr.ib(type=str)
    coefficientsB_serial = attr.ib(type=str)
    coefficientsC_serial = attr.ib(type=str)
    db_id = attr.ib(type=int, init=False)

    # Database information.
    db_table = Table('ode_method', METADATA,
                     Column('db_id', Integer, primary_key=True),
                     Column('name', String, unique=True),
                     Column('stages', Integer),
                     Column('order_', Integer),
                     Column('correctorSteps', Integer),
                     Column('coefficientsA_serial', String),
                     Column('coefficientsB_serial', String),
                     Column('coefficientsC_serial', String),
                     Column('updatedIn', String, default=__version__),
                     Column('updatedOn', DateTime, default=datetime.now, onupdate=datetime.now),
                     Column('updatedBy', String, default=getuser(), onupdate=getuser()),
                     sqlite_autoincrement=True)

    @classmethod
    def from_yaml(cls, yaml_path: Path) -> 'ODEMethod':
        """Construct ODEMethod object from YAML definition.

        Parameters:
        -----------
        yaml_path: Path
            Relative path to this object's YAML file.

        Returns:
        --------
        ODEMethod
            Created ODEMethod object.
        """
        # Load YAML data.
        path = Path(yaml_path)
        yaml = load_yaml(path)
        # Attribute name.
        name = path.stem
        # Attribute stages.
        stages_ = yaml['stages']
        # Attribute order_.
        order_ = yaml['order']
        # Attribute corrector_steps.
        if 'corrector_steps' in yaml:
            # for implicit ODE methods.
            corrector_steps_ = yaml['corrector_steps']
        else:
            # for explicit ODE methods.
            corrector_steps_ = 0
        # Attribute coeff_a.
        coeff_a = yaml['A']
        coeff_a_serial = '; '.join([', '.join(row) for row in yaml['A']])
        # Attribute coeff_b.
        coeff_b = yaml['b']
        coeff_b_serial = ','.join(yaml['b'])
        # Attribute coeff_c.
        coeff_c = yaml['c']
        coeff_c_serial = ','.join(yaml['c'])
        # Create object.
        return cls(name, stages_, order_, corrector_steps_, coeff_a, coeff_b, coeff_c, coeff_a_serial, coeff_b_serial,
                   coeff_c_serial)

    @classmethod
    def from_database(cls, db_session: Session, ode_method_name: str) -> 'ODEMethod':
        """Construct ODEMethod object from database record.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        ode_method_name: str
            Name of the ODE method, which is used as primary key in the database.

        Returns:
        --------
        ODEMethod
            Created ODEMethod object.
        """
        try:
            method: ODEMethod = db_session.query(ODEMethod).filter(ODEMethod.name.like(ode_method_name)).one()
        except NoResultFound:
            raise RuntimeError('Unable to load ODEMethod object from database!')
        except MultipleResultsFound:
            raise RuntimeError('Unable to load ODEMethod object from database!')
        # Attribute coefficientsA.
        method.coefficientsA = list()
        for row in method.coefficientsA_serial.split(';'):
            method.coefficientsA.append(row.split(','))
        # Attribute coefficientsB.
        method.coefficientsB = method.coefficientsB_serial.split(',')
        # Attribute coefficientsC.
        method.coefficientsC = method.coefficientsC_serial.split(',')
        return method

    def to_database(self, db_session: Session) -> 'ODEMethod':
        """Push this ODEMethod object to the database.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        ODEMethod
            Instance of this object connected to database session.
        """
        # Check if database already contains the compiler object.
        method: ODEMethod = db_session.query(ODEMethod).filter(ODEMethod.name.like(self.name)).one_or_none()
        if method:
            # Supplement attributes not saved in database.
            method.coefficientsA = self.coefficientsA
            method.coefficientsB = self.coefficientsB
            method.coefficientsC = self.coefficientsC
            return method
        # Add new object to database.
        insert(db_session, self)
        return self

    @staticmethod
    def database_id(db_session: Session, yaml_path: Path) -> int:
        """Return database ID of an ODE method object given by its YAML description.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        yaml_path: Path
            Relative path to this object's YAML file.

        Returns:
        --------
        int
            Database ID of ODE method or -1 if not in database.
        """
        # Parse ODE method yaml description.
        method = ODEMethod.from_yaml(yaml_path)
        # Query ODE method in database.
        method = db_session.query(ODEMethod).filter(ODEMethod.name.like(method.name)).one_or_none()
        # Return ID of ODE method object.
        return method.db_id if method else -1

    @staticmethod
    def select_all(db_session: Session) -> List['ODEMethod']:
        """Retrieve all ODE method table data records  from the database.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        list of ODEMethod
            List of retrieved data records.
        """
        method_names: List[str] = [x[0] for x in db_session.query(ODEMethod.name).all()]
        return [ODEMethod.from_database(db_session, mname) for mname in method_names]

    @staticmethod
    def fetch_ode_method_name(db_session: Session, db_id: int):
        try:
            name: str = db_session.query(ODEMethod.name).filter(ODEMethod.db_id.like(db_id)).one()[0]
        except NoResultFound:
            raise RuntimeError('Unable to load ODE method object from database!')
        except MultipleResultsFound:
            raise RuntimeError('Unable to load unique ODE method object from database!')
        return name


def parse_method(path: Path) -> ODEMethod:
    """Parse ODE method description YAML file and return ODEMethod object.

    Parameters:
    -----------
    path: Path
        Relative path to file that contains the ODEMethod objects.

    Returns:
    --------
    ODEMethod
        ODEMethod object created.
    """
    if not path.suffix == __ode_method_ext__:
        raise RuntimeError('Invalid extension for ODE method: {}!'.format(path.suffix))
    return ODEMethod.from_yaml(path)


def parse_methods(path: Path) -> List[ODEMethod]:
    """Parse folder or file path for ODEMethod objects.

    Parameters:
    -----------
    path: Path
        Relative path to folder that contains the ODEMethod objects or a particular ODE method object.

    Returns:
    --------
    list of ODEMethod
        ODEMethod objects to be used.
    """
    if path.is_file():
        return [parse_method(path)]
    return [ODEMethod.from_yaml(mf) for mf in (f for f in path.iterdir() if f.suffix == __ode_method_ext__)]


def corrector_steps(method: 'ODEMethod') -> Tuple[str, Union[str, int]]:
    """Return corrector steps variable.

    Parameters:
    -----------
    method : ODEMethod
        Used ODE method.

    Returns:
    --------
    tuple(str, str)
        Corrector step variable.
    """
    return 'm', method.correctorSteps


def stages(method: ODEMethod) -> Tuple[str, Union[str, int]]:
    """Return stages variable.

    Parameters:
    -----------
    method : ODEMethod
        Used ODE method.

    Returns:
    --------
    tuple(str, str)
        Stage variable.
    """
    return 's', method.stages
