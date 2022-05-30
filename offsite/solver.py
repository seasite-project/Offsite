"""@package solver
Definition of class Solver.

@author: Johannes Seiferth
"""

from datetime import datetime
from getpass import getuser
from typing import List

import attr
from sqlalchemy import Column, DateTime, Enum, Integer, String, Table, UniqueConstraint
from sqlalchemy.orm import Session

from offsite import __version__
from offsite.config import SolverType, SolverSpecificTableType
from offsite.database import METADATA, insert


@attr.s
class Solver:
    """Representation of a Solver table database record.

    Attributes:
    -----------
    type: SolverType
        Type of this Solver class.
    db_id: int
        ID of associated solver database table record.
    """
    name = attr.ib(type=str)
    type = attr.ib(type=SolverType, default=SolverType.GENERIC)
    specific_tables = attr.ib(type=List[SolverSpecificTableType], default=[])
    specific_tables_serial = attr.ib(type=str, init=False)
    db_id = attr.ib(type=int, init=False)

    # Database information.
    db_table = Table('solver', METADATA,
                     Column('db_id', Integer, primary_key=True),
                     Column('name', String),
                     Column('type', Enum(SolverType)),
                     Column('specific_tables_serial', String),
                     Column('updatedIn', String, default=__version__),
                     Column('updatedOn', DateTime, default=datetime.now, onupdate=datetime.now),
                     Column('updatedBy', String, default=getuser(), onupdate=getuser()),
                     UniqueConstraint('name'),
                     sqlite_autoincrement=True)

    def to_database(self, db_session: Session) -> 'Solver':
        """Push this Solver record object to the database.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        Solver
            Instance of this object connected to database session.
        """
        # Attribute kernels_serial.
        self.specific_tables_serial = ','.join((x.value for x in self.specific_tables))
        # Check if database already contains this ImplVariant object.
        solver: Solver = db_session.query(Solver).filter(
            Solver.type.is_(self.type), Solver.specific_tables_serial.is_(self.specific_tables_serial)).one_or_none()
        if solver:
            # Supplement attributes not saved in database.
            solver.specific_tables = self.specific_tables
            # Update serialized members.
            solver.specific_tables_serial = ','.join((x.value for x in self.specific_tables))
            return solver
        # Add new object to database.
        # Attribute kernels_serial.
        self.specific_tables_serial = ','.join((x.value for x in self.specific_tables))
        # Insert Solver object.
        insert(db_session, self)
        return self

    @classmethod
    def make_solver(cls, name_str: str, type_str: str):
        if type_str == SolverType.ODE.value:
            return cls(name=name_str, type=SolverType.ODE,
                       specific_tables=[SolverSpecificTableType.IVP, SolverSpecificTableType.ODE_METHOD])
        else:
            return cls(name=name_str, type=SolverType.GENERIC)
