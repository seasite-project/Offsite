"""@package descriptions.machine.compiler
Definition of class Compiler.

@author: Johannes Seiferth
"""

from datetime import datetime
from getpass import getuser
from typing import Dict

import attr
from sqlalchemy import Column, DateTime, Integer, String, Table, UniqueConstraint
from sqlalchemy.exc import NoResultFound, MultipleResultsFound
from sqlalchemy.orm import Session

from offsite import __version__
from offsite.database import METADATA, insert
from offsite.util.process_utils import run_process


@attr.s
class Compiler:
    """Representation of a Compiler object.

    A Compiler object describes a compiler configuration.

    Attributes:
    -----------
    name: str
        Name of the compiler.
    version: str
        Used version of the compiler.
    flags: str
        Used compiler flags.
    db_id: int
        ID of associated Compiler database table record.
    """
    name = attr.ib(type=str)
    version = attr.ib(type=str)
    flags = attr.ib(type=str)
    db_id = attr.ib(type=int, init=False)

    # Database information.
    db_table = Table('compiler', METADATA,
                     Column('db_id', Integer, primary_key=True),
                     Column('name', String),
                     Column('version', String),
                     Column('flags', String),
                     Column('updatedIn', String, default=__version__),
                     Column('updatedOn', DateTime, default=datetime.now, onupdate=datetime.now),
                     Column('updatedBy', String, default=getuser(), onupdate=getuser()),
                     UniqueConstraint('name', 'version', 'flags'),
                     sqlite_autoincrement=True)

    @classmethod
    def from_yaml(cls, yaml: Dict[str, str], used_compiler: str) -> 'Compiler':
        """Construct Compiler object from YAML definition.

        Parameters:
        -----------
        yaml: dict
            YAML object describing this object.
        used_compiler: str
            Name of used compiler.
        Returns:
        --------
        Compiler
            Created Compiler object.
        """
        try:
            item = next(filter(lambda x: x[0] == used_compiler, yaml.items()))
            # Attribute compiler_name.
            name = item[0]
            # Attribute compiler_version.
            version = Compiler.determine_compiler_version(name)
            # Attribute compiler_flags.
            flags = item[1]
        except StopIteration:
            raise RuntimeError('Compiler {} not found in YAML description!'.format(used_compiler))
        return cls(name, version, flags)

    @classmethod
    def from_database(cls, db_session: Session, compiler_id: int) -> 'Compiler':
        """Construct Compiler object from database record.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        compiler_id: int
            Database ID of the requested Compiler object.

        Returns:
        --------
        Compiler
            Created Compiler object.
        """
        try:
            compiler: Compiler = db_session.query(Compiler).filter(Compiler.db_id.is_(compiler_id)).one()
        except NoResultFound:
            raise RuntimeError('Unable to load Compiler object from database!')
        except MultipleResultsFound:
            raise RuntimeError('Unable to load Compiler object from database!')
        return compiler

    def to_database(self, db_session: Session) -> 'Compiler':
        """Push this compiler object to the database.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        Compiler
            Instance of this object connected to database session.
        """
        # Check if database already contains the compiler object.
        compiler: Compiler = db_session.query(Compiler).filter(Compiler.name.like(self.name),
                                                               Compiler.version.like(self.version),
                                                               Compiler.flags.like(self.flags)).one_or_none()
        if compiler:
            return compiler
        # Add new object to database.
        insert(db_session, self)
        return self

    @staticmethod
    def determine_compiler_version(compiler: str) -> str:
        """Determine the version of the compiler used and return the version string.

        Parameters:
        -----------
        compiler: str
            Used compiler command.

        Returns:
        --------
        str
            Version string of the used compiler.
        """
        version = run_process([compiler, '-dumpversion'])
        # Strip trailing EOL characters.
        return version.rstrip()
