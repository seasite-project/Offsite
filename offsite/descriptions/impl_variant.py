"""@package impl_variant
Definition of class ImplVariant.
"""

from datetime import datetime
from getpass import getuser
from typing import List

import attr
from sqlalchemy import Column, DateTime, ForeignKey, Integer, String, Table, UniqueConstraint
from sqlalchemy.orm import Session

from offsite import __version__
from offsite.db import METADATA
from offsite.db.db import insert
from offsite.descriptions.parser_utils import serialize_obj


@attr.s
class ImplVariant:
    """Representation of an ImplVariant table database record.

    Attributes:
    -----------
    skeleton : int
        Used impl skeleton.
    kernels : list of int
        Used machine.
    db_id : int
        ID of associated impl variant database table record.
    """
    skeleton = attr.ib(type=int)
    kernels = attr.ib(type=List[int])
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

    def to_database(self, db_session: Session):
        """Push this ECM record object to the database.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        ImplVariant
            Instance of this object connected to database session.
        """
        # Attribute kernels_serial.
        self.kernels_serial = ','.join(map(str, self.kernels))
        # Check if database already contains this ImplVariant object.
        variant = db_session.query(ImplVariant).filter(
            ImplVariant.skeleton.is_(self.skeleton), ImplVariant.kernels_serial.is_(self.kernels_serial)).one_or_none()
        if variant:
            # Supplement attributes not saved in database.
            variant.kernels = self.kernels
            # Update serialized members.
            variant.kernels_serial = serialize_obj(variant.kernels)
            return variant
        # Add new object to database.
        # Attribute kernels_serial.
        self.kernels_serial = ','.join(map(str, self.kernels))
        # Insert ImplVariantRecord object.
        insert(db_session, self)
        return self

    @staticmethod
    def select_all(db_session: Session) -> List['ImplVariant']:
        """Retrieve all ImplVariant table data record(s) from the database.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        list of ImplVariant
            Retrieved list of data records.
        """
        data = db_session.query(ImplVariant).all()
        return data

    @staticmethod
    def select(db_session: Session, variant_ids: List[int]) -> List['ImplVariant']:
        """
        Retrieve the ImplVariant table data record(s) from the database that match the given implementation variant IDs.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        variant_ids : list of int
            IDs of the impl variants requested.

        Returns:
        --------
        list of ImplVariant
            Retrieved list of data records.
        """
        variants = db_session.query(ImplVariant).filter(ImplVariant.db_id.in_(variant_ids)).all()
        if len(variants) != len(variant_ids):
            raise RuntimeError('Unable to select all requested variants!')
        # Deserialize attributes...
        for variant in variants:
            # Attribute kernels.
            variant.kernels = variant.kernels_serial.split(',')
        return variants
