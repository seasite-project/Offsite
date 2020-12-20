"""@package database
Definition of database functions.
"""

from typing import Any, List

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, Session

from offsite.db import METADATA


def open_db(db_name: str) -> Session:
    """Open connection to database.

    Parameters:
    -----------
    db_name: str
        Name of the used database.

    Returns:
    --------
    sqlalchemy.orm.session.Session
        Opened database connection handle.
    """
    try:
        # Create engine.
        engine = create_engine('sqlite:///{}'.format(db_name), echo=False, connect_args={'timeout': 15})
        # Create tables.
        METADATA.create_all(engine)
        # Open session.
        maker = sessionmaker(bind=engine, expire_on_commit=False)
        session = maker()
    except:
        raise RuntimeError('Failed to connect to database {}!'.format(db_name))
    return session


def close(session: Session):
    """Close connection to database.

    Parameters:
    -----------
    session: sqlalchemy.orm.session.Session
        Used database session.

    Returns:
    --------
    -
    """
    session.close()


def commit(session: Session):
    """Commit changes to the database.

    Parameters:
    -----------
    session: sqlalchemy.orm.session.Session
        Used database session.

    Returns:
    --------
    -
    """
    try:
        session.commit()
    except:
        raise RuntimeError('Failed to commit database changes!')


def rollback(session: Session):
    """Rollback unsaved changes to database.

    Parameters:
    -----------
    session: sqlalchemy.orm.session.Session
        Used database session.

    Returns:
    --------
    -
    """
    try:
        session.rollback()
    except:
        raise RuntimeError('Failed to rollback database changes!')


def insert(session: Session, record: Any):
    """Insert new data record into database.

    Parameters:
    -----------
    session: sqlalchemy.orm.session.Session
        Used database session.
    record: Object
        Record to be inserted.

    Returns:
    --------
    -
    """
    try:
        session.add(record)
    except:
        raise RuntimeError('Failed to insert record into database!')


def bulk_insert(session: Session, records: List[Any]):
    """Insert new data record into database.

    Parameters:
    -----------
    session: sqlalchemy.orm.session.Session
        Used database session.
    record: Object
        Record to be inserted.

    Returns:
    --------
    -
    """
    try:
        session.bulk_save_objects(records)
    except:
        raise RuntimeError('Failed to insert record into database!')
