"""@package database.db
Definition of database functions.

@author: Johannes Seiferth
"""

from typing import Any, List

from sqlalchemy.orm import Session


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
