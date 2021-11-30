"""@package database
Database functions and mapping.
"""

from sqlalchemy import MetaData, create_engine
from sqlalchemy.orm import sessionmaker, Session

from offsite.database.db import bulk_insert, close, commit, insert, rollback

# Setup sqlalchemy meta data class.
METADATA = MetaData()


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
