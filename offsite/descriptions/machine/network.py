"""@package descriptions.machine.network
Definition of class NetworkConfiguration.

@author: Hana Shatri, Johannes Seiferth
"""

from datetime import datetime
from getpass import getuser

import attr
from sqlalchemy import Column, DateTime, Integer, String, Table

from offsite import __version__
from offsite.database import METADATA


@attr.s
class NetworkConfig:
    """
        Representation of a network configuration table database record.
    """
    # clock = attr.ib(type=float, default=None) # TODO shouldn't be needed?!
    network_interface = attr.ib(type=str)
    broadcast = attr.ib(type=int)
    multicast = attr.ib(type=int)
    mtu = attr.ib(type=int)
    qdisc = attr.ib(type=str)
    qlen = attr.ib(type=int)
    conn_speed = attr.ib(type=str)
    db_id = attr.ib(type=int, init=False)

    # Database information.
    db_table = Table('network_config', METADATA,
                     Column('db_id', Integer, primary_key=True),
                     # Column('clock', Float),
                     Column('network_interface', String),
                     Column('broadcast', Integer),
                     Column('multicast', Integer),
                     Column('mtu', Integer),
                     Column('qdisc', String),
                     Column('qlen', Integer),
                     Column('conn_speed', String),
                     Column('updatedIn', String, default=__version__),
                     Column('updatedOn', DateTime, default=datetime.now, onupdate=datetime.now),
                     Column('updatedBy', String, default=getuser(), onupdate=getuser()),
                     sqlite_autoincrement=True)
