"""@package offsite_db2name
Main script of the offsite_db2name application.
"""

from argparse import ArgumentParser, Namespace
from enum import Enum
from typing import List

from sqlalchemy.orm import Session
from sqlalchemy.orm.exc import NoResultFound, MultipleResultsFound

from offsite import __version__
from offsite.codegen.codegen_util import create_variant_name
from offsite.db.db import open_db, close
from offsite.db.db_mapping import mapping
from offsite.descriptions.impl_skeleton import ImplSkeleton
from offsite.descriptions.impl_variant import ImplVariant
from offsite.descriptions.ivp import IVP
from offsite.descriptions.kernel_template import Kernel
from offsite.descriptions.ode_method import ODEMethod


class DBRecordType(Enum):
    """Defines what type of database record is queried.

    - KERNEL
        Query Kernel record information.
    - IMPL
        Query ImplVariant record information.
    - IVP
        Query IVP record information.
    - METHOD
        Query ODEMethod record information.
    - SKELETON
        Query ImplSkeleton record information.

    """
    KERNEL = 'KERNEL'
    IMPL = 'IMPL'
    IVP = 'IVP'
    METHOD = 'METHOD'
    SKELETON = 'SKELETON'


def parse_program_args_app_db2name() -> Namespace:
    """Parse the available program arguments of the offsite_db2name application.

    Parameters:
    -----------
    -

    Returns:
    --------
    argparse.Namespace
        Parsed program arguments.
    """
    # Create argument parser object.
    parser = ArgumentParser(description='Print name identifier of a specific user-given database record.')
    # Available general options.
    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(__version__),
                        help='Print program version and exit.')
    parser.add_argument('--db', action='store', required=True, help='Path to database.')
    parser.add_argument('--record', action='store', required=True, type=DBRecordType,
                        help='Specify what type of record is queried. Possible values: IMPL, IVP, KERNEL, METHOD,'
                             'SKELETON.')
    parser.add_argument('--id', action='store', required=True, type=int,
                        help='Database ID of the database record requested.')
    # Parse program arguments.
    return parser.parse_args()


def db2name(args: Namespace):
    """Run the offsite_db2name application.

    Parameters:
    -----------
    args: argparse.Namespace
        Program arguments.

    Returns:
    --------
    -
    """
    # Open database connection.
    session: Session = open_db(args.db)
    # Print information on the given ...
    if args.record == DBRecordType.IMPL:
        name: str = fetch_impl_variant_name(session, args.id)
        print('Record \'{}\' belongs to ImplVariant \'{}\'.'.format(args.id, name))
    elif args.record == DBRecordType.IVP:
        name: str = fetch_ivp_name(session, args.id)
        print('Record \'{}\' belongs to IVP \'{}\'.'.format(args.id, name))
    elif args.record == DBRecordType.KERNEL:
        name: str = fetch_kernel_name(session, args.id)
        print('Record \'{}\' belongs to Kernel \'{}\'.'.format(args.id, name))
    elif args.record == DBRecordType.METHOD:
        name: str = fetch_ode_method_name(session, args.id)
        print('Record \'{}\' belongs to ODE method \'{}\'.'.format(args.id, name))
    elif args.record == DBRecordType.SKELETON:
        name: str = fetch_impl_skeleton_name(session, args.id)
        print('Record \'{}\' belongs to ImplSkeleton \'{}\'.'.format(args.id, name))
    else:
        raise RuntimeError('Unsupported record type {}!'.format(args.type))
    # Close database session.
    close(session)


def fetch_impl_variant_name(db_session: Session, db_id: int) -> str:
    # Get implementation variant record first.
    variant: ImplVariant = ImplVariant.select(db_session, [db_id])[0]
    # Next, get records of all kernels used in the variant as well as the name of the variant's skeleton.
    try:
        skeleton: str = db_session.query(ImplSkeleton.name).filter(ImplSkeleton.db_id.is_(variant.skeleton)).one()[0]
    except NoResultFound:
        raise RuntimeError('Unable to load IVP object from database!')
    except MultipleResultsFound:
        raise RuntimeError('Unable to load unique IVP object from database!')
    # ... kernel records.
    kernels: List[Kernel] = [Kernel.select(db_session, kid) for kid in variant.kernels]
    # Derive implementation name from given information.
    name: str = create_variant_name(kernels, skeleton)
    return name


def fetch_ivp_name(db_session: Session, db_id: int) -> str:
    try:
        name: str = db_session.query(IVP.name).filter(IVP.db_id.like(db_id)).one()[0]
    except NoResultFound:
        raise RuntimeError('Unable to load IVP object from database!')
    except MultipleResultsFound:
        raise RuntimeError('Unable to load unique IVP object from database!')
    return name


def fetch_kernel_name(db_session: Session, db_id: int) -> str:
    try:
        name: str = db_session.query(Kernel.name).filter(Kernel.db_id.like(db_id)).one()[0]
    except NoResultFound:
        raise RuntimeError('Unable to load kernel object from database!')
    except MultipleResultsFound:
        raise RuntimeError('Unable to load unique kernel object from database!')
    return name


def fetch_ode_method_name(db_session: Session, db_id: int):
    try:
        name: str = db_session.query(ODEMethod.name).filter(ODEMethod.db_id.like(db_id)).one()[0]
    except NoResultFound:
        raise RuntimeError('Unable to load ODE method object from database!')
    except MultipleResultsFound:
        raise RuntimeError('Unable to load unique ODE method object from database!')
    return name


def fetch_impl_skeleton_name(db_session: Session, db_id: int):
    try:
        name: str = db_session.query(ImplSkeleton.name).filter(ImplSkeleton.db_id.like(db_id)).one()[0]
    except NoResultFound:
        raise RuntimeError('Unable to load ImplSkeleton object from database!')
    except MultipleResultsFound:
        raise RuntimeError('Unable to load unique ImplSkeleton object from database!')
    return name


def run():
    """Run command line interface.

    Parameters:
    -----------
    -

    Returns:
    -
    """
    # Map database.
    mapping()
    # Create parser and parse arguments.
    args: Namespace = parse_program_args_app_db2name()
    # Run db2name app.
    db2name(args)
