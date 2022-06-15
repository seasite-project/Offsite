"""@package apps.auxiliary.offsite_db2name
Main script of the offsite_db2name application.

@author: Johannes Seiferth
"""

from enum import Enum

from argparse import ArgumentParser, Namespace
from sqlalchemy.orm import Session

from offsite import __version__
from offsite.database import close, open_db
from offsite.database.db_mapping import mapping
from offsite.descriptions.impl.impl_skeleton import ImplSkeleton
from offsite.descriptions.impl.kernel_template import Kernel
from offsite.descriptions.ode import IVP, ODEMethod
from offsite.train.impl_variant import ImplVariant


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
    parser.add_argument('--db', action='store', required=True, help='Path to used database.')
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
        name: str = ImplVariant.fetch_impl_variant_name(session, args.id)
        print('Record \'{}\' belongs to ImplVariant \'{}\'.'.format(args.id, name))
    elif args.record == DBRecordType.IVP:
        name: str = IVP.fetch_ivp_name(session, args.id)
        print('Record \'{}\' belongs to IVP \'{}\'.'.format(args.id, name))
    elif args.record == DBRecordType.KERNEL:
        name: str = Kernel.fetch_kernel_name(session, args.id)
        print('Record \'{}\' belongs to Kernel \'{}\'.'.format(args.id, name))
    elif args.record == DBRecordType.METHOD:
        name: str = ODEMethod.fetch_ode_method_name(session, args.id)
        print('Record \'{}\' belongs to ODE method \'{}\'.'.format(args.id, name))
    elif args.record == DBRecordType.SKELETON:
        name: str = ImplSkeleton.fetch_impl_skeleton_name(session, args.id)
        print('Record \'{}\' belongs to ImplSkeleton \'{}\'.'.format(args.id, name))
    else:
        raise RuntimeError('Unsupported record type {}!'.format(args.type))
    # Close database session.
    close(session)


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
    args: 'Namespace' = parse_program_args_app_db2name()
    # Run db2name app.
    db2name(args)
