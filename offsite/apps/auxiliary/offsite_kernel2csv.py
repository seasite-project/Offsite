"""@package apps.auxiliary.offsite_kernel2csv
Main script of the offsite_kernel2csv application.

@author: Johannes Seiferth
"""

from argparse import ArgumentParser, Namespace

from pandas import read_sql_query, DataFrame
from sqlalchemy.orm import Session

from offsite import __version__
from offsite.database import close, open_db
from offsite.database.db_mapping import mapping
from offsite.descriptions.impl.kernel_template import Kernel
from offsite.descriptions.ode import IVP, ODEMethod, ivp_system_size, ivp_grid_size
from offsite.util.math_utils import eval_math_expr, solve_equation


def parse_program_args_app_kernel2csv() -> Namespace:
    """Parse the available program arguments of the offsite_kernel2csv application.

    Parameters:
    -----------
    -

    Returns:
    --------
    argparse.Namespace
        Parsed program arguments.
    """
    # Create argument parser object.
    parser = ArgumentParser(description='Write kernel prediction data from the database to CSV.')
    # Available general options.
    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(__version__),
                        help='Print program version and exit.')
    parser.add_argument('--db', action='store', required=True, help='Path to used database.')
    parser.add_argument('--machine', action='store', required=True, type=int, help='Database ID of used machine.')
    parser.add_argument('--compiler', action='store', required=True, type=int, help='Database ID of used compiler.')
    parser.add_argument('--cores', action='store', required=True, type=int,
                        help='Plot data for this number of CPU cores.')
    parser.add_argument('--frequency', action='store', required=True, type=float,
                        help='Plot data for this CPU frequency.')
    parser.add_argument('--method', action='store', required=True, type=str, help='Name of the used ODE method.')
    parser.add_argument('--ivp', action='store', required=True, type=str, default=None, help='Name of the used IVP.')
    parser.add_argument('--kernel', action='store', required=True, type=str, help='Database ID of kernel')
    parser.add_argument('--N', action='store', required=True, type=str,
                        help='Range of considered system sizes N. Used syntax: \'[first_N]:[last_N]:[incr_N]\'.')
    # Parse program arguments.
    return parser.parse_args()


def kernel2csv(args: Namespace):
    """Run the offsite_kernel2csv application.

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
    # Write implementation prediction data to csv.
    write_kernel_prediction_data(args, session)
    # Close database session.
    close(session)


def construct_kernel_prediction_query(args: Namespace, ivp: int, method: int) -> str:
    query = "SELECT kernel, first, last, prediction FROM kernel_prediction WHERE"
    if args.kernel == 'all':
        pass
    else:
        kernels = [x.strip() for x in args.kernel.split(',')]
        if len(kernels) == 1:
            query += " kernel='{}' AND".format(str(args.kernel))
        else:
            query += " (kernel='" + "' OR kernel'".join(kernels) + "') AND"
    query += " machine='{}'".format(str(args.machine))
    query += " AND compiler='{}'".format(str(args.compiler))
    query += " AND method='{}'".format(str(method))
    query += " AND (ivp='{}' OR ivp='-1')".format(str(ivp))
    query += " AND cores='{}'".format(str(args.cores))
    query += " AND frequency='{}'".format(str(args.frequency))
    query += " ORDER BY kernel, first, last"
    return query


def find_prediction(predictions: DataFrame, n: int) -> str:
    # Search in the passed data frame for the particular prediction data that fits the given ODE system size (n).
    if n == 0:
        return str(0.0)
    try:
        row: DataFrame = predictions.loc[(predictions['first'] <= n) & (predictions['last'] >= n)]
    except IndexError:
        raise RuntimeError('Failed to find prediction data for n={}:\n{}'.format(n, predictions))
    return row.prediction.iat[0]


def construct_range_of_N(N_expr: str) -> range:
    # Parse N range argument given by syntax 'first:last:increment' to determine the range of N values considered.
    split = N_expr.split(':')
    if len(split) != 3:
        raise RuntimeError('Failed to parse range of Ns \'{}\'! '.format(N_expr) +
                           'Supported syntax \'[first]:[last]:[increment]\'')
    first = int(split[0])
    last = int(split[1])
    incr = int(split[2])
    # Validate range.
    assert first >= 0
    assert last >= first
    assert incr > 0
    # last + 1 since range's argument 'stop' is exclusive.
    return range(first, last + 1, incr)


def write_kernel_prediction_data(args: Namespace, db_session: Session):
    # Fetch all objects required to determine the correct kernel prediction data, from the database.
    # ... used IVP.
    ivp = IVP.from_database(db_session, args.ivp)
    # ... used ODE method.
    method = ODEMethod.from_database(db_session, args.method)

    # Query fitting kernel prediction data ...
    query = construct_kernel_prediction_query(args, ivp.db_id, method.db_id)
    data = read_sql_query(query, db_session.bind, index_col='kernel')
    if data.empty:
        raise RuntimeWarning('Failed to find fitting data! Check passed database IDs!')
    # ... split these data by kernel ID.
    data = {idx: data.loc[idx] for idx in data.index.unique().values}

    # Evaluate prediction data for all 'N' values given and write results to CSV.
    for kernel_id, kernel_data in data.items():
        df = DataFrame(columns=('N', 't_perComponent', 't_timestep', 'MLUPs'))
        cur_row_idx = 0
        for N in construct_range_of_N(args.N):
            # Determine ODE system size 'n' from 'N' by solving the grid size expression of the IVP.
            # Note: Internally 'g' is used instead of 'N' since 'N' is already reserved in sympy.
            # E.g. g = sqrt(n) for Heat2D
            lhs = eval_math_expr('g', [ivp_grid_size(N)], cast_to=str)
            n: int = int(solve_equation(lhs, ivp.gridSize, 'n')[0])
            # Find the fitting prediction expression string.
            pred_expr: str = find_prediction(kernel_data, n)
            # Evaluate prediction data ...
            # ... predicted time per component.
            t_component: float = eval_math_expr(pred_expr, [ivp_system_size(1), ('x', n)], cast_to=float)
            # ... predicted time per timestep.
            t_timestep: float = eval_math_expr(pred_expr, [ivp_system_size(n), ('x', n)], cast_to=float)
            # ... predicted obtained MLUPs.
            mlups = float(n * 1e-6) / float(t_timestep)
            # Write evaluated prediction data to dataframe.
            df.loc[cur_row_idx] = [N] + [t_component] + [t_timestep] + [mlups]
            cur_row_idx += 1
        # Write data frame to CSV.
        df.to_csv('kernel_{}_{}.csv'.format(kernel_id, Kernel.fetch_kernel_name(db_session, int(kernel_id))))


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
    args: Namespace = parse_program_args_app_kernel2csv()
    # Run plot app.
    kernel2csv(args)
