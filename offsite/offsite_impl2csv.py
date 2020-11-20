"""@package offsite_plot
Main script of the offsite_impl2csv application.
"""

from argparse import ArgumentParser, Namespace
from datetime import datetime

from matplotlib import pyplot
from pandas import read_sql_query, DataFrame

from offsite import __version__
from offsite.db.db import open_db, close
from offsite.db.db_mapping import mapping
from offsite.descriptions.ivp import IVP
from offsite.descriptions.ode_method import ODEMethod
from offsite.evaluation.math_utils import eval_math_expr, ivp_grid_size, ivp_system_size, solve_equation


def parse_program_args_app_impl2csv() -> Namespace:
    """Parse the available program arguments of the offsite_impl_to_csv application.

    Parameters:
    -----------
    -

    Returns:
    --------
    argparse.Namespace
        Parsed program arguments.
    """
    # Create argument parser object.
    parser = ArgumentParser(description='Plot data retrieved from tuning database.')
    # Available general options.
    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(__version__),
                        help='Print program version and exit.')
    parser.add_argument('--verbose', action='store_true', default=False, help='Print further information on this run.')
    parser.add_argument('--db', action='store', default='tune.db', help='Path to database. Default value: tune.db.')
    parser.add_argument('--machine', action='store', required=True, type=int, help='Database ID of used machine.')
    parser.add_argument('--compiler', action='store', required=True, type=int, help='Database ID of used compiler.')
    parser.add_argument('--cores', action='store', required=True, type=int,
                        help='Plot data for this number of CPU cores.')
    parser.add_argument('--frequency', action='store', required=True, type=float,
                        help='Plot data for this CPU frequency.')
    parser.add_argument('--method', action='store', type=str, help='Name of the used ODE method.')
    parser.add_argument('--ivp', action='store', type=str, default=None, help='Name of the used IVP.')
    parser.add_argument('--impl', action='store', type=str, help='Database ID of impl variant.')
    parser.add_argument('--N', action='store', type=str,
                        help='Range of plotted system sizes N. Used syntax: \'[first_N]:[last_N]:[incr_N]\'.')
    # Specific options to tweak appearance of plots.
    parser.add_argument('-tperComponent', action='store_true', default=False, help='')
    parser.add_argument('--no-plot', action='store_true', default=False, help='')

    # Parse program arguments.
    return parser.parse_args()


def plot(args):
    """Run the offsite_plot application.

    Parameters:
    -----------
    args: argparse.Namespace
        Program arguments.

    Returns:
    --------
    -
    """
    # Open database connection.
    session = open_db(args.db)
    # Plot data.
    if not args.method:
        raise RuntimeError('Misses required argument \'--method\'!')
    if not args.ivp:
        raise RuntimeError('Misses required argument \'--ivp\'!')
    if not args.impl:
        raise RuntimeError('Misses required argument \'--impl\'!')
    if not args.N:
        raise RuntimeError('Misses required argument \'--N\'!')
    plot_impl_variant_prediction(args, session)
    # Close database session.
    close(session)


def construct_impl_variant_prediction_query(args, ivp: int, method: int) -> str:
    query = "SELECT impl, first, last, prediction FROM impl_variant_prediction WHERE"
    if args.impl == 'all':
        pass
    else:
        impls = [x.strip() for x in args.impl.split(',')]
        if len(impls) == 1:
            query += " impl='{}' AND".format(str(args.impl))
        else:
            query += " (impl='" + "' OR impl='".join(impls) + "') AND"
    query += " machine='{}'".format(str(args.machine))
    query += " AND compiler='{}'".format(str(args.compiler))
    query += " AND method='{}'".format(str(method))
    query += " AND ivp='{}'".format(str(ivp))
    query += " AND cores='{}'".format(str(args.cores))
    query += " AND frequency='{}'".format(str(args.frequency))
    query += " ORDER BY impl, first, last"
    return query


def find_prediction(predictions: DataFrame, n: int) -> DataFrame:
    # Find prediction data for the the given ODE system size (n).
    if n == 0:
        return 0.0
    try:
        row = predictions.loc[(predictions['first'] <= n) & (predictions['last'] >= n)]
    except IndexError:
        raise RuntimeError('Failed to find prediction data for n={}:\n{}'.format(n, predictions))
    return row


def range_of_N(N_str: str) -> range:
    # Parse N range argument given by syntax 'first:last:increment'.
    split = N_str.split(':')
    if len(split) != 3:
        raise RuntimeError('Failed to parse range of Ns \'{}\'! '.format(N_str) +
                           'Supported syntax \'[first]:[last]:[increment]\'')
    first_N = int(split[0])
    assert first_N >= 0
    last_N = int(split[1])
    assert last_N >= first_N
    incr_N = int(split[2])
    assert incr_N > 0
    return range(first_N, last_N + 1, incr_N)


def plot_impl_variant_prediction(args, db_session):
    # Get data from database ...
    # ... IVP.
    ivp = IVP.from_database(db_session, args.ivp)
    # ... ODE method.
    method = ODEMethod.from_database(db_session, args.method)
    # ... implementation predictions.
    query = construct_impl_variant_prediction_query(args, ivp.db_id, method.db_id)
    data = read_sql_query(query, db_session.bind, index_col='impl')
    if data.empty:
        raise RuntimeWarning('Failed to find fitting data! Check passed database IDs!')

    # Split implementation prediction data by implementation variant ID.
    split_data = {idx: data.loc[idx] for idx in data.index.unique().values}

    plot_legend = list()
    for impl_id, impl_data in split_data.items():
        df = DataFrame(columns=('N', 'prediction'))

        idx = 0
        for N in range_of_N(args.N):
            # Determine corresponding ODE system size 'n' from IVP's grid size expression.
            n = int(solve_equation(eval_math_expr('g', [ivp_grid_size(N)], cast_to=str), ivp.gridSize, 'n')[0])
            # Find fitting prediction data.
            row = find_prediction(impl_data, n)
            # Write prediction data to data frame.
            if args.tperComponent:
                prediction: float = eval_math_expr(
                    str(row.prediction.iat[0]), [ivp_system_size(1), ('x', n)], cast_to=float)
            else:
                prediction: float = eval_math_expr(
                    str(row.prediction.iat[0]), [ivp_system_size(n), ('x', n)], cast_to=float)
            df.loc[idx] = [N] + [prediction]

            idx += 1
        # Plot implementation variant.
        pyplot.plot(df['N'], df['prediction'])
        # Save implementation variant name for legend.
        plot_legend.append('Variant {}'.format(impl_id))
        df.to_csv('impl_{}.csv'.format(impl_id))
    if not args.no_plot:
        # Plot data.
        pyplot.title('Implementation variant predictions')
        pyplot.legend(plot_legend)
        pyplot.xlabel('N')
        if args.tperComponent:
            pyplot.ylabel('time per component [s]')
        else:
            pyplot.ylabel('time [s]')
        # Save plot.
        pyplot.savefig('impl_variant_prediction_{}.png'.format(datetime.now().strftime('%Y-%m-%d_%H:%M:%S')))
        # Show plot.
        pyplot.show()


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
    args: Namespace = parse_program_args_app_impl2csv()
    # Run plot app.
    plot(args)
