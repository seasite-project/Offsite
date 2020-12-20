"""@package offsite_plot
Main script of the offsite_plot application.
"""

from argparse import ArgumentParser, Namespace
from datetime import datetime
from enum import Enum

from matplotlib import pyplot
from pandas import read_sql_query, DataFrame
from sqlalchemy.orm import Session

from offsite import __version__
from offsite.apps.auxiliary.offsite_db2name import fetch_impl_variant_name, fetch_kernel_name
from offsite.db.db import open_db, close
from offsite.db.db_mapping import mapping
from offsite.evaluation.math_utils import eval_math_expr, ivp_system_size


class PlotType(Enum):
    """Defines what type of data is plotted.

    - BENCHMARK
        Plot communication benchmark data.
    - IMPL
        Plot Implementation Variant prediction data.
    - KERNEL
        Plot Kernel prediction data.

    """
    BENCHMARK = 'BENCHMARK'
    KERNEL = 'KERNEL'
    IMPL = 'IMPL'


def parse_program_args_app_plot() -> Namespace:
    """Parse the available program arguments of the offsite_plot application.

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
    parser.add_argument('--db', action='store', required=True, help='Path to database. Default value: tune.db.')
    parser.add_argument('--plot', action='store', required=True, type=PlotType,
                        help='Specify what is plotted. Possible values: BENCHMARK, IMPL, KERNEL.')
    parser.add_argument('--machine', action='store', required=True, type=int, help='Database ID of used machine.')
    parser.add_argument('--compiler', action='store', required=True, type=int, help='Database ID of used compiler.')
    parser.add_argument('--cores', action='store', required=True, type=int,
                        help='Plot data for this number of CPU cores.')
    parser.add_argument('--frequency', action='store', required=True, type=float,
                        help='Plot data for this CPU frequency.')
    # Specific options required by plot type 'BENCHMARK'.
    parser.add_argument('--benchmark', action='store', type=str, default=None,
                        help='Plot data for this communication benchmark.')
    # Specific options required by plot types 'KERNEL' and 'IMPL'.
    parser.add_argument('--method', action='store', type=int, help='Database ID of used ODE method.')
    parser.add_argument('--ivp', action='store', type=int, default=None, help='Database ID of used IVP.')
    # Specific options required by plot type 'KERNEL'.
    parser.add_argument('--kernel', action='store', type=str, help='Database ID of used (pmodel) kernel.')
    # Specific options required by plot type 'IMPL'.
    parser.add_argument('--impl', action='store', type=str, help='Database ID of impl variant.')
    # Specific options to tweak appearance of plots 'IMPL' and 'KERNEL'.
    parser.add_argument('-tperComponent', action='store_true', default=False, help='')
    # Parse program arguments.
    return parser.parse_args()


def plot(args: Namespace):
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
    if args.plot == PlotType.BENCHMARK:
        if not args.benchmark:
            raise RuntimeError('Plot type {} misses required argument \'--benchmark\'!'.format(args.plot))
        plot_benchmark(args, session)
    else:
        if not args.method:
            raise RuntimeError('Plot type {} misses required argument \'--method\'!'.format(args.plot))
        if args.plot == PlotType.KERNEL:
            if not args.kernel:
                raise RuntimeError('Plot type {} misses required argument \'--kernel\'!'.format(args.plot))
            plot_kernel_prediction(args, session)
        elif args.plot == PlotType.IMPL:
            if not args.ivp:
                raise RuntimeError('Plot type {} misses required argument \'--ivp\'!'.format(args.plot))
            if not args.impl:
                raise RuntimeError('Plot type {} misses required argument \'--impl\'!'.format(args.plot))
            plot_impl_variant_prediction(args, session)
        else:
            raise RuntimeError('Unsupported plot type {}!'.format(args.plot))
    # Close database session.
    close(session)


def plot_benchmark(args: Namespace, db_session: Session):
    # Construct SQL query.
    sql_query = "SELECT cores, data FROM benchmark_result WHERE"
    sql_query += " name='{}'".format(args.benchmark)
    sql_query += " AND machine='{}'".format(str(args.machine))
    sql_query += " AND compiler='{}'".format(str(args.compiler))
    sql_query += " AND frequency='{}'".format(str(args.frequency))
    sql_query += " ORDER BY cores"
    # Get data from database.
    data = read_sql_query(sql_query, db_session.bind)
    if data.empty:
        raise RuntimeWarning('Failed to find fitting data! Check passed database IDs!')
    # Plot data.
    data.plot(x='cores', y='data', title='Benchmark: ' + args.benchmark)
    # Configure plot.
    pyplot.xlabel('# cores')
    pyplot.ylabel('time [s]')
    pyplot.savefig('benchmark_{}_{}.png'.format(args.benchmark, datetime.now().strftime('%Y-%m-%d_%H:%M:%S')))
    # Show plot.
    pyplot.show()


def construct_kernel_prediction_query(args: Namespace) -> str:
    query = "SELECT kernel, first, last, prediction FROM kernel_prediction WHERE"
    if args.kernel == 'all':
        pass
    else:
        kernels = [x.strip() for x in args.kernel.split(',')]
        if len(kernels) == 1:
            query += " kernel='{}' AND".format(str(args.kernel))
        else:
            query += " (kernel='" + "' OR kernel='".join(kernels) + "') AND"
    query += " machine='{}'".format(str(args.machine))
    query += " AND compiler='{}'".format(str(args.compiler))
    query += " AND method='{}'".format(str(args.method))
    if args.ivp is not None:
        query += " AND (ivp='{}' OR ivp='-1')".format(str(args.ivp))
    else:
        query += " AND ivp='-1'"
    query += " AND cores='{}'".format(str(args.cores))
    query += " AND frequency='{}'".format(str(args.frequency))
    query += " ORDER BY kernel, first, last"
    return query


def plot_kernel_prediction(args: Namespace, db_session: Session):
    # Construct SQL query.
    query = construct_kernel_prediction_query(args)
    # Get data from database.
    data = read_sql_query(query, db_session.bind, index_col='kernel')
    if data.empty:
        raise RuntimeWarning('Failed to find fitting data! Check passed database IDs!')
    if args.verbose:
        print(data)
    # Split data by kernel ID
    split_data = {idx: data.loc[idx] for idx in data.index.unique().values}
    # Format data.
    legend = list()
    for kernel_id, kernel_data in split_data.items():
        # Fetch kernel name from database.
        kernel = fetch_kernel_name(db_session, int(kernel_id))
        #
        df = DataFrame(columns=('n', 'prediction'))
        #
        cur_row_idx = 0
        for rid, row in kernel_data.iterrows():
            first = row.loc['first']
            last = row.loc['last']
            prediction = row.loc['prediction']
            #
            if args.tperComponent:
                size = ivp_system_size(1)
                #
                df.loc[cur_row_idx] = [first] + [eval_math_expr(prediction, [size, ('x', first)], cast_to=float)]
                df.loc[cur_row_idx + 1] = [last] + [eval_math_expr(prediction, [size, ('x', last)], cast_to=float)]
            else:
                #
                pred_first: float = eval_math_expr(prediction, [ivp_system_size(first), ('x', first)], cast_to=float)
                df.loc[cur_row_idx] = [first] + [pred_first]
                #
                pred_last: float = eval_math_expr(prediction, [ivp_system_size(last), ('x', last)], cast_to=float)
                df.loc[cur_row_idx + 1] = [last] + [pred_last]
            cur_row_idx = cur_row_idx + 2
        # Plot implementation variant.
        pyplot.plot(df['n'], df['prediction'])
        # Save implementation variant name for legend.
        legend.append('{} (id={})'.format(kernel, kernel_id))
    # Configure plot.
    pyplot.title('Kernel predictions')
    pyplot.legend(legend)
    pyplot.xlabel('n')
    if args.tperComponent:
        pyplot.ylabel('time per component [s]')
    else:
        pyplot.ylabel('time [s]')
    # Save plot.
    pyplot.savefig('kernel_prediction_{}.png'.format(datetime.now().strftime('%Y-%m-%d_%H:%M:%S')))
    # Show plot.
    pyplot.show()


def construct_impl_variant_prediction_query(args):
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
    query += " AND method='{}'".format(str(args.method))
    query += " AND ivp='{}'".format(str(args.ivp))
    query += " AND cores='{}'".format(str(args.cores))
    query += " AND frequency='{}'".format(str(args.frequency))
    query += " ORDER BY impl, first, last"
    return query


def plot_impl_variant_prediction(args, db_session):
    # Construct SQL query.
    query = construct_impl_variant_prediction_query(args)
    # Get data from database.
    data = read_sql_query(query, db_session.bind, index_col='impl')
    if data.empty:
        raise RuntimeWarning('Failed to find fitting data! Check passed database IDs!')
    if args.verbose:
        print(data)
    # Split data by implementation variant ID
    split_data = {idx: data.loc[idx] for idx in data.index.unique().values}
    # Format data.
    legend = list()
    #
    for impl_id, impl_data in split_data.items():
        # Fetch impl variant name from database.
        impl = fetch_impl_variant_name(db_session, int(impl_id))
        #
        df = DataFrame(columns=('n', 'prediction'))
        cur_row_idx = 0
        for rid, row in impl_data.iterrows():
            first = row.loc['first']
            last = row.loc['last']
            prediction = row.loc['prediction']
            #
            if args.tperComponent:
                size = ivp_system_size(1)
                #
                df.loc[cur_row_idx] = [first] + [eval_math_expr(prediction, [size, ('x', first)], cast_to=float)]
                df.loc[cur_row_idx + 1] = [last] + [eval_math_expr(prediction, [size, ('x', last)], cast_to=float)]
            else:
                #
                first_pred: float = eval_math_expr(prediction, [ivp_system_size(first), ('x', first)], cast_to=float)
                df.loc[cur_row_idx] = [first] + [first_pred]
                #
                last_pred: float = eval_math_expr(prediction, [ivp_system_size(last), ('x', last)], cast_to=float)
                df.loc[cur_row_idx + 1] = [last] + [last_pred]
            cur_row_idx += 2
        # Plot implementation variant.
        pyplot.plot(df['n'], df['prediction'])
        # Save implementation variant name for legend.
        legend.append('{} (id={})'.format(impl, impl_id))
        df.to_csv('impl_{}.csv'.format(impl_id))
    # Configure plot.
    pyplot.title('Implementation variant predictions')
    pyplot.legend(legend)
    pyplot.xlabel('n')
    if args.tperComponent:
        pyplot.ylabel('time per component [s]')
    else:
        pyplot.ylabel('time [s]')
    # Save plot.
    pyplot.savefig('impl_variant_prediction_{}.png'.format(datetime.now().strftime('%Y-%m-%d_%H:%M:%S')))
    # Show plot.
    pyplot.show()


def evaluate_prediction(predictions: DataFrame, n: int) -> float:
    # Evaluate prediction data by substituting symbol 'n' in the prediction formula with the given ODE system size (n).
    if n == 0:
        return 0.0
    try:
        row = predictions.loc[(predictions['first'] <= n) & (predictions['last'] >= n)]
        prediction = eval_math_expr(str(row.prediction.iat[0]), [ivp_system_size(n)], cast_to=float)
    except IndexError:
        raise RuntimeError('Failed to find prediction data for n={}:\n{}'.format(n, predictions))
    return prediction


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
    args: Namespace = parse_program_args_app_plot()
    # Run plot app.
    plot(args)
