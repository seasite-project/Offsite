"""@package ranking
Definitions of ranking functions.
"""

from multiprocessing import cpu_count, Pool
from sys import maxsize as sys_maxsize
from traceback import print_exc
from typing import List, Tuple, Set

from pandas import DataFrame
from sqlalchemy.orm import Session

import offsite.config
from offsite.descriptions.impl_variant import ImplVariant
from offsite.descriptions.ivp import IVP
from offsite.descriptions.machine import Machine
from offsite.descriptions.ode_method import ODEMethod
from offsite.evaluation.math_utils import eval_math_expr, ivp_system_size, percent_deviation
from offsite.evaluation.performance_model import SampleInterval
from offsite.evaluation.records import ImplVariantRecord, RankingRecord

RankingData = List[Tuple[SampleInterval, Set[int]]]


def deduce_best_variants(data: 'pandas.DataFrame', system_size: float) -> Set[int]:
    """
    Rank implementation variants by ascending runtime prediction and return the set of ids of those implementation
    variants that are within a given tolerance from the best predicted implementation variant.

    Parameters:
    -----------
    data : dict (key: impl variant id, value: ImplVariantRecord)
        Implementation variant prediction data.
    system_size : float
        Rank variants for this system size.

    Returns:
    -------
    Set of int
        Set of implementation variant ids ranked by ascending runtime prediction.
    """
    config = offsite.config.offsiteConfig
    # Evaluate predictions for the given ODE system size.
    constants = [ivp_system_size(system_size), ('x', system_size)]
    for idx, row in data.iterrows():
        data.at[idx, 'prediction'] = eval_math_expr(row[2], constants, cast_to=float)
    # Rank implementation variants by ascending runtime prediction.
    ranking = data.sort_values(by=['prediction'])
    # Return all implementation variants within the tolerance.
    index_best_impl = ranking.head(1).index.values[0]
    return set([idx for idx, impl in ranking.iterrows() if (
        percent_deviation(impl['prediction'], ranking.at[index_best_impl, 'prediction'])) <= config.args.tolerance])


def fuse_equal_rankings(rankings: RankingData) -> RankingData:
    """
    Reduce total number of rankings by combining neighbouring rankings. Adjacent rankings can be combined if they
    contain the same implementation variants.

    Parameters:
    -----------
    rankings : list of tuple(SampleInterval, set of impl variant ids)
        Rankings of implementation variants for ascending intervals.

    Returns:
    --------
    list of tuple(SampleInterval, set of impl variant ids)
        Rankings of implementation variants for ascending intervals after fusing.
    """
    if len(rankings) <= 1:
        return rankings
    # Iterate over the rankings and combine those that contain the same implementation variants.
    combined = list()
    r_prev = set()
    intv = None
    for ranking in rankings:
        i_cur = ranking[0]
        r_cur = ranking[1]
        if r_prev == r_cur:
            # Increase upper bound of current interval.
            intv.last = i_cur.last
        else:
            # Save interval.
            if intv:
                combined.append((intv, r_prev))
            # Start new interval.
            r_prev = r_cur
            intv = SampleInterval(i_cur.first, i_cur.last)
    # Save final interval.
    combined.append((intv, r_prev))
    return combined


def rank(db_session: Session, machine: Machine, methods: List[ODEMethod], ivps: List[IVP]):
    """
    Rank implementation variants by ascending runtime prediction and return the set of ids of those implementation
    variants that are within a given tolerance from the best predicted implementation variant.

    Parameters:
    -----------
    db_session : sqlalchemy.orm.session.Session
        Used database session.
    machine : Machine
        Used Machine.
    methods : List of ODEMethod.
        Used ODE methods.
    ivps : List of IVP
        Used IVPs.

    Returns:
    --------
    -
    """
    config = offsite.config.offsiteConfig
    # Machine frequency used.
    frequency = machine.clock
    # Max number of CPU cores used.
    max_cores = machine.coresPerSocket + 1
    # Retrieve implementation variants from database.
    impl_variants = ImplVariant.select_all(db_session)
    # Initialize worker thread pools.
    pool = Pool(cpu_count())
    # Rank variants for ...
    records = list()
    errors = list()
    # ... all ODE methods
    for method in methods:
        # ... all IVPS:
        for ivp in ivps:
            if config.args.tool is not ivp.modelTool:
                continue
            # Remove old, (possibly obsolete) ranking records from the database.
            # Fixed ODE system size.
            if config.args.ode_size:
                RankingRecord.remove_record(
                    db_session, machine.db_id, machine.compiler.db_id, method.db_id, ivp.db_id, frequency,
                    [*range(1, max_cores)], config.args.ode_size, config.args.ode_size)
            else:
                RankingRecord.remove_records(db_session, machine.db_id, machine.compiler.db_id, method.db_id, ivp.db_id,
                                             frequency, [*range(1, max_cores)])
            # ... all core counts.
            for cores in range(1, max_cores):
                prediction_data = ImplVariantRecord.select(
                    db_session, [impl.db_id for impl in impl_variants], machine.db_id, machine.compiler.db_id,
                    method.db_id, ivp.db_id, frequency, cores)
                # Push rank task to worker pool.
                pool.apply_async(rank_variants, args=(machine, method, ivp, cores, prediction_data),
                                 callback=records.extend, error_callback=errors.append)
    # Wait for all threads and collect results.
    pool.close()
    pool.join()
    # Raise error if ranking failed.
    if errors:
        db_session.rollback()
        raise RuntimeError('Failed to rank implementation variants: Error in worker threads.')
    # Train database with ranking records.
    RankingRecord.update(db_session, records)


def rank_variants(machine: Machine, method: ODEMethod, ivp: IVP, cores: int, data: 'pandas.DataFrame'):
    config = offsite.config.offsiteConfig
    records = list()
    try:
        # Fixed ODE system size.
        if config.args.ode_size:
            ranking = (SampleInterval(config.args.ode_size, config.args.ode_size),
                       deduce_best_variants(data, config.args.ode_size))
            # Create ranking records.
            records.append(RankingRecord(machine.db_id, machine.compiler.db_id, method.db_id, ivp.db_id, cores,
                                         machine.clock, ranking[0], ranking[1]))
        else:
            # ... split prediction data by implementation variant ID.
            prediction_data = {idx: data.loc[idx] for idx in data.index.unique().values}
            #
            ranking_data = list()
            start = 1
            end = 0
            while any(len(impl) > 1 for impl in prediction_data.values()):
                end = sys_maxsize
                # Determine lowest end value of current samples.
                low_sample = None
                for impl in prediction_data.values():
                    first_row = [x for x in impl.head(1).iterrows()][0][1]
                    start = first_row['first']
                    low_sample = first_row if first_row['last'] < end else low_sample
                    end = min(first_row['last'], end)
                    if start == end:
                        break
                # Add sample interval with this end value to 'rank_data'.
                sample = SampleInterval(start, end)
                #
                current_predictions = dict()
                for idx, impl in prediction_data.items():
                    current_predictions[idx] = [x for x in impl.head(1).iterrows()][0][1]
                current_predictions = DataFrame.from_dict(current_predictions, 'index')
                # Rank variants.
                ranking_data.append((sample, deduce_best_variants(current_predictions, end)))
                # Update end of all current samples.
                for impl in prediction_data.values():
                    impl.iat[0, 0] = end + 1
                # Remove all current samples already represented in the list of significant samples.
                for impl_id, impl in prediction_data.items():
                    first_row = [x for x in impl.head(1).iterrows()][0][1]
                    if len(impl) > 1 and first_row['first'] > first_row['last']:
                        prediction_data[impl_id] = impl[1:]
            # Add interval with 'inf' end.
            sample = SampleInterval(end + 1, sys_maxsize)
            #
            current_predictions = dict()
            for idx, impl in prediction_data.items():
                current_predictions[idx] = [x for x in impl.head(1).iterrows()][0][1]
            current_predictions = DataFrame.from_dict(current_predictions, 'index')
            # Rank variants.
            ranking_data.append((sample, deduce_best_variants(current_predictions, end)))
            # Reduce number of rankings by fusing adjacent rankings if possible.
            ranking_data = fuse_equal_rankings(ranking_data)
            # Create ranking records.
            for ranking in ranking_data:
                records.append(RankingRecord(machine.db_id, machine.compiler.db_id, method.db_id, ivp.db_id, cores,
                                             machine.clock, ranking[0], ranking[1]))
        return records
    except Exception as e:
        # Print stack trace of the executing worker thread.
        print_exc()
        print('')
        raise e
