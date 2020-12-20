"""@package ranking
Definitions of ranking functions.
"""

from multiprocessing import cpu_count, Pool
from sys import maxsize as sys_maxsize
from traceback import print_exc
from typing import List, Optional, Tuple, Union

from pandas import DataFrame
from sqlalchemy.orm import Session

from offsite.descriptions.impl_variant import ImplVariant
from offsite.descriptions.ivp import IVP
from offsite.descriptions.machine import Machine
from offsite.descriptions.ode_method import ODEMethod
from offsite.descriptions.ranking_task import RankDeviationTask, RankOrderTask
from offsite.evaluation.performance_model import SampleInterval
from offsite.evaluation.records import ImplVariantRecord, RankingRecord

RankingData = List[Tuple[SampleInterval, List[int]]]


def fuse_equal_rankings(rankings: RankingData) -> RankingData:
    """
    Reduce total number of rankings by combining neighbouring rankings. Adjacent rankings can be combined if they
    contain the same implementation variants.

    Parameters:
    -----------
    rankings: list of tuple(SampleInterval, set of impl variant ids)
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


def create_rankings(db_session: Session, machine: Machine, methods: List[ODEMethod], ivps: List[IVP],
                    rank_tasks: List[Union[RankDeviationTask, RankOrderTask]], ode_size: Optional[int] = None):
    """Rank implementation variants according to the given ranking tasks and store all rankings obtained in database.

    Parameters:
    -----------
    db_session: sqlalchemy.orm.session.Session
        Used database session.
    machine: Machine
        Used Machine.
    methods: List of ODEMethod.
        Used ODE methods.
    ivps: List of IVP
        Used IVPs.
    rank_tasks: List of RankTask
        Used ranking tasks for which ranking will be created and stored in database.
    ode_size: int
        Fixed ODE system size if given.

    Returns:
    --------
    -
    """
    # Machine frequency used.
    frequency: float = machine.clock
    # Max number of CPU cores used.
    max_cores: int = machine.coresPerSocket + 1
    # Retrieve implementation variants from database.
    impl_variants: List[ImplVariant] = ImplVariant.select_all(db_session)
    # Initialize worker thread pools.
    pool = Pool(cpu_count())
    # Rank variants for ...
    records: List[RankingRecord] = list()
    errors = list()
    # ... all ODE methods
    for method in methods:
        # ... all IVPS:
        for ivp in ivps:
            for task in rank_tasks:
                # Remove old, (possibly obsolete) ranking records from the database.
                if ode_size is not None:
                    # Fixed ODE system size.
                    RankingRecord.remove_record(
                        db_session, machine.db_id, machine.compiler.db_id, method.db_id, ivp.db_id, frequency,
                        [*range(1, max_cores)], ode_size, ode_size, task.type, task.cutoff_criteria, task.cutoff_value)
                else:
                    RankingRecord.remove_records(
                        db_session, machine.db_id, machine.compiler.db_id, method.db_id, ivp.db_id, frequency,
                        [*range(1, max_cores)], task.type, task.cutoff_criteria, task.cutoff_value)
                # ... all core counts.
                for cores in range(1, max_cores):
                    prediction_data = ImplVariantRecord.select(
                        db_session, [impl.db_id for impl in impl_variants], machine.db_id, machine.compiler.db_id,
                        method.db_id, ivp.db_id, frequency, cores, ode_size)
                    if prediction_data.empty:
                        pool.close()
                        pool.join()
                        raise RuntimeError('No fitting implementation variant records found in database.')
                    # Push rank task to worker pool.
                    pool.apply_async(rank_variants, args=(machine, method, ivp, cores, prediction_data, task, ode_size),
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


def rank_variants(machine: Machine, method: ODEMethod, ivp: IVP, cores: int, data: DataFrame,
                  rank_task: Union[RankDeviationTask, RankOrderTask], ode_size: Optional[int] = None):
    records = list()
    try:
        # Fixed ODE system size.
        if ode_size is not None:
            intv = SampleInterval(ode_size, ode_size)
            ranking = (intv, rank_task(data, ode_size))
            # Create ranking records.
            records.append(RankingRecord(machine.db_id, machine.compiler.db_id, method.db_id, ivp.db_id, cores,
                                         machine.clock, ranking[0], ranking[1], rank_task.type,
                                         rank_task.cutoff_criteria, rank_task.cutoff_value))
        else:
            # ... split prediction data by implementation variant ID.
            # Pass list [idx] here to ensure DataFrame is used instead of Series.
            prediction_data = {idx: data.loc[[idx]] for idx in data.index.unique().values}
            #
            ranking_data = list()
            start = 1
            end = 0
            while any(len(impl.index) > 1 for impl in prediction_data.values()):
                end = sys_maxsize
                # Determine lowest end value of current samples.
                low_sample = None
                for impl in prediction_data.values():
                    entry = impl.iloc[0]
                    start = entry['first']
                    end = min(entry['last'], end)
                    if start == end:
                        break
                # Add sample interval with this end value to 'rank_data'.
                sample = SampleInterval(start, end)
                #
                current_predictions = {idx: impl.iloc[0]['prediction'] for idx, impl in prediction_data.items()}
                current_predictions = DataFrame.from_dict(current_predictions, 'index', columns=['prediction'])
                # Rank variants.
                ranking_data.append((sample, rank_task(current_predictions, end)))
                # Update end of all current samples.
                for impl in prediction_data.values():
                    impl.iat[0, 0] = end + 1
                # Remove all current samples already represented in the list of significant samples.
                for impl_id, impl in prediction_data.items():
                    entry = impl.iloc[0]
                    if len(impl.index) > 1 and entry['first'] > entry['last']:
                        prediction_data[impl_id] = impl[1:]
            # Add interval with 'inf' end.
            sample = SampleInterval(end + 1, sys_maxsize)
            #
            current_predictions = {idx: impl.iloc[0]['prediction'] for idx, impl in prediction_data.items()}
            current_predictions = DataFrame.from_dict(current_predictions, 'index', columns=['prediction'])

            # Rank variants.
            ranking_data.append((sample, rank_task(current_predictions, end)))
            # Reduce number of rankings by fusing adjacent rankings if possible.
            ranking_data = fuse_equal_rankings(ranking_data)
            # Create ranking records.
            for ranking in ranking_data:
                records.append(RankingRecord(machine.db_id, machine.compiler.db_id, method.db_id, ivp.db_id, cores,
                                             machine.clock, ranking[0], ranking[1], rank_task.type,
                                             rank_task.cutoff_criteria, rank_task.cutoff_value))
        return records
    except Exception as e:
        # Print stack trace of the executing worker thread.
        print_exc()
        print('')
        raise e
