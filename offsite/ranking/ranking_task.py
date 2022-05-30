"""@package ranking.ranking_task
Definitions of classes RankTask, RankOrderTask and RankPercentTask.

@author: Johannes Seiferth
"""

from enum import Enum
from typing import List, Optional

import attr
from pandas import DataFrame

import offsite.config
from offsite.descriptions.ode import ivp_system_size
from offsite.util.math_utils import eval_math_expr, percent_deviation


class RankingCriteriaType(Enum):
    """
    Defines whether implementation variants are ranked by ascending runtime or performance deviation from the best
    rated implementation variant.

    - ORDER
        Rank implementation variants by their runtime.
    - DEVIATION
        Rank implementation variants by their performance deviation from the best rated variant.
    """
    ORDER = 'ORDER'
    DEVIATION = 'DEVIATION'


class RankingCutoffType(Enum):
    """Defines what kind of cutoffValue is used to determine which implementation variants are included in the ranking.

    - TOP
        Include only a fixed number of all implementation variants:
         - e.g. top 20 --> 20 best variants
    - PERCENT
        Include only a fixed percentage of all implementation variants:
         - e.g. 20% --> best 20 % of all variants for cutoffCriteria ORDER.
         - e.g. 20% --> all variants within 20 % deviation of the best variants for cutoffCriteria DEVIATION.
    """
    TOP = 'TOP'
    PERCENT = 'PERCENT'


@attr.s
class RankTask:
    """Representation of a ranking task.

    Attributes:
    -----------
    cutoff_criteria: RankingCutoffType
        Used cutoff criteria.
    cutoff_value: float
        Used cutoff value.
    """
    cutoff_criteria = attr.ib(type=RankingCutoffType)
    cutoff_value = attr.ib(type=float)


@attr.s
class RankOrderTask(RankTask):
    """Representation of a "ranking by variant order" task.

    Attributes:
    -----------
    cutoff_criteria: RankingCutoffType
        Used cutoff criteria.
    cutoff_value: float
        Used cutoff value.
    type: RankingCriteriaType
        Used ranking criteria.
    """
    type = attr.ib(type=RankingCriteriaType, init=False, default=RankingCriteriaType.ORDER)

    def __call__(self, data: DataFrame, system_size: int) -> List[int]:
        if offsite.config.offsiteConfig.args.verbose:
            print('Creating ranking by order using {} for {}{}\n'.format(
                self.cutoff_criteria, self.cutoff_value,
                '%' if self.cutoff_criteria == RankingCutoffType.PERCENT else ''))
        ranking: DataFrame = sort_variants_by_performance(data, system_size)
        if self.cutoff_criteria == RankingCutoffType.TOP:
            # Return top [cutoffValue] of all implementation variants.
            num_variants = int(self.cutoff_value)
        elif self.cutoff_criteria == RankingCutoffType.PERCENT:
            # Return the [cutoffValue] percent of all implementation variants.
            num_variants = int(len(ranking) * (self.cutoff_value / 100.0))
        else:
            assert False
        return [idx for idx, impl in ranking.head(num_variants).iterrows()]


@attr.s
class RankDeviationTask(RankTask):
    """Representation of a "ranking by variant deviation from the best variant" task.

    Attributes:
    -----------
    cutoff_criteria: RankingCutoffType
        Used cutoff criteria.
    cutoff_value: float
        Used cutoff value.
    type: RankingCriteriaType
        Used ranking criteria.
    """
    type = attr.ib(type=RankingCriteriaType, init=False, default=RankingCriteriaType.DEVIATION)

    def __call__(self, data: DataFrame, system_size: int) -> List[int]:
        if offsite.config.offsiteConfig.args.verbose:
            print('Creating ranking by deviation using {} for {}{}\n'.format(
                self.cutoff_criteria, self.cutoff_value,
                '%' if self.cutoff_criteria == RankingCutoffType.PERCENT else ''))
        ranking: DataFrame = sort_variants_by_performance(data, system_size)
        if self.cutoff_criteria == RankingCutoffType.TOP:
            # Return top [cutoffValue] implementation variants.
            return [idx for idx, impl in ranking.head(int(self.cutoff_value)).iterrows()]
        elif self.cutoff_criteria == RankingCutoffType.PERCENT:
            # Return all implementation variants within the tolerance.
            index_best_impl = ranking.head(1).index.values[0]
            return [idx for idx, impl in ranking.iterrows() if (
                percent_deviation(impl['prediction'], ranking.at[index_best_impl, 'prediction'])) <= self.cutoff_value]
        else:
            assert False


def sort_variants_by_performance(data: DataFrame, system_size: int) -> DataFrame:
    """Sort implementation variants by ascending runtime prediction.

    Parameters:
    -----------
    data: DataFrame
        Implementation variant prediction data.
    system_size: int
        Rank variants for this system size.

    Returns:
    -------
    DataFrame
        Implementation variants ranked by ascending runtime prediction.
    """
    # Evaluate predictions for the given ODE system size.
    constants = [ivp_system_size(system_size), ('x', system_size)]
    for idx, row in data.iterrows():
        data.at[idx, 'prediction'] = eval_math_expr(row['prediction'], constants, cast_to=float)
    # Rank implementation variants by ascending runtime prediction.
    ranking = data.sort_values(by=['prediction'])
    return ranking


def parse_ranking_task(task: str) -> Optional[RankTask]:
    """Parse the ranking task string provided and create corresponding RankTask object.

    Parameters:
    -----------
    tasks_str: str
        Used ranking task string.

    Returns:
    --------
    RankTask
        Created ranking task.
    """
    verbose: bool = offsite.config.offsiteConfig.args.verbose
    # Parse task string.
    task_args: List[str] = task.split(':')
    if len(task_args) != 3:
        if verbose:
            print('Skipping invalid ranking task: \'{}\'!'.format(task))
        return None
    # Check cutoff criteria.
    if task_args[1] not in ['t', 'p']:
        if verbose:
            print('Skipping invalid ranking task: \'{}\'! Unknown cutoff criteria \'{}\'!'.format(task, task_args[1]))
        return None
    cutoff_crit: RankingCutoffType = RankingCutoffType.TOP if task_args[1] == 't' else RankingCutoffType.PERCENT
    # Check cutoff value.
    cutoff_val: float = 0.0
    try:
        cutoff_val = float(task_args[2])
    except ValueError:
        if verbose:
            print('Skipping invalid ranking task: \'{}\'! Invalid cutoff value \'{}\'!'.format(task, cutoff_val))
        return None
    # Create ranking task and append to list.
    rank_crit: str = task_args[0]
    if rank_crit not in ['o', 'd']:
        if verbose:
            print('Skipping invalid ranking task: \'{}\'! Unknown ranking criteria \'{}\'!'.format(task, rank_crit))
        return None
    return RankOrderTask(cutoff_crit, cutoff_val) if rank_crit == 'o' else RankDeviationTask(cutoff_crit, cutoff_val)


def parse_ranking_tasks(rank_task_str: str) -> List[RankTask]:
    """Parse the ranking task string provided and create corresponding RankTask objects.

    Parameters:
    -----------
    tasks_str: str
        Used ranking task string.

    Returns:
    --------
    list of RankTask
        Created ranking tasks.
    """
    # Parse single ranking task and add to ranking task list.
    tasks: List[RankTask] = list()
    for task_str in rank_task_str.split(','):
        task = parse_ranking_task(task_str)
        if task is not None:
            tasks.append(task)
    return tasks
