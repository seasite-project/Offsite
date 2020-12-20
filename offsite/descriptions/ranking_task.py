"""@package ranking_task
Definitions of classes RankTask, RankOrderTask and RankPercentTask.
"""

from typing import List

import attr
from pandas import DataFrame

import offsite.config
from offsite.config import RankingCriteriaType, RankingCutoffType
from offsite.evaluation.math_utils import eval_math_expr, ivp_system_size
from offsite.evaluation.math_utils import percent_deviation


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
    """Representation of a ranking by variant order task.

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
    """Representation of a ranking by variant deviation from best variant task.

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
