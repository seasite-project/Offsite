"""@package util.sample_interval
Definition of class SampleInterval.
"""

from enum import Enum
from sys import maxsize as sys_maxsize
from typing import List, Optional, Tuple

import attr

import offsite.config
from offsite.config import Config, ProgramModeType
from offsite.descriptions.ode import IVP, ivp_grid_size, ivp_system_size
from offsite.util.math_utils import eval_math_expr


class SampleType(Enum):
    """
    Defines what kinda of sample is if the sample lies in the border region of two adjacent intervals or not.

    - INNER
        Sample interval in the inner region of an interval.
    - BORDER
        Sample interval in the border region of two adjacent intervals.
    """
    MODEL_INNER = 'MODEL_INNER'
    MODEL_BORDER = 'MODEL_BORDER'
    BENCH_INNER = 'BENCH_INNER'
    BENCH_BORDER = 'BENCH_BORDER'
    ONLINE_FEEDBACK = 'ONLINE_FEEDBACK'


@attr.s(hash=True)
class SampleInterval:
    """Representation of a SampleInterval object.

    Attributes:
    -----------
    first: int
        First value included in the sample interval.
    last: int
        First value included in the sample interval.
    sample: int
        Actual value used to sample the interval.
    type: SampleType
        Type of the interval. See enum class SampleType.
    """
    first = attr.ib(type=int)
    last = attr.ib(type=int)
    type = attr.ib(type=SampleType)
    sample = attr.ib(type=int, default=None)

    def __attrs_post_init__(self):
        self.first = int(self.first)
        self.last = int(self.last)
        if self.sample is not None:
            self.sample = int(self.sample)

    def __hash__(self):
        return hash((self.first, self.last))

    def __eq__(self, other):
        return (self.first, self.last) == (other.first, other.last)

    def __ne__(self, other):
        return not self == other

    def median(self, ivp: Optional[IVP] = None) -> int:
        """Return median value of this interval.

        If an IVP is passed not the exact median is returned but the nearest value that satisfies the constraints of
        the IVP regarding possible system sizes(e.g. square system size for Heat2D, ...).

        Parameters:
        -----------
        ivp: IVP
            Used IVP.

        Returns:
        --------
        int
            Median value of this interval.
        """
        median = int(round(self.first + (self.last - self.first) / 2))
        # Consider restraints of the grid size of the IVP constraining possible but also sampleable, median values.
        if ivp is not None:
            # Round up to previous ...
            grid_size = int(eval_math_expr(ivp.gridSize, [ivp_system_size(median)]))
            median_low = int(eval_math_expr(ivp.ivp_size(), [ivp_grid_size(grid_size)]))
            if median_low < self.first or median_low > self.last:
                # Round up to next ...
                grid_size = int(round(eval_math_expr(ivp.gridSize, [ivp_system_size(median)])))
                median_high = int(round(eval_math_expr(ivp.ivp_size(), [ivp_grid_size(grid_size)])))
                if median_high < self.first or median_high > self.last:
                    median = None
                else:
                    median = median_high
            else:
                median = median_low
        return median


def create_samples_lower_border_working_set(
        wset_start: int, wset_end: int, num_samples: int, step: float, ivp: IVP) -> Tuple[List[SampleInterval], int]:
    """Create sample intervals in the upper border region of a working set.

    Parameters:
    -----------
    wset_start: int
        Smallest system size of the working set.
    wset_end: int
        Largest system size of the working set.
    num_samples: int
        Number of sample intervals created.
    step: float
        Step between created sample intervals.
    ivp: IVP
        Used IVP.

    Returns:
    --------
    list of SampleInterval
        List of sample intervals.
    int
        Highest system size above border region.
    """
    config: Config = offsite.config.offsiteConfig
    sample_type: SampleType
    if config.args.mode == ProgramModeType.MODEL:
        sample_type = SampleType.MODEL_BORDER
    elif config.args.mode == ProgramModeType.RUN:
        sample_type = SampleType.BENCH_BORDER
    else:
        assert False
    samples: List[SampleInterval] = list()
    # Ignore lower border for first working set. Too low values!
    start: int = wset_start
    if start == 1:
        return samples, start
    # Create samples in border region to previous working set.
    for i in range(1, num_samples + 1):
        end: int = min(int(wset_start * (1 + (i / step))), wset_end)
        point: int = SampleInterval(start, end, sample_type).median(ivp)
        if point:
            samples.append(SampleInterval(start, end, sample_type, point))
            # Switch to start of next interval.
            start = end + 1
        # Reached end of working set.
        if start > wset_end:
            break
    return samples, start


def create_samples_upper_border_working_set(
        wset_start: int, wset_end: int, num_samples: int, step: float, ivp: IVP) -> Tuple[List[SampleInterval], int]:
    """Create sample intervals and points in the upper border region of a working set.

    Parameters:
    -----------
    wset_start: int
        Smallest system size of the working set.
    wset_end: int
        Largest system size of the working set.
    num_samples: int
        Number of sample intervals created.
    step: float
        Step between created sample intervals.
    ivp: IVP
        Used IVP.

    Returns:
    --------
    list of SampleInterval
        List of sample intervals.
    int
        Highest system size below border region.
    """
    config: Config = offsite.config.offsiteConfig
    sample_type: SampleType
    if config.args.mode == ProgramModeType.MODEL:
        sample_type = SampleType.MODEL_BORDER
    elif config.args.mode == ProgramModeType.RUN:
        sample_type = SampleType.BENCH_BORDER
    else:
        assert False
    samples: List[SampleInterval] = list()
    end: int = wset_end
    # Create samples in border region to next working set.
    for i in range(1, num_samples + 1):
        start: int = max(int(wset_end * (1 - (i / step))), wset_start)
        point: int = SampleInterval(start, end, sample_type).median(ivp)
        if point:
            samples.append(SampleInterval(start, end, sample_type, point))
            # Switch to end of previous interval.
            end = start - 1
        # Switch to next working set once the end of the last interval from the lower border region was reached.
        if start < wset_start:
            break
    return samples, end


def create_samples_memory_lvl(
        start: int, num_samples: int, sample_offset: int, ivp: Optional[IVP] = None) -> List[SampleInterval]:
    """Create sample intervals and points in the memory level region.

    Parameters:
    -----------
    start: int
        Smallest system size in the memory region.
    num_samples: int
        Number of samples created.
    sample_offset: int
        Offset between samples.
    ivp: IVP
        Used IVP.

    Returns:
    --------
    list of SampleInterval
        List of sample intervals.
    """
    config: Config = offsite.config.offsiteConfig
    sample_type: SampleType
    if config.args.mode == ProgramModeType.MODEL:
        sample_type = SampleType.MODEL_INNER
    elif config.args.mode == ProgramModeType.RUN:
        sample_type = SampleType.BENCH_INNER
    else:
        assert False
    samples: List[SampleInterval] = list()
    # Create samples in memory region.
    intv_length: int = (start * sample_offset) - start + 1
    for _ in range(1, num_samples):
        # Determine end of sample interval.
        end: int = start + intv_length
        # Determine sample point.
        point = SampleInterval(start, end, sample_type).median(ivp)
        if point:
            samples.append(SampleInterval(start, end, sample_type, point))
        else:
            raise ValueError('Failed to create sample interval in main memory range!')
        # Switch to next sample.
        start = end + 1
    # Last sample in memory region.
    point = SampleInterval(start, start + intv_length, sample_type).median(ivp)
    if point:
        samples.append(SampleInterval(start, sys_maxsize, sample_type, point))
    else:
        raise ValueError('Failed to create sample interval in main memory range!')
    return samples
