"""@package util.sample_interval
Definition of class SampleInterval.

@author: Johannes Seiferth
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
        wset_start: int, wset_end: int, num_samples: int, ivp: IVP) -> Tuple[List[SampleInterval], int]:
    """Create sample intervals in the lower 20 percent border region of a working set.

    Parameters:
    -----------
    wset_start: int
        The smallest system size of the working set.
    wset_end: int
        The largest system size of the working set.
    num_samples: int
        Number of sample intervals created.
    ivp: IVP
        Used IVP.

    Returns:
    --------
    list of SampleInterval
        List of sample intervals.
    int
        The highest system size above border region.
    """
    assert wset_start <= wset_end
    config: Config = offsite.config.offsiteConfig
    sample_type: SampleType
    if config.args.mode == ProgramModeType.MODEL:
        sample_type = SampleType.MODEL_BORDER
    elif config.args.mode == ProgramModeType.RUN:
        sample_type = SampleType.BENCH_BORDER
    else:
        assert False
    # Determine the end of the lower 20 percent interval and the length of the interval.
    intv_20p: int = int(wset_start + 0.2 * (wset_end - wset_start + 1))
    len_lower_20p: int = intv_20p - wset_start + 1

    # The number of samples can not be bigger than the length of the interval to be sampled.
    num_samples: int = min(num_samples, len_lower_20p)
    # Determine the length of a sample interval ...
    len_intv: int = int(len_lower_20p / num_samples)
    # ... and the threshold to increase the length for the last 'longer_intv' samples.
    longer_intv = len_lower_20p - num_samples * len_intv

    # Determine 'start' and 'end' of the first sample interval...
    end: int = intv_20p - 1
    start: int = intv_20p - len_intv
    # ... and create all intervals.
    samples: List[SampleInterval] = list()
    for i in reversed(range(1, num_samples + 1)):
        # Create sample interval and sample point.
        point: int = SampleInterval(start, end, sample_type).median(ivp)
        if point:
            samples.append(SampleInterval(start, end, sample_type, point))
            # Switch to next interval.
            end = start - 1
        # Increase the length of the next interval.
        if i == longer_intv:
            len_intv += 1
        start -= len_intv
    samples.reverse()
    return samples, intv_20p


def create_samples_upper_border_working_set(
        wset_start: int, wset_end: int, num_samples: int, ivp: IVP) -> Tuple[List[SampleInterval], int]:
    """Create sample intervals and points in the upper 20 percent border region of a working set.

    Parameters:
    -----------
    wset_start: int
        The smallest system size of the working set.
    wset_end: int
        The largest system size of the working set.
    num_samples: int
        Number of sample intervals created.
    ivp: IVP
        Used IVP.

    Returns:
    --------
    list of SampleInterval
        List of sample intervals.
    int
        The highest system size below border region.
    """
    assert wset_start <= wset_end
    config: Config = offsite.config.offsiteConfig
    sample_type: SampleType
    if config.args.mode == ProgramModeType.MODEL:
        sample_type = SampleType.MODEL_BORDER
    elif config.args.mode == ProgramModeType.RUN:
        sample_type = SampleType.BENCH_BORDER
    else:
        assert False
    # Determine the start of the upper 20 percent interval and the length of the interval.
    intv_80p: int = int(wset_start + 0.8 * (wset_end - wset_start + 1))
    len_upper_20perc: int = wset_end - intv_80p + 1

    # The number of samples can not be bigger than the length of the interval to be sampled.
    num_samples: int = min(num_samples, len_upper_20perc)
    # Determine the length of a sample interval ...
    len_intv: int = int(len_upper_20perc / num_samples)
    # ... and the threshold to increase the length for the last 'longer_intv' samples.
    longer_intv = len_upper_20perc - num_samples * len_intv

    # Determine 'start' and 'end' of the first sample interval...
    start: int = intv_80p + 1
    end: int = intv_80p + len_intv
    # ... and create all intervals.
    samples: List[SampleInterval] = list()
    for i in reversed(range(1, num_samples + 1)):
        # Create sample interval and sample point.
        point: int = SampleInterval(start, end, sample_type).median(ivp)
        if point:
            samples.append(SampleInterval(start, end, sample_type, point))
            # Switch to next interval.
            start = end + 1
        # Increase the length of the next interval.
        if i == longer_intv:
            len_intv += 1
        end += len_intv
    return samples, intv_80p


def create_samples_border_memory_lvl(
        wset_start: int, wset_end: int, num_samples: int, ivp: IVP) -> Tuple[List[SampleInterval], int]:
    """Create sample intervals in the lower 20 percent border region of a working set.

    Parameters:
    -----------
    wset_start: int
        The smallest system size of the working set.
    wset_end: int
        The largest system size of the working set.
    num_samples: int
        Number of sample intervals created.
    ivp: IVP
        Used IVP.

    Returns:
    --------
    list of SampleInterval
        List of sample intervals.
    int
        The highest system size above border region.
    """
    assert wset_start <= wset_end
    config: Config = offsite.config.offsiteConfig
    sample_type: SampleType
    if config.args.mode == ProgramModeType.MODEL:
        sample_type = SampleType.MODEL_BORDER
    elif config.args.mode == ProgramModeType.RUN:
        sample_type = SampleType.BENCH_BORDER
    else:
        assert False
    # Determine the end of the lower 20 percent interval and the length of the interval.
    end_border_region: int = int(wset_start * 1.5)
    if end_border_region > 0.05 * wset_end - wset_start + 1:
        # Border region should not be bigger than 5% of the memory level.
        end_border_region: int = int(0.05 * wset_end)
    len_border_region = end_border_region - wset_start + 1

    # The number of samples can not be bigger than the length of the interval to be sampled.
    num_samples: int = min(num_samples, len_border_region)
    # Determine the length of a sample interval ...
    len_intv: int = int(len_border_region / num_samples)
    # ... and the threshold to increase the length for the last 'longer_intv' samples.
    longer_intv = len_border_region - num_samples * len_intv

    # Determine 'start' and 'end' of the first sample interval...
    end: int = end_border_region - 1
    start: int = end_border_region - len_intv
    # ... and create all intervals.
    samples: List[SampleInterval] = list()
    for i in reversed(range(1, num_samples + 1)):
        # Create sample interval and sample point.
        point: int = SampleInterval(start, end, sample_type).median(ivp)
        if point:
            samples.append(SampleInterval(start, end, sample_type, point))
            # Switch to next interval.
            end = start - 1
        # Increase the length of the next interval.
        if i == longer_intv:
            len_intv += 1
        start -= len_intv
    samples.reverse()
    return samples, end_border_region


def create_samples_memory_lvl(
        start: int, num_samples: int, sample_offset: int, ivp: Optional[IVP] = None) -> List[SampleInterval]:
    """Create sample intervals and points in the memory level region.

    Parameters:
    -----------
    start: int
        The smallest system size in the memory region.
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


def derive_samples_from_range_expr(intv_range: str, mode: ProgramModeType):
    if mode == ProgramModeType.MODEL:
        sample_type = SampleType.MODEL_INNER
    elif mode == ProgramModeType.RUN:
        sample_type = SampleType.BENCH_INNER
    else:
        assert False

    intervals = list()
    # Split interval range string.
    intv_range = intv_range.split(':')
    assert (len(intv_range) == 3)
    cur_sample = int(intv_range[0])
    last_sample = int(intv_range[1])
    incr = int(intv_range[2])
    # First interval starts with first sample.
    e = cur_sample + int(incr / 2)
    intervals.append(SampleInterval(cur_sample, e, sample_type, cur_sample))
    # Add inner intervals.
    cur_sample += incr
    s = e + 1
    e += incr
    while cur_sample < last_sample:
        intervals.append(SampleInterval(s, e, sample_type, cur_sample))
        cur_sample += incr
        s = e + 1
        e = min(e + incr, last_sample)
    # Add last interval.
    if s <= last_sample:
        intervals.append(SampleInterval(s, last_sample, sample_type, last_sample))
    return intervals
