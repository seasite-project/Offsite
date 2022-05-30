"""@package util.time
Util time functions.

@author: Johannes Seiferth
"""

from datetime import timedelta
from time import time


def start_timer() -> float:
    return time()


def stop_timer(start_time: float) -> timedelta:
    return timedelta(seconds=round(time() - start_time, 0))
