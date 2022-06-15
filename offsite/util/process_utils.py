"""@package util.process_utils
Util functions to start subprocesses.

@author: Johannes Seiferth
"""

from subprocess import run, PIPE, CalledProcessError
from typing import List


def run_process(cmd: List[str]) -> str:
    try:
        output = run(cmd, check=True, stdout=PIPE, encoding='utf-8').stdout
    except CalledProcessError as error:
        raise RuntimeWarning('Process exited with error message: {}'.format(error))
    return output
