"""@package train.communication
Training functions for communication patterns.
"""

from offsite.train.communication.openmp.omp_barrier import OmpBarrierBenchmark

# Available predefined benchmark functions.
AVAIL_BENCHMARKS = {
    'omp_barrier': OmpBarrierBenchmark()
}
