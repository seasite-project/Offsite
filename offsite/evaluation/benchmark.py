"""@package benchmark
Definition of class OmpBarrierBenchmark.
"""

from pathlib import Path
from subprocess import run, PIPE, CalledProcessError
from typing import List

import attr

from offsite.descriptions.machine import Machine
from offsite.evaluation.math_utils import remove_outliers
from offsite.evaluation.records import BenchmarkRecord


@attr.s(hash=True)
class OmpBarrierBenchmark:
    """Representation of the OmpBarrier benchmark.

    Attributes:
    -----------
    name: str
        Name of this benchmark.
    folder: Path
        Folder to store source code in.
    source: Path
        Path to source code file associated with this benchmark.
    parameters: set of str
        Parameters of this benchmark.
    binary: Path
        Path to binary object file of this benchmark compiled for a given
        configuration of machine and compiler.
    compiled_for: Machine
        Machine this benchmark was compiled for last.

    """
    name = attr.ib(type=str, default='omp_barrier')
    folder = attr.ib(type=Path, default=Path('tmp'))
    source = attr.ib(type=Path)
    binary = attr.ib(type=Path, default=Path())
    compiled_for = attr.ib(type=Machine, default=None)

    def __del__(self):
        """Destructor that removes all files created during benchmarking.

        Parameters:
        -----------
        -

        Returns:
        --------
        -
        """
        if self.source.exists() and self.source.is_file():
            self.source.unlink()
        if self.binary.exists() and self.binary.is_file():
            self.binary.unlink()

    def run(self, machine: Machine, repetitions: int, save_in_db=True) -> List[BenchmarkRecord]:
        """Run the omp_barrier benchmark on a given machine.

        Parameters:
        -----------
        machine: Machine
            Machine the benchmark is ran on.
        repetitions: int
            Number of times this benchmark is executed.
        save_in_db: bool
            If True create valid BenchmarkRecord database entries.

        Returns:
        --------
        list of BenchmarkRecord
            Benchmark data obtained on the given machine.
        """
        data: List[BenchmarkRecord] = list()
        # Compile benchmark for given machine.
        if self.compiled_for is not machine or not self.binary_exists():
            self.compile(machine)
        # Run benchmark for all core numbers.
        for cores in range(1, machine.coresPerSocket + 1):
            times: List[float] = list()
            try:
                for _ in range(1, repetitions + 1):
                    cmd = ['./{}'.format(self.binary), str(cores), str(repetitions)]
                    # Run benchmark for current configuration.
                    runtime = run(cmd, check=True, encoding='utf-8', stdout=PIPE).stdout
                    times.append(float(runtime))
            except CalledProcessError as error:
                print('OmpBarrierBenchmark failed: {}'.format(error))
            # Remove outliers.
            times: List[float] = remove_outliers(times)
            # Runtime as mean of all times measured.
            runtime: float = sum(times) / len(times)
            # Save benchmark record.
            record = BenchmarkRecord(self.name, machine.db_id if save_in_db else -1,
                                     machine.compiler.db_id if save_in_db else -1, runtime, machine.clock, cores)
            data.append(record)
        return data

    def compile(self, machine: Machine):
        """Compile source file of the omp_barrier benchmark on a given machine.

        Parameters:
        -----------
        machine: Machine
            Machine the benchmark is compiled for.

        Returns:
        --------
        -
        """
        self.binary = self.folder / Path('{}.out'.format(self.name))
        self.compiled_for = machine
        try:
            # Construct compiler call.
            args = [machine.compiler.name, str(self.source), '-fopenmp', '-o{}'.format(self.binary)]
            # Add compiler flags to call.
            args.extend((flag for flag in machine.compiler.flags.split(' ')))
            # Compile benchmark.
            run(args, check=True)
        except CalledProcessError as error:
            raise RuntimeError('Failed to compile {}: {}'.format(self.name, error))

    def binary_exists(self) -> bool:
        """Check if binary file of the omp_barrier benchmark exists.

        Parameters:
        -----------
        -

        Returns:
        --------
        bool
            True if binary file exists else False.
        """
        return False if not self.binary else self.binary.is_file()

    @source.default
    def write_code_to_file(self) -> Path:
        """Write source code of the omp_barrier benchmark to file.

        Parameters:
        -----------
        -

        Returns:
        --------
        pathlib.Path
            Relative path to source file created.
        """
        code = '\
#include <omp.h>\n\
#include <stdlib.h>\n\
#include <stdio.h>\n\
#include <sys/time.h>\n\
double times;\n\
double repetitions;\n\
int p, k;\n\
int threads;\n\
double *timer;\n\
#define start_timer(t) (t = gtod())\n\
#define stop_timer(t) (t = gtod() - t)\n\
const int micro_to_sec = 1E6;\n\
double gtod()\n\
{\n\
  struct timeval tv;\n\
  gettimeofday(&tv, NULL);\n\
  return (double) tv.tv_sec + (double) tv.tv_usec / micro_to_sec;\n\
}\n\
int compare(const void *a, const void *b)\n\
{\n\
  if (*(double*)a > *(double*)b) return 1;\n\
  else if (*(double*)a < *(double*)b) return -1;\n\
  else return 0;\n\
}\n\
int main(int argc, char *argv[])\n\
{\n\
  int p = 1;\n\
  double max_timer = -1.0;\n\
  if (argc < 3) {\n\
    return -1;\n\
  }\n\
  char *end;\n\
  threads = strtol(argv[1], &end, 10);\n\
  repetitions = strtol(argv[2], &end, 10);\n\
  omp_set_dynamic(0);\n\
  omp_set_num_threads(threads);\n\
#pragma omp parallel\n\
  {\n\
#pragma omp single\n\
    {\n\
      p = omp_get_num_threads();\n\
    }\n\
  }\n\
  timer = (double *) malloc(p * sizeof(double));\n\
#pragma omp parallel\n\
  {\n\
    int i, j;\n\
    int me = omp_get_thread_num();\n\
    double t;\n\
    start_timer(t);\n\
    /// MEASURED REGION START ///\n\
    for (j = 0; j < repetitions; ++j)\n\
    {\n\
#pragma omp barrier\n\
      ;\n\
    }\n\
    stop_timer(t);\n\
    /// MEASURED REGION END ///\n\
    timer[me] = t;\n\
  }\n\
  for (k = 0; k < p; ++k)\n\
  {\n\
    if (max_timer < timer[k])\n\
      max_timer = timer[k];\n\
  }\n\
  max_timer /= (double) repetitions;\n\
  printf("%e", max_timer );\n\
  free(timer);\n\
  return 0;\n\
}\n'
        # Create folder if it does not yet exist.
        if self.folder and not self.folder.exists():
            self.folder.mkdir(parents=True)
        # Write source code file.
        path = Path('{}.{}'.format(self.name, 'c'))
        path = self.folder / path
        with path.open('w') as file_handle:
            file_handle.write(code)
            return Path(file_handle.name)


# Available predefined benchmark functions.
AVAIL_BENCHMARKS = {
    'omp_barrier': OmpBarrierBenchmark()
}
