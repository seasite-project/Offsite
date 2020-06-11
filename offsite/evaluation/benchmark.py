"""@package benchmark
Definition of class OmpBarrierBenchmark.
"""

from datetime import datetime
from getpass import getuser
from pathlib import Path
from subprocess import run, PIPE, CalledProcessError
from sys import version_info
from typing import List

import attr
from pandas import read_sql_query
from sqlalchemy import Column, DateTime, Float, ForeignKey, Integer, String, Table
from sqlalchemy.orm import Session
from sqlalchemy.orm.exc import MultipleResultsFound, NoResultFound

from offsite import __version__
from offsite.db import METADATA
from offsite.db.db import insert
from offsite.descriptions.machine import Machine
from offsite.evaluation.math_utils import remove_outliers
from offsite.evaluation.performance_model import KernelRecord


@attr.s(hash=True)
class OmpBarrierBenchmark:
    """Representation of the OmpBarrier benchmark.

    Attributes:
    -----------
    name : str
        Name of this benchmark.
    folder : Path
        Folder to store source code in.
    source : Path
        Path to source code file associated with this benchmark.
    parameters : set of str
        Parameters of this benchmark.
    binary : Path
        Path to binary object file of this benchmark compiled for a given
        configuration of machine and compiler.
    compiled_for : Machine
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

    def run(self, machine: Machine, repetitions: int, save_in_db=True) -> List['BenchmarkRecord']:
        """Run the omp_barrier benchmark on a given machine.

        Parameters:
        -----------
        machine : Machine
            Machine the benchmark is ran on.
        repetitions : int
            Number of times this benchmark is executed.
        save_in_db : bool
            If True create valid BenchmarkRecord database entries.

        Returns:
        --------
        list of BenchmarkRecord
            Benchmark data obtained on the given machine.
        """
        data = list()
        # Compile benchmark for given machine.
        if self.compiled_for is not machine or not self.binary_exists():
            self.compile(machine)
        # Run benchmark for all core numbers.
        for cores in range(1, machine.coresPerSocket + 1):
            times = list()
            try:
                for _ in range(1, repetitions + 1):
                    if version_info[1] > 5:
                        runtime = run(['./{}'.format(self.binary), str(cores), str(repetitions)], check=True,
                                      encoding='utf-8', stdout=PIPE).stdout
                    else:
                        runtime = run(['./{}'.format(self.binary), str(cores), str(repetitions)], check=True,
                                      stdout=PIPE).stdout
                        runtime = runtime.decode("utf-8")
                    times.append(float(runtime))
            except CalledProcessError as error:
                print('OmpBarrierBenchmark failed: {}'.format(error))
            # Remove outliers.
            times = remove_outliers(times)
            # Runtime as mean of all times measured.
            runtime = sum(times) / len(times)
            # Save benchmark record.
            record = BenchmarkRecord(self.name, machine.db_id if save_in_db else -1,
                                     machine.compiler.db_id if save_in_db else -1, runtime, machine.clock, cores)
            data.append(record)
        return data

    def compile(self, machine: Machine):
        """Compile source file of the omp_barrier benchmark on a given machine.

        Parameters:
        -----------
        machine : Machine
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


@attr.s
class BenchmarkRecord:
    """Representation of a benchmark table database record.

    Attributes:
    -----------
    name : str
        Name of this object.
    machine : int
        Used machine.
    compiler : int
        Used compiler.
    data : float
        Benchmark result.
    frequency : float
        Used CPU frequency.
    cores : int
        Used number of CPU cores.
    db_id : int
        ID of associated benchmark result database table record.
    """
    name = attr.ib(type=str)
    machine = attr.ib(type=int)
    compiler = attr.ib(type=int)
    data = attr.ib(type=float)
    frequency = attr.ib(type=float)
    cores = attr.ib(type=int)
    db_id = attr.ib(type=int, init=False)

    # Database information.
    db_table = Table('benchmark_result', METADATA,
                     Column('db_id', Integer, primary_key=True),
                     Column('name', String),
                     Column('machine', Integer, ForeignKey('machine.db_id')),
                     Column('compiler', Integer, ForeignKey('compiler.db_id')),
                     Column('data', Float),
                     Column('frequency', Float),
                     Column('cores', Integer),
                     Column('upadtedIn', String, default=__version__),
                     Column('updatedOn', DateTime, default=datetime.now, onupdate=datetime.now),
                     Column('updatedBy', String, default=getuser(), onupdate=getuser()),
                     sqlite_autoincrement=True)

    def to_database(self, db_session: Session):
        """Push this benchmark record object to the database.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        -
        """
        insert(db_session, self)

    @staticmethod
    def contains(db_session: Session, machine: Machine, benchmark: 'OmpBarrierBenchmark') -> bool:
        """
        Check if the Benchmark table already contains data of a particular benchmark for a given configuration of
        machine and compiler.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.
        machine : Machine
            Used Machine.
        benchmark : Benchmark
            Used Benchmark.

        Returns:
        --------
        boolean
            True if table contains fitting data False else.
        """
        data = db_session.query(BenchmarkRecord).filter(
            BenchmarkRecord.name.like(benchmark.name), BenchmarkRecord.machine.is_(machine.db_id),
            BenchmarkRecord.compiler.is_(machine.compiler.db_id)).first()
        return bool(data)

    @staticmethod
    def update(db_session: Session, machine: Machine, benchmark: 'OmpBarrierBenchmark',
               benchmark_records: List['BenchmarkRecord']):
        """Update data records in Benchmark table with new benchmark data.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.
        machine : Machine
            Used Machine.
        benchmark : Benchmark
            Used Benchmark.
        benchmark_records : list of BenchmarkRecord
            Results of the executed benchmark as list of mathematical expression strings.

        Returns:
        --------
        -
        """
        for record in benchmark_records:
            # Select and remove all already included data records.
            try:
                queried_record = db_session.query(BenchmarkRecord).filter(
                    BenchmarkRecord.name.like(benchmark.name),
                    BenchmarkRecord.machine.is_(machine.db_id),
                    BenchmarkRecord.compiler.is_(machine.compiler.db_id),
                    BenchmarkRecord.frequency.is_(record.frequency),
                    BenchmarkRecord.cores.is_(record.cores)).one()
                # Update data record.
                queried_record.data = (record.data + queried_record.data) / 2
            except NoResultFound:
                # Insert new data record.
                record.to_database(db_session)
            except MultipleResultsFound:
                raise RuntimeError('Unable to update benchmark record!')

    @staticmethod
    def select(db_session: Session, machine: Machine, benchmarks: List[str]) -> List[KernelRecord]:
        """Retrieve BenchmarkRecord table data record(s) from the database.

        Return all records that match the provided configuration of machine, compiler and benchmark name(s).

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        machine: Machine
            Used Machine.
        benchmarks : Benchmark
            Names of used benchmarks.

        Returns:
        --------
        pandas.DataFrame
            Retrieved list of data records.
        """
        query = db_session.query(BenchmarkRecord.cores, BenchmarkRecord.name, BenchmarkRecord.data,
                                 BenchmarkRecord.frequency).filter(BenchmarkRecord.machine.is_(machine.db_id),
                                                                   BenchmarkRecord.compiler.is_(machine.compiler.db_id),
                                                                   BenchmarkRecord.name.in_(benchmarks)).order_by(
            BenchmarkRecord.cores, BenchmarkRecord.name)
        data = read_sql_query(query.statement, db_session.bind, index_col='cores')
        return data
