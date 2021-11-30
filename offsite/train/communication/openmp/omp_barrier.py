"""@package train.communication.openmp.omp_barrier
Definition of class OmpBarrierBenchmark.
"""

from subprocess import run, PIPE, CalledProcessError
from typing import Any, Dict, List, Tuple

import attr
from pandas import read_sql_query, DataFrame
from sqlalchemy.orm import Query, Session
from sqlalchemy.orm.exc import MultipleResultsFound, NoResultFound

import offsite.config
from offsite.config import Config, BenchType
from offsite.descriptions.machine import MachineState
from offsite.descriptions.parser import deserialize_obj, serialize_obj
from offsite.train.communication.benchmark import Benchmark, BenchmarkRecord
from offsite.util.math_utils import remove_outliers


@attr.s(hash=True)
class OmpBarrierBenchmark(Benchmark):
    """Representation of the OmpBarrier benchmark.

    Attributes:
    -----------
    type: BenchType
        Type of this benchmark.
    source: Path
        Path to source code file associated with this benchmark.
    binary: Path
        Path to binary object file of this benchmark compiled for a given
        configuration of machine and compiler.
    compiled_for: MachineState
        MachineState this benchmark was compiled for last.
    folder: Path
        Folder to store source code in.

    """
    type = attr.ib(type=BenchType, default=BenchType.OMP_BARRIER, init=False)

    def run_(self, machine: MachineState, run_config: Dict[str, Any], save_in_db: bool) -> List[
        BenchmarkRecord]:
        """Run the OmpBarrier benchmark on a given machine.

        Parameters:
        -----------
        machine: MachineState
            MachineState the benchmark is ran on.
        run_config: dict (key=str, value=Any)
            Run parameters used to executed this benchmark.
        save_in_db: bool
            If True create valid BenchRecord database entries.

        Returns:
        --------
        list of BenchRecord
            Benchmark data obtained on the given machine.
        """
        config: Config = offsite.config.offsiteConfig
        # Read run parameters.
        try:
            repetitions: int = run_config['reps']
        except KeyError:
            raise RuntimeError('')
        # Compile benchmark for given machine.
        if self.compiled_for is not machine or not self.binary_exists():
            self.compile(machine)
        # Pin CPU frequency for benchmarking.
        if config.args.verbose:
            print('\nPinning CPU frequency to {} Hz for benchmark \'{}\'.'.format(machine.clock, self.type.value))
        MachineState.pin_cpu_frequency(machine.clock)
        # Run benchmark for all core numbers.
        data: List[OmpBarrierRecord] = list()
        for cores in range(1, machine.coresPerSocket + 1):
            run_config = {'freq': machine.clock, 'cores': cores}
            times: List[float] = list()
            try:
                for _ in range(1, repetitions + 1):
                    cmd = ['./{}'.format(self.binary), str(cores), str(repetitions)]
                    # Run benchmark for current configuration.
                    runtime: str = run(cmd, check=True, encoding='utf-8', stdout=PIPE).stdout
                    times.append(float(runtime))
            except CalledProcessError as error:
                print('OmpBarrierBenchmark failed: {}'.format(error))
            # Remove outliers.
            times: List[float] = remove_outliers(times)
            # Runtime as mean of all times measured.
            runtime: float = sum(times) / len(times)
            # Save benchmark record.
            record = OmpBarrierRecord(self.type.value, machine.db_id if save_in_db else -1,
                                      machine.compiler.db_id if save_in_db else -1, runtime, run_config)
            data.append(record)
        # Reset CPU frequency limits after benchmarking.
        if config.args.verbose:
            print('\nResetting CPU frequency limits after benchmark \'{}\'.'.format(self.type.value))
        MachineState.reset_cpu_frequency_limits()
        return data

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

    def contains(self, db_session: Session, machine: MachineState) -> bool:
        """
        Check if the Benchmark table already contains data of a particular benchmark for a given configuration of
        machine and compiler.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        machine: MachineState
            Used MachineState.

        Returns:
        --------
        boolean
            True if table contains fitting data False else.
        """
        return bool(db_session.query(OmpBarrierRecord).filter(
            OmpBarrierRecord.name.like(self.type.value), OmpBarrierRecord.machine.is_(machine.db_id),
            OmpBarrierRecord.compiler.is_(machine.compiler.db_id)).first())

    def update(self, db_session: Session, machine: MachineState, benchmark_records: List['BenchmarkRecord']):
        """Update data records in benchmark table with new benchmark data.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        machine: MachineState
            Used MachineState.
        benchmark_records: list of BenchmarkRecord
            Results of the executed benchmark as list of mathematical expression strings.

        Returns:
        --------
        -
        """
        for record in benchmark_records:
            # Select and remove all already included data records.
            try:
                queried_record: OmpBarrierRecord = db_session.query(OmpBarrierRecord).filter(
                    OmpBarrierRecord.name.like(self.type.value),
                    OmpBarrierRecord.machine.is_(machine.db_id),
                    OmpBarrierRecord.compiler.is_(machine.compiler.db_id),
                    OmpBarrierRecord.run_config_serial.is_(serialize_obj(record.run_config))).one()
                # Update data record.
                queried_record.data = (record.data + queried_record.data) / 2
            except NoResultFound:
                # Insert new data record.
                record.to_database(db_session)
            except MultipleResultsFound:
                raise RuntimeError('Unable to update benchmark record!')

    def select(self, db_session: Session, machine: MachineState) -> DataFrame:
        """Retrieve BenchmarkRecord table data record(s) from the database.

        Return all records that match the provided configuration of machine and compiler.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        machine: MachineState
            Used MachineState.
        benchmarks: Benchmark
            Names of used benchmarks.

        Returns:
        --------
        pandas.DataFrame
            Retrieved list of data records.
        """
        query: Query = db_session.query(
            OmpBarrierRecord.name, OmpBarrierRecord.data, OmpBarrierRecord.run_config_serial).filter(
            OmpBarrierRecord.machine.is_(machine.db_id), OmpBarrierRecord.compiler.is_(machine.compiler.db_id),
            OmpBarrierRecord.name.is_(self.type.value))
        data: DataFrame = read_sql_query(query.statement, db_session.bind)
        #
        data['cores'] = ''
        data['freq'] = ''
        for row in data.itertuples():
            rcfg = deserialize_obj(row.run_config_serial)
            data.at[row.Index, 'freq'] = rcfg['freq']
            data.at[row.Index, 'cores'] = rcfg['cores']
        data.drop('run_config_serial', axis=1, inplace=True)
        data.sort_values(by=['cores'], inplace=True)
        return data


@attr.s
class OmpBarrierRecord(BenchmarkRecord):
    """Representation of a OpenMP barrier benchmark table database record.

    Attributes:
    -----------
    frequency: float
        Used CPU frequency.
    cores: int
        Used number of CPU cores.
    db_id: int
        ID of associated benchmark result database table record.
    """
    # Benchmark specific runtime parameters derived from attribute run_config.
    frequency = attr.ib(type=float, init=False)
    cores = attr.ib(type=int, init=False)

    def __attrs_post_init__(self):
        self.run_config_serial = serialize_obj(self.run_config)
        self.from_run_config_dict()

    def from_run_config_dict(self):
        try:
            self.frequency = self.run_config['freq']
            self.cores = self.run_config['cores']
        except KeyError:
            print('TODO')

    def to_run_config_dict(self):
        self.run_config['freq'] = self.frequency
        self.run_config['cores'] = self.cores
        self.run_config_serial = serialize_obj(self.run_config)

    @staticmethod
    def read_bench_file(yaml: Dict, machine: MachineState) -> Tuple[str, List[BenchmarkRecord]]:
        # Convert data to Benchmark records.
        bench_name = yaml['benchmark']
        bench_data = [OmpBarrierRecord(bench_name, machine.db_id, machine.compiler.db_id, row[1],
                                       {'freq': machine.clock, 'cores': row[0]}) for row in yaml['data']]
        return bench_name, bench_data

    @staticmethod
    def check_bench_file(data, machine: MachineState) -> bool:
        errors: List[Tuple[str, str]] = list()
        # Check machine name.
        if data['machine'] != machine.name:
            errors.append(('machine \'{}\'!'.format(data['machine']), '\'{}\'!'.format(machine.name)))
        # Check compiler.
        name, version = data['compiler'].split(' ')
        if name != machine.compiler.name:
            errors.append(('compiler \'{}\'!'.format(name), '\'{}\'!'.format(machine.compiler.name)))
        if version != machine.compiler.version:
            errors.append(('compiler version \'{}\'!'.format(version), '\'{}\'!'.format(machine.compiler.version)))
        flags = data['flags']
        if flags != machine.compiler.flags:
            errors.append(('compiler flags \'{}\'!'.format(flags), '\'{}\'!'.format(machine.compiler.flags)))
        # Check CPU frequency.
        if data['frequency'] != machine.clock:
            errors.append(('CPU frequency \'{} Hz\'!'.format(data['frequency']), '\'{} Hz\'!'.format(machine.clock)))
        # Print all errors.
        is_valid: bool = True
        if len(errors) > 0:
            is_valid = False
            print('\nRerun required! Passed omp_barrier benchmark data are invalid:')
            for raised, required in errors:
                print('\t* raised for {} but requires {}.'.format(raised, required))
        return is_valid
