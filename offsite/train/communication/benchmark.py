"""@package train.communication.benchmark
Definition of classes Benchmark (abstract), BenchmarkRecord and BenchmarkType.

@author: Johannes Seiferth
"""

from abc import ABC, abstractmethod
from datetime import datetime
from enum import Enum
from getpass import getuser
from subprocess import run, CalledProcessError
from typing import Any, Dict, List

import attr
from pandas import DataFrame
from pathlib2 import Path
from sqlalchemy import Column, DateTime, Float, ForeignKey, Integer, String, Table
from sqlalchemy.orm import Session

from offsite import __version__
from offsite.database import METADATA, insert
from offsite.descriptions.machine import MachineState


class BenchmarkType(Enum):
    """Defines what type of benchmark is executed.

    - OMP_BARRIER
        OpenMP barrier benchmark.
    """
    OMP_BARRIER = 'omp_barrier'


@attr.s(hash=True)
class Benchmark(ABC):
    """Representation of the abstract benchmark base class from which all actual benchmark classes derive from.

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
        Machine this benchmark was compiled for last.
    folder: Path
        Folder to store source code in.
    """
    type = attr.ib(type=BenchmarkType, init=False)
    code = attr.ib(type=str, init=False)
    folder = attr.ib(type=Path, default=Path('tmp'), init=False)
    source = attr.ib(type=Path, init=False)
    compiled_for = attr.ib(type=MachineState, default=None, init=False)
    binary = attr.ib(type=Path, default=Path(), init=False)

    def __attrs_post_init__(self):
        # Write source code file ...
        # ... create folder if it does not yet exist.
        if self.folder and not self.folder.exists():
            self.folder.mkdir(parents=True)
        # ... and write the file.
        path = Path('{}.{}'.format(self.type.value, 'c'))
        path = self.folder / path
        with path.open('w') as file_handle:
            file_handle.write(self.code)
            self.source = Path(file_handle.name)

    def __del__(self):
        """Remove all files created during benchmarking.

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

    def compile(self, machine: MachineState):
        """Compile the benchmark's source file on a given machine.

        Parameters:
        -----------
        machine: MachineState
            Machine the benchmark is compiled for.

        Returns:
        --------
        -
        """
        self.binary = self.folder / Path('{}.out'.format(self.type.value))
        self.compiled_for = machine
        try:
            # Construct compiler call.
            args = [machine.compiler.name, str(self.source), '-fopenmp', '-o{}'.format(self.binary)]
            # Add compiler flags to call.
            args.extend((flag for flag in machine.compiler.flags.split(' ')))
            # Compile benchmark.
            run(args, check=True)
        except CalledProcessError as error:
            raise RuntimeError('Failed to compile {}: {}'.format(self.type.value, error))

    def binary_exists(self) -> bool:
        """Check if the benchmark's binary file exists.

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
        """Write source code of the benchmark to file.

        Parameters:
        -----------
        -

        Returns:
        --------
        pathlib.Path
            Relative path to source file created.
        """

    def run(self, machine: MachineState, run_config: Dict[str, Any], save_in_db: bool = True,
            write_to_file: bool = False) -> List['BenchRecord']:
        """Run the benchmark on a given machine.

        Parameters:
        -----------
        machine: MachineState
            Machine the benchmark is run on.
        run_config: dict (key=str, value=Any)
            Run parameters used to execute this benchmark.
        save_in_db: bool
            If True create valid BenchRecord database entries.

        Returns:
        --------
        list of BenchRecord
            Benchmark data obtained on the given machine.
        """
        # Compile benchmark for the given machine.
        if self.compiled_for is not machine or not self.binary_exists():
            self.compile(machine)
        # Run benchmark.
        data: List['BenchRecord'] = self.run_(machine, run_config, save_in_db)
        # ... and optionally write the results to a file.
        if write_to_file:
            self._write2file(data)
        return data

    @abstractmethod
    def run_(self, machine: MachineState, run_config: Dict[str, Any], save_in_db: bool) -> List['BenchRecord']:
        """Specific function to run a particular benchmark on a given machine.

        Has to be implemented by each benchmark.

        Parameters:
        -----------
        machine: MachineState
            Machine the benchmark is ran on.
        run_config: dict (key=str, value=Any)
            Run parameters used to execute this benchmark.
        save_in_db: bool
            If True create valid BenchRecord database entries.

        Returns:
        --------
        list of BenchRecord
            Benchmark data obtained on the given machine.
        """
        pass

    @abstractmethod
    def contains(self, db_session: Session, machine: MachineState) -> bool:
        """
        Check if the Benchmark table already contains data of this benchmark for a given configuration of machine state
        and compiler.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        machine: MachineState
            Used machine.

        Returns:
        --------
        boolean
            True if table contains fitting data False else.
        """
        pass

    @abstractmethod
    def update(self, db_session: Session, machine: MachineState, benchmark_records: List['BenchmarkRecord']):
        """Update data records in benchmark table with new benchmark data.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        machine: MachineState
            Used machine.
        benchmark_records: list of BenchmarkRecord
            Results of the executed benchmark as list of mathematical expression strings.

        Returns:
        --------
        -
        """
        pass

    @abstractmethod
    def select(self, db_session: Session, machine: MachineState) -> DataFrame:
        """Retrieve BenchmarkRecord table data record(s) from the database.

        Return all records that match the provided configuration of machine state, compiler and benchmark name(s).

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        machine: MachineState
            Used machine.

        Returns:
        --------
        pandas.DataFrame
            Retrieved list of data records.
        """
        pass

    @abstractmethod
    def _write2file(self, data) -> Path:
        """Write benchmark results to file.

        Parameters:
        -----------
        data:
            Obtained benchmark results.

        Returns:
        Path
            Path of the file created.
        """
        pass


@attr.s
class BenchmarkRecord(ABC):
    """Representation of a benchmark table database record.

    Attributes:
    -----------
    name: str
        Name of this object.
    machine: int
        Used machine.
    compiler: int
        Used compiler.
    data: float
        Benchmark result.
    run_config: dict (key=str, value=Any)
        Used run parameters.
    db_id: int
        ID of associated benchmark result database table record.
    """
    name = attr.ib(type=str)
    machine = attr.ib(type=int)
    compiler = attr.ib(type=int)
    data = attr.ib(type=float)
    run_config = attr.ib(type=Dict[str, Any])
    run_config_serial = attr.ib(type=str, init=False)
    db_id = attr.ib(type=int, init=False)

    # Database information.
    db_table = Table('benchmark_result', METADATA,
                     Column('db_id', Integer, primary_key=True),
                     Column('name', String),
                     Column('machine', Integer, ForeignKey('machine.db_id')),
                     Column('compiler', Integer, ForeignKey('compiler.db_id')),
                     Column('data', Float),
                     Column('run_config_serial', String),
                     Column('updatedIn', String, default=__version__),
                     Column('updatedOn', DateTime, default=datetime.now, onupdate=datetime.now),
                     Column('updatedBy', String, default=getuser(), onupdate=getuser()),
                     sqlite_autoincrement=True)

    def to_database(self, db_session: Session):
        """Push this benchmark record object to the database.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        -
        """
        insert(db_session, self)

    @abstractmethod
    def from_run_config_dict(self):
        pass

    @abstractmethod
    def to_run_config_dict(self):
        pass
