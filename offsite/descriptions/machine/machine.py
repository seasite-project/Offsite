"""@package descriptions.machine.machine
Definition of class MachineState, Compiler, Machine, NetworkConfiguration.
"""

from ctypes import sizeof, c_double
from datetime import datetime
from getpass import getuser
from pathlib import Path
from subprocess import run, CalledProcessError

import attr
# noinspection PyUnresolvedReferences
from kerncraft.prefixedunit import PrefixedUnit
from sqlalchemy import Column, DateTime, Float, ForeignKey, Integer, String, Table, UniqueConstraint
from sqlalchemy.orm import Session
from sqlalchemy.orm.exc import NoResultFound, MultipleResultsFound

import offsite.config
from offsite import __version__
from offsite.config import Config
from offsite.database import METADATA, insert
from offsite.descriptions.machine.compiler import Compiler
from offsite.descriptions.parser import load_yaml


@attr.s
class MachineState:
    """Representation of a MachineState object.

    A MachineState object describes the characteristics of the hardware and machine architecture of a single hardware
    platform (e.g. cluster node).

    Attributes:
    -----------
    name: str
        Name of this object.
    db_id: int
        ID of associated MachineState database table record.
    """
    name = attr.ib(type=str)
    path = attr.ib(type=Path)
    microArchitecture = attr.ib(type=str)
    modelType = attr.ib(type=str)
    modelName = attr.ib(type=str)
    clock = attr.ib(type=float)
    coresPerSocket = attr.ib(type=int)
    coresPerNumaDomain = attr.ib(type=int)
    numaDomainsPerSocket = attr.ib(type=int)
    cachelineSize = attr.ib(type=float)
    elements_per_cacheline = attr.ib(type=float)
    l1Cache = attr.ib(type=float)
    l1CacheElements = attr.ib(type=float)
    l2Cache = attr.ib(type=float)
    l2CacheElements = attr.ib(type=float)
    l3Cache = attr.ib(type=float)
    l3CacheElements = attr.ib(type=float)
    flopsDpAdd = attr.ib(type=int)
    flopsDpFma = attr.ib(type=int)
    flopsDpMul = attr.ib(type=int)
    flopsSpAdd = attr.ib(type=int)
    flopsSpFma = attr.ib(type=int)
    flopsSpMul = attr.ib(type=int)
    compiler = attr.ib(type=Compiler)
    db_id = attr.ib(type=int, init=False)

    # Database information.
    db_table = Table('machine_state', METADATA,
                     Column('db_id', Integer, primary_key=True),
                     Column('name', String, unique=True),
                     Column('microArchitecture', String),
                     Column('modelType', String),
                     Column('modelName', String),
                     Column('clock', Float),
                     Column('coresPerSocket', Integer),
                     Column('coresPerNumaDomain', Integer),
                     Column('numaDomainsPerSocket', Integer),
                     Column('cachelineSize', Float),
                     Column('l1Cache', Float),
                     Column('l2Cache', Float),
                     Column('l3Cache', Float),
                     Column('flopsDpAdd', Integer),
                     Column('flopsDpFma', Integer),
                     Column('flopsDpMul', Integer),
                     Column('flopsSpAdd', Integer),
                     Column('flopsSpFma', Integer),
                     Column('flopsSpMul', Integer),
                     Column('updatedIn', String, default=__version__),
                     Column('updatedOn', DateTime, default=datetime.now, onupdate=datetime.now),
                     Column('updatedBy', String, default=getuser(), onupdate=getuser()),
                     sqlite_autoincrement=True)

    @classmethod
    def from_yaml(cls, yaml_path: Path, used_compiler: str) -> 'MachineState':
        """Construct MachineState object from YAML definition.

        Parameters:
        -----------
        yaml_path: Path
            Relative path to this object's YAML file.
        used_compiler: str
            Name of used compiler.
        Returns:
        --------
        MachineState
            Created MachineState object.
        """
        config: Config = offsite.config.offsiteConfig
        # Load YAML data.
        path = Path(yaml_path)
        yaml = load_yaml(path)
        # Attribute name.
        name = path.stem
        # Attribute micro_architecture.
        if 'IACA' in yaml['in-core model']:
            micro_architecture = yaml['in-core model']['IACA']
        elif 'OSACA' in yaml['in-core model']:
            micro_architecture = yaml['in-core model']['OSACA']
        elif 'LLVM-MCA' in yaml['in-core model']:
            micro_architecture = yaml['in-core model']['LLVM-MCA']
        else:
            assert False
        # Attribute incore model.
        incore_model = yaml['in-core model']
        # .. check if this machine supports the selected incore model tool.
        if config.pred_incore_tool.value not in incore_model:
            raise RuntimeError(
                'MachineState \'{}\' does not support incore model \'{}\'!'.format(path, config.pred_incore_tool.value))
        # Attribute model_type.
        model_type = yaml['model type']
        # Attribute model_name.
        model_name = yaml['model name']
        # Attribute clock.
        clock = float(yaml['clock'])
        # Attribute cores_per_socket.
        cores_per_socket = yaml['cores per socket']
        # Attribute cores_per_numa_domain.
        cores_per_numa_domain = yaml['cores per NUMA domain']
        # Attribute numa_domains_per_socket.
        numa_domains_per_socket = yaml['NUMA domains per socket']
        # Attribute cacheline_size [in bytes].
        cacheline_size = yaml['cacheline size'].base_value()
        # Elements per cache line.
        elements_per_cacheline = cacheline_size / sizeof(c_double)
        # Cache information...
        l1_cache = -1
        l2_cache = -1
        l3_cache = -1
        l1_elements = -1
        l2_elements = -1
        l3_elements = -1
        for level in yaml['memory hierarchy']:
            level_name = level['level']
            if level_name == 'MEM':
                continue
            cpg = level['cache per group']
            size = cpg['sets'] * cpg['ways'] * cpg['cl_size']
            if level_name == 'L1':
                # Attribute l1_cache.
                l1_cache = size
                # Attribute l1_elements.
                l1_elements = size / elements_per_cacheline
            elif level_name == 'L2':
                # Attribute l2_cache.
                l2_cache = size
                # Attribute l2_elements.
                l2_elements = size / elements_per_cacheline
            elif level_name == 'L3':
                # Attribute l3_cache.
                l3_cache = size
                # Attribute l3_elements.
                l3_elements = size / elements_per_cacheline
            else:
                raise RuntimeError('Unsupported memory hierarchy level {}!'.format(name))
        # FlOPs per cycle information...
        flops_dp_add = yaml['FLOPs per cycle']['DP']['ADD']
        flops_sp_add = yaml['FLOPs per cycle']['SP']['ADD']
        flops_dp_mul = yaml['FLOPs per cycle']['DP']['MUL']
        flops_sp_mul = yaml['FLOPs per cycle']['SP']['MUL']
        try:
            flops_dp_fma = yaml['FLOPs per cycle']['DP']['FMA']
        except KeyError:
            flops_dp_fma = 0
        try:
            flops_sp_fma = yaml['FLOPs per cycle']['SP']['FMA']
        except KeyError:
            flops_sp_fma = 0
        # Attribute compiler.
        compiler = Compiler.from_yaml(yaml['compiler'], used_compiler)
        # Create object.
        return cls(name, path, micro_architecture, model_type, model_name, clock, cores_per_socket,
                   cores_per_numa_domain, numa_domains_per_socket, cacheline_size, elements_per_cacheline, l1_cache,
                   l1_elements, l2_cache, l2_elements, l3_cache, l3_elements, flops_dp_add, flops_dp_fma, flops_dp_mul,
                   flops_sp_add, flops_sp_fma, flops_sp_mul, compiler)

    @classmethod
    def from_database(cls, db_session: Session, machine_id: int, compiler_id: int) -> 'MachineState':
        """Construct MachineState object from database record.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        machine_id: int
            Database ID of the requested MachineState object.
        compiler_id: int
            Database ID of the associated Compiler object.

        Returns:
        --------
        MachineState
            Created MachineState object.
        """
        try:
            machine: MachineState = db_session.query(MachineState).filter(
                MachineState.db_id.is_(machine_id)).one()
        except NoResultFound:
            raise RuntimeError('Unable to load MachineState object from database!')
        except MultipleResultsFound:
            raise RuntimeError('Unable to load MachineState object from database!')
        # Elements per cache line.
        machine.elements_per_cacheline = machine.cachelineSize / sizeof(c_double)
        # Cache information...
        machine.l1CacheElements = machine.l1Cache / machine.elements_per_cacheline
        machine.l2CacheElements = machine.l2Cache / machine.elements_per_cacheline
        machine.l3CacheElements = machine.l3Cache / machine.elements_per_cacheline
        # Attribute compiler.
        machine.compiler = Compiler.from_database(db_session, compiler_id)
        return machine

    def to_database(self, db_session: Session) -> 'MachineState':
        """Push this MachineState object to the database.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        MachineState
            Instance of this object connected to database session.
        """
        # Check if database already contains the machine object.
        machine: MachineState = db_session.query(MachineState).filter(
            MachineState.name.like(self.name)).one_or_none()
        if machine:
            # Supplement attributes not saved in database.
            machine.path = self.path
            # machine.sockets = self.sockets
            machine.elements_per_cacheline = self.elements_per_cacheline
            machine.l1CacheElements = self.l1CacheElements
            machine.l2CacheElements = self.l2CacheElements
            machine.l3CacheElements = self.l3CacheElements
            machine.compiler = self.compiler.to_database(db_session)
            return machine
        # Add new object to database.
        insert(db_session, self)
        self.compiler = self.compiler.to_database(db_session)
        return self

    @staticmethod
    def select(db_session: Session, machine_id: int, compiler_id: int) -> 'MachineState':
        """
        Retrieve the MachineState table data record and Compiler table record from the database that match the given
        database IDs.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        machine_id: int
            ID of the requested MachineState object.
        compiler_id: int
            ID of the requested Compiler object.

        Returns:
        --------
        MachineState
            MachineState object.
        """
        return MachineState.from_database(db_session, machine_id, compiler_id)

    @staticmethod
    def set_cpu_frequency_limits(lower_limit: float, upper_limit: float):
        config: Config = offsite.config.offsiteConfig
        assert lower_limit <= upper_limit
        if lower_limit == upper_limit:
            return MachineState.pin_cpu_frequency(lower_limit)
        try:
            # Set minimal CPU frequency.
            cmd = [config.likwid_set_frequencies, '-x', str(lower_limit)]
            run(cmd, check=True)
            # Set maximal CPU frequency.
            cmd = [config.likwid_set_frequencies, '-y', str(upper_limit)]
            run(cmd, check=True)
        except CalledProcessError as err:
            raise RuntimeWarning('Failed to limit CPU frequency to range \'{}\' to \'{}\': {}'.format(
                lower_limit, upper_limit, err))

    @staticmethod
    def reset_cpu_frequency_limits():
        config: Config = offsite.config.offsiteConfig
        try:
            cmd = [config.likwid_set_frequencies, '-reset']
            run(cmd, check=True)
        except CalledProcessError as err:
            raise RuntimeWarning('Failed to reset CPU frequency to its default range: {}'.format(err))

    @staticmethod
    def pin_cpu_frequency(frequency: float):
        config: Config = offsite.config.offsiteConfig
        try:
            cmd = [config.likwid_set_frequencies, '-f', str(frequency)]
            run(cmd, check=True)
        except CalledProcessError as err:
            raise RuntimeWarning('Failed to pin CPU frequency to \'{}\': {}'.format(frequency, err))


def parse_machine_state(machine_file: Path, compiler: str) -> MachineState:
    """Parse machine description YAML file and return MachineState object.

    Parameters:
    -----------
    machine_file: pathlib.Path
        Relative path to machine's description file.
    compiler: str
        Name of the compiler used.

    Returns:
    --------
    MachineState
        MachineState object created.
    """
    return MachineState.from_yaml(machine_file, compiler)


@attr.s
class Machine:
    """
        Representation of a machine table database record.
    """
    name = attr.ib(type=str)
    machine_cfg = attr.ib(type=int)  # group_id = attr.ib(type=int, default=None)
    network_cfg = attr.ib(type=int)

    # Database information.
    db_table = Table('machine', METADATA,
                     Column('db_id', Integer, primary_key=True),
                     Column('name', String),
                     Column('machine_cfg', Integer, ForeignKey('machine_state.db_id')),
                     Column('network_cfg', Integer, ForeignKey('network_config.db_id')),
                     Column('updatedOn', DateTime, default=datetime.now, onupdate=datetime.now),
                     Column('updatedBy', String, default=getuser(), onupdate=getuser()),
                     UniqueConstraint('name', 'machine_cfg', 'network_cfg'), sqlite_autoincrement=True)
