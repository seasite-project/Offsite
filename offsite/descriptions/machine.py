"""@package machine
Definition of class Machine.
"""

from ctypes import sizeof, c_double
from datetime import datetime
from getpass import getuser
from pathlib import Path
from subprocess import PIPE, run
from sys import version_info
from typing import Dict

import attr
from kerncraft.prefixedunit import PrefixedUnit
from sqlalchemy import Column, DateTime, Float, Integer, String, Table, \
    UniqueConstraint
from sqlalchemy.orm import Session

from offsite import __version__
from offsite.db import METADATA
from offsite.db.db import insert
from offsite.descriptions.parser_utils import load_yaml


@attr.s
class Compiler:
    """Representation of a Compiler object.

    A Compiler object describes a compiler configuration.

    Attributes:
    -----------
    name : str
        Name of the compiler.
    version : str
        Used version of the compiler.
    flags : str
        Used compiler flags.
    db_id : int
        ID of associated Compiler database table record.
    """
    name = attr.ib(type=str)
    version = attr.ib(type=str)
    flags = attr.ib(type=str)
    db_id = attr.ib(type=int, init=False)

    # Database information.
    db_table = Table('compiler', METADATA,
                     Column('db_id', Integer, primary_key=True),
                     Column('name', String),
                     Column('version', String),
                     Column('flags', String),
                     Column('createdBy', String, default=getuser()),
                     Column('createdIn', String, default=__version__),
                     Column('createdOn', DateTime, default=datetime.now),
                     Column('updatedOn', DateTime, default=datetime.now, onupdate=datetime.now),
                     Column('updatedBy', String, default=getuser(), onupdate=getuser()),
                     UniqueConstraint('name', 'version', 'flags'),
                     sqlite_autoincrement=True)

    @classmethod
    def from_yaml(cls, yaml: Dict[str, str], used_compiler: str) -> 'Compiler':
        """Construct Compiler object from YAML definition.

        Parameters:
        -----------
        yaml : dict
            YAML object describing this object.
        used_compiler : str
            Name of used compiler.
        Returns:
        --------
        Compiler
            Created Compiler object.
        """
        try:
            item = next(filter(lambda x: x[0] == used_compiler, yaml.items()))
            # Attribute compiler_name.
            name = item[0]
            # Attribute compiler_version.
            version = Compiler.determine_compiler_version(name)
            # Attribute compiler_flags.
            flags = item[1]
        except StopIteration:
            raise RuntimeError('Compiler {} not found in YAML description!'.format(used_compiler))
        return cls(name, version, flags)

    def to_database(self, db_session: Session):
        """Push this Machine object to the database.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        Compiler
            Instance of this object connected to database session.
        """
        # Check if database already contains the compiler object.
        compiler = db_session.query(Compiler).filter(Compiler.name.like(self.name), Compiler.version.like(self.version),
                                                     Compiler.flags.like(self.flags)).one_or_none()
        if compiler:
            return compiler
        # Add new object to database.
        insert(db_session, self)
        return self

    @staticmethod
    def determine_compiler_version(compiler: str) -> str:
        """Determine the version of the compiler used and return the version string.

        Parameters:
        -----------
        compiler : str
            Used compiler command.

        Returns:
        --------
        str
            Version string of the used compiler.
        """
        if version_info[1] > 5:
            version = run([compiler, '-dumpversion'], check=True, encoding='utf-8', stdout=PIPE).stdout
        else:
            version = run([compiler, '-dumpversion'], check=True, stdout=PIPE).stdout.decode("utf-8")
        # Strip trailing EOL characters.
        return version.rstrip()


@attr.s
class Machine:
    """Representation of a Machine object.

    A Machine object describes the characteristics of the hardware and machine architecture of a single hardware
    platform / cluster node.

    Attributes:
    -----------
    name : str
        Name of this object.
    db_id : int
        ID of associated Machine database table record.
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
    sockets = attr.ib(type=int)
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
    db_table = Table('machine', METADATA,
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
                     Column('createdBy', String, default=getuser()),
                     Column('createdIn', String, default=__version__),
                     Column('createdOn', DateTime, default=datetime.now),
                     Column('updatedOn', DateTime, default=datetime.now, onupdate=datetime.now),
                     Column('updatedBy', String, default=getuser(), onupdate=getuser()),
                     sqlite_autoincrement=True)

    @classmethod
    def from_yaml(cls, yaml_path: str, used_compiler: str) -> 'Machine':
        """Construct Machine object from YAML definition.

        Parameters:
        -----------
        yaml_path : str
            Relative path to this object's YAML file.
        used_compiler : str
            Name of used compiler.
        Returns:
        --------
        Machine
            Created Machine object.
        """
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
        # Attribute sockets.
        sockets = yaml['sockets']
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
                raise RuntimeError(
                    'Unsupported memory hierarchy level {}!'.format(name))
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
                   cores_per_numa_domain, numa_domains_per_socket, sockets, cacheline_size, elements_per_cacheline,
                   l1_cache, l1_elements, l2_cache, l2_elements, l3_cache, l3_elements, flops_dp_add, flops_dp_fma,
                   flops_dp_mul, flops_sp_add, flops_sp_fma, flops_sp_mul, compiler)

    def to_database(self, db_session: Session):
        """Push this Machine object to the database.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        Machine
            Instance of this object connected to database session.
        """
        # Check if database already contains the machine object.
        machine = db_session.query(Machine).filter(Machine.name.like(self.name)).one_or_none()
        if machine:
            # Supplement attributes not saved in database.
            machine.path = self.path
            machine.sockets = self.sockets
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
