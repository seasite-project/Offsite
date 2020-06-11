"""@package impl_skeleton
Definition of class ImplSkeleton.
"""

from copy import deepcopy
from datetime import datetime
from getpass import getuser
from pathlib import Path
from typing import Dict, List, Tuple

import attr
from lark import Lark
from sqlalchemy import Column, DateTime, Enum, ForeignKey, Integer, String, Table, UniqueConstraint
from sqlalchemy.orm import Session
from sqlalchemy.orm.exc import NoResultFound, MultipleResultsFound

from offsite import __version__
from offsite.codegen.code_tree import CodeTree, CodeTreeGenerator
from offsite.config import ModelToolType
from offsite.db import METADATA
from offsite.db.db import insert
from offsite.descriptions.ivp import IVP
from offsite.descriptions.kernel_template import KernelTemplate
from offsite.descriptions.ode_method import ODEMethod
from offsite.descriptions.parser_utils import load_yaml, deserialize_obj, serialize_obj
from offsite.evaluation.math_utils import corrector_steps, ivp_grid_size, stages, eval_math_expr


@attr.s
class ImplSkeleton:
    """Representation of an ImplSkeleton object.

    Attributes:
    -----------
    name : str
        Name of this object.
    code_tree : CodeTree
        Tree representation of the kernel code.
    communicationOperations : dict
        Communication operations needed when executing a single iteration step.
    loops : dict
        Loops used in the implementation code.
    kernels : dict
        Kernels used in the implementation code.
    connected_templates : list
        List of KernelTemplate objects associated with this object.
    db_id : int
        ID of associated ImplSkeleton database table record.
    """
    name = attr.ib(type=str)
    code_tree = attr.ib(type=CodeTree, hash=False)
    communicationOperations = attr.ib(type=Dict)
    loops = attr.ib(type=Dict)
    kernels = attr.ib(type=Dict)
    connected_templates = attr.ib(type=List['KernelTemplate'])
    codegen = attr.ib(type=Dict, default=list())
    modelTool = attr.ib(type=str, default=ModelToolType.KERNCRAFT)
    isIVPdependent = attr.ib(type=bool, default=False)
    communication_operations_serial = attr.ib(type=str, init=False)
    loops_serial = attr.ib(type=str, init=False)
    kernels_serial = attr.ib(type=str, init=False)
    connected_templates_ids = attr.ib(type=str, init=False)
    code_tree_serial = attr.ib(type=str, init=False)
    codegen_serial = attr.ib(type=str, init=False)
    db_id = attr.ib(type=int, init=False)

    # Database information.
    db_table = Table('impl_skeleton', METADATA,
                     Column('db_id', Integer, primary_key=True),
                     Column('name', String, unique=True),
                     Column('modelTool', Enum(ModelToolType)),
                     Column('communication_operations_serial', String),
                     Column('loops_serial', String),
                     Column('kernels_serial', String),
                     Column('connected_templates_ids', String),
                     Column('code_tree_serial', String),
                     Column('codegen_serial', String),
                     Column('updatedIn', String, default=__version__),
                     Column('updatedOn', DateTime, default=datetime.now, onupdate=datetime.now),
                     Column('updatedBy', String, default=getuser(), onupdate=getuser()),
                     sqlite_autoincrement=True)

    @classmethod
    def from_yaml(cls, yaml_path: Path, kernel_templates: List[KernelTemplate]) -> 'ImplSkeleton':
        """Construct ImplSkeleton object from YAML definition.

        Parameters:
        -----------
        yaml_path : Path
            Relative path to this object's YAML file.
        kernel_templates : list of KernelTemplate
            List of available KernelTemplate objects.

        Returns:
        --------
        ImplSkeleton
            Created ImplSkeleton object.
        """
        # Load YAML data.
        yaml = load_yaml(yaml_path)
        # Attribute name.
        name = yaml_path.stem
        # Attribute code_tree.
        with open('offsite/codegen/code_dsl.lark') as f:
            grammar = f.read()
        parser = Lark(grammar, parser='lalr')
        parsed_code = parser.parse(yaml['code'])
        code_tree = CodeTreeGenerator().generate(parsed_code)
        code_tree = deepcopy(code_tree)
        # Parse code to garner attributes communication operations, loops and kernels.
        communication_operations, loops, kernels = ImplSkeleton.parse_code(yaml['code'])
        # Attribute connected_templates.
        # .. used to temporarily store list of kernel templates for post init.
        connected_templates = kernel_templates
        # Attribute codegen.
        codegen = dict()
        if 'codegen' in yaml:
            for key, value in yaml['codegen'].items():
                if key == 'loop splits':
                    codegen[key] = [v for v in value]
                else:
                    codegen[key] = value
        # Create object.
        return cls(name, code_tree, communication_operations, loops, kernels, connected_templates, codegen)

    @classmethod
    def from_database(cls, db_session: Session, impl_skeleton_name: str,
                      kernel_templates: List[KernelTemplate]) -> 'ImplSkeleton':
        """Construct ImplSkeleton object from database record.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.
        impl_skeleton_name : str
            Name of the impl skeleton, which is used as primary key in the database.
        kernel_templates : list of KernelTemplate
            List of available KernelTemplate objects.

        Returns:
        --------
        ImplSkeleton
            Created ImplSkeleton object.
        """
        try:
            skeleton = db_session.query(ImplSkeleton).filter(ImplSkeleton.name.like(impl_skeleton_name)).one()
        except NoResultFound:
            raise RuntimeError('Unable to load ImplSkeleton object from database!')
        except MultipleResultsFound:
            raise RuntimeError('Unable to load ImplSkeleton object from database!')
        # Attribute code_tree.
        skeleton.code_tree = deserialize_obj(skeleton.code_tree_serial)
        # Attribute communicationOperations.
        skeleton.communicationOperations = deserialize_obj(skeleton.communication_operations_serial)
        # Attribute loops.
        skeleton.loops = deserialize_obj(skeleton.loops_serial)
        # Attribute kernels.
        skeleton.kernels = deserialize_obj(skeleton.kernels_serial)
        # Attribute connected_templates.
        skeleton.connected_templates = list()
        for name, executions in skeleton.kernels.items():
            template = [s for s in kernel_templates if s.name == name]
            if not template:
                raise RuntimeError('Unable to find required KernelTemplate {}!'.format(name))
            if len(template) > 1:
                raise RuntimeError('Unable to unambiguously identify KernelTemplate {}!'.format(name))
            skeleton.connected_templates.append((template[0], executions))
            if template[0].isIVPdependent:
                skeleton.isIVPdependent = True
        # Attribute codegen.
        skeleton.codegen = deserialize_obj(skeleton.codegen_serial)
        return skeleton

    def __attrs_post_init__(self):
        """Create the Kernel objects associated with this object.

        Parameters:
        -----------
        -

        Returns:
        --------
        -
        """
        avail_templates = self.connected_templates
        self.connected_templates = list()
        for name, executions in self.kernels.items():
            template = [s for s in avail_templates if s.name == name]
            if not template:
                raise RuntimeError('Unable to find required KernelTemplate {}!'.format(name))
            if len(template) > 1:
                raise RuntimeError('Unable to unambiguously identify KernelTemplate {}!'.format(name))
            self.connected_templates.append((template[0], executions))
            # Adjust isIVPdependent
            if template[0].isIVPdependent:
                self.isIVPdependent = True
            # Adjust model tool to YASKSITE if one of the kernel templates associated uses it.
            if self.modelTool == ModelToolType.KERNCRAFT and template[0].modelTool != ModelToolType.KERNCRAFT:
                self.modelTool = template[0].modelTool

    def to_database(self, db_session: Session):
        """Push this impl skeleton object to the database.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        ImplSkeleton
            Instance of this object connected to database session.
        """
        # Check if database already contains this KernelTemplate object.
        skeleton = db_session.query(ImplSkeleton).filter(ImplSkeleton.name.like(self.name)).one_or_none()
        if skeleton:
            # Supplement attributes not saved in database.
            skeleton.communicationOperations = self.communicationOperations
            skeleton.loops = self.loops
            skeleton.kernels = self.kernels
            skeleton.connected_templates = self.connected_templates
            for template in skeleton.connected_templates:
                if template[0].isIVPdependent:
                    skeleton.isIVPdependent = True
                    break
            skeleton.code_tree = self.code_tree
            skeleton.codegen = self.codegen
            return skeleton
        # Add new object to database.
        # Serialize data.
        self.communication_operations_serial = serialize_obj(self.communicationOperations)
        self.loops_serial = serialize_obj(self.loops)
        self.kernels_serial = serialize_obj(self.kernels)
        self.code_tree_serial = serialize_obj(self.code_tree)
        self.codegen_serial = serialize_obj(self.codegen)
        # IDs of connected templates.
        self.connected_templates_ids = ','.join(sorted([str(connect[0].db_id) for connect in self.connected_templates]))
        # ImplSkeleton object.
        insert(db_session, self)
        return self

    @staticmethod
    def parse_code(yaml: str) -> Tuple[str, str, str]:
        """Parse code of an implementation skeleton for communication operations, loops and kernels.

        Parameters:
        -----------
        yaml : str
            YAML object describing this object.

        Returns:
        --------
        str
            Communication operations used.
        str
            Loops used.
        str
            Kernels used.
        """
        communication_operations = dict()
        loops = dict()
        kernels = dict()
        # Parse code.
        zero_expr = eval_math_expr('0')
        execution_count = eval_math_expr('1')
        loop_stack = list()
        for line in yaml.split('%'):
            line = line.lstrip()
            line = line.rstrip()
            if line.startswith('COM'):
                op_name = line.split(' ')[1]
                if op_name not in communication_operations:
                    communication_operations[op_name] = zero_expr
                communication_operations[op_name] += execution_count
            elif line.startswith('KERNEL'):
                kernel_name = line.split(' ')[1]
                if kernel_name not in kernels:
                    kernels[kernel_name] = zero_expr
                kernels[kernel_name] += execution_count
            elif line.startswith('LOOP_START'):
                loop_name = line.split(' ')[1]
                loop_iterations = line.split(' ')[2]
                loop_stack.append(loop_name)
                execution_count = execution_count * eval_math_expr(loop_iterations)
                if loop_name in loops:
                    assert False
                loops[loop_name] = eval_math_expr(execution_count)
            elif line.startswith('LOOP_END'):
                loop_name = line.split(' ')[1]
                assert loop_name == loop_stack[-1]
                loop_stack.pop()
                loop_iterations = loops[loop_name]
                execution_count = execution_count / eval_math_expr(loop_iterations)
        return communication_operations, loops, kernels

    def compute_communication_costs(
            self, bench_data: 'pandas.DataFrame', method: ODEMethod, ivp: IVP) -> Dict[int, str]:
        """
        Compute the communication costs of a single iteration step of this object for the given configuration of ODE
        method and IVP.

        Parameters:
        -----------
        bench_data: pandas.DataFrame
            Available benchmark data sorted by number of cores and benchmark name.
        method: ODEMethod
            Used ODE method.
        ivp: IVP
            Used IVP.

        Returns:
        dict(key: number of cores)
            Communication costs of a single iteration step of this impl skeleton sorted by number of cores.
        """
        constants = [corrector_steps(method), stages(method), ivp_grid_size(ivp.gridSize)]
        costs = dict()
        # Compute for all number of cores...
        for cores, row in bench_data.iterrows():
            costs[cores] = 0.0
            operation = row['name']
            frequency = row['frequency']
            # Select measured benchmark result for this communication operation and multiply it with the number of
            # repetitions per iteration step of this operation.
            if operation in self.communicationOperations:
                costs[cores] += row['data'] * eval_math_expr(
                    self.communicationOperations[operation], constants, cast_to=float)
        return costs

    @staticmethod
    def select(db_session: Session, skeleton_id: int) -> 'ImplSkeleton':
        """Retrieve the ImplSkeleton table data record from the database.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        skeleton_id : int
            ID of the impl skeleton requested.

        Returns:
        --------
        ImplSkeleton
            Retrieved data record.
        """
        skeleton = db_session.query(ImplSkeleton).filter(ImplSkeleton.db_id.is_(skeleton_id)).one()
        # Deserialize attributes...
        # Attribute communicationOperations.
        skeleton.communicationOperations = deserialize_obj(skeleton.communication_operations_serial)
        # Attribute loops.
        skeleton.loops = deserialize_obj(skeleton.loops_serial)
        # Attribute kernels.
        skeleton.kernels = deserialize_obj(skeleton.kernels_serial)
        # Attribute connected_templates.
        # skeleton.connected_templates = [KernelTemplateConnect(
        #    name, skeleton, execs, kernel_templates)
        #    for name, execs in skeleton.kernels.items()]
        # Attribute isIVPdependent.
        skeleton.isIVPdependent = True  # TODO
        # Attribute code_tree.
        skeleton.code_tree = deserialize_obj(skeleton.code_tree_serial)
        # Attribute codegen.
        skeleton.codegen = deserialize_obj(skeleton.codegen_serial)
        return skeleton


@attr.s
class ImplVariant:
    """Representation of an ImplVariant table database record.

    Attributes:
    -----------
    skeleton : int
        Used impl skeleton.
    kernels : list of int
        Used machine.
    db_id : int
        ID of associated impl variant database table record.
    """
    skeleton = attr.ib(type=int)
    kernels = attr.ib(type=List[int])
    kernels_serial = attr.ib(type=int, init=False)
    db_id = attr.ib(type=int, init=False)

    # Database information.
    db_table = Table('impl_variant', METADATA,
                     Column('db_id', Integer, primary_key=True),
                     Column('skeleton', Integer, ForeignKey('impl_skeleton.db_id')),
                     Column('kernels_serial', Integer),
                     Column('updatedIn', String, default=__version__),
                     Column('updatedOn', DateTime, default=datetime.now, onupdate=datetime.now),
                     Column('updatedBy', String, default=getuser(), onupdate=getuser()),
                     UniqueConstraint('skeleton', 'kernels_serial'),
                     sqlite_autoincrement=True)

    def to_database(self, db_session: Session):
        """Push this ECM record object to the database.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        ImplVariant
            Instance of this object connected to database session.
        """
        # Attribute kernels_serial.
        self.kernels_serial = ','.join(map(str, self.kernels))
        # Check if database already contains this ImplVariant object.
        variant = db_session.query(ImplVariant).filter(
            ImplVariant.skeleton.is_(self.skeleton), ImplVariant.kernels_serial.is_(self.kernels_serial)).one_or_none()
        if variant:
            # Supplement attributes not saved in database.
            variant.kernels = self.kernels
            return variant
        # Add new object to database.
        # Attribute kernels_serial.
        self.kernels_serial = ','.join(map(str, self.kernels))
        # Insert ImplVariantRecord object.
        insert(db_session, self)
        return self

    @staticmethod
    def select_all(db_session: Session) -> List['ImplVariant']:
        """Retrieve all ImplVariant table data record(s) from the database.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        list of ImplVariant
            Retrieved list of data records.
        """
        data = db_session.query(ImplVariant).all()
        return data

    @staticmethod
    def select(db_session: Session, variant_ids: List[int]) -> List['ImplVariant']:
        """
        Retrieve the ImplVariant table data record(s) from the database that match the given implementation variant IDs.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        variant_ids : list of int
            IDs of the impl variants requested.

        Returns:
        --------
        list of ImplVariant
            Retrieved list of data records.
        """
        variants = db_session.query(ImplVariant).filter(ImplVariant.db_id.in_(variant_ids)).all()
        if len(variants) != len(variant_ids):
            raise RuntimeError('Unable to select all requested variants!')
        # Deserialize attributes...
        for variant in variants:
            # Attribute kernels.
            variant.kernels = variant.kernels_serial.split(',')
        return variants
