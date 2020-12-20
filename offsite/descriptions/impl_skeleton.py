"""@package impl_skeleton
Definition of class ImplSkeleton.
"""

from copy import deepcopy
from datetime import datetime
from getpass import getuser
from pathlib import Path
from typing import Dict, List

import attr
from pandas import DataFrame
from sqlalchemy import Column, DateTime, Enum, Integer, String, Table
from sqlalchemy.orm import Session
from sqlalchemy.orm.exc import NoResultFound, MultipleResultsFound

from offsite import __version__
from offsite.codegen.code_dsl import parse_lark_grammar
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
    name: str
        Name of this object.
    code_tree: CodeTree
        Tree representation of the kernel code.
    communicationOperations : dict
        Communication operations needed when executing a single iteration step.
    kernels: dict
        Kernels used in the implementation code.
    connected_templates : list
        List of KernelTemplate objects associated with this object.
    db_id: int
        ID of associated ImplSkeleton database table record.
    """
    name = attr.ib(type=str)
    code_tree = attr.ib(type=CodeTree, hash=False)
    communicationOperations = attr.ib(type=Dict)
    kernels = attr.ib(type=Dict)
    connected_templates = attr.ib(type=List['KernelTemplate'])
    codegen = attr.ib(type=Dict, default=list())
    modelTool = attr.ib(type=str, default=ModelToolType.KERNCRAFT)
    isIVPdependent = attr.ib(type=bool, default=False)
    connected_templates_ids = attr.ib(type=str, init=False)
    code_tree_serial = attr.ib(type=str, init=False)
    codegen_serial = attr.ib(type=str, init=False)
    db_id = attr.ib(type=int, init=False)

    # Database information.
    db_table = Table('impl_skeleton', METADATA,
                     Column('db_id', Integer, primary_key=True),
                     Column('name', String, unique=True),
                     Column('modelTool', Enum(ModelToolType)),
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
        yaml_path: Path
            Relative path to this object's YAML file.
        kernel_templates: list of KernelTemplate
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
        # add end
        # Attribute code_tree.
        lark_code = parse_lark_grammar(yaml['code'])
        code_tree = CodeTreeGenerator().generate(lark_code)
        code_tree = deepcopy(code_tree)
        # Parse code tree to collect information for attributes communication_operations and kernels.
        # .. communication operations.
        communication_operations = dict()
        communication_operations = CodeTree.count_communication(code_tree.root, communication_operations)
        # .. kernels.
        kernels = dict()
        kernels = CodeTree.count_kernel(code_tree.root, kernels)
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
        return cls(name, code_tree, communication_operations, kernels, connected_templates, codegen)

    @classmethod
    def from_database(cls, db_session: Session, impl_skeleton_name: str,
                      kernel_templates: List[KernelTemplate]) -> 'ImplSkeleton':
        """Construct ImplSkeleton object from database record.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        impl_skeleton_name: str
            Name of the impl skeleton, which is used as primary key in the database.
        kernel_templates: list of KernelTemplate
            List of available KernelTemplate objects.

        Returns:
        --------
        ImplSkeleton
            Created ImplSkeleton object.
        """
        try:
            skeleton: ImplSkeleton = db_session.query(
                ImplSkeleton).filter(ImplSkeleton.name.like(impl_skeleton_name)).one()
        except NoResultFound:
            raise RuntimeError('Unable to load ImplSkeleton object from database!')
        except MultipleResultsFound:
            raise RuntimeError('Unable to load ImplSkeleton object from database!')
        # Attribute code_tree.
        skeleton.code_tree = deserialize_obj(skeleton.code_tree_serial)
        # Attribute communicationOperations.
        communication_operations = dict()
        skeleton.communicationOperations = CodeTree.count_communication(
            skeleton.code_tree.root, communication_operations)
        # Attribute kernels.
        kernels = dict()
        skeleton.kernels = CodeTree.count_kernel(skeleton.code_tree.root, kernels)
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

    def to_database(self, db_session: Session) -> 'ImplSkeleton':
        """Push this impl skeleton object to the database.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        ImplSkeleton
            Instance of this object connected to database session.
        """
        # Check if database already contains this KernelTemplate object.
        skeleton: ImplSkeleton = db_session.query(ImplSkeleton).filter(ImplSkeleton.name.like(self.name)).one_or_none()
        if skeleton:
            # Supplement attributes not saved in database.
            skeleton.connected_templates = self.connected_templates
            for template in skeleton.connected_templates:
                if template[0].isIVPdependent:
                    skeleton.isIVPdependent = True
                    break
            skeleton.code_tree = self.code_tree
            skeleton.codegen = self.codegen
            # ... communication operations.
            communication_operations = dict()
            skeleton.communicationOperations = CodeTree.count_communication(
                skeleton.code_tree.root, communication_operations)
            # .. kernels.
            kernels = dict()
            skeleton.kernels = CodeTree.count_kernel(skeleton.code_tree.root, kernels)
            # Update serialized members.
            skeleton.code_tree_serial = serialize_obj(skeleton.code_tree)
            skeleton.codegen_serial = serialize_obj(skeleton.codegen)
            return skeleton
        # Add new object to database.
        # Serialize data.
        self.code_tree_serial = serialize_obj(self.code_tree)
        self.codegen_serial = serialize_obj(self.codegen)
        # IDs of connected templates.
        self.connected_templates_ids = ','.join(sorted([str(connect[0].db_id) for connect in self.connected_templates]))
        # ImplSkeleton object.
        insert(db_session, self)
        return self

    def compute_communication_costs(self, bench_data: DataFrame, method: ODEMethod, ivp: IVP) -> Dict[int, float]:
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
        costs: Dict[int, float] = dict()
        # Compute for all number of cores...
        for cores, row in bench_data.iterrows():
            costs[cores] = 0.0
            operation = row['name']
            # frequency = row['frequency']
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
        skeleton_id: int
            ID of the impl skeleton requested.

        Returns:
        --------
        ImplSkeleton
            Retrieved data record.
        """
        skeleton: ImplSkeleton = db_session.query(ImplSkeleton).filter(ImplSkeleton.db_id.is_(skeleton_id)).one()
        # Deserialize attributes...
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
        # Parse code tree to collect information for attributes communication_operations and kernels.
        # .. communication operations.
        communication_operations = dict()
        skeleton.communicationOperations = CodeTree.count_communication(
            skeleton.code_tree.root, communication_operations)
        # .. kernels.
        kernels = dict()
        skeleton.kernels = CodeTree.count_kernel(skeleton.code_tree.root, kernels)
        return skeleton
