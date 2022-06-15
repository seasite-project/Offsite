"""@package descriptions.impl.impl_skeleton
Definition of class ImplSkeleton.

@author: Johannes Seiferth
"""

from copy import deepcopy
from datetime import datetime
from getpass import getuser
from typing import Dict, List

import attr
from pathlib2 import Path
from sqlalchemy import Boolean, Column, DateTime, Enum, Integer, String, Table
from sqlalchemy.exc import NoResultFound, MultipleResultsFound
from sqlalchemy.orm import Session

from offsite import __version__
from offsite.codegen.code_dsl.code_dsl import parse_lark_grammar
from offsite.codegen.code_dsl.code_tree import CodeTree, CodeTreeGenerator
from offsite.config import ModelToolType
from offsite.database import METADATA, insert
from offsite.descriptions.impl.kernel_template import KernelTemplate
from offsite.descriptions.parser import load_yaml, deserialize_obj, serialize_obj, ComputationDict
from offsite.util.file_extensions import __impl_skeleton_ext__


@attr.s
class ImplSkeleton:
    """Representation of an ImplSkeleton object.

    Attributes:
    -----------
    name: str
        Name of this object.
    family: str
        Family name of this object. Can be used to group skeletons (E.g. group compatible skeletons for load balancing.)
    code_tree: CodeTree
        Tree representation of the kernel code.
    communicationOperationsClusterLvl : dict
        Cluster-level communication operations needed when executing a single iteration step.
    communicationOperationsNodeLvl : dict
        Node-level communication operations needed when executing a single iteration step.
    kernels: dict
        Kernels used in the impl code.
    connected_templates : list
        List of KernelTemplate objects associated with this object.
    db_id: int
        ID of associated ImplSkeleton database table record.
    """
    name = attr.ib(type=str)
    family = attr.ib(type=str)
    code_tree = attr.ib(type=CodeTree, hash=False)
    communicationOperationsClusterLvl = attr.ib(type=Dict)
    communicationOperationsNodeLvl = attr.ib(type=Dict)
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
                     Column('family', String),
                     Column('modelTool', Enum(ModelToolType)),
                     Column('connected_templates_ids', String),
                     Column('code_tree_serial', String),
                     Column('codegen_serial', String),
                     Column('isIVPdependent', Boolean),
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
        # Attribute family
        family = yaml['family'] if 'family' in yaml else 'default'
        # add end
        # Attribute code_tree.
        lark_code = parse_lark_grammar(yaml['code'])
        code_tree = CodeTreeGenerator().generate(lark_code, ComputationDict())
        code_tree = deepcopy(code_tree)
        # Parse code tree to collect information for attributes communication_operations and kernels.
        # .. communication operations.
        comm_ops_clust_lvl = dict()
        comm_ops_clust_lvl = CodeTree.count_clust_lvl_communication(code_tree.root, comm_ops_clust_lvl)
        comm_ops_node_lvl = dict()
        comm_ops_node_lvl = CodeTree.count_node_lvl_communication(code_tree.root, comm_ops_node_lvl)
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
        return cls(
            name, family, code_tree, comm_ops_clust_lvl, comm_ops_node_lvl, kernels, connected_templates, codegen)

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
            skeleton: ImplSkeleton = db_session.query(ImplSkeleton).filter(
                ImplSkeleton.name.like(impl_skeleton_name)).one()
        except NoResultFound:
            raise RuntimeError('Unable to load ImplSkeleton object from database!')
        except MultipleResultsFound:
            raise RuntimeError('Unable to load ImplSkeleton object from database!')
        # Attribute code_tree.
        skeleton.code_tree = deserialize_obj(skeleton.code_tree_serial)
        # Attribute communicationOperations.
        comm_ops_clust_lvl = dict()
        skeleton.communicationOperationsClusterLvl = CodeTree.count_clust_lvl_communication(
            skeleton.code_tree.root, comm_ops_clust_lvl)
        comm_ops_node_lvl = dict()
        skeleton.communicationOperationsNodeLvl = CodeTree.count_node_lvl_communication(
            skeleton.code_tree.root, comm_ops_node_lvl)
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
            comm_ops_clust_lvl = dict()
            skeleton.communicationOperationsClusterLvl = CodeTree.count_clust_lvl_communication(
                skeleton.code_tree.root, comm_ops_clust_lvl)
            comm_ops_node_lvl = dict()
            skeleton.communicationOperationsNodeLvl = CodeTree.count_node_lvl_communication(
                skeleton.code_tree.root, comm_ops_node_lvl)
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
        # Attribute code_tree.
        skeleton.code_tree = deserialize_obj(skeleton.code_tree_serial)
        # Attribute codegen.
        skeleton.codegen = deserialize_obj(skeleton.codegen_serial)
        # Parse code tree to collect information for attributes communication_operations and kernels.
        # .. communication operations.
        comm_ops_clust_lvl = dict()
        skeleton.communicationOperationsClusterLvl = CodeTree.count_clust_lvl_communication(
            skeleton.code_tree.root, comm_ops_clust_lvl)
        comm_ops_node_lvl = dict()
        skeleton.communicationOperationsNodeLvl = CodeTree.count_node_lvl_communication(
            skeleton.code_tree.root, comm_ops_node_lvl)
        # .. kernels.
        kernels = dict()
        skeleton.kernels = CodeTree.count_kernel(skeleton.code_tree.root, kernels)
        # Attribute connected_templates.
        # TODO: Currently not retrieved since not needed by caller of select
        # TODO read from database like kernel.select
        return skeleton

    @staticmethod
    def select_all(db_session: Session) -> List['ImplSkeleton']:
        """Retrieve all ImplSkeleton table data records  from the database.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        skeleton_id: int

        Returns:
        --------
        list of ImplSkeletons
            List of retrieved data records.
        """
        templates = KernelTemplate.select_all(db_session)
        skeleton_names: List[str] = [x[0] for x in db_session.query(ImplSkeleton.name).all()]
        return [ImplSkeleton.from_database(db_session, name, templates) for name in skeleton_names]

    @staticmethod
    def fetch_impl_skeleton_name(db_session: Session, db_id: int):
        try:
            name: str = db_session.query(ImplSkeleton.name).filter(ImplSkeleton.db_id.like(db_id)).one()[0]
        except NoResultFound:
            raise RuntimeError('Unable to load ImplSkeleton object from database!')
        except MultipleResultsFound:
            raise RuntimeError('Unable to load unique ImplSkeleton object from database!')
        return name


def parse_impl_skeletons(folder: Path, templates: List[KernelTemplate], tool: ModelToolType) -> List[ImplSkeleton]:
    """Parse folder for ImplSkeleton objects.

    Parameters:
    -----------
    folder: Path
        Relative path to folder that contains the ImplSkeleton objects.
    templates: list of KernelTemplate
        List of available KernelTemplate objects.
    tool: ModelToolType
        Limit parsing to ImplSkeleton objects that use this performance modelling tool.

    Returns:
    list of ImplSkeleton
        ImplSkeleton objects to be used.
    """
    # Filter impl skeletons.
    if tool:
        skeletons = list()
        for f in (f for f in folder.iterdir() if f.suffix == __impl_skeleton_ext__):
            try:
                skeleton = ImplSkeleton.from_yaml(f, templates)
                if skeleton.modelTool == tool:
                    skeletons.append(skeleton)
            except RuntimeError:
                pass
        return skeletons
    return [ImplSkeleton.from_yaml(f, templates) for f in (f for f in folder.iterdir()
                                                           if f.suffix == __impl_skeleton_ext__)]
