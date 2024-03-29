"""@package descriptions.impl.kernel_template
Definitions of classes KernelTemplate, Kernel and PModelKernel.

@author: Johannes Seiferth
"""

from copy import deepcopy
from datetime import datetime
from getpass import getuser
from itertools import product
from sys import maxsize as sys_maxsize
from typing import Dict, List, Optional, Set, Tuple

import attr
from pandas import read_sql_query, DataFrame
from pathlib2 import Path
from sqlalchemy import Boolean, Column, DateTime, Enum, Float, ForeignKey, Integer, String, Table, UniqueConstraint
from sqlalchemy.exc import NoResultFound, MultipleResultsFound
from sqlalchemy.orm import Query, Session

import offsite.config
from offsite import __version__
from offsite.codegen.code_dsl.code_dsl import parse_lark_grammar
from offsite.codegen.code_dsl.code_node import CodeNode, CodeNodeType, RootNode
from offsite.codegen.code_dsl.code_tree import CodeTree, CodeTreeGenerator
from offsite.codegen.codegen_util import write_codes_to_file
from offsite.codegen.generator.kerncraft_generator import KerncraftCodeGenerator
from offsite.codegen.generator.kernel_bench_generator import KernelBenchCodeGenerator, KernelBenchFiles
from offsite.codegen.generator.yasksite_generator import YasksiteCodeGenerator
from offsite.config import Config, ModelToolType, ProgramModeType
from offsite.database import METADATA, insert
from offsite.descriptions.machine import MachineState
from offsite.descriptions.ode import ODEMethod, corrector_steps, stages
from offsite.descriptions.ode.ivp import IVP, ivp_grid_size
from offsite.descriptions.parser import load_yaml, ComputationDict, DatastructDict, DatastructDesc, \
    DatastructType, deserialize_obj, serialize_obj
from offsite.util.file_extensions import __kernel_template_ext__
from offsite.util.math_utils import eval_math_expr, solve_equation
from offsite.util.sample_interval import SampleInterval, SampleType, create_samples_lower_border_working_set, \
    create_samples_upper_border_working_set, create_samples_border_memory_lvl, create_samples_memory_lvl
from offsite.util.working_sets.count_references import CountReferences

# Extension of kernel files."""
__kernel_ext__ = '.c'


@attr.s
class KernelTemplate:
    """Representation of a KernelTemplate object.

    Attributes:
    -----------
    name: str
        Name of this object.
    variants: list of Kernel
        Kernel variants of this template.
    isIVPdependent: boolean
        If true this template contains IVP calls.
    modelTool: ModelToolType
        Performance model tool used to obtain performance data of this object.
    datastructs: dict
        Datastructs used by this object.
    computations: dict
        Computations used by this object.
    codegen: dict
        Code generation parameters used to generate code of this object.
    db_id: int
        ID of associated KernelTemplate database table record.
    """
    name = attr.ib(type=str)
    variants = attr.ib(type=List['Kernel'])
    isIVPdependent = attr.ib(type=bool)
    modelTool = attr.ib(type=ModelToolType)
    datastructs = attr.ib(type=DatastructDict)
    computations = attr.ib(type=ComputationDict)
    codegen = attr.ib(type=Dict, default=dict())
    datastructs_serial = attr.ib(type=str, init=False)
    computations_serial = attr.ib(type=str, init=False)
    codegen_serial = attr.ib(type=str, init=False)
    db_id = attr.ib(type=int, init=False)

    # Database information.
    db_table = Table('kernel_template', METADATA,
                     Column('db_id', Integer, primary_key=True),
                     Column('name', String, unique=True),
                     Column('isIVPdependent', Boolean),
                     Column('modelTool', Enum(ModelToolType)),
                     Column('datastructs_serial', String),
                     Column('computations_serial', String),
                     Column('codegen_serial', String),
                     Column('updatedIn', String, default=__version__),
                     Column('updatedOn', DateTime, default=datetime.now, onupdate=datetime.now),
                     Column('updatedBy', String, default=getuser(), onupdate=getuser()),
                     sqlite_autoincrement=True)

    @classmethod
    def from_yaml(cls, yaml_path: Path) -> 'KernelTemplate':
        """Construct KernelTemplate object from YAML definition.

        Parameters:
        -----------
        yaml_path: Path
            Relative path to this object's YAML file.

        Returns:
        --------
        KernelTemplate
            Created KernelTemplate object.
        """
        # Load YAML data.
        yaml = load_yaml(yaml_path)
        # Attribute name.
        name = yaml_path.stem
        # Attribute variants.
        # .. used to temporarily store yaml data for post init.
        variants = yaml['variants']
        # Attribute model_tool.
        tool_str = yaml['model_tool']
        if tool_str == ModelToolType.KERNCRAFT.value:
            model_tool = ModelToolType.KERNCRAFT
        elif tool_str == ModelToolType.YASKSITE.value:
            model_tool = ModelToolType.YASKSITE
        else:
            assert False
        # Attribute computations.
        computations = ComputationDict.from_data(yaml['computations'])
        # Attribute is_ivp_dependent.
        if 'ivp_dependent' in yaml:
            is_ivp_dependent = yaml['ivp_dependent']
        else:
            is_ivp_dependent = any('%RHS' in c.computation for c in computations.values())
        # Attribute datastructs.
        datastructs = DatastructDict.from_data(yaml['datastructs'])
        # ... implicitly add additional datastructures required by IVP-dependent kernels.
        if is_ivp_dependent and model_tool == ModelToolType.KERNCRAFT:
            datastructs['h'] = DatastructDesc('double', DatastructType.scalar, '1', False)
            datastructs['t'] = DatastructDesc('double', DatastructType.scalar, '1', False)
            # datastructs['g'] = DatastructDesc('double', DatastructType.scalar, '1', False)
            datastructs['c'] = DatastructDesc('double', DatastructType.array1D, ['s'], False)
        # Attribute codegen.
        codegen = dict()
        if 'codegen' in yaml:
            for key, value in yaml['codegen'].items():
                if key == 'loop splits':
                    codegen[key] = [v for v in value]
                else:
                    codegen[key] = value
        # Create object.
        return cls(name, variants, is_ivp_dependent, model_tool, datastructs, computations, codegen)

    @classmethod
    def from_database(cls, db_session: Session, kernel_template_name: str) -> 'KernelTemplate':
        """Construct KernelTemplate object from database record.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        kernel_template_name: str
            Name of the kernel template, which is used as primary key in the database.

        Returns:
        --------
        KernelTemplate
            Created KernelTemplate object.
        """
        try:
            kernel_template: KernelTemplate = db_session.query(KernelTemplate).filter(
                KernelTemplate.name.like(kernel_template_name)).one()
        except NoResultFound:
            raise RuntimeError('Unable to load KernelTemplate object from database!')
        except MultipleResultsFound:
            raise RuntimeError('Unable to load KernelTemplate object from database!')
        # Attribute datastructs.
        kernel_template.datastructs = deserialize_obj(kernel_template.datastructs_serial)
        # Attribute computations.
        kernel_template.computations = deserialize_obj(kernel_template.computations_serial)
        # Attribute codegen.
        kernel_template.codegen = deserialize_obj(kernel_template.codegen_serial)
        # Attribute variants.
        try:
            tid = kernel_template.db_id
            knames = [x[0] for x in db_session.query(Kernel.name).filter(Kernel.template_db_id.like(tid)).all()]
            kernel_template.variants = [Kernel.from_database(db_session, name, tid) for name in knames]
        except NoResultFound:
            raise RuntimeError('Unable to load Kernel objects from database!')
        return kernel_template

    def __attrs_post_init__(self):
        """Create the Kernel objects associated with this object.

        Parameters:
        -----------
        -

        Returns:
        --------
        -
        """
        yaml = self.variants
        self.variants = [Kernel.from_yaml(variant, self) for variant in yaml]
        # Create kernel variants for all available Yasksite optimization parameter permutations.
        if self.modelTool == ModelToolType.YASKSITE:
            config: Config = offsite.config.offsiteConfig
            variants = list()
            for variant in self.variants:
                # Add specialized kernel versions that apply the given optimization parameters.
                optimizations = [p for p in
                                 product(config.args.scenario.ys_blockings, config.args.scenario.ys_foldings)]
                for blocking, folding in optimizations:
                    # Deepcopy original variant object.
                    variant_copy = deepcopy(variant)
                    # Adjust kernel name and specify its optimization parameters.
                    variant_copy.name = variant_copy.name + '_' + blocking
                    if folding:
                        variant_copy.name = variant_copy.name + '_folding{}'.format(folding)
                    variant_copy.optimization_parameters['ys_blocking'] = blocking
                    variant_copy.optimization_parameters['ys_folding'] = folding
                    variant_copy.template = self
                    variants.append(variant_copy)
            self.variants = variants

    def to_database(self, db_session: Session) -> 'KernelTemplate':
        """Push this kernel template object to the database.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        KernelTemplate
            Instance of this object connected to database session.
        """
        # Check if database already contains this KernelTemplate object.
        template: KernelTemplate = db_session.query(KernelTemplate).filter(
            KernelTemplate.name.like(self.name)).one_or_none()
        if template:
            # Supplement attributes not saved in database.
            template.datastructs = self.datastructs
            template.computations = self.computations
            template.codegen = self.codegen
            for variant in self.variants:
                variant.template_db_id = template.db_id
            template.variants = [variant.to_database(db_session) for variant in self.variants]
            # Update serialized members.
            template.datastructs_serial = serialize_obj(template.datastructs)
            template.computations_serial = serialize_obj(template.computations)
            template.codegen_serial = serialize_obj(template.codegen)
            return template
        # Add new object to database.
        # Serialize data.
        self.datastructs_serial = self.datastructs.serialize()
        self.computations_serial = self.computations.serialize()
        self.codegen_serial = serialize_obj(self.codegen)
        # KernelTemplate object.
        insert(db_session, self)
        # Kernel objects.
        for variant in self.variants:
            variant.template_db_id = self.db_id
            variant.to_database(db_session)
        return self

    @staticmethod
    def select_all(db_session: Session) -> List['KernelTemplate']:
        """Retrieve all kernel template table data records  from the database.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        list of KernelTemplate
            List of retrieved data records.
        """
        template_names: List[str] = [x[0] for x in db_session.query(KernelTemplate.name).all()]
        return [KernelTemplate.from_database(db_session, name) for name in template_names]


def parse_kernel_templates(folder: Path, tool: ModelToolType) -> List[KernelTemplate]:
    """Parse folder for KernelTemplate objects.

    Parameters:
    -----------
    folder: Path
        Relative path to folder that contains the KernelTemplate objects.
    tool: ModelToolType
        Limit parsing to KernelTemplate objects that use this performance modelling tool.

    Returns:
    --------
    list of KernelTemplate
        KernelTemplate objects to be used.
    """
    templates = [KernelTemplate.from_yaml(kf) for kf in [f for f in folder.iterdir()
                                                         if f.suffix == __kernel_template_ext__]]
    # Filter kernel templates.
    if tool:
        templates = [template for template in templates if template.modelTool == tool]
    return templates


@attr.s(hash=True)
class Kernel:
    """Representation of a Kernel object.

    Attributes:
    -----------
    name: str
        Name of this object.
    template: KernelTemplate
        Reference to associated KernelTemplate object.
    pmodel_kernels: list of PmodelKernel
        PmodelKernel objects associated with this object.
    code_tree: CodeTree
        Tree representation of the kernel code.
    communicationOperationsClusterLvl : dict
        Cluster-level communication operations needed when executing a single iteration step.
    communicationOperationsNodeLvl : dict
        Node-level communication operations needed when executing a single iteration step.
    codegen: dict
        Code generation parameters used to generate code of this object.
    optimization_parameters: dict (key=str, value=str)
        Applied YaskSite optimization parameters.
    db_id: int
        ID of associated Kernel database table record.
    """
    name = attr.ib(type=str)
    template = attr.ib(type=KernelTemplate, hash=False)
    pmodel_kernels = attr.ib(type=List['PModelKernel'], hash=False)
    code_tree = attr.ib(type=CodeTree, hash=False)
    communicationOperationsClusterLvl = attr.ib(type=Dict)
    communicationOperationsNodeLvl = attr.ib(type=Dict)
    used_datastructs = attr.ib(type=Set)
    codegen = attr.ib(type=Dict, default=dict())
    optimization_parameters = attr.ib(type=Dict[str, str], default=dict(), hash=False)
    code_tree_serial = attr.ib(type=str, init=False)
    codegen_serial = attr.ib(type=str, init=False)
    optimization_parameters_serial = attr.ib(type=str, init=False, hash=False)
    template_db_id = attr.ib(type=int, init=False, hash=False)
    db_id = attr.ib(type=int, init=False)

    # Database information.
    db_table = Table('kernel', METADATA,
                     Column('db_id', Integer, primary_key=True),
                     Column('name', String),
                     Column('template_db_id', Integer, ForeignKey('kernel_template.db_id')),
                     Column('code_tree_serial', String),
                     Column('codegen_serial', String),
                     Column('optimization_parameters_serial', String),
                     Column('updatedIn', String, default=__version__),
                     Column('updatedOn', DateTime, default=datetime.now, onupdate=datetime.now),
                     Column('updatedBy', String, default=getuser(), onupdate=getuser()),
                     UniqueConstraint('name', 'template_db_id'),
                     sqlite_autoincrement=True)

    @classmethod
    def from_yaml(cls, yaml: dict, template: KernelTemplate) -> 'Kernel':
        """Construct Kernel object from YAML definition.

        Parameters:
        -----------
        yaml: dict
            YAML object describing this object.
        template: KernelTemplate
            Reference to associated KernelTemplate object.

        Returns:
        --------
        Kernel
            Created Kernel object.
        """
        # Attribute name.
        name = yaml['name']
        # Attribute code_tree.
        lark_code = parse_lark_grammar(yaml['code'])
        code_tree = CodeTreeGenerator().generate(lark_code, template.computations)
        code_tree.visit(code_tree.root)
        code_tree = deepcopy(code_tree)
        # Attribute pmodel_kernels.
        # .. used to temporarily store yaml data for post init.
        pmodel_kernels = yaml
        # Parse code tree to collect information for attributes ...
        # ... communication_operations.
        comm_ops_clust_lvl = dict()
        comm_ops_clust_lvl = CodeTree.count_clust_lvl_communication(code_tree.root, comm_ops_clust_lvl)
        comm_ops_node_lvl = dict()
        comm_ops_node_lvl = CodeTree.count_node_lvl_communication(code_tree.root, comm_ops_node_lvl)
        # ... used_datastructs.
        used_datastructs = set()
        used_datastructs = CodeTree.collect_referenced_datastructs(code_tree.root, used_datastructs)
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
            name, template, pmodel_kernels, code_tree, comm_ops_clust_lvl, comm_ops_node_lvl, used_datastructs, codegen)

    @classmethod
    def from_database(cls, db_session: Session, kernel_name: str, kernel_template_id: int) -> 'Kernel':
        """Construct Kernel object from database record.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        kernel_name: str
            Name of the kernel, which is used as primary key in the database.
        template_id: IVP
            Database ID of the associated KernelTemplate object.

        Returns:
        --------
        Kernel
            Created Kernel object.
        """
        try:
            kernel: Kernel = db_session.query(Kernel).filter(
                Kernel.name.is_(kernel_name), Kernel.template_db_id.is_(kernel_template_id)).one()
        except NoResultFound:
            raise RuntimeError('Unable to load Kernel object from database!')
        except MultipleResultsFound:
            raise RuntimeError('Unable to load Kernel object from database!')
        # Attribute code_tree.
        kernel.code_tree = deserialize_obj(kernel.code_tree_serial)
        # Attribute communicationOperations.
        comm_ops_clust_lvl = dict()
        kernel.communicationOperationsClusterLvl = CodeTree.count_clust_lvl_communication(
            kernel.code_tree.root, comm_ops_clust_lvl)
        comm_ops_node_lvl = dict()
        kernel.communicationOperationsNodeLvl = CodeTree.count_node_lvl_communication(
            kernel.code_tree.root, comm_ops_node_lvl)
        # Attribute used_datastructs.
        used_datastructs = set()
        kernel.used_datastructs = CodeTree.collect_referenced_datastructs(kernel.code_tree.root, used_datastructs)
        # Attribute optimization_parameters.
        kernel.optimization_parameters = deserialize_obj(kernel.optimization_parameters_serial)
        # Attribute codegen.
        kernel.codegen = deserialize_obj(kernel.codegen_serial)
        # Attribute template.
        kernel.template = db_session.query(KernelTemplate).filter(KernelTemplate.db_id.is_(kernel_template_id)).one()
        # Attribute pmodel_kernels.
        try:
            kernel.pmodel_kernels = db_session.query(PModelKernel).filter(
                PModelKernel.kernel_db_id.like(kernel.db_id)).all()
        except NoResultFound:
            raise RuntimeError('Unable to load Kernel object from database!')
        return kernel

    def __attrs_post_init__(self):
        """Create the PmodelKernel objects associated with this object.

        Parameters:
        -----------
        -

        Returns:
        --------
        -
        """
        yaml = self.pmodel_kernels
        # Number of pmodel kernels.
        num_pmodel_kernels = 1 + CodeTree.count_pmodel(self.code_tree.root)
        # No need for name affix if there is only one PModelKernel.
        if num_pmodel_kernels == 1:
            self.pmodel_kernels = [PModelKernel.from_yaml(yaml, self, self.code_tree.root)]
        else:
            self.pmodel_kernels = list()
            for idx in range(num_pmodel_kernels):
                if idx == 0:
                    node = self.code_tree.root
                else:
                    node = CodeTree.find_pmodel_node(self.code_tree.root, idx)[0]
                assert node is not None
                self.pmodel_kernels.append(PModelKernel.from_yaml(yaml, self, node, idx))

    def to_database(self, db_session: Session) -> 'Kernel':
        """Push this kernel object to the database.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        Kernel
            Instance of this object connected to database session.
        """
        # Check if database already contains this Kernel object.
        kernel: Kernel = db_session.query(Kernel).filter(
            Kernel.name.like(self.name), Kernel.template_db_id.is_(self.template_db_id)).one_or_none()
        if kernel:
            # Supplement attributes not saved in database.
            kernel.optimization_parameters = self.optimization_parameters
            kernel.code_tree = self.code_tree
            kernel.codegen = self.codegen
            kernel.template = self.template
            for pmodel in self.pmodel_kernels:
                pmodel.kernel_db_id = kernel.db_id
            kernel.pmodel_kernels = [pmodel.to_database(db_session) for pmodel in self.pmodel_kernels]
            # ... communication operations.
            comm_ops_clust_lvl = dict()
            kernel.communicationOperationsClusterLvl = CodeTree.count_clust_lvl_communication(
                kernel.code_tree.root, comm_ops_clust_lvl)
            comm_ops_node_lvl = dict()
            kernel.communicationOperationsNodeLvl = CodeTree.count_node_lvl_communication(
                kernel.code_tree.root, comm_ops_node_lvl)
            # ... used_datastructures.
            used_datastructs = set()
            kernel.used_datastructs = CodeTree.collect_referenced_datastructs(kernel.code_tree.root, used_datastructs)
            # Update serialized members.
            kernel.optimization_parameters_serial = serialize_obj(kernel.optimization_parameters)
            kernel.code_tree_serial = serialize_obj(kernel.code_tree)
            kernel.codegen_serial = serialize_obj(kernel.codegen)
            return kernel
        # Add new object to database.
        insert(db_session, self)
        # Attribute optimization_parameters_serial.
        self.optimization_parameters_serial = serialize_obj(self.optimization_parameters)
        # Attribute code_tree.
        self.code_tree_serial = serialize_obj(self.code_tree)
        # Attribute codegen.
        self.codegen_serial = serialize_obj(self.codegen)
        # Attribute template_db_id.
        self.template_db_id = self.template.db_id
        # PModelKernel objects.
        for pmodel in self.pmodel_kernels:
            pmodel.kernel_db_id = self.db_id
            pmodel.to_database(db_session)
        return self

    def generate_pmodel_code(
            self, method: Optional[ODEMethod] = None, ivp: Optional[IVP] = None, system_size: Optional[int] = None):
        """Generate code for all pmodel kernels associated with this object.

        For each PModelKernel, either a generalized version of the code or if required a specialised version (unrolled
        loops, ...) is generated for the given configuration of ODEMethod and IVP.

        Parameters:
        -----------
        method: ODEMethod
            Used ODE method.
        ivp: IVP
            Used IVP.
        system_size: int
            Used ODE system size.

        Returns:
        --------
        -
        """
        # Generate code for kerncraft tool.
        if self.template.modelTool == ModelToolType.KERNCRAFT:
            # Generate pmodel codes.
            codes = KerncraftCodeGenerator().generate(self, method, ivp, system_size)
            codes = write_codes_to_file(codes, Path('tmp'), num_workers=1, format_code=False)
            # Link code files to PModelKernel objects.
            for pmodel in self.pmodel_kernels:
                pmodel_suffix = '_{}'.format(pmodel.name) if pmodel.name else ''
                ivp_suffix = '_{}'.format(ivp.name) if (ivp and self.template.isIVPdependent) else ''
                if self.template.isIVPdependent and len(ivp.components) > 1:
                    pmodel.code_path = dict()
                    for idx, name in enumerate(ivp.components):
                        pmodel.code_path[name] = Path(
                            'tmp/{}{}{}_{}.c'.format(self.name, pmodel_suffix, ivp_suffix, idx + 1))
                        if pmodel.code_path[name] not in codes:
                            raise RuntimeError('Unknown pmodel kernel: {}'.format(pmodel.code_path))
                else:
                    pmodel.code_path = Path('tmp/{}{}{}.c'.format(self.name, pmodel_suffix, ivp_suffix))
                    if pmodel.code_path not in codes:
                        raise RuntimeError('Unknown pmodel kernel: {}'.format(pmodel.code_path))
        # Generate code for yasksite tool.
        elif self.template.modelTool == ModelToolType.YASKSITE:
            if ivp and not ivp.characteristic.isStencil:
                raise RuntimeError('Can not generate YASKSITE code for {}: Not a stencil!'.format(ivp.name))
            # Generate Yasksite code.
            codes = YasksiteCodeGenerator().generate(self, method, ivp)
            codes = write_codes_to_file(codes, Path('tmp'), num_workers=1, format_code=False)
            # Link code files to PModelKernel objects.
            for pmodel in self.pmodel_kernels:
                pmodel_suffix = '_{}'.format(pmodel.name) if pmodel.name else ''
                ivp_suffix = '_{}'.format(ivp.name) if (ivp and self.template.isIVPdependent) else ''
                pmodel.code_path = Path('tmp/{}{}{}.c'.format(self.name, pmodel_suffix, ivp_suffix))
                if pmodel.code_path not in codes:
                    raise RuntimeError('Unknown pmodel kernel: {}'.format(pmodel.code_path))
        else:
            assert False

    def generate_rhs_bench_code(self, method: Optional[ODEMethod] = None,
                                ivp: Optional[IVP] = None) -> KernelBenchFiles:
        """Generate benchmarking code for this kernel object.

        Parameters:
        -----------
        method: ODEMethod
            Used ODE method.
        ivp: IVP
            Used IVP.
        system_size: int
            Used ODE system size.

        Returns:
        --------
        KernelBenchFiles
            Path to generated code files.
        """
        assert self.template.modelTool == ModelToolType.KERNCRAFT
        return KernelBenchCodeGenerator().generate(self, method, ivp)

    def deduce_relevant_samples(
            self, machine: MachineState, method: Optional[ODEMethod], ivp: Optional[IVP]) -> List[SampleInterval]:
        """
        Deduce a set of significant sample intervals from the working sets of the pmodel kernels associated with this
        kernel. The shape of the set depends on the working set sizes as well as the characteristics of the given
        machine, ODE method and IVP.

        Parameters:
        -----------
        machine: MachineState
            Trained machine.
        method: ODE method
            Trained ODE method.
        ivp: IVP
            Trained IVP.

        Returns:
        --------
        list of SampleInterval
            List of significant sample intervals.
        """
        # For each PModelKernel, determine the sample intervals required to give prediction for each of the
        # PModelKernel's working sets. These intervals are sorted by their start value in ascending order.
        if self.template.modelTool == ModelToolType.YASKSITE and ivp is not None and ivp.characteristic.isStencil:
            samples_pmk = [pmk.determine_samples_for_yasksite_mode(machine, ivp) for pmk in self.pmodel_kernels]
        else:
            samples_pmk = [pmk.determine_samples_from_working_sets(machine, method, ivp) for pmk in self.pmodel_kernels]
        if not samples_pmk:
            raise RuntimeError('ERROR: Kernel \'{}\' has no PModelKernel!'.format(self.name))
        # Deduce the set of sample intervals required to give prediction for this kernel object from the sample
        # intervals of its associated PModelKernel objects. When having multiple PModelKernel, additional intervals
        # might be necessary to reflect all possible combinations of PModelKernel working sets.
        samples = list()
        # Only a single PModelKernel, we can simply use those intervals.
        if len(samples_pmk) == 1:
            samples = samples_pmk.pop()
        # Multiple PModel, we might have to add addition intervals.
        # 1) Iterate over samples_pmk and determine the interval with the lowest end value.
        # 2) Add that interval to the list of significant sample intervals.
        # 3) Update start values of the other intervals to 'lowest end + 1'.
        # 4) If the start value of an interval is greater than its end value after the update, remove that interval.
        else:
            sample_type = None
            while any(len(le) > 1 for le in samples_pmk):
                end = sys_maxsize
                # Determine lowest end value of current samples.
                low_s = None
                for spk in samples_pmk:
                    low_s = spk[0] if spk[0].last < end else low_s
                    sample_type = spk[0].type if spk[0].last < end else sample_type
                    end = min(spk[0].last, end)
                # Add sample with this end value to list.
                med = low_s.median(ivp)
                if med is None:
                    med = low_s.first
                samples.append(SampleInterval(low_s.first, low_s.last, sample_type, med))
                # Update end of all current samples.
                for spk in samples_pmk:
                    spk[0].first = end + 1
                # Remove all current samples already covered by the list.
                for spk in samples_pmk:
                    if len(spk) > 1 and spk[0].first > spk[0].last:
                        spk.pop(0)
            # Add interval with 'inf' end.
            samples.append(samples_pmk[0][0])
        for sample in samples:
            assert (sample.first <= sample.last)
        samples.sort(key=lambda x: x.first)
        return samples

    @staticmethod
    def select(db_session: Session, kernel_id: int) -> 'Kernel':
        """Retrieve the Kernel table data record from the database that matches the given database ID.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        kernel_id: int
            ID of the requested Kernel object.

        Returns:
        --------
        Kernel
            Kernel object.
        """
        data: Kernel = db_session.query(Kernel).filter(Kernel.db_id.is_(kernel_id)).one()
        return Kernel.from_database(db_session, data.name, data.template_db_id)

    @staticmethod
    def select_kernels(db_session: Session, template: int) -> List[int]:
        """Retrieve IDs and optimization parameters of all kernel of the given KernelTemplate objects.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        template: int
            ID of the used KernelTemplate object.

        Returns:
        --------
        list of tuple of tuples (int, dict (key=str, value=str))
            Lists of tuples that contain a kernel ID and optimization parameters of a kernel.
        """
        data: List[Kernel] = db_session.query(Kernel).filter(Kernel.template_db_id.is_(template)).all()
        return [(kernel.db_id, deserialize_obj(kernel.optimization_parameters_serial)) for kernel in data]

    @staticmethod
    def fetch_kernel_name(db_session: Session, db_id: int) -> str:
        try:
            name: str = db_session.query(Kernel.name).filter(Kernel.db_id.like(db_id)).one()[0]
        except NoResultFound:
            raise RuntimeError('Unable to load kernel object from database!')
        except MultipleResultsFound:
            raise RuntimeError('Unable to load unique kernel object from database!')
        return name


@attr.s
class KernelRecord:
    """Representation of a kernel runtime prediction table database record.

    Attributes:
    -----------
    kernel: int
        Used kernel.
    machine: int
        Used machine.
    compiler: int
        Used compiler.
    method: int
        Used ODE method.
    ivp: int
        Used IVP.
    frequency: float
        Used CPU frequency.
    cores: int
        Used number of CPU cores.
    sample: SampleInterval
        Used sample interval.
    prediction: float
        Kernel runtime prediction.
    mode: str
        Prediction obtained in RUN or MODEL mode.
    db_id: int
        ID of associated kernel prediction database table record.
    """
    kernel = attr.ib(type=int)
    machine = attr.ib(type=int)
    compiler = attr.ib(type=int)
    method = attr.ib(type=int)
    ivp = attr.ib(type=int)
    frequency = attr.ib(type=float)
    cores = attr.ib(type=int)
    sample = attr.ib(type=SampleInterval)
    prediction = attr.ib(type=str)
    mode = attr.ib(type=str)
    weight = attr.ib(type=float, init=False, default=1.0)
    first = attr.ib(type=int, init=False)
    last = attr.ib(type=int, init=False)
    sample_type = attr.ib(type=SampleType, init=False)
    db_id = attr.ib(type=int, init=False)

    # Database information.
    db_table = Table('kernel_prediction', METADATA,
                     Column('db_id', Integer, primary_key=True),
                     Column('kernel', Integer, ForeignKey('kernel.db_id')),
                     Column('machine', Integer, ForeignKey('machine.db_id')),
                     Column('compiler', Integer, ForeignKey('compiler.db_id')),
                     Column('method', Integer, ForeignKey('ode_method.db_id')),
                     Column('ivp', Integer, ForeignKey('ivp.db_id')),
                     Column('frequency', Integer),
                     Column('cores', Integer),
                     Column('first', Integer),
                     Column('last', Integer),
                     Column('sample_type', Enum(SampleType)),
                     Column('prediction', String),
                     Column('mode', String),
                     Column('weight', Float),
                     Column('updatedIn', String, default=__version__),
                     Column('updatedOn', DateTime, default=datetime.now, onupdate=datetime.now),
                     Column('updatedBy', String, default=getuser(), onupdate=getuser()),
                     sqlite_autoincrement=True)

    def __attrs_post_init__(self):
        """Create this ECM Record object.

        Parameters:
        -----------
        -

        Returns:
        --------
        -
        """
        self.first = self.sample.first
        self.last = self.sample.last
        self.sample_type = self.sample.type

    def to_database(self, db_session: Session):
        """Push this ECM record object to the database.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        -
        """
        # Attribute prediction.
        self.prediction = str(self.prediction)
        # Insert KernelRecord object.
        insert(db_session, self)

    @staticmethod
    def update(db_session: Session, kernel: int, machine: int, compiler: int, method: int, ivp: int, cores: int,
               frequency: float, records: List[Tuple[SampleInterval, str]], mode: str):
        """Update data records in PModelRecord table with new data.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        machine: int
            Used machine.
        compiler: int
            Used compiler.
        method: int
            Used ODE method.
        ivp: int
            Used IVP.
        cores: int
            Used CPU cores.
        frequency: float
            Used CPU frequency.
        interval: SampleInterval
            Used sample interval.
        records: list of tuple (SampleInterval, kernel prediction)
            Kernel prediction data.
        mode: str
            Application mode used to obtain data.

        Returns:
        --------
        -
        """
        for record in records:
            interval: SampleInterval = record[0]
            prediction: str = record[1]
            # Select and update all already included data records.
            try:
                queried_record: KernelRecord = db_session.query(KernelRecord).filter(
                    KernelRecord.kernel.is_(kernel), KernelRecord.machine.is_(machine),
                    KernelRecord.compiler.is_(compiler), KernelRecord.method.is_(method), KernelRecord.ivp.is_(ivp),
                    KernelRecord.cores.is_(cores), KernelRecord.frequency.is_(frequency),
                    KernelRecord.first.is_(interval.first), KernelRecord.last.is_(interval.last),
                    KernelRecord.sample_type.is_(interval.type), KernelRecord.mode.is_(mode)).one()
                # Update data record.
                queried_record.prediction = str(prediction)
            except NoResultFound:
                # Insert new data record.
                record = KernelRecord(
                    kernel, machine, compiler, method, ivp, frequency, cores, interval, prediction, mode)
                record.to_database(db_session)
            except MultipleResultsFound:
                raise RuntimeError('Unable to update Kernel record!')

    @staticmethod
    def select(db_session: Session, machine: int, compiler: int, method: int, ivp: int, cores: int, frequency: float,
               mode: Optional[str] = None, ode_size: Optional[int] = None) -> DataFrame:
        """Retrieve KernelRecord table data record(s) from the database.

        Return all records that match the provided configuration of machine, compiler, ODE method, IVP, CPU frequency
        and number of cores.

        If a ODE system size 'ode_size' is passed only records for that fixed 'ode_size' are returned.

        If a prediction mode is passed only records obtained using that mode are returned.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        machine: int
            Used machine.
        compiler: int
            Used compiler.
        method: int
            Used ODE method.
        ivp: int
            Used IVP.
        cores: int
            Used number of cores.
        frequency: float
            Used CPU frequency.
        mode: str
            Application mode used to obtain data.
        ode_size: int
            Used fixed ODE system size.

        Returns:
        --------
        pandas.DataFrame
            Retrieved list of data records.
        """
        query: Query
        # Required ivp.db_id = -1 to include all none IVP-dependent kernels.
        if mode is None:
            if ode_size is None:
                query: Query = db_session.query(KernelRecord).filter(
                    KernelRecord.machine.is_(machine), KernelRecord.compiler.is_(compiler),
                    KernelRecord.method.is_(method), KernelRecord.ivp.in_([ivp, -1]),
                    KernelRecord.frequency.is_(frequency), KernelRecord.cores.is_(cores)).order_by(
                    KernelRecord.kernel.asc())
            else:
                query: Query = db_session.query(KernelRecord).filter(
                    KernelRecord.machine.is_(machine), KernelRecord.compiler.is_(compiler),
                    KernelRecord.method.is_(method), KernelRecord.ivp.in_([ivp, -1]),
                    KernelRecord.frequency.is_(frequency), KernelRecord.cores.is_(cores),
                    KernelRecord.first <= ode_size, KernelRecord.last >= ode_size).order_by(KernelRecord.kernel.asc())
        else:
            if ode_size is None:
                query: Query = db_session.query(KernelRecord).filter(
                    KernelRecord.machine.is_(machine), KernelRecord.compiler.is_(compiler),
                    KernelRecord.method.is_(method), KernelRecord.ivp.in_([ivp, -1]),
                    KernelRecord.frequency.is_(frequency), KernelRecord.cores.is_(cores),
                    KernelRecord.mode.is_(mode)).order_by(KernelRecord.kernel.asc())
            else:
                query: Query = db_session.query(KernelRecord).filter(
                    KernelRecord.machine.is_(machine), KernelRecord.compiler.is_(compiler),
                    KernelRecord.method.is_(method), KernelRecord.ivp.in_([ivp, -1]),
                    KernelRecord.frequency.is_(frequency), KernelRecord.cores.is_(cores), KernelRecord.mode.is_(mode),
                    KernelRecord.first <= ode_size, KernelRecord.last >= ode_size).order_by(KernelRecord.kernel.asc())
        return read_sql_query(query.statement, db_session.bind)

    @staticmethod
    def contains(db_session: Session, kernel: int, machine: int, compiler: int, method: int, ivp: int, frequency: float,
                 max_cores: int, mode: str, ode_size: int, sample_types: Optional[List[SampleType]]) -> bool:
        """
        Check if the KernelRecord table already contains data of a particular kernel for a given configuration of
        machine, compiler, ODE method, IVP, CPU frequency, and number of CPU cores.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        kernel: int
            Used kernel.
        machine: int
            Used machine.
        compiler: int
            Used compiler.
        method: int
            Used ODE method.
        ivp: int
            Used IVP.
        frequency: float
            Used CPU frequency.
        max_cores: int
            Maximum number of cores. Check for core counts up to this value.
        mode: str
            Application mode used to obtain data.
        ode_size: int
            Used ODE size or None.
        sample_type: List of SampleType
            Only include data obtained using a specific set of sample types.

        Returns:
        --------
        boolean
            True if table contains fitting data False else.
        """
        # ODE size is fixed.
        if ode_size:
            # Include data for all sample types.
            if sample_types is None:
                data = db_session.query(KernelRecord).filter(
                    KernelRecord.kernel.is_(kernel), KernelRecord.machine.is_(machine),
                    KernelRecord.compiler.is_(compiler), KernelRecord.method.is_(method), KernelRecord.ivp.is_(ivp),
                    KernelRecord.frequency.is_(frequency), KernelRecord.cores.in_((c for c in range(1, max_cores + 1))),
                    KernelRecord.first.is_(ode_size), KernelRecord.last.is_(ode_size),
                    KernelRecord.mode.is_(mode)).all()
                return bool(len(data) == max_cores)
            # Include only data for the specified sample types.
            data = db_session.query(KernelRecord).filter(
                KernelRecord.kernel.is_(kernel), KernelRecord.machine.is_(machine), KernelRecord.compiler.is_(compiler),
                KernelRecord.method.is_(method), KernelRecord.ivp.is_(ivp), KernelRecord.frequency.is_(frequency),
                KernelRecord.cores.in_((c for c in range(1, max_cores + 1))), KernelRecord.first.is_(ode_size),
                KernelRecord.last.is_(ode_size), KernelRecord.mode.is_(mode),
                KernelRecord.sample_type.in_(sample_types)).all()
            return bool(len(data) == max_cores)
        # ODE size is not fixed.
        if sample_types is None:
            # Include data for all sample types.
            data = db_session.query(KernelRecord).filter(
                KernelRecord.kernel.is_(kernel), KernelRecord.machine.is_(machine), KernelRecord.compiler.is_(compiler),
                KernelRecord.method.is_(method), KernelRecord.ivp.is_(ivp), KernelRecord.frequency.is_(frequency),
                KernelRecord.cores.in_((c for c in range(1, max_cores + 1))), KernelRecord.mode.is_(mode)).order_by(
                KernelRecord.cores.asc(), KernelRecord.first.asc()).all()
        else:
            # Include only data for the specified sample types.
            data = db_session.query(KernelRecord).filter(
                KernelRecord.kernel.is_(kernel), KernelRecord.machine.is_(machine), KernelRecord.compiler.is_(compiler),
                KernelRecord.method.is_(method), KernelRecord.ivp.is_(ivp), KernelRecord.frequency.is_(frequency),
                KernelRecord.cores.in_((c for c in range(1, max_cores + 1))), KernelRecord.mode.is_(mode),
                KernelRecord.sample_type.in_(sample_types)).order_by(KernelRecord.cores.asc(),
                                                                     KernelRecord.first.asc(), ).all()
        if not data:
            return False
        # Check if data for all possible ODE sizes is available.
        cur_core = 0
        cur_last = sys_maxsize
        for record in data:
            if record.cores != cur_core:
                if cur_last != sys_maxsize:
                    return False
                cur_core = record.cores
                if record.first != 1:
                    return False
                cur_last = record.last
            else:
                if (cur_last + 1) != record.first:
                    return False
                cur_last = record.last
        assert cur_core == max_cores
        return bool(cur_last == sys_maxsize)


@attr.s(hash=True)
class PModelKernel:
    """Representation of a PModelKernel object.

    Attributes:
    -----------
    name: str
        Name affix added to code file of this object.
    kernel: Kernel
        Reference to associated Kernel object.
    iterations: str
        Arithmetic expression that describes this object's total iteration
        count.
    workingSets: list of str
        Working sets of this object.
    code_path: Path
        Relative path to generated code file.
    db_id: int
       ID of associated PModelKernel database table record.
    """
    kernel = attr.ib(type=Kernel)
    iterations = attr.ib(type=str, hash=False)
    workingSets = attr.ib(type=List[str], hash=False)
    name = attr.ib(type=str, default='')
    working_sets_serial = attr.ib(type=str, init=False, hash=False)
    kernel_db_id = attr.ib(type=int, init=False, hash=False)
    code_path = attr.ib(type=Path, init=False, hash=False)
    db_id = attr.ib(type=int, init=False)

    # Database information.
    db_table = Table('pmodel_kernel', METADATA,
                     Column('db_id', Integer, primary_key=True),
                     Column('name', String),
                     Column('kernel_db_id', Integer, ForeignKey('kernel.db_id')),
                     Column('iterations', String),
                     Column('working_sets_serial', String),
                     Column('updatedIn', String, default=__version__),
                     Column('updatedOn', DateTime, default=datetime.now, onupdate=datetime.now),
                     Column('updatedBy', String, default=getuser(), onupdate=getuser()),
                     UniqueConstraint('name', 'kernel_db_id'),
                     sqlite_autoincrement=True)

    @classmethod
    def from_yaml(cls, yaml: Dict, kernel: Kernel, pmodel_node: CodeNode, idx: Optional[int] = None) -> 'PModelKernel':
        """Construct PModelKernel object from YAML definition.

        Parameters:
        -----------
        yaml: dict
            YAML object describing this object.
        kernel: Kernel
            Reference to associated Kernel object.
        pmodel_node: CodeNode
            Reference to associated CodeNode object.
        name: str
            Name of this object.

        Returns:
        --------
        PModelKernel
            Created PModelObject object.
        """
        # Attribute iterations.
        iterations: str
        # Compute iteration count.
        if kernel.template.modelTool == ModelToolType.KERNCRAFT:
            if pmodel_node.type == CodeNodeType.PMODEL:
                iterations = CodeTree.iteration_count(pmodel_node.next)
                #
                prev = pmodel_node.prev
                assert prev
                while True:
                    if prev.prev is None:
                        break
                    prev = prev.prev
                if pmodel_node.prev.parent is not None:
                    iterations = CodeTree.iteration_count_up_to_root(pmodel_node.prev.parent, iterations)
            else:
                iterations = CodeTree.iteration_count(pmodel_node)
        elif kernel.template.modelTool == ModelToolType.YASKSITE:
            iterations = CodeTree.iteration_count(pmodel_node)
            # Yasksite code contains an implicit outer loop over the ODE system size 'n'.
            iterations = '{} * n'.format(iterations) if iterations != '' else 'n'
        else:
            assert False
        iterations = eval_math_expr(iterations, cast_to=str)

        # Attribute working_sets.
        # ... first define helper function that adds path to the root node.
        def _nodes_to_root(cur_root):
            assert cur_root is not None
            cur = cur_root
            if cur.prev is not None:
                cur = cur.prev
                assert cur
                while True:
                    if cur.prev is None:
                        break
                    cur = cur.prev
            new_root = deepcopy(cur.parent)
            if new_root is None:
                if cur_root.type is not CodeNodeType.ROOT:
                    new_root = RootNode(0, CodeNodeType.ROOT, None, None, None, cur_root, 'new_root')
                    cur_root.parent = new_root
                    return new_root
                return cur_root
            new_root.next = None
            new_root.child = cur_root
            cur_root.parent = new_root
            return _nodes_to_root(new_root)

        # ...
        name: str = '' if idx is None else str(idx + 1)
        idx: int = 0 if idx is None else idx
        #  ... compute working sets next.
        try:
            # We only support automatic working set derivation for kerncraft kernels at this point.
            # Otherwise, raise an index error exception to continue.
            if kernel.template.modelTool != ModelToolType.KERNCRAFT:
                raise IndexError('Warning: Automatic derivation of working sets is only supported for kernels that'
                                 'use model tool \'{}\'. Next, try to read working sets from the kernel\'s YAML'
                                 'description file...')
            # Prepend nodes on path from pmodel node to root node to include all required loops.
            root = deepcopy(pmodel_node)
            root = _nodes_to_root(root)
            # Compute working sets.
            cr = CountReferences()
            ws = cr.compute_working_sets(root, kernel.template.codegen)
            working_sets: List[str] = list(set(ws))
        except IndexError as err_idx:
            try:
                working_sets: List[str] = [ws for ws in yaml['pmodel'][idx]['working sets']]
            except KeyError:
                print(err_idx)
                raise RuntimeError('Error: Unable to load working set from the kernel\'s YAML description: {}.'.format(
                    kernel.name))
        # Create object.
        return cls(kernel, iterations, working_sets, name)

    @classmethod
    def from_database(cls, db_session: Session, pmodel_name: str, kernel_id: int) -> 'PModelKernel':
        """Construct PModelKernel object from database record.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        pmodel_name: str
            Name of the pmodel kernel, which is used as primary key in the database.
        kernel_id: IVP
            Database ID of the associated Kernels object.

        Returns:
        --------
        Kernel
            Created PModelKernel object.
        """
        try:
            pmodel_kernel: PModelKernel = db_session.query(PModelKernel).filter(
                PModelKernel.name.like(pmodel_name), PModelKernel.kernel_db_id.like(kernel_id)).one()
        except NoResultFound:
            raise RuntimeError('Unable to load PModelKernel object from database!')
        except MultipleResultsFound:
            raise RuntimeError('Unable to load PModelKernel object from database!')
        # Attribute working_sets.
        pmodel_kernel.workingSets = [workset for workset in pmodel_kernel.working_sets_serial.split(',') if workset]
        return pmodel_kernel

    def to_database(self, db_session: Session) -> 'PModelKernel':
        """Push this pmodel kernel object to the database.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        PModelKernel
            Instance of this object connected to database session.
        """
        # Check if database already contains this PModelKernel object.
        pmodel: PModelKernel = db_session.query(PModelKernel).filter(
            PModelKernel.name.like(self.name), PModelKernel.kernel_db_id.is_(self.kernel_db_id)).one_or_none()
        if pmodel:
            # Supplement attributes not saved in database.
            pmodel.workingSets = self.workingSets
            pmodel.kernel = self.kernel
            return pmodel
        # Add new object to database.
        # Attribute working_sets_serial.
        self.working_sets_serial = ','.join(self.workingSets)
        # Attribute kernel_db_id.
        self.kernel_db_id = self.kernel.db_id
        # Insert PModelKernel object.
        insert(db_session, self)
        return self

    def iteration_count(self, method: Optional[ODEMethod] = None, ivp: Optional[IVP] = None):
        """
        Determine number of iterations performed when executing this pmodel kernel object for a given configuration of
        ODE method and IVP and IVP system size 'n'.

        Parameters:
        -----------
        method: ODE method
            Trained ODE method.
        ivp: IVP
            Trained IVP.

        Returns:
        --------
        int
            Iteration count of this pmodel kernel for the given configuration.
        """
        constants = list()
        if method is not None:
            constants.extend([corrector_steps(method), stages(method)])
        if ivp is not None:
            constants.append(ivp_grid_size(ivp.gridSize))
        return eval_math_expr(self.iterations, constants)

    def determine_samples_from_working_sets(self, machine: MachineState, method: Optional[ODEMethod] = None,
                                            ivp: Optional[IVP] = None) -> List[SampleInterval]:
        """
        Determine a list of significant IVP system size sample intervals and points for this object. These samples are
        then used to obtain performance predictions for this object. The shape of the list depends on the objects
        working set sizes and the characteristics of the given run configuration of machine, compiler, ODE method and
        IVP.

        Parameters:
        -----------
        machine: MachineState
            Trained machine.
        method: ODE method
            Trained ODE method.
        ivp: IVP
            Trained IVP.

        Returns:
        --------
        list of SampleInterval
            List of sample intervals.
        """
        config: Config = offsite.config.offsiteConfig
        sample_type: SampleType
        if config.args.mode == ProgramModeType.MODEL:
            sample_type = SampleType.MODEL_INNER
        elif config.args.mode == ProgramModeType.RUN:
            sample_type = SampleType.BENCH_INNER
        else:
            assert False
        assert config.memory_lvl_sample_offset >= 1.0
        # Raise error if no working sets were defined.
        if not self.workingSets:
            raise RuntimeError('No working sets defined for PModelKernel {}{}!'.format(self.kernel.name, self.name))
        # Create SampleInterval objects.
        samples = list()
        wset_start = 1
        wset_end = 1
        # Add samples for each working set.
        for wset_end in self.determine_max_ns_of_working_sets(machine, method):
            # Add samples in border region to previous working set.
            samples_border, start = create_samples_lower_border_working_set(
                wset_start, wset_end, config.samples_per_border, ivp)
            samples.extend(samples_border)
            # Switch to next working set once the end of the current set was reached.
            if start > wset_end:
                break
            # Add samples in border region to next working set.
            samples_border, end = create_samples_upper_border_working_set(
                start, wset_end, config.samples_per_border, ivp)
            samples.extend(samples_border)
            # Switch to next working set once the end of the last interval from the lower border region was reached.
            if end + 1 < start:
                break
            #
            size = end - start + 1
            step = int(size / config.samples_per_interval)
            remain = size % config.samples_per_interval
            #
            if size >= 10 * config.samples_per_interval:
                nothing_added = True
                my_start = start
                my_end = my_start + step
                for i in range(1, config.samples_per_interval + 1):
                    my_end = min(my_end + 1 + i if remain < i else my_end, end)
                    #
                    point = SampleInterval(my_start, my_end, sample_type).median(ivp)
                    if point:
                        nothing_added = False
                        #
                        samples.append(SampleInterval(my_start, my_end, sample_type, point))
                        # Switch to next.
                        my_start = my_end + 1
                        my_end = my_start + step
                    else:
                        # Unable to determine sample point for current interval. Hence, we increase the end value.
                        my_end = my_end + step
                if nothing_added:
                    # Unable to determine any sample point at all for the given interval! Hence, instead incorporate
                    # current working set in next working set.
                    wset_start = start
                else:
                    # Switch to start of next working set.
                    wset_start = wset_end + 1
            else:
                # Add median value of remaining range - the one that wasn't yet considered in the SampleInterval
                # objects of the transition area of the working set as sample.
                point = SampleInterval(start, end, sample_type).median(ivp)
                if point:
                    samples.append(SampleInterval(start, end, sample_type, point))
                    # Switch to start of next working set.
                    wset_start = wset_end + 1
                else:
                    # Unable to determine sample point for the given interval! Hence, instead incorporate current
                    # working set in next working set.
                    wset_start = start
        # Add samples in border region of data from memory and largest cache working set size.
        samples_border, start = create_samples_border_memory_lvl(
            wset_end + 1, sys_maxsize, config.samples_border_region_memory_lvl, ivp)
        samples.extend(samples_border)
        # Add samples in memory region.
        samples.extend(
            create_samples_memory_lvl(start, config.samples_memory_lvl, config.memory_lvl_sample_offset, ivp))
        # Sort samples by start value.
        for sample in samples:
            assert (sample.first <= sample.last)
        samples.sort(key=lambda x: x.first)
        return samples

    def determine_samples_for_yasksite_mode(
            self, machine: MachineState, ivp: Optional[IVP] = None) -> List[SampleInterval]:
        """
        Determine a list of significant IVP system size sample intervals and points for this object. These samples are
        then used to obtain performance predictions for this object.

        Parameters:
        -----------
        machine: MachineState
            Trained machine.
        ivp: IVP
            Trained IVP.

        Returns:
        --------
        list of SampleInterval
            List of sample intervals.
        """
        config: Config = offsite.config.offsiteConfig
        assert config.memory_lvl_sample_offset >= 1.0
        # Create SampleInterval objects.
        samples = list()
        # Define cut-off sample size and sample multiplier.
        cut_off_n = 50
        cut_off_n = eval_math_expr(ivp.ivp_size(), [ivp_grid_size(cut_off_n)], cast_to=int)
        multiplier = 1.1 ** ivp.characteristic.stencil_dim
        # Compute layer of the stencil.
        stencil_layer = 2 * ivp.characteristic.stencil_radius + 2
        # Determine all samples.
        end_size = int(2 * config.memory_lvl_sample_offset * ((machine.l3CacheElements / stencil_layer) ** (
                ivp.characteristic.stencil_dim / (ivp.characteristic.stencil_dim - 1))))
        cur_size = cut_off_n
        prev_size = 0
        while cur_size < end_size:
            point = SampleInterval(prev_size + 1, cur_size, SampleType.MODEL_INNER).median(ivp)
            if point:
                samples.append(SampleInterval(prev_size + 1, cur_size, SampleType.MODEL_INNER, point))
                # Switch to start of next sample.
                prev_size = cur_size
                cur_size = int(cur_size * multiplier)
            else:
                # Unable to determine sample point for the given interval! Hence, instead incorporate current
                # sample into next sample.
                cur_size = int(cur_size * multiplier)
        # Sort samples by start value.
        samples.sort(key=lambda x: x.first)
        return samples

    def determine_max_ns_of_working_sets(self, machine: MachineState, method: ODEMethod) -> List[int]:
        """
        Determine the maximum IVP size 'n' that still fit into the cache level and working sets given by the run
        configuration and this object.

        Parameters:
        -----------
        machine: MachineState
            Trained machine.
        method: ODE method
            Trained ODE method.

        Returns:
        --------
        list of int
            Maximal system sizes fitting into the particular cache levels and working sets.
        """
        # For each cache level and working set permutation, determine the maximum IVP system size 'n' that still fits
        # into that level and working set.
        max_ns: List[int] = list()
        for cache in [machine.l1CacheElements, machine.l2CacheElements, machine.l3CacheElements]:
            max_ns.extend(self.determine_max_n_of_working_sets_for_cache_lvl(method, cache))
        # Sort working set ends by size to iterate list in order when creating the SampleInterval objects.
        max_ns.sort()
        return max_ns

    def determine_max_n_of_working_sets_for_cache_lvl(self, method: ODEMethod, cache_elements: int):
        """
        Determine the maximum system size 'n' that still fit into the particular given cache level and working sets
        given by this object.

        Parameters:
        -----------
        machine: MachineState
            Used machine.
        method: ODE method
            Used ODE method.
        cache_elements: int
            Maximum number of elements fitting into the particular cache level considered.

        Returns:
        --------
        list of int
            Maximal system sizes fitting into the particular cache levels and working sets.
        """
        avail_cache_size: float = offsite.config.offsiteConfig.available_cache_size
        assert avail_cache_size <= 1.0
        # Define some constants that might be part of the arithmetical working set expressions.
        constants = list()
        if method is not None:
            constants.extend([corrector_steps(method), stages(method)])
        # Determine the maximum system size 'n' that still fits into that level and working set.
        max_ns = list()
        for wset in self.workingSets:
            if 'n' not in wset:
                continue
            # Solve equation 'wset = avail_cache * cache'.
            solution: str = solve_equation(wset, '{}*{}'.format(avail_cache_size, cache_elements), 'n', constants)
            # Get end from the solution of the equation
            max_ns.append(int(solution[0]))
        return max_ns


@attr.s
class PModelRecord:
    """Representation of a pmodel kernel ECM prediction table database record.

    Attributes:
    -----------
    pmodel_kernel: int
        Used pmodel kernel.
    machine: int
        Used machine.
    compiler: int
        Used compiler.
    method: int
        Used ODE method.
    ivp: int
        Used IVP.
    cores: int
        Used number of CPU cores.
    sample: SampleInterval
        Used sample interval.
    prediction: float
        Kernel runtime prediction.
    mode: str
        Prediction obtained in IACA or OSACA mode.
    db_id: int
        ID of associated pmodel kernel ECM prediction database table record.
    """
    pmodel_kernel = attr.ib(type=int)
    machine = attr.ib(type=int)
    compiler = attr.ib(type=int)
    method = attr.ib(type=int)
    ivp = attr.ib(type=int)
    cores = attr.ib(type=int)
    sample = attr.ib(type=SampleInterval)
    prediction = attr.ib(type=str)
    mode = attr.ib(type=str)
    first = attr.ib(type=int, init=False)
    last = attr.ib(type=int, init=False)
    db_id = attr.ib(type=int, init=False)

    # Database information.
    db_table = Table('pmodel_kernel_prediction', METADATA,
                     Column('db_id', Integer, primary_key=True),
                     Column('pmodel_kernel', Integer, ForeignKey('pmodel_kernel.db_id')),
                     Column('machine', Integer, ForeignKey('machine.db_id')),
                     Column('compiler', Integer, ForeignKey('compiler.db_id')),
                     Column('method', Integer, ForeignKey('ode_method.db_id')),
                     Column('ivp', Integer, ForeignKey('ivp.db_id')),
                     Column('cores', Integer),
                     Column('first', Integer),
                     Column('last', Integer),
                     Column('prediction', String),
                     Column('mode', String),
                     Column('updatedIn', String, default=__version__),
                     Column('updatedOn', DateTime, default=datetime.now, onupdate=datetime.now),
                     Column('updatedBy', String, default=getuser(), onupdate=getuser()),
                     sqlite_autoincrement=True)

    def __attrs_post_init__(self):
        """Create this PModelRecord object.

        Parameters:
        -----------
        -

        Returns:
        --------
        -
        """
        self.first = self.sample.first
        self.last = self.sample.last

    def to_database(self, db_session: Session):
        """Push this PModelRecord object to the database.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        -
        """
        # Attribute prediction.
        self.prediction = str(self.prediction)
        # Insert KernelRecord object.
        insert(db_session, self)

    @staticmethod
    def update(db_session: Session, pmodel_kernel: int, machine: int, compiler: int, method: int, ivp: int,
               interval: SampleInterval, records: List[Tuple[SampleInterval, str]], mode: str):
        """Update data records in PModelRecord table with new data.

        Parameters:
        -----------
        db_session: sqlalchemy.orm.session.Session
            Used database session.
        pmodel_kernel: int
            Used pmodel kernel.
        machine: int
            Used machine.
        compiler: int
            Used compiler.
        method: int
            Used ODE method.
        ivp: int
            Used IVP.
        interval: SampleInterval
            Used sample interval.
        records: list of tuple (SampleInterval, ECM prediction)
            PModelKernel ECM prediction data.
        mode: str
            Incore model mode used to obtain data.

        Returns:
        --------
        -
        """
        for record in records.items():
            cores: int = record[0]
            prediction: str = record[1]
            # Select and update all already included data records.
            try:
                queried_record: PModelRecord = db_session.query(PModelRecord).filter(
                    PModelRecord.pmodel_kernel.is_(pmodel_kernel), PModelRecord.machine.is_(machine),
                    PModelRecord.compiler.is_(compiler), PModelRecord.method.is_(method), PModelRecord.ivp.is_(ivp),
                    PModelRecord.cores.is_(cores), PModelRecord.first.is_(interval.first),
                    PModelRecord.last.is_(interval.last), PModelRecord.mode.is_(mode)).one()
                # Update data record.
                queried_record.prediction = str(prediction)
            except NoResultFound:
                # Insert new data record.
                record = PModelRecord(pmodel_kernel, machine, compiler, method, ivp, cores, interval, prediction, mode)
                record.to_database(db_session)
            except MultipleResultsFound:
                raise RuntimeError('Unable to update PModelKernel record!')
