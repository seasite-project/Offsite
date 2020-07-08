"""@package ivp
Definitions of classes IVP and IVPCharacteristic.
"""

from datetime import datetime
from getpass import getuser
from pathlib import Path

import attr
from sqlalchemy import Boolean, Column, DateTime, Enum, ForeignKey, Integer, String, Table
from sqlalchemy.orm import Session
from sqlalchemy.orm.exc import NoResultFound, MultipleResultsFound

import offsite.config
from offsite import __version__
from offsite.config import ModelToolType
from offsite.db.db import METADATA
from offsite.db.db import insert
from offsite.descriptions.parser_utils import load_yaml, deserialize_obj, serialize_obj, ComponentDict, ConstantDict
from offsite.evaluation.math_utils import solve_equation


@attr.s
class IVP:
    """Representation of an initial value problem (IVP) object.

    Attributes:
    -----------
    name : str
        Name of this IVP object.
    gridSize : str
        Arithmetic expression that describes the IVP's grid size.
    characteristic : IVPCharacteristic
        Characteristics describing the shape of this IVP object.
    components : list of IVPComponent
        List of components connected to this IVP object.
    constants : list of IVPConstant
        List of constants connected to this IVP object.
    code_eval_range : str
        DSL representation of the IVP function eval_range.
    code_eval_component : str
        DSL representation of the IVP function eval_component.
    code_initial_values : str
        DSL representation of the IVP function initial_values.
    modelTool : ModelToolType
        Performance model tool used to obtain performance data of this object.
    db_id : int
        ID of associated IVP table record in the database.
    """
    name = attr.ib(type=str)
    gridSize = attr.ib(type=str)
    characteristic = attr.ib(type='IVPCharacteristic')
    components = attr.ib(type='ComponentDict')
    constants = attr.ib(type='ConstantDict')
    code_eval_range = attr.ib(type=str)
    code_eval_component = attr.ib(type=str)
    code_stencil_path = attr.ib(type=str)
    code_initial_values = attr.ib(type=str)
    code_required_indices = attr.ib(type=str)
    modelTool = attr.ib(type=ModelToolType)
    components_serial = attr.ib(type=str, init=False)
    constants_serial = attr.ib(type=str, init=False)
    code_eval_range_serial = attr.ib(type=str, init=False)
    code_eval_component_serial = attr.ib(type=str, init=False)
    code_initial_values_serial = attr.ib(type=str, init=False)
    code_required_indices_serial = attr.ib(type=str, init=False)
    db_id = attr.ib(type=int, init=False)

    # Database information.
    db_table = Table('ivp', METADATA,
                     Column('db_id', Integer, primary_key=True),
                     Column('name', String, unique=True),
                     Column('gridSize', String),
                     Column('modelTool', Enum(ModelToolType)),
                     Column('components_serial', String),
                     Column('constants_serial', String),
                     Column('code_eval_range_serial', String),
                     Column('code_eval_component_serial', String),
                     Column('code_stencil_path', String),
                     Column('code_initial_values_serial', String),
                     Column('code_required_indices_serial', String),
                     Column('updatedIn', String, default=__version__),
                     Column('updatedOn', DateTime, default=datetime.now, onupdate=datetime.now),
                     Column('updatedBy', String, default=getuser(), onupdate=getuser()),
                     sqlite_autoincrement=True)

    @classmethod
    def from_yaml(cls, yaml_path: Path) -> 'IVP':
        """Construct IVP object from YAML definition.

        Parameters:
        -----------
        yaml_path : Path
            Relative path to this object's YAML file.

        Returns:
        --------
        IVP
            Created IVP object.
        """
        config = offsite.config.offsiteConfig
        # Load YAML data.
        yaml = load_yaml(yaml_path)
        # Attribute name.
        name = yaml_path.stem
        # Attribute grid_size.
        grid_size = yaml['ODE grid dim']
        # Attribute model_tool.
        tool_str = yaml['model_tool']
        if tool_str == ModelToolType.KERNCRAFT.value:
            model_tool = ModelToolType.KERNCRAFT
        elif tool_str == ModelToolType.YASKSITE.value:
            model_tool = ModelToolType.YASKSITE
        else:
            assert False
        # Attribute characteristic.
        # .. used to temporarily store yaml data for post init.
        characteristic = yaml['characteristics']
        # Attribute components.
        components = ComponentDict.from_data(yaml['components'])
        # Attribute constants.
        constants = ConstantDict.from_data(yaml['constants'])
        if model_tool == ModelToolType.KERNCRAFT:
            # Attribute code_eval_range.
            code_eval_range = yaml['codegen']['eval_range']
            # Attribute code_eval_component.
            code_eval_component = yaml['codegen']['eval_component']
            # Attribute code_required_indices.
            code_required_indices = yaml['codegen']['required_indices']
            # Add empty string for unused codegen properties.
            code_stencil_path = ''
        elif model_tool == ModelToolType.YASKSITE:
            code_stencil_path = yaml['codegen']['stencil']
            # Check if path exists.
            if 'YASKSITE_STENCIL_DIR' in code_stencil_path:
                resolved_path = Path(code_stencil_path.replace('YASKSITE_STENCIL_DIR', config.yasksite_stencil_dir))
            if not resolved_path.exists():
                raise RuntimeError('Stencil file of IVP \'{}\' not found: \'{}\''.format(name, resolved_path))
            # Add empty string for unused codegen properties.
            code_eval_range = ''
            code_eval_component = ''
            code_required_indices = ''
        # Attribute code_initial_values.
        code_initial_values = yaml['codegen']['initial_values']
        # Create object.
        return cls(name, grid_size, characteristic, components, constants, code_eval_range, code_eval_component,
                   code_stencil_path, code_initial_values, code_required_indices, model_tool)

    @classmethod
    def from_database(cls, db_session: Session, ivp_name: str) -> 'IVP':
        """Construct IVP object from database record.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.
        ivp_name : str
            Name of the IVP, which is used as primary key in the database.

        Returns:
        --------
        IVP
            Created IVP object.
        """
        try:
            ivp = db_session.query(IVP).filter(IVP.name.like(ivp_name)).one()
        except NoResultFound:
            raise RuntimeError('Unable to load IVP object from database!')
        except MultipleResultsFound:
            raise RuntimeError('Unable to load IVP object from database!')
        # Attribute ivp_components.
        ivp.components = deserialize_obj(ivp.components_serial)
        # Attribute ivp_characteristic.
        try:
            ivp.characteristic = db_session.query(IVPCharacteristic).filter(
                IVPCharacteristic.ivp_db_id.like(ivp.db_id)).one()
        except NoResultFound:
            raise RuntimeError('Unable to load IVPCharacteristic object from database!')
        except MultipleResultsFound:
            raise RuntimeError('Unable to load IVPComponent objects from database!')
        # Attribute constants.
        ivp.constants = deserialize_obj(ivp.constants_serial)
        # Attribute code_eval_range.
        ivp.code_eval_range = deserialize_obj(ivp.code_eval_range_serial)
        # Attribute code_eval_component.
        ivp.code_eval_component = deserialize_obj(ivp.code_eval_component_serial)
        # Attribute code_initial_values.
        ivp.code_initial_values = deserialize_obj(ivp.code_initial_values_serial)
        # Attribute code_initial_values.
        ivp.code_required_indices = deserialize_obj(ivp.code_required_indices_serial)
        return ivp

    def __attrs_post_init__(self):
        """Create the IVP objects associated with this object.

        Parameters:
        -----------
        -

        Returns:
        --------
        -
        """
        yaml = self.characteristic
        self.characteristic = IVPCharacteristic.from_yaml(yaml, self)

    def to_database(self, db_session: Session):
        """Push this IVP object to the database.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        IVP
            Instance of this object connected to database session.
        """
        # Check if database already contains this IVP object.
        ivp = db_session.query(IVP).filter(IVP.name.like(self.name)).one_or_none()
        if ivp:
            # Supplement attributes not saved in database.
            # IVPComponent objects.
            ivp.components = self.components
            # IVPCharacteristic object.
            self.characteristic.ivp_db_id = ivp.db_id
            ivp.characteristic = self.characteristic.to_database(db_session)
            # IVPConstant objects.
            ivp.constants = self.constants
            # Eval_range function.
            ivp.code_eval_range = self.code_eval_range
            # Eval_component function.
            ivp.code_eval_component = self.code_eval_component
            # Initial_value function.
            ivp.code_initial_values = self.code_initial_values
            # Required_indices function.
            ivp.code_required_indices = self.code_required_indices
            # Update serialized members.
            ivp.components_serial = serialize_obj(ivp.components)
            ivp.constants_serial = serialize_obj(ivp.constants)
            ivp.code_eval_range_serial = serialize_obj(ivp.code_eval_range)
            ivp.code_eval_component_serial = serialize_obj(ivp.code_eval_component_serial)
            ivp.code_initial_values_serial = serialize_obj(ivp.code_initial_values)
            ivp.code_required_indices_serial = serialize_obj(ivp.code_required_indices)
            return ivp
        # Add new object to database.
        # Serialize data.
        self.components_serial = self.components.serialize()
        self.constants_serial = self.constants.serialize()
        self.code_eval_range_serial = serialize_obj(self.code_eval_range)
        self.code_eval_component_serial = serialize_obj(self.code_eval_component)
        self.code_initial_values_serial = serialize_obj(self.code_initial_values)
        self.code_required_indices_serial = serialize_obj(self.code_required_indices)
        # IVP object.
        insert(db_session, self)
        # IVPCharacteristic object.
        self.characteristic.ivp_db_id = self.db_id
        self.characteristic.to_database(db_session)
        return self

    def ivp_size(self) -> str:
        """IVP size as a function of the IVP's grid size.

        Parameters:
        -----------
        -

        Returns:
        --------
        - str
            Arithmetic expression that describes the IVP's size.
        """
        return str(solve_equation('g', self.gridSize, 'n')[0])

    @staticmethod
    def database_id(db_session: Session, yaml_path: str) -> int:
        """Return database ID ob an IVP object given by its YAML description.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.
        yaml_path : str
            Relative path to this object's YAML file.

        Returns:
        --------
        int
            Database ID of IVP or -1 if not in database.
        """
        # Parse IVP yaml description.
        ivp = IVP.from_yaml(yaml_path)
        # Query IVP in database.
        ivp = db_session.query(IVP).filter(IVP.name.like(ivp.name)).one_or_none()
        # Return ID of IVP object.
        return ivp.db_id if ivp else -1


class IVPComputationType(Enum):
    """Characterizes behaviour/boundaries of the computation(s) of an IVP.
    """
    UNKNOWN = 'unknown'
    MIX = 'mixed-bound'
    COMPUTE = 'compute-bound'
    MEMORY = 'memory-bound'


@attr.s
class IVPCharacteristic:
    """Representation of an IVPCharacteristic object.

    Attributes:
    -----------
    ivp : IVP
        Reference to associated IVP object.
    isSparse : bool
        True if referenced IVP is a sparse problem.
    isStencil : bool
        True if referenced IVP is a stencil-like problem.
    computationType : IVPComputationType
        Computation
    db_id : int
         ID of associated IVPCharacteristic table record in the database.
     """
    ivp = attr.ib(type=IVP)
    isSparse = attr.ib(type=bool)
    isStencil = attr.ib(type=bool)
    access_distance = attr.ib(type=str)
    computationType = attr.ib(type=IVPComputationType, default=IVPComputationType.UNKNOWN)
    stencil_radius = attr.ib(type=Integer, default=None)
    stencil_dim = attr.ib(type=Integer, default=None)
    ivp_db_id = attr.ib(type=Integer, init=False)
    db_id = attr.ib(type=int, init=False)

    # Database information.
    db_table = Table('ivp_characteristic', METADATA,
                     Column('db_id', Integer, primary_key=True),
                     Column('ivp_db_id', Integer, ForeignKey('ivp.db_id'), unique=True),
                     Column('isSparse', Boolean),
                     Column('isStencil', Boolean),
                     Column('computationType', String),
                     Column('access_distance', String),
                     Column('stencil_radius', Integer),
                     Column('stencil_dim', Integer),
                     Column('updatedIn', String, default=__version__),
                     Column('updatedOn', DateTime, default=datetime.now, onupdate=datetime.now),
                     Column('updatedBy', String, default=getuser(), onupdate=getuser()),
                     sqlite_autoincrement=True)

    @classmethod
    def from_yaml(cls, yaml: str, ivp: IVP) -> 'IVPCharacteristic':
        """Construct IVPCharacteristic object from YAML definition.

        Parameters:
        -----------
        yaml : str
            YAML object describing this object.
        ivp : IVP
            Reference to associated IVP object.

        Returns:
        --------
        IVPCharacteristic
            Created IVPCharacteristic object.
        """
        # Attribute is_sparse.
        is_sparse = 'sparse' in yaml
        # Attribute is_stencil.
        is_stencil = 'stencil' in yaml
        # Attribute computation_type.
        # set proper value...
        ctype = IVPComputationType.UNKNOWN
        for component in ivp.components.values():
            if ctype is IVPComputationType.UNKNOWN:
                ctype = component.component_type
            elif component.component_type == IVPComputationType.MIX or ctype != component.component_type:
                ctype = 'mixed-bound'
                break
        # Attribute access distance.
        access_distance = None
        for entry in yaml:
            if isinstance(entry, dict):
                if 'access_distance' in entry:
                    access_distance = str(entry['access_distance'])
                    break
        assert access_distance is not None
        # Stencil-specific attributes.
        stencil_radius = None
        stencil_dim = None
        if is_stencil:
            for entry in yaml:
                if isinstance(entry, dict):
                    # Attribute stencil_dim.
                    if 'stencil_dim' in entry:
                        stencil_dim = int(entry['stencil_dim'])
                    # Attribute stencil_radius.
                    elif 'stencil_radius' in entry:
                        stencil_radius = int(entry['stencil_radius'])
            if not stencil_dim:
                raise RuntimeError('Stencil-like IVP \'{}\' requires attribute \'stencil_dim\''.format(ivp.name))
            if not stencil_radius:
                raise RuntimeError('Stencil-like IVP \'{}\' requires attribute \'stencil_radius\''.format(ivp.name))
        # Create object.
        return cls(ivp, is_sparse, is_stencil, access_distance, ctype, stencil_radius, stencil_dim)

    @classmethod
    def from_database(cls, db_session: Session, ivp_id: int) -> 'IVPCharacteristic':
        """Construct IVPCharacteristic object from database record.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.
        ivp_id : IVP
            Database ID of the associated IVP object.

        Returns:
        --------
        IVPCharacteristic
            Created IVPCharacteristic object.
        """
        try:
            ivp_characteristic = db_session.query(IVPCharacteristic).filter(
                IVPCharacteristic.ivp_db_id.like(ivp_id)).one()
        except NoResultFound:
            raise RuntimeError('Unable to load IVPCharacteristic object from database!')
        except MultipleResultsFound:
            raise RuntimeError('Unable to load IVPCharacteristic object from database!')
        return ivp_characteristic

    def to_database(self, db_session: Session):
        """Push this IVPCharacteristic object to the database.

        Parameters:
        -----------
        db_session : sqlalchemy.orm.session.Session
            Used database session.

        Returns:
        --------
        IVPCharacteristic
            Instance of this object connected to database session.
        """
        # Check if database already contains this IVPCharacteristic object.
        characteristic = db_session.query(IVPCharacteristic).filter(
            IVPCharacteristic.ivp_db_id.is_(self.ivp_db_id)).one_or_none()
        if characteristic:
            # Supplement attributes not saved in database.
            characteristic.ivp = self.ivp
            return characteristic
        # Add new object to database.
        # Attribute ivp_db_id.
        self.ivp_db_id = self.ivp.db_id
        # Insert IVPCharacteristic object.
        insert(db_session, self)
        return self
