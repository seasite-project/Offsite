"""@package descriptions.parser
Functions to parse YAML descriptions as well as different utility functions and classes used during YAML parsing.

@author: Johannes Seiferth
"""

from codecs import decode, encode
from enum import IntEnum
from pickle import dumps, loads
from re import findall, split as re_split, S
from typing import Any, Dict, List, Tuple

import attr
from pathlib2 import Path
from ruamel import yaml


def load_yaml(path_to_yaml: Path) -> Dict:
    """Load a YAML file.

    Parameters:
    -----------
    path_to_yaml: Path
        Path to loaded YAML file.

    Returns:
    --------
    dict
        Dictionary containing the loaded YAML data.
    """
    with path_to_yaml.open('r') as file_handle:
        return yaml.load(file_handle, Loader=yaml.Loader)


@attr.s
class ComponentDesc:
    """Representation of a single component of an IVP.
    """
    component_type = attr.ib(type=str)
    code = attr.ib(type=str)
    equations = attr.ib(type=str)
    first = attr.ib(type=str)

    @classmethod
    def from_yaml(cls, yaml_data: Dict[str, str]) -> 'ComponentDesc':
        """Construct ComponentDesc object from YAML definition.

        Parameters:
        -----------
        yaml_data: dict
            YAML object describing this object.

        Returns:
        --------
        ComponentDesc
            Created ComponentDesc object.
        """
        # Attribute component_type.
        component_type = yaml_data['type']
        # Attribute code.
        code = yaml_data['code']
        # Attribute equations.
        equations = yaml_data['size']
        # Attribute first.
        first = yaml_data['first']
        # Create object.
        return cls(component_type, code, equations, first)


@attr.s
class ComponentDict(dict):
    """Dictionary to store the components of a kernel.
    """

    @classmethod
    def from_data(cls, components: Dict[str, str]) -> 'ComponentDict':
        """Construct initialized ComponentDict object.

        Parameters:
        -----------
        components : dict
            Dict of component descriptions to be stored.

        Returns:
        --------
        ComponentDict
            Created ComponentDict object.
        """
        dict_obj = cls.__new__(cls)
        # Insert components.
        for component in components:
            name = component['name'] if 'name' in component else ''
            if dict_obj.__contains__(name):  # Ignore duplicate keys.
                pass
            else:
                dict_obj[name] = ComponentDesc.from_yaml(component)
        return dict_obj

    def get_component(self, name: str) -> ComponentDesc:
        """Return description of a computation in the dictionary.

        Searches for component 'name' in the dictionary and returns its component description.

        Parameters:
        -----------
        name: str
            Identifier of the component requested.

        Returns:
        --------
        ComponentDesc
            Description of the component requested.
        """
        try:
            return next(filter(lambda x: x[0] == name, self.items()))[1]
        except StopIteration:
            raise RuntimeError('Component \'{}\' not in dict!'.format(name))

    def serialize(self) -> str:
        """Serialize dictionary data in base64 encoding.

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Pickled representation of the dictionary data in base64 encoding.
        """
        return serialize_obj(self)


@attr.s
class ComputationDesc:
    """Representation of a single computation of a kernel.
    """
    computation = attr.ib(type=str)


@attr.s
class ComputationDict(dict):
    """Dictionary to store the computations of a kernel.
    """

    @classmethod
    def from_data(cls, computations: Dict[str, str]) -> 'ComputationDict':
        """Construct initialized ComputationDict object.

        Parameters:
        -----------
        computations: dict
            Dict of computation descriptions to be stored.

        Returns:
        --------
        ComputationDict
            Created ComputationDict object.
        """
        dict_obj = cls.__new__(cls)
        # Insert computations.
        for name, computation in computations.items():
            if dict_obj.__contains__(name):  # Ignore duplicate keys.
                pass
            else:
                dict_obj[name] = ComputationDesc(computation)
        return dict_obj

    def get_computation(self, name: str) -> ComputationDesc:
        """Return description of a computation in the dictionary.

        Searches for computation 'name' in the dictionary and returns its computation description.

        Parameters:
        -----------
        name: str
            Identifier of the computation requested.

        Returns:
        --------
        ComputationDesc
            Description of the computation requested.
        """
        try:
            return next(filter(lambda x: x[0] == name, self.items()))[1]
        except StopIteration:
            raise RuntimeError(
                'Computation \'{}\' not in dict!'.format(name))

    def serialize(self) -> str:
        """Serialize dictionary data in base64 encoding.

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Pickled representation of the dictionary data in base64 encoding.
        """
        return serialize_obj(self)


@attr.s
class ConstantDesc:
    """Representation of a single constant of an IVP."""
    datatype = attr.ib(type=str)
    value = attr.ib(type=str)


@attr.s
class ConstantDict(dict):
    """Dictionary to store the constants of an IVP."""

    @classmethod
    def from_data(cls, constants: List[str]) -> 'ConstantDict':
        """Construct initialized ConstantDict object.

        Parameters:
        -----------
        constants: list of str
            List of constant descriptions to be stored.

        Returns:
        --------
        ConstantDict
            Created ConstantDict object.
        """
        dict_obj = cls.__new__(cls)
        # Insert datastruct.
        for constant in constants:
            constant = constant.strip()
            parts = re_split(r'=', constant)
            assert len(parts) == 2
            # ... value
            value = parts[1].strip()
            # Parse lhs of constant expression.
            lhs = parts[0].strip()
            parts = re_split(r'\s+', lhs)
            assert len(parts) == 2
            # ... datatype
            datatype = parts[0]
            # ... name
            name = parts[1]
            if dict_obj.__contains__(name):  # Ignore duplicate keys.
                pass
            else:
                dict_obj[name] = ConstantDesc(datatype, value)
        return dict_obj

    def serialize(self) -> str:
        """Serialize dictionary data in base64 encoding.

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Pickled representation of the dictionary data in base64 encoding.
        """
        return serialize_obj(self)

    def get_constant(self, name: str) -> ConstantDesc:
        """Return description of a constant in the dictionary.

        Searches for constant 'name' in the dictionary and returns its constant description.

        Parameters:
        -----------
        name: str
            Identifier of the constant requested.

        Returns:
        --------
        ConstantDesc
            Description of the constant requested.
        """
        try:
            return next(filter(lambda x: x[0] == name, self.items()))[1]
        except StopIteration:
            raise RuntimeError('Constant \'{}\' not in dict!'.format(name))

    def as_tuple(self) -> List[Tuple[str, str]]:
        """Return descriptions as pairs  of a constant in the dictionary.

        Parameters:
        -----------
        -

        Returns:
        --------
        List of Tuple of str
            Descriptions as tuples (name, value) pairs.
        """
        return [(key, desc.value) for key, desc in self.items()]


class DatastructType(IntEnum):
    """Datastruct types supported by the parser.
    """
    scalar = 0
    array1D = 1
    array2D = 2
    array3D = 3


@attr.s
class DatastructDesc:
    """Representation of a single datastruct of a kernel."""
    datatype = attr.ib(type=str)
    struct_type = attr.ib(type=DatastructType)
    size = attr.ib(type=List[str])
    isYasksiteParam = attr.ib(type=bool)

    @classmethod
    def from_yaml(cls, dtype: str, dimensions: str, is_yasksite_param: bool) -> 'DatastructDesc':
        """Construct DatastructDesc object from YAML definition.

        Parameters:
        -----------
        datatype: dict
            YAML object describing this object.

        Returns:
        --------
        DatastructDesc
            Created DatastructDesc object.
        """
        # Attribute datatype.
        datatype = dtype
        # Attribute type.
        struct_type = DatastructType(dimensions.count('['))
        # Attribute size.
        dimensions = findall(r'\w+', dimensions, S)
        if struct_type == DatastructType.scalar:
            size = '1'
        elif struct_type in (DatastructType.array1D, DatastructType.array2D, DatastructType.array3D):
            size = [str(x) for x in dimensions]
        else:
            assert False
        # Create object.
        return cls(datatype, struct_type, size, is_yasksite_param)

    def dimensions_to_string(self) -> str:
        """Convert dimension sizes to static C array dimension declaration.

        scalar --> do nothing
        1D array --> [size_x]
        2D array --> [size_x][size_y]
        3D array --> [size_x][size_y][size_z]

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Static array dimension declaration.
        """
        if self.struct_type == DatastructType.scalar:
            dim_to_str = ''
        elif self.struct_type == DatastructType.array1D:
            dim_to_str = '[{}]'.format(self.size[0])
        elif self.struct_type == DatastructType.array2D:
            dim_to_str = '[{}][{}]'.format(self.size[0], self.size[1])
        elif self.struct_type == DatastructType.array3D:
            dim_to_str = '[{}][{}][{}]'.format(
                self.size[0], self.size[1], self.size[2])
        else:
            assert False
        return dim_to_str

    def dimensions_to_string_yasksite(self) -> str:
        """Convert dimension sizes to static C array dimension declaration.

        scalar --> do nothing
        1D array --> [size_x]
        2D array --> [size_x][size_y]
        3D array --> [size_x][size_y][size_z]

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Static array dimension declaration.
        """
        if self.struct_type == DatastructType.scalar:
            dim_to_str = ''
        elif self.struct_type == DatastructType.array1D:
            if self.size[0] == 'n':
                dim_to_str = ''
            else:
                dim_to_str = '[{}]'.format(self.size[0])
        elif self.struct_type == DatastructType.array2D:
            if self.size[0] == 'n' and self.size[1] == 'n':
                dim_to_str = ''
            elif self.size[0] == 'n':
                dim_to_str = '[{}]'.format(self.size[1])
            elif self.size[1] == 'n':
                dim_to_str = '[{}]'.format(self.size[0])
            else:
                dim_to_str = '[{}][{}]'.format(self.size[0], self.size[1])
        elif self.struct_type == DatastructType.array3D:
            if self.size[0] == 'n' and self.size[1] == 'n' and self.size[2] == 'n':
                dim_to_str = ''
            elif self.size[0] == 'n' and self.size[1] == 'n':
                dim_to_str = '[{}]'.format(self.size[2])
            elif self.size[0] == 'n' and self.size[2] == 'n':
                dim_to_str = '[{}]'.format(self.size[1])
            elif self.size[1] == 'n' and self.size[2] == 'n':
                dim_to_str = '[{}]'.format(self.size[0])
            elif self.size[0] == 'n':
                dim_to_str = '[{}][{}]'.format(self.size[1], self.size[2])
            elif self.size[1] == 'n':
                dim_to_str = '[{}][{}]'.format(self.size[0], self.size[2])
            elif self.size[2] == 'n':
                dim_to_str = '[{}][{}]'.format(self.size[0], self.size[1])
            else:
                dim_to_str = '[{}][{}][{}]'.format(self.size[0], self.size[1], self.size[2])
        else:
            assert False
        return dim_to_str

    def dimensions_1d_to_string(self, idx_values: List[str] = None) -> str:
        """Convert dimension sizes to flattened 1D C array statement using the index values given.

        scalar --> do nothing
        1D array --> [idx_val_x]
        2D array --> [(idx_val_x*size_x) + idx_val_y]
        3D array --> [(idx_val_x*size_x) + (idx_val_y*size_y) + idx_val_z]

        Parameters:
        -----------
        idx_values : list of str
            Indices used to access the dimensions.

        Returns:
        --------
        str
            Array statement.
        """
        if self.struct_type == DatastructType.scalar:
            dim_to_str = ''
        elif self.struct_type == DatastructType.array1D:
            assert len(idx_values) == 1
            dim_to_str = '[{}]'.format(idx_values[0])
        elif self.struct_type == DatastructType.array2D:
            assert len(idx_values) == 2
            dim_to_str = '[{}+({}*{})]'.format(idx_values[1], idx_values[0], self.size[0])
        elif self.struct_type == DatastructType.array3D:
            assert len(idx_values) == 3
            dim_to_str = '[{}+({}*{})+({}*{})]'.format(
                idx_values[2], idx_values[1], self.size[1], idx_values[0], self.size[0])
        else:
            assert False
        return dim_to_str


@attr.s
class DatastructDict(dict):
    """Dictionary to store the datastructs of a kernel."""

    @classmethod
    def from_data(cls, datastructs: List[str]) -> 'DatastructDict':
        """Construct initialized DatastructDict object.

        Parameters:
        -----------
        datastruct: list of str
            List of datastruct descriptions to be stored.

        Returns:
        --------
        DatastructDict
            Created DatastructDict object.
        """
        dict_obj = cls.__new__(cls)
        # Insert datastruct.
        for datastruct in datastructs:
            datastruct = datastruct.strip()
            # Parse datastruct string using regular expressions...
            parts = re_split(r'\s+', datastruct)
            assert len(parts) in [2, 3]
            # ... check if datatype is a parameter in yasksite mode
            if '@' in parts[0]:
                is_yasksite_param = True
                parts[0] = parts[0].replace('@', '')
            else:
                is_yasksite_param = False
            # ... datatype
            datatype = parts[0].strip()
            if len(parts) == 3:
                # ... name
                name = parts[1].strip()
                # ... dimension
                size = parts[2].strip()
            else:
                parts = parts[1].strip()
                parts = re_split(r'\[', parts, maxsplit=1)
                # ... name
                name = parts[0].strip()
                # ... dimension
                size = '' if len(parts) == 1 else '[' + parts[1].strip()
            if dict_obj.__contains__(name):  # Ignore duplicate keys.
                pass
            else:
                dict_obj[name] = DatastructDesc.from_yaml(datatype, size, is_yasksite_param)
        return dict_obj

    def serialize(self) -> str:
        """Serialize dictionary data in base64 encoding.

        Parameters:
        -----------
        -

        Returns:
        --------
        str
            Pickled representation of the dictionary data in base64 encoding.
        """
        return serialize_obj(self)

    def get_datastruct(self, name: str) -> DatastructDesc:
        """Return description of a datastruct in the dictionary.

        Searches for datastruct 'name' in the dictionary and returns its datastruct description.

        Parameters:
        -----------
        name: str
            Identifier of the datastruct requested.

        Returns:
        --------
        DatastructDesc
            Description of the datastruct requested.
        """
        try:
            return next(filter(lambda x: x[0] == name, self.items()))[1]
        except StopIteration:
            raise RuntimeError('Datastruct \'{}\' not in dict!'.format(name))


def serialize_obj(obj) -> str:
    """Serialize object data in base64 encoding.

    Parameters:
    -----------
    -

    Returns:
    --------
    str
        Pickled representation of the object in base64 encoding.
    """
    return encode(dumps(obj, 4), 'base64').decode()


def deserialize_obj(serialized_obj: str) -> Any:
    """Deserialize object data encoded in base64.

    Parameters:
    -----------
    serialized_obj: str
        Serialized object data.

    Returns:
    --------
    Object
        Unpickled object data.
    """
    return loads(decode(serialized_obj.encode(), "base64"))
