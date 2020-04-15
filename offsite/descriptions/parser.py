"""@package parser
Functions to parse YAML descriptions and create corresponding Python objects.
"""

from pathlib import Path
from typing import List

from sqlalchemy.orm import Session

import offsite.config
from offsite.config import ModelToolType, __impl_skeleton_ext__, __ivp_ext__, __kernel_template_ext__, \
    __ode_method_ext__
from offsite.descriptions.impl_skeleton import ImplSkeleton
from offsite.descriptions.ivp import IVP
from offsite.descriptions.kernel_template import KernelTemplate
from offsite.descriptions.machine import Machine
from offsite.descriptions.ode_method import ODEMethod


def parse_machine(machine_file: Path, compiler: str) -> Machine:
    """Parse machine description YAML file and return Machine object.

    Parameters:
    -----------
    machine_file : pathlib.Path
        Relative path to machine's description file.
    compiler : str
        Name of the compiler used.

    Returns:
    --------
    Machine
        Machine object created.
    """
    return Machine.from_yaml(machine_file, compiler)


def parse_kernel_templates(folder: Path, tool: ModelToolType) -> List[KernelTemplate]:
    """Parse folder for KernelTemplate objects.

    Parameters:
    -----------
    folder : Path
        Relative path to folder that contains the KernelTemplate objects.
    tool : ModelToolType
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


def parse_impl_skeletons(folder: Path, templates: List[KernelTemplate], tool: ModelToolType) -> List[ImplSkeleton]:
    """Parse folder for ImplSkeleton objects.

    Parameters:
    -----------
    folder : Path
        Relative path to folder that contains the ImplSkeleton objects.
    templates : list of KernelTemplate
        List of available KernelTemplate objects.
    tool : ModelToolType
        Limit parsing to ImplSkeleton objects that use this performance modelling tool.

    Returns:
    list of ImplSkeleton
        ImplSkeleton objects to be used.
    """
    # Filter implementation skeletons.
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


def parse_method(method_file: Path) -> ODEMethod:
    """Parse ODE method description YAML file and return ODEMethod object.

    Parameters:
    -----------
    method_file : Path
        Relative path to ODE method's description file.

    Returns:
    --------
    ODEMethod
        ODEMethod object used.
    """
    return ODEMethod.from_yaml(method_file)


def parse_methods(path: Path) -> List[ODEMethod]:
    """Parse folder or file path for ODEMethod objects.

    Parameters:
    -----------
    path : Path
        Relative path to folder that contains the ODEMethod objects or a particular ODE method object.

    Returns:
    --------
    list of ODEMethod
        ODEMethod objects to be used.
    """
    if path.is_file():
        if not path.suffix == __ode_method_ext__:
            raise RuntimeError('Invalid extension for ODE method: {}!'.format(path.suffix))
        return [parse_method(path)]
    return [parse_method(mf) for mf in (f for f in path.iterdir() if f.suffix == __ode_method_ext__)]


def parse_ivp(ivp_file: Path) -> IVP:
    """Parse IVP description YAML file and return IVP object.

    Parameters:
    -----------
    ivp_file : Path
        Relative path to IVP's description file.

    Returns:
    --------
    IVP
        IVP object used.
    """
    return IVP.from_yaml(ivp_file)


def parse_ivps(path: Path) -> List[IVP]:
    """Parse folder or file path for IVP objects.

    path : str
        Relative path to folder that contains the IVP objects or a particular IVP object.

    Returns:
    --------
    list of IVP
        IVP objects created.
    """
    if path.is_file():
        if not path.suffix == __ivp_ext__:
            raise RuntimeError('Invalid extension for IVP: {}!'.format(path.suffix))
        return [parse_ivp(path)]
    return [parse_ivp(ivpf) for ivpf in (f for f in path.iterdir() if f.suffix == __ivp_ext__)]


def parse_verify_yaml_desc(db_session: Session):
    """Parse and verify the YAML description files given by the user arguments.

    Parses and verifies the YAML descriptions given by the program arguments. For each description file corresponding
    Python class objects are created (Machine, ImplSkeleton, KernelTemplate, ...) and database records are created.

    Parameters:
    -----------
    db_session : sqlalchemy.orm.session.Session
        Used database session.

    Returns:
    --------
    Machine
        Description object of the machine used.
    List of ImplSkeleton
        List of description objects of the implementation skeletons used.
    List of KernelTemplate
        List of description objects of the kernel templates used.
    ODEMethod
        List of description objects of the ODE methods used.
    List of IVP
        List of description objects of the IVPs used.
    """
    args = offsite.config.offsiteConfig.args
    # Parse the machine description.
    machine = parse_machine(args.machine, args.compiler)
    machine = machine.to_database(db_session)
    # Parse the kernel template descriptions.
    templates = parse_kernel_templates(args.kernel, args.tool)
    templates = [template.to_database(db_session) for template in templates]
    # Parse the implementation skeleton descriptions.
    skeletons = parse_impl_skeletons(args.impl, templates, args.tool)
    skeletons = [skeleton.to_database(db_session) for skeleton in skeletons]
    # Parse the ODE method descriptions.
    methods = parse_methods(args.method)
    methods = [method.to_database(db_session) for method in methods]
    # Parse the IVP descriptions.
    ivps = parse_ivps(args.ivp)
    ivps = [ivp.to_database(db_session) for ivp in ivps]
    # Return parsed objects.
    return machine, skeletons, templates, methods, ivps


def print_yaml_desc(machine: Machine, skeletons: List[ImplSkeleton], templates: List[ImplSkeleton],
                    methods: List[ODEMethod], ivps: List[IVP]):
    """Print information on the parsed YAML descriptions.

    Parameters:
    -----------
    Machine
        Description object of the machine used.
    List of ImplSkeleton
        List of description objects of the implementation skeletons used.
    List of KernelTemplate
        List of description objects of the kernel templates used.
    ODEMethod
        List of description objects of the ODE methods used.
    List of IVP
        List of description objects of the IVPs used.

    Returns:
    --------
    -
    """
    print('')
    print('Machine: {}'.format(machine.path))
    print('--------')
    print('* {}'.format(machine.modelName))
    print('* {}'.format(machine.modelType))
    print('')
    print('Compiler:')
    print('---------')
    print('* {} {}'.format(machine.compiler.name, machine.compiler.version))
    print('')
    print('Implementation skeletons:')
    print('-------------------------')
    for skeleton in skeletons:
        print('* {} (Kernels: {})'.format(skeleton.name, ', '.join(skeleton.kernels)))
    print('')
    print('Kernel templates:')
    print('-----------------')
    for temp in templates:
        print('* {}'.format(temp.name))
    if methods:
        print('')
        print('ODE methods:')
        print('------------')
        for method in methods:
            print('* {} (stages={}, corrector steps={})'.format(method.name, method.stages, method.correctorSteps))
    if ivps:
        print('')
        print('IVPs:')
        print('-----')
        for ivp in ivps:
            print('* {} (grid dim = {})'.format(ivp.name, ivp.gridSize))
    print('')
