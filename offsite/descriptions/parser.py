"""@package parser
Functions to parse YAML descriptions and create corresponding Python objects.
"""

from pathlib import Path
from typing import List, Optional

from sqlalchemy.orm import Session

import offsite.config
from offsite.config import ModelToolType, RankingCutoffType, __impl_skeleton_ext__, __ivp_ext__, \
    __kernel_template_ext__, __ode_method_ext__
from offsite.descriptions.impl_skeleton import ImplSkeleton
from offsite.descriptions.ivp import IVP
from offsite.descriptions.kernel_template import KernelTemplate
from offsite.descriptions.machine import Machine
from offsite.descriptions.ode_method import ODEMethod
from offsite.descriptions.ranking_task import RankTask, RankDeviationTask, RankOrderTask


def parse_machine(machine_file: Path, compiler: str) -> Machine:
    """Parse machine description YAML file and return Machine object.

    Parameters:
    -----------
    machine_file: pathlib.Path
        Relative path to machine's description file.
    compiler: str
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


def parse_method(path: Path) -> ODEMethod:
    """Parse ODE method description YAML file and return ODEMethod object.

    Parameters:
    -----------
    path: Path
        Relative path to file that contains the ODEMethod objects.

    Returns:
    --------
    ODEMethod
        ODEMethod object created.
    """
    if not path.suffix == __ode_method_ext__:
        raise RuntimeError('Invalid extension for ODE method: {}!'.format(path.suffix))
    return ODEMethod.from_yaml(path)


def parse_methods(path: Path) -> List[ODEMethod]:
    """Parse folder or file path for ODEMethod objects.

    Parameters:
    -----------
    path: Path
        Relative path to folder that contains the ODEMethod objects or a particular ODE method object.

    Returns:
    --------
    list of ODEMethod
        ODEMethod objects to be used.
    """
    if path.is_file():
        return [parse_method(path)]
    return [ODEMethod.from_yaml(mf) for mf in (f for f in path.iterdir() if f.suffix == __ode_method_ext__)]


def parse_ivp(path: Path, tool: ModelToolType) -> IVP:
    """Parse IVP description YAML file and return IVP object.

    Parameters:
    -----------
    path: str
        Relative path to file that contains the IVP object.
    tool: ModelToolType
        Limit parsing to IVP objects that use this performance modelling tool.

    Returns:
    --------
    IVP
        IVP object created.
    """
    if not path.suffix == __ivp_ext__:
        raise RuntimeError('Invalid extension for IVP: {}!'.format(path.suffix))
    ivp = IVP.from_yaml(path)
    if tool and ivp.modelTool is not tool:
        raise RuntimeError('Invalid model tool for IVP: {}!'.format(path.stem))
    return ivp


def parse_ivps(path: Path, tool: Optional[ModelToolType]) -> List[IVP]:
    """Parse folder or file path for IVP objects.

    path: str
        Relative path to folder that contains the IVP objects or a particular IVP object.
    tool: ModelToolType
        Limit parsing to IVP objects that use this performance modelling tool.

    Returns:
    --------
    list of IVP
        IVP objects created.
    """
    # Single IVP.
    if path.is_file():
        return [parse_ivp(path, tool)]
    # Filter IVP folder.
    if tool:
        ivps = list()
        for f in (f for f in path.iterdir() if f.suffix == __ivp_ext__):
            try:
                ivp = IVP.from_yaml(f)
                if ivp.modelTool == tool:
                    ivps.append(ivp)
            except RuntimeError:
                pass
        return ivps
    # Parse IVP folder.
    return [IVP.from_yaml(ivpf) for ivpf in (f for f in path.iterdir() if f.suffix == __ivp_ext__)]


def parse_verify_yaml_desc(db_session: Session):
    """Parse and verify the YAML description files given by the user arguments.

    Parses and verifies the YAML descriptions given by the program arguments. For each description file corresponding
    Python class objects are created (Machine, ImplSkeleton, KernelTemplate, ...) and database records are created.

    Parameters:
    -----------
    db_session: sqlalchemy.orm.session.Session
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
    if not templates:
        raise RuntimeError('No valid kernel templates found: \'{}\''.format(args.kernel))
    # Parse the implementation skeleton descriptions.
    skeletons = parse_impl_skeletons(args.impl, templates, args.tool)
    skeletons = [skeleton.to_database(db_session) for skeleton in skeletons]
    if not skeletons:
        raise RuntimeError('No valid implementation skeletons found: \'{}\''.format(args.impl))
    # Parse the ODE method descriptions.
    methods = parse_methods(args.method)
    methods = [method.to_database(db_session) for method in methods]
    if not methods:
        raise RuntimeError('No valid ODE methods found: \'{}\''.format(args.method))
    # Parse the IVP descriptions.
    ivps = parse_ivps(args.ivp, args.tool)
    ivps = [ivp.to_database(db_session) for ivp in ivps]
    if not ivps:
        raise RuntimeError('No valid IVPs found: \'{}\''.format(args.ivp))
    # Return parsed objects.
    return machine, skeletons, templates, methods, ivps


def print_yaml_desc(machine: Machine, skeletons: Optional[List[ImplSkeleton]], templates: Optional[List[ImplSkeleton]],
                    methods: Optional[List[ODEMethod]], ivps: Optional[List[IVP]]):
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
    if skeletons is not None:
        print('')
        print('Implementation skeletons:')
        print('-------------------------')
        for skeleton in skeletons:
            print('* {} (Kernels: {})'.format(skeleton.name, ', '.join(skeleton.kernels)))
    if templates is not None:
        print('')
        print('Kernel templates:')
        print('-----------------')
        for temp in templates:
            print('* {}'.format(temp.name))
    if methods is not None:
        print('')
        print('ODE methods:')
        print('------------')
        for method in methods:
            print('* {} (stages={}, corrector steps={})'.format(method.name, method.stages, method.correctorSteps))
    if ivps is not None:
        print('')
        print('IVPs:')
        print('-----')
        for ivp in ivps:
            print('* {} (grid dim = {})'.format(ivp.name, ivp.gridSize))
    print('')


def parse_ranking_tasks(rank_task_str: str) -> List[RankTask]:
    """Parse the ranking task string provided and create corresponding RankTask objects.

    Parameters:
    -----------
    tasks_str: str
        Used ranking task string.

    Returns:
    --------
    list of RankTask
        Created ranking tasks.
    """
    # Parse single ranking task and add to ranking task list.
    tasks: List[RankTask] = list()
    for task in rank_task_str.split(','):
        task_args: List[str] = task.split(':')
        if len(task_args) != 3:
            print('Skipping invalid ranking task: \'{}\'!'.format(task))
            continue
        # Check cutoff criteria.
        if task_args[1] not in ['t', 'p']:
            print(
                'Skipping invalid ranking task: \'{}\'! Unknown cutoff criteria \'{}\'!'.format(task, task_args[1]))
            continue
        cutoff_criteria: RankingCutoffType = RankingCutoffType.TOP if task_args[1] == 't' else RankingCutoffType.PERCENT
        # Check cutoff value.
        cutoff_val: float = 0.0
        try:
            cutoff_val = float(task_args[2])
        except ValueError:
            print('Skipping invalid ranking task: \'{}\'! Invalid cutoff value \'{}\'!'.format(task, cutoff_val))
            continue
        # Create ranking task and append to list.
        rank_criteria: str = task_args[0]
        if rank_criteria not in ['o', 'd']:
            print('Skipping invalid ranking task: \'{}\'! Unknown ranking criteria \'{}\'!'.format(task, rank_criteria))
            continue
        tasks.append(
            RankOrderTask(cutoff_criteria, cutoff_val) if rank_criteria == 'o' else RankDeviationTask(cutoff_criteria,
                                                                                                      cutoff_val))
    return tasks
