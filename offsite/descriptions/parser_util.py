"""@package descriptions.parser_util
Functions to parse YAML description as well as different utility functions and classes used during YAML parsing.

@author: Johannes Seiferth
"""

from typing import Dict, List, Optional

from sqlalchemy.orm import Session

import offsite.config
from offsite.config import __impl_skeleton_ext__, __ivp_ext__, __kernel_template_ext__, __ode_method_ext__, Config
from offsite.descriptions.impl.impl_skeleton import ImplSkeleton, parse_impl_skeletons
from offsite.descriptions.impl.kernel_template import KernelTemplate, parse_kernel_templates
from offsite.descriptions.machine import MachineState, parse_machine_state
from offsite.descriptions.ode import IVP, ODEMethod, parse_ivp, parse_ivps, parse_methods
from offsite.solver import SolverSpecificTableType
from offsite.tuning_scenario import TuningScenario


def parse_verify_yaml_desc(db_session: Session, scenario: TuningScenario):
    """Parse and verify the YAML description files given by the user arguments.

    Parses and verifies the YAML descriptions given by the tuning scenario passed. For each description file
    corresponding Python class objects are created (MachineState, ImplSkeleton, KernelTemplate, ...) and database records
    are created.

    Parameters:
    -----------
    db_session: sqlalchemy.orm.session.Session
        Used database session.
    scenario: TuningScenario
        Used tuning scenario.

    Returns:
    --------
    MachineState
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
    config: Config = offsite.config.offsiteConfig
    # Parse the machine description.
    machine = parse_machine_state(scenario.machine, scenario.compiler)
    machine = machine.to_database(db_session)
    # Parse the kernel template descriptions.
    templates_dict: Dict[str, KernelTemplate] = dict()
    if scenario.template_path is not None:
        for p in (p for p in scenario.template_path if p.suffix == __kernel_template_ext__):
            t = KernelTemplate.from_yaml(p)
            templates_dict[t.name] = t
    if scenario.template_dir is not None:
        for d in scenario.template_dir:
            for t in parse_kernel_templates(d, config.pred_model_tool):
                templates_dict[t.name] = t
    templates: List[KernelTemplate] = [template.to_database(db_session) for template in templates_dict.values()]
    if not templates:
        raise RuntimeError('No valid kernel templates found: \'{}\''.format(config.args.kernel))
    # Parse the implementation skeleton descriptions if required by solver.
    skeletons_dict: Dict[str, ImplSkeleton] = dict()
    if scenario.skeleton_path is not None:
        for p in (p for p in scenario.skeleton_path if p.suffix == __impl_skeleton_ext__):
            s = ImplSkeleton.from_yaml(p, templates)
            skeletons_dict[s.name] = s
    if scenario.skeleton_dir is not None:
        for d in scenario.skeleton_dir:
            for s in parse_impl_skeletons(d, templates, config.pred_model_tool):
                skeletons_dict[s.name] = s
    skeletons: List[ImplSkeleton] = [skeleton.to_database(db_session) for skeleton in skeletons_dict.values()]
    if not skeletons:
        raise RuntimeError('No valid implementation skeletons found.')
    # Filter out all kernel templates which are not actually used by any of the implementation skeletons given.
    # We need to check here if the arg exists since this it is not available in all apps.
    if 'include_all_kernels' in vars(config.args) and not config.args.include_all_kernels:
        tmp_list: List[KernelTemplate] = list()
        for skeleton in skeletons:
            for kernel in skeleton.kernels.keys():
                for template in templates:
                    if kernel == template.name:
                        tmp_list.append(template)
                        templates.remove(template)
        templates = tmp_list
    # Parse the ODE method descriptions if required by solver.
    methods: List[ODEMethod] = list()
    if SolverSpecificTableType.ODE_METHOD in scenario.solver.specific_tables:
        methods_dict: Dict[str, ODEMethod] = dict()
        if scenario.method_path is not None:
            for p in (p for p in scenario.method_path if p.suffix == __ode_method_ext__):
                m = ODEMethod.from_yaml(p)
                methods_dict[m.name] = m
        if scenario.method_dir is not None:
            for d in scenario.method_dir:
                for m in parse_methods(d):
                    methods_dict[m.name] = m
        methods = [method.to_database(db_session) for method in methods_dict.values()]
        if not methods:
            raise RuntimeError('No valid ODE methods found.')
    # Parse the IVP descriptions if required by solver.
    ivps: List[IVP] = list()
    if SolverSpecificTableType.IVP in scenario.solver.specific_tables:
        ivps_dict: Dict[str, IVP] = dict()
        if scenario.ivp_path is not None:
            for p in (p for p in scenario.ivp_path if p.suffix == __ivp_ext__):
                i = parse_ivp(p, config.pred_model_tool)
                ivps_dict[i.name] = i
        if scenario.ivp_dir is not None:
            for d in scenario.ivp_dir:
                for i in parse_ivps(d, config.pred_model_tool):
                    if i.name not in methods:
                        ivps_dict[i.name] = i
        ivps = [ivp.to_database(db_session) for ivp in ivps_dict.values()]
        if not ivps:
            raise RuntimeError('No valid IVPs found.')
    # Return parsed objects.
    return machine, skeletons, templates, methods, ivps


def print_yaml_desc(
        machine: MachineState, skeletons: Optional[List[ImplSkeleton]], templates: Optional[List[ImplSkeleton]],
        methods: Optional[List[ODEMethod]], ivps: Optional[List[IVP]]):
    """Print information on the parsed YAML descriptions.

    Parameters:
    -----------
    MachineState
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
