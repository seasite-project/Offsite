"""@package apps.offsite_codegen_yaml
Main script of the offsite_codegen_from_yaml application.
"""

from argparse import ArgumentParser, Namespace
from os import remove
from pathlib import Path

import offsite.config
from offsite import __version__
from offsite.codegen.codegen_util import write_codes_to_file
from offsite.codegen.generator.impl_generator_c import ImplCodeGenerator
from offsite.codegen.generator.impl_generator_cpp import ImplCodeGeneratorCPP
from offsite.codegen.generator.impl_generator_cpp_mpi import ImplCodeGeneratorCppMPI
from offsite.config import __config_ext__, __tuning_scenario_ext__, init_config, Config, GeneratedCodeLanguageType, \
    ModelToolType
from offsite.database import commit, open_db
from offsite.database.db_mapping import mapping
from offsite.descriptions.ode import IVP, ODEMethod
from offsite.descriptions.parser_util import parse_verify_yaml_desc
from offsite.train.train_utils import deduce_available_impl_variants
from offsite.tuning_scenario import TuningScenario


def parse_program_args_app_codegen() -> Namespace:
    """Create argument parser for the offsite_codegen application.

    Parameters:
    -----------
    -

    Returns:
    --------
    ArgumentParser
        Created argument parser object.
    """
    # Create argument parser object.
    parser = ArgumentParser(description='Handle code generation requests for impl variants.')
    # Available options.
    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(__version__),
                        help='Print program version and exit.')
    parser.add_argument('--verbose', action='store_true', default=False, help='Print further information on this run.')
    parser.add_argument('--folder', action='store', default='tmp/variants', type=Path,
                        help='Path to the location of the generated code files.')
    parser.add_argument('--folder_ds', action='store', default=None, type=Path,
                        help='Path to the location of the generated datastructure code files.'
                             'Prioritized over option \'--folder\'.')
    parser.add_argument('--folder_impl', action='store', default=None, type=Path,
                        help='Path to the location of the generated impl variant code files.'
                             'Prioritized over option \'--folder\'.')
    parser.add_argument('--folder_ivp', action='store', default=None, type=Path,
                        help='Path to the location of the generated ODE problem code files.'
                             'Prioritized over option \'--folder\'.')
    parser.add_argument('--folder_method', action='store', default=None, type=Path,
                        help='Path to the location of the generated ODE method code files.'
                             'Prioritized over option \'--folder\'.')
    parser.add_argument('--db', action='store', default=None, help='Path to used database.')
    parser.add_argument('--mode', action='store', required=True, type=GeneratedCodeLanguageType,
                        help='Supported program languages: C (default), C++, MPI')
    parser.add_argument('--tile', action='store_true', default=False,
                        help='If set to true an extra tiling version is generated for each impl variant.')
    parser.add_argument('--config', action='store', type=Path,
                        help='Tweak code generation process by passing a custom configuration ({}) file'.format(
                            __config_ext__))
    # Required options to specify from which YAML descriptions impl variant code is derived.
    parser.add_argument('--scenario', action='store', required=True, type=Path,
                        help='Path to the used tuning scenario ({}) file which contains all required information on '
                             'the used machine, compiler, impl skeleton, kernel templates, ODE methods and '
                             'IVPs.'.format(__tuning_scenario_ext__))
    return parser.parse_args()


def run_code_generation():
    """Run the code generation application for use-case data from YAML.

    Parameters:
    -----------
    -

    Returns:
    -
    """
    config: Config = offsite.config.offsiteConfig
    config.pred_model_tool = ModelToolType.KERNCRAFT  # Currently only codegen for KERNCRAFT is supported
    # Connect to database.
    if config.args.db is None:
        db_session = open_db('tmp_db_codegen_yaml_config.db')
    else:
        db_session = open_db(config.args.db)
    # Read tuning scenario.
    scenario: TuningScenario = TuningScenario.from_file(config.args.scenario)
    # Parse YAML files and store in database.
    machine, skeletons, templates, methods, ivps = parse_verify_yaml_desc(db_session, scenario)
    # Select ODE method if tuning scenario contains multiple ODE methods.
    if len(methods) == 1:
        method: ODEMethod = methods[0]
    else:
        method_ids = list()
        select_str = 'Select ODE method: (specify number)\n'
        for method in methods:
            select_str += ' *  ({})  {}\n'.format(method.db_id, method.name)
            method_ids.append(method.db_id)
        # Select the used ODE method.
        method_id = -1
        while method_id not in method_ids:
            ipt = input(select_str)
            try:
                ipt_id = int(ipt)
                if ipt_id not in method_ids:
                    print('Input out of valid range!')
            except:
                print('Invalid input! Not a valid number!')
        method: ODEMethod = next(x for x in ivps if x.db_id == ipt_id)
    # Select IVP if tuning scenario contains multiple IVPs.
    if len(ivps) == 1:
        ivp: IVP = ivps[0]
    else:
        ivp_ids = list()
        select_str = 'Select IVP: (specify number)\n'
        for ivp in ivps:
            select_str += ' *  ({})  {}\n'.format(ivp.db_id, ivp.name)
            ivp_ids.append(ivp.db_id)
        # Select the used IVP.
        ipt_id = -1
        while ipt_id not in ivp_ids:
            ipt = input(select_str)
            try:
                ipt_id = int(ipt)
                if ipt_id not in ivp_ids:
                    print('Input out of valid range!')
            except:
                print('Invalid input! Not a valid number!')
        ivp: IVP = next(x for x in ivps if x.db_id == ipt_id)
    # Folders used to store generated code files.
    folder_ds = config.args.folder_ds if config.args.folder_ds is not None else config.args.folder
    folder_impl = config.args.folder_impl if config.args.folder_impl is not None else config.args.folder
    folder_ivp = config.args.folder_ivp if config.args.folder_ivp is not None else config.args.folder
    folder_method = config.args.folder_method if config.args.folder_method is not None else config.args.folder
    # Derive all available impl variants ...
    for skeleton in skeletons:
        available_variants, _ = deduce_available_impl_variants(db_session, skeleton, False)
        available_variants = [impl.to_database(db_session) for impl in available_variants]
        commit(db_session)
        # ... and generate impl variant codes.
        impl_variants = [(impl.db_id, impl.kernels) for impl in available_variants]
        if config.args.mode == GeneratedCodeLanguageType.C:
            codes = ImplCodeGenerator(db_session, folder_impl, folder_ivp, folder_method).generate(
                skeleton, impl_variants, ivp, method)
            write_codes_to_file(codes, suffix='')
        elif config.args.mode == GeneratedCodeLanguageType.CPP:
            codes = ImplCodeGeneratorCPP(db_session, folder_ds, folder_impl, folder_ivp, folder_method).generate(
                skeleton, impl_variants, ivp, method)
            write_codes_to_file(codes, suffix='')
        elif config.args.mode == GeneratedCodeLanguageType.CPP_MPI:
            codes = ImplCodeGeneratorCppMPI(db_session, folder_ds, folder_impl, folder_ivp, folder_method).generate(
                skeleton, impl_variants, ivp, method)
            write_codes_to_file(codes, suffix='')
        else:
            assert False
    if config.args.db is None:
        remove('tmp_db_codegen_yaml_config.db')


def run():
    """Run command line interface used for use-case 'data from YAML'.

    Here the data needed to create impl variants are read from YAML files. Afterwards all read YAML
    information are stored in database.

    Parameters:
    -----------
    -

    Returns:
    -
    """
    # Map database.
    mapping()
    # Create parser and parse arguments.
    args: Namespace = parse_program_args_app_codegen()
    # Create custom or default configuration.
    init_config(args)
    # Run code generator.
    run_code_generation()
