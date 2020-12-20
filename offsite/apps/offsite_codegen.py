"""@package offsite_codegen
Main script of the offsite_codegen application.
"""

from argparse import ArgumentParser, Namespace
from enum import Enum
from pathlib import Path

import offsite.config
from offsite import __version__
from offsite.codegen.impl_generator import ImplCodeGenerator
from offsite.codegen.impl_generator_cpp import ImplCodeGeneratorCpp
from offsite.config import __config_ext__, ModelToolType, init_config, __impl_skeleton_ext__, __ivp_ext__, \
    __kernel_template_ext__, __ode_method_ext__
from offsite.db.db import open_db, commit
from offsite.db.db_mapping import mapping
from offsite.descriptions.impl_skeleton import ImplSkeleton
from offsite.descriptions.impl_variant import ImplVariant
from offsite.descriptions.ivp import IVP
from offsite.descriptions.ode_method import ODEMethod
from offsite.descriptions.parser import parse_impl_skeletons, parse_ivp, parse_kernel_templates, parse_method, \
    parse_machine, parse_ranking_tasks
from offsite.evaluation.ranking import create_rankings
from offsite.evaluation.records import RankingRecord
from offsite.train.train_impl import deduce_available_impl_variants


class LanguageType(Enum):
    """Defines what type of code is generated.

    - C
        C style code.
    - CPP
        C++ style code.

    """
    C = 'C'
    CPP = 'CPP'


def create_args_parser_app_codegen() -> ArgumentParser:
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
    parser = ArgumentParser(description='Handle code generation requests for implementation variants.')
    # Available options.
    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(__version__),
                        help='Print program version and exit.')
    parser.add_argument('--verbose', action='store_true', default=False, help='Print further information on this run.')
    parser.add_argument('--folder', action='store', default='tmp/variants', type=Path,
                        help='Path to the location of the generated code files.')
    parser.add_argument('--folder_impl', action='store', default=None, type=Path,
                        help='Path to the location of the generated implementation variant code files.'
                             'Prioritized over option \'--folder\'.')
    parser.add_argument('--folder_ivp', action='store', default=None, type=Path,
                        help='Path to the location of the generated ODE problem code files.'
                             'Prioritized over option \'--folder\'.')
    parser.add_argument('--folder_method', action='store', default=None, type=Path,
                        help='Path to the location of the generated ODE method code files.'
                             'Prioritized over option \'--folder\'.')
    parser.add_argument('--db', action='store', default='tune.db', help='Path to database. Default: tune.db.')
    parser.add_argument('--mode', action='store', required=True, type=LanguageType,
                        help='Supported program languages: C (default), C++')
    parser.add_argument('--tile', action='store_true', default=False,
                        help='If set to true an extra tiling version is generated for each implementation variant.')
    parser.add_argument('--config', action='store', type=Path,
                        help='Tweak code generation process by passing a custom configuration ({}) file'.format(
                            __config_ext__))
    return parser


def parse_program_args_app_codegen_db() -> Namespace:
    """Parse the available program arguments of the offsite_codegen application for use-case data from database.

    Parameters:
    -----------
    -

    Returns:
    --------
    argparse.Namespace
        Parsed program arguments.
    """
    # Create argument parser object.
    parser = create_args_parser_app_codegen()
    # Required options to specify which data are selected from the database.
    parser.add_argument('--ivp', action='store', required=True, type=str,
                        help='Name of the IVP for which specialized code is generated for.')
    parser.add_argument('--method', action='store', required=True, type=str,
                        help='Name of the ODE method for which specialized code is generated for.')
    parser.add_argument('--variants', action='store', type=str, default=None,
                        help='Comma-separated list of variants whose code will be generated.')
    parser.add_argument('--rank', action='store', type=str, default=None,
                        help='TODO')
    #
    parser.add_argument('--machine', action='store', type=Path, default=None,
                        help='Path to YAML machine description (.yaml) file. Required and used only by argument --rank')
    parser.add_argument('--compiler', action='store', required=False, default=None,
                        help='Name of the compiler. Required and used only by argument --rank')
    parser.add_argument('-n', '--ode-size', action='store', type=int, default=False,
                        help='Tune for a fixed ODE system size only. If argument is not passed the working set model '
                             'is applied to determine the set of tested ODE system sizes. Optional argument that would '
                             'only be required and used by argument --rank')
    # Parse program arguments.
    return parser.parse_args()


def verify_program_args_app_codegen_db():
    config = offsite.config.offsiteConfig
    if config.args.variants is None and config.args.rank is None:
        raise RuntimeError('Missing program argument. Either --variant or --rank are required!')
    if config.args.variants is not None and config.args.rank is not None:
        raise RuntimeError('Ambiguous variant selection argument. Must either use --variant or --rank but not both!')
    if config.args.rank is not None:
        if config.args.machine is None:
            raise RuntimeError('Missing program argument --machine required by --rank!')
        if config.args.compiler is None:
            raise RuntimeError('Missing program argument --compiler required by --rank!')


def parse_program_args_app_codegen_yaml() -> Namespace:
    """Parse the available program arguments of the offsite_codegen application for use-case data from YAML.

    Parameters:
    -----------
    -

    Returns:
    --------
    argparse.Namespace
        Parsed program arguments.
    """
    # Create argument parser object.
    parser = create_args_parser_app_codegen()
    # Required options to specify from which YAML descriptions implementation variant code is derived.
    parser.add_argument('--kernel', action='store', required=True, type=Path,
                        help='Path to folder containing YAML kernel template description ({}) files.'.format(
                            __kernel_template_ext__))
    parser.add_argument('--impl', action='store', required=True, type=Path,
                        help='Path to folder containing YAML implementation skeleton description ({}) files.'.format(
                            __impl_skeleton_ext__))
    parser.add_argument('--method', action='store', required=True, type=Path,
                        help='Path to a ODE method description ({}) file'.format(__ode_method_ext__))
    parser.add_argument('--ivp', action='store', required=True, type=Path,
                        help='Path to a IVP description ({}) file.'.format(__ivp_ext__))
    # Parse program arguments.
    return parser.parse_args()


def run_code_generation_db():
    """Run the code generation application for use-case data from database.

    Parameters:
    -----------
    -

    Returns:
    -
    """
    config = offsite.config.offsiteConfig
    # Connect to database.
    db_session = open_db(config.args.db)
    # Retrieve IVP from database.
    ivp = IVP.from_database(db_session, config.args.ivp)
    if ivp.modelTool is not ModelToolType.KERNCRAFT:
        raise RuntimeError(
            'IVP \'{}\' uses model tool \'{}\'! Application offsite_codegen supports only model tool \'{}\'.!'.format(
                ivp.name, ivp.modelTool.value, ModelToolType.KERNCRAFT.value))
    # Retrieve ODE method from database.
    method = ODEMethod.from_database(db_session, config.args.method)
    # Retrieve implementation variant IDs.
    if config.args.variants is not None:
        variant_ids = list()
        for vid in config.args.variants.split(','):
            # Parse ID ranges given by 'first:last'.
            if ':' in vid:
                vid_split = vid.split(':')
                if len(vid_split) != 2:
                    raise RuntimeError('Failed to parse implementation variant ID range \'{}\'! '.format(vid) +
                                       'Supported syntax \'[first]:[last]\'')
                first = int(vid_split[0])
                last = int(vid_split[1])
                variant_ids.extend([x for x in range(first, last + 1)])
            else:
                variant_ids.append(int(vid))
        # Remove duplicate IDs.
        variant_ids = list(set(variant_ids))
        # Sort variant IDS.
        variant_ids = sorted(variant_ids)
    elif config.args.rank is not None:
        task = parse_ranking_tasks(config.args.rank)
        print(task)
        if len(task) == 0:
            raise RuntimeError('Codegen failed: no valid rank task specified in \'{}\'!'.format(config.args.rank))
        if len(task) != 1:
            raise RuntimeError('Codegen failed: using multiple rank tasks in one run is not supported!')
        task = task[0]
        # Parse machine and compiler information.
        machine = parse_machine(config.args.machine, config.args.compiler)
        machine = machine.to_database(db_session)
        # Select ranking data from database.
        df_ranking = RankingRecord.select(db_session, machine, method.db_id, ivp.db_id, task, config.args.ode_size)
        # No fitting data found in database. Create required ranking and select again.
        if len(df_ranking) == 0:
            # Create ranking.
            create_rankings(db_session, machine, [method], [ivp], [task])
            # Save new ranking in database.
            commit(db_session)
            # Select ranking data from database.
            df_ranking = RankingRecord.select(db_session, machine, method.db_id, ivp.db_id, task, config.args.ode_size)
            if len(df_ranking) == 0:
                raise RuntimeError('Unable to select the required ranking data from database.')
        # Returned ranking data entries might contain multiple implementation variant IDs (as comma-separated string).
        # For further processing, we require the set of unique IDs. Thus, we need to split these strings here.
        variant_ids = set()
        for entry in df_ranking.variants_serial.unique().tolist():
            entries = (int(x) for x in entry.split(','))
            variant_ids.update(entries)
        variant_ids = sorted(variant_ids)
    else:
        assert False
    # Retrieve implementation variants from database.
    skeletons = dict()
    for variant in ImplVariant.select(db_session, variant_ids):
        sid = variant.skeleton
        if sid in skeletons:
            skeletons[sid].append((variant.db_id, variant.kernels))
        else:
            skeletons[sid] = [(variant.db_id, variant.kernels)]
    # Retrieve implementation skeletons from database.
    for sid in skeletons:
        variants = skeletons[sid]
        skeletons[sid] = (ImplSkeleton.select(db_session, sid), variants)
        skeleton = skeletons[sid][0]
        if skeleton.modelTool is not ModelToolType.KERNCRAFT:
            raise RuntimeError(
                'Implementation skeleton \'{}\' uses model tool \'{}\'! Application offsite_codegen supports only ' \
                'model tool \'{}\'.!'.format(skeleton.name, skeleton.modelTool.value, ModelToolType.KERNCRAFT.value))
    # Folders used to store generated code files.
    folder_impl = config.args.folder_impl if config.args.folder_impl is not None else config.args.folder
    folder_ivp = config.args.folder_ivp if config.args.folder_ivp is not None else config.args.folder
    folder_method = config.args.folder_method if config.args.folder_method is not None else config.args.folder
    # Generate implementation variant codes.
    for skeleton in skeletons.values():
        if config.args.mode == LanguageType.C:
            codes = ImplCodeGenerator(db_session, folder_impl, folder_ivp, folder_method).generate(
                skeleton[0], skeleton[1], ivp, method)
            ImplCodeGenerator.write_codes_to_file(codes)
        elif config.args.mode == LanguageType.CPP:
            codes = ImplCodeGeneratorCpp(db_session, folder_impl, folder_ivp, folder_method).generate(
                skeleton[0], skeleton[1], ivp, method)
            ImplCodeGeneratorCpp.write_codes_to_file(codes)
        else:
            assert False


def run_code_generation_yaml():
    """Run the code generation application for use-case data from YAML.

    Parameters:
    -----------
    -

    Returns:
    -
    """
    config = offsite.config.offsiteConfig
    args = config.args
    tool = ModelToolType.KERNCRAFT  # currently only codegen for KERNCRAFT is supported
    # Connect to database.
    db_session = open_db(config.args.db)
    # Parse YAML files and store in database.
    # .. parse the kernel template descriptions.
    templates = parse_kernel_templates(args.kernel, tool)
    templates = [template.to_database(db_session) for template in templates]
    # .. parse the implementation skeleton descriptions.
    skeletons = parse_impl_skeletons(args.impl, templates, tool)
    skeletons = [skeleton.to_database(db_session) for skeleton in skeletons]
    # .. parse the ODE method descriptions.
    method = parse_method(args.method)
    method = method.to_database(db_session)
    # .. parse the IVP descriptions.
    ivp = parse_ivp(args.ivp, tool)
    ivp = ivp.to_database(db_session)
    # Folders used to store generated code files.
    folder_impl = config.args.folder_impl if config.args.folder_impl is not None else config.args.folder
    folder_ivp = config.args.folder_ivp if config.args.folder_ivp is not None else config.args.folder
    folder_method = config.args.folder_method if config.args.folder_method is not None else config.args.folder
    # Derive all available implementation variants ...
    for skeleton in skeletons:
        available_variants, _ = deduce_available_impl_variants(db_session, skeleton)
        available_variants = [impl.to_database(db_session) for impl in available_variants]
        commit(db_session)
        # ... and generate implementation variant codes.
        impl_variants = [(impl.db_id, impl.kernels) for impl in available_variants]
        if config.args.mode == LanguageType.C:
            codes = ImplCodeGenerator(db_session, folder_impl, folder_ivp, folder_method).generate(
                skeleton, impl_variants, ivp, method)
            ImplCodeGenerator.write_codes_to_file(codes)
        elif config.args.mode == LanguageType.CPP:
            codes = ImplCodeGeneratorCpp(db_session, folder_impl, folder_ivp, folder_method).generate(
                skeleton, impl_variants, ivp, method)
            ImplCodeGeneratorCpp.write_codes_to_file(codes)
        else:
            assert False


def run_db():
    """Run command line interface used for use-case data from database.

    Here implementation variants are created from database records.

    Parameters:
    -----------
    -

    Returns:
    -
    """
    # Map database.
    mapping()
    # Create parser and parse arguments.
    args: Namespace = parse_program_args_app_codegen_db()
    # Create custom or default configuration.
    init_config(args)
    verify_program_args_app_codegen_db()
    # Run code generator.
    run_code_generation_db()


def run_yaml():
    """Run command line interface used for use-case data from YAML.

    Here the data needed to create implementation variants are read from YAML files. Afterwards all read YAML
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
    args: Namespace = parse_program_args_app_codegen_yaml()
    # Create custom or default configuration.
    init_config(args)
    # Run code generator.
    run_code_generation_yaml()
