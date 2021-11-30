"""@package apps.offsite_codegen_db
Main script of the offsite_codegen application.
"""

from argparse import ArgumentParser, Namespace
from pathlib import Path

import offsite.config
from offsite import __version__
from offsite.codegen.codegen_util import write_codes_to_file
from offsite.codegen.generator.impl_generator_c import ImplCodeGenerator
from offsite.codegen.generator.impl_generator_cpp import ImplCodeGeneratorCPP
from offsite.codegen.generator.impl_generator_cpp_mpi import ImplCodeGeneratorCppMPI
from offsite.config import __config_ext__, init_config, Config, GeneratedCodeLanguageType, ModelToolType
from offsite.database import commit, open_db
from offsite.database.db_mapping import mapping
from offsite.descriptions.impl.impl_skeleton import ImplSkeleton
from offsite.descriptions.machine.machine import parse_machine_state
from offsite.descriptions.ode import IVP, ODEMethod
from offsite.ranking.ranking import create_rankings_ode, RankingRecord
from offsite.ranking.ranking_task import parse_ranking_tasks
from offsite.train.impl_variant import ImplVariant


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
    parser.add_argument('--db', action='store', required=True, help='Path to used database.')
    parser.add_argument('--mode', action='store', required=True, type=GeneratedCodeLanguageType,
                        help='Supported program languages: C (default), CPP, MPI')
    parser.add_argument('--tile', action='store_true', default=False,
                        help='If set to true an extra tiling version is generated for each impl variant.')
    parser.add_argument('--config', action='store', type=Path,
                        help='Tweak code generation process by passing a custom configuration ({}) file'.format(
                            __config_ext__))
    # Required options to specify which data are selected from the database.
    parser.add_argument('--ivp', action='store', required=True, type=str,
                        help='Name of the IVP for which specialized code is generated for.')
    parser.add_argument('--method', action='store', required=True, type=str,
                        help='Name of the ODE method for which specialized code is generated for.')
    parser.add_argument('--variants', action='store', type=str, default=None,
                        help='Comma-separated list of variants whose code will be generated.')
    parser.add_argument('--rank', action='store', type=str, default=None,
                        help='Create ranking of impl variants first and only generate code for all ranked '
                             'variants.')
    # Required additional options when using --rank.
    parser.add_argument('--machine', action='store', type=Path, default=None,
                        help='Path to YAML machine state description (.yaml) file. Required and used only by'
                             'argument --rank')
    parser.add_argument('--compiler', action='store', required=False, default=None,
                        help='Name of the compiler. Required and used only by argument --rank')
    parser.add_argument('-n', '--ode-size', action='store', type=int, default=False,
                        help='Create ranking for this particular ODE system size. Optional argument that would only be '
                             'required and used by argument --rank')
    return parser.parse_args()


def verify_program_args_app_codegen():
    config: Config = offsite.config.offsiteConfig
    if config.args.variants is None and config.args.rank is None:
        raise RuntimeError('Missing program argument. Either --variant or --rank are required!')
    if config.args.variants is not None and config.args.rank is not None:
        raise RuntimeError('Ambiguous variant selection argument. Must either use --variant or --rank but not both!')
    if config.args.rank is not None:
        if config.args.machine is None:
            raise RuntimeError('Missing program argument --machine required by --rank!')
        if config.args.compiler is None:
            raise RuntimeError('Missing program argument --compiler required by --rank!')


def run_code_generation():
    """Run the code generation application for use-case data from database.

    Parameters:
    -----------
    -

    Returns:
    -
    """
    config: Config = offsite.config.offsiteConfig
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
    # Retrieve impl variant IDs.
    if config.args.variants is not None:
        variant_ids = list()
        for vid in config.args.variants.split(','):
            # Parse ID ranges given by 'first:last'.
            if ':' in vid:
                vid_split = vid.split(':')
                if len(vid_split) != 2:
                    raise RuntimeError('Failed to parse impl variant ID range \'{}\'! '.format(vid) +
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
        # Parse machine state and compiler information.
        machine = parse_machine_state(config.args.machine, config.args.compiler)
        machine = machine.to_database(db_session)
        # Select ranking data from database.
        df_ranking = RankingRecord.select(db_session, machine, method.db_id, ivp.db_id, task, config.args.ode_size)
        # No fitting data found in database. Create required ranking and select again.
        if len(df_ranking) == 0:
            # Create ranking.
            create_rankings_ode(db_session, machine, [method], [ivp], [task])
            # Save new ranking in database.
            commit(db_session)
            # Select ranking data from database.
            df_ranking = RankingRecord.select(db_session, machine, method.db_id, ivp.db_id, task, config.args.ode_size)
            if len(df_ranking) == 0:
                raise RuntimeError('Unable to select the required ranking data from database.')
        # Returned ranking data entries might contain multiple impl variant IDs (as comma-separated string).
        # For further processing, we require the set of unique IDs. Thus, we need to split these strings here.
        variant_ids = set()
        for entry in df_ranking.variants_serial.unique().tolist():
            entries = (int(x) for x in entry.split(','))
            variant_ids.update(entries)
        variant_ids = sorted(variant_ids)
    else:
        assert False
    # Retrieve impl variants from database.
    skeletons = dict()
    for variant in ImplVariant.select(db_session, variant_ids):
        sid = variant.skeleton
        if sid in skeletons:
            skeletons[sid].append((variant.db_id, variant.kernels))
        else:
            skeletons[sid] = [(variant.db_id, variant.kernels)]
    # Retrieve impl skeletons from database.
    for sid in skeletons:
        variants = skeletons[sid]
        skeletons[sid] = (ImplSkeleton.select(db_session, sid), variants)
        skeleton = skeletons[sid][0]
        if skeleton.modelTool is not ModelToolType.KERNCRAFT:
            raise RuntimeError(
                'Implementation skeleton \'{}\' uses model tool \'{}\'! Application offsite_codegen supports only '
                'model tool \'{}\'.!'.format(skeleton.name, skeleton.modelTool.value, ModelToolType.KERNCRAFT.value))
    # Folders used to store generated code files.
    folder_ds = config.args.folder_ds if config.args.folder_ds is not None else config.args.folder
    folder_impl = config.args.folder_impl if config.args.folder_impl is not None else config.args.folder
    folder_ivp = config.args.folder_ivp if config.args.folder_ivp is not None else config.args.folder
    folder_method = config.args.folder_method if config.args.folder_method is not None else config.args.folder
    # Generate impl variant codes.
    for skeleton in skeletons.values():
        if config.args.mode == GeneratedCodeLanguageType.C:
            codes = ImplCodeGenerator(db_session, folder_impl, folder_ivp, folder_method).generate(
                skeleton[0], skeleton[1], ivp, method)
            write_codes_to_file(codes, suffix='')
        elif config.args.mode == GeneratedCodeLanguageType.CPP:
            codes = ImplCodeGeneratorCPP(db_session, folder_ds, folder_impl, folder_ivp, folder_method).generate(
                skeleton[0], skeleton[1], ivp, method)
            write_codes_to_file(codes, suffix='')
        elif config.args.mode == GeneratedCodeLanguageType.CPP_MPI:
            codes = ImplCodeGeneratorCppMPI(db_session, folder_ds, folder_impl, folder_ivp, folder_method).generate(
                skeleton[0], skeleton[1], ivp, method)
            write_codes_to_file(codes, suffix='')
        else:
            assert False


def run():
    """Run command line interface used for use-case 'data from database'.

    Here impl variants are created from database records.

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
    verify_program_args_app_codegen()
    # Run code generator.
    run_code_generation()
