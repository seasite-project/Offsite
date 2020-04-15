"""@package offsite_codegen
Main script of the offsite_codegen application.
"""

from argparse import ArgumentParser
from pathlib import Path

from offsite import __version__
from offsite.codegen.impl_generator import ImplCodeGenerator
from offsite.config import Config, __config_ext__, ModelToolType
from offsite.db.db import open_db
from offsite.db.db_mapping import mapping
from offsite.descriptions.impl_skeleton import ImplSkeleton, ImplVariant
from offsite.descriptions.ivp import IVP
from offsite.descriptions.ode_method import ODEMethod


def parse_program_args_app_codegen() -> 'argparse.Namespace':
    """Parse the available program arguments of the offsite_code_server application.

    Parameters:
    -----------
    -

    Returns:
    --------
    argparse.Namespace
        Parsed program arguments.
    """
    # Create argument parser object.
    parser = ArgumentParser(
        description='Start code server to handle code generation requests for implementation variants.')
    # Available options.
    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(__version__),
                        help='Print program version and exit.')
    parser.add_argument('--db', action='store', default='tune.db', help='Path to database. Default: tune.db.')
    parser.add_argument('--folder', action='store', default='tmp/variants', type=Path,
                        help='Path to the location of the generated code files.')
    parser.add_argument('--ivp', action='store', required=True, type=str,
                        help='Name of the IVP for which specialized code is generated for.')
    parser.add_argument('--method', action='store', required=True, type=str,
                        help='Name of the ODE method for which specialized code is generated for.')
    parser.add_argument('--variants', action='store', required=True, type=str,
                        help='Path to folder containing YAML kernel template description ({}) files.')
    parser.add_argument('--config', action='store', type=Path,
                        help='Tweak code generation process by passing a custom configuration ({}) file'.format(
                            __config_ext__))
    # Parse program arguments.
    return parser.parse_args()


def run_code_generation(config: 'Config'):
    """Run the code generation application.

    Parameters:
    -----------
    config: Config
        Program configuration.

    Returns:
    -
    """
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
    # Retrieve implementation variants from database.
    skeletons = dict()
    # Retrieve implementation variant IDs.
    variant_ids = config.args.variants.split(',')
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
                'Implementation skeleton \'{}\' uses model tool \'{}\'! Application offsite_codegen supports only model tool \'{}\'.!'.format(
                skeleton.name, skeleton.modelTool.value, ModelToolType.KERNCRAFT.value))
    # Generate implementation variant codes.
    impl_generator = ImplCodeGenerator(db_session=db_session, folder=config.args.folder, config=config)
    for skeleton in skeletons.values():
        generated_codes = impl_generator.generate(skeleton[0], skeleton[1], ivp, method)
        impl_generator.write_to_file(generated_codes)


def run():
    """Run command line interface.

    Parameters:
    -----------
    -

    Returns:
    -
    """
    # Map database.
    mapping()
    # Create parser and parse arguments.
    args = parse_program_args_app_codegen()
    # Parse custom configuration file if present else use default configuration.
    config = Config.from_file(args) if args.config else Config(args)
    # Run code generator.
    run_code_generation(config)
