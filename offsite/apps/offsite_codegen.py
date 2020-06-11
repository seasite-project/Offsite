"""@package offsite_codegen
Main script of the offsite_codegen application.
"""

from argparse import ArgumentParser
from enum import Enum
from pathlib import Path

import offsite.config
from offsite import __version__
from offsite.codegen.impl_generator import ImplCodeGenerator
from offsite.config import __config_ext__, ModelToolType, init_config, __impl_skeleton_ext__, __ivp_ext__, \
    __kernel_template_ext__, __ode_method_ext__
from offsite.db.db import open_db, commit
from offsite.db.db_mapping import mapping
from offsite.descriptions.impl_skeleton import ImplSkeleton, ImplVariant
from offsite.descriptions.ivp import IVP
from offsite.descriptions.ode_method import ODEMethod
from offsite.descriptions.parser import parse_impl_skeletons, parse_ivp, parse_kernel_templates, parse_method
from offsite.train.train_impl import deduce_available_impl_variants


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
    parser.add_argument('--folder', action='store', default='tmp/variants', type=Path,
                        help='Path to the location of the generated code files.')
    parser.add_argument('--db', action='store', default='tune.db', help='Path to database. Default: tune.db.')
    parser.add_argument('--config', action='store', type=Path,
                        help='Tweak code generation process by passing a custom configuration ({}) file'.format(
                            __config_ext__))
    return parser


def parse_program_args_app_codegen_db() -> 'argparse.Namespace':
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
    parser.add_argument('--variants', action='store', required=True, type=str,
                        help='Comma-seperated list of variants whose code will be generated.')
    # Parse program arguments.
    return parser.parse_args()


def parse_program_args_app_codegen_yaml() -> 'argparse.Namespace':
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
                'Implementation skeleton \'{}\' uses model tool \'{}\'! Application offsite_codegen supports only ' \
                'model tool \'{}\'.!'.format(skeleton.name, skeleton.modelTool.value, ModelToolType.KERNCRAFT.value))
    # Generate implementation variant codes.
    for skeleton in skeletons.values():
        codes = ImplCodeGenerator(db_session, config.args.folder).generate(skeleton[0], skeleton[1], ivp, method)
        ImplCodeGenerator.write_codes_to_file(codes)


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
    # Derive all available implementation variants ...
    for skeleton in skeletons:
        available_variants, _ = deduce_available_impl_variants(db_session, skeleton)
        available_variants = [impl.to_database(db_session) for impl in available_variants]
        commit(db_session)
        # ... and generate implementation variant codes.
        impl_variants = [(impl.db_id, impl.kernels) for impl in available_variants]
        codes = ImplCodeGenerator(db_session, config.args.folder).generate(skeleton, impl_variants, ivp, method)
        ImplCodeGenerator.write_codes_to_file(codes)


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
    args = parse_program_args_app_codegen_db()
    # Create custom or default configuration.
    init_config(args)
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
    args = parse_program_args_app_codegen_yaml()
    # Create custom or default configuration.
    init_config(args)
    # Run code generator.
    run_code_generation_yaml()
