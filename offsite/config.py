"""@package config
Configuration options of the off_tune application.

@author: Johannes Seiferth
"""

from configparser import ConfigParser
from enum import Enum
from sys import version_info

import attr
from argparse import Namespace

from offsite.tuning_scenario import TuningScenario


class ProgramModeType(Enum):
    """Defines how the performance of implementation variants is rated.

    - MODEL
        Performance prediction using the ECM model.
    - RUN
        Runtime sampling of variants.
    """
    MODEL = 'MODEL'
    RUN = 'RUN'


class ModelToolType(Enum):
    """Defines which tool is used for performance modelling.

    - KERNCRAFT
        Performance prediction using the KERNCRAFT tool.
    - YASKSITE
        Performance prediction using the YASKSITE tool.
    """
    KERNCRAFT = 'kerncraft'
    YASKSITE = 'yasksite'


class IncoreToolType(Enum):
    """Defines what tool is used to determine the in-core model.

    - IACA
        IACA tool.
    - OSACA
        OSACA tool.
    - LLVM-MCA
        LLVM-MCA tool.
    """
    IACA = 'IACA'
    OSACA = 'OSACA'
    LLVMMCA = 'LLVM-MCA'


@attr.s
class Config:
    """Representation of a Config object.

    Attributes:
    -----------
    """
    # Program arguments
    args = attr.ib(type=Namespace)
    scenario = attr.ib(type=TuningScenario, default=None)
    # Configuration options ...
    # ... for working sets
    available_cache_size = attr.ib(type=float, default=0.99)
    samples_per_border = attr.ib(type=int, default=2)
    samples_per_interval = attr.ib(type=int, default=2)
    samples_border_region_memory_lvl = attr.ib(type=int, default=5)
    samples_memory_lvl = attr.ib(type=int, default=5)
    memory_lvl_sample_offset = attr.ib(type=float, default=1.015)
    # ... for layer condition analysis
    layer_condition_safety_margin = attr.ib(type=float, default=2.0)
    # ... for benchmarking
    repetitions_communication_operations = attr.ib(type=int, default=100)
    repetitions_implementation_variants = attr.ib(type=int, default=10)
    # ... for code generation
    var_idx = attr.ib(type=str, default='j')
    yasksite_stencil_dir = attr.ib(type=str, default='examples/ivps/yasksite_stencils')
    # ... for LIKWID
    likwid_set_frequencies = attr.ib(type=str, default='likwid-setFrequencies')
    # ... for predictions
    pred_model_tool = attr.ib(type=ModelToolType, default=ModelToolType.KERNCRAFT)
    pred_incore_tool = attr.ib(type=IncoreToolType, default=IncoreToolType.IACA)

    @classmethod
    def from_file(cls, args: Namespace) -> 'Config':
        """Construct Config object from configuration file.

        Parameters:
        -----------
        args: argparse.Namespace
            Passed program arguments.

        Returns:
        --------
        Config
            Created Config file.
        """
        # Create default Config object.
        config_obj: Config = cls(args)
        # Read config file.
        path = args.config
        parser = ConfigParser()
        # Check kerncraft version.
        if version_info[1] > 6:
            parser.read(path)
        else:
            parser.read(str(path))
        # Program arguments
        config_obj.args = args
        if 'scenario' in vars(args):
            config_obj.scenario = TuningScenario.from_file(args.scenario)
        # Parse configuration options for...
        # ... working sets.
        tag_ws = 'WORKING SETS'
        if tag_ws in parser:
            tag_cache = 'AvailableCacheSize'
            if tag_cache in parser[tag_ws]:
                config_obj.available_cache_size = float(parser[tag_ws][tag_cache])
            tag_samples = 'SamplesPerBorderRegion'
            if tag_samples in parser[tag_ws]:
                config_obj.samples_per_border = int(parser[tag_ws][tag_samples])
            tag_interval = 'SamplesPerInterval'
            if tag_interval in parser[tag_ws]:
                config_obj.samples_per_interval = int(parser[tag_ws][tag_interval])
            tag_memory_offset = 'MemoryLvlSampleOffset'
            if tag_memory_offset in parser[tag_ws]:
                config_obj.memory_lvl_sample_offset = float(parser[tag_ws][tag_memory_offset])
            tag_memory_border_samples = 'SamplesBorderRegionToMemoryLvl'
            if tag_memory_border_samples in parser[tag_ws]:
                config_obj.samples_border_region_memory_lvl = int(parser[tag_ws][tag_memory_border_samples])
            tag_memory_samples = 'SamplesMemoryLvl'
            if tag_memory_samples in parser[tag_ws]:
                config_obj.samples_memory_lvl = int(parser[tag_ws][tag_memory_samples])
        # ... layer condition analysis.
        tag_lc = 'LAYER CONDITION'
        if tag_lc in parser:
            tag_margin = 'SafetyMargin'
            if tag_margin in parser[tag_lc]:
                config_obj.layer_condition_safety_margin = float(parser[tag_lc][tag_margin])
        # ... benchmarking.
        tag_bench = 'BENCHMARK'
        if tag_bench in parser:
            tag_comm = 'RepetitionsCommunicationOperations'
            if tag_comm in parser[tag_bench]:
                config_obj.repetitions_communication_operations = int(parser[tag_bench][tag_comm])
            tag_impl = 'RepetitionsImplementationVariants'
            if tag_impl in parser[tag_bench]:
                config_obj.repetitions_implementation_variants = int(parser[tag_bench][tag_impl])
        # ... code generation.
        tag_codegen = 'CODEGEN'
        if tag_codegen in parser:
            tag_idx = 'idx'
            if tag_idx in parser[tag_codegen]:
                config_obj.var_idx = parser[tag_codegen][tag_idx]
            tag_yasksite_stencil_dir = 'yasksite_stencil_dir'
            if tag_yasksite_stencil_dir in parser[tag_codegen]:
                config_obj.yasksite_stencil_dir = parser[tag_codegen][tag_yasksite_stencil_dir]
        # ... LIKWID.
        tag_likwid = 'LIKWID'
        if tag_likwid in parser:
            tag_freq = 'set-frequencies'
            if tag_freq in parser[tag_freq]:
                config_obj.likwid_set_frequencies = parser[tag_likwid][tag_freq]
        # .. prediction tools.
        tag_pred = 'PREDICTION'
        if tag_pred in parser:
            tag_tool = 'ModelingTool'
            if tag_tool in parser[tag_pred]:
                tool_str = parser[tag_pred][tag_tool]
                config_obj.pred_model_tool = ModelToolType[tool_str]
            tag_incore = 'IncoreModel'
            if tag_incore in parser[tag_pred]:
                incore_str = parser[tag_pred][tag_incore]
                config_obj.pred_incore_tool = IncoreToolType[incore_str]
        # Return created Config object.
        return config_obj


offsiteConfig = None


def init_config(args):
    global offsiteConfig
    if args and args.config:
        offsiteConfig = Config.from_file(args)
    else:
        offsiteConfig = Config(args)
