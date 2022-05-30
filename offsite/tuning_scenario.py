"""@package tuning_scenario
Tuning scenario passed to the different offsite applications.

@author: Johannes Seiferth
"""

from collections import OrderedDict
from configparser import RawConfigParser
from sys import version_info
from typing import List, Set

import attr
from pathlib2 import Path

from offsite.config import __impl_skeleton_ext__, __ivp_ext__, __kernel_template_ext__, __ode_method_ext__, \
    __tuning_scenario_ext__
from offsite.solver import Solver, SolverType


class MultiOrderedDict(OrderedDict):
    def __setitem__(self, key, value):
        if isinstance(value, list) and key in self:
            self[key].extend(value)
        else:
            super().__setitem__(key, value)


@attr.s
class TuningScenario:
    """Representation of a TuningScenario object.

    Attributes:
    -----------
    """
    # Tuning scenario parameters and options ...
    # ... for machine configuration.
    machine = attr.ib(type=Path, init=False)
    compiler = attr.ib(type=str, init=False)
    # ... for benchmarks.
    bench_omp_barrier = attr.ib(type=Path, init=False, default=None)
    # ... for solver
    solver = attr.ib(type=Solver, default=Solver('GENERIC_SOLVER'))
    # ... for impl skeletons.
    skeleton_dir = attr.ib(type=Set[Path], default=None)
    skeleton_path = attr.ib(type=Set[Path], default=None)
    # ... for kernel templates.
    template_dir = attr.ib(type=Set[Path], default=None)
    template_path = attr.ib(type=Set[Path], default=None)
    # ... for ODE methods.
    method_dir = attr.ib(type=Set[Path], default=None)
    method_path = attr.ib(type=Set[Path], default=None)
    # ... for IVPs.
    ivp_dir = attr.ib(type=Set[Path], default=None)
    ivp_path = attr.ib(type=Set[Path], default=None)
    # ... for YaskSite.
    ys_foldings = attr.ib(type=List[str], default=list(['']))
    ys_blockings = attr.ib(type=List[str], default=list(['plain']))

    @classmethod
    def from_file(cls, path: Path) -> 'TuningScenario':
        """Construct TuningScenario object from configuration file.

        Parameters:
        -----------
        path: Path
            Path to configuration file. Required file ending: '.scenario'.

        Returns:
        --------
        Config
            Created Config file.
        """
        # Create default Config object.
        scenario_obj: TuningScenario = cls()
        # Check file ending of tuning scenario file.
        if path.suffix != __tuning_scenario_ext__:
            raise RuntimeError('Expected file ending: \'{}\' for tuning scenario! Found ending \'{}\' instead'.format(
                __tuning_scenario_ext__, path.suffix))
        # Read config file.
        parser = RawConfigParser(dict_type=MultiOrderedDict, strict=False)
        # Check kerncraft version.
        if version_info[1] > 6:
            parser.read(path)
        else:
            parser.read(str(path))
        # Parse configuration options for...
        tag_dir = 'dir'
        tag_path = 'path'
        # ... machine configuration.
        tag_machine = 'MACHINE'
        if tag_machine in parser:
            if tag_path in parser[tag_machine]:
                scenario_obj.machine = Path(parser[tag_machine][tag_path])
                if scenario_obj.machine.suffix not in ['.yml', '.yaml']:
                    raise RuntimeError(
                        'Expected file ending: \'{}\' for machine! Found ending \'{}\' instead.'.format(
                            '.yml', scenario_obj.machine.suffix))
            else:
                raise RuntimeError('Missing entry:\n[MACHINE]\npath = MY_CUSTOM_PATH')
            tag_compiler = 'compiler'
            if tag_compiler in parser[tag_machine]:
                scenario_obj.compiler = parser[tag_machine][tag_compiler]
            else:
                raise RuntimeError('Missing entry:\n[MACHINE]\ncompiler = MY_CUSTOM_COMPILER')
        else:
            raise RuntimeError('Missing entry:\n[MACHINE]\npath = MY_CUSTOM_PATH\ncompiler = MY_CUSTOM_COMPILER')
        # ... benchmark.
        tag_benchmark = 'BENCHMARK'
        if tag_benchmark in parser:
            tag_omp_barrier = 'omp_barrier'
            if tag_omp_barrier in parser[tag_benchmark]:
                scenario_obj.bench_omp_barrier = Path(parser[tag_benchmark][tag_omp_barrier])
        # ... solver.
        tag_solver = 'SOLVER'
        if tag_solver in parser:
            tag_name = 'name'
            tag_type = 'type'
            if bool(tag_name in parser[tag_solver]) != bool(tag_type in parser[tag_solver]):
                raise RuntimeError('Missing entry. Both entries are required:\n[SOLVER]\nname = MY_CUSTOM_SOLVER_NAME'
                                   '\ntype = MY_CUSTOM_SOLVER_TYPE')
            solver_name: str = 'GENERIC_SOLVER'
            if tag_name in parser[tag_solver]:
                solver_name = parser[tag_solver][tag_name]
            if tag_type in parser[tag_solver]:
                scenario_obj.solver = Solver.make_solver(solver_name, parser[tag_solver][tag_type])
        # ... impl skeleton.
        tag_skeleton = 'IMPLEMENTATION SKELETON'
        if tag_skeleton in parser:
            if tag_dir in parser[tag_skeleton]:
                scenario_obj.skeleton_dir = set(Path(x) for x in parser[tag_skeleton][tag_dir].split('\n'))
                if any((not x.is_dir() for x in scenario_obj.skeleton_dir)):
                    raise RuntimeError('Implementation skeleton directory does not exist!')
            if tag_path in parser[tag_skeleton]:
                scenario_obj.skeleton_path = set(Path(x) for x in parser[tag_skeleton][tag_path].split('\n'))
                if any((x.suffix != __impl_skeleton_ext__ for x in scenario_obj.skeleton_path)):
                    raise RuntimeError(
                        'Expected file ending: \'{}\' for impl skeleton!'.format(__impl_skeleton_ext__))
        if scenario_obj.skeleton_dir is None and scenario_obj.skeleton_path is None:
            raise RuntimeError('Section [IMPLEMENTATION SKELETON] requires at least one of the following entries:\n'
                               '\"dir = MY_CUSTOM_PATH\" or \"path = \"MY_CUSTOM_PATH\"')
        # ... kernel template.
        tag_kernel = 'KERNEL TEMPLATE'
        if tag_kernel in parser:
            if tag_dir in parser[tag_kernel]:
                scenario_obj.template_dir = set(Path(x) for x in parser[tag_kernel][tag_dir].split('\n'))
                if any((not x.is_dir() for x in scenario_obj.template_dir)):
                    raise RuntimeError('Kernel template directory does not exist!')
            if tag_path in parser[tag_kernel]:
                scenario_obj.template_path = set(Path(x) for x in parser[tag_kernel][tag_path].split('\n'))
                if any((x.suffix != __kernel_template_ext__ for x in scenario_obj.template_path)):
                    raise RuntimeError(
                        'Expected file ending: \'{}\' for kernel template!'.format(__kernel_template_ext__))
        if scenario_obj.template_dir is None and scenario_obj.template_path is None:
            raise RuntimeError('Section [KERNEL TEMPLATE] requires at least one of the following entries:\n'
                               '\"dir = MY_CUSTOM_PATH\" or \"path = \"MY_CUSTOM_PATH\"')
        # ... ODE method.
        tag_method = 'ODE METHOD'
        if tag_method in parser:
            if tag_dir in parser[tag_method]:
                scenario_obj.method_dir = set(Path(x) for x in parser[tag_method][tag_dir].split('\n'))
                if any((not x.is_dir() for x in scenario_obj.method_dir)):
                    raise RuntimeError('ODE method directory does not exist!')
            if tag_path in parser[tag_method]:
                scenario_obj.method_path = set(Path(x) for x in parser[tag_method][tag_path].split('\n'))
                if any((x.suffix != __ode_method_ext__ for x in scenario_obj.method_path)):
                    raise RuntimeError('Expected file ending: \'{}\' for ODE method!'.format(__ode_method_ext__))
        if scenario_obj.solver.type == SolverType.ODE and scenario_obj.method_dir is None and scenario_obj.method_path is None:
            raise RuntimeError('Section [ODE METHOD] requires at least one of the following entries:\n'
                               '\"dir = MY_CUSTOM_PATH\" or \"path = \"MY_CUSTOM_PATH\"')
        # ... IVP.
        tag_ivp = 'IVP'
        if tag_ivp in parser:
            if tag_dir in parser[tag_ivp]:
                scenario_obj.ivp_dir = set(Path(x) for x in parser[tag_ivp][tag_dir].split('\n'))
                if any((not x.is_dir() for x in scenario_obj.ivp_dir)):
                    raise RuntimeError('IVP directory does not exist!')
            if tag_path in parser[tag_ivp]:
                scenario_obj.ivp_path = set(Path(x) for x in parser[tag_ivp][tag_path].split('\n'))
                if any((x.suffix != __ivp_ext__ for x in scenario_obj.ivp_path)):
                    raise RuntimeError('Expected file ending: \'{}\' for IVP!'.format(__ivp_ext__))
        if scenario_obj.solver.type == SolverType.ODE and scenario_obj.ivp_dir is None and scenario_obj.ivp_path is None:
            raise RuntimeError('Section [IVP] requires at least one of the following entries:\n'
                               '\"dir = MY_CUSTOM_PATH\" or \"path = \"MY_CUSTOM_PATH\"')
        # ... YaskSite.
        tag_ys = 'YASKSITE'
        if tag_ys in parser:
            tag_folding = 'folding'
            if tag_folding in parser[tag_ys]:
                fold_str = parser[tag_ys][tag_folding]
                scenario_obj.ys_foldings = [f.replace("''", "").strip() for f in fold_str.split(';') if f]
            tag_blocking = 'blocking'
            if tag_blocking in parser[tag_ys]:
                block_str = parser[tag_ys][tag_blocking]
                scenario_obj.ys_blockings = [b.replace("''", "").strip() for b in block_str.split(';') if b]
        # Return created tuning scenario object.
        return scenario_obj
