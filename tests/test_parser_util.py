#!/usr/bin/env python3

"""@package test_parser_util
Tests for function of package parser_util.
"""

from pathlib import Path
from unittest import TestCase

from offsite.config import ModelToolType
from offsite.descriptions.impl.impl_skeleton import ImplSkeleton, parse_impl_skeletons
from offsite.descriptions.impl.kernel_template import KernelTemplate, parse_kernel_templates
from offsite.descriptions.machine import MachineState, parse_machine_state
from offsite.descriptions.ode.ivp import IVP, parse_ivps
from offsite.descriptions.ode.ode_method import ODEMethod, parse_methods
from offsite.descriptions.parser_util import print_yaml_desc


class TestParserUtil(TestCase):
    def setUp(self):
        self.kernel_templates = [
            KernelTemplate.from_yaml(Path('examples/kernels/pirk/Approx.kernel')),
            KernelTemplate.from_yaml(Path('examples/kernels/pirk/LC.kernel')),
            KernelTemplate.from_yaml(Path('examples/kernels/pirk/RHS.kernel')),
            KernelTemplate.from_yaml(Path('examples/kernels/pirk/RHS_predictor.kernel')),
            KernelTemplate.from_yaml(Path('examples/kernels/pirk/Update.kernel'))]
        self.impl_skeleton = ImplSkeleton.from_yaml(Path('examples/impls/pirk/openmp/A.impl'), self.kernel_templates)
        self.ivp = IVP.from_yaml(Path('examples/ivps/Heat2D.ivp'))
        self.machine = MachineState.from_yaml(Path('examples/machines/node17.yml'), 'icc')
        self.method = ODEMethod.from_yaml(Path('examples/methods/implicit/radauIIA7.ode'))

    def test_types(self):
        # Test function parse_machine.
        parse_machine_state(Path('examples/machines/Emmy.yml'), 'gcc')
        with self.assertRaises(RuntimeError):
            parse_machine_state(Path('examples/machines/Emmy.yml'), 'imaginaryCompiler')
        # Test function parse_kernel_templates.
        templates = parse_kernel_templates(Path('examples/kernels/pirk'), ModelToolType.KERNCRAFT)
        # Test function parse_impl_skeletons.
        parse_impl_skeletons(Path('examples/impls/pirk/openmp'), templates, ModelToolType.KERNCRAFT)
        # Test function parse_methods.
        parse_methods(Path('examples/methods'))
        # Test function parse_ivps.
        parse_ivps(Path('examples/ivps'), None)
        # Test function print_yaml_desc.
        print_yaml_desc(self.machine, [self.impl_skeleton], self.kernel_templates, [self.method], [self.ivp])
