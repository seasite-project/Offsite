#!/usr/bin/env python3

"""@package test_kerncraft
Tests for kerncraft's ECM modes for kernel prediction.
"""

from pathlib import Path
from subprocess import run
from unittest import TestCase

import attr

import offsite.config
from offsite.codegen.codegen_util import write_codes_to_file
from offsite.codegen.generator.kerncraft_generator import KerncraftCodeGenerator
from offsite.config import IncoreToolType, ModelToolType, init_config
from offsite.descriptions.impl.kernel_template import KernelTemplate
from offsite.descriptions.machine import MachineState
from offsite.descriptions.ode.ivp import IVP
from offsite.descriptions.ode.ode_method import ODEMethod
from offsite.train.node.util.kerncraft_utils import execute_kerncraft_ecm_mode, parse_kerncraft_output_ecm_mode, \
    execute_kerncraft_bench_mode, parse_kerncraft_output_bench_mode, execute_kerncraft_lc_mode, \
    parse_kerncraft_output_lc_mode
from offsite.tuning_scenario import TuningScenario


@attr.s
class Helper:
    """Helper class to mimique args.Namespace."""
    incore = attr.ib(type=IncoreToolType, default=IncoreToolType.OSACA)
    config = attr.ib(type=Path, default=None)

class TestKerncraft(TestCase):
    def setUp(self):
        # Create custom or default configuration.
        args = Helper()
        init_config(args)
        self.config = offsite.config.offsiteConfig
        self.config.scenario = TuningScenario.from_file(Path('tests/data/test_kerncraft_prediction.scenario'))
        self.machine = MachineState.from_yaml(self.config.scenario.machine, self.config.scenario.compiler)
        self.ivps = [IVP.from_yaml(p) for p in self.config.scenario.ivp_path]
        self.methods = [ODEMethod.from_yaml(p) for p in self.config.scenario.method_path]
        self.templates = [KernelTemplate.from_yaml(p) for p in self.config.scenario.template_path]

    def test_types(self):
        folder = Path('tests/tmp')
        for template in self.templates:
            for method in self.methods:
                if template.isIVPdependent:
                    for ivp in self.ivps:
                        for kernel in template.variants:
                            paths = list()
                            if template.modelTool == ModelToolType.KERNCRAFT:
                                if ivp.name in ['InverterChain', 'Medakzo']:
                                    codes = KerncraftCodeGenerator().generate(kernel, method, ivp)
                                    paths = write_codes_to_file(codes, folder)
                            for path in paths:
                                # Test ECM mode.
                                expected_results = self.results_ecm[str(path)][method.name][ivp.name]
                                self._test_ecm_mode(path, method, ivp, expected_results)
                                # Test bench mode.
                                # expected_results = self.results_bench[str(path)][method.name][ivp.name]
                                # self._test_bench_mode(path, method, ivp, expected_results)
                                # Test LC mode.
                                expected_results = self.results_lc[str(path)][method.name][ivp.name]
                                self._test_lc_mode(path, method, ivp, expected_results)
                else:
                    for kernel in template.variants:
                        paths = list()
                        if template.modelTool == ModelToolType.KERNCRAFT:
                            codes = KerncraftCodeGenerator().generate(kernel, method, None)
                            paths = write_codes_to_file(codes, folder)
                        for path in paths:
                            # Test ECM mode.
                            expected_results = self.results_ecm[str(path)][method.name][None]
                            self._test_ecm_mode(path, method, None, expected_results)
                            # Test bench mode.
                            expected_results = self.results_bench[str(path)][method.name][None]
                            # self._test_bench_mode(path, method, None, expected_results)
                            # Test LC mode.
                            expected_results = self.results_lc[str(path)][method.name][None]
                            self._test_lc_mode(path, method, None, expected_results)

    results_ecm = {
        'tests/tmp/Update_j.c': {'radauIIA7': {None: {
            IncoreToolType.IACA: {1: 12.0, 2: 6.0, 3: 4.0, 4: 3.0, 5: 2.4, 6: 2.0, 7: 1.7, 8: 1.5},
            IncoreToolType.OSACA: {1: 12.0, 2: 6.0, 3: 4.0, 4: 3.0, 5: 2.4, 6: 2.0, 7: 1.7, 8: 1.5}}}},
        'tests/tmp/RHS_jl_InverterChain.c': {'radauIIA7': {'InverterChain': {
            IncoreToolType.IACA: {1: 281.5, 2: 142.1, 3: 97.4, 4: 76.3, 5: 64.3, 6: 56.8, 7: 51.8, 8: 48.1},
            IncoreToolType.OSACA: {1: 224.0, 2: 113.7, 3: 79.2, 4: 63.3, 5: 54.5, 6: 49.0, 7: 45.4, 8: 42.8}}}},
        'tests/tmp/RHS_lj_InverterChain.c': {'radauIIA7': {'InverterChain': {
            IncoreToolType.IACA: {1: 56.0, 2: 28.4, 3: 19.8, 4: 15.8, 5: 13.6, 6: 12.3, 7: 11.3, 8: 10.7},
            IncoreToolType.OSACA: {1: 56.0, 2: 28.4, 3: 19.8, 4: 15.8, 5: 13.6, 6: 12.3, 7: 11.3, 8: 10.7}}}}}

    results_bench = {
        'tests/tmp/Update_j.c': {'radauIIA7': {None: {IncoreToolType.IACA: 12.0, IncoreToolType.OSACA: 1.5}}},
        'tests/tmp/RHS_jl_InverterChain.c': {
            'radauIIA7': {'InverterChain': {IncoreToolType.IACA: 281.5, IncoreToolType.OSACA: 48.1}}},
        'tests/tmp/RHS_lj_InverterChain.c': {
            'radauIIA7': {'InverterChain': {IncoreToolType.IACA: 56.0, IncoreToolType.OSACA: 10.7}}}}

    results_lc = {
        'tests/tmp/Update_j.c': {'radauIIA7': {None: {
            IncoreToolType.IACA: {'L1': [('n', 2048)], 'L2': [('n', 16384)], 'L3': [('n', 4194304)]},
            IncoreToolType.OSACA: {'L1': [('n', 2048)], 'L2': [('n', 16384)], 'L3': [('n', 4194304)]}}}},
        'tests/tmp/RHS_jl_InverterChain.c': {'radauIIA7': {'InverterChain': {
            IncoreToolType.IACA: {
                'L1': [('2*n*s + s', 4096), ('8*n', 4095), ('2*n', 1025)],
                'L2': [('2*n*s + s', 32768), ('8*n', 32767), ('2*n', 8193)],
                'L3': [('2*n*s + s', 8388608), ('8*n', 8388607), ('2*n', 2097153)]},
            IncoreToolType.OSACA: {
                'L1': [('2*n*s + s', 4096), ('8*n', 4095), ('2*n', 1025)],
                'L2': [('2*n*s + s', 32768), ('8*n', 32767), ('2*n', 8193)],
                'L3': [('2*n*s + s', 8388608), ('8*n', 8388607), ('2*n', 2097153)]}}}},
        'tests/tmp/RHS_lj_InverterChain.c': {'radauIIA7': {'InverterChain': {
            IncoreToolType.IACA: {
                'L1': [('2*n*s + s', 4096)], 'L2': [('2*n*s + s', 32768)], 'L3': [('2*n*s + s', 8388608)]},
            IncoreToolType.OSACA: {
                'L1': [('2*n*s + s', 4096)], 'L2': [('2*n*s + s', 32768)], 'L3': [('2*n*s + s', 8388608)]}}}}}

    def _test_ecm_mode(self, path, method, ivp, expected_results):
        self.config.pred_incore_tool = IncoreToolType.IACA
        out_iaca = execute_kerncraft_ecm_mode(path, self.machine, method, ivp, 2000000)
        ecm_iaca = parse_kerncraft_output_ecm_mode(out_iaca)
        self.assertDictEqual(ecm_iaca, expected_results[self.config.pred_incore_tool])
        self.config.pred_incore_tool = IncoreToolType.OSACA
        out_osaca = execute_kerncraft_ecm_mode(path, self.machine, method, ivp, 2000000)
        ecm_osaca = parse_kerncraft_output_ecm_mode(out_osaca)
        self.assertDictEqual(ecm_osaca, expected_results[self.config.pred_incore_tool])
        # if ecm_iaca != ecm_osaca:
        #    print('\n#####################################################################\n')
        #    print(path)
        #    print('-------------- IACA:\n')
        #    print(out_iaca)
        #    print('-------------- OSACA:\n')
        #    print(out_osaca)

    def _test_bench_mode(self, path, method, ivp, expected_results):
        self.config.pred_incore_tool = IncoreToolType.IACA
        out_iaca = execute_kerncraft_bench_mode(path, self.machine, method, ivp, 2000000, self.machine.coresPerSocket)
        ecm_iaca = parse_kerncraft_output_bench_mode(out_iaca)
        self.assertDictEqual(ecm_iaca, expected_results[self.config.pred_incore_tool])
        self.config.pred_incore_tool = IncoreToolType.OSACA
        out_osaca = execute_kerncraft_bench_mode(path, self.machine, method, ivp, 2000000, self.machine.coresPerSocket)
        ecm_osaca = parse_kerncraft_output_bench_mode(out_osaca)
        self.assertDictEqual(ecm_osaca, expected_results[self.config.pred_incore_tool])
        # if ecm_iaca != ecm_osaca:
        #    print('\n#####################################################################\n')
        #    print(path)
        #    print('-------------- IACA:\n')
        #    print(out_iaca)
        #    print('-------------- OSACA:\n')
        #    print(out_osaca)

    def _test_lc_mode(self, path, method, ivp, expected_results):
        self.config.pred_incore_tool = IncoreToolType.IACA
        out_iaca = execute_kerncraft_lc_mode(path, self.machine, method, ivp)
        ecm_iaca = parse_kerncraft_output_lc_mode(out_iaca)
        self.assertDictEqual(ecm_iaca, expected_results[self.config.pred_incore_tool])
        self.config.pred_incore_tool = IncoreToolType.OSACA
        out_osaca = execute_kerncraft_lc_mode(path, self.machine, method, ivp)
        ecm_osaca = parse_kerncraft_output_lc_mode(out_osaca)
        self.assertDictEqual(ecm_osaca, expected_results[self.config.pred_incore_tool])

    def tearDown(self):
        run(['rm', '-rf', 'tests/tmp/'])
