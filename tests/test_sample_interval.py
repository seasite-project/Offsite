#!/usr/bin/env python3

"""@package test_sample_interval
Tests for SampleInterval and its utility functions.
"""

from pathlib import Path
from typing import List
from unittest import TestCase

import attr

import offsite.config
from offsite.config import ProgramModeType, init_config
from offsite.descriptions.ode.ivp import IVP
from offsite.tuning_scenario import TuningScenario
from offsite.util.sample_interval import SampleInterval, SampleType, create_samples_lower_border_working_set, \
    create_samples_upper_border_working_set, create_samples_memory_lvl


@attr.s
class Helper:
    """Helper class to mimique args.Namespace."""
    # incore = attr.ib(type=IncoreToolType, default=IncoreToolType.OSACA)
    mode = attr.ib(type=ProgramModeType, default=ProgramModeType.MODEL)
    config = attr.ib(type=Path, default=None)


class TestSampleInterval(TestCase):
    def setUp(self):
        # Create custom or default configuration.
        args = Helper()
        init_config(args)
        self.config = offsite.config.offsiteConfig
        self.config.scenario = TuningScenario.from_file(Path('tests/data/test_kerncraft_prediction.scenario'))
        self.ivps = [IVP.from_yaml(p) for p in self.config.scenario.ivp_path]

    def test_types(self):
        intv1: SampleInterval = SampleInterval(0, 100, SampleType.MODEL_INNER, 40)
        median1: int = intv1.median()
        self.assertEqual(median1, 50)
        intv2: SampleInterval = SampleInterval(0, 77, SampleType.MODEL_INNER, 40)
        median2: int = intv2.median()
        self.assertEqual(median2, 38)
        # Test create_samples_lower_border_working_set.
        lb_ws, next_point = create_samples_lower_border_working_set(400, 1000, 15, 2.3, None)
        self.assertEqual(next_point, 1001)
        points: List[float] = [s.sample for s in lb_ws]
        self.assertListEqual(points, [486, 660, 834, 961])
        lb_ws_ivp, next_point_ivp = create_samples_lower_border_working_set(400, 1000, 15, 2.3, self.ivps[0])
        self.assertEqual(next_point_ivp, 1001)
        points_ivp: List[float] = [s.sample for s in lb_ws_ivp]
        self.assertListEqual(points_ivp, [486, 660, 834, 961])
        # Test create_samples_upper_border_working_set.
        ub_ws, prev_point = create_samples_upper_border_working_set(400, 1000, 5, 3.3, None)
        self.assertEqual(prev_point, 399)
        points2: List[float] = [s.sample for s in ub_ws]
        self.assertListEqual(points2, [486, 660, 834, 961])
        ub_ws_ivp, prev_point_ivp = create_samples_upper_border_working_set(400, 1000, 5, 2.3, self.ivps[0])
        self.assertEqual(prev_point_ivp, 399)
        points_ivp2: List[float] = [s.sample for s in ub_ws_ivp]
        self.assertListEqual(points_ivp2, [486, 660, 834, 961])
        # Test create_samples_memory_lvl.
        create_samples_memory_lvl()
