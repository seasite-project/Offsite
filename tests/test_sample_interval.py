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
    create_samples_upper_border_working_set, create_samples_border_memory_lvl, create_samples_memory_lvl


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
        lb_ws, next_point = create_samples_lower_border_working_set(400, 1000, 5, None)
        self.assertEqual(next_point, 520)
        points: List[float] = [s.sample for s in lb_ws]
        self.assertListEqual(points, [412, 436, 460, 484, 508])
        lb_ws_ivp, next_point_ivp = create_samples_lower_border_working_set(400, 1000, 7, self.ivps[0])
        self.assertEqual(next_point_ivp, 520)
        points_ivp: List[float] = [s.sample for s in lb_ws_ivp]
        self.assertListEqual(points_ivp, [408, 426, 443, 460, 477, 494, 511])
        # Test create_samples_upper_border_working_set.
        ub_ws, prev_point = create_samples_upper_border_working_set(400, 1000, 5, None)
        self.assertEqual(prev_point, 880)
        points2: List[float] = [s.sample for s in ub_ws]
        self.assertListEqual(points2, [892, 916, 940, 964, 988])
        ub_ws_ivp, prev_point_ivp = create_samples_upper_border_working_set(400, 1000, 7, self.ivps[0])
        self.assertEqual(prev_point_ivp, 880)
        points_ivp2: List[float] = [s.sample for s in ub_ws_ivp]
        self.assertListEqual(points_ivp2, [889, 906, 923, 940, 957, 974, 992])
        # Test create_samples_border_memory_lvl.
        border_mem_lvl, next_point2 = create_samples_border_memory_lvl(40000, 1000000, 5, None)
        self.assertEqual(next_point2, 50000)
        points3: List[float] = [s.sample for s in border_mem_lvl]
        self.assertListEqual(points3, [41000, 43000, 45000, 47000, 49000])
        border_mem_lvl_ivp, next_point_ivp2 = create_samples_border_memory_lvl(40000, 100000000, 8, self.ivps[0])
        self.assertEqual(next_point_ivp2, 60000)
        points3_ivp: List[float] = [s.sample for s in border_mem_lvl_ivp]
        self.assertListEqual(points3_ivp, [41250, 43750, 46250, 48750, 51250, 53750, 56250, 58750])
        # Test create_samples_memory_lvl.
        mem_lvl = create_samples_memory_lvl(40000, 5, 1500, None)
        points4: List[float] = [s.sample for s in mem_lvl]
        self.assertListEqual(points4, [30020000, 89980002, 149940004, 209900006, 269860008])
        mem_lvl_ivp = create_samples_memory_lvl(40000, 7, 1500, self.ivps[0])
        points_ivp4: List[float] = [s.sample for s in mem_lvl_ivp]
        self.assertListEqual(points_ivp4, [30020000, 89980002, 149940004, 209900006, 269860008, 329820010, 389780012])
