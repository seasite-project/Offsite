#!/usr/bin/env python3

"""@package test_performance_model
Tests for package performance_model.
"""

from unittest import TestCase

from offsite.train.node.util.performance_model import compute_pmodel_kernel_pred, compute_kernel_pred, \
    compute_impl_variant_pred


class TestPerformanceModel(TestCase):
    def test_types(self):
        # Test function compute_kernel_pred.
        self.assertEquals(float(compute_pmodel_kernel_pred(8, 2.0, 8)), 2.0)
        self.assertEquals(float(compute_pmodel_kernel_pred(3, 2.0, 8)), 0.75)
        # Test function compute_kernel_pred.
        self.assertAlmostEquals(float(compute_kernel_pred([3.0], 3400)), 0.000882352941176471)
        self.assertAlmostEquals(float(compute_kernel_pred([4.0, 4.5, 2.3], 34000)), 0.000317647058823529)
        # Test function compute_impl_variant_pred.
        self.assertAlmostEquals(float(compute_impl_variant_pred({1: (4.2, 1.0)}, {1: 400}, 2.14)), 1682.140)
        self.assertAlmostEquals(float(compute_impl_variant_pred(
            {1: (4.7, 1.0), 3: (1.1, 2.1)}, {1: 400, 3: 200}, 22.54)), 2364.5400000)
        with self.assertRaises(RuntimeError):
            compute_impl_variant_pred({1: (4.2, 1.0)}, {1: 400, 3: 200}, 2.14)
        with self.assertRaises(RuntimeError):
            compute_impl_variant_pred({1: (4.2, 1.0), 3: (1.1, 2.1)}, {1: 400}, 2.14)
