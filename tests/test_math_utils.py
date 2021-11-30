#!/usr/bin/env python3

"""@package test_math_utils
Tests for math utility functions.
"""

from typing import List
from unittest import TestCase

from offsite.util.math_utils import eval_math_expr, percent_deviation, remove_outliers, solve_equation


class TestMathUtils(TestCase):
    def test_types(self):
        # Test eval_math_expr.
        t1: int = eval_math_expr('2*4')
        self.assertEqual(t1, 8)
        t2: int = eval_math_expr('2*n', [('n', 4)])
        self.assertEqual(t2, 8)
        t3: float = eval_math_expr('4*n', [('n', 2)], cast_to=float)
        self.assertEqual(t3, 8.0)
        t4: str = eval_math_expr('2*n+s', [('n', 2)], cast_to=str)
        self.assertEqual(t4, 's + 4.0')
        # Test percent_deviation.
        pd: float = percent_deviation(1.2, 1.4)
        self.assertAlmostEqual(pd, -14.285714285714283)
        # Test remove_outliers.
        ro: List[float] = remove_outliers([-29.1, 4.1, 4.3, 4.2, 4.6, 4.9, 4.2, 4.3, 4.1, 4.0, 3.7, 5.2, 39.0, 38.0])
        self.assertEqual(len(ro), 11)
        self.assertListEqual(ro, [4.1, 4.3, 4.2, 4.6, 4.9, 4.2, 4.3, 4.1, 4.0, 3.7, 5.2])
        # Test solve_equation.
        eq1: float = solve_equation('y + 2', '11', 'y')
        self.assertEqual(eq1[0], 9)
        eq2: float = solve_equation('y + 2', 'x', 'y', [('x', 42)])
        self.assertEqual(eq2[0], 40)
        eq3: float = solve_equation('y + 2', '3', 'y', [('y', 42)])
        self.assertEqual(len(eq3), 0)
