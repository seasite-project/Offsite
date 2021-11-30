#!/usr/bin/env python3

"""@package test_database
Delete test database after all tests were executed.
"""

from os import remove
from unittest import TestCase


class TestDbRemoval(TestCase):
    def tearDown(self):
        remove('test_database.db')

    def test_types(self):
        pass
