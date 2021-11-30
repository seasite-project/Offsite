#!/usr/bin/env python3

"""@package test_ranking
Tests for function of package ranking.
"""

from unittest import TestCase

from offsite.database.db_mapping import mapping


class TestDbMapping(TestCase):
    def setUp(self):
        pass

    def test_types(self):
        mapping()
