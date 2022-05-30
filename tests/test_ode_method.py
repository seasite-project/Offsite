#!/usr/bin/env python3

"""@package test_ode_method
Tests for function of package ode_method.
"""

from pathlib import Path
from unittest import TestCase

from offsite.database import close, open_db, rollback
from offsite.descriptions.ode.ode_method import ODEMethod


class TestODEMethod(TestCase):
    def setUp(self):
        self.db_session = open_db('test_database.db')

    def tearDown(self):
        rollback(self.db_session)
        close(self.db_session)

    def test_types(self):
        method = ODEMethod.from_yaml(Path('examples/methods/implicit/implicitEuler.ode'))
        self.assertEquals(method.name, 'implicitEuler')
        self.assertEquals(method.stages, 1)
        self.assertEquals(method.order_, 1)
        self.assertEquals(method.correctorSteps, 0)
        self.assertEquals(method.coefficientsA, [["1.0"]])
        self.assertEquals(method.coefficientsB, ["1.0"])
        self.assertEquals(method.coefficientsC, ["1.0"])
        with self.assertRaises(RuntimeError):
            ODEMethod.from_database(self.db_session, 'implicitEuler')
        method.to_database(self.db_session)
        method2 = ODEMethod.from_database(self.db_session, 'implicitEuler')
        self.assertEquals(method, method2)
        mid = ODEMethod.database_id(
            self.db_session, 'examples/methods/implicit/implicitEuler.ode')
        self.assertEquals(mid, method.db_id)
