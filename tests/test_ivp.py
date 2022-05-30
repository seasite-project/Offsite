#!/usr/bin/env python3

"""@package test_ivp
Tests for function of package ivp.
"""

from pathlib import Path
from unittest import TestCase

from offsite.database import close, open_db, rollback
from offsite.descriptions.ode.ivp import IVP


class TestIVP(TestCase):
    def setUp(self):
        self.db_session = open_db('test_database.db')

    def tearDown(self):
        rollback(self.db_session)
        close(self.db_session)

    def test_types(self):
        ivp = IVP.from_yaml(Path('examples/ivps/Wave1D.ivp'))
        self.assertEquals(ivp.gridSize, 'n')
        with self.assertRaises(RuntimeError):
            IVP.from_database(self.db_session, 'Wave1D')
        ivp.to_database(self.db_session)
        ivp2 = IVP.from_database(self.db_session, 'Wave1D')
        self.assertEquals(ivp, ivp2)
        iid = IVP.database_id(self.db_session, Path('examples/ivps/Wave1D.ivp'))
        self.assertEquals(iid, ivp.db_id)
