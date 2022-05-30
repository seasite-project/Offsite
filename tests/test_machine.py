#!/usr/bin/env python3

"""@package test_machine
Tests for function of package machine.
"""

from pathlib import Path
from unittest import TestCase

from offsite.database import close, open_db, rollback
from offsite.descriptions.machine import MachineState


class TestMachine(TestCase):
    def setUp(self):
        self.db_session = open_db('test_database.db')

    def tearDown(self):
        rollback(self.db_session)
        close(self.db_session)

    def test_types(self):
        machine = MachineState.from_yaml(Path('examples/machines/Emmy.yml'), 'gcc')
        machine.to_database(self.db_session)
