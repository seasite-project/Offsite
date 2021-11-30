#!/usr/bin/env python3

"""@package test_database
Tests for function of package database.
"""

from pathlib import Path
from unittest import TestCase

import attr

import offsite.config
from offsite.database import close, commit, open_db, rollback
from offsite.descriptions.impl.impl_skeleton import ImplSkeleton
from offsite.descriptions.impl.kernel_template import KernelTemplate
from offsite.descriptions.machine import MachineState
from offsite.descriptions.ode.ivp import IVP
from offsite.descriptions.ode.ode_method import ODEMethod


@attr.s
class Helper:
    """Helper class to mimique args.Namespace."""
    verbose = attr.ib(type=bool, default=False)
    config = attr.ib(type=Path, default=None)


class TestDatabase(TestCase):
    def setUp(self):
        # Create custom or default configuration.
        args = Helper()
        offsite.config.init_config(args)
        self.config = offsite.config.offsiteConfig
        self.machine = MachineState.from_yaml(Path('examples/machines/BroadwellEP_E5-2630v4.yml'), 'icc')
        self.machine2 = MachineState.from_yaml(Path('examples/machines/Emmy.yml'), 'icc')
        self.ivp = IVP.from_yaml(Path('examples/ivps/Heat2D.ivp'))
        self.ivp2 = IVP.from_yaml(Path('examples/ivps/InverterChain.ivp'))
        self.method = ODEMethod.from_yaml(Path('examples/methods/implicit/radauIIA7.ode'))
        self.method2 = ODEMethod.from_yaml(Path('examples/methods/implicit/lobattoIIIC8.ode'))
        self.kernel_templates = [
            KernelTemplate.from_yaml(Path('examples/kernels/pirk/Approx.kernel')),
            KernelTemplate.from_yaml(Path('examples/kernels/pirk/LC.kernel')),
            KernelTemplate.from_yaml(Path('examples/kernels/pirk/RHS.kernel')),
            KernelTemplate.from_yaml(Path('examples/kernels/pirk/RHS_predictor.kernel')),
            KernelTemplate.from_yaml(Path('examples/kernels/pirk/Update.kernel'))]
        self.impl_skeleton = ImplSkeleton.from_yaml(
            Path('examples/impls/pirk/openmp/A.impl'), self.kernel_templates)

    def test_types(self):
        # Test function open_db.
        db_session = open_db('test_database.db')
        rollback(db_session)
        # Test function insert.
        self.machine.to_database(db_session)
        self.machine.to_database(db_session)
        self.ivp.to_database(db_session)
        self.ivp2.to_database(db_session)
        self.method.to_database(db_session)
        self.method2.to_database(db_session)
        for template in self.kernel_templates:
            template.to_database(db_session)
        self.impl_skeleton.to_database(db_session)
        # Test function commit.
        commit(db_session)
        # Test function rollback.
        rollback(db_session)
        # Test function from_database.
        ivp3 = IVP.from_database(db_session, 'Heat2D')
        self.assertEquals(self.ivp.name, ivp3.name)
        method3 = ODEMethod.from_database(db_session, 'radauIIA7')
        self.assertEquals(self.method.name, method3.name)
        template3 = KernelTemplate.from_database(db_session, 'Approx')
        self.assertEquals(self.kernel_templates[0].name, template3.name)
        impl_skeleton3 = ImplSkeleton.from_database(db_session, 'A', self.kernel_templates)
        self.assertEquals(self.impl_skeleton.name, impl_skeleton3.name)
        # Test function close.
        close(db_session)
