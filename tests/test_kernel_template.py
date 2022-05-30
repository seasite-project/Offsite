#!/usr/bin/env python3

"""@package test_kernel_template
Tests for function of package kernel_template.
"""

from os import remove
from pathlib import Path
from unittest import TestCase

from offsite.config import Config
from offsite.database import close, open_db
from offsite.descriptions.impl.kernel_template import KernelTemplate
from offsite.descriptions.ode.ivp import IVP
from offsite.descriptions.ode.ode_method import ODEMethod


class TestKernelTemplate(TestCase):
    def setUp(self):
        args = None
        self.config = Config(args)
        self.db_session = open_db('test_kernel_template.db')
        self.ivp = IVP.from_yaml(Path('examples/ivps/Heat2D.ivp'))
        self.kernel_template_lc = KernelTemplate.from_yaml(Path('examples/kernels/pirk/LC.kernel'))
        self.kernel_template_rhs = KernelTemplate.from_yaml(Path('examples/kernels/pirk/RHS.kernel'))
        self.method = ODEMethod.from_yaml(Path('examples/methods/implicit/radauIIA7.ode'))

    def tearDown(self):
        close(self.db_session)
        remove('test_kernel_template.db')

    def test_types(self):
        # Test function generate_pmodel_code.
        for kernel in self.kernel_template_lc.variants:
            kernel.generate_pmodel_code(self.method)
        for kernel in self.kernel_template_rhs.variants:
            kernel.generate_pmodel_code(self.method, self.ivp)
        # Test function to_database.
        self.kernel_template_lc.to_database(self.db_session)
        self.kernel_template_rhs.to_database(self.db_session)
        # Test function from_database.
        KernelTemplate.from_database(self.db_session, 'LC')
