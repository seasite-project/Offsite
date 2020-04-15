"""@package db_mapping
Handles mapping of database tables and classes.
"""

from sqlalchemy.orm import mapper

from offsite.descriptions.impl_skeleton import ImplSkeleton, ImplVariant
from offsite.descriptions.ivp import IVP, IVPCharacteristic
from offsite.descriptions.kernel_template import KernelTemplate, Kernel, PModelKernel
from offsite.descriptions.machine import Machine, Compiler
from offsite.descriptions.ode_method import ODEMethod
from offsite.evaluation.benchmark import BenchmarkRecord
from offsite.evaluation.performance_model import KernelRecord, ImplVariantRecord
from offsite.evaluation.ranking import RankingRecord


def mapping():
    """Map database tables and classes.

    Parameters:
    -----------
    -

    Returns:
    --------
    -
    """
    mapper(Machine, Machine.db_table)
    mapper(Compiler, Compiler.db_table)
    mapper(IVP, IVP.db_table)
    mapper(IVPCharacteristic, IVPCharacteristic.db_table)
    mapper(ImplSkeleton, ImplSkeleton.db_table)
    mapper(KernelTemplate, KernelTemplate.db_table)
    mapper(Kernel, Kernel.db_table)
    mapper(PModelKernel, PModelKernel.db_table)
    mapper(ODEMethod, ODEMethod.db_table)
    mapper(ImplVariant, ImplVariant.db_table)
    mapper(BenchmarkRecord, BenchmarkRecord.db_table)
    mapper(KernelRecord, KernelRecord.db_table)
    mapper(ImplVariantRecord, ImplVariantRecord.db_table)
    mapper(RankingRecord, RankingRecord.db_table)
