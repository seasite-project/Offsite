"""@package database.db_mapping
Handles mapping of database tables and classes.

@author: Johannes Seiferth
"""

from sqlalchemy.orm import mapper

from offsite.descriptions.impl.impl_skeleton import ImplSkeleton
from offsite.descriptions.impl.kernel_template import KernelTemplate, Kernel, PModelKernel, KernelRecord, PModelRecord
from offsite.descriptions.machine import Compiler, Machine, MachineState, NetworkConfig
from offsite.descriptions.ode import IVP, IVPCharacteristic, ODEMethod
from offsite.ranking.ranking import RankingRecord
from offsite.solver import Solver
from offsite.train.communication.benchmark import BenchmarkRecord
from offsite.train.communication.openmp.omp_barrier import OmpBarrierRecord
from offsite.train.impl_variant import ImplVariant
from offsite.train.impl_variant import ImplVariantRecord


def mapping():
    """Map database tables and classes.

    Parameters:
    -----------
    -

    Returns:
    --------
    -
    """
    mapper(MachineState, MachineState.db_table)
    mapper(Compiler, Compiler.db_table)
    mapper(NetworkConfig, NetworkConfig.db_table)
    mapper(Machine, Machine.db_table)
    mapper(Solver, Solver.db_table)
    mapper(IVP, IVP.db_table)
    mapper(IVPCharacteristic, IVPCharacteristic.db_table)
    mapper(ImplSkeleton, ImplSkeleton.db_table)
    mapper(KernelTemplate, KernelTemplate.db_table)
    mapper(Kernel, Kernel.db_table)
    mapper(PModelKernel, PModelKernel.db_table)
    mapper(ODEMethod, ODEMethod.db_table)
    mapper(ImplVariant, ImplVariant.db_table)
    mapper(BenchmarkRecord, BenchmarkRecord.db_table)
    mapper(OmpBarrierRecord, OmpBarrierRecord.db_table)
    mapper(KernelRecord, KernelRecord.db_table)
    mapper(PModelRecord, PModelRecord.db_table)
    mapper(ImplVariantRecord, ImplVariantRecord.db_table)
    mapper(RankingRecord, RankingRecord.db_table)
