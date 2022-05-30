"""@package codegen.generator.impl
Impl variant code generator classes.
"""
from enum import Enum
from typing import Optional

from pathlib2 import Path
from sqlalchemy.orm import Session

from offsite.codegen.generator.impl.impl_generator_c import ImplCodeGeneratorC
from offsite.codegen.generator.impl.impl_generator_cpp import ImplCodeGeneratorCPP
from offsite.codegen.generator.impl.impl_generator_cpp_mpi import ImplCodeGeneratorCppMPI


class GeneratedCodeLanguageType(Enum):
    """Defines what type of code is generated.

    - C
        C code.
    - CPP
        C++ code.
    - CPP_MPI
        C++ MPI code.
    """
    C = 'C'
    CPP = 'CPP'
    CPP_MPI = 'MPI'


def make_impl_code_generator(type_str: GeneratedCodeLanguageType, db_session: Session, folder_impl: Path,
                             folder_ivp: Path, folder_method: Path, folder_ds: Optional[Path] = None):
    if type_str == GeneratedCodeLanguageType.C:
        return ImplCodeGeneratorC(db_session, folder_impl, folder_ivp, folder_method)
    elif type_str == GeneratedCodeLanguageType.CPP:
        return ImplCodeGeneratorCPP(db_session, folder_impl, folder_ivp, folder_method, folder_ds)
    elif type_str == GeneratedCodeLanguageType.CPP_MPI:
        return ImplCodeGeneratorCppMPI(db_session, folder_impl, folder_ivp, folder_method, folder_ds)
    else:
        raise RuntimeError('TODO')  # TODO
