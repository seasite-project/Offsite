#!/usr/bin/env python3

from codecs import open
from os import path
from re import search, M

from setuptools import setup, find_packages

# Versions should comply with PEP440.  For a discussion on single-sourcing
# the version across setup.py and the project code, see
# https://packaging.python.org/en/latest/single_source_version.html

here = path.abspath(path.dirname(__file__))


def read(*names, **kwargs):
    with open(path.join(path.dirname(__file__), *names), encoding=kwargs.get("encoding", "utf8")) as fp:
        return fp.read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_file, M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


setup(
    name='offsite',
    version=find_version('offsite', '__init__.py'),

    description='',
    # long_description=long_description,

    # url='',

    # Author details
    author='Johannes Seiferth',
    author_email='johannes.seiferth@uni-bayreuth.de',

    # Choose your license
    # license='AGPLv3',

    # https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 4 - Beta',

        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Topic :: Software Development',
        'Topic :: Utilities',

        'License :: OSI Approved :: GNU Affero General Public License v3',

        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],

    packages=find_packages(exclude=['tests*']),

    install_requires=[
        'argparse>=1.4.0',
        'attrs>=21.4.0',
        'kerncraft==0.8.14',
        'lark-parser>=0.12.0',
        'matplotlib>=3.5.2',
        'numpy>=1.21.6',
        'osaca==0.4.8',
        'pandas>=1.3.5',
        'pathlib2>=2.3.7',
        'pycparser==2.21',
        'ruamel.yaml>=0.17.21',
        'sqlalchemy>=1.4.36',
        'sortedcontainers>=2.4.0',
        'sympy>=1.9',
    ],
    python_requires='>=3.7,<3.10',

    entry_points={
        'console_scripts': [
            # Offsite tuning.
            'offsite_bench_openmp=offsite.apps.offsite_bench_openmp:run',
            'offsite_tune_node=offsite.apps.offsite_tune_node:run',
            'offsite_rank=offsite.apps.offsite_rank:run',
            # Code generator.
            'offsite_codegen=offsite.apps.offsite_codegen_db:run',
            'offsite_codegen_from_yaml=offsite.apps.offsite_codegen_yaml:run',
            # Online tuning.
            'offsite_gensolver=offsite.apps.offsite_gensolver:run',
            # Auxiliary apps.
            'offsite_db2name=offsite.apps.auxiliary.offsite_db2name:run',
            'offsite_impl2csv=offsite.apps.auxiliary.offsite_impl2csv:run',
            'offsite_plot=offsite.apps.auxiliary.offsite_plot:run',
        ],
    }
)
