#!/usr/bin/env python3

from codecs import open  # To use a consistent encoding
from os import path
from re import search, M

from setuptools import setup, find_packages

# Versions should comply with PEP440.  For a discussion on single-sourcing
# the version across setup.py and the project code, see
# https://packaging.python.org/en/latest/single_source_version.html

here = path.abspath(path.dirname(__file__))


def read(*names, **kwargs):
    with open(path.join(path.dirname(__file__), *names),
              encoding=kwargs.get("encoding", "utf8")) as fp:
        return fp.read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                           version_file, M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


# Get the long description from the relevant file
# with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
#    long_description = f.read()

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
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        # 'Intended Audience :: Developers',
        # 'Intended Audience :: Science/Research',
        # 'Topic :: Scientific/Engineering',
        # 'Topic :: Software Development',
        # 'Topic :: Utilities',

        # Pick your license as you wish (should match "license" above)
        # 'License :: OSI Approved :: GNU Affero General Public License v3',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],

    # What does your project relate to?
    # keywords='',

    packages=find_packages(exclude=['tests*']),

    install_requires=[
        'argparse>=1.4.0',
        'attrs>=20.3.0',
        'datetime>=4.3',
        # 'graphviz>=0.14.0', -> offsite_decision_tree
        'kerncraft==0.8.5',
        'lark-parser>=0.11.1',
        'matplotlib>=3.3.3',
        'numpy>=1.19.4',
        'osaca>=0.3.14',
        'pandas>=1.1.5',
        'pathlib>=1.0.1',
        'ruamel.yaml>=0.16.12',
        # 'scikit-learn>=0.23.0', -> offsite_decision_tree
        # 'sklearn', -> offsite_decision_tree
        'sqlalchemy>=1.3.22',
        'sortedcontainers>=2.2.2',
        'sympy>=1.7.1',
    ],
    python_requires='>=3.6',

    entry_points={
        'console_scripts': [
            # Offsite tuning.
            'offsite_bench=offsite.apps.offsite_bench:run',
            'offsite_tune=offsite.apps.offsite_tune:run',
            'offsite_rank=offsite.apps.offsite_rank:run',
            # Code generator.
            'offsite_codegen=offsite.apps.offsite_codegen:run_db',
            'offsite_codegen_from_yaml=offsite.apps.offsite_codegen:run_yaml',
            # Auxiliary apps.
            'offsite_db2name=offsite.apps.auxiliary.offsite_db2name:run',
            'offsite_impl2csv=offsite.apps.auxiliary.offsite_impl2csv:run',
            'offsite_plot=offsite.apps.auxiliary.offsite_plot:run',
        ],
    }
)
