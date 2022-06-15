#!/usr/bin/env python3

from codecs import open
from os import path
from re import search, M

from setuptools import setup, find_packages

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
    name='offsite_aux',
    version=find_version('offsite_aux', '__init__.py'),

    description='',

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
        'Programming Language :: Python :: 3.10',
    ],

    packages=find_packages(exclude=['tests*']),

    install_requires=[
        # 'graphviz>=0.14.0', -> offsite_decision_tree
        'matplotlib>=3.5.2',
        'offsite>=0.3.3',
        # 'scikit-learn>=0.23.0', -> offsite_decision_tree
        # 'sklearn', -> offsite_decision_tree
    ],
    python_requires='>=3.7',

    entry_points={
        'console_scripts': [
            # Auxiliary apps.
            'offsite_aux_db2name=offsite.apps.auxiliary.offsite_db2name:run',
            'offsite_aux_impl2csv=offsite.apps.auxiliary.offsite_impl2csv:run',
            'offsite_aux_kernel2csv=offsite.apps.auxiliary.offsite_kernel2csv:run',
            'offsite_aux_plot=offsite.apps.auxiliary.offsite_plot:run',
        ],
    }
)
