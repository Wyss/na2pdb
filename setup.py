#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
test with:
    python setup.py install
"""
DESCRIPTION = ("Create PDB files of DNA")
LONG_DESCRIPTION = """
**na2pdb** is a Python package

Create PDB files of DNA based on sequence and apply spatial manipulations to 
them

License is BSD3
RNA is not yet supported
"""

DISTNAME = 'na2pdb'
LICENSE = 'BSD3'
AUTHORS = "Nick Conway"
EMAIL = "nick.conway@wyss.harvard.edu"
URL = ""
DOWNLOAD_URL = ''
CLASSIFIERS = [
    'Development Status :: 1 - Beta',
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3.3',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Cython',
    'Topic :: Scientific/Engineering',
]

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

import os
import sys
import shutil

pjoin = os.path.join
rpath = os.path.relpath

PACKAGE_PATH =      os.path.abspath(os.path.dirname(__file__))
MODULE_PATH =       pjoin(PACKAGE_PATH, 'na2pdb')
DATASETS_PATH =     pjoin(MODULE_PATH, 'data')

# PDB dataset files to include in installation
na2pdb_files = [rpath(pjoin(root, f), MODULE_PATH) for root, _, files in
                 os.walk(DATASETS_PATH) for f in files if '.pdb' in f]

is_py_3 = int(sys.version_info[0] > 2)

setup(
    name=DISTNAME,
    maintainer=AUTHORS,
    packages=['na2pdb'],
    package_data={'na2pdb': na2pdb_files},
    maintainer_email=EMAIL,
    description=DESCRIPTION,
    license=LICENSE,
    url=URL,
    download_url=DOWNLOAD_URL,
    long_description=LONG_DESCRIPTION,
    classifiers=CLASSIFIERS,
    zip_safe=False
)
