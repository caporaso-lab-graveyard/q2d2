#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import re
import ast
from setuptools import find_packages, setup

# version parsing from __init__ pulled from Flask's setup.py
# https://github.com/mitsuhiko/flask/blob/master/setup.py
_version_re = re.compile(r'__version__\s+=\s+(.*)')

with open('q2d2/__init__.py', 'rb') as f:
    hit = _version_re.search(f.read().decode('utf-8')).group(1)
    version = str(ast.literal_eval(hit))

classes = """
    Development Status :: 4 - Beta
    License :: OSI Approved :: BSD License
    Topic :: Software Development :: Libraries
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python
    Programming Language :: Python :: 3
    Programming Language :: Python :: 3.3
    Programming Language :: Python :: 3.4
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classes.split('\n') if s]

description = ('Prototype/experiments for microbiome analyses.')

with open('README.md') as f:
    long_description = f.read()

setup(name='q2d2',
      version=version,
      license='BSD',
      description=description,
      long_description=long_description,
      author="scikit-bio development team",
      author_email="gregcaporaso@gmail.com",
      maintainer="scikit-bio development team",
      maintainer_email="gregcaporaso@gmail.com",
      url='http://caporasolab.us',
      test_suite='nose.collector',
      packages=find_packages(),
      scripts=['scripts/q2d2'],
      package_data={'q2d2': ['q2d2/markdown/*md']},
      setup_requires=['numpy >= 1.9.2'],
      install_requires=[
          'scikit-bio >= 0.4.0',
          'IPython >= 3.2.0',
          'ipymd',
          'click',
          'seaborn'
      ],
      extras_require={'test': ["nose", "pep8", "flake8",
                               "python-dateutil"],
                      'doc': ["Sphinx == 1.2.2", "sphinx-bootstrap-theme"]},
      classifiers=classifiers,
      )
