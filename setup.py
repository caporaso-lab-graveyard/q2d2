#!/usr/bin/env python

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
    Development Status :: 1 - Planning
    License :: OSI Approved :: BSD License
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

description = 'Prototype/experiments for microbiome analyses.'

with open('README.md') as f:
    long_description = f.read()

authors = 'https://github.com/gregcaporaso/q2d2/graphs/contributors'

setup(name='q2d2',
      version=version,
      license='BSD',
      description=description,
      long_description=long_description,
      author=authors,
      author_email="gregcaporaso@gmail.com",
      maintainer=authors,
      maintainer_email="gregcaporaso@gmail.com",
      url='https://github.com/gregcaporaso/q2d2',
      packages=find_packages(),
      scripts=['scripts/q2d2'],
      package_data={'q2d2': ['q2d2/markdown/*md']},
      install_requires=[
          'scikit-bio',
          'IPython < 4.0.0',
          'ipymd',
          'click',
          'seaborn'
      ],
      classifiers=classifiers,
      )
