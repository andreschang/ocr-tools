#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#
#  OCR TOOLS - Making climate data malleable
#  Copyright (C) 2018 Andres Chang
#
###############################################################################

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

data = ['ocrtools/var_lists/cam_vars.csv',
             'ocrtools/var_lists/cice_vars.csv',
             'ocrtools/var_lists/clm_vars.csv',
             'ocrtools/var_lists/pop_vars.csv',
             'ocrtools/images/lambert_cylindrical.gif']

setup(name='ocrtools',
      version='0.1.0',
      description='Tools for interpreting and generating new climate data',
      long_description=long_description,
      long_description_content_type='text/markdown',
      url='https://github.com/andreschang/ocr-tools',
      author='Andres Chang',
      author_email='andresdanielchang@gmail.com',
      license='MIT',
      packages=['ocrtools'],
      package_data={'ocrtools': data},
      python_requires='>=3.0',
      install_requires=['numpy', 'matplotlib', 'pandas', 'xarray', 'geopy'],
      setup_requires=["pytest-runner"],
      tests_require=["pytest"])
