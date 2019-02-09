#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#
#  OCR TOOLS - Making climate data malleable
#  Copyright (C) 2018 Andres Chang
#
###############################################################################

from setuptools import setup, find_packages

setup(name='ocrtools',
      version='0.1.0',
      description='Tools for interpreting and generating new climate data',
      url='https://github.com/andreschang/ocr-tools',
      author='Andres Chang',
      author_email='andresdanielchang@gmail.com',
      license='MIT',
      packages=find_packages(),
      python_requires='>=3',
      install_requires=[
        'numpy', 'geopy', 'matplotlib', 'xarray', 'pandas', 'scipy',
        'tkinter', 'xarray'])
