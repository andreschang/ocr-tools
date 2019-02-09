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

setup(name='ocrtools',
      version='0.1.0',
      description='Tools for interpreting and generating new climate data',
      long_description=long_description,
      long_description_content_type='text/markdown',
      url='https://github.com/andreschang/ocr-tools',
      author='Andres Chang',
      author_email='andresdanielchang@gmail.com',
      license='MIT',
      packages=find_packages(),
      python_requires='>=3',
      dependency_links=[
       'ftp://ftp.tcl.tk/pub/tcl/tcl8_6/tk8.6.8-src.tar.gz#egg=tkinter-8.6.8',
       'https://github.com/nbren12/sklearn-xarray#egg=sklearn_xarray-0.0.0'],
      install_requires=[
        'numpy', 'geopy', 'matplotlib', 'xarray', 'pandas', 'scipy',
        'xarray'],
        extras_require={
          'tk_selector': ['tkinter>8.6'],
          'build_regression': ['sklearn_xarray==0.0.0']}
      )
