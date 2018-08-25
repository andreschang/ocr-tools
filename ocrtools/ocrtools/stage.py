#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#  
#  OCRTOOLS - Making climate data malleable
#  Copyright (C) 2018 Andres Chang
#
###############################################################################

import numpy as np
import os, os.path
import errno
from datetime import datetime
from netCDF4 import Dataset, num2date

class stage(object):

  ndivs = {'daily': 365, 'monthly':12}
  cice_vars = ['aice', 'hi', 'flwdn', 'fswdn']
  cam_vars = ['TS', 'PRECT']
  now = datetime.now()
  scratchId = now.strftime("%Y%m%d%H%M")
  
  def __init__(self, preset = 'andres_local', time_as = 'sequence'):

    """
    stage is used to manage file organization, naming, and output in conjunction with query and build

    A few challenges of working with climate data:
    - Raw data files can be very large (and slow)
    - Filenames are hard to read and can look very similar
    - Hard to remember the details of a specific analysis or summary plot once the Python script
    has been closed (ex. what are the exact coordinates that were used in a spatial average?)

    stage addesses these by allowing you to define a directory structure and naming convention for
    both the source data and output of ocrtools. Most query functions rely on this

    Args:
    * preset (str) [optional]: specifies folder structure for existing data and new output. Please add your 
      own preset and make it the default!
    * time_as (str) [optional]: 'sequence' (ex. base1920.34310, as in 34310 days after 1920-000 [2014]) or
      'date' (ex. 2014-000)
    """

    if preset == 'andres_local':
      self.directories = {"cesm-raw": "/Volumes/Samsung_T5/Open_Climate_Research-Projects/data/raw/cesm", \
      "cesm-reformatted": "/Volumes/Samsung_T5/Open_Climate_Research-Projects/data/reformatted/cesm", \
      "other-reformatted": "/Volumes/Samsung_T5/Open_Climate_Research-Projects/data/reformatted", \
      "scratch": "/Volumes/Samsung_T5/Open_Climate_Research-Projects/data/scratch", \
      "plot": "/Volumes/Samsung_T5/Open_Climate_Research-Projects/data/test_plots", \
      "csv": "/Volumes/Samsung_T5/Open_Climate_Research-Projects/data/streams"}
      self.subfolders = {"cesm-raw": ['dt', 'var_name'], "cesm-reformatted": ["var_name", "yr_range"],
      "other-reformatted": ["var_name", "yr_range"], "scratch": [], "plot": ["function"]}
    elif preset == 'ucar':
      self.directories = {"cesm-raw": "/glade/p/cesmLE/CESM-CAM5-BGC-LE", \
      "cesm-reformatted": "/glade/p/work/andresc/OCR/reformatted/cesm", \
      "other-reformatted": "/glade/p/work/andresc/OCR/reformatted", \
      "scratch": "/glade/p/work/andresc/OCR/scratch"}
      self.subfolders = {"cesm-raw": ['model', 'proc', 'tseries', 'dt', 'src_var_name'], \
      "cesm-reformatted": ["var_name", "yr_range"], "cesm-reformatted": ["var_name", "yr_range"],
      "scratch": []}
      import matplotlib
      matplotlib.use('Agg')
    elif preset == 'demo':
      self.directories = {"cesm-raw": "data/raw", "cesm-reformatted": "data/reformatted", \
      "other-reformatted": "data/reformatted", "scratch": "data/scratch", "plot":"test_plots", \
      "csv": "data/streams"}
    self.naming = {"cesm-reformatted": ["var_name", "dt", "mem", "time_slice"], \
    "other-reformatted": ["var_name", "dt", "time_slice"], "spatial_average": ["var_name", \
    "dt", "time_slice", "mean", "scratchId"], "cesm_report": ["var_name", "dt", "mem"], \
    "other-report": ["var_name", "dt"]}
    self.time_as = time_as
    self.base_yr = 1920

