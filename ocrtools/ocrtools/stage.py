import numpy as np
import os, os.path
import errno
from datetime import datetime
from netCDF4 import Dataset, num2date

## Tools to improve malleability of climate data for 
## research and visualization purposes.
## 

## Global variables

version = "4.3"
ndivs = {'daily': 365, 'monthly':12}
cice_vars = ['aice', 'hi', 'flwdn', 'fswdn']
cam_vars = ['TS', 'PRECT']
now = datetime.now()
scratchId = now.strftime("%Y%m%d%H%M")

class stage(object):

  def __init__(self, preset = 'andres_local', time_as = 'sequence'):
  ## Preset specifies folder schema and location of CESM data
  # preset = 'andres_local'
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

  def add_directory(self, add_directory_name, add_directory_path):
    self.directories[add_directory_name] = add_directory_path

