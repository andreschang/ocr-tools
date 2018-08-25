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

ndivs = {'daily': 365, 'monthly':12}
cice_vars = ['aice', 'hi', 'flwdn', 'fswdn']
cam_vars = ['TS', 'PRECT']
now = datetime.now()
scratchId = now.strftime("%Y%m%d%H%M")

class query(object):
  def __init__(self, stage = None, verbose = True, **kwargs):
    """
    Initializes a query. This is the core analysis feature of ocrtools.
    Should always be followed by set_params function
    
    Args:
    * stage (stage class object) [optional]: see stage.py for documentation
    * verbose (bool)

    Kwargs:
    * file (str): path to data addressed by query. If not specified, this will be
    assigned by the user or automatically when the set_params function is executed
    * src (str): helps ocrtools autofill parameters, if the source is known (ex. CESM)
    """

    self.verbose = verbose
    try:
      self.file = [kwargs["file"]]
    except:
      pass
    try:
      self.src = kwargs["src"]
    except:
      pass
    try:
      self.directories = stage.directories
    except:
      self.directories = {"cesm-raw": "", "cesm-reformatted" : "", "other-reformatted": "", \
      "scratch": "", "plot": "", "csv": ""}
    try:
      self.subfolders = stage.subfolders
    except:
      self.subfolders = {"cesm-raw": "", "cesm-reformatted" : "", "other-reformatted": "", \
      "scratch": "", "plot": "", "csv": ""}
    try:
      self.naming = stage.naming
    except:
      self.naming = {"cesm-reformatted": ["var_name", "dt", "mem", "time_slice"], \
    "other-reformatted": ["var_name", "dt", "time_slice"], "spatial-average": ["var_name", \
    "dt", "yr_range", "mean", "scratchId"], "cesm-report": ["var_name", "dt", "mem"], \
    "other-report": ["var_name", "dt"]}
    try:
      self.time_as = stage.time_as
    except:
      self.time_as = "date"
    try:
      self.base_yr = stage.base_yr
    except:
      self.base_yr = 1920

  ## Primary functions  ##
  ########################

  def reduce_data(self, **kwargs):
    """
    Slices large datasets into manageable pieces that can later be accessed
    to save time and memory
    
    Kwargs:
    * yr0 (int) [optional]: First year that will be used to make chunked data
      files
    * yrf (int) [optional]: Last year to be used for chunked data files
    * custom_tags (str) [optional]: String added to output filenames
    * print_report (bool) [optional]: If True, adds a report with conversion
      info to the reduced data directory
    * add_to_report (str) [optional]: Text to be added to the report
    """
    try:
      self.yr0 = kwargs["yr0"]
    except:
      pass
    try:
      self.yrf = kwargs["yrf"]
    except:
      pass

    custom_tags, print_report, add_to_report = self.assign_reporting_vars(kwargs)
    ts_var, offset = self.timesync()

    print(ts_var.shape)
    report, report_tag = self.reduce_report()

    for i in range(int(self.nt/120)):
      div_in_var = i*120
      data = ts_var[div_in_var: div_in_var+120]
      div_since_base_yr = div_in_var+offset + (self.data_yr0[0]-self.base_yr)*self.ndiv
      div_since_zero = div_in_var +offset+ self.data_yr0[0]*self.ndiv

      fYr0, fDiv0 = divmod(div_since_zero, self.ndiv)
      fYrF, fDivF = divmod(div_since_zero+120, self.ndiv)

      div_format = str(int(div_since_base_yr))

      f_out = self.outfile(mode = 'reduce', div_since_base_yr = div_format, \
        fYr0 = fYr0, fDiv0 = fDiv0, fYrF = fYrF, fDivF = fDivF, custom_tags = custom_tags)

      dataset = Dataset(f_out, 'w', format='NETCDF4_CLASSIC')
      fdims = {}
      for index, dimension in enumerate(self.dim):
        print(dimension)
        fdims[dimension] = dataset.createDimension(dimension, data.shape[index])

      fvars = {}

      if len(self.lat0.shape) > 1:
        if self.axes[self.lat_name] < self.axes[self.lon_name]:
          grid_axes = ('lat', 'lon')
        else:
          grid_axes = ('lon', 'lat')

      for index, dimension in enumerate(self.dim):
        if (((dimension == self.lat_name) or (dimension == self.lon_name)) and (len(self.lat0.shape) > 1)):
          fvars[dimension] = dataset.createVariable(dimension, np.float32, grid_axes)
        else:
          fvars[dimension] = dataset.createVariable(dimension, np.float32, (dimension,))

      print(data.shape)
      print(fvars)

      if self.dt == 'daily':
        fYr0, fDiv0 = '{:04d}'.format(fYr0), '{:03d}'.format(fDiv0)
        fYrF, fDivF = '{:04d}'.format(fYrF), '{:03d}'.format(fDivF)
        print('\n formatting ... yr ' + str(fYr0) + ' day ' + str(fDiv0) + 
          '\n to yr ' + str(fYrF) + ' day ' + str(fDivF))
        fvars['time'].units = 'days since '+str(self.base_yr)+'-00-00'

      elif self.dt == 'monthly':
        fYr0, fDiv0 = '{:04d}'.format(fYr0), '{:02d}'.format(fDiv0)
        fYrF, fDivF = '{:04d}'.format(fYrF), '{:02d}'.format(fDivF)
        print('\n formatting ... yr ' + str(fYr0) + ' month ' + str(fDiv0) + 
          '\n to yr ' + str(fYrF) + ' month ' + str(fDivF))
        fvars['time'].units = 'months since '+str(self.base_yr)+'-00'

      print(div_since_base_yr)
      print(div_since_base_yr+120)
      print(np.arange(div_since_base_yr,div_since_base_yr+120))
      fvars['time'][:] = np.arange(div_since_base_yr,div_since_base_yr+120)
      fvars['lat'].units = self.units['lat']
      fvars['lon'].units = self.units['lon']
      fvars['lat'][:] = self.lat0
      fvars['lon'][:] = self.lon0

      if 'level' in self.dim:
        fvars['level'][:] = self.level0
        fvars['level'].units = self.units['level']
      if len(self.lat0.shape) > 1:
        fvars['cellarea'] = dataset.createVariable('cellarea', np.float32, grid_axes)
        fvars['cellarea'][:] = self.cellarea0
        fvars['cellarea'].units = self.units['cellarea']

      saveVar = dataset.createVariable(self.var_name,  np.float64, tuple(self.dim))
      saveVar.units = self.units['var']
      dataset.description = report
      saveVar[:] = data
      dataset.close()

    ## 3. PRINT REPORT
    ##
    if print_report == True:
      report += '\n\n'+add_to_report
      report_out = self.outfile(mode = 'report', custom_tags = report_tag)

      np.savetxt(report_out, [report], fmt = '%s')

  def spatial_average(self, lat_bounds = [-89., 89.], lon_bounds = [-179., 179.], lon_convert = True, output = 'plot', **kwargs):
    """
    Calculates the weighted average of a 3-dimensional variable (2D space plus time)
    and outputs a graph or list of the timeseries. Must be run after set_params

    Args:
    * lat_bounds (list) [optional]: 2-item list of lower and upper latitude bounds
    * lon_bounds (list) [optional]: 2-item list of longitude bounds (use '-' for longitudes
      west of the prime meridian and positive for east [-180 to 180] OR set values in the range
      [0-360] and set lon_convert to FALSE)
    * lon_convert (bool) [optional]: most climate datasets (perhaps not all) use longitude
      between 0-360, instead of the more intuitive -180 to 180. This will convert the longitudinal
      range by default
    * output (str) [optional]: output timeseries as a 'plot', 'csv', or 'list'. Plot and csv are
      saved as specified by the stage. List does not save external to python script.

    Kwargs:
    * custom_tags (str) [optional]: String added to output filenames
    * print_report (bool) [optional]: If True, get_latlon_indices() will output a map
      off data points corresponding to the spatial average
    """

    custom_tags, print_report, add_to_report = self.assign_reporting_vars(kwargs)
    ts_var, offset = self.timesync()
    ts_var = ts_var[:self.nt]
    print(ts_var.shape)
    self.function = 'spatial_average'
    mvar = []

    ## 2. CRUNCH IT UP!
    ## Only works with 2d variables currently
    ##
    if output == 'list':
      lat_indices, lon_indices = self.get_latlon_indices(lat_bounds, lon_bounds, print_report = False)
    else:
      lat_indices, lon_indices = self.get_latlon_indices(lat_bounds, lon_bounds, print_report = print_report, \
      directory = output)

    if len(self.lat0.shape) == 1:
      bounded_var = self.bound_var(ts_var, lat_indices, lon_indices, lat_axis = self.axes['lat'], lon_axis = self.axes['lon'])
      mm_var = np.mean(bounded_var, axis = self.axes['lon'])

      lat_range = [self.lat0[k] for k in lat_indices]
      wgt = self.reg_wgt(lat_range[0], lat_range[-1], len(lat_range))
      mvar0 = 0
      for k in range(len(lat_range)):
        mvar0 += mm_var[:,k]*wgt[k]

    elif len(self.lat0.shape) == 2:
      total_area = 0.
      sum_var = 0.
      for i in range(len(lat_indices)):
        if (np.isnan(ts_var[0, lat_indices[i], lon_indices[i]]) == False and \
          ts_var[0, lat_indices[i], lon_indices[i]]) < 3e10:
          sum_var += ts_var[:, lat_indices[i], lon_indices[i]]*self.cellarea0[lat_indices[i], lon_indices[i]]
          total_area += self.cellarea0[lat_indices[i], lon_indices[i]]
      # print(sum_var)
      # print(total_area)
      mvar0 = sum_var/total_area
      # print(mvar0)

    roundEnd = (self.yrf-self.yr0+1)*self.ndiv
    mvar += list(mvar0)[:roundEnd]
    xtime = np.arange(self.yr0, self.yrf+1, 1./self.ndiv)
    print('Spatial average created with length (time), '+str(len(mvar)))

    if output == 'plot':
      import matplotlib.pyplot as plt
      print(xtime.shape, len(mvar))
      plt.plot(xtime, mvar)
      f1 = self.outfile(mode = 'spatial_average-plot', custom_tags = custom_tags)
      plt.savefig(f1, dpi = 200)
      if print_report == True:
        self.name_sync(f1, ['mean'], ['bounded_area'])

    elif output == 'csv':
      f1 = self.outfile(mode = 'spatial_average-csv', custom_tags = custom_tags)
      print(f1)
      np.savetxt(f1, mvar, delimiter = ',')
      if self.verbose == True:
        name_sync(f1, ['mean','csv'], ['bounded_area','png'])

    elif output == 'list':
      return mvar

  def set_params(self, **kwargs):
    """
    Sets the parameters that enable ocrtools to read the data file.
    This function is highly flexible.

    A key feature of ocrtools is that is reads metadata and autofills as
    many parameters as possible. If src is unknown, ocrtools will prompt the user
    to set parameters interactively (and offer lots of help). 

    If src is known and supported (only 'cesm' in v0.1), many parameters (including 
    file path) are set automatically; however, the user still must specify some info
    about the analysis (ex. yr0 and yrf), so that ocrtools knows which file(s) to use

    Parameters do NOT have to be set interactively; they can also be included
    as kwargs. See below

    Kwargs:
    * src (str) [optional]: source of raw data file. Overrides any src set in query
      init. If src has not been set at all (or src is not recognized), src is 
      treated as 'unknown'

    Kwargs (src = 'unknown'):
    * file (str) [optional]: overrides any filepath set in query init. If both left
      blank, user must input interactively
    * var (str) [optional]: main variable examined by query. If left blank,
      ocrtools will run first_look() and prompt user to input interactively
    * lat_name (str) [optional]: name of lat dimension. If left blank, ocrtools will
      autofill or prompt user to input interactively
    * lon_name (str) [optional]: name of lon dimension. If left blank, ocrtools will
      autofill or prompt user to input interactively
    * time_name (str) [optional]: name of time dimension. If left blank, ocrtools will
      autofill or prompt user to input interactively
    * data_yr0 (int) [optional]: first year of data. ALL DATA ASSUMED TO START AT 
      BEGINNING OF CALENDAR YEAR. If left blank, ocrtools will provide help so the 
      user can input interactively
    * dt (str) [optional]: "monthly" or "daily". If left blank, ocrtools will provide
      help so the user can input interactively
    * dim (list) [optional]: list of all main variable dimensions (strings) in order. 
      If left blank, ocrtools will provide help so the user can input interactively
    * cellarea_name (str) [optional]: only needed for non-rectilinear grids (i.e. grid 
      where each lat/lon point has a 2d reference pair). If left blank and grid is 
      non-rectilinear, will autofill or prompt user to input interactively
    * level_name (str) [optional]: only needed for 4D data (ex. space, time, and depth
      like atmospheric or ocean level). If left blank and data is 4D, will autofill or
      prompt user to input interactively
    * yr0 (int) [optional]: calendar year 0 examined by query. Can be set later with
      set_yr0() or if never specified, query will set yr0 = data_yr0
    * yrf (int) [optional]: calendar year F (inclusive) examined by query. Can be set 
      later with set_yrf or if never specified, query will set yrf as the final year 
      of data

    Kwargs (src = 'cesm'):
    * var (str) [optional]: name of main variable examined by query. Does not 
      need to include timestep tag (ex. daily TS should be var = 'TS', 
      not var = 'TS_d'). If left blank, user must input interactively
    * dt (str) [optional]: "monthly" or "daily". If left blank, user must 
      input interactively
    * mem (int) [optional]: CESM member 2-34 or 1850 (control run). If left blank,
      user must input interactively
    * hemisphere (str) [optional]: only for cice model variables "nh" or "sh"
    * yr0 (int) [optional]: calendar year 0 examined by query. If left blank,
      user must input interactively
    * yrf (int) [optional]: calendar year F (inclusive) examined by query. If left
      blank, user must input interactively
    * ndim (int) [optional]: number of dimensions in main variable (3 or 4). If left
      blank, user must input interactively
    """

    try:
      src = kwargs["src"]
    except KeyError:
      try:
        src = self.src
      except AttributeError:
        src = 'unknown'
    self.src = src

    if src == 'cesm':
      try:
        var_name0 = kwargs["var"]
      except KeyError:
        var_name0 = input("Please enter variable name: ")
      if var_name0[-2:len(var_name0)] == "_d":
        self.var_name = var_name0[:(len(var_name0)-2)]
      else:
        self.var_name = var_name0

      try:
        self.dt = kwargs["dt"]
      except KeyError:
        self.dt = input("Please enter dt (monthly or daily): ")

      try:
        self.mem = kwargs["mem"]
      except KeyError:
        self.mem = int(input("Please enter CESM member (2-34 or 1850): "))

      if self.var_name in cice_vars:
        try:
          self.hemisphere = kwargs["hemisphere"]
        except KeyError:
          self.hemisphere = input("Please enter CESM hemisphere (nh or sh): ")

      try:
        self.yr0 = int(kwargs["yr0"])
      except KeyError:
        self.yr0 = int(input("Please enter yr0 of analysis: "))

      try:
        self.yrf = int(kwargs["yrf"])
      except KeyError:
        self.yrf = int(input("Please enter yrf of analysis (inclusive): "))

      try:
        self.dim = kwargs["dim"]
      except KeyError:
        try:
          self.ndim = kwargs["ndim"] 
        except KeyError:
          self.ndim = int(input("Please enter the number of dimensions in your CESM variable (3 or 4): "))

      self.fname = self.fill_cesm_params()
      self.f_open = None

    else:
      print("")
      try:
        self.file = [kwargs["file"]]
      except KeyError:
        try:
          self.file
        except AttributeError:
          self.file = [input('Enter file: ')]

      try:
        self.var_name = kwargs["var"]
        self.f_open = Dataset(self.file[0])
        self.fvars = self.f_open.variables
      except KeyError:        
        self.first_look()
        self.var_name = input('Enter main variable: ')

      self.src_var_name = self.var_name
      fdim = list(self.fvars[self.src_var_name].dimensions)

      ## SET LON AND LAT NAMES
      ##
      try:
        self.lat_name, self.lon_name = kwargs["lat_name"], kwargs["lon_name"]
      except KeyError:
        manual_entry = False
        if(('lon' in self.fvars) and ('lat' in self.fvars)):
          self.lon_name, self.lat_name = 'lon', 'lat'
        elif(('LON' in self.fvars) and ('LAT' in self.fvars)):
          self.lon_name, self.lat_name = 'LON', 'LAT'
        elif(('TLON' in self.fvars) and ('TLON' in fdim) and \
          ('TLAT' in self.fvars) and ('TLAT' in fdim)):
          self.lon_name, self.lat_name = 'TLON', 'TLAT'
        elif(('ULON' in self.fvars) and ('ULON' in fdim) and\
          ('ULAT' in self.fvars) and ('ULAT' in fdim)):
          self.lon_name, self.lat_name = 'ULON', 'ULAT'
        else:
          print("\nCould not set lat and lon names automatically.\nShowing all variable metadata ... \n")
          print(self.fvars[self.var_name])
          self.lat_name = input('\nEnter lat_name: ')
          self.lon_name = input('Enter lon_name: ')
          manual_entry = True

        if(manual_entry == False):
          print("lat_name set to "+self.lat_name)
          print("lon_name set to "+self.lon_name)

      ## SET CELLAREA NAME
      ##
      try:
        self.cellarea_name = kwargs["cellarea_name"]
      except KeyError:
        if(len(self.fvars[self.lat_name].shape)> 1):
          manual_entry = False
          if(self.lat_name == 'TLAT'):
            if('tarea' in self.fvars):
              self.cellarea_name = 'tarea'
            elif('TAREA' in self.fvars):
              self.cellarea_name = 'TAREA'
          elif(self.lat_name == 'ULAT'):
            if('uarea' in self.fvars):
              self.cellarea_name = 'uarea'
            elif('UAREA' in self.fvars):
              self.cellarea_name = 'UAREA'
          else:          
            self.cellarea_name = input("\nIt looks like this data is spread across a curvilinear grid\n"+
            "If there is a cellarea variable, please enter it here (or NONE): ")
            manual_entry = True

          if(manual_entry == False):
            print("cellarea_name set to "+self.cellarea_name)
        else:
          self.cellarea_name = None

      ## SET TIME NAME
      ##
      try:
        self.time_name = kwargs["time_name"] 
      except KeyError:
        if('time' in self.fvars):
          self.time_name = 'time'
        elif('TIME' in self.fvars):
          self.time_name = 'TIME'
        else:
          print("\nCould not set time name automatically.\nShowing all variable metadata ... \n")
          print(self.fvars[self.var_name])
          self.time_name = input('\nEnter time_name: ')

      ## SET LEVEL NAME
      ##
      try:
        self.level_name = kwargs["level_name"]
      except KeyError:
        if(len(fdim)> 3):
          for j in range(len(fdim)):
            if((fdim[j] != self.time_name) and (fdim[j] != self.lat_name) and 
              (fdim[j] != self.lon_name)):
              if(fdim[j] in self.fvars):
                self.level_name = fdim[j]
          print("level_name set to "+self.level_name)
        else:
          self.level_name = None

      ## SET DATA YR0
      ##
      try:
        self.data_yr0 = [kwargs["data_yr0"]]
      except KeyError:
        data_yr0 = input('Enter data yr0 (or HELP): ')
        if(data_yr0 == "HELP"):
          time_attrs = self.fvars[self.time_name].ncattrs()
          if('units' in time_attrs):
            units = self.fvars[self.time_name].units
            if('calendar' in time_attrs):
              calendar = self.fvars[self.time_name].calendar
            else:
              calendar = 'standard'

            times = num2date(self.fvars[self.time_name][:], units = units, calendar = calendar)
            print("Best guess for time[0]: "+times[0].strftime("%B %d, %Y"))
            show_time_meta = input("Enter SHOW to show more time metadata: ")
            if show_time_meta == "SHOW":
              print("\nShowing time metadata ...\n")
              print(self.fvars[self.time_name])
              print("\nShowing first values in time variable...\n")
              print(self.fvars[self.time_name][0:4])
          else:
            print("\nShowing time metadata ...\n")
            print(self.fvars[self.time_name])
            print("Showing first values in time variable...\n")
            print(self.fvars[self.time_name][0:4])

          data_yr0 = input('Enter data yr0: ')

        self.data_yr0 = [int(data_yr0)]
        print(self.data_yr0)

      ## SET DT
      ##
      try:
        self.dt = kwargs["dt"]
      except KeyError:
        dt = input('Enter dt (monthly, daily, or HELP): ')
        if dt == "HELP":
          if('units' in self.fvars[self.time_name].ncattrs()):
            units = self.fvars[self.time_name].units
            step = self.fvars[self.time_name][:][1]- self.fvars[self.time_name][:][0]
            if('days since' in units):
              step_unit = 'days'
            elif('hours since'):
              step_unit = 'hours'
            else:
              step_unit = 'units'
            print(str(step)+" "+step_unit+" in first time step")
          else:
            print("\nShowing time metadata ...\n")
            print(self.fvars[self.time_name])
            print("\nShowing first values in time variable...\n")
            print(self.fvars[self.time_name][0:4])
          dt = input('Enter dt (monthly or daily): ')
        self.dt = dt

      ## SET DIM ORDER FOR VARIABLE IN RAW DATA
      ##
      try:
        self.dim = kwargs["dim"]
      except KeyError:
        dim = input('Enter list of '+ self.var_name+' dimensions ordered by axis (or HELP): ')
        if dim == "HELP":
          print("\nHere, we are specifying which axes of "+self.var_name+", correspond to which dimensions."+
            "\n\nShowing variable dims and shape from metadata")
          print("["+", ".join([str(fdim[k]) + " : " + str(self.fvars[self.var_name][:].shape[k]) for k in \
            range(len(self.fvars[self.var_name][:].shape))])+"]")
          alldims = [self.time_name, self.lat_name, self.lon_name]
          if self.level_name != None:
            alldims.append(self.level_name)
          dim = input('\nEnter list of '+ self.var_name+' dimensions ordered by axis (like this: '+
            '["'+'", "'.join(alldims)+'"], but maybe in a different order): ')
        dim = dim.replace("[", "").replace("]", "").replace("'", "").replace('"', '')
        dim = dim.split(",")
        dim = [dim0.strip() for dim0 in dim]
        self.dim = dim

      if "yr0" in kwargs:
        self.yr0 = int(kwargs["yr0"])
      if "yrf" in kwargs:
        self.yrf = int(kwargs["yrf"])

  def first_look(self, **kwargs):
    """
    Opens a netcdf file and prints a list of variables
    or metadata of a specified variable

    Kwargs:
    * file (string) [optional]: path to netcdf data file. If left blank,
      user must input interactively
    * var (string) [optional]: name of main variable to be examined. If
      left blank, ocrtools will just print all variable keys (skipping 
      metadata of the main variable)
    """

    try:
      self.f_open
    except:
      try:
        self.file
      except:
        try:
          self.file = [kwargs["file"]]
        except:
          self.file = [input('Enter file: ')]

      self.f_open = Dataset(self.file[0])

    self.fvars = self.f_open.variables

    try:
      var_look = kwargs["var"]
      try:
        print("Showing metadata for variable, "+kwargs["var"]+"\n")
        print(self.f_open.variables[var_look])
      except:
        print("Variable not in file")
    except:
      print("\nShowing variables ...\n")
      print(', '.join(self.f_open.variables.keys())+'\n')

  ## Core data-wrangling utilities ##
  ###################################

  def extract_data(self):
    """
    Binds data to query class attributes once
    params have been set
    """
    
    ## 1. LOAD DATA
    start_index = 0
    if self.f_open == None:
      try:
        self.file, start_index, end_index = self.grab_reduced_data()
        self.data_yr0 = [self.yr0]
      except:
        pass
      self.f_open = Dataset(self.file[0])

    self.axes = {}
    self.axes['time']  = self.dim.index(self.time_name)
    self.var = self.f_open.variables[self.src_var_name][start_index:]
    print(self.var.shape)

    for ff in np.arange(1, len(self.file)):
      print('Concatenating netcdf '+str(ff))
      add_f_open = Dataset(self.file[ff])
      self.var = np.concatenate((self.var, add_f_open.variables[self.src_var_name][:]), axis = self.axes['time'])
      add_f_open.close()

    print(self.var.shape)
    self.setup_axes()
    self.units = {}

    self.lat0 = self.f_open.variables[self.lat_name][:]
    self.lon0 = self.f_open.variables[self.lon_name][:]
    self.units['var'] = self.f_open.variables[self.src_var_name].units
    self.units['lat'] = self.f_open.variables[self.lat_name].units
    self.units['lon'] = self.f_open.variables[self.lon_name].units

    if len(self.lat0.shape) > 1:
      self.cellarea0 = self.f_open.variables[self.cellarea_name][:]
      self.units['cellarea'] = self.f_open.variables[self.cellarea_name].units
    else:
      self.cellarea0 = False 
    if 'level' in self.dim:
      self.level0 = self.f_open.variables[self.level_name][:]
      self.units['level'] = self.f_open.variables[self.level_name].units
    else:
      self.level0 = False

    self.f_open.close()
    self.f_open = None

  def grab_reduced_data(self):
    """
    Checks to see if climate data has been chunked by reduce_data
    for the query time period under study. If so, ocrtools will use
    these data files to save time and memory
    """

    if self.src == 'cesm':
      folder_type = 'cesm-reformatted'
    else:
      folder_type = 'other-reformatted'

    if "yr_range" in self.subfolders[folder_type]:
      yr_range = str(self.yr0)+'-'+str(self.yrf)
      ideal_folder = self.folder_path(folder_type)
      end_index = ideal_folder.find(yr_range)
      f0 = ideal_folder[0:end_index]
      print("Checking for reduced data in "+f0)

      subdirs = self.get_subs(f0)
      yr_ranges = [[int(x.split('-')[0]), int(x.split('-')[1])] for x in subdirs]

      for span in yr_ranges:
        if ((span[0] <= self.yr0) and (span[1]-1 >= self.yrf)):
          target_folder = f0+('-').join(map(str, span))
          print('Fetching reduced data from '+target_folder)
          rf_yr0 = span[0]
          nstartfile, start_index = divmod((self.yr0-rf_yr0)*ndivs[self.dt], 120)
          endfile, end_index = divmod((self.yrf-rf_yr0+1)*ndivs[self.dt], 120)
          all_ncs = self.get_ncs(target_folder)
          target_ncs = [target_folder+'/'+file for file in all_ncs[nstartfile:endfile+1]]
          self.simple_params()

          return target_ncs, start_index, end_index
        else:
          return False
      return False

  def get_latlon_indices(self, lat_bounds = [-89., 89.], lon_bounds = [-179., 179.], \
    lon_bound_convert = True, data_wrap_lon = 0., print_report = True, directory = 'scratch'):
    """
    Returns a tuple of indice lists for area bounded by lat and lon bounds.
    If print_report is True, the function will also produce a map
    of the data points used for reference

    Args:
    * lat_bounds (list) [optional]: 2-item list of lower and upper latitude bounds
    * lon_bounds (list) [optional]: 2-item list of longitude bounds (use '-' for longitudes
      west of the prime meridian and positive for east [-180 to 180] OR set values in the range
      [0-360] and set lon_convert to FALSE)
    * lon_convert (bool) [optional]: most climate datasets (perhaps not all) use longitude
      between 0-360, instead of the more intuitive -180 to 180. This will convert the longitudinal
      range by default
    * data_wrap_lon (float) [optional]: The longitude value where data "wraps" around the globe.
      For 0-360 datasets (0 at prime meridian), this would be 0 or 360 (equivalent)
    * print_report (bool) [optional]: if True, function draws a map of data points and saves it
      in the specified directory
    * directory (str) [optional]: Stage directory where plot is saved

    """

    resubmit = False
    nsubmit = 1
    mlon_bounds = lon_bounds
    if lon_bound_convert == True:
      lon_bounds = [(bound+360. if bound < 0 else bound) for bound in lon_bounds]
      if self.verbose == True:
        print("Converting lon bounds from 180W-180E to 0-360")
        print("New lon bounds: "+",".join(map(str, lon_bounds)))

    if lon_bounds[0] > lon_bounds[1]:
      if ((data_wrap_lon == 0.) or (data_wrap_lon == 360.)):
        wrap_points = [360., 0.]
      elif ((data_wrap_lon == -180.) or (data_wrap_lon == 180.)):
        wrap_points = [180., -180.]
      resubmit = True
      nsubmit = 2
      lon_bounds2 = [wrap_points[1], lon_bounds[1]]
      lon_bounds = [lon_bounds[0], wrap_points[0]]

    if print_report == True:
      from mpl_toolkits.basemap import Basemap
      x, y = [], []
      llcrnrlon = mlon_bounds[0]-15 if mlon_bounds[0] > -165 else -180
      urcrnrlon = mlon_bounds[1]+15 if mlon_bounds[1] < 165 else 180
      llcrnrlat = lat_bounds[0]-10 if lat_bounds[0] > -80 else -90
      urcrnrlat = lat_bounds[1]+10 if lat_bounds[1] < 80 else 90
      m = Basemap(projection='cyl',llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat,\
          llcrnrlon=llcrnrlon,urcrnrlon=urcrnrlon,resolution='c')
      m.drawparallels(np.arange(-75,75,25),labels=[1,0,0,0])
      m.drawmeridians(np.arange(0,360,40), labels = [0,0,1,0])
      m.fillcontinents(zorder = 1, color = 'wheat', alpha = 0.7)

    while nsubmit > 0:
      if len(self.lat0.shape) == 1:
        lat_minmax = [self.find_nearest(self.lat0, bound) for bound in lat_bounds]
        lon_minmax = [self.find_nearest(self.lon0, bound) for bound in lon_bounds]

        if lon_minmax[1][0] < lon_minmax[0][0]:
          lon_minmax = list(reversed(lon_minmax))
        if lat_minmax[1][0] < lat_minmax[0][0]:
          lat_minmax = list(reversed(lat_minmax))

        lat_indices = np.arange(lat_minmax[0][0], lat_minmax[1][0]+1)
        lon_indices = np.arange(lon_minmax[0][0], lon_minmax[1][0]+1)

        if self.verbose == True:
          print("Area spans "+str(lat_minmax[0][1])+ " to "+str(lat_minmax[1][1]) + " and " + \
            str(lon_minmax[0][1]) + " to " + str(lon_minmax[1][1]))
        if print_report == True:
          for i in lat_indices:
            for j in lon_indices:
              x.append(self.lon0[j])
              y.append(self.lat0[i])

      elif len(self.lat0.shape) == 2:
        curv_grid = self.make_curv_grid()
        lat_indices = []
        lon_indices = []
        for i in range(self.lat0.shape[0]):
          for j in range(self.lat0.shape[1]):
            pair = curv_grid[i,j]
            if ((pair[0] >= lat_bounds[0]) and (pair[0] <= lat_bounds[1]) and \
              (pair[1] >= lon_bounds[0]) and (pair[1] <= lon_bounds[1])):
              lat_indices.append(i)
              lon_indices.append(j)
              if print_report == True:
                x.append(pair[1])
                y.append(pair[0])

      if nsubmit == 2:
        lat_indices1, lon_indices1 = lat_indices, lon_indices
        lon_bounds = lon_bounds2
      elif ((nsubmit == 1) and (resubmit == True)):
        lat_indices2, lon_indices2 = lat_indices, lon_indices
        lat_indices= np.concatenate((lat_indices1, lat_indices2))
        lon_indices= np.concatenate((lon_indices1, lon_indices2))

      nsubmit += -1

    if(print_report == True):
      import matplotlib.pyplot as plt
      m.scatter(x, y, s = 1.5,color='b', latlon = True, zorder = 10)
      f1 = self.folder_path(directory)+'bounded_area.'+scratchId+'.png'
      print("Bounded area on grid printed as map")
      plt.savefig(f1, dpi = 200)
      plt.close()

    return lat_indices, lon_indices

  def setup_axes(self):
    """
    Reorders variable axes so that time axis is first
    """

    if self.axes['time'] > 0:
      varF = np.moveaxis(self.var, self.axes['time'], 0)
      self.dim.remove(self.time_name)
      dimF = [self.time_name]+self.dim
    elif self.axes['time'] == 0:
      varF = self.var
      dimF = self.dim

    self.axes['time'] = dimF.index('time')
    self.axes['lat'], self.axes['lon'] = dimF.index('lat'), dimF.index('lon')

    if len(self.dim) > 3:
      self.axes['level'] = dimF.index('level')
    else:
      self.axes['level'] = False

  def make_curv_grid(self):
    """
    Returns a curvilinear grid from lat(i,j) and lon(i,j) pairs.
    """

    ni = self.lat0.shape[0]
    nj = self.lat0.shape[1]
    curv_grid = np.zeros((ni, nj, 2))

    for i in range(ni):
      for j in range(nj):
        curv_grid[i,j] = [self.lat0[i, j], self.lon0[i,j]]

    return curv_grid

  def timesync(self, **kwargs):
    """
    Returns a tuple of the main variable array starting from yr0
    and the offset index where the variable has been sliced

    Kwargs:
    * yr0 (int) [optional]: If yr0 has already been assigned to the query
      instance, yr0 kwarg will be ignored. If not, timesync() will look for
      a kwarg. If not specified, timesync() will use data_yr0 as yr0
    * yrf (int) [optional]: If yrf has already been assigned, yrf kwarg 
      will be ignores. If no yrf specified here or elsewhere, yrf set to None
      and ocrtools will generally just work with yr0 to the end of the datafile
    """

    try:
      self.yr0
    except AttributeError:
      try:
        self.yr0 = kwargs['yr0']
      except KeyError:
        self.yr0 = self.data_yr0[0]
        print(self.yr0)
    try:
      self.yrf
    except AttributeError:
      try:
        self.yrf = kwargs['yrf']
      except KeyError:
        self.yrf = None

    try:
      self.var
    except:
      self.extract_data()

    self.ndiv = ndivs[self.dt]
    print(self.yr0)
    print(self.data_yr0)
    offset = (self.yr0-self.data_yr0[0])*self.ndiv
    try:
      self.nt = (self.yrf-self.yr0+1)*self.ndiv
    except TypeError:
      self.nt = ((self.var.shape[0]-(self.yr0-self.data_yr0[0])*self.ndiv)-1)
      self.yrf = self.yr0-1+int(self.nt/self.ndiv)

    print(self.yrf)
    print('offset')
    print(offset)
    print('nt')
    print(self.nt)
    try:
      return self.var[offset:], offset
    except:
      raise

  def fill_cesm_params(self, **kwargs):

    if ((self.mem != 1850) and (self.yr0 < 2006) and (self.yrf < 2006)):
      self.data_yr0, self.data_yrf = [1920], [2005]
    elif ((self.mem != 1850) and (self.yr0 >= 2006) and (self.yrf >= 2006)):
      self.data_yr0, self.data_yrf = [2006], [2080]
    elif ((self.mem == 1850) and (self.yr0 < 2000) and (self.yrf < 2000)):
      self.data_yr0, self.data_yrf = [1900], [1999]
    elif ((self.mem == 1850) and (self.yr0 >= 2000) and (self.yrf >= 2000)):
      self.data_yr0, self.data_yrf = [2000], [2099]
    elif ((self.mem == 1850) and (self.yr0 < 2000) and (self.yrf >= 2000)):
      self.data_yr0, self.data_yrf = [1900, 2000], [1999, 2099]
    elif ((self.mem != 1850) and (self.yr0 < 2006) and (self.yrf >= 2006)):
      self.data_yr0, self.data_yrf = [1920, 2006], [2005, 2080]

    print(self.data_yr0)
    print(self.data_yrf)
    self.file = []
    for m in range(len(self.data_yr0)):

      self.time_name = 'time'
      if '_d' in self.var_name:
        self.var_name = self.var_name.replace('_d', '')

      if self.var_name in cice_vars:
        self.cesm_grid, self.model = 'cice', 'ice'
        self.lon_name = 'TLON'
        self.lat_name = 'TLAT'
        self.cellarea_name = 'tarea'
        f_hem = "_"+self.hemisphere
        if self.dt == 'daily':
          self.src_var_name = self.var_name+'_d'
          f_h = "h1"
          f_date = str(self.data_yr0[m]) + '0101-' + str(self.data_yrf[m]) + '1231'
        elif self.dt == 'monthly':
          self.src_var_name = self.var_name
          f_h = "h"
          f_date = str(self.data_yr0[m]) + '01-' + str(self.data_yrf[m]) + '12'

      elif self.var_name in cam_vars:
        self.cesm_grid, self.model = 'cam', 'atm'
        self.lon_name = 'lon'
        self.lat_name = 'lat'
        f_hem = ""
        self.src_var_name = self.var_name
        self.level_name = 'lev'
        if self.dt == 'daily':
          f_h = "h1"
          f_date = str(self.data_yr0[m]) + '0101-' + str(self.data_yrf[m]) + '1231'
        elif self.dt == "monthly":
          f_h = "h0"
          f_date = str(self.data_yr0[m]) + '01-' + str(self.data_yrf[m]) + '12'

      if (self.mem == 1850):
        f0 = 'b.e11.B1850C5CN.f09_g16.005.' + self.cesm_grid

      else:
        fmem = '{:03d}'.format(self.mem)
        if self.data_yr0[m] < 2006:
          f0 = 'b.e11.B20TRC5CNBDRD.f09_g16.' + fmem + '.' + self.cesm_grid
        elif self.data_yr0[m] >= 2006:
          f0 = 'b.e11.BRCP85C5CNBDRD.f09_g16.' + fmem + '.' + self.cesm_grid

      try:
        self.dim
      except AttributeError:
        self.dim = ['time', 'lat', 'lon'] if self.ndim == 3 else ['time', 'level', 'lat', 'lon']

      f2 = [f0, f_h, self.src_var_name+f_hem, f_date, 'nc']
      f = '.'.join(f2)
      folder_path = self.folder_path("cesm-raw")
      self.file.append(folder_path+f)
      print(self.file)

  def set_yr0(self, yr0):
    """
    Can be used to set yr0 of query object at any time

    Args:
    yr0 (int): calendar year 0 of query analysis
    """

    try:
      self.yr0 = yr0
    except:
      raise

  def set_yrf(self, yrf):
    """
    Can be used to set yrf of query object at any time

    Args:
    yrf (int): calendar year F (inclusive) of query analysis
    """
    try:
      self.yrf = yrf
    except:
      raise

  ## More data-wrangling utilities ##
  ###################################

  def bound_var(self, var, lat_indices, lon_indices, lat_axis = 1, lon_axis = 2):
    ## trims variable down
    ## to area covered by lon and lat indices

    lon_bounded_var = np.take(var, lon_indices, axis = lon_axis)
    bounded_var = np.take(lon_bounded_var, lat_indices, axis = lat_axis)
    return bounded_var

  def reg_wgt(self, latmin, latmax, nlat):
    ## returns a list of weighted values (sum = 1)
    ## based on zonal parameters

    dy = (latmax-latmin)/(nlat-1)
    lats = np.arange(latmin,latmax+dy,dy)[:nlat]
    wgts = np.zeros(nlat)
    if latmin == -90:
      for i in range(nlat):
        if ((i != 0) & (i != nlat-1)):
          wgts[i] = np.abs(np.sin(np.deg2rad(lats[i]+(dy/2))) - np.sin(np.deg2rad(lats[i]-(dy/2))))
        else:
          wgts[i] = 1-np.abs(np.sin(np.deg2rad(lats[i]+(dy/2))))
    else:
      for i in range(nlat):
        wgts[i] = np.abs(np.sin(np.deg2rad(lats[i]+(dy/2))) - np.sin(np.deg2rad(lats[i]-(dy/2))))

    wgts = wgts/(np.sum(wgts))
    return(wgts)

  def assign_reporting_vars(self, kwargs0):
    try:
      custom_tags = kwargs0['custom_tags']
    except KeyError:
      custom_tags = []
    try:
      print_report = kwargs0['print_report']
    except:
      print_report = True
    try:
      add_to_report = kwargs0['add_to_report']
    except:
      add_to_report = ''
    return custom_tags, print_report, add_to_report

  def find_nearest(self, array,value):
    idx = (np.abs(array-value)).argmin()
    return (idx, array[idx])

  ## Organizational tools ##
  ##########################

  def folder_path(self, directory):
    if self.directories != None:

      subs = {"dt": self.dt, "var_name": self.var_name, "proc": "proc", \
      "tseries": "teseries", "src_var_name": self.src_var_name}
      try:
        subs["yr_range"] = str(self.yr0)+'-'+str(self.yrf)
      except:
        pass
      try:
        subs["function"] = self.function
      except:
        pass
      if self.src == 'cesm':
        subs["mem"] = self.mem
        subs["model"] = self.model
        if self.var_name in cice_vars:
          subs["hemisphere"] = self.hemisphere

      subfolders = [subs[n] for n in self.subfolders[directory]]
      fpath = '/'.join(([self.directories[directory]]+subfolders))+'/'
      return fpath
    else:
      return ''

  def outfile(self, mode = '', custom_tags = [], **kwargs):
    abbrev = {'daily': 'd', 'monthly': 'mon'}
    exts = {'reduce': 'nc', 'report': 'txt', 'spatial_average-plot': 'png',
    'spatial_average-csv': 'csv'}

    subs = {'var_name': self.var_name, "dt": abbrev[self.dt], "mean": "mean", "scratchId": scratchId}

    if ((mode == 'reduce') or (mode == 'report')):
      if mode == 'reduce':
        if self.time_as == 'sequence':
          subs["time_slice"] = "base"+str(self.base_yr)+"."+ str(kwargs["div_since_base_yr"])
        elif self.time_as == 'date':
          if self.dt == 'monthly':
            subs["time_slice"] = '{:04d}'.format(kwargs["fYr0"]) + '{:02d}'.format(kwargs["fDiv0"]) + \
            '-' + '{:04d}'.format(kwargs["fYrF"]) + '{:02d}'.format(kwargs["fDivF"])
          elif self.dt == 'daily':
            subs["time_slice"] = '{:04d}'.format(kwargs["fYr0"]) + '{:03d}'.format(kwargs["fDiv0"]) + \
            '-' + '{:04d}'.format(kwargs["fYrF"]) + '{:03d}'.format(kwargs["fDivF"])

      if ((self.src == 'cesm')):
        name_type, folder_type = "cesm-reformatted", "cesm-reformatted"
        subs["mem"] = 'm'+'{:03d}'.format(self.mem)
      else:
        name_type, folder_type = "other-reformatted", "other-reformatted"

      if mode == 'report':
        name_type = "cesm_report" if self.src == 'cesm' else "other-report"

    elif (mode == 'spatial_average-plot' or mode == 'spatial_average-csv'):
      subs["time_slice"] =str(self.yr0)+'-'+str(self.yrf)
      name_type = 'spatial_average'
      if self.src == 'cesm':
        self.naming[name_type].insert(2, 'mem')
        subs["mem"] = 'm'+'{:03d}'.format(self.mem)

      if mode == 'spatial_average-plot':
        folder_type = 'plot'
      elif mode == 'spatial_average-csv':
        folder_type = 'csv'


    outfile_name = '.'.join([subs[seg] for seg in self.naming[name_type]]+custom_tags+[exts[mode]])
    outfile_folder = self.folder_path(folder_type)
    self.mkdir_p(outfile_folder)

    return outfile_folder+outfile_name

  def reduce_report(self):
    report_yr00, report_div00 = self.yr0, 0
    nfiles = int((self.yrf-self.yr0+1)*self.ndiv/120)
    yrf_div = nfiles*120+self.yr0*self.ndiv
    report_yrf00, report_divf00 = divmod(yrf_div, self.ndiv)

    report_yr0, report_yrf = '{:04d}'.format(report_yr00), '{:04d}'.format(int(report_yrf00))
    if self.ndiv == 365:
      report_div0, report_divf = '{:03d}'.format(report_div00), '{:03d}'.format(report_divf00)
    elif self.ndiv == 12:
      report_div0, report_divf = '{:02d}'.format(report_div00), '{:02d}'.format(report_divf00)

    date_time = now.strftime("%Y-%m-%d %H:%M")
    report = '\nData reformatting executed by the OCR Tools preanalysis module, version '+version
    report += ', on '+date_time+'. \n\nOutput '+str(nfiles) + ' files'
    report += ' covering the timespan '+report_yr0+'-'+report_div0 + ' to '+ report_yrf + '-' + report_divf + '.'
    report += '\n\nMain variable - ' + self.var_name + ' - as a function of ' +  \
    ', '.join([x for x in self.dim]) + '.'
    report_tag = ['_' +  report_yr0+ report_div0+ '-' + report_yrf+ report_divf]
    return report, report_tag

  def mkdir_p(self, path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

  def name_sync(self, source_file, keywords, replacements):
    ## Utility that takes scratch files produced with
    ## scratchId and syncs their names to the main
    ## output figure
    parent_folder = os.path.dirname(source_file)
    for file in os.listdir(parent_folder):
      if ((scratchId in file) and (replacements[0] in file)):
        target_file = parent_folder+'/'+file
        print(target_file)
    for i in range(len(keywords)):
      source_file = source_file.replace(keywords[i], replacements[i])

    new_name = source_file
    os.rename(target_file, new_name)

  def get_subs(self, a_dir):
    if os.path.isdir(a_dir):
      return [name for name in os.listdir(a_dir)
              if os.path.isdir(os.path.join(a_dir, name))]
    else:
      return []

  def get_ncs(self, a_dir):
    all_ncs = []
    for file in os.listdir(a_dir):
      if file.endswith(".nc"):
        all_ncs.append(file)
    return all_ncs

  def simple_params(self):
    self.lat_name, self.lon_name, self.time_name = 'lat', 'lon', 'time'
    self.cellarea_name = 'cellarea'
    self.dim = ['time', 'lat', 'lon']

    try:
      self.src_var_name = self.var_name
    except:
      pass

