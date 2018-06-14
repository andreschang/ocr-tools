import numpy as np
import os, os.path
import errno
from datetime import datetime
from netCDF4 import Dataset

## Tools to improve malleability of climate data for 
## research and visualization purposes.
## 

## Global variables

version = "4.1"
ndivs = {'daily': 365, 'monthly':12}
cice_vars = ['aice', 'hi', 'flwdn', 'fswdn']
cam_vars = ['TS', 'PRECT']
now = datetime.now()
scratchId = now.strftime("%Y%m%d%H%M")

class reformat(object):

  ## Primary functions  ##
  ########################

  def __init__(self, stage = None):
    ## Preset specifies folder schema and location of CESM data
    # preset = 'andres_local'
    if preset == 'andres_local':
      self.directories = {"cesm-raw": "/Volumes/Samsung_T5/Open_Climate_Research-Projects/data/raw/cesm/", \
      "cesm-reformatted": "/Volumes/Samsung_T5/Open_Climate_Research-Projects/data/reformatted/cesm/", \
      "other-reformatted": "/Volumes/Samsung_T5/Open_Climate_Research-Projects/data/reformatted/", \
      "scratch": "/Volumes/Samsung_T5/Open_Climate_Research-Projects/data/scratch/"}
      self.cesm_raw_subfolders = ['dt', 'variable']
    elif preset == 'ucar':
      self.directories = {"cesm-raw": "/glade/p/cesmLE/CESM-CAM5-BGC-LE/", \
      "cesm-reformatted": "/glade/p/work/andresc/OCR/reformatted/cesm/", \
      "other-reformatted": "/glade/p/work/andresc/OCR/reformatted/", \
      "scratch": "/glade/p/work/andresc/OCR/scratch/"}
      self.cesm_raw_subfolders = ['model', 'proc', 'tseries', 'dt', 'var_name']
      import matplotlib
      matplotlib.use('Agg')
    

  def reduce_data(self, variable, yr0, yrf, dt = 'daily', mode = 'cesm', cesm_mem = 2, cesm_hemisphere = 'nh',
    file = '', data_yr0 = -1, lat_name = 'lat', lon_name = 'lon', level_name = 'level', cellarea_name = 'tarea',
    dim = ['time', 'lat', 'lon'], saveto = '', output = 'sequence', base_yr = 1920,
    custom_tags = [], report = 'yes', add_to_report = ''):

    ## modes: 'cesm', 'custom'
    ##
    ## Data from netCDF is sliced into manageable segments
    ## in files that cover the time period yr0 to yrf.
    ## Can be output as files with a sequential title (in days since base_yr)
    ## or as files with date range.
    ## This module assumes that input files begin on January 1, if dt = daily, and
    ## the month of January, if dt = monthly. Assumes 365 days/year if dt = daily.
    ## Items in dim must be named 'time', 'lat', 'lon', and 'level' (optional).
    ## These are used as uniform output dimensions, as well (with units preserved).
    ## Time and level assumed to be 1D arrays.
    ## Lat and lon can be 1D arrays or curvilinear. 
    ## 
    ## Details of 'cesm' mode:
    ## Required arguments are variable, yr0, yrf, dt, cesm_mem, and cesm_hemisphere (cice only).
    ## Variable should not include dt specifiers (ex. use 'aice', not 'aice_d').
    ## Passes arguments to grab_cesm_data, which returns necessary info.

    ## 0. EXTRACT RELEVANT DATA FROM FILE(S)

    if ((mode == 'cesm') and (variable in cice_vars)):
      custom_tags.append(cesm_hemisphere)

    var, units, dim, axes, lat0, lon0, level0, cellarea0, data_yr0 = self.extract_data(variable, \
      yr0, yrf, dt = dt, mode = mode, cesm_mem = cesm_mem, \
      cesm_hemisphere = cesm_hemisphere, file = file, data_yr0 = data_yr0, \
      lat_name = lat_name, lon_name = lon_name, level_name = level_name, cellarea_name = cellarea_name, \
      dim = ['time', 'lat', 'lon'])

    ndiv = ndivs[dt]
    nt, offset = self.timesync(yr0, yrf, dt = dt, data_yr0 = data_yr0)
    var = var[offset:]

    ## 2. MAKE SIMPLIFIED, SEGMENT FILES
    ##
    summary, sum_tag = self.reformat_summary(variable, yr0, yrf, dim, ndiv = ndivs[dt])

    for i in range(int(nt/120)):
      div_in_var = i*120+offset
      data = var[div_in_var: div_in_var+120]
      div_since_base_yr = div_in_var + (data_yr0-base_yr)*ndiv
      div_since_zero = div_in_var + data_yr0*ndiv

      fYr0, fDiv0 = divmod(div_since_zero, ndiv)
      fYrF, fDivF = divmod(div_since_zero+120, ndiv)

      div_format = str(int(div_since_base_yr))

      f_out = self.outfile(variable, output = output, dt = dt, mode = mode, mem = cesm_mem,
        base_yr = base_yr, div_since_base_yr = div_format, fYr0 = fYr0, fDiv0 = fDiv0, fYrF = fYrF,
        fDivF = fDivF, custom_tags = custom_tags, yr0 = yr0, yrf = yrf)

      dataset = Dataset(f_out, 'w', format='NETCDF4_CLASSIC')
      fdims = {}
      for index, dimension in enumerate(dim):
        fdims[dimension] = dataset.createDimension(dimension, data.shape[index])

      fvars = {}

      if len(lat0.shape) > 1:
        if axes['lat'] < axes['lon']:
          grid_axes = ('lat', 'lon')
        else:
          grid_axes = ('lon', 'lat')

      for index, dimension in enumerate(dim):
        if (((dimension == 'lat') or (dimension == 'lon')) and (len(lat0.shape) > 1)):
          fvars[dimension] = dataset.createVariable(dimension, np.float32, grid_axes)
        else:
          fvars[dimension] = dataset.createVariable(dimension, np.float32, (dimension,))

      if dt == 'daily':
        fYr0, fDiv0 = '{:04d}'.format(fYr0), '{:03d}'.format(fDiv0)
        fYrF, fDivF = '{:04d}'.format(fYrF), '{:03d}'.format(fDivF)
        print('\n formatting ... yr ' + str(fYr0) + ' day ' + str(fDiv0) + 
          '\n to yr ' + str(fYrF) + ' day ' + str(fDivF))
        fvars['time'].units = 'days since '+str(base_yr)+'-00-00'

      elif dt == 'monthly':
        fYr0, fDiv0 = '{:04d}'.format(fYr0), '{:02d}'.format(fDiv0)
        fYrF, fDivF = '{:04d}'.format(fYrF), '{:02d}'.format(fDivF)
        print('\n formatting ... yr ' + str(fYr0) + ' month ' + str(fDiv0) + 
          '\n to yr ' + str(fYrF) + ' month ' + str(fDivF))
        fvars['time'].units = 'months since '+str(base_yr)+'-00'

      fvars['time'][:] = np.arange(div_since_base_yr,div_since_base_yr+120)
      fvars['lat'].units = units['lat']
      fvars['lon'].units = units['lon']
      fvars['lat'][:] = lat0
      fvars['lon'][:] = lon0

      if 'level' in dim:
        fvars['level'][:] = level0
        fvars['level'].units = units['level']
      if len(lat0.shape) > 1:
        fvars['cellarea'] = dataset.createVariable('cellarea', np.float32, grid_axes)
        fvars['cellarea'][:] = cellarea0
        fvars['cellarea'].units = units['cellarea']

      saveVar = dataset.createVariable(variable,  np.float64, tuple(dim))
      saveVar.units = units['var']
      dataset.description = summary
      saveVar[:] = data
      dataset.close()

    ## 3. PRINT SUMMARY
    ##

    summary += '\n\n'+add_to_report
    sum_out = self.outfile(variable, output = 'summary', dt = dt, mode = mode, mem = cesm_mem,
      base_yr = base_yr, yr0 = yr0, yrf = yrf, custom_tags = sum_tag)

    np.savetxt(sum_out, [summary], fmt = '%s')

  def spatial_average(self, variable, yr0, yrf, dt = 'daily', mode = 'cesm', cesm_mem = 2, cesm_hemisphere = 'nh',
    file = '', data_yr0 = -1,  lat_name = 'lat', lon_name = 'lon', level_name = 'level', cellarea_name = 'tarea',
    dim = ['time', 'lat', 'lon'], lat_bounds = [-89., 89.], lon_bounds = [-179., 179.], saveto = '', custom_tags = [],
    lon_bound_convert = True, output = 'plot', print_summary = True):

    ## modes: 'cesm', 'custom'
    ## output: 'plot', 'list'
    ##
    ## Data from netCDF is averaged over space
    ## and output as a single list of values
    ## along a time axis.
    ## Lat indices can run -90 to 90
    ## Lon indices are input as -180 to 180 by default
    ## and converted to 0-360 for climate data.
    ## 
    ## Details of 'cesm' mode:
    ## Required arguments are variable, yr0, yrf, dt, cesm_mem, and cesm_hemisphere (cice only).
    ## Variable should not include dt specifiers (ex. use 'aice', not 'aice_d').
    ## Passes arguments to grab_cesm_data, which returns necessary info.

    ## 0. EXTRACT RELEVANT DATA FROM FILES
    ## 

    reformatted = False
    nsubmit = 1
    ndiv = ndivs[dt]
    xtime = np.arange(0, (ndiv)*(yrf-yr0+1))
    mvar = []

    if self.grab_reduced_data(variable, yr0, yrf, dt = dt, mode = mode, cesm_mem = cesm_mem, \
    cesm_hemisphere = cesm_hemisphere):
      target_ncs, start_index, end_index  = self.grab_reduced_data(variable, yr0, yrf, dt = dt, \
        mode = mode, cesm_mem = cesm_mem, cesm_hemisphere = cesm_hemisphere, print_summary = True)
      reformatted = True
      nsubmit = len(target_ncs)

    for i in range(nsubmit):

      if print_summary == True:
        print_summary0 = True if i == 0 else False
      else:
        print_summary0 = False

      if reformatted == True:
        file = target_ncs[i]
        var, units, dim, axes, lat0, lon0, level0, cellarea0, data_yr0 = self.extract_data(variable, \
          yr0, yrf, dt = dt, mode = 'custom', file = file, lat_name = 'lat', lon_name = 'lon', \
          level_name = 'level', cellarea_name = 'cellarea', dim = ['time', 'lat', 'lon'], \
          print_summary = print_summary0)
        if i == 0:
          nt, offset = 120-start_index, start_index
        elif i < (nsubmit-1):
          nt, offset = 120, 0
        elif i == (nsubmit-1):
          end_index = 120 if end_index == 0 else end_index
          nt, offset = end_index, 0

      else:
        var, units, dim, axes, lat0, lon0, level0, cellarea0, data_yr0 = self.extract_data(variable, \
        yr0, yrf, dt = dt, mode = mode, cesm_mem = cesm_mem, cesm_hemisphere = cesm_hemisphere, \
        file = file, data_yr0 = data_yr0, lat_name = lat_name, lon_name = lon_name, \
        level_name = level_name, cellarea_name = cellarea_name, dim = ['time', 'lat', 'lon'])

        nt, offset = self.timesync(yr0, yrf, dt = dt, data_yr0 = data_yr0)

      var = var[offset:offset+nt]

      ## 2. CRUNCH IT UP!
      ## Only works with 2d variables rn
      ##

      lat_indices, lon_indices = self.get_latlon_indices(lat0, lon0, lat_bounds, lon_bounds, print_summary = print_summary0)

      if len(lat0.shape) == 1:
        bounded_var = self.bound_var(var, lat_indices, lon_indices, lat_axis = axes['lat'], lon_axis = axes['lon'])
        mm_var = np.mean(bounded_var, axis = axes['lon'])

        lat_range = [lat0[k] for k in lat_indices]
        wgt = self.reg_wgt(lat_range[0], lat_range[-1], len(lat_range))
        mvar0 = 0
        for k in range(len(lat_range)):
          mvar0 += mm_var[:,k]*wgt[k]

      elif len(lat0.shape) == 2:
        total_area = 0.
        sum_var = 0.
        for i in range(len(lat_indices)):
          if np.isnan(var[0, lat_indices[i], lon_indices[i]]) == False:
            sum_var += var[:, lat_indices[i], lon_indices[i]]*cellarea0[lat_indices[i], lon_indices[i]]
            total_area += cellarea0[lat_indices[i], lon_indices[i]]
        mvar0 = sum_var/total_area

      mvar += list(mvar0)

    if output == 'plot':
      import matplotlib.pyplot as plt
      print(xtime.shape, len(mvar))
      plt.plot(xtime, mvar)
      f1 = outfile(variable, output = 'spatial_average-plot', dt = dt, mode = mode, mem = cesm_mem, yr0 = yr0,\
       yrf = yrf, custom_tags = custom_tags)
      plt.savefig(f1, dpi = 200)
      if print_summary == True:
        name_sync(f1, ['mean'], ['bounded_area'])

    elif output == 'csv':
      f1 = outfile(variable, output = 'spatial_average-csv', dt = dt, mode = mode, mem = cesm_mem, yr0 = yr0,\
       yrf = yrf, custom_tags = custom_tags)
      print(f1)
      np.savetxt(f1, mvar, delimiter = ',')
      if print_summary == True:
        name_sync(f1, ['mean','csv'], ['bounded_area','png'])

    elif output == 'list':
      return mvar

  ## Core data wrangling utilities ##
  ###################################
  def extract_data(self, variable, yr0, yrf, dt = 'daily', mode = 'cesm', cesm_mem = 2, cesm_hemisphere = 'nh',
    file = '', data_yr0 = -1, lat_name = 'lat', lon_name = 'lon', level_name = 'level', cellarea_name = 'cellarea',
    dim = ['time', 'lat', 'lon'], print_summary = True):

    resubmit = False
    if mode == 'cesm':
      cellarea_name = 'tarea'
      if ((cesm_mem != 1850) and (yr0 < 2006) and (yrf >= 2006)):
        resubmit = True
        yr0, yrf2 = yr0, yrf
        yrf, yr02 = 2005, 2006
      elif ((cesm_mem == 1850) and (yr0 < 2000) and (yrf >= 2000)):
        resubmit = True
        yr0, yrf2 = yr0, yrf
        yrf, yr02 = 1999, 2000

      file, data_yr0, var_name, lat_name, lon_name = grab_cesm_data(variable, yr0, yrf, \
        mem = cesm_mem, dt = dt, hemisphere = cesm_hemisphere, mode = 'reduce_data', \
        subfolders = self.cesm_raw_subfolders)

      if resubmit == True:
        file2, data_yr02, var_name, lat_name, lon_name = grab_cesm_data(variable, yr02, yrf2, \
          mem = cesm_mem, dt = dt, hemisphere = cesm_hemisphere, mode = 'reduce_data', \
          subfolders = self.cesm_raw_subfolders)
        yrf = yrf2

    else:
      var_name = variable
      if data_yr0 == -1:
        data_yr0 = yr0

    ## 1. LOAD DATA
    f_open = Dataset(file)
    if print_summary == True:
      print("\nShowing metadata ...")
      print(f_open.variables[var_name])
    axes = {}
    axes['time']  = dim.index('time')

    if resubmit == False:
      var = f_open.variables[var_name][:]
    elif resubmit == True:
      f_open2 = Dataset(file2)
      print("\nOpening two files...")
      var = np.concatenate((f_open.variables[var_name][:], f_open2.variables[var_name][:]), axis = axes['time'])

    var, dim, axes = setup_axes(var, dim)
    units = {}

    lat0 = f_open.variables[lat_name][:]
    lon0 = f_open.variables[lon_name][:]
    units['var'] = f_open.variables[var_name].units
    units['lat'] = f_open.variables[lon_name].units
    units['lon'] = f_open.variables[lat_name].units

    if len(lat0.shape) > 1:
      cellarea0 = f_open.variables[cellarea_name][:]
      units['cellarea'] = f_open.variables[cellarea_name].units
    else:
      cellarea0 = False 
    if 'level' in dim:
      level0 = f_open.variables[level_name][:]
      units['level'] = f_open.variables[level_name].units
    else:
      level0 = False

    f_open.close()

    return var, units, dim, axes, lat0, lon0, level0, cellarea0, data_yr0

  def grab_cesm_data(self, variable, yr0, yrf, mem = 2, dt = 'daily', hemisphere = 'nh',
    data_folder = "", subfolders = [], mode = 'getfile'):

    data_folder = self.directories['cesm-raw']

    if ((mem != 1850) and (yr0 < 2006) and (yrf < 2006)):
      data_yr0, data_yrf = 1920, 2005
    elif ((mem != 1850) and (yr0 >= 2006) and (yrf >= 2006)):
      data_yr0, data_yrf = 2006, 2080
    elif ((mem == 1850) and (yr0 < 2000) and (yrf < 2000)):
      data_yr0, data_yrf = 1900, 1999
    elif ((mem == 1850) and (yr0 >= 2000) and (yrf >= 2000)):
      data_yr0, data_yrf = 2000, 2099

    if variable in cice_vars:
      grid, model = 'cice', 'ice'
      lon_name = 'TLON'
      lat_name = 'TLAT'
      f_hem = "_"+hemisphere
      if dt == 'daily':
        var_name = variable+'_d'
        f_h = "h1"
        f_date = str(data_yr0) + '0101-' + str(data_yrf) + '1231'
      elif dt == 'monthly':
        var_name = variable
        f_h = "h"
        f_date = str(data_yr0) + '01-' + str(data_yrf) + '12'

    elif variable in cam_vars:
      grid, model = 'cam', 'atm'
      lon_name = 'lon'
      lat_name = 'lat'
      f_hem = ""
      var_name = variable
      if dt == 'daily':
        f_h = "h1"
        f_date = str(data_yr0) + '0101-' + str(data_yrf) + '1231'
      elif dt == "monthly":
        f_h = "h0"
        f_date = str(data_yr0) + '01-' + str(data_yrf) + '12'

    if (mem == 1850):
      f1 = 'b.e11.B1850C5CN.f09_g16.005.' + grid

    else:
      fmem = '{:03d}'.format(mem)
      if data_yr0 < 2006:
        f1 = 'b.e11.B20TRC5CNBDRD.f09_g16.' + fmem + '.' + grid
      elif data_yr0 >= 2006:
        f1 = 'b.e11.BRCP85C5CNBDRD.f09_g16.' + fmem + '.' + grid


    f0 = data_folder
    subs = {"dt": dt, "variable": variable, "hemisphere": hemisphere, "mem": str(mem), \
    "model": model, "proc": "proc", "tseries": "tseries", "var_name": var_name}
    for sub in subfolders:
      f0 += subs[sub] + '/'
    f = f0 + f1 + '.' + f_h + '.' + var_name + f_hem + '.' + f_date + '.nc'

    print("\nGrabbing "+dt+" "+variable+" data for "+str(yr0)+"-"+str(yrf) + " ... \n"+f)

    if mode == 'reduce_data':
      return f, data_yr0, var_name, lat_name, lon_name
    elif mode == 'getfile':
      return f

  def grab_reduced_data(self, variable, yr0, yrf, dt = 'daily', mode = 'cesm', cesm_mem = 'any',
    cesm_hemisphere = 'nh', print_summary = False):

    cesm_mem = 2

    if mode == 'cesm':
      f0 = self.directories["cesm-reformatted"]+variable+'/'
    else:
      f0 = self.directories["other-reformatted"]+variable+'/'

    # subdirs = [x for x in os.walk(f0)]
    subdirs = get_subs(f0)
    yr_ranges = [[int(x.split('-')[0]), int(x.split('-')[1])] for x in subdirs]
    for span in yr_ranges:
      if ((span[0] <= yr0) and (span[1]-1 >= yrf)):
        target_folder = f0+('-').join(map(str, span))
        rf_yr0 = span[0]
        if dt == 'daily':
          ndt = 365
        elif dt == 'monthly':
          ndt = 12
        nstartfile, start_index = divmod((yr0-rf_yr0)*ndt, 120)
        endfile, end_index = divmod((yrf-rf_yr0+1)*ndt, 120)
        all_ncs = get_ncs(target_folder)
        target_ncs = [target_folder+'/'+file for file in all_ncs[nstartfile:endfile+1]]
        if print_summary == True:
          print('Fetching reformatted data from '+target_folder)
        return target_ncs, start_index, end_index
      else:
        return False
    return False

  def setup_axes(self, var0, dim0):
    ## reorders var axes so that time axis = 0
    ## returns var, an ordered list of dimensions,
    ## and axis of each independent variable

    axes = {}
    axes['time']  = dim0.index('time')
    if axes['time'] > 0:
      varF = np.moveaxis(var0, axes['time'], 0)
      dim0.remove('time')
      dimF = ['time']+dim0
    elif axes['time'] == 0:
      varF = var0
      dimF = dim0

    axes['time'], axes['lat'], axes['lon'] = dimF.index('time'), dimF.index('lat'), dimF.index('lon')

    if len(dim0) > 3:
      axes['level'] = dimF.index('level')
    else:
      axes['level'] = False

    return varF, dimF, axes

  def get_latlon_indices(self, lat, lon, lat_bounds = [-89., 89.], lon_bounds = [-179., 179.], \
    lon_bound_convert = True, print_summary = False, data_wrap_lon = 0):
    ## returns lists of indices in area bounded by
    ## lat and lon bounds
    ##

    resubmit = False
    nsubmit = 1
    mlon_bounds = lon_bounds
    if lon_bound_convert == True:
      lon_bounds = [(bound+360. if bound < 0 else bound) for bound in lon_bounds]
      if print_summary == True:
        print("Converting lon bounds from 180W-180E to 0-360")
        print("New lon bounds: "+",".join(map(str, lon_bounds)))

    if lon_bounds[0] > lon_bounds[1]:
      if ((data_wrap_lon == 0) or (data_wrap_lon == 360)):
        wrap_points = [360., 0.]
      elif ((data_wrap_lon == -180) or (data_wrap_lon == 180)):
        wrap_points = [180., -180.]
      resubmit = True
      nsubmit = 2
      lon_bounds2 = [wrap_points[1], lon_bounds[1]]
      lon_bounds = [lon_bounds[0], wrap_points[0]]

    if print_summary == True:
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
      if len(lat.shape) == 1:
        lat_minmax = [find_nearest(lat, bound) for bound in lat_bounds]
        lon_minmax = [find_nearest(lon, bound) for bound in lon_bounds]

        if lon_minmax[1][0] < lon_minmax[0][0]:
          lon_minmax = list(reversed(lon_minmax))
        if lat_minmax[1][0] < lat_minmax[0][0]:
          lat_minmax = list(reversed(lat_minmax))

        lat_indices = np.arange(lat_minmax[0][0], lat_minmax[1][0]+1)
        lon_indices = np.arange(lon_minmax[0][0], lon_minmax[1][0]+1)

        if print_summary == True:
          print("Area spans "+str(lat_minmax[0][1])+ " to "+str(lat_minmax[1][1]) + " and " + \
            str(lon_minmax[0][1]) + " to " + str(lon_minmax[1][1]))
          for i in lat_indices:
            for j in lon_indices:
              x.append(lon[j])
              y.append(lat[i])

      elif len(lat.shape) == 2:
        curv_grid = make_curv_grid(lat, lon)
        lat_indices = []
        lon_indices = []
        for i in range(lat.shape[0]):
          for j in range(lat.shape[1]):
            pair = curv_grid[i,j]
            if ((pair[0] >= lat_bounds[0]) and (pair[0] <= lat_bounds[1]) and \
              (pair[1] >= lon_bounds[0]) and (pair[1] <= lon_bounds[1])):
              lat_indices.append(i)
              lon_indices.append(j)
              if print_summary == True:
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

    if(print_summary == True):
      import matplotlib.pyplot as plt
      m.scatter(x, y, s = 1.5,color='b', latlon = True, zorder = 10)
      f1 = self.directories['scratch']+'bounded_area.'+scratchId+'.png'
      print("Bounded area on grid printed as map")
      plt.savefig(f1, dpi = 200)
      plt.close()

    return lat_indices, lon_indices

  def make_curv_grid(self, lat, lon):
    ni = lat.shape[0]
    nj = lat.shape[1]
    curv_grid = np.zeros((ni, nj, 2))

    for i in range(ni):
      for j in range(nj):
        curv_grid[i,j] = [lat[i, j], lon[i,j]]

    return curv_grid

  ## More data wrangling utilities ##
  #####################################
  def timesync(self, yr0, yrf, dt = 'daily', data_yr0 = -1):
    ndiv = ndivs[dt]
    offset = (yr0-data_yr0)*ndiv
    nt = (yrf-yr0+1)*ndiv
    return nt, offset

  def bound_var(self, var, lat_indices, lon_indices, lat_axis = 1, lon_axis = 2):
    ## trims variable down
    ## to area covered by lon and lat indices

    lon_bounded_var = np.take(var, lon_indices, axis = lon_axis)
    bounded_var = np.take(lon_bounded_var, lat_indices, axis = lat_axis)
    return bounded_var

  def reg_wgt(latmin, latmax, nlat):
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

  def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return (idx, array[idx])

  def first_look(file):
    f_open = Dataset(file)
    print("\nShowing metadata ...")
    print(f_open.variables)

  ## Naming and documentation utilities ##
  ########################################
  def outfile(variable, output = 'sequence', dt = 'daily', mode = 'custom', mem = -1, 
    base_yr = 1920, div_since_base_yr = -1, fYr0 = -1, fDiv0 = -1, fYrF = -1, fDivF = -1, custom_tags = [],
    yr0 = -1, yrf = -1):

    abbrev = {'daily': 'd', 'monthly': 'mon'}
    exts = {'sequence': 'nc', 'date': 'nc', 'summary': 'txt', 'spatial_average-plot': 'png',
    'spatial_average-csv': 'csv'}
    f1 = [variable, abbrev[dt]]
    if mode == 'cesm':
      f0 = self.directories["cesm-reformatted"]
      f1.append('m'+'{:03d}'.format(mem))
    else:
      f0 = self.directories["other-reformatted"]
    f0 += variable + '/' + str(yr0)+'-'+str(yrf)+'/'

    if output == 'sequence':
      f1.append('base'+str(base_yr))
      f1.append(div_since_base_yr) 
    elif output == 'date':
      f1.append(str(fYr0) + str(fDiv0) + '-' + str(fYrF) + str(fDivF))
    elif('spatial_average' in output):
      f0 = self.directories["scratch"]
      f1.append(str(yr0)+'-'+str(yrf))
      f1.append('mean')
      f1.append(scratchId)

    if custom_tags != []:
      f1.extend(custom_tags)

    f1.append(exts[output])
    f_out = f0+('.').join(f1)
    mkdir_p(f0)

    return f_out

  def reformat_summary(variable, yr0, yrf, dim, ndiv = 365):
    sum_yr00, sum_div00 = yr0, 0
    nfiles = int((yrf-yr0+1)*ndiv/120)
    yrf_div = nfiles*120+yr0*ndiv
    sum_yrf00, sum_divf00 = divmod(yrf_div, ndiv)

    sum_yr0, sum_yrf = '{:04d}'.format(sum_yr00), '{:04d}'.format(int(sum_yrf00))
    if ndiv == 365:
      sum_div0, sum_divf = '{:03d}'.format(sum_div00), '{:03d}'.format(sum_divf00)
    elif ndiv == 12:
      sum_div0, sum_divf = '{:02d}'.format(sum_div00), '{:02d}'.format(sum_divf00)

    date_time = now.strftime("%Y-%m-%d %H:%M")
    summary = '\nData reformatting executed by the OCR Tools preanalysis module, version '+version
    summary += ', on '+date_time+'. \n\nOutput '+str(nfiles) + ' files'
    summary += ' covering the timespan '+sum_yr0+'-'+sum_div0 + ' to '+ sum_yrf + '-' + sum_divf + '.'
    summary += '\n\nMain variable - ' + variable + ' - as a function of ' +  \
    ', '.join([x for x in dim]) + '.'
    sum_tag = ['_' +  sum_yr0+ sum_div0+ '-' + sum_yrf+ sum_divf]
    return summary, sum_tag

  ## File organization utilities ##
  #################################
  def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

  def get_subs(a_dir):
    if os.path.isdir(a_dir):
      return [name for name in os.listdir(a_dir)
              if os.path.isdir(os.path.join(a_dir, name))]
    else:
      return []

  def get_ncs(a_dir):
    all_ncs = []
    for file in os.listdir(a_dir):
      if file.endswith(".nc"):
        all_ncs.append(file)
    return all_ncs

  def name_sync(source_file, keywords, replacements):
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


    