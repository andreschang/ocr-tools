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
from scipy import stats
from scipy import signal
import matplotlib.pyplot as plt
import random as rand

## Global variables

ndivs = {'daily': 365, 'monthly':12}
cice_vars = ['aice', 'hi', 'flwdn', 'fswdn']
cam_vars = ['TS', 'PRECT']
now = datetime.now()
scratchId = now.strftime("%Y%m%d%H%M")

class build(object):

  def __init__(self, query, base_list, stage = None, verbose = True, dt = 'monthly', **kwargs):

    self.base_list = base_list
    self.verbose = verbose
    self.dt = dt
    self.ndiv = ndivs[dt]
    try:
      self.dt = kwargs['dt']
    except:
      try:
        self.dt = query.dt
      except InputError:
        print("Could not autofill dt. Please specify dt in build object")

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
    try:
      self.combine_steps = kwargs['combine_steps']
    except:
      pass

    ## Modify the input query class for naming
    self.query = query
    self.query.src = 'build'

  def annual_cycle(self, **kwargs):
    try:
      list0 = kwargs['list']
    except:
      list0 = self.base_list
    try:
      combine_steps = kwargs['combine_steps']
    except:
      try:
        combine_steps = self.combine_steps
      except:
        combine_steps = 1

    all_step0 = np.arange(0, self.ndiv, combine_steps)
    pad_list = np.concatenate((list0, list0[-self.ndiv:]))
    all_vals, average_steps, div_average_steps = [], [], []

    for i in all_step0:
      # print('step')
      # print(i)
      ## ind0: list of starting segment indices (ex. january year 0 OR january+february year 0)
      ## indf: list of ending segment indices (ex. february year f OR march+april year f)
      ind0, indf, all_ind, nstep = self.get_build_indices(list0, i, combine_steps)
      all_vals.append([pad_list[int(k)] for k in all_ind])

    for j in range(len(all_vals)-1):
      average_steps.append(np.mean(all_vals[j+1])-np.mean(all_vals[j]))
    average_steps.append(np.mean(all_vals[0])-np.mean(all_vals[-1]))

    for i in range(len(all_vals)):
      for j in range(combine_steps):
        div_average_steps.append(average_steps[i])

    new_list = [list0[0]]
    for i in range(len(list0)-1):
      year, div = divmod(i, self.ndiv)
      which_step = div
      step_i = div_average_steps[which_step]+new_list[i]
      new_list.append(step_i)

    print('Annual cycle of length '+str(len(new_list))+' calculated from list with combine_steps = '+str(combine_steps))

    return new_list

  def get_variance(self, **kwargs):
    try:
      list0 = kwargs['list']
    except:
      list0 = self.base_list
    try:
      combine_steps = kwargs['combine_steps']
    except:
      try:
        combine_steps = self.combine_steps
      except:
        combine_steps = 1

    self.ndiv = ndivs[self.dt]
    all_step0 = np.arange(0, self.ndiv, combine_steps)
    pad_list = np.concatenate((list0, list0[-self.ndiv:]))
    all_blur_vars, all_step_vars  = [], []

    for i in all_step0:
      ind0, indf, all_ind, nstep = self.get_build_indices(list0, i, combine_steps)
      div_vals = [pad_list[int(k)] for k in all_ind]

      ## get blur variances for each div
      all_blur_var = []
      all_blur_means = []

      for j in range(nstep+1):
        select_blur_range = np.arange(j*combine_steps, (j+1)*combine_steps)
        # print(select_blur_range)
        # print([div_vals[z] for z in select_blur_range])
        step_vals = [div_vals[z] for z in select_blur_range]
        all_blur_var.append(np.var(step_vals))
        all_blur_means.append(np.mean(step_vals))

      all_blur_vars.append(all_blur_var)
      all_step_vars.append(np.var(signal.detrend(all_blur_means)))
      # print(all_blur_vars)
      # print(np.var(signal.detrend(all_blur_means)))

    
    print('Blur var of length '+str(len(all_blur_vars))+', with sublist of length '+str(len(all_blur_vars[0])) +\
      ', and step var of length '+str(len(all_step_vars))+' calculated')


    return all_step_vars, all_blur_vars

  def new(self, list2, **kwargs):

    try:
      combine_steps = kwargs['combine_steps']
    except:
      try:
        combine_steps = self.combine_steps
      except:
        combine_steps = 1
    try:
      print_report = kwargs["print_report"]
    except:
      print_report = False
    try:
      step_var_a = kwargs["step_var_a"]
    except:
      step_var_a = 10.
    try:
      blur_var_a = kwargs["blur_var_a"]
    except:
      blur_var_a = 10.
    try:
      var_min = kwargs["var_min"]
    except:
      var_min = None
    try:
      var_max = kwargs["var_max"]
    except:
      var_max = None
    try:
      a = kwargs['a']
    except:
      a = 1
    try:
      snap = kwargs["snap"]
    except:
      snap = 20.

    list1 = self.base_list

    list1_base = self.annual_cycle(list = list1, combine_steps = combine_steps)
    list2_base = self.annual_cycle(list = list2, combine_steps = combine_steps)
    ## step var is detrended variance of each div (ex. january or mean(january, february))
    ## blur var is variance of each "combine_steps" segment (ex. january (0) or var(january, february) year 0)

    list1_step_var, list1_blur_var = self.get_variance(list = list1, combine_steps = combine_steps)
    list2_step_var, list2_blur_var = self.get_variance(list = list2, combine_steps = combine_steps)

    ## Make an average of the list1 and list2 bases
    list0_base = np.mean([list1_base, list2_base], axis = 0)
    max1, max2 = np.amax(list1), np.amax(list2)
    min1, min2 = np.amin(list1), np.amin(list2)
    max0 = max1 if max1>max2 else max2
    min0 = min1 if min1<min2 else min2
    data_range = max0-min0

    if print_report == True:
      print("List 1 step_var: ")
      print(list1_step_var)
      print("List 1 blur_var: ")
      print(list2_blur_var)
      print("List 2 step_var: ")
      print(list2_step_var)
      print("List 2 blur_var: ")
      print(list2_blur_var)

      x = np.arange(0, len(list1_base))
      plt.plot(x, list1_base, label = 'list 1')
      plt.plot(x, list2_base, label = 'list 2')
      plt.legend(loc = 'best')
      f1 = self.query.outfile(mode = 'spatial_average-plot', custom_tags = ['annual_cycle'])
      plt.savefig(f1, dpi = 200)
      plt.close()

    all_step0 = np.arange(0, self.ndiv, combine_steps)
    pad_list1 = np.concatenate((list1, list1[-self.ndiv:]))
    pad_list2 = np.concatenate((list2, list2[-self.ndiv:]))

    full_steps, opt_steps = [], []
    div_opt_steps = []

    for i in all_step0:
      for j in range(combine_steps):
        ## ind0: list of starting segment indices (ex. january year 0 OR january+february year 0)
        ## indf: list of ending segment indices (ex. february year f OR march+april year f)
        ind0, indf, all_ind, nstep = self.get_build_indices(list1, i, combine_steps)
        all_val1, all_val2 = [pad_list1[int(k)] for k in all_ind], [pad_list2[int(h)] for h in all_ind]
        full_step = a*(np.mean([pad_list2[k] for k in indf])-np.mean([pad_list1[h] for h in ind0]))
        opt_step = full_step/(len(list1)-1)

        full_steps.append(full_step)
        opt_steps.append(opt_step)
        if (print_report == True and j == 0):
          print('ind0, indf, nstep')
          print(ind0, indf, nstep)
          print('Build indices:')
          print(all_ind)
          print('Full step')
          print(full_step)
          print('opt_step')
          print(opt_step)

    # print('step_var, blur_var')
    # print(list1_step_var, list1_blur_var)

    new_list = [list1[0]]
    ghost_list = [list1[0]]
    list2_blur_var.append(list2_blur_var[0])

    for i in range(len(list1)-1):
      which_year, which_div = divmod(i, self.ndiv)
      which_blur = int(which_div/combine_steps)+1
      if which_blur == len(list2_step_var):
        which_blur = 0
      # print(which_year)
      # print(which_blur)

      step_dev = (list2_step_var[which_blur])**0.5
      blur_dev = ((list2_blur_var[which_blur][which_year])**0.5+(list2_blur_var[which_blur+1][which_year])**0.5)/2
      # print(div_opt_steps[which_div])
      step_0 = list0_base[i+1]
      # step_i = (list1_base[i+1]-list1_base[i])+opt_steps[which_div]+new_list[i]
      step_i = (list1_base[i+1]-list1_base[i])+opt_steps[which_div]+ghost_list[i]
      step_i = step_i + np.random.normal(0,step_dev)*(step_var_a/50.) + np.random.normal(0,blur_dev) * (blur_var_a/50.)
      if var_min != None:
        if step_i < var_min:
          step_i = var_min
      if var_max != None:
        if step_i > var_max:
          step_i = var_max

      ghost_list.append(step_i)

      snap_c0 = step_dev/data_range
      if snap_c0 >= snap/100.:
        snap_c = 1.
      else:
        snap_c = snap_c0/(snap/100.)

      step_i = (1-snap_c)*list0_base[i+1]+(snap_c)*step_i

      ## do some displacement here
      new_list.append(step_i)

    return new_list

  def magnet(self, list0, listM, **kwargs):
    try:
      magnet_a= kwargs["magnet_a"]
    except:
      magnet_a = 50.
    try:
      which_divs = kwargs["divs"]
    except:
      which_divs = np.arange(0, self.ndiv)

    new_list = list0

    for i in range(len(list0)-1):
      which_year, which_div = divmod(i, self.ndiv)
      if which_div in which_divs:

        new_list[i] = (list0[i]*(100.-magnet_a)/100.) + (listM[i]*magnet_a/100.)

    return new_list



  def get_build_indices(self, list0, div, combine_steps):
    yrf, yrf_div = divmod(len(list0), self.ndiv)

    ind0 = np.arange(div, div+combine_steps)
    ## -1 or -2 ???
    divf = (yrf-1)*self.ndiv+(div+combine_steps)
    # divf = (yrf-1)*self.ndiv+(div+combine_steps)
    if yrf_div >= (div+2*combine_steps):
      divf = divf+self.ndiv
      # print(divf)
    indf = np.arange(divf, divf+combine_steps)
    nstep = int((divf-div-1)/self.ndiv)
    # print(nstep)
    # print(len(ind0))
    all_ind = []
    for j in range(nstep+1):
      for g in range(combine_steps):
        all_ind.append(div+self.ndiv*j+g)
    return ind0, indf, all_ind, nstep



