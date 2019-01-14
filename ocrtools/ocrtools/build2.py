#!/usr/bin/env python
# -*- coding: utf-8 -*-

########################################
#
#  OCRTOOLS - Making climate data malleable
#  Copyright (C) 2018 Andres Chang
#
########################################

import xarray as xr
import pandas as pd
from scipy import signal
import numpy as np


def annual_cycle(data, dt, combine_steps=1):
    """
    Returns average annual cycle (i.e. climatology) of input dataset

    Args:
    * data (xarray dataset or data_array)
    * dt (str): monthly or daily
    * combine_steps (int): how many consecutive datapoints should
    be used for the analysis (ex. January 1-5, combine_steps = 5)
    """

    by, fby = get_groupings(dt)
    acycle = (data.resample(time=str(combine_steps)+fby, keep_attrs=True)
                  .mean()
                  .groupby(by)
                  .mean('time'))

    return(acycle)


def get_variance(data, dt, combine_steps=1):
    """
    Returns step_var - variance of each equivalent div (ex. each January) -
    & blur_var - variance within each div
    (only non-zero if combine_steps > 0)

    Args:
    * data (xarray dataset or data_array)
    * dt (str): monthly or daily
    * combine_steps (int): how many consecutive datapoints should
    be used for the analysis (ex. January 1-5, combine_steps = 5)
    """

    by, fby = get_groupings(dt)
    blur_var = (data.resample(time=str(combine_steps)+fby, keep_attrs=True)
                    .std()
                    .dropna('time', how='all'))**2

    def step_var(dataset):
        xr_vars = [x for x in dataset.data_vars]

        da = dataset.to_array().dropna('time', how='all')
        t_axis = da.get_axis_num('time')
        ddata = np.var(signal.detrend(da, axis=t_axis), axis=t_axis)
        dims = [x for x in dataset.dims if x != 'time']
        ddata_array = xr.Dataset(
            {xr_vars[0]: (dims, np.take(ddata, 0, axis=0))})
        for i in range(len(xr_vars) - 1):
            ddata_array[xr_vars[i + 1]] = xr.Dataset(
                {xr_vars[i + 1]: (dims, np.take(ddata, i + 1, axis=0))})
        return(ddata_array)

    step_var = (data.resample(time=str(combine_steps)+fby, keep_attrs=True)
                    .mean()
                    .groupby(by)
                    .apply(step_var))

    return(step_var, blur_var)


class build(object):

    def __init__(self, data1, dt, data2=None, **kwargs):
        """
        Makes new climate data based on existing modeled data and user options
        using a step-wise approach (i.e. OCR Tools calculates each timestep in
        sequence)
        """

        ac1 = annual_cycle(data1, dt, **kwargs)
        d1_step_var, d1_blur_var = get_variance(data1, dt, **kwargs)

        if data2 is not None:
            ac2 = annual_cycle(data2, dt, **kwargs)
            d2_step_var, d2_blur_var = get_variance(data2, dt, **kwargs)
            d12 = xr.concat([data1, data2], 'scenario')
            max0, min0 = np.amax(d12), np.amin(d12)
            self.base_data = d12.mean('scenario')
        else:
            self.base_data = ac1
            max0, min0 = np.amax(data1), np.amin(data1)
        d_range = max0 - min0




def get_groupings(dt):
    if dt == 'monthly':
        by = 'time.month'
        fby = 'MS'
    elif dt == 'daily':
        by = 'time.day'
        fby = 'D'
    return(by, fby)
