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
from ocrtools.tk_selector import plt
from ocrtools.explore import spatial_average
import numpy as np


def annual_cycle(data, dt, combine_steps=1, **kwargs):
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
                  .mean('time')
                  .groupby(by)
                  .mean('time'))

    return(acycle)


def get_stds(data, dt, combine_steps=1, **kwargs):
    """
    Returns step_std - standard dev of each equivalent div (ex. each January) -
    & blur_std - standard dev within each div
    (only non-zero if combine_steps > 0)

    Args:
    * data (xarray dataset or data_array)
    * dt (str): monthly or daily
    * combine_steps (int): how many consecutive datapoints should
    be used for the analysis (ex. January 1-5, combine_steps = 5)
    """

    by, fby = get_groupings(dt)
    blur_std = (data.resample(time=str(combine_steps)+fby, keep_attrs=True)
                    .std('time')
                    .dropna('time', how='all'))

    def step_std(dataset):
        xr_vars = [x for x in dataset.data_vars]

        da = dataset.to_array().dropna('time', how='all')
        t_axis = da.get_axis_num('time')
        ddata = np.std(signal.detrend(da, axis=t_axis), axis=t_axis)
        dims = [x for x in dataset.dims if x != 'time']
        ddata_array = xr.Dataset(
            {xr_vars[0]: (dims, np.take(ddata, 0, axis=0))})
        for i in range(len(xr_vars) - 1):
            ddata_array[xr_vars[i + 1]] = xr.Dataset(
                {xr_vars[i + 1]: (dims, np.take(ddata, i + 1, axis=0))})
        return(ddata_array)

    step_std = (data.resample(time=str(combine_steps)+fby, keep_attrs=True)
                    .mean('time')
                    .groupby(by)
                    .apply(step_std))

    return(step_std, blur_std)


class build(object):

    def __init__(self, data1, dt, data2=None, **kwargs):
        """
        Makes new climate data based on existing modeled data and user options
        using a step-wise approach (i.e. OCR Tools calculates each timestep in
        sequence)
        """

        # Set parameters
        try:
            combine_steps = kwargs["combine_steps"]
        except KeyError:
            combine_steps = 1
        try:
            step_std_a = kwargs["step_std_a"]
        except KeyError:
            step_std_a = 10.
        try:
            blur_std_a = kwargs["blur_std_a"]
        except KeyError:
            blur_std_a = 10.
        try:
            var_min = kwargs["var_min"]
        except KeyError:
            var_min = None
        try:
            var_max = kwargs["var_max"]
        except KeyError:
            var_max = None
        try:
            snap = kwargs["snap"]
            if snap > 15:
                snap = 15
        except KeyError:
            snap = 10.
        try:
            snap_atten = kwargs["snap_atten"]
            if snap_atten > 50:
                snap_atten = 50
        except KeyError:
            snap_atten = 20
        try:
            head = kwargs["head"]
        except KeyError:
            head = 1
        try:
            tail = kwargs["tail"]
        except KeyError:
            tail = 1
        try:
            load_rand = kwargs["load_rand"]
        except KeyError:
            load_rand = None
        try:
            savgol_window = kwargs["savgol_window"]
        except KeyError:
            savgol_window = 0
        try:
            hist_stretch = kwargs["hist_stretch"]
        except KeyError:
            hist_stretch = False
        try:
            c_dist = kwargs["hist_dist"]
        except KeyError:
            def c_dist(x):
                return(x)
        try:
            plot = kwargs["plot"]
        except KeyError:
            plot = False

        main_vars = data1.attrs['main_vars']
        nvars = len(main_vars)
        by, fby = get_groupings(dt)
        t_dim = by.split('.')[1]

        # Make annual cycle plots of data1 and data2 for diagnostics
        if plot:
            fig0, ax0 = plt.subplots(
                ncols=1, nrows=nvars, figsize=(8.5, 3 * nvars))

        def plot_sa(data, n=1):
            if len([x for x in data.dims]) > 1:
                sa = spatial_average(data, **kwargs)
            else:
                sa = data
            if nvars == 1:
                sa[main_vars[0]].plot(ax=ax0, label='Scenario '+str(n))
            else:
                for vi in range(nvars):
                    sa[main_vars[vi]].plot(ax=ax0[vi])

        # Get annual cycle and standard dev params of data1
        # calculate standard devs:
        # step std is detrended standard dev of each div
        # (ex. january or mean(january, february))
        # blur std is standard dev of each "combine_steps" segment
        # (ex. january (0) or std(january, february) year 0)
        ac1 = (data1.resample(time=str(combine_steps)+fby, keep_attrs=True)
                   .mean('time').dropna('time', how='all'))
        d1_step_std, d1_blur_std = get_stds(data1, dt, **kwargs)
        if plot:
            plot_sa(ac1, 1)

        # Do the same for data2 if given, otherwise refer always to data1
        if data2 is not None:
            ac2 = (data2.resample(time=str(combine_steps)+fby, keep_attrs=True)
                   .mean('time').dropna('time', how='all'))
            d2_step_std, d2_blur_std = get_stds(data2, dt, **kwargs)
            if plot: plot_sa(ac2, 2)

            d12 = xr.concat([ac1, ac2], 'scenario')
            max0, min0 = np.amax(d12), np.amin(d12)
            step_std = (d1_step_std + d2_step_std)/2
            blur_std = (d1_blur_std + d2_blur_std)/2
            self.base_data = d12.mean('scenario')
        else:
            self.base_data = ac1
            max0, min0 = np.amax(data1), np.amin(data1)
            step_std, blur_std = d1_step_std, d1_blur_std
        base_ac = self.base_data.groupby(by).mean('time')

        # Calculate full range of input data
        d_range = max0 - min0

        # Absolute 'full_step' is calculated, which is the change between
        # div(i~), yr0 and div(i~+1), yrf. Where i~ here means the full
        # 'combine_steps' div. Ex. January 1-5 1995, January 6-10 2000.
        # Opt steps is the full step divided by number of years

        time_range = data1.coords['time'].to_index()
        yr0, yrf = np.amin(time_range).year, np.amax(time_range).year
        d_yr0 = annual_cycle(
            data1.sel(time=slice("{:04d}".format(yr0) + '-01-01',
                                 "{:04d}".format(yr0 + head) + '-12-31')),
            dt, **kwargs)

        if data2 is not None:
            d_yrf = annual_cycle(
                data2.sel(time=slice("{:04d}".format(yrf - tail) + '-01-01',
                                     "{:04d}".format(yrf) + '-12-31')),
                dt, **kwargs)
        else:
            d_yrf = annual_cycle(
                data1.sel(time=slice("{:04d}".format(yrf - tail) + '-01-01',
                                     "{:04d}".format(yrf) + '-12-31')),
                dt, **kwargs)

        full_steps = d_yrf - d_yr0.roll({t_dim: 1}, roll_coords=True)
        opt_steps = full_steps/(yrf-yr0)
        self.rand = []
        self.snap_list = []
        out_vars = []

        for vi in range(nvars):
            var_i = main_vars[vi]
            new_var = (
                data1[var_i].copy(deep=True)
                            .resample(time=str(combine_steps) + fby,
                                      keep_attrs=True)
                            .mean('time')
                            .dropna('time', how='all')
                        )

            ctime = new_var.coords['time'].to_index()

            for i in range(len(new_var['time']) - 1):
                if dt == 'daily':
                    t_attr = 'dayofyear'
                elif dt == 'monthly':
                    t_attr = 'month'

                # Step is calculated based on previous plus opt_step
                new_step = (new_var.isel(time=i) + opt_steps[var_i].sel(
                        {t_dim: getattr(ctime[i], t_attr)}))
                # Some standard deviation added based on average of std for
                # step i and step i+1

                add_step_std = (step_std.sel(
                    {t_dim: [getattr(ctime[i], t_attr),
                                  getattr(ctime[i+1], t_attr)]})
                                        .mean(t_dim))
                add_blur_std = (blur_std.sel(
                    {'time': slice(ctime[i], ctime[i+1])})
                                        .mean('time'))

                # produce some randomness (or load it from previous)
                if load_rand is not None:
                    r1, r2 = load_rand[i]
                else:
                    if vi == 0:
                        r1, r2 = np.random.normal(), np.random.normal()
                        self.rand.append((r1, r2))
                    else:
                        r1, r2 = self.rand[i]

                new_step = (new_step + r1 * add_step_std * step_std_a/50. +
                             r2 * add_blur_std * blur_std_a/50.)

                # Normalize snapping envelope to 5 standard deviations
                dev_from_base = (new_step - base_ac[var_i].sel(
                    {t_dim: getattr(ctime[i+1], t_attr)})).apply(
                        np.abs)/(add_step_std * (15-snap))
                dev_from_base = dev_from_base.where(dev_from_base > 1, 1)
                snap_amt = dev_from_base**((50 - snap_atten)/5)
                self.snap_list.append(snap_amt)

                new_step = (
                    snap_amt * self.base_data[var_i].sel({'time': ctime[i+1]}) +
                    (1 - snap_amt) * new_step)

                # Replace variables that exceed min or max with min/max values
                new_step = filter_minmax(new_step, var_i, var_min, var_max)

                # modify the array!
                new_var = new_var.where(new_step != new_var, new_step)

            # enhance contrast
            if hist_stretch:
                in_mean = new_var.mean('time')
                in_std = new_var.std('time')
                contrast_lims = [in_mean-in_std, in_mean+in_std]
                c_range = contrast_lims[1]-contrast_lims[0]
                print(c_range)

                if var_min is not None:
                    out_min = var_min
                else:
                    out_min = min0
                if var_max is not None:
                    out_max = var_max
                else:
                    out_max = max0

                hist_var = (((new_var - contrast_lims[0]) /
                            c_range).apply(c_range) * (out_max - out_min)
                            + out_min)

                new_var = new_var.where(
                    new_var > contrast_lims[0] and new_var < contrast_lims[1],
                    hist_var)

            if combine_steps > 1:
                new_var = new_var.resample(time=fby).interpolate()

            new_var = new_var.to_array(name=var_i).squeeze()
            print(new_var)
            
            out_vars.append(new_var)
            if plot:
                if nvars == 1:    
                    spatial_average(new_var).plot(ax=ax0, label='new')
                    ax0.legend(loc='best')
                else:
                    spatial_average(new_var).plot(ax=ax0[vi], label='new')

        new_ds = xr.merge(out_vars)
        # plt.show()
        # plt.clf()
        self.new = new_ds


def filter_minmax(dataset, var_i, var_min, var_max):
    if var_min is not None:
        if isinstance(var_min, dict):
            try:
                dataset = dataset.where(
                    dataset < var_min[var_i], var_min[var_i])
            except KeyError:
                pass
        else:
            dataset = dataset.where(
                dataset < var_min, var_min)

    if var_max is not None:
        if isinstance(var_max, dict):
            try:
                if dataset < var_max[var_i]:
                    dataset = var_max
            except KeyError:
                pass
        else:
            if dataset < var_max:
                dataset = var_max

    return(dataset)


def get_groupings(dt):
    if dt == 'monthly':
        by = 'time.month'
        fby = 'MS'
    elif dt == 'daily':
        by = 'time.day'
        fby = 'D'
    return(by, fby)
