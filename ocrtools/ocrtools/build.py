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
            ddata_array[xr_vars[i+1]] = np.take(ddata, i+1, axis=0)
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
        main_vars = data1.attrs['main_vars']
        nvars = len(main_vars)
        by, fby = get_groupings(dt)
        t_dim = by.split('.')[1]

        try:
            combine_steps = kwargs["combine_steps"]
        except KeyError:
            combine_steps = 1
        try:
            a = kwargs["a"]
        except KeyError:
            a = 1.
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
            plot = kwargs["plot"]
        except KeyError:
            plot = False
        try:
            debug = kwargs["debug"]
        except KeyError:
            debug = False
        try:
            debug_step = kwargs["debug_step"]
        except KeyError:
            debug_step = 0

        step_std_a = try_dict(kwargs, 'step_std_a', 10., main_vars)
        blur_std_a = try_dict(kwargs, 'blur_std_a', 10., main_vars)
        var_min = try_dict(kwargs, 'var_min', None, main_vars)
        var_max = try_dict(kwargs, 'var_max', None, main_vars)
        snap = try_dict(kwargs, 'snap', None, main_vars)
        snap = {k: 15 if v > 15 else v for k, v in snap.items()}
        snap_atten = try_dict(kwargs, 'snap_atten', 20, main_vars)
        snap_atten = {k: 50 if v > 50 else v for k, v in snap.items()}
        hist_stretch = try_dict(kwargs, 'hist_stretch', False, main_vars)
        hist_args = try_dict(kwargs, 'hist_args', False, {})
        savgol_window = try_dict(kwargs, 'savgol_window', 1, main_vars)

        def none_dist(x):
            return(x)

        c_dist = try_dict(kwargs, 'hist_dist', none_dist, main_vars)


        # Make plots of data1 and data2 for comparison
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
                ax0.legend(loc='best')
                # ax0.set_title('OCR Build '+main_vars[0])
            else:
                for vi in range(nvars):
                    sa[main_vars[vi]].plot(ax=ax0[vi], label='Scenario '+str(n))
                    ax0[vi].legend(loc='best')
                    # ax0[vi].set_title('OCR Build '+main_vars[vi])


        def apply_savgol(dataset):
            for var in main_vars:
                dataset[var] = xr.DataArray(signal.savgol_filter(
                    dataset[var], savgol_window[var], 2,
                    axis=dataset[var].get_axis_num('time')),
                    coords=dataset[var].coords, dims=dataset[var].dims)
            return(dataset)


        # Get annual cycle and standard dev params of data1
        # calculate standard devs:
        # step_std is detrended standard dev of each div
        # (ex. january or mean(january, february))
        # blur_std is standard dev of each "combine_steps" segment
        # (ex. january (0) or std(january, february) year 0)
        d1 = (data1.resample(time=str(combine_steps)+fby, keep_attrs=True)
                   .mean('time').dropna('time', how='all'))
        d1 = apply_savgol(d1)

        d1_step_std, d1_blur_std = get_stds(data1, dt, **kwargs)
        if plot:
            plot_sa(d1, 1)

        # Do the same for data2 if given, otherwise refer always to data1
        if data2 is not None:
            d2 = (data2.resample(time=str(combine_steps)+fby, keep_attrs=True)
                       .mean('time').dropna('time', how='all'))
            d2 = apply_savgol(d2)
            d2_step_std, d2_blur_std = get_stds(data2, dt, **kwargs)
            if plot: plot_sa(d2, 2)

            d12 = xr.concat([d1, d2], 'scenario')
            max0, min0 = np.amax(d12), np.amin(d12)
            step_std = (d1_step_std + d2_step_std)/2
            blur_std = (d1_blur_std + d2_blur_std)/2
            self.base_data = d12.mean('scenario')
        else:
            self.base_data = d1
            max0, min0 = np.amax(data1), np.amin(data1)
            step_std, blur_std = d1_step_std, d1_blur_std
        base_ac = self.base_data.groupby(by).mean('time')
        d1_group = d1.groupby(by).mean('time')

        # Calculate full range of input data
        d_range = max0 - min0

        # opt_step is the variable delta between timesteps based on 1) standard
        # annual cycle of data1 and 2) affect of scenario choice on overall
        # variable between data1 and data12

        time_range = data1.coords['time'].to_index()
        yr0, yrf = np.amin(time_range).year, np.amax(time_range).year
        d_yr0 = annual_cycle(
            data1.sel(time=slice("{:04d}".format(yr0) + '-01-01',
                                 "{:04d}".format(yr0 + head-1) + '-12-31')),
            dt, **kwargs)

        if data2 is not None:
            d_yrf = annual_cycle(
                data2.sel(time=slice("{:04d}".format(yrf - tail+1) + '-01-01',
                                     "{:04d}".format(yrf) + '-12-31')),
                dt, **kwargs)
        else:
            d_yrf = annual_cycle(
                data1.sel(time=slice("{:04d}".format(yrf - tail+1) + '-01-01',
                                     "{:04d}".format(yrf) + '-12-31')),
                dt, **kwargs)

        d1_steps = d1_group.roll({t_dim: -1}, roll_coords=False) - d1_group
        full_steps = d_yrf.roll({t_dim: -1}, roll_coords=False) - d_yr0
        opt_steps = d1_steps + a * (full_steps - d1_steps)/(yrf - yr0 + 1)

        self.rand = []
        self.snap_list = []
        out_vars = []

        if debug:
            print('New loop', opt_steps)

        for vi in range(nvars):
            var = main_vars[vi]
            new_var = (
                data1[var].copy(deep=True).resample(time=str(combine_steps) + fby,
                                      keep_attrs=True)
                            .mean('time')
                            .dropna('time', how='all'))

            ctime = new_var.coords['time'].to_index()

            for i in range(len(new_var['time']) - 1):
                if dt == 'daily':
                    t_attr = 'dayofyear'
                elif dt == 'monthly':
                    t_attr = 'month'

                # Step is calculated based on previous plus opt_step

                new_step = ((new_var.isel(time=i).copy(deep=True) + opt_steps[var].sel(
                        {t_dim: getattr(ctime[i], t_attr)}, method='nearest'))
                                   .assign_coords(time=xr.DataArray(ctime[i+1])))

                if debug and i == debug_step:
                    print('\n[OCR debug] Date of loop: '+str(ctime[i]) +
                          ', Date of replacement '+str(ctime[i+1]))
                    print('[OCR debug] New_step from opt_step sum')
                    print(new_step)
                    print('from sum of')
                    print((new_var.isel(time=i).copy(deep=True)))
                    print('and')
                    print(opt_steps[var].sel(
                        {t_dim: getattr(ctime[i], t_attr)}, method='nearest'))

                # Some standard deviation added based on average of std for
                # step i and step i+1

                add_step_std = (step_std.sel(
                    {t_dim: [getattr(ctime[i], t_attr),
                             getattr(ctime[i+1], t_attr)]})
                                        .mean(t_dim))[var]

                add_blur_std = (blur_std.sel(
                    {'time': slice(ctime[i], ctime[i+1])})
                                        .mean('time'))[var]
                # produce some randomness (or load it from previous)
                if load_rand is not None:
                    r1, r2 = load_rand[i]
                else:
                    if vi == 0:
                        r1, r2 = np.random.normal(), np.random.normal()
                    else:
                        r1, r2 = self.rand[i]
                self.rand.append((r1, r2))

                new_step = (new_step + r1 * add_step_std * step_std_a[var]/50. +
                            r2 * add_blur_std * blur_std_a[var]/50.)

                if debug and i == debug_step:
                    print('[OCR debug] New_step after the addition of randomness')
                    print(new_step)

                # Normalize snapping envelope to 5 standard deviations default
                # unless snap is 0
                if snap[var] == 0:
                    snap_amt = 0
                else:
                    dev_from_base = (np.abs(
                        new_step - base_ac[var].sel(
                            {t_dim: getattr(ctime[i+1], t_attr)})) /
                        (add_step_std * (15-snap[var])))

                    dev_from_base = xr.where(dev_from_base > 1, 1, dev_from_base)
                    dev_from_base = xr.where(np.isnan(dev_from_base), 1, dev_from_base)
                    snap_amt = dev_from_base**((50 - snap_atten[var])/5)
                    new_step = (
                        snap_amt * self.base_data[var].sel({'time': ctime[i+1]}) +
                        (1 - snap_amt) * new_step)

                self.snap_list.append(snap_amt)

                if debug and i == debug_step:
                    print('[OCR debug] New_step after snapping')
                    print(new_step)
                    print(' by amount ')
                    print(snap_amt)

                # Replace variables that exceed min or max with min/max values
                new_step = filter_minmax(new_step, var_min[var], var_max[var])

                if debug and i == debug_step:
                    print('[OCR debug] New_step after filter_minmax')
                    print(new_step)

                # modify the array!
                mask = (new_var.coords['time'] == xr.DataArray(ctime[i+1]))

                new_var = xr.where(mask, new_step, new_var)

                if debug and i == debug_step:
                    print('Replaced new step looks like: ')
                    print(new_var)

            # enhance contrast
            if hist_stretch[var]:
                # Histogram stretching only supported for spatially averaged
                # data at the moment
                in_mean = new_var.mean('time')
                in_std = new_var.std('time')
                contrast_lims = [in_mean-in_std, in_mean+in_std]
                c_range = contrast_lims[1]-contrast_lims[0]

                if var_min[var] is not None:
                    out_min = var_min[var]
                else:
                    out_min = min0[var].item()
                if var_max[var] is not None:
                    out_max = var_max[var]
                else:
                    out_max = max0[var]

                if debug:
                    print('\n[OCR debug] Starting hist_stretch. Min/max: ')
                    print(out_min, out_max)

                hist_var = (c_dist[var](
                    ((new_var - contrast_lims[0]) /
                     c_range), **hist_args[var]) * (out_max - out_min) + out_min)

                cond = xr.DataArray((new_var > contrast_lims[0]) &
                                    (new_var < contrast_lims[1]), dims=['time'])
                new_var = xr.where(cond, hist_var, new_var)
                new_var = filter_minmax(new_var, var_min[var], var_max[var])

            if combine_steps > 1:
                new_var = new_var.resample(time=fby).interpolate()

            new_var.name = var
            out_vars.append(new_var)

        # print(out_vars)
        new_ds = xr.merge(out_vars)
        new_ds.attrs = data1.attrs
        if plot:
            plot_sa(new_ds, 'new')
        self.new = new_ds


def filter_minmax(dataset, var_min, var_max):

    if var_min is not None:
        # mina = np.ones(len(dataset['time'])) * var_min
        dataset = xr.where(
            dataset < var_min, var_min,
            dataset)

    if var_max is not None:
        # maxa = np.ones(len(dataset['time'])) * var_max
        dataset = xr.where(
            dataset > var_max, var_max,
            dataset)

    return(dataset)


def try_dict(kwargs, kw, except_val, main_vars):

    # If kwarg given, try to use it to fill params
    try:
        ndict = kwargs[kw]
        # If kwarg given, but its not a dictionary set it to the except_val
        if not isinstance(ndict, dict):
            except_val = ndict
            ndict = {}
    # Otherwise, set them all according to except_val
    except KeyError:
        ndict = {}

    # Check to see if all vars in ndict, otherwise fill with except_val
    for vi in main_vars:
        try:
            ndict[vi]
        except KeyError:
            ndict[vi] = except_val
    return(ndict)


def get_groupings(dt):
    if dt == 'monthly':
        by = 'time.month'
        fby = 'MS'
    elif dt == 'daily':
        by = 'time.dayofyear'
        fby = 'D'
    return(by, fby)
