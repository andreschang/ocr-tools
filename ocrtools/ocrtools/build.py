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
            a = kwargs["a"]
        except KeyError:
            a = 1
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
        try:
            debug = kwargs["debug"]
        except KeyError:
            debug = False
        try:
            debug_step = kwargs["debug_step"]
        except KeyError:
            debug_step = 0

        main_vars = data1.attrs['main_vars']
        nvars = len(main_vars)
        by, fby = get_groupings(dt)
        t_dim = by.split('.')[1]

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
                ax0.set_title('OCR Build '+main_vars[0])
            else:
                for vi in range(nvars):
                    sa[main_vars[vi]].plot(ax=ax0[vi])
                    ax0[vi].legend(loc='best')
                    ax0[vi].title('OCR Build '+main_vars[vi])


        # Get annual cycle and standard dev params of data1
        # calculate standard devs:
        # step_std is detrended standard dev of each div
        # (ex. january or mean(january, february))
        # blur_std is standard dev of each "combine_steps" segment
        # (ex. january (0) or std(january, february) year 0)
        d1 = (data1.resample(time=str(combine_steps)+fby, keep_attrs=True)
                   .mean('time').dropna('time', how='all'))
        d1_step_std, d1_blur_std = get_stds(data1, dt, **kwargs)
        if plot:
            plot_sa(d1, 1)

        # Do the same for data2 if given, otherwise refer always to data1
        if data2 is not None:
            d2 = (data2.resample(time=str(combine_steps)+fby, keep_attrs=True)
                       .mean('time').dropna('time', how='all'))
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

        d1_steps = d1_group.roll({t_dim: -1}, roll_coords=False) - d1_group
        full_steps = d_yrf.roll({t_dim: -1}, roll_coords=False) - d_yr0
        opt_steps = d1_steps + a * (full_steps - d1_steps)/(yrf - yr0 + 1)
        self.rand = []
        self.snap_list = []
        out_vars = []

        if debug:
            print('New loop', opt_steps)

        for vi in range(nvars):
            var_i = main_vars[vi]
            new_var = (
                data1[var_i].resample(time=str(combine_steps) + fby,
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
                new_step = ((new_var.isel(time=i).copy(deep=True) + opt_steps[var_i].sel(
                        {t_dim: getattr(ctime[i], t_attr)}))
                                   .assign_coords(time=xr.DataArray(ctime[i+1])))

                if debug and i == debug_step:
                    print('\n[OCR debug] Date of loop: '+str(ctime[i]) +
                          ', Date of replacement '+str(ctime[i+1]))
                    print('[OCR debug] New_step from opt_step sum')
                    print(new_step)
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
                    else:
                        r1, r2 = self.rand[i]
                self.rand.append((r1, r2))

                new_step = (new_step + r1 * add_step_std * step_std_a/50. +
                             r2 * add_blur_std * blur_std_a/50.)

                if debug and i == debug_step:
                    print('[OCR debug] New_step after the addition of randomness')
                    print(new_step)

                # Normalize snapping envelope to 5 standard deviations by default
                dev_from_base = (new_step - base_ac[var_i].sel(
                    {t_dim: getattr(ctime[i+1], t_attr)})).apply(
                        np.abs)/(add_step_std * (15-snap))
                dev_from_base = xr.where(dev_from_base > 1, 1, dev_from_base)
                snap_amt = dev_from_base**((50 - snap_atten)/5)
                self.snap_list.append(snap_amt)

                new_step = (
                    snap_amt * self.base_data[var_i].sel({'time': ctime[i+1]}) +
                    (1 - snap_amt) * new_step)

                if debug and i == debug_step:
                    print('[OCR debug] New_step after snapping')
                    print(new_step)

                # Replace variables that exceed min or max with min/max values
                new_step = filter_minmax(new_step, var_i, var_min, var_max)

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
            if hist_stretch:
                in_mean = new_var.mean('time')
                in_std = new_var.std('time')
                contrast_lims = [in_mean-in_std, in_mean+in_std]
                c_range = contrast_lims[1]-contrast_lims[0]

                if var_min is not None:
                    out_min = var_min
                else:
                    out_min = min0
                if var_max is not None:
                    out_max = var_max
                else:
                    out_max = max0

                # print('\n\n\nHIST CALC')
                hist_in = ((new_var - contrast_lims[0]) /
                           c_range)
                hist_var = (c_dist(
                    ((new_var - contrast_lims[0]) /
                     c_range)) * (out_max - out_min) + out_min)

                new_var = xr.where(
                    new_var > contrast_lims[0] and new_var < contrast_lims[1],
                    hist_var, new_var)

            if combine_steps > 1:
                new_var = new_var.resample(time=fby).interpolate()

            # new_var = new_var.to_array(name=var_i).squeeze()
            # print(new_var)
            
            out_vars.append(new_var)

        new_ds = xr.merge(out_vars)
        new_ds.attrs = data1.attrs
        if plot:
            plot_sa(new_ds, 'new')
        self.new = new_ds


def filter_minmax(dataset, var_i, var_min, var_max):

    if var_min is not None:
        if isinstance(var_min, dict):
            try:
                dataset = xr.where(
                    dataset < var_min[var_i], var_min[var_i],
                    dataset)
            except KeyError:
                pass
        else:
            dataset = xr.where(
                dataset < var_min, var_min,
                dataset)

    if var_max is not None:
        if isinstance(var_max, dict):
            try:
                dataset = xr.where(
                    dataset > var_max[var_i], var_max[var_i],
                    dataset)
            except KeyError:
                pass
        else:
            dataset = xr.where(
                dataset > var_max, var_max,
                dataset)

    return(dataset)


def get_groupings(dt):
    if dt == 'monthly':
        by = 'time.month'
        fby = 'MS'
    elif dt == 'daily':
        by = 'time.day'
        fby = 'D'
    return(by, fby)
