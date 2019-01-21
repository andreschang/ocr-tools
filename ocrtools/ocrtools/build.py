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
from ocrtools.build_params import alabels, conversions, ylim


class build(object):

    def __init__(self, data1, var_lims, data2=None, combine_steps=1,
                 head=1, tail=1, savgol_window=0):
        """
        Makes new climate data based on existing modeled data and user options
        using a step-wise approach (i.e. OCR Tools calculates each timestep in
        sequence)
        """

        def resample_each(data, combine_steps):
            if self.dt == 'monthly':
                return (
                    data.resample(time=str(combine_steps) + self.fby, keep_attrs=True)
                        .mean('time'))
            elif self.dt == 'daily':
                grouped_data = data.groupby('time.year')
                return xr.concat(
                    [n.resample(time=str(combine_steps) + self.fby)
                      .mean('time')
                      .isel(time=slice(0, int(365/combine_steps)))
                     for m, n in grouped_data], 'time')

        def apply_savgol(dataset):
            for var in self.vars:
                if(savgol_window < 3):
                    pass
                else:
                    dataset[var] = xr.DataArray(signal.savgol_filter(
                        dataset[var], savgol_window[var], 2,
                        axis=dataset[var].get_axis_num('time')),
                        coords=dataset[var].coords, dims=dataset[var].dims)
            return(dataset)

        # Set parameters ======================================================
        self.vars = list(data1.data_vars)
        dt0 = int(data1.time[1]-data1.time[0])
        self.dt = 'daily' if dt0 < 2.3e15 else 'monthly'
        self.nvars = len(self.vars)
        self.by, self.fby = get_groupings(self.dt)
        self.t_dim = self.by.split('.')[1]
        self.yr0 = np.amin(data1.time.to_index()).year
        self.yrf = np.amax(data1.time.to_index()).year
        self.var_lims = xr.Dataset(
            {k: ('bound', var_lims[k]) for k in self.vars},
            coords={'bound': ['min', 'max']})

        self.v1 = resample_each(data1, combine_steps)
        self.v1 = apply_savgol(self.v1)
        s1 = self.v1.roll({'time': -1}, roll_coords=False) - self.v1
        self.V1 = self.v1.groupby(self.by).mean('time')
        self.S1 = s1.groupby(self.by).mean('time')
        self.DV1 = self.v1.groupby(self.by).std('time')
        self.DS1 = s1.groupby(self.by).std('time')
        self.max0, self.min0 = np.amax(self.v1), np.amin(self.v1)

        # Do the same for data2 if given, otherwise refer always to data1
        if data2 is not None:
            self.v2 = resample_each(data2, combine_steps)
            self.v2 = apply_savgol(self.v2)
            s2 = self.v2.roll({'time': -1}, roll_coords=False) - self.v2
            self.V2 = self.v2.groupby(self.by).mean('time')
            # self.S2 = s2.groupby(self.by).mean('time')
            self.DV2 = self.v2.groupby(self.by).std('time')
            self.DS2 = s2.groupby(self.by).std('time')

            self.max0 = np.amax(self.v2) if np.amax(self.v2) > self.max0 else self.max0
            self.min0 = np.amin(self.v2) if np.amin(self.v2) < self.min0 else self.min0
            v2_tail = (self.v2.sel(
                time=slice("{:04d}".format(self.yrf - tail + 1) + '-01-01',
                           "{:04d}".format(self.yrf) + '-12-31'))
                       .groupby(self.by).mean('time'))
            v1_head = (self.v2.sel(
                time=slice("{:04d}".format(self.yr0) + '-01-01',
                           "{:04d}".format(self.yr0 + head - 1) + '-12-31'))
                       .groupby(self.by).mean('time'))
            self.F = v2_tail.roll({self.t_dim: -1}, roll_coords=False) - v1_head

        else:
            self.V2, self.S2 = None, None
            self.DV2, self.DS2 = None, None
            self.F = None

    def mix(self, a):

        def scenario_merge(data1, data2):
            if data2 is None:
                return data1
            else:
                return (1 - a) * data1 + a * data2

        self.V = scenario_merge(self.V1, self.V2)
        self.OS = scenario_merge(self.S1, self.F)
        self.DV = scenario_merge(self.DV1, self.DV2)
        self.DS = scenario_merge(self.DS1, self.DS2)

    def regression(self, y_var, x_vars):
        from sklearn.linear_model import LinearRegression
        from sklearn.pipeline import make_pipeline, make_union
        from sklearn_xarray import Stacker, Select
        from scipy.stats import pearsonr

        vn = [self.v1]
        if self.V2 is not None:
            vn.append(self.v1)

        all_c, all_int, all_corr = [], [], []
        for vi in vn:
            x = vi[x_vars]
            y = vi[y_var]

            x_in = make_union(
                    *[make_pipeline(Select(xi), Stacker())
                      for xi in x_vars])
            mod = make_pipeline(x_in, LinearRegression())
            y_np = Stacker().fit_transform(y)
            mod.fit(x, y_np)
            lm = mod.named_steps['linearregression']
            coefs = list(lm.coef_.flat)
            intercept = lm.intercept_[0]
            y_pred = intercept + sum(
                [x[x_vars[i]] * coefs[i] for i in range(len(x_vars))])
            corr = pearsonr(y.values.ravel(), y_pred.values.ravel())[0]
            all_c.append(coefs), all_int.append(intercept)
            all_corr.append(corr)
        c = xr.DataArray(
            np.mean(np.array(all_c), axis=0), coords={
                'variable': x_vars}, dims='variable')
        corr, intercept = np.mean(all_corr), np.mean(all_int)
        return(c, corr, intercept)


    def new(self, y_vars=[], x_vars=[], a=1):

        self.mix(a)
        if x_vars == []:
            x_vars = self.vars
        else:
            pass

        # Get regression variables
        if len(y_vars) > 0:
            c, corr, intercept = [], [], []
            for yi in y_vars:
                ci, corr_i, intercept_i = self.regression(yi, x_vars)
                c.append(ci), corr.append(corr_i), intercept.append(intercept_i)

            # Matrices of regression coeffs, corrs, and intercepts
            c = xr.concat(c, 'y_var').assign_coords(y_var=y_vars)
            intercept = xr.DataArray(
                intercept, dims=['y_var'], coords={'y_var': y_vars})
            corr = xr.DataArray(corr, dims=['variable'], coords={'variable': y_vars})

        # Build new data for x_vars (i.e. vars that do not require regression)
        # Ux = self.v1[x_vars].to_array('x_var')
        # Ux = self.v1[x_vars]
        Ux = self.v1
        ctime = Ux.time.to_index()

        for i in range(len(ctime) - 1):

            t_sel = {self.t_dim: getattr(ctime[i], self.t_dim)}

            # Calculate adjusted center of normal distribution for next step
            norm_center = (((self.V.sel(t_sel) - Ux.isel(time=i)) /
                           self.DV.sel(t_sel))).to_array()
            norm_center = xr.where(np.isfinite(norm_center), norm_center, 0)
            r_norm = xr.DataArray(
                np.random.normal(norm_center),
                dims=['variable'], coords={'variable': self.vars})

            # Add a step amount based on random draw from normal distribution
            # multiplied by standard deviation added to optimum step
            new_step = (Ux.isel(time=i) + self.OS.sel(t_sel) +
                        (self.DS.sel(t_sel) * r_norm.to_dataset('variable')))

            # Calculate y_vars using regression and update random var
            # assignment proportionally to the correlation var
            y_pred = (((new_step[x_vars].to_array() * c).sum('variable') + intercept))
            new_step = new_step.update(
                y_pred.to_dataset('y_var') * corr +
                new_step[y_vars] * (1-corr)).squeeze()

            # Replace values that exceed min/max with min/max values
            new_step = xr.where(new_step < self.var_lims.sel(bound='min'),
                                self.var_lims.sel(bound='min'), new_step)
            new_step = xr.where(new_step > self.var_lims.sel(bound='max'),
                                self.var_lims.sel(bound='max'), new_step)
            Ux = xr.where(Ux.time == xr.DataArray(ctime[i + 1]), new_step, Ux)

        return(Ux)


def get_groupings(dt):
    if dt == 'monthly':
        by = 'time.month'
        fby = 'MS'
    elif dt == 'daily':
        by = 'time.dayofyear'
        fby = 'D'
    return(by, fby)
