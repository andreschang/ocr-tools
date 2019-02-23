#!/usr/bin/env python
# -*- coding: utf-8 -*-

########################################
#
#  OCRTOOLS - Making climate data malleable
#  Copyright (C) 2018 Andres Chang
#
########################################

import xarray as xr
from scipy import signal
import numpy as np


def mps_2_cmpday(data):
    t_day = 60.*60.*24.   # seconds * minutes * hours
    lconvert = data * t_day * 100
    return lconvert


def K_2_F(data):
    t_f = data * (9./5.) - 459.67
    return t_f


def mmH2O_2_inSNO(data):
    H2O_row = 997  # kg/m^3
    SNO_row = 100  # kg/m^3
    mm_2_inch = 0.0393701
    mmSNO = (mm_2_inch * H2O_row/SNO_row) * data
    return(mmSNO)


def mmps_2_cmday(data):
    return mps_2_cmpday(data/1000.)


def percentify(data): return 100 * data


def inverse(data): return 100 * (1-data)


# Conversion and plotting settings
ylabels = {'TS': 'Temperature (F)', 'PRECT': 'Precipitation (cm/day)',
           'RAIN': 'Precipitation (cm/day)', 'H2OSNO': 'Snow cover (in)',
           'TREFHT': 'Temperature (F)', 'RELHUM': 'Relative humidity (%)',
           'FSNSCLOUD': "Cloudiness (% SW energy blocked)"}
conversions = {'PRECT': mps_2_cmpday, 'TS': K_2_F, 'RAIN': mmps_2_cmday,
               'H2OSNO': mmH2O_2_inSNO, 'TREFHT': K_2_F, 'RELHUM': percentify,
               'FSNSCLOUD': inverse}
ylim = {'PRECT': [0, 20], 'TS': [0, 100], 'RAIN': [0, 20], 'H2OSNO': [0, 5],
        'TREFHT': [0, 100], 'FSNSCLOUD': [0, 80], 'RELHUM': [60, 100]}


# Build settings
var_lims = {'PRECT': [0., 1000.], 'TS': [0., 1000.], 'RAIN': [0., 1000.],
            'H2OSNO': [0., 1000.], 'FSNSCLOUD': [0., 100.],
            'TREFHT': [0., 1000.], 'RELHUM': [0., 100.]}


class build(object):

    def __init__(self, data1, data2=None, var_lims=var_lims, combine_steps=1,
                 head=1, tail=1, savgol_window=0, y_vars=[], x_vars=[],
                 debug=False):
        """
        Makes new climate data based on existing modeled data and user options
        using a step-wise approach (i.e. OCR Tools calculates each timestep in
        sequence)
        """

        def resample_each(data, combine_steps):
            if self.dt == 'monthly':
                return (
                    data.resample(time=str(combine_steps) + self.fby,
                                  keep_attrs=True)
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
        self.debug = debug
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
        self.nt = self.V1.dims[self.t_dim]
        self.nyrs = self.yrf-self.yr0+1
        self.v1_head = (self.v1.sel(
            time=slice("{:04d}".format(self.yr0) + '-01-01',
                       "{:04d}".format(self.yr0 + head - 1) + '-12-31'))
                   .groupby(self.by).mean('time'))

        # Do the same for data2 if given, otherwise refer always to data1
        if data2 is not None:
            self.v2 = resample_each(data2, combine_steps)
            self.v2 = apply_savgol(self.v2)
            s2 = self.v2.roll({'time': -1}, roll_coords=False) - self.v2
            self.V2 = self.v2.groupby(self.by).mean('time')
            self.S2 = s2.groupby(self.by).mean('time')
            self.DV2 = self.v2.groupby(self.by).std('time')
            self.DS2 = s2.groupby(self.by).std('time')

            self.max0 = (np.amax(self.v2) if np.amax(self.v2) > self.max0
                         else self.max0)
            self.min0 = (np.amin(self.v2) if np.amin(self.v2) < self.min0
                         else self.min0)
            self.v2_tail = (self.v2.sel(
                time=slice("{:04d}".format(self.yrf - tail + 1) + '-01-01',
                           "{:04d}".format(self.yrf) + '-12-31'))
                       .groupby(self.by).mean('time'))

        else:
            self.V2, self.S2 = None, None
            self.DV2, self.DS2 = None, None
            self.v2_tail = (self.v1.sel(
                time=slice("{:04d}".format(self.yrf - tail + 1) + '-01-01',
                           "{:04d}".format(self.yrf) + '-12-31'))
                       .groupby(self.by).mean('time'))

        # Get regression variables
        self.y_vars = y_vars
        if y_vars == []:
            self.x_vars = self.vars
        else:
            self.x_vars = x_vars

        if len(self.y_vars) > 0:
            c, corr, intercept = [], [], []
            for yi in y_vars:
                ci, corr_i, intercept_i = self.regression(yi, x_vars)
                c.append(ci), corr.append(corr_i)
                intercept.append(intercept_i)

            # Matrices of regression coeffs, corrs, and intercepts
            self.c = xr.concat(c, 'y_var').assign_coords(y_var=y_vars)
            self.intercept = xr.DataArray(
                intercept, dims=['y_var'], coords={'y_var': y_vars})
            self.corr = xr.DataArray(
                corr, dims=['y_var'], coords={'y_var': y_vars})

    def apply_var_lims(self, data):
        # Replace values that exceed min/max with min/max values
        data = xr.where(data < self.var_lims.sel(bound='min'),
                        self.var_lims.sel(bound='min'), data)
        data = xr.where(data > self.var_lims.sel(bound='max'),
                        self.var_lims.sel(bound='max'), data)
        return(data)

    def mix(self, a):

        def scenario_merge(data1, data2):
            if data2 is None:
                return data1
            else:
                am = 1 if a > 1 else a
                return (1 - am) * data1 + am * data2

        self.DV = scenario_merge(self.DV1, self.DV2)
        self.DS = scenario_merge(self.DS1, self.DS2)

        # Interpolation function
        def rs(group):
            g = group.resample(time='1YS').interpolate('linear')
            return(g)

        # Softmax <3 normalizes all data to sum to one
        def softmax(data, dim):
            e_data = np.exp(data - np.max(data))
            return e_data / e_data.sum(dim=dim)

        # weirdo function that ensures total sum of year steps, before F, = 0
        def normo(data):
            excess = data.sum()
            pos, neg = data.where(data >= 0), data.where(data < 0)
            data = xr.where(
                data >= 0, data - excess/2 * softmax(pos, 'time'),
                data - excess/2 * softmax(neg, 'time'))
            return(data)

        # Get annual step cycle of S1 head and S2 tail
        # Step: delta from e.g. Jan to Feb, with coordinate aligned to Jan
        SH1 = ((self.v1_head.roll({self.t_dim: -1}, roll_coords=False)
               - self.v1_head))
        ST2 = ((self.v2_tail.roll({self.t_dim: -1}, roll_coords=False)
               - self.v2_tail))
        ones = xr.DataArray(np.ones((self.nyrs, self.nt)),
                            dims=['year', self.t_dim])

        # Interpolate step cycle for all years between S1 and S2
        OS = (xr.concat([SH1, ST2], 'time')
                .assign_coords(time=self.v1.time[[0, -self.nt]])
                .groupby(self.t_dim).apply(rs)
                .rename({'time': 'year'})
                .stack(time=('year', self.t_dim))
                .assign_coords(time=self.v1.time))
        OS = OS.groupby('time.year').apply(normo)

        # Add an amount that sums to the total, mean delta between v2_tail and
        # v1_head weighted by step in year and distributed across n yrs evenly
        step_delt = (self.v2_tail - self.v1_head)
        mean_delt = self.v2_tail.mean() - self.v1_head.mean()
        F = (a * softmax(step_delt, self.t_dim) * mean_delt) / (self.nyrs)

        # Add F to OS, but shift backward by one month because this is a step
        # (i.e. Jan to Feb), not an amount added to each month
        self.OS = (OS + (
            F.roll({self.t_dim: -1}, roll_coords=False) * ones)
            .stack(time=('year', self.t_dim)).assign_coords(time=OS.time))

        # Make a scenario from v1 with all steps added. Just cumulative sum of
        # all opt_steps where time0 is equal to v1
        OS_V = xr.where(
            self.OS.time == np.amax(self.OS.time),
            self.v1.isel(time=0), self.OS).roll({'time': 1}, roll_coords=False)
        self.V = self.apply_var_lims(
            np.cumsum(OS_V, axis=0).assign_coords(time=self.v1.time))

        if self.debug:
            from ocrtools import plt
            for vi in self.vars:
                self.V[vi].plot()
                plt.show()
                plt.close()
                self.v1_head[vi].plot(label='v1 head')
                self.v2_tail[vi].plot(label='v2 tail')
                self.V1[vi].plot(label='V1')
                self.V2[vi].plot(label='V2')
                plt.legend()
                plt.show()
                plt.close()

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

    def new(self, a=1, unravel=0, debug=False):

        self.mix(a)
        Ux = self.V.copy().assign_coords(
            lat=self.V1.coords['lat'],
            lon=self.V1.coords['lon'])
        ctime = Ux.time.to_index()

        for i in range(len(ctime) - 1):

            t_sel = {self.t_dim: getattr(ctime[i], self.t_dim)}

            # Calculate adjusted center of normal distribution for next step
            norm_center = (((self.V.isel(time=i) - Ux.isel(time=i)) /
                           self.DV.sel(t_sel))).to_array() * (1 - unravel)
            norm_center = xr.where(np.isfinite(norm_center), norm_center, 0)

            r_norm = xr.DataArray(
                np.random.normal(norm_center),
                dims=list(norm_center.dims),
                coords=norm_center.coords)

            # Add a step amount based on random draw from normal distribution
            # multiplied by standard deviation added to optimum step
            new_step = (Ux.isel(time=i) + self.OS.isel(time=i) +
                        (self.DS.sel(t_sel) * r_norm.to_dataset('variable')))
            # Calculate y_vars using regression and update random var
            # assignment proportionally to the correlation var
            if self.y_vars != []:
                y_pred = (((new_step[self.x_vars].to_array() * self.c)
                          .sum('variable') + self.intercept))
                new_step = new_step.update(
                    (y_pred * self.corr +
                     new_step[self.y_vars].to_array('y_var') *
                     (1-self.corr)).to_dataset('y_var'))

            # Replace values that exceed min/max with min/max values
            new_step = self.apply_var_lims(new_step)
            Ux = xr.where(
                Ux.time == xr.DataArray(ctime[i + 1]),
                new_step, Ux)

        return(Ux)


def get_groupings(dt):
    if dt == 'monthly':
        by = 'time.month'
        fby = 'MS'
    elif dt == 'daily':
        by = 'time.dayofyear'
        fby = 'D'
    return(by, fby)
