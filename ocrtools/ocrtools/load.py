#!/usr/bin/env python
# -*- coding: utf-8 -*-

########################################
#
#  OCRTOOLS - Making climate data malleable
#  Copyright (C) 2018 Andres Chang
#
########################################

import numpy as np
import xarray as xr


def xr_load(path):
    """
    Loads a netCDF and fixes date convention for CESM-LE or other Gregorian
    calendar data

    Args:
    * path (str): path to data
    """

    try:
        d = xr.open_dataset(path)
    except ValueError:
        import cf_units
        d0 = xr.open_dataset(path, decode_times=False)
        d = d0.assign_coords(time=cf_units.num2date(d0.time-1, d0.time.units,
                             d0.time.calendar))
    return(d)


def standardize_coords(dataset, center=True):
    """
    Standardize coordinates of an xarray dataset (i.e. make sure that there
    are coordinates named lat, lon, and time, which can be specified from
    other namings)

    Args:
    * dataset (xarray.Dataset)
    * center (bool): True by default - uses TLAT/TLON if multiple coordinate
    systems. If False, uses ULAT/ULON (case insensitive)
    """

    for ci in ['lat', 'lon']:
        if ci not in dataset.coords:
            c_options = [x for x in dataset.coords if ci in x.lower()]

            if len(c_options) == 1:
                c_key = c_options[0]
            elif len(c_options) > 1:
                try_tc = [x for x in c_options if 't' in
                          x.lower().replace('lat', 'la')]
                try_uc = [x for x in c_options if 'u' in str.lower(x)]
                if center and len(try_tc) == 1:
                    c_key = try_tc[0]
                elif center is False and len(try_uc) == 1:
                    c_key = try_uc[0]

            try:
                c_key
            except NameError:
                print(dataset.coords)
                c_key = input('Please enter ' + ci +
                              ', as shown in the original coordinate list: ')

            dataset.coords[ci] = dataset.coords[c_key]


    return(dataset)


def load(data, var=[], interactive=True, drop=True, **kwargs):
    """
    Return a dataset with standardized coordinate names for easy querying

    Args:
    * data (str or xarray.Dataset): either an xarray dataset or a path to
    netCDF that will be read in by xr_load function
    * var (string): name of main variable (optional)
    * interactive (bool): instructions for the less-experienced programmer :)
    * drop (bool): if True and var(s) specified, other vars are removed from
    the dataset
    """

    # Load xarray dataset
    if isinstance(data, str):
        ds = xr_load(data)
    elif isinstance(data, xr.core.dataset.Dataset):
        ds = data
    else:
        raise TypeError(
            "data argument must be a path to a netCDF or xarray dataset")

    # Add an attribute for main_vars
    if isinstance(var, str):
        var = [var]
    elif isinstance(var, list):
        var = var
    else:
        raise TypeError(
            "var argument must be a string or list indicating the main" +
            " variable, if any")

    # Specify in command line, if interactive set to True
    if interactive and var == []:
        print("Variables in dataset:\n", str(ds.data_vars))
        v0 = input('Please enter a main variable or list of main' +
                   ' variables separated by commas: ')
        if ',' in v0:
            var = str.strip(v0.split(','))
        elif isinstance(v0, str):
            var = [v0]

    if all(vi in ds for vi in var):
        ds.attrs['main_vars'] = var
    else:
        raise ValueError('main_vars must specify data variables in the xarray')

    # Standardize coordinate names
    try:
        ds = standardize_coords(ds, center=kwargs['center'])
    except KeyError:
        ds = standardize_coords(ds)

    if drop:
        ds = ds[ds.main_vars]

    return(ds)


def piece(dataset, scope):
    """
    Returns a subsetted dataset based on scope
    Args:
    * dataset (xarray.Dataset)
    * scope (ocrtools.scope)
    """

    scope_keys = [x for x in scope.__dict__.keys()]

    d_subset = dataset.sel(time=slice(scope.yr0+'-01-01', scope.yrf+'-12-31'))
    d_subset = d_subset.where(
        (dataset.lat >= scope.lat_min) &
        (dataset.lat <= scope.lat_max) &
        (dataset.lon >= scope.lon_min) &
        (dataset.lon <= scope.lon_max))

    if 'z_max' in scope_keys:
        d_subset = d_subset.where(d_subset.z <= scope.z_max)
    if 'z_min' in scope_keys:
        d_subset = d_subset.where(d_subset.z >= scope.z_min)

    return(d_subset)


class scope(object):
    """
    Class that is used to subset data in space and time

    Args:
    * interactive (bool): instructions for the less-experienced programmer :)
    * kwargs: lat_min, lat_max, lon_min, lon_max, yr0, yrf are always defined
    to some extent. If z_min and z_max are included as keywords, they will also
    be used added as attributes to the scope object
    """

    def __init__(self, interactive=True, **kwargs):

        scopes = ['yr0', 'yrf', 'lat_min', 'lat_max', 'lon_min', 'lon_max']
        none_vals = {'lat_min': -91, 'lat_max': 91, 'lon_min': -361,
                     'lon_max': 361, 'yr0': 1, 'yrf': 10000}

        try:
            if(kwargs['z_min']):
                scopes = scopes + ['z_min']
        except KeyError:
            pass
        try:
            if(kwargs['z_max']):
                scopes = scopes + ['z_max']
        except KeyError:
            pass

        for ai in scopes:
            try:
                setattr(self, ai, kwargs[ai])
            except KeyError:
                if interactive:
                    lim0 = input('Enter ' + ai + ': ')
                    if lim0 == '':
                        lim0 = none_vals[ai]

                    if 'yr' in ai:
                        setattr(self, ai, str(lim0))
                    else:
                        setattr(self, ai, float(lim0))

                else:
                    setattr(self, ai, none_vals[ai])
