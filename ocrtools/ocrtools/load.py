#!/usr/bin/env python
# -*- coding: utf-8 -*-

########################################
#
#  OCRTOOLS - Making climate data malleable
#  Copyright (C) 2018 Andres Chang
#
########################################

import numpy as np
import pandas as pd
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

            if(ci == 'lat'):
                dataset = dataset.assign_coords(
                    lat=dataset.coords[c_key])
            elif(ci == 'lon'):
                dataset = dataset.assign_coords(
                    lon=dataset.coords[c_key])

    return(dataset)


def load(data, var=[], interactive=True, **kwargs):
    """
    Return a dataset with standardized coordinate names for easy querying

    Args:
    * data (str or dataset): either an xarray dataset or a path to netCDF that
    will be read in by xr_load function
    * var (string): name of main variable. If '?', will be set interactively
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

    ds = standardize_coords(ds)

    return(ds)
