#!/usr/bin/env python
# -*- coding: utf-8 -*-

########################################
#
#  OCR TOOLS - Making climate data malleable
#  Copyright (C) 2018 Andres Chang
#
########################################

import numpy as np
import xarray as xr
from ocrtools.tk_selector import plt
import cf_units
from datetime import date, datetime

def spatial_average(data, cell_area=None):
    """
    Returns approximate spatial average of 3D data

    Args:
    * data (xarray dataset or data_array)
    * cell_area (data_array): array of cell areas
    """

    try:
        d0 = data.to_array()
        dataset_in = True
    except AttributeError:
        d0 = data
        dataset_in = False

    if len(d0['lat'].shape) == 1:
        lat_weights = xr.DataArray(
            reg_wgt(np.min(d0['lat']), np.max(d0['lat']),
                    d0['lat'].shape[0]), coords=[('lat', d0['lat'])])

        # If dataset contains NaN values, calculate average of cell value,
        # weighted by latitude
        if np.sum(np.isnan(d0)) > 0:
            cell_weights = lat_weights * d0['lat'].shape[0] / np.sum(
                np.isfinite(d0.isel(variable=0, time=0)))
            weighted_mean = (d0 * cell_weights).sum(dim=['lat', 'lon'])

        # Otherwise, zonal mean and then weighted meridional mean (faster)
        else:
            zonal_mean = d0.mean(dim='lon')
            weighted_mean = zonal_mean.dot(lat_weights)

    elif len(d0['lat'].shape) == 2:
        if cell_area is None:
            raise ValueError('Cell area must be defined as a function arg',
                             ' to calculate spatial average of a variable',
                             ' across a curvilinear grid')
        else:
            weighted_mean = (
                d0 * cell_area/cell_area.sum()).sum(
                dim=[n for n in d0.dims if n != 'time' and n != 'variable'])

    if dataset_in:
        weighted_mean = weighted_mean.to_dataset(dim='variable')

    return(weighted_mean)


def reg_wgt(latmin, latmax, nlat):
    """
    Returns a list of weighted values (sum = 1) based on
    latitudinal range for calculating spatial averages.
    Assumes that lat bands are evenly spaced (ex. 10N, 20N, 30N...)

    Args:
    * latmin (int): Minimum lat of query area
    * latmax (int): Maximum lat of query area
    * nlat (int): Number of latitudinal bands
    """
    if latmin == latmax:
        return([1])
    else:
        dy = (latmax-latmin)/(nlat-1)
        lats = np.arange(latmin, latmax+dy, dy)[:nlat]
        wgts = np.zeros(nlat)
        if latmin == -90:
            for i in range(nlat):
                if ((i != 0) & (i != nlat-1)):
                    wgts[i] = np.abs(np.sin(np.deg2rad(lats[i]+(dy/2))) -
                                     np.sin(np.deg2rad(lats[i]-(dy/2))))
                else:
                    wgts[i] = 1-np.abs(np.sin(np.deg2rad(lats[i]+(dy/2))))
        else:
            for i in range(nlat):
                wgts[i] = np.abs(np.sin(np.deg2rad(lats[i]+(dy/2))) -
                                 np.sin(np.deg2rad(lats[i]-(dy/2))))

        wgts = wgts/(np.sum(wgts))
        return(wgts)

