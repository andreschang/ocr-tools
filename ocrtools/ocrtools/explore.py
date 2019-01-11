#!/usr/bin/env python
# -*- coding: utf-8 -*-

########################################
#
#  OCR TOOLS - Making climate data malleable
#  Copyright (C) 2018 Andres Chang
#
########################################

import numpy as np


def spatial_average(dataset):
    zonal_mean = dataset.mean(dim='lon')
    print(zonal_mean)


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

