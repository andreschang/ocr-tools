#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#
#  OCR TOOLS - Making climate data malleable
#  Copyright (C) 2018 Andres Chang
#
###############################################################################

"""
OCR Tools is an interactive program, so you can visualize and transform climate
data with just a few basic functions and minimal syntax

PART 1: LOADING DATA
Climate data are most commonly stored in netCDF files -- a great format, but
one that is recognized by almost no one outside the geosciences community
"""

import ocrtools as ocr

f_temp_obs = 'data/raw/air.mon.mean.nc'
temp_obs = ocr.load(f_temp_obs)

# ocr.load() initializes a dataset from a netCDF. Some netCDFs contain multiple
# variables; in these cases, the variable(s) of interests needs to be specified
f_temp_simulated = 'data/raw/b.e11.BRCP85C5CNBDRD.f09_g16.002.cam.h0.TS.200601-208012.nc'
# temp_simulated = ocr.load(f_temp_simulated)

# Instead of specifying the main variable interactively, you can also include
# it as a function argument, var. The result is identical
temp_simulated = ocr.load(f_temp_simulated, var='TS')


"""
PART 2: SUBSETTING THE DATA
One challenge of working with climate data is that the files can be huge!
Monthly data, like we have here, isn't so bad, but daily or hourly data can
be a few gigabytes
"""

# In OCR TOOLS, 'scope' is an object that represents the geographic bounds
# and/or year range of inquiry. It can be used to subset any dataset
# map_selection = ocr.scope(
#     yr0=2008, yrf=2010, lat_min=0, lat_max=60, lon_min=-8.5, lon_max=-0.5)
# subsetted_temp_obs = ocr.subset(temp_obs, map_selection)
# subsetted_temp_sim = ocr.subset(temp_simulated, map_selection)


"""
PART 3: SOME BASIC ANALYSIS
Many basic transformations can be applied to xarray datasets
Let's take a look at the average temperature across areas that we selected
"""

# TESTING 1D WEIGHTED MEAN.. IT WORKS :)

# gwtest = ocr.subset(ocr.load(f_temp_simulated, var=['TS', 'gw']), map_selection).to_array()
# gw = gwtest.sel(variable='gw').isel(lon=0, time=0)
# gw = gw/gw.sum()
# mean = (gwtest.sel(variable='TS')
#               .mean(dim='lon')
#               .dot(gw))

# print(ocr.spatial_average(subsetted_temp_sim)['TS'])
# print(mean)

# print(mean)
# print(average_temp_sim)
# print(average_temp_obs)

# TESTING 2D
faice = 'data/raw/b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.aice_nh.192001-200512.nc'
aice_data = ocr.load(faice, var=['aice', 'tarea'])
map_selection = ocr.scope(
    yr0=2008, yrf=2010, lat_min=70, lat_max=90, lon_min=-20, lon_max=20)
subsetted_temp_obs = ocr.subset(aice_data, map_selection)
subsetted_temp_sim = ocr.subset(aice_data, map_selection)

ocr.spatial_average(aice_data, aice_data['tarea'])

# APPENDIX: SETTING THE SCOPE MANUALLY
# iberian_peninsula_2008_2010 = ocr.scope(
#     yr0=2008, yrf=2010, lat_min=37.5, lat_max=43, lon_min=-8.5, lon_max=-0.5)
# iberian_temp_obs = ocr.subset(temp_obs, iberian_peninsula_2008_2010)
# iberian_temp_sim = ocr.subset(temp_simulated, iberian_peninsula_2008_2010)
