#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#  
#  OCRTOOLS - Making climate data malleable
#  Copyright (C) 2018 Andres Chang
#
###############################################################################

import ocrtools.stage as st
import ocrtools.query as query

"""
File directory structure and naming conventions are specified
by passing a stage object to the query class.
"""
demo = st(preset = 'demo', time_as = 'date')
query.stage = demo

"""
ocrtools.query is used to explore climate data, set parameters,
and transform the data.

First, call a query object with file (or src)
Not sure what is in the file? Call first_look
"""
print("\nMYSTERY FILE ANALYSIS")
mystery_file = query(file = 'data/raw/air.mon.mean.nc' )
mystery_file.first_look()


"""
Ok, so we can tell that there are 4 variables: lat, lon, time, and air
Let's find out more about the 'air' var.
Uncomment the hashed lines to keep going!
"""

# mystery_file.first_look(var = 'air')

"""
Still quite a few unknowns but we keep on moving. 
The set_params() function is interactive and provides us with all kinds
of helpful feedback. I wonder if we can learn about air temperature in NYC...
Let's make a new pointer to the query instance called nyc_temp
"""

# nyc_temp = mystery_file
# nyc_temp.set_params()


"""
set_params takes a look at the data, autofills the parameters that it can figure out, 
and offers suggestions for some of the trickier ones like data_yr0.
Now we are ready to look at temperature over NYC. Uhh where is NYC?
Maybe someday we will be able to select the analysis area on a map, but for now...
That is a question for Google. Fill in the lat and lon ranges when you've got it figured out
"""

# nyc_temp.spatial_average(lat_bounds = [40,42], lon_bounds = [-75, -73])

"""
FYI, we can also set all of the parameters before submitting the script to save time 
once we know what we're doing...

nyc_temp.set_params(var = 'air', dt = 'monthly', lat_name = 'lat', lon_name = 'lon', \
dim = ["time", "lat", "lon"], data_yr0 = 1948)

Next, we will look at some CESM data. CESM is a world-class climate model developed by NCAR.
These data are from a series of runs called CESM-LE that simulates Earth's climate from
around 1920-2100. One version uses a "business as usual scenario", 
and the other is a "control run" with no humans.

ocr-tools is well-tuned to CESM data, so a lot of the parameters are filled automatically.
In fact, if you know the variable you're looking for, all you need to enter is the year-range.
"""

# print("\nCESM ANALYSIS")
# cesm_aice = query(src = 'cesm')
# cesm_aice.set_params(yr0 = 1980, yrf = 2050, mem = 2, dt = 'monthly', var = 'aice', hemisphere = 'nh')
# cesm_aice.spatial_average(lat_bounds = [80, 89], lon_bounds = [-179, 179])
