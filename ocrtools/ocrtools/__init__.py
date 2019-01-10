#!/usr/bin/env python
# -*- coding: utf-8 -*-

########################################
#
#  OCRTOOLS - Making climate data malleable
#  Copyright (C) 2018 Andres Chang
#
########################################

from load import *

f2 = '../data/raw/b.e11.BRCP85C5CNBDRD.f09_g16.002.cice.h.aice_nh.200601-208012.nc'
f0 = '../data/raw/air.mon.mean.nc'

aice_mon = load(f2, 'aice')
# print(aice_mon.where((aice_mon.lat > 50) & (aice_mon.lat < 54)))

scope0 = scope()

aice_lim = piece(aice_mon, scope0)
print(aice_lim)
