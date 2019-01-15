#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#
#  OCR TOOLS - Making climate data malleable
#  Copyright (C) 2018 Andres Chang
#
###############################################################################

import ocrtools as ocr
from ocrtools import plt
import xarray as xr

dt = 'monthly'
TS2 = ocr.load_cesmLE('TS', dt, yr0=2030, yrf=2050, mem=2)
TS5 = ocr.load_cesmLE('TS', dt, yr0=2030, yrf=2050, mem=5)

TS_build = ocr.build(TS2, dt, data2=TS5, combine_steps=2, plot=True).new

TS_mean = ocr.spatial_average(TS_build)

# fig, ax = plt.subplots(nrows=1, ncols=1)

# plt.show()
plt.clf()

ocr.spatial_average(TS2)['TS'].plot(label='02')
ocr.spatial_average(TS5)['TS'].plot(label='05')
TS_mean['TS'].plot(label='new')
plt.legend(loc='best')
plt.show()