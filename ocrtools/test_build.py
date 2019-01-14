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

TS = ocr.load_cesmLE('TS', dt='monthly', yr0=1990, yrf=2000, mem=2)
# print(TS['TS'])
TS_mean = ocr.annual_cycle(TS, dt='monthly', combine_steps=2)['TS']
step_var, blur_var = ocr.get_variance(TS, dt='monthly', combine_steps=1)

# TS_mean4 = ocr.annual_cycle(TS, dt='monthly', combine_steps=3)['TS']
# TS_mean3.mean(dim=['lat', 'lon']).plot()
# TS_mean4.mean(dim=['lat', 'lon']).plot()
step_var['TS'].mean(dim=['lat', 'lon']).plot()
# blur_var['TS'].mean(dim=['lat', 'lon']).plot()
plt.show()