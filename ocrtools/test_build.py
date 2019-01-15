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
fig_dpi = 200
TS2 = ocr.load_cesmLE('TS', dt, yr0=2008, yrf=2020, mem=2)
TS5 = ocr.load_cesmLE('TS', dt, yr0=2008, yrf=2020, mem=5)
scope = ocr.scope()
TS2 = ocr.spatial_average(ocr.subset(TS2, scope))
TS5 = ocr.spatial_average(ocr.subset(TS5, scope))

for i in range(5):
    if i == 0:
        debug = True
    else:
        debug = False
    TS_build = ocr.build(TS5, dt, TS2, combine_steps=2, plot=True, snap=0, debug=debug).new
    # TS2['TS'].mean(['lat', 'lon']).plot(label='mean_1')
    # TS5['TS'].mean(['lat', 'lon']).plot(label='mean_2')
    # ocr.spatial_average(TS_build['TS']).plot()
    # TS_build['TS'].plot()
    # plt.show()
    plt.savefig('TS_base_test_peru'+str(i)+'.png', dpi=fig_dpi)
    plt.clf()

# TS_mean = ocr.spatial_average(TS_build)

# figb, axb = plt.subplots(nrows=1, ncols=1, figsize=(9, 3.5))
# TS_mean['TS'].plot()

# plt.savefig('TS_build_test_peru.png', dpi=fig_dpi)