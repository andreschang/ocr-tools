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
TS2 = ocr.load_cesmLE('TS', dt, mem=2, yr0=2000, yrf=2004)
TS5 = ocr.load_cesmLE('TS', dt, mem=5, yr0=2000, yrf=2004)
scope = ocr.scope(location='Lima, Peru')
TS2 = ocr.spatial_average(ocr.subset(TS2, scope))
TS5 = ocr.spatial_average(ocr.subset(TS5, scope))
print(TS2)

# ocr.save_reformatted(TS2, dt=dt, pre='test', id='002')

# for i in range(5):
#     for j in range(4):
#         if i == 4 and j == 1:
#             debug = True
#         else:
#             debug = False
#         snap = (i) * 3
#         cs = j + 1

#         TS_build = ocr.build(
#             TS5, dt, TS2, combine_steps=cs, plot=True, snap=snap, debug=debug, debug_step=15,
#             head=2, tail=2).new
#         plt.savefig('TS_base_test_peru.snap'+str(snap)+'.cs'+str(cs)+'.png', dpi=fig_dpi)
#         plt.clf()
