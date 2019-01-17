#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#
#  OCR TOOLS - Making climate data malleable
#  Copyright (C) 2018 Andres Chang
#
###############################################################################

import numpy as np
import ocrtools as ocr
from ocrtools import plt
import xarray as xr
from ocrviz import colormap
import matplotlib as mpl
from datetime import datetime

steps = [2]

# Project settings
project = 'EWYHF'
location = 'Wynnewood, PA, USA'
yr0, yrf = 1958, 1981
which_sims = [1850, 2]
all_vars = ['TS', 'H2OSNO']

# Build settings
n_build = 10
a_range = [1, 2.5]
dt = 'daily'
ocr.build_kwargs['combine_steps'] = 5
ocr.build_kwargs['plot'] = True
ocr.build_kwargs['debug'] = False
ocr.build_kwargs['head'] = 5
ocr.build_kwargs['tail'] = 5
ocr.build_kwargs['hist_args']['H2OSNO'] = {
    'steepness': 30, 'p': 0.025, 'knee': 0.96}

# File settings
fig_dpi = 200
show_snap = True
ddir = ('/Volumes/Samsung_T5/Open_Climate_Research-Projects/prj-EWYHF_Exact' +
        'ly_Where_You_Had_Fallen/data')
pdir = ('/Volumes/Samsung_T5/Open_Climate_Research-Projects/prj-EWYHF_Exact' +
        'ly_Where_You_Had_Fallen/scratch_plots')
now = datetime.now()
scratchId = now.strftime("%Y%m%d%H%M")

# Style settings
colors = ["#94B4DE", "#1E44D9", "#C71303", "#710F3F"]
build_colors = colormap().make_grad(
    ["#94B4DE", "#1E44D9", "#710F3F", "#C71303"], nsteps=2*(n_build+1))
style = 'pf_test'

# Initialize filenaming variables
fname_dict = {'yr_range': '-'.join(["{:04d}".format(x) for x in [yr0, yrf]]),
              'dt': dt,
              'all_vars': '_'.join(all_vars)}
which_sims = ["{:04d}".format(x) for x in which_sims]

# Make reformatted, combined datasets in project subfolder
if 1 in steps:
    scope = ocr.scope(location=location, interactive=False)

    for si in which_sims:
        merge_datasets = []
        for vi in all_vars:
            data0 = ocr.subset(
                ocr.load_cesmLE(vi, dt, yr0, yrf, int(si)), scope)
            data0 = data0.assign_coords(lat=scope.lat, lon=scope.lon)
            merge_datasets.append(data0)
            # ocr.save_reformatted(data0, dt, directory=project,
            #                      post=si)
        dataf = xr.merge(merge_datasets)
        dataf.attrs['main_vars'] = all_vars
        ocr.save_reformatted(dataf, dt, directory=ddir, post=si)

# Make new data
if 2 in steps:
    # Set style
    mpl.style.use(style)
    mpl.rcParams['lines.linewidth'] = 0.5

    # Generate reformatted data filenames and add to two-item list
    build_data = []
    fname0 = fname_dict
    fname0['ending'] = 'nc'

    for si in which_sims:
        fname0['post'] = si
        build_data.append(
            ocr.gen_path(ocr.formatted_fname, fname0, '.', top=ddir))
    # Load data
    build_data = [ocr.load(x, all_vars, interactive=False, average=False)
                  for x in build_data]

    # Build n times with function args
    for ni in range(n_build):
        try:
            a = a_range[0] + (a_range[1]-a_range[0]) * ni / (n_build-1)
        except ZeroDivisionError:
            a = a_range[0]

        ocr.build_kwargs['a'] = a
        ocr.build_kwargs['step_std_a']['H2OSNO'] = 20 + 3 * ni
        ocr.build_kwargs['blur_std_a']['H2OSNO'] = 20 + 3 * ni

        new = ocr.build(build_data[0], dt, build_data[1], **ocr.build_kwargs)
        # new = ocr.build(build_data[0], dt, **ocr.build_kwargs)

        # save figure
        fname0['post'] = scratchId + "-" + "{:02d}".format(ni)
        fname0['ending'] = 'png'
        outf = ocr.gen_path(ocr.formatted_fname, fname0, '.', top=pdir)
        plt.savefig(outf, dpi=fig_dpi)
        plt.clf()

        # save settings
        fname_txt = '.'.join(outf.split('.')[0:-1])+'.txt'
        report_f = str(ocr.build_kwargs)
        np.savetxt(fname_txt, [report_f], fmt='%s')

