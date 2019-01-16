#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#
#  OCR TOOLS - Making climate data malleable
#  Copyright (C) 2018 Andres Chang
#
###############################################################################

import ocrtools as ocr
from ocrtools import plt, mps_2_cmpday, K_2_F, p2
import xarray as xr
from ocrviz import colormap
import matplotlib as mpl

steps = [2]

# Project setting
project = 'EWYHF'
location = 'Wynnewood, PA, USA'
yr0, yrf = 1958, 1981
which_sims = [1850, 2]
all_vars = ['TS', 'RAIN']

# Plotting settings
alabels = {'TS': 'Temperature (F)', 'PRECT': 'Precipitation (cm/day)',
           'RAIN': 'Precipitation (cm/day)'}
conversions = {'PRECT': mps_2_cmpday, 'TS': K_2_F, 'RAIN': mps_2_cmpday}
var_ranges = {'PRECT': [0., None], 'TS': [None, None], 'RAIN': [0., None]}
dt = 'monthly'
fig_dpi = 200
show_snap = True

# Build settings
n_build = 5
step_std_a = {'PRECT': 8, 'TS': 13, 'RAIN': 8}
blur_std_a = {'PRECT': 8, 'TS': 22, 'RAIN': 8}
combine_steps = 1
a_range = [1, 2.5]
snap = {'PRECT': 11.5, 'TS': 9.5, 'RAIN': 8}
snap_curve = {'PRECT': 8, 'TS': 10, 'RAIN': 8}
savgol_window = {'PRECT': 5, 'TS': 3, 'RAIN': 5}
head = 3
tail = 3
combine_plots = False
hist_stretch = {'PRECT': True, 'TS': False, 'RAIN': True}
hist_dist = {'PRECT': p2, 'TS': False, 'RAIN': p2}

# Style settings
colors = ["#94B4DE", "#1E44D9", "#C71303", "#710F3F"]
b0 = colormap()
build_colors = b0.make_grad(
    ["#94B4DE", "#1E44D9", "#710F3F", "#C71303"], nsteps=2*(n_build+1))
style = 'pf_test'
ylim = {'PRECT': [0, 20], 'TS': [0, 100], 'RAIN': [0, 20]}

# Initialize a naming dict
fname_dict = {'yr_range': '-'.join(["{:04d}".format(x) for x in [yr0, yrf]]),
              'dt': dt}
which_sims = ["{:04d}".format(x) for x in which_sims]

# Make reformatted, combined datasets in project subfolder
if 1 in steps:
    scope = ocr.scope(location=location, interactive=False)

    for si in which_sims:
        merge_datasets = []
        for vi in all_vars:
            data0 = ocr.subset(ocr.load_cesmLE(vi, dt, yr0, yrf, int(si)), scope)
            ocr.save_reformatted(data0, dt, directory=project,
                                 post=si)

# Make new data
if 2 in steps:
    # Get reformatted data
    build_sims = {which_sims[0]: {},
                  which_sims[1]: {}}

    for si in which_sims:
        for var in all_vars:
            fname0 = fname_dict
            fname0['all_vars'] = var
            fname0['ending'] = 'nc'
            fname0['post'] = si

            fpath = ocr.gen_path(ocr.formatted_fname, fname0, '.', top=project)
            build_sims[si][var] = ocr.load(fpath, average=False)

    # For now, build applies the same parameters to each var in a dataset
    # so run the build script separately for each variable
    for ni in range(n_build):
        a = a_range[0] + (a_range[1]-a_range[0]) * ni / (n_build-1)

        for var in all_vars:
            if var == all_vars[0]:
                load_rand = None

            new = ocr.build(
                build_sims[which_sims[0]][var], dt, build_sims[which_sims[1]][var],
                combine_steps=combine_steps, step_std_a=step_std_a[var],
                blur_std_a=blur_std_a[var], a=a, var_min=var_ranges[var][0],
                var_max=var_ranges[var][1], snap=snap[var],
                snap_atten=snap_curve[var], head=head, tail=tail,
                hist_stretch=hist_stretch[var], hist_dist=hist_dist[var],
                load_rand=load_rand, debug=True, plot=True)

            fname0 = fname_dict
            fname0['all_vars'] = var
            fname0['post'] = "{:02d}".format(ni)
            fname0['ending'] = 'png'
            plt.savefig(ocr.gen_path(ocr.formatted_fname, fname0, '.',
                                     top=project), dpi=fig_dpi)
            plt.clf()

            if var == all_vars[0]:
                load_rand = new.rand


