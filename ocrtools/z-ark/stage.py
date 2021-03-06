#!/usr/bin/env python
# -*- coding: utf-8 -*-

########################################
#
#  OCRTOOLS - Making climate data malleable
#  Copyright (C) 2018 Andres Chang
#
########################################

import os
import os.path
from datetime import datetime


class stage(object):

    ndivs = {'daily': 365, 'monthly': 12}

    def rcsv(fname):
        import csv
        with open(fname) as f:
            reader = csv.reader(f)
            next(reader)  # skip header
            data = [r for r in reader]

        return data

    cice = rcsv(os.path.join(os.path.dirname(__file__),
                'var_lists/cice_vars.csv'))
    cice_vars = [itm[3] for itm in cice]
    cam = rcsv(os.path.join(os.path.dirname(__file__),
               'var_lists/cam_vars.csv'))
    cam_vars = [itm[3] for itm in cam]
    clm = rcsv(os.path.join(os.path.dirname(__file__),
               'var_lists/clm_vars.csv'))
    clm_vars = [itm[3] for itm in clm]
    pop = rcsv(os.path.join(os.path.dirname(__file__),
               'var_lists/pop_vars.csv'))
    pop_vars = [itm[3] for itm in pop]

    now = datetime.now()
    scratchId = now.strftime("%Y%m%d%H%M")

    def __init__(self, preset='andres_local', time_as='sequence'):

        """
        stage is used to manage file organization, naming, and output in
        conjunction with query and build

        A few challenges of working with climate data:
        - Raw data files can be very large (and slow)
        - Filenames are hard to read and can look very similar
        - Hard to remember the details of a specific analysis or summary plot
        once the Python script has been closed (ex. what are the exact
        coordinates that were used in a spatial average?)

        stage addesses these by allowing you to define a directory structure
        and naming convention for both the source data and output of ocrtools.
        Most query functions rely on this

        Args:
        * preset (str) [optional]: specifies folder structure for existing data
        and new output. Please add your own preset and make it the default!
        * time_as (str) [optional]: 'sequence' (ex. base1920.34310, as in
        34310 days after 1920-000 [2014]) or 'date' (ex. 2014-000)
        """

        if preset == 'andres_local':
            path0 = "/Volumes/Samsung_T5/Open_Climate_Research-Projects"
            self.directories = {
                "cesm-raw": path0 + "/data/raw/cesm",
                "cesm-reformatted": path0 + "/data/reformatted/cesm",
                "other-reformatted": path0 + "/data/reformatted",
                "scratch": path0 + "/data/scratch",
                "plot": path0 + "/data/test_plots",
                "csv": path0 + "/data/streams"}
            self.subfolders = {
                "cesm-raw": ['dt', 'var_name'], "cesm-reformatted":
                ["var_name", "yr_range"], "other-reformatted":
                ["var_name", "yr_range"], "scratch": [], "plot": ["function"]}

        elif preset == 'ucar':
            self.directories = {
                "cesm-raw": "/glade/p_old/cesmLE/CESM-CAM5-BGC-LE",
                "cesm-reformatted": "/glade/p/work/andresc/OCR/reformatted/\
                cesm",
                "other-reformatted": "/glade/p/work/andresc/OCR/reformatted",
                "scratch": "/glade/p/work/andresc/OCR/scratch"}
            self.subfolders = {
                "cesm-raw": ['model', 'proc', 'tseries', 'dt', 'src_var_name'],
                "cesm-reformatted": ["var_name", "yr_range"], "scratch": []}
            import matplotlib
            matplotlib.use('Agg')

        elif preset == 'demo':
            self.directories = {
                "cesm-raw": "data/raw", "cesm-reformatted": "data/reformatted",
                "other-reformatted": "data/reformatted", "scratch":
                "data/scratch", "plot": "test_plots", "csv": "data/streams"}
            self.subfolders = {
                "cesm-raw": [], "cesm-reformatted": [],
                "other-reformatted": [], "plot": [], "scratch": []}

        elif preset == 'default':
            self.directories = None
            self.subfolders = {
                "cesm-raw": [], "cesm-reformatted": [], "other-reformatted":
                [], "plot": [], "scratch": []}

        self.naming = {
            "cesm-reformatted": ["var_name", "dt", "mem", "time_slice"],
            "other-reformatted": ["var_name", "dt", "time_slice"],
            "spatial_average":
            ["var_name",  "dt", "time_slice", "mean", "scratchId"],
            "cesm_report": ["var_name", "dt", "mem"],
            "other-report": ["var_name", "dt"]}

        self.time_as = time_as
        self.base_yr = 1920
