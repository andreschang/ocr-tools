#!/usr/bin/env python
# -*- coding: utf-8 -*-

########################################
#
#  OCR TOOLS - Making climate data malleable
#  Copyright (C) 2018 Andres Chang
#
########################################

import os
import numpy as np
import xarray as xr
from ocrtools.load import load

directory_map = ['type', 'src', 'dt', 'var', 'file']
cesm_fname = ['compset', 'code_base', 'compset_short', 'res_short', 'desc',
              'nnn', 'scomp', 'type', 'string', 'date', 'ending']
cesmLE_map = {'compset': 'b', 'code_base': 'e11', 'res_short': 'f09_g16',
              'ending': 'nc'}
top_directory = '/Volumes/Samsung_T5/Open_Climate_Research-Projects/data'


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


def gen_path(path_map, path_info, join='/', top=''):
    """
    Returns a path
    Args:
    * path_map (list): Mapping order
    * path_info (dict): Dictionary with entries for each relevant map item
    """
    if top is None:
        try:
            top = top_directory + '/'
        except NameError:
            top = ''
    elif top is '':
        top = ''
    else:
        top = top + '/'

    return(
        top +
        join.join([path_info[i] for i in path_map if i in path_info.keys()]))


def cesmLE_fname(var, dt, yr0, mem=0, hem=''):
    """
    Generates cesmLE filename
    Args:
    * var (str): name of variable
    * dt (str): monthly or daily
    * yr0 (numeric): first year of data file
    * mem (numeric): member number (if 0, defaults to 002 for BAU runs; and
    all control runs are 005)
    * hem (str): cice variables need to be specified as 'nh' or 'sh'
    """

    cesm_d = cesmLE_map

    if '{:04d}'.format(yr0)[2:4] == '00':
        cesm_d['compset_short'] = 'B1850C5CN'
        yrf = yr0 + 99
        if mem == 0:
            mem = 5
    else:
        if mem == 0:
            mem = 2
        if yr0 == 1920:
            cesm_d['compset_short'] = 'B20TRC5CNBDRD'
            yrf = 2005
        else:
            cesm_d['compset_short'] = 'BRCP85C5CNBDRD'
            if yr0 == 2006:
                yrf = 2080
            else:
                yrf = 2100

    cesm_d['nnn'] = '{:03d}'.format(mem)

    if dt == 'daily':
        cesm_d['date'] = ('{:04d}'.format(yr0) + '0101' + '-' +
                          '{:04d}'.format(yrf) + '1231')
    elif dt == 'monthly':
        cesm_d['date'] = ('{:04d}'.format(yr0) + '01' + '-' +
                          '{:04d}'.format(yrf) + '12')

    if var in cice_vars:
        cesm_d['scomp'] = 'cice'
        if dt == 'daily':
            cesm_d['type'] = 'h1'
            var = var + '_d' if not var.endswith('_d') else var
        elif dt == 'monthly':
            cesm_d['type'] = 'h'
        if hem != 'nh' and hem != 'sh':
            hem = input('\n[OCR] Which hemisphere (nh or sh)? ')

        cesm_d['string'] = var + '_' + hem

    else:
        cesm_d['string'] = var
        if var in cam_vars or var in clm_vars:
            if var in cam_vars:
                cesm_d['scomp'] = 'cam'
            else:
                cesm_d['scomp'] = 'clm2'
            if dt == 'daily':
                cesm_d['type'] = 'h1'
            elif dt == 'monthly':
                cesm_d['type'] = 'h0'

        elif var in pop_vars:
            cesm_d['scomp'] = 'pop'
            if dt == 'daily':
                cesm_d['type'] = 'h.nday1'
            elif dt == 'monthly':
                cesm_d['type'] = 'h0'

    return(gen_path(cesm_fname, cesm_d, '.'))


def mkdir_p(path):
    """
    Makes a new directory, if does not already exist
    """
    import errno

    try:
            os.makedirs(path)
    except OSError as exc:  # Python >2.5
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                    pass
            else:
                raise


def get_ncs(f_dir):
    """
    Returns list of netcdf filenames in directory
    """
    all_ncs = []
    for file in os.listdir(f_dir):
        if file.endswith(".nc"):
            all_ncs.append(file)
    return all_ncs


def find_nearest(array, value, index=False):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    if index:
        return idx
    else:
        return array[idx]


def load_cesmLE(var, dt, yr0, yrf, mem, **kwargs):
    """
    Convenience function that loads cesm LE data, appending multiple
    files if needed.

    Args:
    * var (str): name of variable
    * dt (str): monthly or daily
    * yr0 (numeric): first year of inquiry
    * yrf (numeric): last year of inquiry
    """
    raw_cesm = {'type': 'raw', 'src': 'cesm', 'dt': dt, 'var': var}
    f0 = gen_path(directory_map, raw_cesm, top=None)
    raw_ncs = get_ncs(f0)
    raw_yr0 = [int(x.split('/')[-1]
                    .split('.')[-2]
                    .split('-')[0][0:4]) for x in raw_ncs if
               x.split('/')[-1]
                .split('.')[4] == '{:03d}'.format(mem)]

    raw_yrf = [int(x.split('/')[-1]
                    .split('.')[-2]
                    .split('-')[1][0:4]) for x in raw_ncs if
               x.split('/')[-1]
                .split('.')[4] == '{:03d}'.format(mem)]

    f_cesm = []
    f_yr0, f_yrf = yr0, -1

    while f_yrf < yrf:
        opts = [(x, y) for x, y in zip(raw_yr0, raw_yrf) if x <= f_yr0]
        if len(opts) == 0:
            if mem == 5:
                guess_yr0 = round(yr0-50, -2)
            else:
                if yr0 < 2006:
                    guess_yr0 = 1920
                else:
                    guess_yr0 = 2006
            raise ValueError(
                '\n[OCR]File not found in expected location, ' +
                f0 + '. Best guess for missing filename is ' +
                cesmLE_fname(var, dt, guess_yr0, mem, **kwargs))

        ni = find_nearest([x[0] for x in opts], f_yr0, index=True)
        f_yr0, f_yrf = opts[ni]
        fname = cesmLE_fname(var, dt, f_yr0, mem, **kwargs)
        f_cesm.append(fname)
        f_yr0 = f_yrf + 1

    print('\n[OCR] Requested CESM LE files were successfully found in ' + f0)
    load_cesm = [load(f0 + '/' + x, var=x.split('.')[7]
                                         .replace('_nh', '')
                                         .replace('_sh', '')) for x in f_cesm]
    t0 = '{:04d}'.format(yr0) + '-01-01'
    tf = '{:04d}'.format(yrf) + '-12-31'

    return(xr.concat(load_cesm, 'time').sel(time=slice(t0, tf)))


# Additional functions to be written -- save_formatted
