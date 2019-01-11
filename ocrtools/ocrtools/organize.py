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

directory_map = ['type', 'src', 'dt', 'file']
cesm_fname = ['compset', 'code_base', 'compset_short', 'res_short', 'desc',
              'nnn', 'scomp', 'type', 'string', 'date', 'ending']
cesmLE_map = {'compset': 'b', 'code_base': 'e11', 'res_short': 'f09_g16',
              'ending': 'nc'}


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


def gen_path(path_map, path_info, join='/'):
    """
    Returns a path
    Args:
    * path_map (list): Mapping order
    * path_info (dict): Dictionary with entries for each relevant map item
    """
    return(
        join.join([path_info[i] for i in path_map if i in path_info.keys()]))


def cesmLE_fname(var, dt, yr0, mem=0, hem=''):

    cesm_d = cesmLE_map

    if '{:04d}'.format(yr0)[2:4] == '00':
        cesm_d['compset_short'] = 'B1850C5CN'
        yrf = yr0 + 99
        if mem == 0:
            mem = 5
    else:
        if mem == 0:
            mem = 2
        if yr0 < 2006:
            cesm_d['compset_short'] = 'B20TRC5CNBDRD'
            yrf = 1999
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
        elif dt == 'monthly':
            cesm_d['type'] = 'h'
        if hem != 'nh' and hem != 'sh':
            hem = input('\n[OCR] Which hemisphere (nh or sh)? ')

        cesm_d['string'] = var + '_' + hem

    elif var in cam_vars or var in clm_vars:
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
