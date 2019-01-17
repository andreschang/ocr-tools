
import numpy as np
import xarray as xr
from ocrtools import plt


# Dist functions
def p0(x, **kwargs):
    if x < 0.6:
        return(x**2)
    else:
        return(0.36+((x-0.6)/0.4))


def p1(x, **kwargs):
    x0 = 0.7

    if x < 0.7:
        return 0
    else:
        return ((1/(1-x0))*(x-x0))**3


def p2(x, knee=0.7, steepness=5, p=0.18, dice_a=1, **kwargs):

    d_mult = dice_a/p

    try:
        roll_dice = np.random.rand(x.dims['time'])
    except TypeError:
        roll_dice = np.random.rand(len(x['time']))

    x_else = ((1/(1-knee))*(x-knee))**steepness
    x_roll_under = roll_dice * d_mult
    x_roll_over = 0

    x_roll = xr.where(roll_dice < p, x_roll_under, x_roll_over)
    xf = xr.where(x < knee, x_roll, x_else)

    return(xf)


def print_dist(dist):
    import matplotlib.pyplot as plt
    xa = xr.DataArray(np.arange(0, 1, 0.01), dims=['time']) 
    x1 = dist(xa)
    plt.plot(xa, x1)
    plt.savefig('test_dist.png')


# Build settings
step_std_a = {'PRECT': 10, 'TS': 13, 'RAIN': 15, 'H2OSNO': 40}
blur_std_a = {'PRECT': 10, 'TS': 22, 'RAIN': 15, 'H2OSNO': 40}
snap = {'PRECT': 1.2, 'TS': 1.2, 'RAIN': 0.8, 'H2OSNO': 0.4}
snap_atten = {'PRECT': 1, 'TS': 8, 'RAIN': 1, 'H2OSNO': 0}
head = 2
tail = 2
combine_plots = False
hist_stretch = {'PRECT': True, 'TS': False, 'RAIN': True, 'H2OSNO': True}
hist_dist = {'PRECT': p2, 'TS': False, 'RAIN': p2, 'H2OSNO': p2}
hist_args = {'PRECT': {}, 'TS': {}, 'RAIN': {}, 'H2OSNO': {}}
var_min = {'PRECT': 0., 'TS': None, 'RAIN': 0., 'H2OSNO': 0.}
var_max = {'PRECT': None, 'TS': None, 'RAIN': None, 'H2OSNO': None}
savgol_window = {'PRECT': 5, 'TS': 3, 'H2OSNO': 3, 'RAIN': 3}


build_kwargs = {
    'step_std_a': step_std_a,
    'blur_std_a': blur_std_a,
    'snap': snap,
    'snap_atten': snap_atten,
    'hist_stretch': hist_stretch,
    'hist_dist': hist_dist,
    'hist_args': hist_args,
    'var_min': var_min,
    'var_max': var_max,
    'head': head,
    'tail': tail,
    'savgol_window': savgol_window
}


# Conversion functions

def mps_2_cmpday(data):
    t_day = 60.*60.*24.   # seconds * minutes * hours
    lconvert = data * t_day * 100
    return lconvert


def K_2_F(data):
    t_f = data * 9./5 - 459.67
    return t_f


def mmH2O_2_inSNO(data):
    # http://webarchiv.ethz.ch/arolla/Arolla_Data/SnowConditions/depth_to_swe.pdf
    # https://bit.ly/2Cpc4nR
    H2O_row = 997  # kg/m^3
    SNO_row = 100  # kg/m^3
    mm_2_inch = 0.0393701
    mmSNO = (mm_2_inch * H2O_row/SNO_row) * data
    return(mmSNO)


def mmps_2_cmday(data):
    return mps_2_cmpday(data/1000.)

# Conversion and plotting settings
alabels = {'TS': 'Temperature (F)', 'PRECT': 'Precipitation (cm/day)',
           'RAIN': 'Precipitation (cm/day)', 'H2OSNO': 'Snow cover (in)'}
conversions = {'PRECT': mps_2_cmpday, 'TS': K_2_F, 'RAIN': mmps_2_cmday,
               'H2OSNO': mmH2O_2_inSNO}
ylim = {'PRECT': [0, 20], 'TS': [0, 100], 'RAIN': [0, 20], 'H2O_SNO': [0, 16]}

