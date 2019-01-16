
import numpy as np
import xarray as xr


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


def p2(x, x0=0.7, steepness=5, **kwargs):
    dice0 = 0.18
    dice_a = 4.4
    d_mult = dice_a/dice0

    # print(x)
    # print(x.dims)
    # print(len(x['time']))

    try:
        roll_dice = np.random.rand(x.dims['time'])
    except TypeError:
        roll_dice = np.random.rand(len(x['time']))

    x_else = ((1/(1-x0))*(x-x0))**steepness
    x_roll_under = roll_dice * dice_a
    x_roll_over = 0

    x_roll = xr.where(roll_dice < dice0, x_roll_under, x_roll_over)
    xf = xr.where(x < x0, x_roll, x_else)
    return(xf)


def print_dist(dist):
    import matplotlib.pyplot as plt
    xa = np.arange(0, 1, 0.01)
    x1 = [dist(x) for x in xa]
    plt.plot(xa, x1)
    plt.savefig('test_dist.png')


# Build settings
step_std_a = {'PRECT': 8, 'TS': 13, 'RAIN': 8, 'H2OSNO': 8}
blur_std_a = {'PRECT': 8, 'TS': 22, 'RAIN': 8, 'H2OSNO': 8}
snap = {'PRECT': 11.5, 'TS': 9.5, 'RAIN': 8, 'H2OSNO': 8}
snap_atten = {'PRECT': 8, 'TS': 10, 'RAIN': 8, 'H2OSNO': 8}
head = 3
tail = 3
combine_plots = False
hist_stretch = {'PRECT': True, 'TS': False, 'RAIN': True, 'H2OSNO': True}
hist_dist = {'PRECT': p2, 'TS': False, 'RAIN': p2, 'H2OSNO': p2}
var_min = {'PRECT': 0., 'TS': None, 'RAIN': 0., 'H2OSNO': 0.}
var_max = {'PRECT': None, 'TS': None, 'RAIN': None, 'H2OSNO': 0.}

build_kwargs = {
    'step_std_a': step_std_a,
    'blur_std_a': blur_std_a,
    'snap': snap,
    'snap_atten': snap_atten,
    'hist_stretch': hist_stretch,
    'hist_dist': hist_dist,
    'var_min': var_min,
    'var_max': var_max,
    'head': head,
    'tail': tail
}


# Conversion functions

def mps_2_cmpday(input_list):
    t_day = 60.*60.*24.   # seconds * minutes * hours
    lconvert = [n * t_day * 100 for n in input_list]
    return lconvert


def K_2_F(input_list):
    t_f = [t * 9./5 - 459.67 for t in input_list]
    return t_f


# Conversion and plotting settings
alabels = {'TS': 'Temperature (F)', 'PRECT': 'Precipitation (cm/day)',
           'RAIN': 'Precipitation (cm/day)'}
conversions = {'PRECT': mps_2_cmpday, 'TS': K_2_F, 'RAIN': mps_2_cmpday}
ylim = {'PRECT': [0, 20], 'TS': [0, 100], 'RAIN': [0, 20]}

