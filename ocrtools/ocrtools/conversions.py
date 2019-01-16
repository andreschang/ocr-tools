import numpy as np
import xarray as xr

def mps_2_cmpday(input_list):
    t_day = 60.*60.*24.   # seconds * minutes * hours
    lconvert = [n * t_day * 100 for n in input_list]
    return lconvert


def K_2_F(input_list):
    t_f = [t * 9./5 - 459.67 for t in input_list]
    return t_f


def p0(x):
    if x < 0.6:
        return(x**2)
    else:
        return(0.36+((x-0.6)/0.4))


def p1(x):
    x0 = 0.7

    if x < 0.7:
        return 0
    else:
        return ((1/(1-x0))*(x-x0))**3


def p2(x, x0=0.7, steepness=5):
    dice0 = 0.18
    dice_a = 4.4
    d_mult = dice_a/dice0

    roll_dice = np.random.rand(x.dims['time'])
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


