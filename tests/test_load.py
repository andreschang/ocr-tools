import ocrtools as ocr
import numpy as np

def test_load_var():
    f1 = 'b.e11.B1850C5CN.f09_g16.005.cice.h.aice_nh.190001-199912.nc'
    load1 = ocr.load('../data/raw/' + f1, 'aice')
    assert ((list(load1.data_vars) == ['aice']) and
            'lat' in list(load1.coords.keys()) and
            'lon' in list(load1.coords.keys()))

def test_scope_location():
    scope = ocr.scope(location='Lima, Peru', interactive=False)
    assert round(scope.lat) == -12 and round(scope.lon) == 103


def test_scope_lat_min():
    scope = ocr.scope(lat_min=20, lon_min=20, lon_max=40, interactive=False)
    assert scope.lat_min == 20 and scope.lat_max == 90


def test_subset():
    f1 = 'b.e11.B1850C5CN.f09_g16.005.cice.h.aice_nh.190001-199912.nc'
    load1 = ocr.load('../data/raw/' + f1, 'aice')
    scope = ocr.scope(lat_min=70., lon_min=-20., lon_max=20., interactive=False)
    subset = ocr.subset(load1, scope)
    # In this example with 2D lat/lon grid, some dim pairs exit outside of scope
    # range and are assigned NaN
    assert (np.amin(subset.lat.values) > 65)
