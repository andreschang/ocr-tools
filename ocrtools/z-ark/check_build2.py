import ocrtools as ocr
from ocrtools import plt, conversions
import xarray as xr
import numpy as np
import pandas as pd

dt = 'monthly'
location = 'Wynnewood, PA USA'
d0 = ocr.load(
    '/Volumes/Samsung_T5/Open_Climate_Research-Projects/prj-EWYHF_Exactly_Where_You_Had_Fallen/data/WynnewoodPA_data/TREFHT_H2OSNO.1958-1981.daily.0002.nc',
    standardize=False, var=['TREFHT', 'H2OSNO'])
d1 = ocr.load(
    '/Volumes/Samsung_T5/Open_Climate_Research-Projects/prj-EWYHF_Exactly_Where_You_Had_Fallen/data/WynnewoodPA_data/TREFHT_H2OSNO.1958-1981.daily.1850.nc',
    standardize=False, var=['TREFHT', 'H2OSNO'])

var_lims = {'TREFHT': [0, 500], 'H2OSNO': [0, 100]}
b = ocr.build(d0, var_lims, d1, combine_steps=5)
new = b.new(x_vars=['TREFHT'], y_vars=['H2OSNO'])
new['TREFHT'].plot()
d0['TREFHT'].plot()
plt.show()
plt.close()
new['H2OSNO'].fillna(0).plot()
d0['H2OSNO'].plot()
plt.show()
plt.close()
