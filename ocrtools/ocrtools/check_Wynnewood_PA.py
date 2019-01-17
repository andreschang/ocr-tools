import ocrtools as ocr
from ocrtools import plt, conversions


wpa = ocr.scope(location='Wynnewood, PA USA')
TS_daily = ocr.load_cesmLE('TS', 'daily', 1960, 1980, 2)
TS_wpa = ocr.subset(TS_daily, wpa)


f = conversions['K_2_F']
f(TS_wpa['TS']).plot()
plt.show()