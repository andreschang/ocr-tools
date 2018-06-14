import ocrtools.stage as st
import ocrtools.query as query


## Directory structure and naming conventions are set
## with ocrtools.stage.

demo_stage = st(preset = 'demo', time_as = 'date')



# ice = qy(stage = stage)

# # ice.set_params('/Volumes/Samsung_T5/Open_Climate_Research-Projects/data/raw/cesm/monthly/aice/b.e11.B1850C5CN.f09_g16.005.cice.h.aice_nh.200001-209912.nc')
# ice.set_params('/Volumes/Samsung_T5/Open_Climate_Research-Projects/data/raw/air.mon.mean.nc', dt = 'monthly', \
#   data_yr0 = 1948, var = 'air', dim = ['time', 'lat', 'lon'])
# # ice.extract_data()
# ice.reduce_data()

# ice.set_params(src = 'cesm', cesm_mem = 2, cesm_hemisphere = 'nh', variable = 'aice', dt = 'monthly')

ts = query(stage = stage)
ts.set_params(var = 'TS', src = 'cesm', cesm_mem = 2, dt = 'monthly')
# ts.reduce_data()
ts.spatial_average(lat_bounds = [-20, 20], lon_bounds = [50, 60], yr0 = 2010, yrf = 2015)