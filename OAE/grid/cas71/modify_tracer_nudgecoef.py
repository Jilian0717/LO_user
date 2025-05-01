# modify nudge_coef for tracer
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import scipy.io
from lo_tools import plotting_functions as pfun
import pickle
from lo_tools import zrfun

# old nudge file (generate from pgrid)
ds = xr.open_dataset('nudgcoef.nc')
#M2_NudgeCoef = ds.M2_NudgeCoef.values
#M3_NudgeCoef = ds.M3_NudgeCoef.values
#tracer_NudgeCoef = ds.tracer_NudgeCoef.values
#temp_NudgeCoef = ds.temp_NudgeCoef.values
#salt_NudgeCoef = ds.salt_NudgeCoef.values

alk_nudge_val = ds.tracer_NudgeCoef.values * 100
TIC_nudge_val = alk_nudge_val.copy()

vn = 'alkalinity_NudgeCoef'
dims = ('s_rho','eta_rho','xi_rho')
ds[vn] = (dims, alk_nudge_val)
ds[vn].attrs['long_name'] = "alkalinity inverse nudging coefficients"
ds[vn].attrs['units'] = 'day-1'

vn = 'TIC_NudgeCoef'
dims = ('s_rho','eta_rho','xi_rho')
ds[vn] = (dims, TIC_nudge_val)
ds[vn].attrs['long_name'] = "TIC inverse nudging coefficients"
ds[vn].attrs['units'] = 'day-1'

Enc_dict = {vn:zrfun.enc_dict for vn in ds.data_vars}

ds.to_netcdf('nudgcoef.nc.alk_TIC', encoding=Enc_dict)

#