# calculate hourly volume, hourly salt*volume and EminusP
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from lo_tools import zrfun, Lfun
import pickle
from datetime import datetime, timedelta
from time import time
import sys
import pandas as pd

tt0 = time()

#% load salish sea j,i
seg_name = 'seg_info_dict_cas7_c2_noriv.p'
seg_df = pd.read_pickle(seg_name)
ji_list = seg_df['sog6_m']['ji_list']
jj = [x[0] for x in ji_list]
ii = [x[1] for x in ji_list]

#
Ldir = Lfun.Lstart()
Ldir['roms_out'] = Ldir['roms_out2']
Ldir['gtagex'] = 'cas7_t0_x4b'

ds0 = '2017.01.01'
ds1 = '2017.01.31'
Ldir['ds0'] = ds0
in_dir = Ldir['roms_out'] / Ldir['gtagex']
G, S, T = zrfun.get_basic_info(in_dir / ('f' + Ldir['ds0']) / 'ocean_his_0002.nc')

fn0 = xr.open_dataset(in_dir / ('f' + Ldir['ds0']) / 'ocean_his_0002.nc')
dx = 1/fn0.pm.values
dy = 1/fn0.pn.values
area = dx * dy

dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
dt1 = datetime.strptime(ds1, Lfun.ds_fmt)
dt00 = dt0

t = []
SV = [] # sum(salt*vol) in the whole domain, hourly
vol_hrly = []  # total volume in the domain, hourly
surf_s_flux = []  # EminusP: m s-1, see /ROMS/External/varinfo.ymal for more descriptions

while dt00 <= dt1:
    print(dt00)
    sys.stdout.flush()
    ds00 = dt00.strftime(Lfun.ds_fmt)
    fn_list = Lfun.get_fn_list('hourly', Ldir, ds00, ds00)
    
    for fn in fn_list[0:-1]: 
        ds_his = xr.open_dataset(fn)
        EminusP = ds_his.EminusP.values.squeeze()  # EminusP
        salt_surf = ds_his.salt.values.squeeze()[-1,:,:]  # surface salinity
        #print(fn)
        h = ds_his.h.values      
        zeta = ds_his.zeta.values.squeeze()
        zw = zrfun.get_z(h, zeta, S, only_w=True)
        dz = np.diff(zw, axis=0)
        vol = dx*dy*dz
            
        #-------- salt*vol --------
        salt_tmp = ds_his.salt.values.squeeze()
        SV.append(np.nansum(salt_tmp[:,jj,ii] * vol[:,jj,ii]))  # salt*vol
        vol_hrly.append(np.nansum(vol[:,jj,ii]))
        tmp = salt_surf*EminusP*area  # also see Parker's github /LO/extract/tef2/tracer_budget.py line 177
        surf_s_flux.append(np.nansum(tmp[jj,ii]))
        t.append(ds_his.ocean_time.values)
        
    dt00 = dt00 + timedelta(days=1)
    
dict_tmp = {'t': t,
            'salt_vol_sum_hrly': SV,
            'vol_hrly': vol_hrly,
            'surf_s_flux': surf_s_flux
           }
pickle.dump(dict_tmp, open("vol_SV_hrly_"+ds0+'_'+Ldir['gtagex']+'.p','wb'))
