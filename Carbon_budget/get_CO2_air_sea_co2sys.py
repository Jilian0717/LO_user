"""
Get air-sea CO2 flux based on saved history files
"""
import xarray as xr
import numpy as np
from lo_tools import zfun, zrfun, Lfun
from datetime import datetime, timedelta
import sys
import PyCO2SYS as pyco2  # https://pyco2sys.readthedocs.io/en/latest/co2sys_nd/
import gsw
import pandas as pd

Ldir = Lfun.Lstart()
Ldir['roms_out'] = Ldir['roms_out2']
Ldir['gtagex'] = 'cas7_t0_x4b'

ds0 = '2014.01.01'
ds1 = '2014.01.31'
Ldir['ds0'] = ds0
in_dir = Ldir['roms_out'] / Ldir['gtagex']
G, S, T = zrfun.get_basic_info(in_dir / ('f' + Ldir['ds0']) / 'ocean_his_0002.nc')

fn0 = xr.open_dataset(in_dir / ('f' + Ldir['ds0']) / 'ocean_his_0002.nc')
dx = 1/fn0.pm.values
dy = 1/fn0.pn.values
area = dx * dy
NX, NY = dx.shape

dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
dt1 = datetime.strptime(ds1, Lfun.ds_fmt)
dt00 = dt0

A_CO2 = 2073.1      # Schmidt
B_CO2 = 125.62      # number
C_CO2 = 3.6276      # transfer
D_CO2 = 0.043219    # coefficients
E_CO2 = 0.0

A1 = -60.2409       # surface
A2 = 93.4517        # CO2
A3 = 23.3585        # solubility
B1 = 0.023517       # coefficients
B2 = -0.023656
B3 = 0.0047036

pi2 = 6.2831853071796

D0 = 282.6          # coefficients
D1 = 0.125          # to calculate
D2 = -7.18          # secular trend in
D3 = 0.86           # atmospheric pCO2
D4 = -0.99
D5 = 0.28
D6 = -0.80
D7 = 0.06

# initialize
t = []
CO2_flux = np.zeros((24*32, NX, NY)) # CO2 air-sea flux
pCO2     = np.zeros((24*32, NX, NY)) # pCO2 (in the surface layer)
TA       = np.zeros((24*32, NX, NY)) # total alkalinity in umol/kg

dtdays = 3600 * (1.0/86400) / 1
cff2   = dtdays * 0.31 * 24.0/100.0

cnt = 0
   
while dt00 <= dt1:  # loop every day and each history file
    print(dt00)
    sys.stdout.flush()
    ds00 = dt00.strftime(Lfun.ds_fmt)
    year = dt00.year
    yday = dt00 - datetime(year,1,1)
    yday = yday.days
    fn_list = Lfun.get_fn_list('hourly', Ldir, ds00, ds00)
    
    for fn in fn_list[0:-1]: 
        ds = xr.open_dataset(fn)
        zeta = ds.zeta.values.squeeze()
        h = ds.h.values
        lon = ds.lon_rho.values
        lat = ds.lat_rho.values
        mask_rho = ds.mask_rho.values # 1 = water, 0 = land
     
        z_w = zrfun.get_z(h, zeta, S, only_rho=False, only_w=True)
        vol = np.diff(z_w, axis=0) * area  # grid cell volume
                
        #---------- air-sea CO2 flux ----------
        Uwind = ds.Uwind.values.squeeze()
        Vwind = ds.Vwind.values.squeeze()
        temp_surf = ds.temp.values[0,-1,:,:] # surface temp
        salt_surf = ds.salt.values[0,-1,:,:] # surface salt
        TAlk_surf = ds.alkalinity.values[0,-1,:,:] # milliequivalents m-3
        TIC_surf = ds.TIC.values[0,-1,:,:]  #  millimole_carbon m-3
        
        NI, NJ = TIC_surf.shape
        
        # Compute CO2 transfer velocity: u10squared (u10 in m/s)
        u10squ = Uwind * Uwind + Vwind * Vwind  # ifdef BULK_FLUXES
        SchmidtN = A_CO2 - temp_surf*(B_CO2 - temp_surf*(C_CO2 - temp_surf*(D_CO2 - temp_surf*E_CO2)))
        cff3 = cff2 * u10squ * np.sqrt(660.0/SchmidtN)  # meter?
        
        #  Calculate CO2 solubility [mol/(kg.atm)] using Weiss (1974) formula.
        TempK = 0.01 * (temp_surf + 273.15)
        CO2_sol = np.exp(A1 + A2/TempK + A3*np.log(TempK) + salt_surf*(B1 + TempK*(B2 + B3*TempK)))

        # pCO2air_secular
        pmonth = year - 1951.0 + yday/365.0  # months since Jan 1951
        #if defined PCO2AIR_SECULAR
        pCO2air_secular = D0 + D1*pmonth*12.0 + D2*np.sin(pi2*pmonth + D3) + D4*np.sin(pi2*pmonth + D5) + D6*np.sin(pi2*pmonth + D7)
                
        zz = 0*h # layer depth
        pres = gsw.p_from_z(zz, lat) # pressure [dbar]
        SA = gsw.SA_from_SP(salt_surf, pres, lon, lat) # absolute salinity
        CT = gsw.CT_from_pt(SA, temp_surf) # conservative temperature
        rho = gsw.rho(SA, CT, pres) # in situ density
        temp = gsw.t_from_CT(SA, CT, pres) # in situ temperature
        # convert from umol/L to umol/kg using in situ density
        alkalinity = 1000 * TAlk_surf / rho # umol/kg
        alkalinity[alkalinity < 100] = np.nan
        TIC = 1000 * TIC_surf / rho # umol/kg
        TIC[TIC < 100] = np.nan
        # See LPM/co2sys_test/test0.py for info.
        CO2dict = pyco2.sys(par1=alkalinity, par2=TIC, par1_type=1, par2_type=2, salinity=salt_surf, temperature=temp_surf, pressure=pres, total_silicate=50, total_phosphate=2, opt_pH_scale=1, opt_k_carbonic=10, opt_k_bisulfate=1)
        pCO2[cnt,:,:] = CO2dict['pCO2'] # uatm
        CO2_flux[cnt,:,:] = cff3 * CO2_sol * (pCO2air_secular - pCO2[cnt,:,:])
        TA[cnt,:,:] = alkalinity # umol/kg
        cnt += 1
        t.append(ds.ocean_time.values)
        ds.close()       
    dt00 = dt00 + timedelta(days=1)

# save netcdf  - only save CO2_flux
from netCDF4 import Dataset
nc = Dataset('co2sys_CO2_air_sea_'+ds0+'.nc','w')
time = nc.createDimension('time', len(t))
eta_rho = nc.createDimension('eta_rho', NX)
xi_rho  = nc.createDimension('xi_rho', NY)
s_rho   = nc.createDimension('s_rho', 30)

times = nc.createVariable('time','f8',('time',))
times.units = 'seconds*1e9 since 1970-01-01 00:00:00'
CO2_flux_tmp = nc.createVariable('CO2_flux','f4', ('time',),compression='zlib',complevel=9)
CO2_flux_tmp.units = 'mmol CO2/hr'

times[:] = t
#% load salish sea index
seg_name = 'seg_info_dict_cas7_c2_noriv.p'
seg_df = pd.read_pickle(seg_name)
ji_list = seg_df['sog6_m']['ji_list']
jj = [x[0] for x in ji_list]
ii = [x[1] for x in ji_list]
tmp = CO2_flux[0:len(t),jj,ii] * area[jj,ii]
CO2_flux_tmp[:] = np.nansum(tmp,axis=1)  # only save integrated CO2 flux over the whole Salish Sea domain

nc.close()
