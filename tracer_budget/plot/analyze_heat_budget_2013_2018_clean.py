import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from lo_tools import zrfun, zfun
import pandas as pd
import pickle, scipy, get_two_layer
from time import time
from datetime import datetime
from pathlib import Path
from lo_tools import plotting_functions as pfun
from matplotlib.ticker import ScalarFormatter

Cp = 3985 # Joules/kg/degC
rho0 = 1025 # kg/m3

yr1 = 2013
yr2 = 2018

#%% tef2 extracted tracer flux  -- jdf1
dir0 = Path('/Users/jilian/Desktop/LiveOcean/LO_output/extract/cas7_t0_x4b/tef2/bulk_'\
            + str(yr1) + '_' + str(yr2) + '/');
fn_list = sorted(dir0.glob('jdf1*.nc'))
bulk = xr.open_mfdataset(fn_list, concat_dim="time", combine='nested')
tef_df, vn_list, vec_list = get_two_layer.get_two_layer(bulk)
temp_p_jdf1 = tef_df['temp_p'] # Tout
temp_m_jdf1 = tef_df['temp_m'] # Tin
temp_p_flux_jdf1 = (tef_df['temp_p'])*tef_df['q_p'] * Cp*rho0  # m3/s degC to J/s or Watts
temp_m_flux_jdf1 = (tef_df['temp_m'])*tef_df['q_m'] * Cp*rho0

#%----------------------------- sog6 -------------------------------
fn_list = sorted(dir0.glob('sog6*.nc'))
bulk = xr.open_mfdataset(fn_list, concat_dim="time", combine='nested')
tef_df, vn_list, vec_list = get_two_layer.get_two_layer(bulk)
temp_p_sog6 = tef_df['temp_p'] # Tout
temp_m_sog6 = tef_df['temp_m'] # Tin
temp_p_flux_sog6 = tef_df['temp_p']*tef_df['q_p'] * Cp*rho0  # m3/s degC to J/s or Watts
temp_m_flux_sog6 = tef_df['temp_m']*tef_df['q_m'] * Cp*rho0

# add jdf1 and sog6 together
temp_p_flux = temp_p_flux_jdf1 + temp_p_flux_sog6
temp_m_flux = temp_m_flux_jdf1 + temp_m_flux_sog6

QT_in  = -temp_m_flux  
QT_out = -temp_p_flux

#%% riverine heat input
seg_name = '/Users/jilian/Desktop/LiveOcean/LO_output/extract/tef2/seg_info_dict_cas7_c2_trapsV00.p'
seg_df = pd.read_pickle(seg_name)
ji_list = seg_df['sog6_m']['ji_list']
jj = [x[0] for x in ji_list]
ii = [x[1] for x in ji_list]
# all rivers inside the Salish segment
riv_list = seg_df['sog6_m']['riv_list']

# all rivers inside LO domain
fn = '/Users/jilian/Desktop/LiveOcean/LO_output/forcing/cas7/f2013.01.01/traps00/rivers.nc'
ds = xr.open_dataset(fn)
river_name_all = ds.river_name.values

# river data
fn_riv0 = '/Users/jilian/Desktop/LiveOcean/LO_output/pre/river1/cas7_traps00/Data_roms/extraction_'\
        + str(2013) + '.01.01_' + str(2015) + '.12.31.nc'
fn_riv1 = '/Users/jilian/Desktop/LiveOcean/LO_output/pre/river1/cas7_trapsF00/Data_roms/extraction_'\
        + str(2016) + '.01.01_' + str(2018) + '.12.31.nc'
fn_riv = np.array([fn_riv0, fn_riv1])
ds_riv = xr.open_mfdataset(fn_riv)
t_riv = ds_riv.time.values # daily
temp_riv = ds_riv.temp.values # river temperature
Q_riv = ds_riv.transport.values # river discharge
temp_riv_noon = dict()

for riv in riv_list:
    ix=np.where(riv==river_name_all)[0]
    temp_riv_noon[riv] = temp_riv[:,ix[0]] * Q_riv[:,ix[0]] * Cp*rho0 # riverine heat flux
# add all riverine heat input together
QT_riv = np.zeros(len(temp_riv_noon[riv])) 
for riv in riv_list:
    QT_riv += temp_riv_noon[riv]
    
#%% load air-sea heat flux and Temp*Vol (hourly)
fn_list = []
for yr in range(yr1,yr2+1):
    dir0 = Path('./data/TNTempVol_Denitri_AirSeaHeat_'+str(yr)+'/'); 
    fn_list0 = sorted(dir0.glob('*.p'))
    fn_list += fn_list0

t = []; 
shflux_sum = []; swrad_sum = []; lwrad_sum =[]; latent_sum=[]; sensible_sum=[]
temp_vol = []
for fn in fn_list:
    tmp = pickle.load(open(fn, 'rb'))
    shflux_sum += tmp['shflux_sum'] # watts or J/s
    swrad_sum  += tmp['swrad_sum']
    lwrad_sum  += tmp['lwrad_sum']
    latent_sum += tmp['latent_sum']
    sensible_sum += tmp['sensible_sum']   
    temp_vol += tmp['temp_vol']
    t += tmp['t']
# change data type
shflux_sum = np.array(shflux_sum)
swrad_sum = np.array(swrad_sum)
lwrad_sum = np.array(lwrad_sum)
sensible_sum = np.array(sensible_sum)
latent_sum = np.array(latent_sum)
temp_vol = np.array(temp_vol)  # Temperature*vol
temp_rate = np.diff(temp_vol)/3600 * Cp*rho0 # m3/s degC to Watts
# low pass using godin filter, then re-sampled to daily at noon
shflux_sum_lp = zfun.lowpass(shflux_sum, f='godin')[36:-34:24]
swrad_sum_lp = zfun.lowpass(swrad_sum, f='godin')[36:-34:24]
lwrad_sum_lp = zfun.lowpass(lwrad_sum, f='godin')[36:-34:24]
latent_sum_lp = zfun.lowpass(latent_sum, f='godin')[36:-34:24]
sensible_sum_lp = zfun.lowpass(sensible_sum, f='godin')[36:-34:24]
temp_rate_lp = zfun.lowpass(temp_rate, f='godin')[36:-34:24]
t_lp = t[36:-34:24]

#%%
                                        #---- plot ----#
Qnet = QT_in + QT_out + QT_riv + shflux_sum_lp
Qnet_E = QT_in + QT_out # net exchange flow term

fig = plt.figure(figsize=(14, 8));  fz=12; 
# wind
ds = xr.open_dataset('/Users/jilian/Desktop/LiveOcean/LO_output/extract/cas7_t0_x4b/moor/wind_stress/shelf_'+str(yr1)+'.01.01_'+str(yr2)+'.12.31.nc')
sustr = ds.sustr.values
svstr = ds.svstr.values
t_wind = ds.ocean_time.values
angle = 20*np.pi/180
along = -sustr * np.sin(angle) + svstr * np.cos(angle)
cross = sustr * np.cos(angle) + svstr * np.sin(angle)
along_filt = zfun.filt_AB8d(along)
yy = np.zeros(len(t_wind))
ax = fig.add_axes([0.08, 0.75, 0.86, 0.18]) # wind stress
ax.fill_between(t_wind, along_filt, yy, color='g', alpha=.5, lw=0.1, where=(along_filt>=yy))
ax.fill_between(t_wind, along_filt, yy, color='k', alpha=.5, lw=0.1, where=(along_filt<yy))
ax.set_ylim([-0.4, 0.4])
for yrr in range(yr1, yr2+1):
    plt.vlines(x=datetime(yrr,1,1), ymin=-90, ymax=90, color='k', ls='--', lw=0.5)
plt.xlim(datetime(yr1+1,1,1), datetime(yr2+1,1,1))
plt.ylabel('Wind Stress (Pa)', fontsize=12)
ax.tick_params(axis='y', labelsize=12)
plt.xticks([])

#---- heat budget component
ax = fig.add_axes([0.08, 0.4, 0.86, 0.35]); N = 10 # day
plt.plot(t_lp, zfun.lowpass(Qnet_E.values/1e11, f='hanning', n=N), 'b',
         label='QT$_{exchange}$', lw=2)
plt.plot(t_lp, zfun.lowpass(shflux_sum_lp/1e11, f='hanning', n=N), 'peru', 
         label='QT$_{air-sea}$', lw=2)
plt.plot(t_riv, zfun.lowpass(QT_riv/1e11, f='hanning', n=N), 'lightseagreen', 
         label='QT$_{river}$', lw=2)
plt.plot(t_lp, zfun.lowpass(temp_rate_lp/1e11, f='hanning', n=N), 'pink', 
         label='d(T*V)/dt', lw=1.5)
error = zfun.lowpass((Qnet.values-temp_rate_lp)/1e11, f='hanning', n=N)
plt.plot(t_lp, error, '--', c='gray', lw=1.5, label='Error')
for yrr in range(yr1+1, yr2+1):
    plt.vlines(x=datetime(yrr,1,1), ymin=-90, ymax=90, color='k', ls='--', lw=0.5)
plt.ylabel('Heat flux (10$^{11}$W)', fontsize=fz)
plt.xlim(datetime(yr1+1,1,1), datetime(yr2+1,1,1)); 
plt.ylim(-70, 70)
plt.legend(loc='upper left', fontsize=10, ncol=3, frameon=False)
plt.xticks(fontsize=fz); plt.yticks(fontsize=fz)
plt.xticks([])

#% vol-avg temperature
fn_list = []
for yr in range(yr1, yr2+1):
    dir0 = Path('./data/vol_weighted_tracer_conc_'+str(yr)+'/');
    fn_list0 = sorted(dir0.glob('*.p'))
    fn_list += fn_list0
t1=[] 
temp_all_depth = []
temp_shallow = []
temp_deep = []
for fn in fn_list:
    tmp = pickle.load(open(fn, 'rb'))
    t1 += tmp['t']
    temp_all_depth += tmp['temp_all_depth']
    temp_shallow += tmp['temp_shallow']
    temp_deep += tmp['temp_deep']
temp_all_depth = np.array(temp_all_depth)
temp_shallow = np.array(temp_shallow)
temp_deep = np.array(temp_deep)
temp_all_depth_lp = zfun.lowpass(temp_all_depth, f='godin')[36:-30:24]
temp_shallow_lp = zfun.lowpass(temp_shallow, f='godin')[36:-30:24]
temp_deep_lp = zfun.lowpass(temp_deep, f='godin')[36:-30:24]
t1_noon = t1[36:-30:24]

ax = fig.add_axes([0.08, 0.05, 0.86, 0.35])
plt.plot(t1_noon, zfun.lowpass(temp_shallow_lp, f='hanning', n=N), 'g', 
         label='Shallow (<20m)')
plt.plot(t1_noon, zfun.lowpass(temp_deep_lp, f='hanning', n=N), 'k', 
         label='Deep (>20m)')
plt.plot(t1_noon, zfun.lowpass(temp_all_depth_lp, f='hanning', n=N), 'peru', 
         label='All depth')
for yrr in range(yr1+1, yr2+1):
    plt.vlines(x=datetime(yrr,1,1), ymin=6, ymax=16, color='k', ls='--', lw=0.5)
plt.legend(frameon=False, loc='upper left')
plt.ylim([6,14])
plt.xlim(datetime(yr1+1,1,1), datetime(yr2+1,1,1))
plt.ylabel('Vol-weigthed temperature (˚C)', fontsize=12)
ax.tick_params(axis='y', labelsize=12)
ax.tick_params(axis='x', labelsize=12)
plt.tight_layout() 
plt.show()

