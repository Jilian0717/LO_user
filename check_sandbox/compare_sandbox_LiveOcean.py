"""
compare NOAA Sandbox and LiveOcean 
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pinfo

def dar(ax):
    """
    Fixes the plot aspect ratio to be locally Cartesian.
    """
    yl = ax.get_ylim()
    yav = (yl[0] + yl[1])/2
    ax.set_aspect(1/np.cos(np.pi*yav/180))

var = 'oxygen'
# Sandbox 
fn0 = './f2023.07.04/ocean_his_0002.nc'
ds0 = xr.open_dataset(fn0)
if var == 'oxygen':
    vn0 = ds0[var].values.squeeze() *32/1000 # mmol/m3 to mg/L
else:
    vn0 = ds0[var].values.squeeze()
lonr0 = ds0.lon_rho.values
latr0 = ds0.lat_rho.values

# LiveOcean
fn1 = '/Users/jilian/Desktop/LiveOcean/LO_roms/cas7_t0_x4b/f2023.07.04/ocean_his_0002.nc'
ds1 = xr.open_dataset(fn1)
if var == 'oxygen':
    vn1 = ds1[var].values.squeeze() *32/1000 # # mmol/m3 to mg/L
else:
    vn1 = ds1[var].values.squeeze()
lonr1 = ds1.lon_rho.values
latr1 = ds1.lat_rho.values

#%% plot
fig = plt.figure(figsize=(12,12))
plt.rcParams['font.size'] = 12
vmin = pinfo.vlims_dict[var][0]; vmax=pinfo.vlims_dict[var][1]
vmin_diff = pinfo.vlims_dict_diff[var][0]; vmax_diff=pinfo.vlims_dict_diff[var][1]

# surface layer
ax = fig.add_subplot(231)
plt.pcolormesh(lonr0, latr0, vn0[-1,:,:], vmin=vmin, vmax=vmax, cmap='jet')
dar(ax)
plt.title('Sandbox'); 
plt.text(-0.6, 0.5, 'Surface\nlayer', c='r', transform=ax.transAxes)
plt.colorbar()

ax = fig.add_subplot(232)
plt.pcolormesh(lonr1, latr1, vn1[-1,:,:], vmin=vmin, vmax=vmax, cmap='jet')
dar(ax)
plt.colorbar() 
plt.title('LiveOcean')

ax = fig.add_subplot(233)
plt.pcolormesh(lonr1, latr1, vn0[-1,:,:] - vn1[-1,:,:],cmap='RdBu', 
               vmin=vmin_diff, vmax=vmax_diff)
plt.colorbar()
dar(ax); plt.title('Sandbox - LiveOcean')

# bottom layer
ax = fig.add_subplot(234)
plt.pcolormesh(lonr0, latr0, vn0[0,:,:], vmin=vmin, vmax=vmax, cmap='jet')
dar(ax); 
plt.title('Sandbox')
plt.colorbar()
plt.text(-0.6, 0.5, 'Bottom\nlayer', c='r', transform=ax.transAxes)

ax = fig.add_subplot(235)
plt.pcolormesh(lonr1, latr1, vn1[0,:,:], vmin=vmin, vmax=vmax, cmap='jet')
dar(ax)
plt.title('LiveOcean')
plt.colorbar(); 

ax = fig.add_subplot(236)
plt.pcolormesh(lonr1, latr1, vn0[0,:,:] - vn1[0,:,:],cmap='RdBu', 
               vmin=vmin_diff, vmax=vmax_diff)
dar(ax)
plt.title('Sandbox - LiveOcean')
plt.colorbar()


dt_str = pd.to_datetime(ds0.ocean_time.values).strftime('%Y-%m-%d %H:%M:%S').item()
plt.suptitle(var + pinfo.units_dict[var] + ', ' + dt_str)
