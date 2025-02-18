"""
calculate bottom pressure anomaly
refer to: https://github.com/parkermac/LPM/tree/main/bpress
Qs: model is smoother than obs
"""
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from datetime import datetime, timedelta
from mat4py import loadmat
from lo_tools import zfun
import pandas as pd
import gsw

sn_list = ['CP01', 'CP015', 'CP02', 'CP03']
fig = plt.figure(figsize=(14,10))
plt.rcParams['font.size'] = 12
ii = 1

LO_cas7 = False
LO_cas6 = False

for sn in sn_list:
    # data from Matt
    data = loadmat('./Matt_processed_data/'+sn+'.mat')
    data = data[sn]
    data_lp = loadmat('./Matt_processed_data/'+sn+'lp.mat')
    data_lp = data_lp[sn+'lp']

    # pressure
    prs   = np.array(data['prs']) #detided, de-drifted pressure
    prsdd = np.array(data['prsdd']).squeeze() # decimal days from 1/1/2017
    prs_lp   = np.array(data_lp['prs'])
    prsdd_lp = np.array(data_lp['prsdd']).squeeze() 

    prs_time    = []
    prs_time_lp = []
    for i in range(len(prsdd)):
        prs_time.append(datetime(2017,1,1) + timedelta(days=prsdd[i])); 

    for i in range(len(prsdd_lp)):
        prs_time_lp.append(datetime(2017,1,1) + timedelta(days=prsdd_lp[i])); 

    #%% prs_lp from Matt is very close to applying a 3-day hanning window to prs
    #fig = plt.figure(figsize=(16,10))
    #ax = fig.add_subplot(411)
    #plt.plot(prs_time, prs*100, 'k', lw=0.1);
    #plt.plot(prs_time, zfun.lowpass(prs, f='hanning', n=24*3)*100, c='r', label='3-day hanning (jx)')
    #plt.ylabel('hPa');
    #plt.title('De-tide, de-drifted pressure');
    # dB to hPa
    #plt.plot(prs_time_lp, prs_lp*100, 'b', lw=1, label='lp-matt')
    #plt.legend()
    #%% measured velocity
    u = np.array(data['u']).squeeze()
    v = np.array(data['v']).squeeze()
    uvdd = np.array(data['uvdd']).squeeze()
    u_lp = np.array(data_lp['u']).squeeze()
    v_lp = np.array(data_lp['v']).squeeze()
    uvdd_lp = np.array(data_lp['uvdd']).squeeze()
    uvdd_time = []
    for i in range(len(uvdd)):
        uvdd_time.append(datetime(2017,1,1) + timedelta(days=uvdd[i])); 
    uvdd_time_lp = []
    for i in range(len(uvdd_lp)):
        uvdd_time_lp.append(datetime(2017,1,1) + timedelta(days=uvdd_lp[i])); 

    #fig = plt.figure(figsize=(16,10))
    #ax = fig.add_subplot(211)
    #plt.plot(uvdd_time, u, 'k', lw=0.1);
    #plt.plot(uvdd_time, zfun.lowpass(u, f='hanning', n=24*3), c='r', label='3-day hanning (jx)')
    #plt.ylabel('u');
    ## dB to hPa
    #plt.plot(uvdd_time_lp, u_lp, 'b', lw=1, label='lp-matt')
    #plt.legend()    
    #%% load hourly model extraction
    if LO_cas7:
        dir0 = '/Users/jilian/Desktop/LiveOcean/LO_output/extract/cas7_t0_x4b/moor/matt_wei/'
        ds = xr.open_dataset(dir0+sn+'_2017.04.01_2017.12.31.nc') #hourly extraction
    elif LO_cas6:
        dir0 = '/Users/jilian/Desktop/LiveOcean/LO_output/extract/cas6_v0_live/moor/matt_wei/'
        ds = xr.open_dataset(dir0+sn+'_2017.04.01_2017.12.31.nc') #hourly extraction
    else:
        dir0='/Users/jilian/Desktop/LiveOcean/LO_output/extract/or2_t0_x0/moor/matt_wei/'
        ds=xr.open_dataset(dir0+sn+'_2017.04.01_2017.12.31.nc')

    # get time axis
    ot = ds.ocean_time.values # an array with dtype='datetime64[ns]'
    dti = pd.to_datetime(ot) # a pandas DatetimeIndex with dtype='datetime64[ns]'
    dt = dti.to_pydatetime() # an array of datetimes
    
    # set constants
    pad = 36 # this trims the ends after the low pass so there are no nan's
    g = 9.81 # gravity [m s-2]
    
    # pull fields from model dataset
    z  = ds.z_rho.values
    zw = ds.z_w.values
    eta  = ds.zeta.values
    salt = ds.salt.values
    temp = ds.temp.values
    lon = ds.lon_rho.values
    lat = ds.lat_rho.values
    NT, N = z.shape
    
    #%%
    # make time-mean z positions so we can isolate baroclinic and SSH
    # contributions to pressure variation
    Z  = np.mean(z, axis=0)
    ZW = np.mean(zw, axis=0)
    DZ = np.diff(ZW)
    # adjust so free surface is at 0
    Z  -= ZW[-1]
    ZW -= ZW[-1]
    
    # Equation of state calculations
    p  = gsw.p_from_z(Z, lat)
    SA = gsw.SA_from_SP(salt, p, lon, lat) # absolute salinity
    CT = gsw.CT_from_pt(SA, temp) # conservative temperature
    if True:
        rho = gsw.rho(SA, CT, p)
        # This is denser than ROMS rho by 0.037 [kg m-3] at the bottom and 0.0046 [kg m-3]
        # (annual averages), and it is the full density, not density minus 1000.
        # There was no visual difference between the pressure time series.
    else:
        rho = ds.rho.values # makes no difference to add 1000 [kg m-3]
            
    # low pass filtered version
    etalp  = zfun.lowpass(eta, f='godin')[pad:-pad:24]
    etalp  = etalp  - np.mean(etalp) # remove mean SSH
    rholp  = zfun.lowpass(rho, f='godin')[pad:-pad:24, :]
    #saltlp = zfun.lowpass(SA, f='godin')[pad:-pad:24, :]
    #templp = zfun.lowpass(CT, f='godin')[pad:-pad:24, :]
    
    # also make associated time vectors
    tlp  = dt[pad:-pad:24]
    tlpf = dt[24:-24:24] # "full" version used for pcolormesh plots
    NTlp = len(etalp)
    # calculate the baroclinic pressure
    plp = np.flip(np.cumsum(np.flip(g * rholp * DZ.reshape((1,N)), axis=1), axis=1), axis=1)

    # calculate the pressure due to SSH
    plp0 = g * 1025 *etalp
    
    # annual means
    #Rho  = np.mean(rholp, axis=0)
    #Salt = np.mean(saltlp, axis=0)
    #Temp = np.mean(templp, axis=0)
    Pm    = np.mean(plp, axis=0)

    # anomalies from the annual mean
    #saltlp_a = saltlp - Salt.reshape((1,N))
    #templp_a = templp - Temp.reshape((1,N))
    #rholp_a  = rholp  - Rho.reshape((1,N))
    plp_a    = plp    - Pm.reshape((1,N))

    # separate contributions of salt and temp to density
    #rho_only_salt = gsw.rho(saltlp, templp - templp_a, p) - Rho.reshape((1,N))
    #rho_only_temp = gsw.rho(saltlp - saltlp_a, templp, p) - Rho.reshape((1,N))

    # make the full pressure anomaly
    plp_aa = plp_a + plp0.reshape((NTlp,1))
    
    ax = fig.add_subplot(4,2,ii)
    nfilt = 1 # days for Hanning filter
    P = zfun.lowpass(plp_aa[:,0]/100, 'hanning', n=nfilt) # Pa to hPa, bottom
    ax.plot(tlp, P, '-b', lw=2, label='Model')
    ax.set_ylim(-7,7)
    ax.set_xlim(datetime(2017,4,1),datetime(2017,12,1))
    ax.axhline(c='k',lw=1.5)
    #ax.grid(True)
    Depth = int(-ZW[0])
    
    # observation
    ax.plot(prs_time_lp, prs_lp*1e2, 'r', lw=2, label='obs') # dB to hPa
 #   ax.plot(prs_time, zfun.lowpass(prs*1e2, f='hanning',n=24*10), 'g', lw=1, label='obs') # dB to hPa
    
    #%% jx's version of total bottom pressure anomaly, close to Parker's method
    add_jx_cal = False
    if add_jx_cal:
        dz_jx = np.diff(zw, axis=1)
        p_jx  = gsw.p_from_z(z, lat)
        SA_jx = gsw.SA_from_SP(salt, p_jx, lon, lat) # absolute salinity
        CT_jx = gsw.CT_from_pt(SA_jx, temp) # conservative temperature
        rho_jx = gsw.rho(SA_jx, CT_jx, p_jx)
        
        ptot = np.flip(np.cumsum(np.flip(g * rho_jx * dz_jx, axis=1), axis=1), axis=1)
        ptot_m = np.mean(ptot, axis=0)
        ptot_a = ptot - ptot_m.reshape((1,30))
        ax.plot(dt, zfun.lowpass(ptot_a[:,0]/100, f='godin'), c='g')
    #%%
    ax.text(.01, .87, '%s (Bottom Depth = %d m)' % (sn, Depth), transform=ax.transAxes, weight='bold')
    plt.ylabel('hPa')
    if ii == 1:
        if LO_cas7:
            plt.suptitle('LiveOcean - cas7', fontsize=15, weight='bold')
        elif LO_cas6:
            plt.suptitle('LiveOcean - cas6', fontsize=15, weight='bold')
        else:
            plt.suptitle('Small domain model', fontsize=15, weight='bold')
        ax.legend(ncols=2)
    if ii == 7:
        ax.set_xlabel('Date')
    else:
        ax.set_xticklabels([])

    #%% velocity
    u_mod = ds.u.values
    v_mod = ds.v.values
    u_mod_lp = zfun.lowpass(u_mod, f='godin')[pad:-pad:24, :]
    u_mod_lp = zfun.lowpass(u_mod_lp, f='hanning', n=nfilt)
    v_mod_lp = zfun.lowpass(v_mod, f='godin')[pad:-pad:24, :]
    v_mod_lp = zfun.lowpass(v_mod_lp, f='hanning', n=nfilt)

    reference_date = pd.to_datetime('2017-01-01')
    tt = (dti[pad:-pad:24]-reference_date)/np.timedelta64(1, 'D')

    ax = fig.add_subplot(4,2,ii+1)
    plt.quiver(uvdd_lp, np.zeros(len(uvdd_lp)), u_lp, v_lp, color='r', 
               scale_units='y', scale=1, width=0.002,)
    plt.quiver(tt, np.zeros(len(tt)), u_mod_lp[:,0]*100, v_mod_lp[:,0]*100, color='b', 
               scale_units='y', scale=1, width=0.002,)
    plt.ylim([-8, 8])
    ax.axhline(y=0, c='k',lw=1.5)
    ax.text(.01, .87,  sn, transform=ax.transAxes, weight='bold')
    if ii == 7:
        ax.set_xlabel('Days since 2017.01.01')
    else:
        ax.set_xticklabels([])
    plt.ylabel('cm/s')    
    ii += 2

plt.tight_layout(h_pad=0.5, w_pad=0.1)
if LO_cas7:
    plt.savefig('LiveOcean-cas7',dpi=300,bbox_inches='tight')
elif LO_cas6:
    plt.savefig('LiveOcean-cas6',dpi=300,bbox_inches='tight')
else:
    plt.savefig('or2',dpi=300,bbox_inches='tight')


