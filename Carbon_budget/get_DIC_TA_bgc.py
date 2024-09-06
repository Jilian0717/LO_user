# calculate DIC & TA production and comsumption from bgc process
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from lo_tools import zfun, zrfun, Lfun
from datetime import datetime, timedelta
import scipy, pickle, sys
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
lonr = fn0.lon_rho.values
latr = fn0.lat_rho.values
area = dx * dy
NX, NY = dx.shape

# We increased AttSW from 0.05 m-1 to 0.15 m-1 inside Salish Sea
AttSW = np.zeros((NX, NY)) + 0.05
AttSW[(lonr > -123.89) & (latr < 50.29) & (latr > 47.02)] = 0.15
AttSW[(lonr > -125.31) & (lonr < -123.89) & (latr < 51.02) & (latr > 49.13)] = 0.15

dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
dt1 = datetime.strptime(ds1, Lfun.ds_fmt)
dt00 = dt0

#% load salish sea index
seg_name = 'seg_info_dict_cas7_c2_noriv.p'
seg_df = pd.read_pickle(seg_name)
ji_list = seg_df['sog6_m']['ji_list']
jj = [x[0] for x in ji_list]
ii = [x[1] for x in ji_list]

#rho0 = 1025.0 # ROMS/Modules/mod_scalars.F
#sec2day = 1.0/86400
#dt = 3600 # sec, hourly history file
#BioIter = 1
#dtdays = dt * sec2day / BioIter
dtdays = 3600 * (1.0/86400) / 1
#cff1 = rho0 * 550.0
#cff1 = 1025.0 * 550.0 

# light attenuation
#PARfrac = 0.43
#Cp = 3985 # Joules/kg/degC
#AttSW = 0.05
#AttChl = 0.012
#Vp = 1.7;
#PhyIS = 0.07 # initial slope of P-I curve [1/(Watts m-2 day)]
#rOxNO3 = 138/16
#rOxNH4 = 106/16
#K_NO3 = 10 # [1/(millimole_N m-3)]
#K_NH4 = 10
#SDeRRN = 0.1 # Small detritus remineralization rate N-fraction [1/day]
#LDeRRN = 0.1
#NitriR = 0.05 # Nitrification rate: oxidation of NH4 to NO3 [1/day]

#Ws_L = 80 # m/d
#Ws_S = 8 # m/d

# initialization
# these variables were integrated over the Salish domain
DIC_vol_sum   = []  # DIC * vol
DIC_photo_sum = []  # DIC consumed by photosynthesis
DIC_remi_sum  = []  # DIC production in water column remineralization
DIC_sed_sum   = []  # DIC production in sediments
TA_vol_sum    = []  # TA * vol
TA_photo_sum  = []  # TA by photosynthesis: TA increase via uptake of NO3, TA decrease via                                               uptake of NH4
TA_remi_sum   = []  # TA production in water column remineralization
TA_nitri_sum  = []  # 2TA decrease during nitrification: NH4 --> NO3
TA_sed_sum    = []  # TA production in sediments

t = []
cnt = 0

while dt00 <= dt1:  # loop each day and every history file
    print(dt00)
    sys.stdout.flush()
    ds00 = dt00.strftime(Lfun.ds_fmt)
    fn_list = Lfun.get_fn_list('hourly', Ldir, ds00, ds00)
    for fn in fn_list[0:-1]: 
        ds    = xr.open_dataset(fn)
        swrad = ds.swrad.values.squeeze()   # w/m2
        chl   = ds.chlorophyll.values.squeeze()
        zeta  = ds.zeta.values.squeeze()
        h     = ds.h.values
        zw   = zrfun.get_z(h, zeta, S, only_w=True)
        dz   = np.diff(zw, axis=0)
        vol  = dz * area # volume of each grid cell
        salt = ds.salt.values.squeeze()
        NH4  = ds.NH4.values.squeeze()
        NO3  = ds.NO3.values.squeeze()
        phy  = ds.phytoplankton.values.squeeze()
        Oxy  = ds.oxygen.values.squeeze()
        SDeN = ds.SdetritusN.values.squeeze()
        LDeN = ds.LdetritusN.values.squeeze()
        
        DIC  = ds.TIC.values.squeeze()
        SDeC = ds.SdetritusC.values.squeeze()
        LDeC = ds.LdetritusC.values.squeeze()
        TA   = ds.alkalinity.values.squeeze()

        # PARsur = PARfrac * swrad #* rho0 * Cp # surface PAR, watts/m2
        Att = np.zeros(salt.shape)
        # ExpAtt = np.zeros(Att.shape)
        nk,ni,nj = Att.shape
        PAR = np.zeros([nk+1, ni, nj])
        PAR[-1,:,:] = 0.43 * swrad   # PAR[-1,:,:] = PARsur
        
        DIC_photo = np.zeros(Att.shape)
        DIC_remi  = np.zeros(Att.shape)
        TA_photo  = np.zeros(Att.shape)
        TA_remi   = np.zeros(Att.shape)
        TA_nitri  = np.zeros(Att.shape)
                
        #if np.nanmin(PARsur) > 0: # doing photosysnthesis
        if np.nanmin(0.43*swrad) > 0: # doing photosysnthesis
            for k in np.arange(29,-1,-1): # surface layer to bottom layer
                # O2 production from photosynthesis, mmol O2/m3/hr ?
                # dz = zw[k+1,:,:] - zw[k,:,:]
                # Att[k,:,:] = (AttSW + AttChl*chl[k,:,:] - 0.0065*(salt[k,:,:]-32))*dz
                Att[k,:,:] = (AttSW + 0.012*chl[k,:,:] - 0.0065*(salt[k,:,:]-32))* (zw[k+1,:,:] - zw[k,:,:])                       
                # ExpAtt[k,:,:] = np.exp(-Att[k,:,:])
                # Itop = PAR[k+1,:,:]
                # PAR[k,:,:] = Itop * (1-ExpAtt[k,:,:])/Att[k,:,:]
                PAR[k,:,:] = PAR[k+1,:,:] * (1-np.exp(-Att[k,:,:]))/Att[k,:,:] # average at cell center
                fac1 = PAR[k,:,:] * 0.07  # fac1 = PAR[k,:,:] * PhyIS
                Epp  = 1.7/np.sqrt(1.7*1.7 + fac1*fac1) # Epp = Vp/np.sqrt(Vp*Vp + fac1*fac1)
                t_PPmax = Epp * fac1
                cff1 = NH4[k,:,:] * 10; cff1[cff1<0] = 0 #cff1 = NH4[k,:,:] * K_NH4;
                cff2 = NO3[k,:,:] * 10; cff2[cff2<0] = 0 #cff2 = NO3[k,:,:] * K_NO3;
                inhNH4 = 1.0/(1.0+cff1)

                # fac1 = dtdays * t_PPmax
                # cff4 = fac1 * K_NO3 * inhNH4/(1.0 + cff2 + 2.0*np.sqrt(cff2)) * phy[k,:,:]
                cff4 = dtdays * t_PPmax * 10 * inhNH4/(1.0 + cff2 + 2.0*np.sqrt(cff2)) * phy[k,:,:]
                # cff5 = fac1 * K_NH4 / (1.0 + cff1 + 2.0 * np.sqrt(cff1)) * phy[k,:,:]
                cff5 = dtdays * t_PPmax * 10 / (1.0 + cff1 + 2.0 * np.sqrt(cff1)) * phy[k,:,:]
                # N_Flux_NewProd = NO3[k,:,:] * cff4
                # N_Flux_RegProd = NH4[k,:,:] * cff5
                # Oxy_pro[k,:,:] = N_Flux_NewProd * rOxNO3 + N_Flux_RegProd * rOxNH4
                # Oxy_pro[k,:,:] = NO3[k,:,:] * cff4 * 138/16 + NH4[k,:,:] * cff5 * 106/16
                
                # O2 comsumption by nitrification in water column, mmol O2/m3/hr ?
                fac2 = Oxy[k,:,:];    fac2[fac2<0] = 0
                fac3 = fac2/(3+fac2); fac3[fac3<0] = 0
                # fac1 = dtdays * NitriR * fac3
                # cff3 = fac1
                # N_Flux_Nitrifi = NH4[k,:,:] * cff3
                # Oxy_nitri[k,:,:] = 2.0 * N_Flux_Nitrifi
                # Oxy_nitri[k,:,:] = 2.0 * NH4[k,:,:] * dtdays * 0.05 * fac3
              
                # PAR[k,:,:] = Itop * ExpAtt[k,:,:]  # !!! light attenuation at the bottom of grid cell !!!!
                PAR[k,:,:] = PAR[k+1,:,:] * np.exp(-Att[k,:,:])  
                
                # DIC uptake during phytoplankton growth
                # PhyCN = 6.625  # mole_C/mole_N
                DIC_photo[k,:,:] = 6.625*(NO3[k,:,:] * cff4 + NH4[k,:,:] * cff5)
                
                # TA: N_Flux_NewProd - N_Flux_RegProd
                TA_photo[k,:,:] = NO3[k,:,:] * cff4 - NH4[k,:,:] * cff5
                # 2 * N_Flux_Nitrification
                TA_nitri[k,:,:] = 2*NH4[k,:,:] * dtdays * 0.05 * fac3
        
        else: # PARsur = 0, nitrification occurs at the maximum rate (NitriR)
            # cff3 = dtdays * NitriR
            # cff3 = dtdays * 0.05
            for k in np.arange(29,-1,-1):
                # N_Flux_Nitrifi = NH4[k,:,:] * cff3
                # Oxy_nitri[k,:,:] = 2.0 * N_Flux_Nitrifi
                # Oxy_nitri[k,:,:] = 2.0 * NH4[k,:,:] * dtdays * 0.05
                TA_nitri[k,:,:] = 2.0 * NH4[k,:,:] * dtdays * 0.05
            
        #DIC production by remineralization in water column, mmol C/m3/hr
        for k in np.arange(0,30):
            # fac1 = 0; fac2 = 1
            # cff1 = dtdays * SDeRRN * fac2
            # cff3 = dtdays * LDeRRN * fac2
            # N_Flux_RemineS = SDeN[k,:,:]*cff1
            # N_Flux_RemineL = LDeN[k,:,:]*cff3
            # Oxy_remi[k,:,:] = (N_Flux_RemineS + N_Flux_RemineL)*rOxNH4
            # Oxy_remi[k,:,:] = (SDeN[k,:,:]*dtdays*0.1*1 + LDeN[k,:,:]*dtdays*0.1*1) * 106/16
            DIC_remi[k,:,:] = SDeC[k,:,:]*dtdays*0.1*1 + LDeC[k,:,:]*dtdays*0.1*1
            # N_Flux_RemieS + N_Flux_RemineL
            TA_remi[k,:,:] = SDeN[k,:,:]*dtdays*0.1*1 + LDeN[k,:,:]*dtdays*0.1*1

        # mmol C/m3/hr to mmol C/hr
        # only save Salish Sea
        DIC_photo = DIC_photo * vol
        DIC_photo_sum.append(np.nansum(DIC_photo[:,jj,ii]))
        DIC_remi = DIC_remi * vol
        DIC_remi_sum.append(np.nansum(DIC_remi[:,jj,ii]))
        DIC_vol_sum.append(np.nansum(DIC[:,jj,ii]*vol[:,jj,ii]))
        
        TA_photo = TA_photo * vol
        TA_photo_sum.append(np.nansum(TA_photo[:,jj,ii]))
        TA_remi = TA_remi * vol
        TA_remi_sum.append(np.nansum(TA_remi[:,jj,ii]))
        TA_nitri = TA_nitri * vol
        TA_nitri_sum.append(np.nansum(TA_nitri[:,jj,ii]))
        TA_vol_sum.append(np.nansum(TA[:,jj,ii]*vol[:,jj,ii]))
        
        #---------- sediment DIC production ----------
        # DIC production in sediments from SDeC, LDeC
        DIC_sed_tmp = (SDeC[0,:,:] * 8)/24  + (LDeC[0,:,:] * 80)/24 # mmol C/m2/hr
        DIC_sed_sum.append(np.nansum(DIC_sed_tmp[jj,ii] * area[jj,ii]))  # mmol/hr
        
        # ---------- sediment TA production -----------
        TA_sed_tmp = (LDeN[0,:,:] * 80)/24 + (SDeN[0,:,:] * 8)/24  # mmol C/m2/hr
        TA_sed_sum.append(np.nansum(TA_sed_tmp[jj,ii] * area[jj,ii]))
        
        cnt += 1
        t.append(ds.ocean_time.values)
        ds.close()       
    dt00 = dt00 + timedelta(days=1)
        
# save netcdf
from netCDF4 import Dataset
nc   = Dataset('DIC_TA_bgc_'+ds0+'.nc','w')
time = nc.createDimension('time', len(t))
eta_rho = nc.createDimension('eta_rho', NX)
xi_rho  = nc.createDimension('xi_rho', NY)
s_rho   = nc.createDimension('s_rho', 30)

times       = nc.createVariable('time','f8',('time',))
times.units = 'seconds*1e9 since 1970-01-01 00:00:00'
DIC_photo_sum_tmp       = nc.createVariable('DIC_photo_sum','f4',('time',),compression='zlib',complevel=9)
DIC_photo_sum_tmp.units = 'mmol C/hr'
DIC_remi_sum_tmp        = nc.createVariable('DIC_remi_sum','f4',('time',),compression='zlib',complevel=9)
DIC_remi_sum_tmp.units  = 'mmol C/hr'
DIC_sed_sum_tmp         = nc.createVariable('DIC_sed_sum','f4',('time',),compression='zlib',complevel=9)
DIC_sed_sum_tmp.units   = 'mmol C/hr'
DIC_vol_sum_tmp         = nc.createVariable('DIC_vol_sum','f4', ('time',),compression='zlib',complevel=9)
DIC_vol_sum_tmp.units   = 'mmol C'

TA_photo_sum_tmp       = nc.createVariable('TA_photo_sum','f4',('time',),compression='zlib',complevel=9)
TA_photo_sum_tmp.units = 'mmol C/hr'
TA_remi_sum_tmp        = nc.createVariable('TA_remi_sum','f4',('time',),compression='zlib',complevel=9)
TA_remi_sum_tmp.units  = 'mmol C/hr'
TA_nitri_sum_tmp       = nc.createVariable('TA_nitri_sum','f4',('time',),compression='zlib',complevel=9)
TA_nitri_sum_tmp.units = 'mmol C/hr'
TA_sed_sum_tmp         = nc.createVariable('TA_sed_sum','f4',('time',),compression='zlib',complevel=9)
TA_sed_sum_tmp.units   = 'mmol C/hr'
TA_vol_sum_tmp         = nc.createVariable('TA_vol_sum','f4', ('time',),compression='zlib',complevel=9)
TA_vol_sum_tmp.units   = 'mmol C'

#
times[:] = t
DIC_photo_sum_tmp[:] = DIC_photo_sum
DIC_remi_sum_tmp[:]  = DIC_remi_sum
DIC_sed_sum_tmp[:]   = DIC_sed_sum
DIC_vol_sum_tmp[:]   = DIC_vol_sum
TA_photo_sum_tmp[:]  = TA_photo_sum
TA_remi_sum_tmp[:]   = TA_remi_sum
TA_nitri_sum_tmp[:]  = TA_nitri_sum
TA_sed_sum_tmp[:]    = TA_sed_sum
TA_vol_sum_tmp[:]    = TA_vol_sum

nc.close()
