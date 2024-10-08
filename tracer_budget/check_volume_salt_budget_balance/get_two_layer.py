"""
 From /LO/extract/tef2/tef_fun.py
"""
import xarray as xr
import numpy as np
import pandas as pd

#def get_two_layer(in_dir, sect_name):

def get_two_layer(bulk):
    """
    Form time series of 2-layer TEF quantities, from the multi-layer bulk values.
    """
   # bulk = xr.open_dataset(in_dir / (sect_name + '.nc'))
    
    # determine which variables to process
    vn_list = []
    vec_list = []
    for vn in bulk.data_vars:
        if ('time' in bulk[vn].coords) and ('layer' in bulk[vn].coords):
            vn_list.append(vn)
        elif ('time' in bulk[vn].coords) and ('layer'  not in bulk[vn].coords):
            vec_list.append(vn)
    vn_list.remove('q') # transport is handled separately from tracers
    
    # separate positive and negative transports
    QQ = bulk.q.values
    QQp = QQ.copy()
    QQp[QQ<=0] = np.nan
    QQm = QQ.copy()
    QQm[QQ>=0] = np.nan
    
    # form two-layer versions of volume and tracer transports
    Qp = np.nansum(QQp, axis=1)
    Qm = np.nansum(QQm, axis=1)
    if True:
        # Mask out times when the transport is too small to
        # use for tracer averaging.
        # Not needed? Seems like a possible source of budget errors.
        Qp[Qp<np.nanmean(Qp)/100] = np.nan
        Qm[Qm>np.nanmean(Qm)/100] = np.nan
    QCp = dict()
    QCm = dict()
    for vn in vn_list:
        QCp[vn] = np.nansum(QQp*(bulk[vn].values), axis=1)
        QCm[vn] = np.nansum(QQm*(bulk[vn].values), axis=1)
    
    # form flux-weighted tracer concentrations
    Cp = dict()
    Cm = dict()
    for vn in vn_list:
        Cp[vn] = QCp[vn]/Qp
        Cm[vn] = QCm[vn]/Qm
        
    # pack results in a DataFrame
    tef_df = pd.DataFrame(index=bulk.time.values)
    tef_df['q_p']=Qp
    tef_df['q_m']=Qm
    for vn in vn_list:
        tef_df[vn+'_p'] = Cp[vn]
        tef_df[vn+'_m'] = Cm[vn]
    # also pass back time series like qprism, for convenience
    for vn in vec_list:
        tef_df[vn] = bulk[vn].values
        
    bulk.close()
            
    return tef_df, vn_list, vec_list
