#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 11:18:51 2024

@author: jilian
"""
from lo_tools import zfun, zrfun, Lfun
from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt

D0 = 282.6;    D1 = 0.125;   D2 = -7.18;   D3 = 0.86
D4 = -0.99;    D5 = 0.28;    D6 = -0.80;   D7 = 0.06

pi2 = np.pi * 2

ds0 = '1951.01.01'
ds0 = '2010.01.01'
ds1 = '2020.01.01'
dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
dt1 = datetime.strptime(ds1, Lfun.ds_fmt)

dt00 = dt0

pCO2air_secular = []
t = []

while dt00 <= dt1:
    t.append(dt00)
    year = dt00.year
    yday = dt00 - datetime(year,1,1)
    yday = yday.days
    
    pmonth = year - 1951.0 + yday/365.0  # months since Jan 1951
    
    tmp = D0 + D1*pmonth*12.0 + D2*np.sin(pi2*pmonth + D3) + \
                      D4*np.sin(pi2*pmonth + D5) + D6*np.sin(pi2*pmonth + D7)
                      
    pCO2air_secular.append(tmp)
    
    dt00 = dt00 + timedelta(days=1)

plt.figure(figsize=(10,6))
plt.plot(t, pCO2air_secular)
plt.title('pCO2air_secular (cas7_t0_x4b, fennel.h)')
plt.grid()
plt.show()


