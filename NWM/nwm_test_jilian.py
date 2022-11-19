"""
test downloading NWM forecast data

IOOS NWM API prototype: https://prototype.ioos.us/nwm

"""
from lo_tools import plotting_functions as pfun
import matplotlib.pyplot as plt
import netCDF4 as nc4
import requests 

#%% I copied this very lon url from NWM website after I selected the variables
# and box domain mannually.   
url='https://prototype.ioos.us/dap/noaa.nwm.short_range.channel_rt/run/latest.nc?lat_lon%5B0:1:2776737%5D%5B0:1:1%5D,time_offset%5B0:1:17%5D,time_run%5B0:1:17%5D,time%5B0:1:17%5D,reference_time%5B0:1:17%5D%5B0:1:0%5D,crs,feature_id%5B0:1:2776737%5D,streamflow%5B0:1:17%5D%5B0:1:2776737%5D,nudge%5B0:1:17%5D%5B0:1:2776737%5D,velocity%5B0:1:17%5D%5B0:1:2776737%5D,qSfcLatRunoff%5B0:1:17%5D%5B0:1:2776737%5D,qBucket%5B0:1:17%5D%5B0:1:2776737%5D,qBtmVertRunoff%5B0:1:17%5D%5B0:1:2776737%5D,bbox(-132.25324630737308,40.93309804899676,-118.92425537109376,52.24041522350549)'   
r = requests.get(url)
out_fn = 'test.nc'
with open(out_fn,'wb') as f:
    f.write(r.content)   
#%% check donwloaded data
nc = nc4.Dataset(out_fn)
streamflow = nc['streamflow']
lat = nc['lat_lon'][:, 0]
lon = nc['lat_lon'][:, 1]
t = nc['time'] #'seconds since 2022-11-18 22:00:00'
ref_t = nc['reference_time'] #minutes since 1970-01-01 00:00:00 UTC
#
plt.close('all')
pfun.start_plot(figsize = (10, 13))
fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(lon, lat,'k.', markersize = 1)
pfun.add_coast(ax, color='r', linewidth = 1.5)
#%% try self-defined bbox and other choices
min_lon = -132
max_lon = -118
min_lat = 42
max_lat = 52
tmp1 = '0:1:2776737'
tmp2 = '0:1:1'
tmp3 = '0:1:17'
tmp4 = '0:1:0'
url_test = ('https://prototype.ioos.us/dap/noaa.nwm.short_range.channel_rt/run/latest.nc?'
       +'lat_lon'+'%5B'+tmp1+'%5D%5B'+tmp2+'%5D,'
       +'time_offset'+'%5B'+tmp3+'%5D,'
       +'time_run'+'%5B'+tmp3+'%5D,'
       +'time+''%5B'+tmp3+'%5D,'
       +'reference_time'+'%5B'+tmp3+'%5D%5B'+tmp4+'%5D,'
       +'crs,'
       +'feature_id'+'%5B'+tmp1+'%5D,'
       +'streamflow'+'%5B'+tmp3+'%5D%5B'+tmp1+'%5D,'
       +'nudge'+'%5B'+tmp3+'%5D%5B'+tmp1+'%5D,'
       +'velocity'+'%5B'+tmp3+'%5D%5B'+tmp1+'%5D,'
       +'qSfcLatRunoff'+'%5B'+tmp3+'%5D%5B'+tmp1+'%5D,'
       +'qBucket'+'%5B'+tmp3+'%5D%5B'+tmp1+'%5D,'
       +'qBtmVertRunoff'+'%5B'+tmp3+'%5D%5B'+tmp1+'%5D,'
       +'bbox('+str(min_lon)+','+str(min_lat)+','+str(max_lon)+','+str(max_lat)+')')

#r_test = requests.get(url_test)
#out_fn = 'test0.nc'
#with open(out_fn,'wb') as f:
#    f.write(r.content)   
