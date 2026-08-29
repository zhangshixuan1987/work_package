import netCDF4
import sys
import os 
import numpy as np
from netCDF4 import Dataset as NetCDFFile
from netCDF4 import num2date,date2num
from netCDF4 import Dataset,num2date,date2num
import datetime
import pandas as pd
from datetime import datetime, timedelta

#Notes: we use AMIP simulation where the sea ice fraction 
#       and surface temperature are prescribed. Therefore 
#       variables like the land, ocean and sea ice fraction 
#       will not depend on specific model configurations for 
#       nuding, and all simulations shoud share these data 

topdir = "/global/cfs/cdirs/e3sm/www/zhan391/darpa_temporary_data_share"
outdir = "/global/cfs/cdirs/e3sm/www/zhan391/SEA_CROGS"
syear  = 2007
eyear  = 2017
tunit  = "hours since 2007-01-01 00:00:00"

#constant values
deg2rad=np.pi/180.  # degree to radians
rad2deg=180./np.pi  # radians to degree

#output path
outpath  = os.path.join(outdir,"New_Training","E3SMv2_Scalar")
if not os.path.exists(outpath):
  os.makedirs(outpath)

var_list = ["area","time", "lat",      "lon"      ]
var_outs = ["area","time", "latitude", "longitude"]

#note: we use this simulation to extract the time information 
#note: to keep the consistenty time format used for other quantity 
exp = "E3SMv2_NDGUVTQ_SRF1_tau6"
datadir = os.path.join(topdir,"SE_PG2","before_nudging")
smdh = "01010000"
emdh = "12312100"
for year in range(syear,eyear+1):
  fname = "{}_3hourly_{:04d}{}-{:04d}{}.nc".format(exp,year,smdh,year,emdh)
  f = NetCDFFile(os.path.join(datadir,fname),'r')
  lat_val = f.variables['lat']
  lon_val = f.variables['lon']
  area_val= f.variables['area']
  tim_val = f.variables['time']
  tim_old = num2date(tim_val[:], tim_val.units, tim_val.calendar)
  tim_new = date2num(tim_old[:], tunit, tim_val.calendar)

  print('latitude range:',  lat_val[:].min(),lat_val[:].max())
  print('longitude range:', lon_val[:].min(),lon_val[:].max())
  print('time range:', tim_new[:].min(),tim_new[:].max())

  ax1 = len(tim_val) # itim2 - itim1 + 1
  ax2 = len(lat_val) # ilat2 - ilat1 + 1
  ax3 = len(lon_val) # ilon2 - ilon1 + 1
  
  for ii in range(len(var_list)):
    var = var_list[ii]
    vou = var_outs[ii]
    print("vars are :{}".format(var))

    data_np = np.zeros(shape=(ax1, ax2))
    if var == "time":
      data_np = np.zeros(shape=(ax1))  
      data_np[:] = tim_new[:]
    else:
      data_np = np.zeros(shape=(ax2))
      if var == "lat":  
        data_np[:] = lat_val[:]  
      elif var == "lon":
        data_np[:] = lon_val[:] 
      else:
        data_np[:] = area_val[:]

    data_np = np.array(data_np)
    outname = "{}_3hourly_{:04d}.npy".format(vou,year)
    fout    = os.path.join(outpath,outname)
    if os.path.exists(fout):
      os.remove(fout)
    np.save(fout,data_np,allow_pickle=True,fix_imports=True)
    del(data_np,outname,fout)
  f.close()  
  del(f,fname,ax1,ax2,ax3,lat_val,lon_val,tim_val,tim_old,tim_new)
del(datadir,var_list,var_outs)
