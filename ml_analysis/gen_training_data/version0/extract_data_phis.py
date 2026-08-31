import netCDF4
import sys
import os 
import numpy as np
from netCDF4 import Dataset as NetCDFFile
from netCDF4 import num2date,date2num

topdir   = "/global/cfs/cdirs/e3sm/www/zhan391/darpa_temporary_data_share"
outdir   = "/global/cfs/cdirs/e3sm/www/zhan391/SEA_CROGS"
exp      = "E3SMv2_NDGUVTQ_SRF1_tau6"

var      = "PHIS"
vou      = "PHIS"
syear    = 2007
eyear    = 2017
smdh     = "01010300"
emdh     = "01010000"
datadir  = os.path.join(topdir,"SE_PG2","nudging_tendency")

#output path
outpath  = os.path.join(outdir,"New_Training","E3SMv2_Scalar")
if not os.path.exists(outpath):
  os.makedirs(outpath)

for year in range(syear,eyear+1):
  print("vars are :{}".format(var))  
  fname = "{}_3hourly_{:04d}{}-{:04d}{}.nc".format(exp,year,smdh,year+1,emdh)
  f = NetCDFFile(os.path.join(datadir,fname),'r')

  lat_val = f.variables['lat'][:]
  lon_val = f.variables['lon'][:]
  lev_val = f.variables['lev'][:]
  tim_val = f.variables['time'][:]
  print('latitude range:',  lat_val.min(),lat_val.max())
  print('longitude range:', lon_val.min(),lon_val.max())
  print('level range:', lev_val.min(),lev_val.max())
  print('time range:', tim_val.min(),tim_val.max())
  
  ax1 = len(tim_val) # itim2 - itim1 + 1
  ax2 = len(lev_val) #
  ax3 = len(lat_val) # ilat2 - ilat1 + 1
  ax4 = len(lon_val) # ilon2 - ilon1 + 1

  data_np = np.zeros(shape=(ax1, ax3))
  data_np[:,:] = np.array(f.variables[var][:,:,:])
  
  outname = "{}_3hourly_{:04d}.npy".format(vou,year)

  fout    = os.path.join(outpath,outname)
  if os.path.exists(fout):
    os.remove(fout)

  np.save(fout, data_np, allow_pickle=True,fix_imports=True)

  f.close()

  del(data_np,ax1,ax2,ax3,ax4,f,fname)
del(smdh,emdh,var_list,var_outs)

