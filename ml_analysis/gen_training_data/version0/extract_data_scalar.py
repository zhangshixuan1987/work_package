import netCDF4
import sys
import os
import xarray as xr 
import numpy as np
from netCDF4 import Dataset as NetCDFFile
from netCDF4 import num2date,date2num
import datetime
import pandas as pd
import pylab as plt
import matplotlib.cm as cm
from matplotlib import colors

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

#process data #
var_list = ["LANDFRAC", "OCNFRAC", "ICEFRAC"]
var_outs = ["LANDFRAC", "OCNFRAC", "ICEFRAC"]

#note: we use this simulation to extract the time information 
#note: to keep the consistenty time format used for other quantity 
exp = "E3SMv2_NDGUVTQ_SRF1_tau6"
datadir = os.path.join(topdir,"New_Training","nc_data",exp)
dataref = os.path.join(topdir,"SE_PG2","before_nudging")
smdh = "01010000"
emdh = "12312100"
for year in range(syear,eyear+1):
  fref    = "{}_3hourly_{:04d}{}-{:04d}{}.nc".format(exp,year,smdh,year,emdh)
  fr      = xr.open_dataset(os.path.join(dataref,fref),decode_times=True,decode_timedelta=True) 
  lat_val = fr['lat'] 
  lon_val = fr['lon'] 
  tim_val = fr['time'].astype('datetime64[ns]') #[:].to_pandas().astype('datetime64[ns]') #to_pandas()
  print('latitude range:',  lat_val[:].min(),lat_val[:].max())
  print('longitude range:', lon_val[:].min(),lon_val[:].max())
  print('time range:', tim_val[0],tim_val[len(tim_val)-1])

  ax1 = len(tim_val) # itim2 - itim1 + 1
  ax2 = len(lat_val) # ilat2 - ilat1 + 1
  ax3 = len(lon_val) # ilon2 - ilon1 + 1

  for ii in range(len(var_list)):
    var      = var_list[ii]
    vou      = var_outs[ii]
    fname    = "{}_Scalar_monthly_{}.nc".format(exp,year)
    df       = xr.open_dataset(os.path.join(datadir,fname),decode_times=True,decode_timedelta=True)
    data     = df[var][:,:]  
    str_time = df['time_bnds'][:,0].astype('datetime64[ns]')
    end_time = df['time_bnds'][:,1].astype('datetime64[ns]')

    data_np = np.zeros(shape=(ax1, ax2))
    for jj in range(len(str_time)):
      mask = (tim_val >= str_time[jj]) & (tim_val < end_time[jj])
      data_np[mask,:] = data[jj,:]
      #print("monthly:", data[jj,:].data)
      #print("3hourly:", data_np[mask,:])
    print(np.min(data[:,:].data),np.max(data[:,:].data))
    print(np.min(data_np[:,:]),np.max(data_np[:,:]))

#    #sanity check
#    #plt.figure()
#    fig, axs = plt.subplots(2,1, figsize=(9, 6))
#    fig.suptitle('{}'.format(var))
#    # create a single norm to be shared across all images
#    #norm = colors.Normalize(vmin=0.0, vmax=1.0)
#    #images = []
#    #images.append(axs[0].imshow(np.array(data.data),cmap=cm.RdYlGn)) #, norm=norm))
#    #images.append(axs[1].imshow(data.data,cmap=cm.RdYlGn)) #, norm=norm))
#    cntr1 = axs[0].contourf(data.data, levels=14, cmap="RdBu_r")
#    cntr2 = axs[1].contourf(data_np, levels=14, cmap="RdBu_r")
#    #fig.colorbar(cntr1[0], ax=axs, orientation='horizontal', fraction=.1)
#    plt.show()
#    #plt.savefig(os.path.join(fig_path, var + ".png"))
#    #plt.close()

    outname = "{}_3hourly_{:04d}.npy".format(vou,year)
    fout    = os.path.join(outpath,outname)
    if os.path.exists(fout):
      os.remove(fout)
    np.save(fout, data_np, allow_pickle=True,fix_imports=True)
    del(outname,fout,data_np)
  fr.close()
  df.close()
  del(ax1,ax2,ax3,fr,df,fname,lat_val,lon_val,tim_val)  
del(datadir,var_list,var_outs)

