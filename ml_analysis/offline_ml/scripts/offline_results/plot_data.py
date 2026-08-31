import os 
import netCDF4 
import numpy as np
import xarray as xr
import pandas as pd
import warnings
from netCDF4 import Dataset,num2date
from scipy import stats 
from datetime import datetime

def main(exp,year,tunit,datadir,mldir,outpath):
   #var_state  = ["U","V","T","Q"]
   var_tend   = ["U_TEND","V_TEND","T_TEND","Q_TEND"]
   var_unit   = ["m/s","m/s","K","g/kg"]
   var_scale  = [1.0,1.0,1.0,1000.0]
   #read information for time, latitude and longitude 
   fpred = os.path.join(mldir,'Test_pred.npy')
   ftrue = os.path.join(mldir,'Test_true.npy')
   vpred = np.load(fpred,allow_pickle=True)
   vtrue = np.load(ftrue,allow_pickle=True)
   nvt   = vpred.shape[0]

   time,lat,lon,area,date = extract_grid_time_info(datadir,tunit)
   print("shape:", time.shape)
   print("shape:", lat.shape)
   print("shape:", lon.shape)
   print("shape:", area.shape)

   #tdate = num2date(time, tunit, calendar="noleap").astype('datetime64[ns]')
   group = "Unet_4var"
   ntime = len(time)
   ncol  = len(lat)
   nlev  = 72
   
   date = num2date(time[0:nvt-1],tunit, calendar="noleap").astype('datetime64[ns]')
   var_names = ["Nudge_U","Nudge_V","Nudge_T","Nudge_Q"]
   var_scals = [x * 86400.0 for x in var_scale]
   var_units = [unit + " day$^{-1}$" for unit in var_unit]

   dsm = construct_data_4var(var_names,var_units,var_scals,date,vpred,lat,lon,area,tunit)
   dsm.fillna(value=1e20)
   dsm.to_netcdf("model_pred_data.nc")

   dso = construct_data_4var(var_names,var_units,var_scals,date,vtrue,lat,lon,area,tunit)
   dso.fillna(value=1e20)
   dso.to_netcdf("model_true_data.nc")
   exit()

   print("plot 2d map distribution")
   print(dsm["Nudge_U"].values.min(),dsm["Nudge_U"].values.max())
   ymdh = "2008_N20"
   ilev = 48
   fgrd = xr.open_dataset(os.path.join(datadir,"E3SMv2_grid_info","ne30pg2_scrip.nc"),decode_times=False)
   plot_2d_map_4var(var_names,var_units,dsm,dso,ymdh,ilev,exp,group,fgrd)
     
   pcorr = np.corrcoef(ds0[vname].values,ds1[vname].values)
   print(pcorr.shape)
   print(pcorr.min(),pcorr.max())
   exit()
   
   #for it in range(0,ntime,1):
   #  #select a time step to plot (first step as example below)
   #  data = ds.isel(time=it,lev=48) #* vscal 
   #  ymdh = date[it] #.strftime("%Y-%m-%d_%HZ")[0] #data.time.dt.strftime("%Y-%m-%d_%HZ").data
   #  varMin, varMax, nlevs = var_min[i],var_max[i],11
   #  plot_2d_map(vname,data,ymdh,exp,group,fgrd,
   #              varMin,varMax,nlevs)
   
   return

def construct_data_4var(var_name,var_unit,var_scale,date,vdata,lat,lon,area,tunit):
   ds0 = xr.Dataset(
           data_vars={
               var_name[0]:(("time", "lev", "ncol"), 
                          np.array(vdata[:,0,:,:] * var_scale[0]),
                          {"units": var_unit[0] + " day$^{-1}$"}),
               var_name[1]:(("time", "lev", "ncol"), 
                          np.array(vdata[:,1,:,:] * var_scale[1]),
                          {"units": var_unit[1] + " day$^{-1}$"}),
               var_name[2]:(("time", "lev", "ncol"), 
                          np.array(vdata[:,2,:,:] * var_scale[2]),
                          {"units": var_unit[2] + " day$^{-1}$"}),
               var_name[3]:(("time", "lev", "ncol"), 
                          np.array(vdata[:,3,:,:] * var_scale[3]),
                          {"units": var_unit[3] + " day$^{-1}$"}),
               },
           coords={"lon": (["ncol"], lon),
                   "lat": (["ncol"], lat),
                   "area": (["ncol"],area),
                   "lev": (["lev"], np.arange(0,72,1)),
                   "reference_time": tunit,
                   },
           attrs = dict(long_name="nudging tendency(truth)"),
           )
   return ds0 

def extract_grid_time_info(datadir,tunit):
   #read information for time, latitude and longitude
   ftim = os.path.join(datadir,'E3SMv2_Scalar',"{}_3hourly_{}.npy".format("time",year))
   flat = os.path.join(datadir,'E3SMv2_Scalar',"{}_3hourly_{}.npy".format("latitude",year))
   flon = os.path.join(datadir,'E3SMv2_Scalar',"{}_3hourly_{}.npy".format("longitude",year))
   farg = os.path.join(datadir,'E3SMv2_Scalar',"{}_3hourly_{}.npy".format("area",year))
   #read and decode them to real date
   time = np.load(ftim,allow_pickle=True)
   lat  = np.load(flat,allow_pickle=True)
   lon  = np.load(flon,allow_pickle=True)
   area = np.load(farg,allow_pickle=True)
   date = num2date(time,tunit, calendar="noleap").astype('datetime64[ns]')
   print("shape:", time.shape)
   print("shape:", lat.shape)
   print("shape:", lon.shape)
   print("shape:", area.shape)
   return time, lat, lon, area, date 

def plot_time_series(var,vdat1,label1,vdat2,label2,year,exp):
   #start to plot data
   import cartopy.crs as ccrs
   import matplotlib.pyplot as plt
   import numpy as np
   import xarray as xr
   time  = vdat1.coords["time"]
   fontz = 16
   plt.rcParams.update({'font.size':fontz})
   f, ax = plt.subplots(1, 1, figsize=(12,6), sharey=True, layout='constrained')
   ax.plot(time, vdat1[var].data, label=label1, color='red',   linewidth=3, linestyle='-')  
   ax.plot(time, vdat2[var].data, label=label2, color='black', linewidth=3, linestyle='--')  
   ax.set_xlabel('Model Time')   # Add an x-label to the Axes.
   ax.set_ylabel('{} ({})'.format(vdat1.long_name,vdat1.units)) # Add a y-label to the Axes.
   ax.set_title("Time seires of {}".format(var))  # Add a title to the Axes.
   ax.legend()  # Add a legend.
   #plt.show()
   #save figure 
   figname = "time_series_{}_{}_{:04d}.png".format(exp,var,year)
   plt.savefig(os.path.join("./diag_figure",figname))
   plt.close()

def plot_2d_hov(var,ds,year,exp): 
   #start to plot data
   import cartopy.crs as ccrs
   import matplotlib.pyplot as plt
   import numpy as np
   import xarray as xr

   pcorr = []
   for j in range(0,len(time)):
      corr = stats.pearsonr(data0[j,:,:],data1[j,:,:],axis=1).statistic.shape
      print(corr.shape)
      exit()
      pcorr.append(corr)
   dam = ds.groupby("time.month").mean("time")
   xvals = np.arange(len(ds[var][0,:]))
   vtim1 = ds.time.values.astype('datetime64[ms]').astype('O')
   vtim2 = dam.month.values.astype('O')
   #start to plot data
   #print(np.min(dam.data),np.max(dam.data))
   fontz = 16
   plt.rcParams.update({'font.size':fontz})
   clevs = np.linspace(np.min(dam[var]),np.max(dam[var]),11)
   fig, axs = plt.subplots(2,1, figsize=(10,12))
   #fig.suptitle('{}'.format(var_out))
   cntr0 = axs[0].contourf(xvals,vtim1, ds[var],   clevs,  cmap=plt.cm.jet)
   axs[0].set_title('3hourly Mean', loc='left', fontsize=fontz)
   axs[0].set_xlabel("Model columns (#)",fontsize=fontz)
   cntr1 = axs[1].contourf(xvals,vtim2, dam[var],  clevs,  cmap=plt.cm.jet)
   axs[1].set_title('3hourly to Monthly Mean', loc='left', fontsize=fontz)
   axs[1].set_xlabel("Model columns (#)", fontsize=fontz)
   axs[1].set_ylabel("Months", fontsize=fontz)
   plt.subplots_adjust(bottom=0.1, right = 0.8, hspace=0.3)  #, right=0.9, top=0.9)
   cbar_ax = fig.add_axes([0.82, 0.10, 0.02, 0.78]) #(left, bottom, width, height)
   cbar = fig.colorbar(cntr1, cax=cbar_ax, shrink=0.6, orientation='vertical',
                       pad=0.1, aspect=10, extendrect=True)
   cbar.set_label('{} ({})'.format(ds.long_name,ds.units),fontsize=fontz)
   figname = "hov_2d_{}_{}_{:04d}.png".format(exp,var,year)
   plt.savefig(os.path.join("./diag_figure",figname))
   plt.close()
   del(vtim1,vtim2,xvals,ds,dam,clevs,fig,axs,cntr0,cntr1,cbar_ax,cbar)
   return

def plot_2d_map_4var(vnames,vunits,dsm,dso,ymdh,ilev,exp,group,dsgrid):
   import time, os
   import xarray as xr
   import numpy as np
   import numpy.ma as ma
   import matplotlib as mpl
   import matplotlib.pyplot as plt
   import matplotlib.tri as tri
   from   matplotlib.collections import PolyCollection
   import cartopy.crs as ccrs
   import cartopy.feature as cfeature
   import matplotlib.ticker as mticker
   from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

   var_min    = [-2.,-2.,-1.,-1.]
   var_max    = [ 2., 2., 1., 1.]
   #select a time step to plot (first step as example below)
   datm = dsm.isel(lev=ilev).mean(dim="time") 
   dato = dso.isel(lev=ilev).mean(dim="time")

   fontz = 16
   plt.rcParams.update({'font.size':fontz})

   #plt.switch_backend('agg')
   t1 = time.time()    #-- retrieve start time
   #only select first time step: 
   lon  = datm.lon.data #    = ((data.lon.data - 180) % 360) - 180
   lat  = datm.lat.data
   #-- get coordinates and (if not in degree) convert radians to degrees
   clon = dsgrid.grid_center_lon.values
   clat = dsgrid.grid_center_lat.values
   clon_vertices = dsgrid.grid_corner_lon.values
   clat_vertices = dsgrid.grid_corner_lat.values
   # sanity check to ensure data grid are consistent with the grid info 
   if np.max(np.abs(clon-lon)) > 1e-5 or np.max(np.abs(clat-lat)) > 1e-5:
     print("data grid is inconsistent with grid info")
     exit("ERROR in plot_2d_map")

   ncells, nv = clon_vertices.shape[0], clon_vertices.shape[1]
   #-- set title string
   title = 'Nudging tendency'
   #-- set projection
   projection = ccrs.PlateCarree()
   for i in range(2):
     for j in range(2):
       k = j*2 + 1 
       varMin, varMax, nlevs = var_min[k],var_max[k],11
       vnam = vnames[k]
       vunt = vunits[k]
       #-- plot the title string
       #plt.title(title)
       #-- define color map
       cmap     = plt.get_cmap('Spectral_r', nlevs)        #-- read the color map
       cmaplist = [i for i in range(cmap.N)]               #-- color bar indices
       ncol     = len(cmaplist)                            #-- number of colors
       colors   = np.ones([ncells,4], np.float32)       #-- assign color array for triangles

       #-- set contour levels, labels
       levels = np.linspace(varMin,varMax,nlevs)
       if np.max(np.abs(levels)) < 1:
         xformat = '%.2f'
         labels = ['{:.2f}'.format(x) for x in levels]
       else:
         xformat = '%0.1f'
         labels = ['{:.1f}'.format(x) for x in levels]

       #-- create figure and axes instances; we need subplots for plot and colorbar
       fig, axs = plt.subplots(1,2,figsize=(20,8), subplot_kw=dict(projection=projection))
       model = ["Truth (Nudging)", "ML prediction (Unet)"]
       gl = []
       for gg in range(2):
         if gg == 0: 
           var = dato[vnam]
         else:
           var = datm[vnam]
         #-- print information to stdout
         print('')
         print('Cells:            %6d ' % clon.size)
         print('Variable min/max: %6.2e ' % np.nanmin(var)+'/'+' %.2e' % np.nanmax(var))
         print('Contour  min/max: %6.2e ' % varMin+'/'+' %.2e' % varMax)
         print('')
         print('levels:      ',levels)
         print('nlevs:       %3d' %nlevs)
         print('ncol:        %3d' %ncol)
         print('')
         #-- set color index of all cells in between levels
         for m in range(0,ncol-1):
           vind = []
           for i in range(0,ncells-2, 1):
             if (var[i] >= levels[m] and var[i] < levels[m+1]):
                colors[i,:] = cmap(cmaplist[m])
                vind.append(i)
           print('set colors: finished level %3d' % m ,
                 ' -- %5d ' % len(vind) ,
                 ' polygons considered')
           del vind

         colors[np.where(var < varMin),:]  = cmap(cmaplist[0])
         colors[np.where(var >= varMax),:] = cmap(cmaplist[ncol-1])
         #-- plot the grid cells 
         clon_vertices = np.where(clon_vertices < -180., clon_vertices + 360., clon_vertices)
         clon_vertices = np.where(clon_vertices >  180., clon_vertices - 360., clon_vertices)
         triangles = np.zeros((ncells, nv, 2), np.float32)
         for i in range(0, ncells, 1):
             triangles[i,:,0] = np.array(clon_vertices[i,:])
             triangles[i,:,1] = np.array(clat_vertices[i,:])
         print('')
         print('--> triangles done')

         axs[gg].tick_params(labeltop=False,labelright=False,direction='out', 
                                  length=6, width=2, colors='r',
                                  grid_color='r', grid_alpha=0.5)

         axs[gg].set_title(model[gg], loc='left', fontsize=fontz)
         axs[gg].set_global()
         #-- plot land areas at last to get rid of the contour lines at land
         axs[gg].coastlines(linewidth=0.5, zorder=2)

         #-- create polygon/triangle collection
         coll = PolyCollection(triangles, array=None,fc=colors, edgecolors='black',
                               linewidth=0.05, transform=ccrs.Geodetic(), zorder=0)
         axs[gg].add_collection(coll)
         print('--> polygon collection done')
  
         gl = axs[gg].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, 
                                     linewidth=0.5,color='dimgray',alpha=0.4,zorder=2)
         gl.xlabels_top = False
         gl.ylabels_right= False

         #-- plot the title string
         #plt.title(title)

         cbar_ax = fig.add_axes([0.13, 0.19, 0.78, 0.02]) #(bottom,right,width,height)
         cb      = plt.cm.ScalarMappable(cmap=cmap,
                       norm=plt.Normalize(vmin=varMin, vmax=varMax))
         cbar    = plt.colorbar(cb,cax=cbar_ax, 
                                orientation='horizontal',
                                ticks=levels,
                                boundaries=levels, 
                                format=xformat,
                                shrink=0.8,
                                pad=0.04,
                                aspect=30,)

         #plt.setp(cbar.ax.get_xticklabels()[::2], visible=False)
         cbar.set_label('{}({})'.format(vnam,vunt))

       #-- maximize and save the PNG file
       figname = os.path.join("./diag_figure",
                              "map_2d_{}_{}_{}.png".format(exp,vnam,ymdh))
       plt.savefig(figname, bbox_inches='tight',dpi=300)
 
   #-- get wallclock time
   t2 = time.time()
   print('Wallclock time:  %0.3f seconds' % (t2-t1))
   print('')
   return 

def triangulate(vertices, x="Longitude", y="Latitude"):
    from scipy.spatial import Delaunay
    """
    Generate a triangular mesh for the given x,y, z vertices, using Delaunay triangulation.
    For large n, typically results in about double the number of triangles as vertices.
    """
    triang = Delaunay(vertices[[x, y]].values)
    print('Given', len(vertices), "vertices, created", len(triang.simplices), 'triangles.')
    return pd.DataFrame(triang.simplices, columns=['v0', 'v1', 'v2'])

def get_var_info (var, grp, datadir):
   import json
   from pprint import pprint
   print(os.path.join(datadir,'variable_dictionary.json'))
   var_dic = json.load(open(os.path.join(wkdir,'variable_dictionary.json'))) 
   print("----------------------")
   print("Inquire variable: ",var)
   pprint(var_dic[grp][var])
   print("----------------------")
   vunit = var_dic[grp][var]["units"]
   vname = var_dic[grp][var]["longname"]
   return vunit,vname 

if __name__ == "__main__":

   datadir    = "/global/cfs/cdirs/e3sm/www/zhan391/SEA_CROGS/New_Training"
   mldir      = "/pscratch/sd/z/zhan391/SEACROGS_project/e3sm_model/machine_learning/offline_results"
   exp        = "Unet_4var"
   time_unit  = "hours since 2009-01-01 00:00:00"

   #output path
   outpath    = "./"

   year = 2009
   main(exp,year,time_unit,datadir,mldir,outpath)

