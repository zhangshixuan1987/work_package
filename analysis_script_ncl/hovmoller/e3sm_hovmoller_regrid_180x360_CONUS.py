'''
Description: Python script using PyNGL Python module

 - Hovmoller diagram for v-wind over Northern hemisphere mide-latitudes
 - This plots illustrate upper-level wave and energy propogation such 
   as the downstream baroclinic development

shixuan.zhang@pnnl.gov
'''
from datetime import datetime
from scipy import stats
from sklearn.metrics import mean_squared_error
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib 
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
import numpy as np
import xarray as xr
import Ngl, Nio, os

#-----------------------------------------------------------------------
#-- Function: main
#-----------------------------------------------------------------------
def main():
    diri    = '/global/cscratch1/sd/zhan391/DARPA_project/post_process/trainning_data/FV_180x360' #-- data directory
    varName = 'V'
    varDesc = 'Meridional wind'
    varUnit = 'm $s^{-1}$'
    expndg  = 'NDGUV' # nudged experiments for plotting 
    expndg  = 'NDGUVT' # nudged experiments for plotting 
    expndg  = 'NDGUVTQ' # nudged experiments for plotting 

    # Create time slice from dates
    season     = 'JJA'
    plev       = 250
    start_time = '2010-06-01'
    end_time   = '2010-08-31'
    #season     = 'DJF'
    #plev       = 250
    #start_time = '2009-12-01'
    #end_time   = '2010-02-28'

    # Create slice variables subset domain
    time_slice = slice(start_time, end_time)
    region     = "CONUS"           
    regstr     = "CONUS (20-55N)"  
    minlat     = 20                
    maxlat     = 55                
    minlon     = 205   
    maxlon     = 315   
    cenlon     = 180
    lat_slice  = slice(minlat, maxlat) #slice(40, 60)
    lon_slice  = slice(minlon, maxlon) #slice(0, 360)

    # minimum and maximum range for plot 
    minval  = -50
    maxval  =  51
    incval  =  5

    # data for the plot 
    expname = []
    filname = []
    expname.append('ERA5')
    filname.append(diri+'/reference_data/ML_REF_ERA5_3hourly_200911302100-201011301800.nc')
    expname.append('CLIM')
    filname.append(diri+'/clim/E3SMv2_CLIM_3hourly_200911302100-201011301800.nc')

    if (expndg == 'NDGUV'): 
        expname.append('NDGUV')
        filname.append(diri+'/after_nudging/ML_IN_E3SMv2_NDGUV_tau6_3hourly_200911302100-201011301800.nc')
    elif (expndg == 'NDGUVT'): 
        expname.append('NDGUVT')
        filname.append(diri+'/after_nudging/ML_IN_E3SMv2_NDGUVT_tau6_3hourly_200911302100-201011301800.nc')
    else: 
        expname.append('NDGUVTQ')
        filname.append(diri+'/after_nudging/ML_IN_E3SMv2_NDGUVTQ_tau6_3hourly_200911302100-201011301800.nc')

    #read the reference data 
    print ('working on ' + expname[0])
    ds  = xr.open_dataset(filname[0])

    # find required vertical levels
    lev = ds.variables["lev"][:]
    dp  = np.abs(lev - plev)
    ik  = np.where(dp == np.min(dp))

    # Get data, selecting time, level, lat/lon slicea
    data = ds[varName].sel(time=time_slice,
              lev=lev[ik], 
              lat=lat_slice,
              lon=lon_slice)

    # Specify longitude values for chosen domain
    lats = data.lat.values
    lons = data.lon.values

    # Get times and make array of datetime objects
    # date_time = datetime.fromtimestamp(data.time.values)
    # vtimes = date_time.strftime("%Y%m%d, %HZ")
    vtimes = data.time.values.astype('datetime64[ms]').astype('O')
    timstr = []
    for k in range(len(vtimes)):
        #{0:%Y-%m-%d %HZ}
        timstr.append('{0:%Y-%m-%d}'.format(vtimes[k]))

    # Read the CLIM and Nudged simulation
    print ('working on ' + expname[1])
    ds1  = xr.open_dataset(filname[1])
    # Get data, selecting time, level, lat/lon slicea
    dat1 = ds1[varName].sel(time=time_slice,
              lev=lev[ik], 
              lat=lat_slice,
              lon=lon_slice)

    print ('working on ' + expname[2])
    ds2  = xr.open_dataset(filname[2])
    # Get data, selecting time, level, lat/lon slicea
    dat2 = ds2[varName].sel(time=time_slice,
              lev=lev[ik], 
              lat=lat_slice,
              lon=lon_slice)

    # Compute weights and take weighted average over latitude dimension
    weights  = np.cos(np.deg2rad(lats))
    ref_var  = (data[:,0,:,:] * weights[None, :, None]).sum(dim='lat') / np.sum(weights)
    clm_var  = (dat1[:,0,:,:] * weights[None, :, None]).sum(dim='lat') / np.sum(weights)
    ndg_var  = (dat2[:,0,:,:] * weights[None, :, None]).sum(dim='lat') / np.sum(weights)

    # Calculate the pattern correlation and RMSD 
    corr = []
    rmsd = []
    for loop in range(2):
        if(loop == 0 ):
            x = np.ravel(clm_var)
            y = np.ravel(ref_var)
        else:
            x = np.ravel(ndg_var)
            y = np.ravel(ref_var)
        rc,_  = stats.pearsonr(x, y)
        mse   = mean_squared_error(x,y)
        re    = np.sqrt(mse)
        corr.append('{0:.2f}'.format(rc))
        rmsd.append('{0:.2f}'.format(re))

    # Start figure
    figname = "fig_hovmoller_"+region+"_"+expndg+"_"+varName+"_"+season+"_"+"{0:.0f}hPa".format(plev)+".png"
    fntsz   = 12
    cmap    = plt.cm.bwr
    fig     = plt.figure(figsize=(10, 12), constrained_layout=True)
    # Use gridspec to help size elements of plot; small top plot and big bottom plot
    gs      = fig.add_gridspec(nrows=4, ncols=3, height_ratios=[1, 4, 1, 1], hspace=0.03)

    # Tick labels
    if( region == "NH_MidLat"): 
        x_tick_labels = [u'0\N{DEGREE SIGN}E', u'90\N{DEGREE SIGN}E',
                u'180\N{DEGREE SIGN}E', u'90\N{DEGREE SIGN}W',
                u'0\N{DEGREE SIGN}E']
        # Top plot for geographic reference (makes small map)
        ax1 = fig.add_subplot(gs[0, :], projection=ccrs.PlateCarree(central_longitude=cenlon))
        ax1.set_extent([0, 357.5, 35, 65], ccrs.PlateCarree(central_longitude=cenlon))
        ax1.set_yticks([40, 60])
        ax1.set_yticklabels([u'40\N{DEGREE SIGN}N', u'60\N{DEGREE SIGN}N'],fontsize=fntsz)
        ax1.set_xticks([-180, -90, 0, 90, 180])
        ax1.set_xticklabels(x_tick_labels,fontsize=fntsz)
        ax1.grid(linestyle='dotted', linewidth=1.0)
    elif ( region == "Tropics"): 
        x_tick_labels = [u'0\N{DEGREE SIGN}E', u'90\N{DEGREE SIGN}E',
                u'180\N{DEGREE SIGN}E', u'90\N{DEGREE SIGN}W',
                u'0\N{DEGREE SIGN}E']
        # Top plot for geographic reference (makes small map)
        ax1 = fig.add_subplot(gs[0, :], projection=ccrs.PlateCarree(central_longitude=cenlon))
        ax1.set_extent([0, 357.5, -20, 20], ccrs.PlateCarree(central_longitude=cenlon))
        ax1.set_yticks([-15, 15])
        ax1.set_yticklabels([u'15\N{DEGREE SIGN}S', u'15\N{DEGREE SIGN}N'],fontsize=fntsz)
        ax1.set_xticks([-180, -90, 0, 90, 180])
        ax1.set_xticklabels(x_tick_labels,fontsize=fntsz)
        ax1.grid(linestyle='dotted', linewidth=1.0)
    else: 
        x_tick_labels = []
        xlon = np.arange(minlon+5,maxlon-5,30)
        xlon = xlon - 360
        for kk in range(len(xlon)):
            x_tick_labels.append(u'{} \N{DEGREE SIGN}W'.format(np.abs(xlon[kk])))
        # Top plot for geographic reference (makes small map)
        ax1 = fig.add_subplot(gs[0, :], projection=ccrs.PlateCarree())
        ax1.set_extent([minlon-360, maxlon-360, minlat-5, maxlat+5], ccrs.PlateCarree())
        ax1.set_yticks([20, 55])
        ax1.set_yticklabels([u'20\N{DEGREE SIGN}N', u'55\N{DEGREE SIGN}N'],fontsize=fntsz)
        ax1.set_xticks(xlon)
        ax1.set_xticklabels(x_tick_labels,fontsize=fntsz)
        ax1.grid(linestyle='dotted', linewidth=1.0)


    # Add geopolitical boundaries for map reference
    # Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
    states_provinces = cfeature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale='50m',
            linewidths=0.1, 
            facecolor='none')
    ax1.add_feature(cfeature.BORDERS.with_scale('50m'),linewidths=0.1)
    ax1.add_feature(cfeature.LAND.with_scale('50m'),linewidths=0.1)
    ax1.add_feature(cfeature.OCEAN.with_scale('50m'),facecolor='lightblue')
    ax1.add_feature(cfeature.COASTLINE.with_scale('50m'),linewidths=0.1)
    ax1.add_feature(cfeature.RIVERS.with_scale('50m'),linewidth=0.1)
    ax1.add_feature(cfeature.LAKES.with_scale('50m'), linewidths=0.1)
    ax1.add_feature(states_provinces, edgecolor='gray')

    # Set some titles
    plt.title(regstr, loc='left')
    #plt.title(expname[i], loc='right')
    # Set some titles
    #plt.title('{0:.1f}-hPa V-wind'.format(float(lev[ik].values)), loc='left', fontsize=fntsz)
    #plt.title('Time Range: {0:%Y%m%d %HZ} - {1:%Y%m%d %HZ}'.format(vtimes[0], vtimes[-1]),
    #        loc='right', fontsize=fntsz)

    # Bottom plot for Hovmoller diagram
    for ip in range(3):
        if(ip == 0):
            avg_data = ref_var
        elif(ip == 1):
            avg_data = clm_var
        else:
            avg_data = ndg_var

        ax2 = fig.add_subplot(gs[1, ip])
        ax2.invert_yaxis()  # Reverse the time order to do oldest first
        
        # Plot of chosen variable averaged over latitude and slightly smoothed
        clevs = np.arange(minval, maxval, incval)
        cf = ax2.contourf(lons, vtimes, mpcalc.smooth_n_point(
            avg_data, 9, 2), clevs, cmap=cmap, extend='both')
        cs = ax2.contour(lons, vtimes, mpcalc.smooth_n_point(
            avg_data, 9, 2), clevs, colors='k', linewidths=0.2)
        
        # Make some ticks and tick labels
        if( region == "NH_MidLat"):
            ax2.set_xticks([0, 90, 180, 270, 357.5])
            ax2.set_xticklabels(x_tick_labels,fontsize = fntsz)
        elif( region == "Tropics"):
            ax2.set_xticks([0, 90, 180, 270, 357.5])
            ax2.set_xticklabels(x_tick_labels,fontsize = fntsz)
        else:
            x_tick_labels = []
            xlon = np.arange(minlon+5,maxlon-5,30)
            for kk in range(len(xlon)):
                x_tick_labels.append(u'{} \N{DEGREE SIGN}W'.format(np.abs(xlon[kk]-360)))
            ax2.set_xticks(xlon)
            ax2.set_xticklabels(x_tick_labels,fontsize = fntsz)

        ax2.set_yticks(vtimes[4::32])
        if (ip == 0): 
            ax2.set_yticklabels(timstr[4::32],fontsize = fntsz)
        else:
            ax2.set_yticklabels("")
        
        # Set some titles
        plt.title(expname[ip], loc='left')
        plt.title('{0:.1f}-hPa V-wind'.format(float(lev[ik].values)), loc='right', fontsize=fntsz)

    # attach a color bar
    cbar_ax2 = fig.add_axes([0.12, 0.24, 0.85, 0.015])
    cb2 = matplotlib.colorbar.ColorbarBase(cbar_ax2, cmap=cmap, ticks=np.arange(0, 1.001, 0.05), orientation='horizontal')
    cb2.set_ticklabels(np.arange(minval, maxval, incval))
    cbar_ax2.tick_params(labelsize=12)
    cbar_ax2.text(0.5, -2.0, varUnit, ha='center', va='top', fontsize=fntsz)

    # attach a table for mean statistics
    ax4 = fig.add_axes([0.22, 0.13, 0.75, 0.015])
    fig.patch.set_visible(False)
    ax4.axis('off')
    ax4.set_title('Hovmoller diagram during {0:%Y%m%d %HZ} - {1:%Y%m%d %HZ}'.format(vtimes[0], vtimes[-1]), pad = 24.0 )
    
    # Set row headers (first column)
    row_text  = [ "Mean metrics", "Corr. (unitless)","RMSE ("+varUnit+")"]
    rowcolors = ['lightgray', 'lightgray', 'lightgray']
    # Set cell text 
    cell_text = [[expname[1], expname[2]], 
                 [corr[0],corr[1]],
                 [rmsd[0],rmsd[1]]]
    colors    = [['lightgray', 'lightgray'],
                 ["White", "White"],
                 ["White", "White"]]
    #draw tables 
    table = ax4.table(rowLabels=row_text,
                 rowColours=rowcolors,
                 rowLoc='center',
                 cellText=cell_text,
                 cellColours=colors,
                 cellLoc='center',
                 loc='center', 
                 fontsize=fntsz*1.2)

    plt.savefig(figname)


if __name__ == '__main__':
    main()


