'''
Description: Python script using PyNGL Python module

 - contour plot on map (rectilinear data)

shixuan.zhang@pnnl.gov
'''
import argparse
import os
from pathlib import Path
import numpy as np
import xarray as xr
import pickle
import pandas as pd
from numpy import unravel_index

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import cartopy.crs as ccrs                 # For plotting maps
import cartopy.feature as cfeature         # For plotting maps
from cartopy.util import add_cyclic_point  # For plotting maps
import plotly.express as px

#Green Function Experiments 
target = "SOM_Patches"
npatch = 108 # total patch number 
nyr    = 150 # note: total number of years
nmo    = 12  # note: total number of month
DEFAULT_OUTPUT_ROOT = Path(
    "/global/cfs/cdirs/m1867/zhan391/paper_work/green_function"
)
diri = None
outdr = None
figures_dir = None
#There are 18 variables in total, the first three dimension in each data file, 
#i.e (q_sign, q_lon, q_lat) provides the info for warm/ cold patches. 
#the coordinates “q_lon” and “q_lat” indicate the location of the logitude and 
#latitude center of the q-flux forcing, while the coordinate “q_sign” indicates 
#the responses are forced with either the warm or the cold patch.
#In addition to (q_sign, q_lon, q_lat), we can group the quantities into three types: 
#The 3-D(time,lat,lon) variables (surface level), time = nyr x nmo months 
varList1 = ['TS','PRECT','FSNT','FLNT','FSNS','FLNS','LHFLX', 'SHFLX']
#The 4-D(time,lev,lat,lon) variables (model level), time = nyr x nmo months
varList2 = ['U','V','T','Q']
#The 2-D(lat,lon) variabels, climatological annual mean (average over all nyr x nmo)
varList2 = ['R_plk','R_alb','R_lr','R_q','R_SWcld','R_LWcld']

def map_plot(lat,lon,var,projection,title,figname):

    var_cyc, lon_cyc = add_cyclic_point(var, coord=lon)
    var_cyc = np.where(abs(var_cyc) < 1.0,np.nan,var_cyc)

    #define contour levels
    clev = np.arange(-12,1, 1)

    plt.figure(figsize=(20,8))

    #Define projection
    if projection == "robinson":
        ax = plt.axes(projection=ccrs.Robinson())
        # make the map global rather than have it zoom in to
        # the extents of any plotted data
        ax.set_global()
        #ax.stock_img()
        #ax.coastlines()
        #plt.plot(-0.08, 51.53, 'o', transform=ccrs.PlateCarree())
        #plt.plot([-0.08, 132], [51.53, 43.17], transform=ccrs.PlateCarree())
        #plt.plot([-0.08, 132], [51.53, 43.17], transform=ccrs.Geodetic())
        #plot the data
        if(len(var_cyc.shape) >2):
            for i in np.arange(len(var_cyc[:,0,0])):
                plt.contourf(lon_cyc,lat,var_cyc[i,:,:],clev,cmap='Spectral_r',transform=ccrs.PlateCarree())
        else:
            plt.contourf(lon_cyc,lat,var_cyc,clev,cmap='Spectral_r',transform=ccrs.PlateCarree())
    else:
        ax = plt.axes(projection=ccrs.PlateCarree())
        #ax.set_extent([0, 360, -90, 90])
        #plot the data
        plt.contourf(lon_cyc,lat,var_cyc,clev,cmap='Spectral_r',extend='both')


    # Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')

    # add coastlines
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(states_provinces, edgecolor='gray')

    #add lat lon grids
    gl = ax.gridlines(draw_labels=True,color='grey', alpha=0.8, linestyle='--')
    gl.top_labels   = False
    gl.right_labels = False

    # Titles
    # Main
    plt.rcParams['font.size'] = 18
    # Set the default text font size
    plt.rc('font', size=16)
    # Set the axes title font size
    plt.rc('axes', titlesize=16)
    # Set the axes labels font size
    plt.rc('axes', labelsize=16)
    # Set the font size for x tick labels
    plt.rc('xtick', labelsize=16)
    # Set the font size for y tick labels
    plt.rc('ytick', labelsize=16)
    # Set the legend font size
    plt.rc('legend', fontsize=18)
    # Set the font size of the figure title
    plt.rc('figure', titlesize=20)

    plt.title(title,fontsize=18)

    # y-axis
    ax.text(-0.08, 0.5, 'Latitude', va='bottom', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax.transAxes)

    # x-axis
    ax.text(0.5, -0.12, 'Longitude', va='bottom', ha='center',
        rotation='horizontal', rotation_mode='anchor',
        transform=ax.transAxes)

    # legend
    ax.text(1.15, 0.5, 'Q-flux', va='bottom', ha='center',
        rotation='vertical', rotation_mode='anchor',
        transform=ax.transAxes)

    plt.colorbar()

    #plt.show()
    plt.savefig(figures_dir / Path(figname).name)

def patch_plot(lats,lons,var):
    # draw map with markers for float locations
    m = Basemap(projection='hammer',lon_0=180)
    x, y = m(lons,lats)
    m.drawmapboundary(fill_color='#99ffff')
    m.fillcontinents(color='#cc9966',lake_color='#99ffff')
    m.scatter(x,y,3,marker='o',color='k')
    plt.title('Locations of %s patches for Green Function Experiment %s and %s' %\
              (len(lats),date1,date2),fontsize=12)
    plt.show()

def xy_plot(x,y,xmin,xmax,xint,ymin,ymax,yint,title,legend,xlabel,ylabel,figname):

    # plot
    plt.style.use('_mpl-gallery')
    plt.rcParams["figure.figsize"] = [7.50, 3.50]
    plt.rcParams["figure.autolayout"] = True
    # Main
    plt.rcParams['font.size'] = 18
    # Set the default text font size
    plt.rc('font', size=16)
    # Set the axes title font size
    plt.rc('axes', titlesize=16)
    # Set the axes labels font size
    plt.rc('axes', labelsize=16)
    # Set the font size for x tick labels
    plt.rc('xtick', labelsize=16)
    # Set the font size for y tick labels
    plt.rc('ytick', labelsize=16)
    # Set the legend font size
    plt.rc('legend', fontsize=18)
    # Set the font size of the figure title
    plt.rc('figure', titlesize=20)

    fig, ax = plt.subplots()
    plt.title(title,fontsize=18)

    ax.plot(x, y[0][:], linewidth=2.0,
            label = legend[0], linestyle="-", color = 'red')
    ax.plot(x, y[1][:], linewidth=2.0,
            label = legend[1], linestyle="-", color = 'blue')

    if abs(ymax) > 10: 
        yubd = ymax + 1
    else:
        yubd = ymax + 0.1

    if abs(xmax) > 10:
        xubd = xmax + 1
    else:
        xubd = xmax + 0.1

    ax.set(xlim=(xmin, xmax), xticks=np.arange(xmin, xubd,xint),
           ylim=(ymin, ymax), yticks=np.arange(ymin, yubd))

    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(start, end, xint))
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))

    ax.tick_params(axis ='both', which ='both', length = 0)
    plt.xticks(np.arange(xmin, xubd, xint))
    plt.yticks(np.arange(ymin, yubd, yint))

    #plt.grid(False)
    plt.axis('on')
    plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    #plt.show()
    plt.savefig(figures_dir / Path(figname).name)


def scatter_plot(x,y,xmin,xmax,xint,ymin,ymax,yint,title,legend,xlabel,ylabel,figname):
    # plot
    plt.style.use('_mpl-gallery')
    plt.rcParams["figure.figsize"] = [7.50, 3.50]
    plt.rcParams["figure.autolayout"] = True
    # Main
    plt.rcParams['font.size'] = 18
    # Set the default text font size
    plt.rc('font', size=16)
    # Set the axes title font size
    plt.rc('axes', titlesize=16)
    # Set the axes labels font size
    plt.rc('axes', labelsize=16)
    # Set the font size for x tick labels
    plt.rc('xtick', labelsize=16)
    # Set the font size for y tick labels
    plt.rc('ytick', labelsize=16)
    # Set the legend font size
    plt.rc('legend', fontsize=18)
    # Set the font size of the figure title
    plt.rc('figure', titlesize=20)

    fig, ax = plt.subplots()
    plt.title(title,fontsize=18)

    ax.scatter(x[0][:], y[0][:], linewidth=2.0,
            label = legend[0], linestyle="-", color = 'red')
    ax.scatter(x[1][:], y[1][:], linewidth=2.0,
            label = legend[1], linestyle="-", color = 'blue')

    if abs(ymax) > 10:
        yubd = ymax + 1
    else:
        yubd = ymax + 0.1

    if abs(xmax) > 10:
        xubd = xmax + 1
    else:
        xubd = xmax + 0.1

    ax.set(xlim=(xmin, xmax), xticks=np.arange(xmin, xubd,xint),
           ylim=(ymin, ymax), yticks=np.arange(ymin, yubd))

    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(start, end, xint))
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))

    ax.tick_params(axis ='both', which ='both', length = 0)
    plt.xticks(np.arange(xmin, xubd, xint))
    plt.yticks(np.arange(ymin, yubd, yint))

    #plt.grid(False)
    plt.axis('on')
    plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    #plt.show()
    plt.savefig(figures_dir / Path(figname).name)

def qflx_patch_mean(lat,lon,qdp,outfile):

    #calculate the weighted sum of each patch
    qdp     = qdp.where(abs(qdp)>0)
    weights = np.cos(np.deg2rad(qdp.lat))
    weights.name = "weights"
    qdp_wgt = qdp.weighted(weights)
    qdp_ave = qdp_wgt.sum(("lat","lon")) #/sum(weights)

    latp = [] 
    lonp = []
    for i in np.arange(qdp.shape[0]):
        #indx = unravel_index(abs(qdp[i,:,:]).argmin(), qdp[i,:,:].shape)
        ilat,ilon = unravel_index(qdp[i,:,:].argmin(), qdp[i,:,:].shape)
        latp.append(lat[ilat].data)
        lonp.append(lon[ilon].data)
        del ilat,ilon
    
    #output data into nc file 
    qdp_pos = qdp_ave.copy()
    qdp_neg = qdp_ave.copy()
    qdp_pos = qdp_pos*-1.0
    qdp_pos = np.nan_to_num(qdp_pos,nan    = -9999.0)
    qdp_neg = np.nan_to_num(qdp_neg,nan    = -9999.0)
    latp    = np.nan_to_num(latp,nan       = -9999.0)
    lonp    = np.nan_to_num(lonp,nan       = -9999.0)
    patch   = np.arange(1,qdp.shape[0]+1,1)
    ds = xr.Dataset(
            data_vars=dict(
                qdp_patch_neg = (["patch"], qdp_pos.data),
                qdp_patch_pos = (["patch"], qdp_neg.data),
                lat_patch     = (["patch"], latp),
                lon_patch     = (["patch"], lonp),
                ),
            coords=dict(
                patch=(["patch"], patch),
                ),
            attrs=dict(description="average qfluxes from Green Function experiments"),)
    #ds.fillna(-9999.0)
    ds.qdp_patch_pos.name                = 'Qflux (+qdp)'
    ds.qdp_patch_pos.attrs['long_name']  = 'Patch mean ocean heat flux convergence (positive patch)'
    ds.qdp_patch_pos.attrs['units']      = 'W/m^2'
    ds.qdp_patch_pos.attrs['_FillValue'] = -9999.0
    ds.qdp_patch_neg.name                = 'Qflux (-qdp)'
    ds.qdp_patch_neg.attrs['long_name']  = 'Patch mean ocean heat flux convergence (negative patch)'
    ds.qdp_patch_neg.attrs['units']      = 'W/m^2'
    ds.qdp_patch_neg.attrs['_FillValue'] = -9999.0
    ds.lat_patch.name                    = 'latitude'
    ds.lat_patch.attrs['long_name']      = 'Latitude at patch center (maximum |qdp|)'
    ds.lat_patch.attrs['units']          = 'degree_north'
    ds.lat_patch.attrs['_FillValue']     = -9999.0
    ds.lon_patch.name                    = 'logitude'
    ds.lon_patch.attrs['long_name']      = 'Longitude at patch center (maximum |qdp|)'
    ds.lon_patch.attrs['units']          = 'degree_east'
    ds.lon_patch.attrs['_FillValue']     = -9999.0
    ds.patch.name                        = 'patch'
    ds.patch.attrs['long_name']          = 'patch number'
    ds.patch.attrs['units']              = '1'

    #save data to ncfile
    ds.to_netcdf(outfile,encoding={'patch': {'dtype': 'i4'}},mode='w')
    del ds

    #plot data
    x = np.arange(0, len(qdp_ave), 1) + 1
    y = []
    y.append(qdp_pos)
    y.append(qdp_neg)
    xmin    = 0
    xmax    = 109
    xint    = 10
    ymin    = -4
    ymax    =  4
    yint    =  1
    title   = 'Q-flux averaged over each patch region'
    figname = 'Qflux_app_patch_time_series.png'
    legend  = ['positive qflux','negative qflux']
    xlabel  = 'patches'
    ylabel  = 'qdp (W/m^2)'
    xy_plot(x,y,xmin,xmax,xint,ymin,ymax,yint,title,legend,xlabel,ylabel,figname)

    title   = 'Central latitude for each patch region (max|qdp|)'
    figname = 'Latitude_app_patch_time_series.png'
    legend  = ['positive qflux','negative qflux']
    xlabel  = 'patches'
    ylabel  = 'latitude (degree)'
    y       = []
    y.append(latp)
    y.append(latp)
    xmin    = 0
    xmax    = 109
    xint    = 10
    ymin    = -90
    ymax    =  90
    yint    =  20 
    xy_plot(x,y,xmin,xmax,xint,ymin,ymax,yint,title,legend,xlabel,ylabel,figname)

    title   = 'Central logitude for each patch region (max|qdp|)'
    figname = 'Longitude_app_patch_time_series.png'
    legend  = ['positive qflux','negative qflux']
    xlabel  = 'patches'
    ylabel  = 'longitude (degree)'
    y       = []
    y.append(lonp)
    y.append(lonp)
    xmin    = 0
    xmax    = 109
    xint    = 10 
    ymin    = 0
    ymax    = 361
    yint    = 60
    xy_plot(x,y,xmin,xmax,xint,ymin,ymax,yint,title,legend,xlabel,ylabel,figname)

    del qdp_pos,qdp_neg,latp,lonp,qdp_ave,patch


def global_mean_2d(qlat,qlon,qsign,lat,lon,TS,infile,outfile):
    df0          = xr.open_dataset(infile,decode_times=False)
    weights      = np.cos(np.deg2rad(TS.lat))
    weights.name = "weights"
    TS_qdp_pos   = []
    TS_qdp_neg   = []
    delta_TS1    = []
    delta_TS2    = []
    fqlat        = []
    fqlon        = []

    for i in np.arange(len(df0.lat_patch)):
        latr  = df0.lat_patch[i].data
        lonr  = df0.lon_patch[i].data
        #print(latr,lonr)
        dlat  = abs(qlat.data - latr) 
        dlon  = abs(qlon.data - lonr) 
        #print(dlat,dlon)
        ilat  = unravel_index(dlat.argmin(), dlat.shape)
        ilon  = unravel_index(dlon.argmin(), dlon.shape)
        #print(ilat,ilon)
        #print(dlat[ilat],dlon[ilon])
        #print(latr,"", lonr," ",qlat[ilat].data," ",qlon[ilon].data)
        #print(lat1[ilat][ilon].data,lon1[ilat][ilon].data)
        x_pos     = TS[1][ilon][ilat][:][:][:].mean(("time"))
        x_neg     = TS[0][ilon][ilat][:][:][:].mean(("time"))
        delx1     = 0.5*(x_pos - x_neg)
        delx2     = 0.5*(x_pos + x_neg)
        x_pos_wgt = x_pos.weighted(weights)
        x_neg_wgt = x_neg.weighted(weights)
        delx1_wgt = delx1.weighted(weights)
        delx2_wgt = delx2.weighted(weights)
        x_pos_ave = x_pos_wgt.mean(("lat","lon")) #/sum(weights)
        x_neg_ave = x_neg_wgt.mean(("lat","lon")) #/sum(weights)
        delx1_ave = delx1_wgt.mean(("lat","lon")) #/sum(weights)
        delx2_ave = delx2_wgt.mean(("lat","lon")) #/sum(weights)
        TS_qdp_pos.append(x_pos_ave.data)
        TS_qdp_neg.append(x_neg_ave.data)
        delta_TS1.append(delx1_ave.data)
        delta_TS2.append(delx2_ave.data)
        fqlat.append(qlat[ilat].data)
        fqlon.append(qlon[ilon].data)
        del latr,lonr,dlat,dlon,ilat,ilon
        del x_pos,x_neg,delx1,delx2,x_pos_wgt,x_neg_wgt,x_pos_ave,x_neg_ave
        del delx1_wgt,delx2_wgt,delx1_ave,delx2_ave

    #scatter plot for patch lat-lon 
    x = [] 
    y = []
    x.append(df0.lon_patch.data)
    x.append(fqlon)
    y.append(df0.lat_patch.data)
    y.append(fqlat)
    xmin    =  0
    xmax    =  360
    xint    =  60
    ymin    = -90.
    ymax    =  90
    yint    =  30
    title   = 'lat-lon location of each patch' 
    figname = 'patch_lat-lon_scatter.png'
    legend  = ['qdp_data','response_data']
    xlabel  = 'Longitude (degree)'
    ylabel  = 'Latitude (degree)'
    scatter_plot(x,y,xmin,xmax,xint,ymin,ymax,yint,title,legend,xlabel,ylabel,figname)

    # remove the point that does not match between qdp and climate response simulation 
    dxdy         = abs(fqlat - df0.lat_patch.data) + abs(fqlon - df0.lon_patch.data)
    delta_TS_pos = np.ma.masked_where(dxdy>=2.0, TS_qdp_pos)
    delta_TS_neg = np.ma.masked_where(dxdy>=2.0, TS_qdp_neg)
    delta_TS_ln  = np.ma.masked_where(dxdy>=2.0, delta_TS1)
    delta_TS_nl  = np.ma.masked_where(dxdy>=2.0, delta_TS2)
    delta_qlat   = np.ma.masked_where(dxdy>=2.0,fqlat)
    delta_qlon   = np.ma.masked_where(dxdy>=2.0,fqlon)
    delta_TS_pos = np.ma.filled(delta_TS_pos,fill_value=np.nan)
    delta_TS_neg = np.ma.filled(delta_TS_neg,fill_value=np.nan)
    delta_TS_ln  = np.ma.filled(delta_TS_ln,fill_value=np.nan)
    delta_TS_nl  = np.ma.filled(delta_TS_nl,fill_value=np.nan)
    delta_qlat   = np.ma.filled(delta_qlat,fill_value=np.nan)
    delta_qlon   = np.ma.filled(delta_qlon,fill_value=np.nan)
    #print(delta_qlat)
    #print(np.count_nonzero(np.isnan(delta_qlat)))
    #print(np.count_nonzero(np.isnan(delta_qlon)))
    #exit()

    x = []
    y = []
    x1 = df0.lon_patch.data.copy() 
    x1 = np.ma.masked_where(dxdy>=2.0, x1)
    x1 = np.ma.filled(x1,fill_value=np.nan) 
    y1 = df0.lat_patch.data.copy()
    y1 = np.ma.masked_where(dxdy>=2.0, y1)
    y1 = np.ma.filled(y1,fill_value=np.nan)
    x.append(x1)
    x.append(delta_qlon)
    y.append(y1)
    y.append(delta_qlat)
    xmin    =  0
    xmax    =  360
    xint    =  60
    ymin    = -90.
    ymax    =  90
    yint    =  30
    title   = 'lat-lon location of each patch'
    figname = 'patch_lat-lon_scatter_miss.png'
    legend  = ['qdp','response']
    xlabel  = 'Longitude (degree)'
    ylabel  = 'Latitude (degree)'
    scatter_plot(x,y,xmin,xmax,xint,ymin,ymax,yint,title,legend,xlabel,ylabel,figname)

    #plot data
    x = np.arange(0, len(df0.lat_patch), 1) + 1
    y = []
    y.append(delta_TS_pos)
    y.append(delta_TS_neg)
    xmin    = 0
    xmax    = 109
    xint    = 10
    ymin    = -0.4
    ymax    =  0.3
    yint    =  0.1
    title   = 'Global mean TS response to qflux'
    figname = 'TS1_response_108_patch_time_series.png'
    legend  = ['TS (positive)','TS (negative)']
    xlabel  = 'patches'
    ylabel  = 'delta TS (K)'
    xy_plot(x,y,xmin,xmax,xint,ymin,ymax,yint,title,legend,xlabel,ylabel,figname)

    x = np.arange(0, len(df0.lat_patch), 1) + 1
    y = []
    y.append(delta_TS_ln)
    y.append(delta_TS_nl)
    xmin    = 0
    xmax    = 109
    xint    = 10
    ymin    = -0.3
    ymax    =  0.1
    yint    =  0.1
    title   = 'Global mean TS response (linear vs nonlinear)'
    figname = 'TS2_response_108_patch_time_series.png'
    legend  = ['TS (linear)','TS (nonlinear)']
    xlabel  = 'patches'
    ylabel  = 'delta TS (K)'
    xy_plot(x,y,xmin,xmax,xint,ymin,ymax,yint,title,legend,xlabel,ylabel,figname)

    #output data into nc file 
    patch = np.arange(1,len(df0.lat_patch)+1,1)
    delta_TS_pos = np.nan_to_num(delta_TS_pos,nan  = -9999.0)
    delta_TS_neg = np.nan_to_num(delta_TS_neg,nan  = -9999.0)
    delta_TS_ln  = np.nan_to_num(delta_TS_ln, nan  = -9999.0)
    delta_TS_nl  = np.nan_to_num(delta_TS_nl, nan  = -9999.0)
    delta_qlat   = np.nan_to_num(delta_qlat,  nan  = -9999.0)
    delta_qlon   = np.nan_to_num(delta_qlon,  nan  = -9999.0)
    patch        = np.arange(1,len(delta_TS1)+1,1)
    ds = xr.Dataset(
            data_vars=dict(
                delta_TS_pos = (["patch"], delta_TS_pos),
                delta_TS_neg = (["patch"], delta_TS_neg),
                delta_TS_ln  = (["patch"], delta_TS_ln),
                delta_TS_nl  = (["patch"], delta_TS_nl),
                qlat         = (["patch"], delta_qlat),
                qlon         = (["patch"], delta_qlon),
                lat_patch    = (["patch"], df0.lat_patch.data),
                lon_patch    = (["patch"], df0.lon_patch.data),
                ),
            coords=dict(
                patch=(["patch"], patch),
                ),
            attrs=dict(description="Climatological mean forced responses from Green Function experiments"),)

    ds.delta_TS_pos.name               = 'TS (+qdp)'
    ds.delta_TS_pos.attrs['long_name'] = 'Forced response of Surface temperature (positive qflux)'
    ds.delta_TS_pos.attrs['units']     = 'K'
    ds.delta_TS_pos.attrs['_FillValue']= -9999.0
    ds.delta_TS_neg.name               = 'TS (-qdp)'
    ds.delta_TS_neg.attrs['long_name'] = 'Forced response of Surface temperature (negative qflux)'
    ds.delta_TS_neg.attrs['units']     = 'K'
    ds.delta_TS_neg.attrs['_FillValue']= -9999.0
    ds.delta_TS_ln.name                = 'TS (linear)'
    ds.delta_TS_ln.attrs['long_name']  = 'Forced response of Surface temperature (linear) 0.5*[TS (+qdp) - TS (-qdp)]'
    ds.delta_TS_ln.attrs['units']      = 'K'
    ds.delta_TS_ln.attrs['_FillValue'] = -9999.0
    ds.delta_TS_nl.name                = 'TS (nonlinear)'
    ds.delta_TS_nl.attrs['long_name']  = 'Forced response of Surface temperature (nonlinear) 0.5*[TS (+qdp) + TS (-qdp)]'
    ds.delta_TS_nl.attrs['units']      = 'K'
    ds.delta_TS_nl.attrs['_FillValue'] = -9999.0
    ds.lat_patch.name                  = 'latitude'
    ds.lat_patch.attrs['long_name']    = 'Latitude at patch center (maximum |qdp|)'
    ds.lat_patch.attrs['units']        = 'degree_north'
    ds.lat_patch.attrs['_FillValue']   = -9999.0
    ds.lon_patch.name                  = 'logitude'
    ds.lon_patch.attrs['long_name']    = 'Longitude at patch center (maximum |qdp|)'
    ds.lon_patch.attrs['units']        = 'degree_east'
    ds.lon_patch.attrs['_FillValue']   = -9999.0

    ds.qlat.name                       = 'latitude'
    ds.qlat.attrs['long_name']         = 'Latitude of q-flux center (in SOM_Patches*.nc)'
    ds.qlat.attrs['units']             = 'degree_north'
    ds.qlat.attrs['_FillValue']        = -9999.0
    ds.qlon.name                       = 'logitude'
    ds.qlon.attrs['long_name']         = 'Longitude of q-flux center (in SOM_Patches*.nc)'
    ds.qlon.attrs['units']             = 'degree_east'
    ds.qlon.attrs['_FillValue']        = -9999.0

    #save data to ncfile
    ds.to_netcdf(outfile,encoding={'patch': {'dtype': 'i4'}},mode='w')
    del ds, TS_qdp_pos, TS_qdp_neg, delta_TS1, delta_TS2

#-----------------------------------------------------------------------
#-- Function: main
#-----------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        required=True,
        help="External directory containing qdp.PatchEXP.108.nc and SOM_Patches_TS.nc.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="External root for generated outputs, figures, and logs.",
    )
    return parser.parse_args()


def main():
    global diri, outdr, figures_dir
    args = parse_args()
    diri = args.data_dir.expanduser().resolve()
    run_root = args.output_root.expanduser() / "GFexp2_patch_qflux"
    outdr = run_root / "outputs"
    figures_dir = run_root / "figures"
    logs_dir = run_root / "logs"
    for directory in (outdr, figures_dir, logs_dir):
        directory.mkdir(parents=True, exist_ok=True)

    fqflx = diri / "qdp.PatchEXP.108.nc"
    fts   = diri / "SOM_Patches_TS.nc"

    # read the qdp and calculate patch weighted sum
    df0   = xr.open_dataset(fqflx,decode_times=False)
    qdp_patch = df0.qdp
    lat_patch = df0.lat
    lon_patch = df0.lon
    wgt       = np.cos(np.deg2rad(lat_patch))
    wgt.name  = "wgt"
    qdp_wgt   = qdp_patch.weighted(wgt)
    qdp_sum   = qdp_wgt.sum(("lat","lon")) / wgt.sum()

    #read the variable TS(q_sign,q_lon,q_lat,time,lat,lon)
    df1       = xr.open_dataset(fts,decode_times=False)
    TS        = df1.TS[:,:,:,0,0,0]
    q_sign    = df1.q_sign
    q_lat     = df1.q_lat
    q_lon     = df1.q_lon

    qdp3d     = df1.TS[:,:,:,0,0,0].copy()
    qdp3d.attrs['units'] = qdp_patch.attrs['units']
    qdp3d.attrs['long_name'] = qdp_patch.attrs['long_name']

    #find corresponding patch qfluxes 
    for i in np.arange(len(q_lat)):
        for j in np.arange(len(q_lon)):
            dxy  = abs(lat_patch - q_lat[i]) + abs(lon_patch - q_lon[j])
            ixy  = unravel_index(dxy.argmin(), dxy.shape)
            print(ixy)
            print(q_lat[i], q_lon[j]) 
            print(lat_patch[ixy], lon_patch[ixy])
            exit()
            qdp3d[0,i,j] = qdp_sum[ixy] 
            print(lat_patch)
            print(qdp_patch)
            exit()
    print(qdp3d)
    exit()
    #define the structure to the output patch integral qdp
    ds = xr.Dataset({"qdp": (("q_sign","q_lon", "q_lat"), qdp3d.data)},
                    coords={"q_sign": q_sign.data,
                            "q_lat": q_lon.data,
                            "q_lon": q_lat.data,
                            },)
    ds.q_lat.attrs['long_name']   = "sign of fluxes (negetive or positive)"
    ds.q_lat.attrs['units']       = "degree_north"
    ds.q_lat.attrs['long_name']   = "latitude"
    ds.q_lon.attrs['units']       = "degree_east"
    ds.q_lon.attrs['long_name']   = "longitude"
    ds.qdp.attrs['units']         = qdp_patch.attrs['units']
    ds.qdp.attrs['long_name']     = qdp_patch.attrs['long_name']
    
    for i in np.arange(len(qdp_sum)):
        x       = qdp_patch[i,:,:]
        x       = np.nan_to_num(x,nan=0)
        ix,iy   = np.where(x!=0)
        ix1,iy1 = np.where(x == np.amin(x))
        print(lat_patch[ix1].values,lon_patch[iy1].values)
        exit()
        minlat  = lat_patch[ix].min().values
        maxlat  = lat_patch[ix].max().values
        minlon  = lon_patch[iy].min().values
        maxlon  = lon_patch[iy].max().values
        print(minlat,maxlat,minlon,maxlon)
        exit()
        wgt     = np.cos(np.deg2rad(lat_patch[ix]))
        x_wgt   = qdp.weighted(weights)
        x_sum   = qdp_wgt.sum(("lat","lon")) #/sum(weights)

        print(lats)
        print(lons)
        exit()
        print(y.lat.data)
        exit()
        print(x)
        exit()

    #copy the structure to the output patch mean qdp
    qdp           = TS.copy()
    qdp.units     = "W/m2"
    qdp.long_name = ""
    print(qdp)
    exit()
    var_cyc, lon_cyc = add_cyclic_point(var, coord=lon)
    var_cyc = np.where(abs(var_cyc) < 1.0,np.nan,var_cyc)

    #define contour levels
    clev = np.arange(-12,1, 1)

    plt.figure(figsize=(20,8))

    ax = plt.axes(projection=ccrs.Robinson())
    ax.set_global()
    #ax.stock_img()
    #ax.coastlines()
    #plt.plot(-0.08, 51.53, 'o', transform=ccrs.PlateCarree())
    #plt.plot([-0.08, 132], [51.53, 43.17], transform=ccrs.PlateCarree())
    #plt.plot([-0.08, 132], [51.53, 43.17], transform=ccrs.Geodetic())
    #plot the data
    plt.contourf(lon_cyc,lat,var_cyc,clev,cmap='Spectral_r',transform=ccrs.PlateCarree())

    plt.show()

    exit()

    df0   = xr.open_dataset(fqflx,decode_times=False)
    df1   = xr.open_dataset(fts,decode_times=False)

   
    #check the map distribution for all patches 
    title   = 'Q-flux for Green Function Experiment (total 108 patch)'
    figname = 'Qflux_app_patch_108.png'
    map_plot(df0.lat,df0.lon,df0.qdp,projection="robinson",title=title,figname=figname)

    #calculate weighted average qflux at each patch 
    outfl_qfl = outdr / "qflux_forcing_108_patches.nc"
    if os.path.isfile(outfl_qfl):
        os.remove(outfl_qfl)
    qflx_patch_mean(df0.lat,df0.lon,df0.qdp,outfile=outfl_qfl)

    #calculate the weighted average TS over globe
    outfl_ts = outdr / "ts_clim_response_to_qflx_108_patches.nc"
    if os.path.isfile(outfl_ts):
        os.remove(outfl_ts)
    global_mean_2d(df1.q_lat,df1.q_lon,df1.q_sign,df1.lat,df1.lon,df1.TS,infile=outfl_qfl,outfile=outfl_ts)
 

    #check the map distribution for non-missing patches
    vqdp    = df0.qdp
    df2     = xr.open_dataset(outfl_ts,decode_times=False)
    delts   = df2.delta_TS_neg
    #print(np.count_nonzero(np.isnan(vqdp.data)))
    for i in np.arange(len(delts)):
        if delts[i].isnull():
            vqdp[i,:,:] = np.nan
    num_nnan = len(delts) - np.count_nonzero(np.isnan(delts.data))
    print(np.count_nonzero(np.isnan(vqdp.data)))
    title   = 'Q-flux for Green Function Experiment (non-missing)'
    figname = 'Qflux_app_patch_'+str(num_nnan)+'.png'
    map_plot(df0.lat,df0.lon,vqdp,projection="robinson",title=title,figname=figname)


if __name__ == '__main__':
    main()
