import numpy as np
import xarray as xr
from cartopy.mpl.gridliner import LatitudeFormatter, LongitudeFormatter
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import math
import datetime
import pandas as pd
import random

from joblib import Parallel, delayed
import multiprocessing
from time import sleep
from parfor import parfor

import geocat.comp as gc
import geocat.datafiles as gdf
import geocat.viz as gv

def read_beta_drift(fl,lats,lons):
    ds   = xr.open_dataset(fl, decode_times=False)
    row,col,ns = ds.beta_drift_mod.shape
    lonb  = lons+np.arange(row)*2.5
    latb  = lats+np.arange(col)*2.5
    spcb  = np.arange(0,ns,1)
    beta  = xr.DataArray(
            data=ds.beta_drift_mod,
            dims=["lon", "lat", "space"],
            coords=dict(lon=(["lon"], lonb),
                        lat=(["lat"], latb),
                        space=(["space"],spcb),),
            attrs=dict(description="Beta drift",
                       units="m/s",),
            )
    return beta

def read_sample_track(fl):
    ds    = xr.open_dataset(fl2, decode_times=False)
    sample = ds.sample
    nrdtrk,nspace= sample.shape
    #print(sample.shape)
    #print(sample)
    return sample,nrdtrk,nspace

def read_model_data(fl,varname,lats,latn,lone,lonw,syear,eyear,freq):
    #specify the index identificator:
    LatIndexer, LonIndexer = 'lat', 'lon'
    ds          = xr.open_dataset(fl, decode_times=True)
    date_range  = pd.date_range(datetime.datetime(syear, 1, 1, 0, 0), datetime.datetime(eyear, 12, 31, 23, 59), freq=freq)
    date_noleap = date_range[(date_range.day != 29) | (date_range.month != 2)]
    ds['time']  = pd.to_datetime(date_noleap)
    var         = ds[varname].loc[{LatIndexer: slice(lats, latn),
                                   LonIndexer: slice(lone, lonw)}]
    return var

def read_land_mask(fl,varname,lats,latn,lone,lonw):
    #specify the index identificator:
    LatIndexer, LonIndexer = 'lat', 'lon'
    ds  = xr.open_dataset(fl, decode_times=False)
    var = ds[varname].loc[{LatIndexer: slice(lats, latn),
                           LonIndexer: slice(lone, lonw)}]
    return var

def degree_per_meter(lat_input, dim_type):
    rlat = lat_input * math.pi /180
    if dim_type == 'lat':
        #Meters per degree Latitude:
        meter_per_degree = 111132.92 - 559.82 * math.cos(2* rlat) + 1.175*math.cos(4*rlat)
        degree_per_meter = 1.0 / meter_per_degree
    elif dim_type == 'lon':
        #Meters per degree Longitude: 
        meter_per_degree = 111412.84 * math.cos(rlat) - 93.5 * math.cos(3*rlat)
        degree_per_meter = 1.0 / meter_per_degree
    else:
        print("unknow dim_type: " + dim_type)
    return degree_per_meter

def yrdoy2date(year,doy):
    time_date = datetime.datetime(year,1,1,0,0) + datetime.timedelta(days=doy-1)
    #TimeTuple = time_date.timetuple()
    #for value in TimeTuple:
    #    print(value)
    return time_date

def sythetic_track_forecst(nvars,max_steps,time_step,sample,land_mask,rglat,rglon,year,u200,v200,u850,v850):

    strks        = np.empty((nvars,max_steps),dtype=float)
    stime        = np.empty((max_steps),dtype='U20')
    strks[:,:]   = np.nan
    stime[:]     = yrdoy2date(1,1)
    #print(stime)
    #print(strks)

    #initial condition
    step         = 0
    storm_lat    = sample[0].data
    storm_lon    = sample[1].data
    storm_time   = sample[2].data
    time_on_land = 0.0

    info = "working on location and time: lat = " \
            +str(storm_lat)+"; lon = " + str(storm_lon) + "; time = " + str(yrdoy2date(year,storm_time))
    print(info)

    lats = rglat[0] + 4
    latn = rglat[1] 
    lone = rglon[0]  
    lonw = rglon[1] - 5

    while ((storm_lat >= lats) and (storm_lat <= latn) and \
           (storm_lon >= lone) and (storm_lon <= lonw) and \
           (step < max_steps) and (time_on_land < 4.0)):

        #interpolate data to storm time
        #using xarray package; 'slinear' for spline interpolation
        time_out = datetime.datetime(year,1,1,0,0) + datetime.timedelta(days=storm_time*1.0) 
        tu200    = u200.interp(time=time_out,method='slinear') #or 'linear'
        tv200    = v200.interp(time=time_out,method='slinear') #or 'linear'
        tu850    = u850.interp(time=time_out,method='slinear') #or 'linear'
        tv850    = v850.interp(time=time_out,method='slinear') #or 'linear'
 
        #interpolate data to storm location 
        #using the geocat package
        #xu200  = gc.interp_multidim(tu200,storm_lat,storm_lon,tu200.lat,tu200.lon,cyclic=False,missing_val=np.nan)
        #xv200  = gc.interp_multidim(tv200,storm_lat,storm_lon,tv200.lat,tv200.lon,cyclic=False,missing_val=np.nan)
        #xu850  = gc.interp_multidim(tu850,storm_lat,storm_lon,tu850.lat,tu850.lon,cyclic=False,missing_val=np.nan)
        #xv850  = gc.interp_multidim(tv850,storm_lat,storm_lon,tv850.lat,tv850.lon,cyclic=False,missing_val=np.nan)

        #interpolate data to storm location 
        #using the xarray package
        xu200  = tu200.interp(lon=storm_lon,lat=storm_lat,method='nearest') 
        xv200  = tv200.interp(lon=storm_lon,lat=storm_lat,method='nearest')
        xu850  = tu850.interp(lon=storm_lon,lat=storm_lat,method='nearest')
        xv850  = tv850.interp(lon=storm_lon,lat=storm_lat,method='nearest')

        #extract beta_u, beta_v at current location
        #note: we use existing beta drift data on 2.5 degree
        #grid within Atlantic Basin constructed by Wenwei Xu
        #lat_index = math.ceil((lat-rglat[0])/2.5)
        #lon_index = math.ceil((lon-rglon[0])/2.5)
        #beta_u = beta_drift[lon_index,lat_index,1].data
        #beta_v = beta_drift[lon_index,lat_index,2].data
        #print(beta_drift)
        beta_uv = beta_drift.interp(lon=storm_lon,lat=storm_lat,method='nearest') 
        beta_u  = beta_uv[0]
        beta_v  = beta_uv[1]

        #calculate sythetic u, v
        alfa = 0.8;
        u = alfa * xu850 + (1-alfa) * xu200 + beta_u # apply beta correction
        v = alfa * xv850 + (1-alfa) * xv200 + beta_v # unit m/s %2.5
        
        #calculate new lat lon and convert from meter to degrees using conversion formula
        storm_lat = storm_lat + v * time_step * 24 * 60 * 60 * degree_per_meter(storm_lat, 'lat')
        storm_lon = storm_lon + u * time_step * 24 * 60 * 60 * degree_per_meter(storm_lat, 'lon')

        #count how long has the track being on land
        land_frac = land_mask.interp(lon=storm_lon, lat=storm_lat, method = 'nearest') 
        if land_frac.data >= 1.0: 
            time_on_land = time_on_land + time_step

        if step == 0:
            print('xu200=',xu200.data)
            print('xv200=',xv200.data)
            print('xu850=',xu850.data)
            print('xv850=',xv850.data)
            print('beta_u=',beta_u.data)
            print('beta_v=',beta_v.data)
            print('storm_lat=',storm_lat)
            print('storm_lon=',storm_lon)
            print('u=',u.data)
            print('v=',v.data)

        #save the data in matrix
        strks[:,step] = [storm_time, storm_lat, storm_lon, u.data, v.data, xu200.data, xv200.data, xu850.data, xv850.data]
        stime[step]   = yrdoy2date(year,storm_time)
        #print(strks[:,i,step])

        #print('step = ', step)
        #print('storm_time = ', storm_time)
        step = step + 1
        storm_time = storm_time + time_step

        del(time_out, u, v, beta_uv, beta_u, beta_v, tu200, tv200, tu850, tv850, xu200, xv200, xu850, xv850)

    return strks,stime

#Main function 
exps  = ["NDGUV"]
nexp  = len(exps)
syear = 2007
eyear = 2007
nyear = eyear - syear + 1 

#output frequency 
freq  = '6H'
if freq == '6H':
    time_step = 6.0/24.0 #frequency of output (6-hr)
elif freq == '3H':
    time_step = 3.0/24.0 #frequency of output (3-hr)
elif freq == '1H':
    time_step = 1.0/24.0 #frequency of output (1-hr)

max_steps = math.ceil(30/time_step) #Maximum time is 30 days, convert to time steps  
ntrk_nsub = 5000 

#TC basin and domain size
basin = "AL" #Atlantic 
rglat = [0,50]
rglon = [260,360]

#Read Beta drift data (from Wenwei Xu)
dir        = '/global/cscratch1/sd/zhan391/DARPA_project/Track_model/data'
fl1        = dir + '/beta_drift_linear_regression.nc'
beta_drift = read_beta_drift(fl1,rglat[0],rglon[0])
#print(beta)
#exit()

#Read sample tracks data (from Wenwei Xu)
fl2    = dir + "/sample_50000.nc"
sample,nrdtrk,nspace = read_sample_track(fl2)
#print(nrdtrk,nspace)
#exit()

#Read land mask data
vnm = "LANDFRAC"
fnm = dir + "/"+vnm+"_E3SM_180x360.nc"
land_mask = read_land_mask(fnm,vnm,rglat[0],rglat[1],rglon[0],rglon[1])

#loop each E3SM experiment for forecast 
for i in np.arange(nexp): 

    #loop each year for forecast 
    for j in np.arange(nyear): 

        year  = int(syear) + j 

        #Read model data
        vnms   = ['U200','V200','U850','V850']
        for vnm in vnms:
            fnm = dir + "/"+exps[i] + "/"+vnm+"_"+str(year)+".nc"
            if vnm == 'U200':
                u200  = read_model_data(fnm,vnm,rglat[0],rglat[1],rglon[0],rglon[1],syear,eyear,freq)
            elif vnm == 'V200':
                v200  = read_model_data(fnm,vnm,rglat[0],rglat[1],rglon[0],rglon[1],syear,eyear,freq)
            elif vnm == 'U850':
                u850  = read_model_data(fnm,vnm,rglat[0],rglat[1],rglon[0],rglon[1],syear,eyear,freq)
            elif vnm == 'V850':
                v850  = read_model_data(fnm,vnm,rglat[0],rglat[1],rglon[0],rglon[1],syear,eyear,freq)

        #select ntrk_nsub samples out of the nrdtrk sample tracks for each year
        track_list = np.arange(0,nrdtrk)
        sub_list   = random.choices(track_list, k=ntrk_nsub)

        # initialize array to track forecast data
        varlist               = ['time', 'lat', 'lon', 'u', 'v','u200','v200','u850','v850']
        nvars                 = len(varlist)
        sythetic_track        = np.empty((nvars,ntrk_nsub,max_steps),dtype=float)
        sythetic_track[:,:,:] = np.nan
        storm_time            = np.empty((ntrk_nsub,max_steps),dtype='U20')
        storm_time[:,:]       = yrdoy2date(1,1)
        #print(sythetic_track)
        #print(stime)

        #do track forecast
        for m in np.arange(0,ntrk_nsub):
            k              = sub_list[m]
            fsamp          = sample[k,:]
            fstrks,fstimes = \
                    sythetic_track_forecst(nvars,max_steps,time_step,fsamp,land_mask,\
                                           rglat,rglon,year,u200,v200,u850,v850)

            sythetic_track[:,m,:] = fstrks
            storm_time[m,:]       = fstimes

        #save data for analysis
        #strks,stime,sname
        nx,ny,nz = sythetic_track.shape
        ds = xr.Dataset({"strks": (("nvars","ntrks","nsteps"), sythetic_track),
                         "stime": (("ntrks","nsteps"),storm_time),
                         "sname": (("nvars"),varlist),
                         "sub_sample":(("ntrks"),sub_list)},
                        coords={
                            "nvars": np.arange(nx),
                            "ntrks": np.arange(ny),
                            "nsteps": np.arange(nz),
                            },
                        )

        file_out = 'sythetic_track_forecast_'+exps[i]+'_'+str(year)+'.nc'
        ds.to_netcdf(file_out)

exit()

