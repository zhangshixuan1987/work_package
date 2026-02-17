import numpy as np
import xarray as xr
from cartopy.mpl.gridliner import LatitudeFormatter, LongitudeFormatter
from scipy import interpolate
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import math
import datetime
import pandas as pd

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

def sythetic_track_forecst(nrdtrk,max_steps,sample,land_mask,rglat,rglon,year,u200,v200,u850,v850):
    sname        = ['time', 'lat', 'lon', 'u', 'v','u200','v200','u850','v850']
    nvars        = len(sname)

    strks        = np.empty((nvars,nrdtrk,max_steps),dtype=float)
    strks[:,:,:] = np.nan
    #print(strks)

    stime        = np.empty((nrdtrk,max_steps),dtype='U20')  
    stime[:,:]   = yrdoy2date(1,1) 
    #print(stime)

    for i in np.arange(nrdtrk):
        #initial condition
        step     = 0
        lat      = sample[i,0].data
        lon      = sample[i,1].data
        time     = sample[i,2].data
        time_on_land = 0.0

        info = "working on location and time: lat = " \
                +str(lat)+"; lon = " + str(lon) + "; time = " + str(yrdoy2date(year,time))
        print(info)

        #time_on_land = 0
        #temp_sythetic_tracks = nan (1,max_steps, 9)
        while ((lat > (rglat[0]+4)) and (lat < rglat[1]) and \
               (lon > rglon[0]) and (lon < (rglon[1]-5)) and \
               (step < max_steps) and (time_on_land < 4.0)):

            #interpolate data to storm location 
            tu200  = gc.interp_multidim(u200,lat, lon, u200.lat,u200.lon,cyclic=False, missing_val=np.nan)
            tv200  = gc.interp_multidim(v200,lat, lon, v200.lat,v200.lon,cyclic=False, missing_val=np.nan)
            tu850  = gc.interp_multidim(u850,lat, lon, u850.lat,u850.lon,cyclic=False, missing_val=np.nan)
            tv850  = gc.interp_multidim(v850,lat, lon, v850.lat,v850.lon,cyclic=False, missing_val=np.nan)

            #interpolate data to storm time 

            time_out = pd.to_datetime(yrdoy2date(year,time))
            fu200    = gc.linint1(tu200,time_out,tu200.)
            exit() 
            fu200    = tu200.interp(time=time_out,method = 'slinear') #or 'linear'
            fv200    = tv200.interp(time=time_out,method = 'slinear') #or 'linear'
            fu850    = tu850.interp(time=time_out,method = 'slinear') #or 'linear'
            fv850    = tv850.interp(time=time_out,method = 'slinear') #or 'linear'

            #extract beta_u, beta_v at current location
            #note: we use existing beta drift data on 2.5 degree
            #grid within Atlantic Basin constructed by Wenwei Xu
            #lat_index = math.ceil((lat-rglat[0])/2.5)
            #lon_index = math.ceil((lon-rglon[0])/2.5)
            #beta_u = beta_drift[lon_index,lat_index,1].data
            #beta_v = beta_drift[lon_index,lat_index,2].data
            print(beta_drift)
            beta_uv = beta_drift.interp(lon=lon, lat=lat,method = 'nearest') 
            beta_u  = beta_uv[0].data
            beta_v  = beta_uv[1].data

            #calculate sythetic u, v
            alfa = 0.8;
            u = alfa * fu850 + (1-alfa) * fu200 + beta_u # apply beta correction
            v = alfa * fv850 + (1-alfa) * fv200 + beta_v # unit m/s %2.5

            #calculate new lat lon and convert from meter to degrees using conversion formula
            lat = lat + v * time_step * 24 * 60 * 60 * degree_per_meter(lat, 'lat')
            lon = lon + u * time_step * 24 * 60 * 60 * degree_per_meter(lat, 'lon')

            #count how long has the track being on land
            land_frac = land_mask.interp(lon=lon, lat=lat, method = 'nearest') 
            if land_frac.data >= 1.0: 
                time_on_land = time_on_land + time_step

            #save the data in matrix 
            strks[:,i,step] = [time, lat, lon, u, v, u200, v200, u850, v850]
            stime[i,step]   = yrdoy2date(year,time)

            step = step + 1
            time = time + time_step

    return strks,stime,sname

#Main function 
exps  = ["CLIM","NDGUV"]
nexp  = len(exps)
syear = 2007
eyear = 2017
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
#print(sample.shape)
#print(sample)
#print(nrdtrk)
#print(nspace)
#exit()

print(min(sample[:,2].data))
print(max(sample[:,2].data))

for i in np.arange(nexp): 

    #Read model data 
    vnms   = ['U200','V200','U850','V850']
    for vnm in vnms:
        fnm = dir + "/"+exps[i] + "/"+vnm+"_"+str(syear)+"-"+str(eyear)+".nc"  
        if vnm == 'U200':
            u200  = read_model_data(fnm,vnm,rglat[0],rglat[1],rglon[0],rglon[1],syear,eyear,freq)
        elif vnm == 'V200':
            v200  = read_model_data(fnm,vnm,rglat[0],rglat[1],rglon[0],rglon[1],syear,eyear,freq)
        elif vnm == 'U850':
            u850  = read_model_data(fnm,vnm,rglat[0],rglat[1],rglon[0],rglon[1],syear,eyear,freq)
        elif vnm == 'V850':
            v850  = read_model_data(fnm,vnm,rglat[0],rglat[1],rglon[0],rglon[1],syear,eyear,freq)

    #Read land mask data
    vnm = "LANDFRAC"
    fnm = dir + "/"+vnm+"_E3SM_180x360.nc"
    land_mask = read_land_mask(fnm,vnm,rglat[0],rglat[1],rglon[0],rglon[1])

    for j in np.arange(nyear): 

        year  = int(syear + j) 

        #select data at a specific year 
        fu200 = u200.sel(time=str(year))
        fv200 = v200.sel(time=str(year))
        fu850 = u850.sel(time=str(year))
        fv850 = v850.sel(time=str(year))
        #print(fu200)

        #do track forecast 
        sythetic_track,vartime,varlist = sythetic_track_forecst(nrdtrk,max_steps,sample,\
                                                                land_mask,rglat,rglon,year,\
                                                                fu200,fv200,fu850,fv850)

        #save data for analysis 
        #strks,stime,sname

exit()

