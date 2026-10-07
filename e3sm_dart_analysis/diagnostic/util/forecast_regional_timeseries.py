"""Cached regional ensemble/observation time series for atmospheric forecasts."""
from __future__ import annotations
import os,re
from pathlib import Path
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from util.bias_metrics import experiment_registry

REGIONS={
 "global":((-90.,90.),(-180.,180.)),
 "CONUS":((25.,50.),(-125.,-95.)),
 "NHMidLat":((25.,50.),(-60.,150.)),
 "Tropics":((-10.,10.),(-90.,60.)),
 "NINO3.4":((-5.,5.),(-170.,-120.)),
 "MAR_CONT":((-10.,10.),(95.,150.)),
}
VARIABLES={
 "PRECT":{"observation":"GPM","obs_dir":"IMERG","obs_file":"PRECT.daily.{year}.nc","units":"mm day-1"},
 "FLUT":{"observation":"NOAA-OLR","obs_dir":"NOAA_IOLR","obs_file":"FLUT.daily.{year}.nc","units":"W m-2"},
 "TS":{"observation":"ERA5","units":"K"},
 "TREFHT":{"observation":"ERA5","units":"K"},
 "PSL":{"observation":"ERA5","units":"hPa"},
}

def _period(period):
 m=re.fullmatch(r"(\d{4})(\d{2})-(\d{4})(\d{2})",period)
 if not m: raise ValueError("period must be YYYYMM-YYYYMM")
 y0,m0,y1,m1=map(int,m.groups()); return y0,m0,y1,m1,range(y0,y1+1)

def _select_period(field,period):
 y0,m0,y1,m1,_=_period(period); key=field.time.dt.year*100+field.time.dt.month
 return field.where((key>=y0*100+m0)&(key<=y1*100+m1),drop=True)

def _normalize(field):
 rename={}
 if "latitude" in field.dims: rename["latitude"]="lat"
 if "longitude" in field.dims: rename["longitude"]="lon"
 if rename: field=field.rename(rename)
 if float(field.lon.max())>180: field=field.assign_coords(lon=((field.lon+180)%360)-180).sortby("lon")
 return field

def _subset(field,region):
 if region not in REGIONS: raise ValueError(f"Unknown region {region}; choose {list(REGIONS)}")
 field=_normalize(field); (la0,la1),(lo0,lo1)=REGIONS[region]
 lat_slice=slice(la0,la1) if float(field.lat[0])<float(field.lat[-1]) else slice(la1,la0)
 return field.sel(lat=lat_slice,lon=slice(lo0,lo1))

def _regional_mean(field,region):
 field=_subset(field,region); weights=np.cos(np.deg2rad(field.lat)).clip(min=0).broadcast_like(field.isel(time=0))
 return field.weighted(weights).mean(("lat","lon"),skipna=True)

def _convert(field,variable,is_obs=False):
 units=str(field.attrs.get("units","")).lower()
 if variable=="PRECT" and ("s-1" in units or units in {"m/s","m s-1"}): field=field*86400*1000
 if variable=="PSL" and "pa" in units and "hpa" not in units: field=field*.01
 field.attrs["units"]=VARIABLES[variable]["units"]; return field

def read_model_regional_series(data_root,experiment,season,run_segment,variable,region,period=None,max_members=None,frequency="daily"):
 registry=experiment_registry(data_root,season,run_segment)
 if experiment not in registry: raise ValueError(f"{experiment} unavailable; choose {list(registry)}")
 meta=registry[experiment]; period=period or meta["period"]; *_,years=_period(period)
 count=meta["nens"] if max_members is None else min(max_members,meta["nens"]); members=[]; opened=[]
 try:
  for number in range(1,count+1):
   member=f"EN{number:02d}"; folder=meta["path"]/member/"archive/post/atm/180x360_aave/ts"/frequency
   files=[folder/f"{variable}.{member}.{year}.nc" for year in years]; missing=[p for p in files if not p.is_file()]
   if missing: raise FileNotFoundError("Missing model files: "+", ".join(map(str,missing)))
   ds=xr.open_mfdataset(files,combine="by_coords",parallel=True,chunks={}); opened.append(ds)
   field=_convert(_select_period(ds[variable],period),variable).assign_coords(time=lambda x:x.time.dt.floor("D"))
   series=_regional_mean(field,region); series.attrs["units"]=VARIABLES[variable]["units"]
   members.append(series.expand_dims(ens=[member]))
  # One dask graph over all members so files are read in parallel.
  ensemble=xr.concat(members,dim="ens").compute()
 finally:
  for ds in opened: ds.close()
 return ensemble,meta

def read_observation_regional_series(reference_root,variable,region,period):
 cfg=VARIABLES[variable]; *_,years=_period(period); folder=Path(reference_root)/cfg["obs_dir"]/"daily"
 files=[folder/cfg["obs_file"].format(year=year) for year in years]
 for path in files:
  if not path.is_file(): raise FileNotFoundError(path)
 ds=xr.open_mfdataset(files,combine="by_coords",parallel=True,chunks={})
 try:
  name=next((x for x in (variable,"OLR","olr","precip") if x in ds),None)
  if name is None: raise KeyError(f"{variable} not found in {files[0]}: {list(ds.data_vars)}")
  field=_convert(_select_period(ds[name],period),variable,is_obs=True).assign_coords(time=lambda x:x.time.dt.floor("D"))
  series=_regional_mean(field,region).compute()
 finally: ds.close()
 series.attrs["units"]=VARIABLES[variable]["units"]
 return series

def timeseries_cache_path(cache_dir,variable,region,season,run_segment,period):
 return Path(cache_dir)/f"{variable}_{region}_{season}_{run_segment}_{period}_regional_timeseries.nc"

def cache_regional_timeseries(path,*,data_root=None,reference_root=None,experiments=None,season=None,run_segment=None,variable=None,region=None,period=None,max_members=None,force_compute=False):
 path=Path(path)
 if path.is_file() and not force_compute: return xr.load_dataset(path)
 if any(x is None for x in (data_root,reference_root,experiments,season,run_segment,variable,region,period)): raise ValueError("raw-data parameters are required for a missing cache")
 model=[]
 for experiment in experiments:
  field,_=read_model_regional_series(data_root,experiment,season,run_segment,variable,region,period,max_members=max_members)
  model.append(field.expand_dims(exp=[experiment]))
 model=xr.concat(model,dim="exp",join="outer"); obs=read_observation_regional_series(reference_root,variable,region,period)
 model,obs=xr.align(model,obs,join="inner"); ds=xr.Dataset({"model_member":model,"observation":obs})
 ds["model_mean"]=model.mean("ens",skipna=True); ds["model_q25"]=model.quantile(.25,"ens",skipna=True).reset_coords("quantile",drop=True); ds["model_q75"]=model.quantile(.75,"ens",skipna=True).reset_coords("quantile",drop=True)
 ds.attrs.update(variable=variable,region=region,season=season,run_segment=run_segment,period=period,observation=VARIABLES[variable]["observation"])
 path.parent.mkdir(parents=True,exist_ok=True); tmp=path.with_suffix(path.suffix+".tmp"); ds.to_netcdf(tmp); os.replace(tmp,path); return ds

def plot_regional_timeseries(ds,label_map=None,title=""):
 fig,ax=plt.subplots(figsize=(10,4.5)); lead=np.arange(ds.sizes["time"])
 for exp in ds.exp.values:
  mean=ds.model_mean.sel(exp=exp); lo=ds.model_q25.sel(exp=exp); hi=ds.model_q75.sel(exp=exp); label=(label_map or {}).get(str(exp),str(exp))
  ax.plot(lead,mean,label=label); ax.fill_between(lead,lo,hi,alpha=.16)
 ax.plot(lead,ds.observation,color="black",lw=2,label=ds.attrs.get("observation","Observation")); ax.set_xlabel("Lead day"); ax.set_ylabel(ds.model_member.attrs.get("units",VARIABLES[ds.attrs["variable"]]["units"])); ax.grid(alpha=.3); ax.legend(ncol=2)
 ax.set_title(title or f"{ds.attrs['variable']} — {ds.attrs['region']}"); fig.tight_layout(); return fig
