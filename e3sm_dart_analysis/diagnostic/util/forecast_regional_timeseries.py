"""Regional boxes and daily regional-mean ensemble series for atmospheric forecasts."""
from __future__ import annotations
import re
import numpy as np
import xarray as xr
from util.bias_metrics import experiment_registry
from util.model_daily import open_model_daily

REGIONS={
 "global":((-90.,90.),(-180.,180.)),
 "60S-60N":((-60.,60.),(-180.,180.)),
 "CONUS":((25.,50.),(-125.,-95.)),
 "NHMidLat":((25.,50.),(-60.,150.)),
 "Tropics":((-10.,10.),(-90.,60.)),
 "NINO3.4":((-5.,5.),(-170.,-120.)),
 "MAR_CONT":((-10.,10.),(95.,150.)),
}
VARIABLES={
 "PRECT":{"units":"mm day-1"},
 "FLUT":{"units":"W m-2"},
 "TS":{"units":"K"},
 "TREFHT":{"units":"K"},
 "PSL":{"units":"hPa"},
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

def _convert(field,variable):
 units=str(field.attrs.get("units","")).lower()
 if variable=="PRECT" and ("s-1" in units or units in {"m/s","m s-1"}): field=field*86400*1000
 if variable=="PSL" and "pa" in units and "hpa" not in units: field=field*.01
 field.attrs["units"]=VARIABLES[variable]["units"]; return field

def read_model_regional_series(data_root,experiment,season,run_segment,variable,region,period=None,max_members=None):
 """Daily regional-mean member series: ts/daily_mean, else ts/6hourly averaged to days (util.model_daily)."""
 registry=experiment_registry(data_root,season,run_segment)
 if experiment not in registry: raise ValueError(f"{experiment} unavailable; choose {list(registry)}")
 meta=registry[experiment]; period=period or meta["period"]; *_,years=_period(period)
 count=meta["nens"] if max_members is None else min(max_members,meta["nens"]); members=[]; opened=[]
 try:
  for number in range(1,count+1):
   member=f"EN{number:02d}"
   daily=open_model_daily(meta["path"]/member/"archive/post/atm/180x360_aave/ts",variable,member,years)
   field=_convert(_select_period(daily,period),variable).assign_coords(time=lambda x:x.time.dt.floor("D"))
   series=_regional_mean(field,region); series.attrs["units"]=VARIABLES[variable]["units"]
   members.append(series.expand_dims(ens=[member]))
  # One dask graph over all members so files are read in parallel.
  ensemble=xr.concat(members,dim="ens").compute()
 finally:
  for ds in opened: ds.close()
 return ensemble,meta
