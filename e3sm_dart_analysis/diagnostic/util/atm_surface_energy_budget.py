"""Cached surface-energy-budget diagnostics for ensemble analysis and forecasts."""
from __future__ import annotations
import os,re
from pathlib import Path
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from util.bias_metrics import experiment_registry

INPUT_VARIABLES=("LHFLX","SHFLX","FSNS","FLNS")
OUTPUT_VARIABLES=("Rn","LE","SH","Residual","BowenRatio")
REGIONS=("global","land","ocean")

def _select_period(data,period):
    m=re.fullmatch(r"(\d{4})(\d{2})-(\d{4})(\d{2})",period)
    if not m: raise ValueError("period must be YYYYMM-YYYYMM")
    y0,m0,y1,m1=map(int,m.groups()); key=data.time.dt.year*100+data.time.dt.month
    return data.where((key>=y0*100+m0)&(key<=y1*100+m1),drop=True)

def _open_member_field(root,member,variable,period,run_segment):
    years=range(int(period[:4]),int(period[7:11])+1)
    choices=("6hourly_da","daily","6hourly") if run_segment=="da" else ("daily","6hourly")
    candidates=("LHFLX","QFLX") if variable=="LHFLX" else (variable,)
    for freq in choices:
        for candidate in candidates:
            folder=root/member/"archive/post/atm/180x360_aave/ts"/freq
            files=[folder/f"{candidate}.{member}.{year}.nc" for year in years]
            if all(p.is_file() for p in files):
                arrays=[]
                for p in files:
                    ds=xr.open_dataset(p)
                    try: arrays.append(ds[candidate].load())
                    finally: ds.close()
                field=xr.concat(arrays,dim="time")
                if candidate=="QFLX": field=field*2.5e6
                if freq.startswith("6hourly"): field=field.resample(time="1D").mean()
                return _select_period(field,period)
    raise FileNotFoundError(f"No {variable} files for {member} under {root}")

def energy_budget_cache_path(cache_dir,experiment,member,period):
    return Path(cache_dir)/experiment/f"surface_energy_budget_{member}_{period}.nc"

def cache_surface_energy_budget(path,*,data_root,experiment,season,run_segment,member,period=None,force_compute=False):
    path=Path(path)
    if path.is_file() and not force_compute: return xr.load_dataset(path)
    registry=experiment_registry(data_root,season,run_segment)
    if experiment not in registry: raise ValueError(f"{experiment} unavailable; choose {list(registry)}")
    meta=registry[experiment]; period=period or meta["period"]
    fields={v:_open_member_field(meta["path"],member,v,period,run_segment) for v in INPUT_VARIABLES}
    ds=xr.Dataset(fields); rn=ds.FSNS-ds.FLNS; residual=rn-(ds.LHFLX+ds.SHFLX)
    out=xr.Dataset({"Rn":rn,"LE":ds.LHFLX,"SH":ds.SHFLX,"Residual":residual,"BowenRatio":xr.where(np.abs(ds.LHFLX)>1e-6,ds.SHFLX/ds.LHFLX,np.nan)})
    out.attrs.update(experiment=experiment,member=member,season=season,run_segment=run_segment,period=period,input_variables=",".join(INPUT_VARIABLES))
    path.parent.mkdir(parents=True,exist_ok=True); tmp=path.with_suffix(path.suffix+".tmp"); out.to_netcdf(tmp); os.replace(tmp,path)
    return out

def _mask_for(field,region,landmask_file):
    if region=="global": return xr.ones_like(field.isel(time=0),dtype=bool)
    ds=xr.open_dataset(landmask_file)
    try:
        name=next((n for n in ("landfrac","LANDFRAC","landmask","mask") if n in ds),None)
        if name is None: name=next(iter(ds.data_vars))
        mask=ds[name].load().squeeze(drop=True)
    finally: ds.close()
    rename={}
    if "latitude" in mask.dims: rename["latitude"]="lat"
    if "longitude" in mask.dims: rename["longitude"]="lon"
    if rename: mask=mask.rename(rename)
    mask=mask.interp(lat=field.lat,lon=field.lon,method="nearest")
    return mask>=0.5 if region=="land" else mask<0.5

def regional_timeseries(datasets,variable,region,landmask_file):
    series=[]
    for member,ds in datasets.items():
        field=ds[variable]; mask=_mask_for(field,region,landmask_file)
        weights=np.cos(np.deg2rad(field.lat)).clip(min=0).broadcast_like(field.isel(time=0)).where(mask)
        ts=field.where(mask).weighted(weights).mean(("lat","lon"),skipna=True)
        series.append(ts.expand_dims(member=[member]))
    return xr.concat(series,dim="member")

def plot_energy_budget_timeseries(experiments,variable,region,landmask_file,label_map=None):
    fig,ax=plt.subplots(figsize=(11,4.5))
    for experiment,members in experiments.items():
        ts=regional_timeseries(members,variable,region,landmask_file)
        mean=ts.mean("member",skipna=True); lower=ts.quantile(.25,"member"); upper=ts.quantile(.75,"member")
        ax.plot(mean.time,mean,label=(label_map or {}).get(experiment,experiment)); ax.fill_between(mean.time.values,lower.values,upper.values,alpha=.18)
    ax.set_title(f"{region.title()} {variable}"); ax.grid(alpha=.3); ax.legend(); ax.set_ylabel("SH/LE" if variable=="BowenRatio" else "W m$^{-2}$"); fig.tight_layout(); return fig
