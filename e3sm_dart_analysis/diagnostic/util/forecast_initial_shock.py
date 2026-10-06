"""Initial-shock decay diagnostics from cached regional ensemble time series."""
from __future__ import annotations
import os
from pathlib import Path
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def _decay(t,a,tau,c): return a*np.exp(-t/tau)+c

def derive_shock_decay(shock,baseline,*,shock_experiment,baseline_experiment,variable,region):
 shock,baseline=xr.align(shock,baseline,join="inner",exclude={"ens"}); reference=baseline.mean("ens",skipna=True); signed=shock.mean("ens",skipna=True)-reference; magnitude=np.abs(signed)
 elapsed=((magnitude.time-magnitude.time.isel(time=0))/np.timedelta64(1,"D")).astype(float); t=np.asarray(elapsed); y=np.asarray(magnitude)
 valid=np.isfinite(t)&np.isfinite(y); params=np.array([np.nan,np.nan,np.nan])
 if valid.sum()>=5:
  guess=(max(float(y[valid][0]-y[valid][-1]),1e-8),max(float(t[valid][-1]-t[valid][0])/3,1.),max(float(y[valid][-1]),0.))
  try: params,_=curve_fit(_decay,t[valid],y[valid],p0=guess,bounds=([0,.05,0],[np.inf,1e3,np.inf]),maxfev=20000)
  except (RuntimeError,ValueError): pass
 fit=xr.full_like(magnitude,np.nan) if not np.isfinite(params).all() else xr.DataArray(_decay(t,*params),dims="time",coords={"time":magnitude.time})
 member_tau=[]
 for member in shock.ens.values:
  yy=np.abs(np.asarray(shock.sel(ens=member)-reference)); vv=np.isfinite(t)&np.isfinite(yy); tau=np.nan
  if vv.sum()>=5:
   try:
    p,_=curve_fit(_decay,t[vv],yy[vv],p0=(max(float(yy[vv][0]-yy[vv][-1]),1e-8),guess[1],max(float(yy[vv][-1]),0.)),bounds=([0,.05,0],[np.inf,1e3,np.inf]),maxfev=20000); tau=p[1]
   except (RuntimeError,ValueError): pass
  member_tau.append(tau)
 ds=xr.Dataset({"signed_anomaly":signed,"shock_magnitude":magnitude,"fitted_magnitude":fit,"member_tau_days":xr.DataArray(member_tau,dims="ens",coords={"ens":shock.ens})})
 ds["tau_days"]=xr.DataArray(params[1]); ds.attrs.update(shock_experiment=shock_experiment,baseline_experiment=baseline_experiment,variable=variable,region=region,fit="abs(shock ensemble mean - baseline ensemble mean) = a*exp(-t/tau)+c"); return ds

def shock_cache_path(cache_dir,variable,region,shock_experiment,baseline_experiment,period): return Path(cache_dir)/f"{variable}_{region}_{shock_experiment}_minus_{baseline_experiment}_{period}_shock_decay.nc"
def cache_shock_decay(path,*,shock=None,baseline=None,force_compute=False,**metadata):
 path=Path(path)
 if path.is_file() and not force_compute:return xr.load_dataset(path)
 if shock is None or baseline is None:raise ValueError("shock and baseline series required for missing cache")
 ds=derive_shock_decay(shock,baseline,**metadata); path.parent.mkdir(parents=True,exist_ok=True); tmp=path.with_suffix(path.suffix+".tmp"); ds.to_netcdf(tmp); os.replace(tmp,path); return ds

def plot_shock_decay(ds):
 fig,ax=plt.subplots(figsize=(8,4)); t=((ds.time-ds.time.isel(time=0))/np.timedelta64(1,"D")).astype(float); ax.plot(t,ds.shock_magnitude,label="|shock-baseline|",lw=2); ax.plot(t,ds.fitted_magnitude,"r--",label=f"fit: tau={float(ds.tau_days):.2f} days"); ax.axhline(0,color="gray",lw=.7); ax.set_xlabel("Lead day"); ax.set_ylabel(ds.attrs["variable"]+" anomaly"); ax.grid(alpha=.3); ax.legend(); ax.set_title(f"{ds.attrs['shock_experiment']} − {ds.attrs['baseline_experiment']} ({ds.attrs['region']})"); fig.tight_layout(); return fig
