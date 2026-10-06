"""Extended ensemble model-evaluation metrics built on the maintained bias readers."""
from __future__ import annotations

import os
from pathlib import Path
import numpy as np
import xarray as xr
import xskillscore as xs
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from util.bias_metrics import (
    ATM_VARIABLES, LND_VARIABLES, experiment_registry, read_model_ensemble,
    read_reference, derive_metric_dataset,
)

MAP_METRICS = (
    "mean", "reference_mean", "bias", "bias_stddev", "rmse", "spread",
    "crps", "coverage",
)
SCALAR_METRICS = (
    "model_global_mean", "reference_global_mean", "model_rms", "reference_rms",
    "rmse_global", "pattern_correlation", "bias_global", "spread_global",
    "crps_global", "coverage_global", "ssr", "bias_squared", "spread_squared",
    "mse", "sem",
)
OPTIONAL_METRICS = ("rank_histogram",)
METRIC_DESCRIPTIONS = {
    "bias": "ensemble-mean model minus reference",
    "rmse": "member RMSE against the reference",
    "spread": "ensemble standard deviation",
    "crps": "continuous ranked probability score",
    "coverage": "fraction of times the reference lies inside the ensemble envelope",
    "pattern_correlation": "cosine-latitude-weighted spatial correlation of climatologies",
    "ssr": "global spread divided by global RMSE",
    "rank_histogram": "reference rank among ensemble members over all valid samples",
}


def _weights(field):
    return np.cos(np.deg2rad(field.lat)).clip(min=0).broadcast_like(field)


def _weighted_mean(field):
    w = _weights(field)
    return field.weighted(w).mean(("lat", "lon"), skipna=True)


def derive_evaluation_dataset(model, reference, *, experiment, variable, component,
                              run_segment, season, period, n_bootstrap=0,
                              random_seed=42, include_rank_histogram=True):
    """Return cached map and scalar ensemble verification diagnostics."""
    model, reference = xr.align(model, reference, join="inner")
    if model.sizes.get("time", 0) == 0:
        raise ValueError("Model and reference have no overlapping times")
    base = derive_metric_dataset(
        model, reference, experiment=experiment, variable=variable,
        component=component, run_segment=run_segment, season=season,
        period=period, n_bootstrap=n_bootstrap, random_seed=random_seed,
    )
    ensemble_mean = model.mean("ens", skipna=True)
    model_clim = ensemble_mean.mean("time", skipna=True)
    reference_clim = reference.mean("time", skipna=True)
    difference = ensemble_mean - reference
    ensemble_spread = model.std("ens", skipna=True)
    coverage = ((reference >= model.min("ens", skipna=True)) &
                (reference <= model.max("ens", skipna=True))).astype("float32")
    base["mean_map"] = model_clim
    base["reference_mean_map"] = reference_clim
    base["bias_stddev_map"] = difference.std("time", skipna=True)
    base["coverage_map"] = coverage.mean("time", skipna=True)
    w = _weights(model_clim)
    model_global = model_clim.weighted(w).mean(("lat", "lon"), skipna=True)
    reference_global = reference_clim.weighted(w).mean(("lat", "lon"), skipna=True)
    model_rms = np.sqrt(((model_clim-model_global)**2).weighted(w).mean(("lat","lon"),skipna=True))
    reference_rms = np.sqrt(((reference_clim-reference_global)**2).weighted(w).mean(("lat","lon"),skipna=True))
    rmse_global = xs.rmse(reference_clim, model_clim, dim=["lat","lon"], weights=w, skipna=True)
    pcor = xs.pearson_r(reference_clim-reference_global, model_clim-model_global,
                        dim=["lat","lon"], weights=w, skipna=True)
    bias_global = _weighted_mean(base["bias_map"])
    spread_global = _weighted_mean(base["spread_map"])
    crps_global = _weighted_mean(base["crps_map"])
    coverage_global = _weighted_mean(base["coverage_map"])
    scalars = {
        "model_global_mean": model_global, "reference_global_mean": reference_global,
        "model_rms": model_rms, "reference_rms": reference_rms,
        "rmse_global": rmse_global, "pattern_correlation": pcor,
        "bias_global": bias_global, "spread_global": spread_global,
        "crps_global": crps_global, "coverage_global": coverage_global,
        "ssr": spread_global/rmse_global, "bias_squared": bias_global**2,
        "spread_squared": spread_global**2, "mse": rmse_global**2,
        "sem": spread_global/np.sqrt(model.sizes["ens"]),
    }
    for name,value in scalars.items(): base[name]=value
    if include_rank_histogram:
        ranks=(model < reference).sum("ens").where(reference.notnull())
        values=np.asarray(ranks.values).ravel(); values=values[np.isfinite(values)].astype(int)
        counts=np.bincount(values,minlength=model.sizes["ens"]+1).astype("float64")
        if counts.sum(): counts/=counts.sum()
        base["rank_histogram"] = xr.DataArray(counts,dims="rank",coords={"rank":np.arange(counts.size)})
    base.attrs["available_map_metrics"] = ",".join(MAP_METRICS)
    base.attrs["available_scalar_metrics"] = ",".join(SCALAR_METRICS)
    return base


def evaluation_cache_path(cache_dir, experiment, variable, period):
    return Path(cache_dir)/f"{experiment}_{variable}_{period}_model_evaluation.nc"


def cache_evaluation(path, *, model=None, reference=None, force_compute=False, **metadata):
    path=Path(path)
    if path.is_file() and not force_compute:
        cached=xr.load_dataset(path)
        requested=int(metadata.get("n_bootstrap",0))
        actual=int(cached.attrs.get("bootstrap_samples",0))
        if requested != actual:
            raise ValueError(f"Cached bootstrap_samples={actual}, requested={requested}; set force_compute=True")
        return cached
    if model is None or reference is None:
        raise ValueError("model and reference are required for a missing cache")
    ds=derive_evaluation_dataset(model,reference,**metadata)
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+".tmp")
    ds.to_netcdf(temporary); os.replace(temporary,path)
    return ds


def plot_evaluation_maps(datasets, metric, variable, *, label_map=None, title=""):
    key=f"{metric}_map"
    absent=[name for name,ds in datasets.items() if key not in ds]
    if absent: raise KeyError(f"{key} missing for {absent}; map metrics are {MAP_METRICS}")
    fields={name:ds[key] for name,ds in datasets.items()}
    finite=np.concatenate([x.values[np.isfinite(x.values)] for x in fields.values()])
    if finite.size==0: raise ValueError(f"No finite values for {variable} {metric}")
    diverging=metric in {"bias","bias_stddev"}
    if diverging:
        limit=float(np.nanpercentile(np.abs(finite),98)); vmin,vmax,cmap=-limit,limit,"RdBu_r"
    else:
        vmin,vmax,cmap=0.0,float(np.nanpercentile(finite,98)),"YlGnBu"
    n=len(fields); cols=min(3,n); rows=int(np.ceil(n/cols)); projection=ccrs.PlateCarree()
    fig,axes=plt.subplots(rows,cols,figsize=(5*cols,3.8*rows),subplot_kw={"projection":projection},squeeze=False)
    image=None
    for i,(name,field) in enumerate(fields.items()):
        ax=axes.flat[i]; image=ax.pcolormesh(field.lon,field.lat,field,shading="auto",transform=projection,cmap=cmap,vmin=vmin,vmax=vmax)
        if metric=="bias" and {"bias_ci_lower","bias_ci_upper"} <= set(datasets[name]):
            sig=xr.where((datasets[name].bias_ci_lower>0)|(datasets[name].bias_ci_upper<0),1.0,np.nan)
            ax.contourf(sig.lon,sig.lat,sig,levels=[0.5,1.5],colors="none",hatches=[".."],transform=projection)
        ax.coastlines(linewidth=.6); ax.set_global(); ax.set_title((label_map or {}).get(name,name))
    for i in range(n,rows*cols): axes.flat[i].axis("off")
    if image is not None: fig.colorbar(image,ax=axes.ravel().tolist(),shrink=.75,label=f"{variable} {metric}")
    if title: fig.suptitle(title)
    return fig
