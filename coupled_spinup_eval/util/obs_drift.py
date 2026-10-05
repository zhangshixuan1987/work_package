"""
Drift of the global-mean surface climate relative to an observed pre-industrial reference.

* Observations: annual global means (cos-latitude weighted) of the e3sm_diags observational time
  series (e.g. NOAA-20C ``<var>_183601_201512.nc``) over a reference period such as 1871-1890.
* Model: annual global means from ``post/time_series/atm`` area means (util.area_mean_drift).
* Metrics over a sliding window of the reference length: mean bias (model window mean - observed
  mean) and the RMSE between the window's annual means and the observed annual means (bias and
  interannual mismatch together).
"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from util.area_mean_drift import DPM


def obs_annual_global_mean(path, var, years, scale=1.0, offset=0.0):
    """Annual (days-in-month weighted) cos(lat)-weighted global mean of ``var`` over years (y0, y1)."""
    with xr.open_dataset(path) as ds:
        da = ds[var].sel(time=slice(f"{years[0]:04d}-01-01", f"{years[1]:04d}-12-31")).load()
    w = np.cos(np.deg2rad(da["lat"]))
    gm = da.weighted(w.fillna(0)).mean(("lat", "lon"))
    dpm = gm["time"].dt.days_in_month
    ann = (gm * dpm).groupby("time.year").sum() / dpm.groupby("time.year").sum()
    return ann * scale + offset


def model_annual_global_mean(area_means, var):
    """Annual (days-in-month weighted) global mean from a read_area_means() DataArray -> (year,)."""
    x = area_means.sel(variable=var, region="Global")
    w = xr.DataArray(DPM, dims="month", coords={"month": x["month"]})
    return (x * w).sum("month", skipna=False) / w.sum()


def sliding_drift(model, obs, window=None):
    """Sliding-window bias and RMSE of model annual means against the observed annual means.

    model: (year,) annual means; obs: (year,) annual means of the reference period.
    The window length is len(obs) unless given; the result is labelled by the window's centre year."""
    o = np.asarray(obs, float)
    n = window or o.size
    o = o[:n]
    m = np.asarray(model, float)
    years = np.asarray(model["year"])
    k = m.size - n + 1
    if k < 1:
        raise ValueError(f"model series ({m.size} yr) shorter than the window ({n} yr)")
    win = np.lib.stride_tricks.sliding_window_view(m, n)          # (k, n)
    bias = win.mean(1) - o.mean()
    rmse = np.sqrt(np.mean((win - o[None, :]) ** 2, axis=1))
    centre = years[: k] + (n - 1) / 2
    return xr.Dataset({"bias": ("year", bias), "rmse": ("year", rmse)}, coords={"year": centre})


def plot_drift(results, fields, *, metric="bias", obs_std=None, markers=(), ncols=2, figsize=None,
               fontz=12.0, panel_letters=True, save=None):
    """
    One panel per field: the sliding-window ``metric`` ("bias" or "rmse") against model year.

    results: {field: {label: dict(ds=sliding_drift() Dataset, color=..., lw=...)}}
    fields:  {field: (unit, title)}
    obs_std: optional {field: interannual std of the observed annual means}; bias panels get a
             ±2σ/√n band (sampling uncertainty of the observed mean) and a zero line.
    markers: model years marked with a vertical dotted line (e.g. recoupling, end of spin-up).
    """
    names = list(fields)
    nrows = int(np.ceil(len(names) / ncols))
    figsize = figsize or (7.0 * ncols, 2.8 * nrows)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False, sharex=True, layout="constrained")
    letters = iter("abcdefghijklmnopqrstuvwxyz")
    for ax, f in zip(axes.flat, names):
        for lab, r in results[f].items():
            ds = r["ds"]
            ax.plot(ds["year"], ds[metric], color=r["color"], lw=r.get("lw", 1.2), label=lab)
        if metric == "bias":
            ax.axhline(0, color="k", lw=0.8)
            if obs_std is not None and f in obs_std:
                ax.axhspan(-2 * obs_std[f]["sem"], 2 * obs_std[f]["sem"], color="0.85", zorder=0)
        for yr in markers:
            ax.axvline(yr, color="0.4", ls=":", lw=1)
        unit, title = fields[f]
        ax.set_title((f"({next(letters)}) " if panel_letters else "") + title, fontsize=fontz, loc="left")
        ax.set_ylabel(f"{'bias' if metric == 'bias' else 'RMSE'} ({unit})", fontsize=0.9 * fontz)
        ax.tick_params(labelsize=0.85 * fontz)
        ax.grid(alpha=0.3)
    for ax in axes.flat[len(names):]:
        ax.set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel("Model year (window centre)", fontsize=0.9 * fontz)
    h, l = axes.flat[0].get_legend_handles_labels()
    fig.legend(h, l, loc="outside lower center", ncol=len(l), fontsize=fontz, frameon=False)
    if save:
        fig.savefig(save, dpi=300, bbox_inches="tight")
    return fig, axes
