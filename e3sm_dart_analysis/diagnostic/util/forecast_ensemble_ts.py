"""Hindcast ensemble time series and initial-shock diagnostics.

Each hindcast is one continuous forecast from the analysis time. Model records
carry ``time_bnds``: the first record has zero width and holds the initial
state; later records are interval means (``time: mean``) or end-of-interval
values (``time: point``). Every record is placed at its lead time (interval
midpoint for means, interval end for point values) and, where available,
verified against ERA5 6-hourly samples averaged over the same interval.
"""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from util.bias_metrics import _period_years
from util.forecast_regional_timeseries import REGIONS as BOX_REGIONS

FREQUENCIES = ("6hourly", "daily")
REGIONS = ("global", "land", "ocean") + tuple(r for r in BOX_REGIONS if r != "global")
# Variables with an ERA5 6-hourly reference on the model's 180x360 grid.
REFERENCE_VARIABLES = ("PSL", "TS", "TREFHT", "PRECT")
_HPA = {"PS", "PSL"}
_PRECIPITATION = {"PRECT", "PRECC", "PRECL"}

DIAGNOSTICS = {
    "member_mean": "regional mean per member",
    "member_min": "regional minimum per member",
    "member_max": "regional maximum per member",
    "member_tendency": "regional mean of |dX/dt| per member (units per day)",
    "ensemble_spread": "regional mean of the gridpoint ensemble standard deviation",
    "reference_mean": "ERA5 regional mean over the same intervals",
    "reference_tendency": "ERA5 regional mean of |dX/dt| (units per day)",
    "ensemble_mean_rmse": "gridpoint RMSE of the ensemble mean against ERA5",
}
SUMMARY_METRICS = {
    "tendency_ratio": "mean |dX/dt| over the shock window / over the settled window",
    "tendency_ratio_era5": "mean |dX/dt| over the shock window / ERA5 over the same leads",
    "adjustment_days": "last lead before the settled window where the 1-day running |dX/dt| exceeds (1 + tolerance) x the settled level",
    "initial_bias": "ensemble-mean minus ERA5 regional mean at lead 0",
    "drift": "bias change from lead 0 to the last day of the shock window",
    "spread_ratio": "ensemble spread at the end of the shock window / at lead 0",
}


def _units(field, variable):
    attrs = dict(field.attrs)
    units = str(attrs.get("units", ""))
    if variable in _HPA and units == "Pa":
        field, units = field * 0.01, "hPa"
    elif variable in _PRECIPITATION and units in {"m/s", "m s-1"}:
        field, units = field * 86400 * 1000, "mm day-1"
    elif variable in _PRECIPITATION and units in {"kg m**-2 s**-1", "kg m-2 s-1", "mm/s"}:
        field, units = field * 86400, "mm day-1"
    field.attrs = {**attrs, "units": units}
    return field


def _normalize(field):
    rename = {k: v for k, v in (("latitude", "lat"), ("longitude", "lon")) if k in field.dims}
    return field.rename(rename) if rename else field


def region_masks(lat, lon, regions, landmask_file):
    """Boolean (region, lat, lon) masks; box regions use -180..180 longitudes."""
    lon180 = ((lon + 180) % 360) - 180
    land = None
    masks = []
    for region in regions:
        if region in ("land", "ocean"):
            if land is None:
                with xr.open_dataset(landmask_file) as ds:
                    land = _normalize(ds["landfrac"].load()).interp(lat=lat, lon=lon, method="nearest") >= 0.5
            mask = land if region == "land" else ~land
        elif region in BOX_REGIONS:
            (la0, la1), (lo0, lo1) = BOX_REGIONS[region]
            mask = (lat >= la0) & (lat <= la1) & (lon180 >= lo0) & (lon180 <= lo1)
        else:
            raise ValueError(f"Unknown region {region!r}; choose from {list(REGIONS)}")
        masks.append(mask.broadcast_like(lat * lon))
    return xr.concat(masks, dim=pd.Index(list(regions), name="region"))


def _lead_coordinates(ds, variable, period):
    """Record selection and lead/interval coordinates from ``time_bnds``."""
    bounds = ds[ds.time.attrs.get("bounds", "time_bnds")].load()
    start, end = bounds.isel(nbnd=0), bounds.isel(nbnd=1)
    y0, m0, y1, m1, _ = _period_years(period)
    key = start.dt.year * 100 + start.dt.month
    keep = ((key >= y0 * 100 + m0) & (key <= y1 * 100 + m1)).values
    start, end = start.values[keep], end.values[keep]
    point = "time: point" in str(ds[variable].attrs.get("cell_methods", ""))
    position = end if point else start + (end - start) / 2
    lead = (position - start[0]) / np.timedelta64(1, "h")
    return keep, lead.astype(float), start, end


def read_member_ensemble(config, variable, period, *, frequency="6hourly", max_members=None):
    """Lazy (member, lead, lat, lon) field for one experiment, in analysis units."""
    if frequency not in FREQUENCIES:
        raise ValueError(f"frequency must be one of {FREQUENCIES}")
    count = config["nens"] if max_members is None else min(max_members, config["nens"])
    *_, years = _period_years(period)
    fields, opened, reference_lead = [], [], None
    try:
        for number in range(1, count + 1):
            member = f"EN{number:02d}"
            folder = config["path"] / member / "archive/post/atm/180x360_aave/ts" / frequency
            files = [folder / f"{variable}.{member}.{year}.nc" for year in years]
            missing = [str(p) for p in files if not p.is_file()]
            if missing:
                raise FileNotFoundError("Missing model files: " + ", ".join(missing))
            ds = xr.open_mfdataset(files, combine="by_coords", chunks={"time": 16})
            opened.append(ds)
            keep, lead, start, end = _lead_coordinates(ds, variable, period)
            if reference_lead is None:
                reference_lead, interval = lead, (start, end)
            elif not np.array_equal(lead, reference_lead):
                raise ValueError(f"{member} lead times differ from EN01")
            field = _normalize(ds[variable].isel(time=keep)).astype("float32")
            field = field.drop_vars("time").rename(time="lead").assign_coords(lead=lead)
            fields.append(field.expand_dims(member=[member]))
    except Exception:
        for ds in opened:
            ds.close()
        raise
    field = _units(xr.concat(fields, dim="member"), variable)
    field = field.assign_coords(interval_start=("lead", interval[0]), interval_end=("lead", interval[1]))
    field.attrs.update(variable=variable, frequency=frequency)
    return field, opened


def read_era5_reference(reference_root, variable, field):
    """ERA5 6-hourly samples averaged (trapezoid) over each model interval.

    Zero-width records take the sample at that time. For 6-hourly model means
    the reference is the mean of the two bounding samples, an approximation of
    the interval mean.
    """
    start = pd.DatetimeIndex(field.interval_start.values)
    end = pd.DatetimeIndex(field.interval_end.values)
    folder = Path(reference_root) / "ERA5" / "6hourly"
    files = [folder / f"ERA5.6hourly.en00.{variable}.{year}01-{year}12.nc"
             for year in range(start.min().year, end.max().year + 1)]
    missing = [str(p) for p in files if not p.is_file()]
    if missing:
        raise FileNotFoundError("Missing ERA5 files: " + ", ".join(missing))
    with xr.open_mfdataset(files, combine="by_coords") as ds:
        era5 = _normalize(ds[variable]).sel(time=slice(start.min(), end.max())).load()
    era5 = _units(era5, variable)
    times = pd.DatetimeIndex(era5.time.values)
    weights = np.zeros((len(start), len(times)))
    for row, (a, b) in enumerate(zip(start, end)):
        inside = (times >= a) & (times <= b)
        if a == b or "time: point" in str(field.attrs.get("cell_methods", "")):
            inside = times == b
        weights[row, inside] = 1.0
        if a != b:
            weights[row, times == a] *= 0.5
            weights[row, times == b] *= 0.5
        if weights[row].sum() == 0:
            weights[row] = np.nan
    weights /= weights.sum(axis=1, keepdims=True)
    weights = xr.DataArray(weights, dims=("lead", "time"), coords={"lead": field.lead})
    reference = xr.dot(weights, era5.fillna(0), dim="time").where(era5.notnull().all("time"))
    reference = reference.assign_coords(lat=field.lat, lon=field.lon)
    reference.attrs["units"] = era5.attrs["units"]
    return reference


def _area_mean(values, weights):
    return values.weighted(weights).mean(("lat", "lon"), skipna=True)


def _abs_tendency(values):
    step = np.diff(values.lead.values) / 24.0
    change = np.abs(values.diff("lead")) / xr.DataArray(step, dims="lead", coords={"lead": values.lead.values[1:]})
    return change.reindex(lead=values.lead)


def derive_ensemble_ts(field, masks, reference=None):
    """Regional member, spread, tendency, and (optional) ERA5 diagnostics."""
    intervals = {name: field[name] for name in ("interval_start", "interval_end")}
    field = field.drop_vars(list(intervals))
    if reference is not None:
        reference = reference.drop_vars(list(intervals), errors="ignore")
    weights = np.cos(np.deg2rad(field.lat)).clip(min=0) * masks
    masked = field.where(masks)
    out = xr.Dataset({
        "member_mean": _area_mean(field, weights),
        "member_min": masked.min(("lat", "lon")),
        "member_max": masked.max(("lat", "lon")),
        "member_tendency": _area_mean(_abs_tendency(field), weights),
        "ensemble_spread": _area_mean(field.std("member", ddof=1), weights),
    })
    if reference is not None:
        error = field.mean("member") - reference
        out["reference_mean"] = _area_mean(reference, weights)
        out["reference_tendency"] = _area_mean(_abs_tendency(reference), weights)
        out["ensemble_mean_rmse"] = np.sqrt(_area_mean(error ** 2, weights))
    out = out.compute().assign_coords(intervals)
    for name in out.data_vars:
        units = field.attrs["units"]
        out[name].attrs["units"] = f"{units} day-1" if name.endswith("tendency") else units
        out[name].attrs["long_name"] = DIAGNOSTICS[name]
    return out


def ensemble_ts_cache_path(cache_dir, experiment, variable, frequency, period):
    return Path(cache_dir) / f"{experiment}_{variable}_{frequency}_{period}_ensemble_ts.nc"


def ensemble_ts_cache_status(path, *, regions, members, use_reference):
    """'' when the cache covers the request, otherwise why it must be rebuilt."""
    path = Path(path)
    if not path.is_file():
        return "missing"
    with xr.open_dataset(path) as cached:
        missing = set(regions) - set(cached.region.values.tolist())
        if missing:
            return f"regions {sorted(missing)} not cached"
        if cached.sizes["member"] != members:
            return f"members: cached {cached.sizes['member']}, requested {members}"
        if use_reference and "reference_mean" not in cached:
            return "no ERA5 reference cached"
    return ""


def cache_ensemble_ts(path, *, config, experiment, variable, period, frequency, regions,
                      landmask_file, reference_root=None, max_members=None, force_compute=False,
                      log=print):
    """Compute or reuse one experiment/variable cache covering ``regions``."""
    path = Path(path)
    members = config["nens"] if max_members is None else min(max_members, config["nens"])
    use_reference = reference_root is not None and variable in REFERENCE_VARIABLES
    reason = "forced" if force_compute else ensemble_ts_cache_status(
        path, regions=regions, members=members, use_reference=use_reference)
    if not reason:
        log(f"[REUSE] {path.name}")
        return xr.load_dataset(path).sel(region=list(regions))
    field, opened = read_member_ensemble(config, variable, period, frequency=frequency,
                                         max_members=max_members)
    try:
        masks = region_masks(field.lat, field.lon, regions, landmask_file)
        reference = read_era5_reference(reference_root, variable, field) if use_reference else None
        ds = derive_ensemble_ts(field, masks, reference)
    finally:
        for opened_ds in opened:
            opened_ds.close()
    ds.attrs.update(experiment=experiment, label=config["label"], variable=variable,
                    period=period, frequency=frequency, ensemble_members=members,
                    reference="ERA5 6-hourly" if use_reference else "",
                    lead_definition="hours from the initial time to the interval midpoint "
                                    "(time: mean) or interval end (time: point)")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".nc.tmp")
    ds.to_netcdf(temporary)
    os.replace(temporary, path)
    log(f"[WRITE] {path.name} ({reason})")
    return ds


def _running_daily_mean(series):
    lead = series.lead.values
    per_day = max(1, int(round(24.0 / np.median(np.diff(lead[1:])))))
    return series.rolling(lead=per_day, center=True, min_periods=1).mean()


def shock_summary(datasets, *, shock_days=2.0, settled_days=(10.0, 30.0), tolerance=0.2):
    """One row per experiment and region with the ``SUMMARY_METRICS``."""
    rows = []
    for experiment, ds in datasets.items():
        days = ds.lead / 24.0
        shock = (days > 0) & (days <= shock_days)
        settled = (days >= settled_days[0]) & (days <= settled_days[1])
        shock_end = (days > shock_days - 1) & (days <= shock_days)
        for region in ds.region.values:
            sel = ds.sel(region=region)
            tendency = sel.member_tendency.mean("member")
            settled_level = float(tendency.where(settled).mean())
            row = {"experiment": experiment, "label": ds.attrs.get("label", experiment),
                   "variable": ds.attrs["variable"], "region": str(region),
                   "units": sel.member_mean.attrs.get("units", ""),
                   "tendency_ratio": float(tendency.where(shock).mean()) / settled_level}
            smooth = _running_daily_mean(tendency)
            above = smooth.where((days > 0) & (days < settled_days[0])) > (1 + tolerance) * settled_level
            exceed = days.values[above.values]
            row["adjustment_days"] = float(exceed.max()) if exceed.size else 0.0
            spread = sel.ensemble_spread
            row["spread_ratio"] = float(spread.where(shock_end).mean()) / float(spread.isel(lead=0))
            if "reference_mean" in sel:
                row["tendency_ratio_era5"] = (float(tendency.where(shock).mean())
                                              / float(sel.reference_tendency.where(shock).mean()))
                bias = sel.member_mean.mean("member") - sel.reference_mean
                row["initial_bias"] = float(bias.isel(lead=0))
                row["drift"] = float(bias.where(shock_end).mean()) - row["initial_bias"]
            rows.append(row)
    columns = ["experiment", "label", "variable", "region", "units", *SUMMARY_METRICS]
    return pd.DataFrame(rows).reindex(columns=columns)


def _ylabel(ds, name):
    return f"{ds.attrs['variable']} ({ds[name].attrs.get('units', '')})"


def plot_ensemble_timeseries(datasets, region, *, labels, colors, zoom_days=10, show_members=True,
                             title=""):
    """Members, ensemble mean, and ERA5 over the full forecast and the first ``zoom_days``."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.5), gridspec_kw={"width_ratios": (3, 2)})
    reference = None
    for experiment, ds in datasets.items():
        sel = ds.sel(region=region)
        days = sel.lead / 24.0
        for ax in axes:
            if show_members:
                for member in sel.member.values:
                    ax.plot(days, sel.member_mean.sel(member=member), color=colors[experiment],
                            lw=0.5, alpha=0.3)
            ax.plot(days, sel.member_mean.mean("member"), color=colors[experiment], lw=2,
                    label=labels[experiment])
        if reference is None and "reference_mean" in sel:
            reference = sel.reference_mean
    for ax in axes:
        if reference is not None:
            ax.plot(reference.lead / 24.0, reference, color="black", lw=1.5, ls="--", label="ERA5")
        ax.set_xlabel("Lead (days)")
        ax.grid(alpha=0.3)
    first = next(iter(datasets.values()))
    axes[0].set_ylabel(_ylabel(first, "member_mean"))
    axes[1].set_xlim(0, zoom_days)
    axes[1].set_title(f"First {zoom_days:g} days")
    axes[0].legend(fontsize=8, ncol=2)
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_shock_diagnostics(datasets, region, *, labels, colors, zoom_days=10, title=""):
    """Tendency, drift, RMSE/spread, and spatial range over the first ``zoom_days``."""
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
    (ax_tend, ax_bias), (ax_err, ax_range) = axes
    reference_drawn = False
    for experiment, ds in datasets.items():
        sel = ds.sel(region=region)
        days = sel.lead / 24.0
        color, label = colors[experiment], labels[experiment]
        ax_tend.plot(days, sel.member_tendency.mean("member"), color=color, lw=1.5, label=label)
        mean = sel.member_mean.mean("member")
        if "reference_mean" in sel:
            ax_bias.plot(days, mean - sel.reference_mean, color=color, lw=1.5)
            ax_err.plot(days, sel.ensemble_mean_rmse, color=color, lw=1.5)
            if not reference_drawn:
                ax_tend.plot(days, sel.reference_tendency, color="black", ls="--", lw=1.2, label="ERA5")
                reference_drawn = True
        else:
            ax_bias.plot(days, mean - mean.isel(lead=0), color=color, lw=1.5)
        ax_err.plot(days, sel.ensemble_spread, color=color, lw=1.5, ls=":")
        ax_range.plot(days, sel.member_max.max("member"), color=color, lw=1.2)
        ax_range.plot(days, sel.member_min.min("member"), color=color, lw=1.2, ls="--")
    first = next(iter(datasets.values()))
    units = first.member_mean.attrs.get("units", "")
    ax_tend.set_ylabel(f"|d{first.attrs['variable']}/dt| ({units} day-1)")
    ax_tend.set_title("Regional mean absolute tendency")
    has_reference = reference_drawn
    ax_bias.set_title("Ensemble-mean bias vs ERA5" if has_reference else "Ensemble-mean change since lead 0")
    ax_bias.axhline(0, color="gray", lw=0.7)
    ax_bias.set_ylabel(units)
    ax_err.set_title("Ensemble-mean RMSE (solid) and spread (dotted)" if has_reference
                     else "Ensemble spread")
    ax_err.set_ylabel(units)
    ax_range.set_title("Ensemble max of spatial max (solid), min of spatial min (dashed)")
    ax_range.set_ylabel(units)
    for ax in axes.flat:
        ax.set_xlim(0, zoom_days)
        ax.grid(alpha=0.3)
    for ax in axes[1]:
        ax.set_xlabel("Lead (days)")
    ax_tend.legend(fontsize=8)
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_shock_summary(summary, metrics, *, colors, title=""):
    """Grouped bars: one panel per metric, variables on the x axis, one bar per experiment."""
    metrics = [m for m in metrics if m in summary and summary[m].notna().any()]
    variables = list(dict.fromkeys(summary["variable"]))
    experiments = list(dict.fromkeys(summary["experiment"]))
    width = 0.8 / len(experiments)
    fig, axes = plt.subplots(1, len(metrics), figsize=(3.6 * len(metrics), 3.8), squeeze=False)
    x = np.arange(len(variables))
    for ax, metric in zip(axes[0], metrics):
        for index, experiment in enumerate(experiments):
            rows = summary[summary["experiment"] == experiment].set_index("variable")
            values = [rows[metric].get(v, np.nan) for v in variables]
            label = rows["label"].iloc[0] if len(rows) else experiment
            ax.bar(x + (index - (len(experiments) - 1) / 2) * width, values, width,
                   color=colors[experiment], label=label)
        if metric in {"tendency_ratio", "tendency_ratio_era5", "spread_ratio"}:
            ax.axhline(1, color="gray", lw=0.8, ls="--")
        else:
            ax.axhline(0, color="gray", lw=0.8)
        ax.set_xticks(x, variables)
        ax.set_title(metric, fontsize=10)
        ax.grid(axis="y", alpha=0.3)
    axes[0][0].legend(fontsize=8)
    fig.suptitle(title)
    fig.tight_layout()
    return fig
