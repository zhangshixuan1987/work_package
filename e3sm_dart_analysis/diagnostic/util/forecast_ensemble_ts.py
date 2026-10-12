"""Hindcast ensemble time series and initial-shock diagnostics.

Each hindcast is one continuous forecast from the analysis time. Model records
carry ``time_bnds``: the first record has zero width and holds the initial
state; later records are interval means (``time: mean``) or end-of-interval
values (``time: point``). Every record is placed at its lead time (interval
midpoint for means, interval end for point values) and, where available,
verified against a reference averaged over the same interval: ERA5 6-hourly
samples, or daily means (GPM IMERG precipitation, NOAA interpolated OLR).
"""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from util.bias_metrics import _period_years, experiment_registry
from util.forecast_regional_timeseries import REGIONS as BOX_REGIONS
from util.model_evaluation import experiment_colors

FREQUENCIES = ("6hourly", "daily")
REGIONS = ("global", "land", "ocean") + tuple(r for r in BOX_REGIONS if r != "global")
# Reference products on the model's 180x360 grid. "samples": 6-hourly values at
# 00/06/12/18Z; "daily": calendar-day means (any stamp within the day).
REFERENCE_PRODUCTS = {
    "ERA5": {"kind": "samples", "file": "ERA5/6hourly/ERA5.6hourly.en00.{variable}.{year}01-{year}12.nc"},
    "GPM": {"kind": "daily", "file": "IMERG/daily/{variable}.daily.{year}.nc"},
    "NOAA OLR": {"kind": "daily", "file": "NOAA_IOLR/daily/{variable}.daily.{year}.nc"},
}
# Products available per variable; the first is the default.
REFERENCES = {
    "PSL": ("ERA5",),
    "TS": ("ERA5",),
    "TREFHT": ("ERA5",),
    "PRECT": ("ERA5", "GPM"),
    "FLUT": ("NOAA OLR",),
}
REFERENCE_VARIABLES = tuple(REFERENCES)
MEMBER_DISPLAYS = ("members", "quartiles", "none")
PLOT_KINDS = ("timeseries", "shock")
_FIGURE_SUFFIX = {"timeseries": "ensemble_ts", "shock": "shock_diagnostics"}
MODE_DEFAULTS = {"references": {}, "member_display": "members", "zoom_days": 10, "plots": PLOT_KINDS}
MODE_OPTIONS = ("frequency", "variables", *MODE_DEFAULTS)
# Bump when cached content changes meaning; older caches are rebuilt.
CACHE_VERSION = 2
_RESOLUTION_HOURS = {"samples": 6.0, "daily": 24.0}
_HPA = {"PS", "PSL"}
_PRECIPITATION = {"PRECT", "PRECC", "PRECL"}

DIAGNOSTICS = {
    "member_mean": "regional mean per member",
    "member_min": "regional minimum per member",
    "member_max": "regional maximum per member",
    "member_tendency": "regional mean of |dX/dt| per member (units per day)",
    "ensemble_spread": "regional mean of the gridpoint ensemble standard deviation",
    "reference_mean": "reference regional mean over the same intervals",
    "reference_tendency": "reference regional mean of |dX/dt| (units per day); missing when the reference is coarser than the model output",
    "ensemble_mean_rmse": "gridpoint RMSE of the ensemble mean against the reference",
}
SUMMARY_METRICS = {
    "tendency_ratio": "mean |dX/dt| over the shock window / over the settled window",
    "tendency_ratio_ref": "mean |dX/dt| over the shock window / the reference's over the same leads",
    "adjustment_days": "last lead before the settled window where the 1-day running |dX/dt| exceeds (1 + tolerance) x the settled level",
    "initial_bias": "ensemble-mean minus reference regional mean at the first lead with a reference (lead 0 for ERA5)",
    "drift": "bias change from that lead to the last day of the shock window",
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


def reference_product(variable, product=None):
    """Validated reference product name for ``variable`` (None -> its default)."""
    choices = REFERENCES.get(variable, ())
    if product is None:
        return choices[0] if choices else None
    if product not in choices:
        raise ValueError(f"{variable}: reference {product!r} unavailable; choose from {list(choices)}")
    return product


def read_reference(reference_root, variable, field, product=None):
    """Reference averaged over each model interval, on the model grid.

    ``samples`` (ERA5 6-hourly): trapezoid mean of the samples in the interval;
    zero-width and instantaneous records take the sample at their time. For
    6-hourly model means this is the mean of the two bounding samples, an
    approximation of the interval mean. ``daily`` (GPM, NOAA OLR): overlap-weighted
    mean of the calendar days the interval covers; zero-width records have no
    reference. Gridpoints the product does not cover stay missing.
    """
    name = reference_product(variable, product)
    spec = REFERENCE_PRODUCTS[name]
    start = pd.DatetimeIndex(field.interval_start.values)
    end = pd.DatetimeIndex(field.interval_end.values)
    files = [Path(reference_root) / spec["file"].format(variable=variable, year=year)
             for year in range(start.min().year, end.max().year + 1)]
    missing = [str(p) for p in files if not p.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing {name} files: " + ", ".join(missing))
    with xr.open_mfdataset(files, combine="by_coords") as ds:
        data = _normalize(ds[variable]).sel(time=slice(start.min().floor("D"), end.max())).load()
    data = _units(data, variable)
    times = pd.DatetimeIndex(data.time.values)
    day_start = times.floor("D").values
    day_end = day_start + np.timedelta64(1, "D")
    weights = np.zeros((len(start), len(times)))
    point = "time: point" in str(field.attrs.get("cell_methods", ""))
    for row, (a, b) in enumerate(zip(start, end)):
        if spec["kind"] == "daily":
            if a < b:
                a64, b64 = a.to_datetime64(), b.to_datetime64()
                overlap = (np.minimum(day_end, b64) - np.maximum(day_start, a64)) / np.timedelta64(1, "h")
                weights[row] = np.clip(overlap, 0, None)
        elif a == b or point:
            weights[row, times == b] = 1.0
        else:
            inside = (times >= a) & (times <= b)
            weights[row, inside] = 1.0
            weights[row, times == a] *= 0.5
            weights[row, times == b] *= 0.5
        if weights[row].sum() == 0:
            weights[row] = np.nan
    weights /= weights.sum(axis=1, keepdims=True)
    weights = xr.DataArray(weights, dims=("lead", "time"), coords={"lead": field.lead})
    reference = xr.dot(weights, data.fillna(0), dim="time").where(data.notnull().all("time"))
    reference = reference.assign_coords(lat=field.lat, lon=field.lon)
    reference.attrs.update(units=data.attrs["units"], name=name,
                           resolution_hours=_RESOLUTION_HOURS[spec["kind"]])
    return reference


def _area_mean(values, weights):
    return values.weighted(weights).mean(("lat", "lon"), skipna=True)


def _abs_tendency(values, point=True):
    """|dX/dt| per day between consecutive records, labelled by the later record.

    For interval means (``point=False``) the step from the zero-width initial
    snapshot to the first mean compares different quantities, so tendencies
    start at the first interval mean.
    """
    series = values if point else values.isel(lead=slice(1, None))
    step = np.diff(series.lead.values) / 24.0
    change = np.abs(series.diff("lead")) / xr.DataArray(step, dims="lead", coords={"lead": series.lead.values[1:]})
    return change.reindex(lead=values.lead)


def derive_ensemble_ts(field, masks, reference=None):
    """Regional member, spread, tendency, and (optional) reference diagnostics."""
    intervals = {name: field[name] for name in ("interval_start", "interval_end")}
    point = "time: point" in str(field.attrs.get("cell_methods", ""))
    field = field.drop_vars(list(intervals))
    if reference is not None:
        reference = reference.drop_vars(list(intervals), errors="ignore")
        # Restrict every diagnostic to gridpoints the reference covers, so model
        # and reference regional means describe the same area.
        masks = masks & reference.notnull().any("lead")
    weights = np.cos(np.deg2rad(field.lat)).clip(min=0) * masks
    masked = field.where(masks)
    out = xr.Dataset({
        "member_mean": _area_mean(field, weights),
        "member_min": masked.min(("lat", "lon")),
        "member_max": masked.max(("lat", "lon")),
        "member_tendency": _area_mean(_abs_tendency(field, point), weights),
        "ensemble_spread": _area_mean(field.std("member", ddof=1), weights),
    })
    if reference is not None:
        error = field.mean("member") - reference
        out["reference_mean"] = _area_mean(reference, weights)
        step = float(np.median(np.diff(field.lead.values[1:])))
        tendency = _area_mean(_abs_tendency(reference, point), weights)
        # A reference coarser than the model output has no within-interval change.
        resolved = reference.attrs["resolution_hours"] <= step * 1.01
        out["reference_tendency"] = tendency if resolved else xr.full_like(tendency, np.nan)
        out["ensemble_mean_rmse"] = np.sqrt(_area_mean(error ** 2, weights))
    out = out.compute().assign_coords(intervals)
    for name in out.data_vars:
        units = field.attrs["units"]
        out[name].attrs["units"] = f"{units} day-1" if name.endswith("tendency") else units
        out[name].attrs["long_name"] = DIAGNOSTICS[name]
    return out


def ensemble_ts_cache_path(cache_dir, experiment, variable, frequency, period):
    return Path(cache_dir) / f"{experiment}_{variable}_{frequency}_{period}_ensemble_ts.nc"


def ensemble_ts_cache_status(path, *, regions, members, reference=None):
    """'' when the cache covers the request, otherwise why it must be rebuilt."""
    path = Path(path)
    if not path.is_file():
        return "missing"
    with xr.open_dataset(path) as cached:
        if int(cached.attrs.get("cache_version", 1)) != CACHE_VERSION:
            return f"cache version {cached.attrs.get('cache_version', 1)}, current {CACHE_VERSION}"
        missing = set(regions) - set(cached.region.values.tolist())
        if missing:
            return f"regions {sorted(missing)} not cached"
        if cached.sizes["member"] != members:
            return f"members: cached {cached.sizes['member']}, requested {members}"
        cached_reference = str(cached.attrs.get("reference", ""))
        if reference and not cached_reference.startswith(reference):
            return f"reference: cached {cached_reference or 'none'}, requested {reference}"
    return ""


def cache_ensemble_ts(path, *, config, experiment, variable, period, frequency, regions,
                      landmask_file, reference_root=None, reference=None, max_members=None,
                      force_compute=False, log=print):
    """Compute or reuse one experiment/variable cache covering ``regions``.

    ``reference`` picks the product for variables with several (see
    ``REFERENCES``); None uses the default. No ``reference_root`` -> model only.
    """
    path = Path(path)
    members = config["nens"] if max_members is None else min(max_members, config["nens"])
    reference = reference_product(variable, reference) if reference_root is not None else None
    reason = "forced" if force_compute else ensemble_ts_cache_status(
        path, regions=regions, members=members, reference=reference)
    if not reason:
        log(f"[REUSE] {path.name}")
        return xr.load_dataset(path).sel(region=list(regions))
    field, opened = read_member_ensemble(config, variable, period, frequency=frequency,
                                         max_members=max_members)
    try:
        masks = region_masks(field.lat, field.lon, regions, landmask_file)
        reference_field = (read_reference(reference_root, variable, field, reference)
                           if reference else None)
        ds = derive_ensemble_ts(field, masks, reference_field)
    finally:
        for opened_ds in opened:
            opened_ds.close()
    ds.attrs.update(experiment=experiment, label=config["label"], variable=variable,
                    period=period, frequency=frequency, ensemble_members=members,
                    reference=reference or "", cache_version=CACHE_VERSION,
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
                row["tendency_ratio_ref"] = (float(tendency.where(shock).mean())
                                             / float(sel.reference_tendency.where(shock).mean()))
                bias = (sel.member_mean.mean("member") - sel.reference_mean).dropna("lead")
                row["initial_bias"] = float(bias.isel(lead=0))
                row["drift"] = float(bias.where(shock_end).mean()) - row["initial_bias"]
            rows.append(row)
    columns = ["experiment", "label", "variable", "region", "units", *SUMMARY_METRICS]
    return pd.DataFrame(rows).reindex(columns=columns)


def validate_mode(mode, modes):
    """Settings of ``modes[mode]`` with defaults filled in, after checking every option."""
    if mode not in modes:
        raise ValueError(f"Unknown mode {mode!r}; choose from {list(modes)}")
    config = {**MODE_DEFAULTS, **modes[mode]}
    missing = {"frequency", "variables"} - set(config)
    if missing:
        raise ValueError(f"Mode {mode!r} lacks {sorted(missing)}")
    unknown = set(config) - set(MODE_OPTIONS)
    if unknown:
        raise ValueError(f"{mode}: unknown options {sorted(unknown)}; choose from {list(MODE_OPTIONS)}")
    if config["frequency"] not in FREQUENCIES:
        raise ValueError(f"{mode}: frequency must be one of {FREQUENCIES}")
    if config["member_display"] not in MEMBER_DISPLAYS:
        raise ValueError(f"{mode}: member_display must be one of {MEMBER_DISPLAYS}")
    unknown = set(config["plots"]) - set(PLOT_KINDS)
    if unknown:
        raise ValueError(f"{mode}: unknown plots {sorted(unknown)}; choose from {PLOT_KINDS}")
    for variable, product in config["references"].items():
        reference_product(variable, product)
    return config


def validate_regions(regions):
    unknown = set(regions) - set(REGIONS)
    if unknown:
        raise ValueError(f"Unknown regions {sorted(unknown)}; choose from {list(REGIONS)}")


def read_or_process_ensemble_ts(cases, variables, *, data_root, cache_dir, frequency, regions,
                                landmask_file, reference_root=None, references=None,
                                max_members=None, force_compute=False, log=print):
    """Reuse or build every case/variable/experiment cache.

    Returns ``(results, registries, failures)`` with
    ``results[season][variable][experiment] -> Dataset``; a failing request is
    logged and skipped.
    """
    references = references or {}
    registries, results, failures = {}, {}, []
    for season, experiments in cases.items():
        registry = registries[season] = experiment_registry(data_root, season, "fc")
        unknown = set(experiments) - set(registry)
        if unknown:
            raise ValueError(f"{season}: {sorted(unknown)} unavailable; choose from {list(registry)}")
        for variable in variables:
            for experiment in experiments:
                meta = registry[experiment]
                target = ensemble_ts_cache_path(cache_dir, experiment, variable, frequency, meta["period"])
                try:
                    ds = cache_ensemble_ts(
                        target, config=meta, experiment=experiment, variable=variable,
                        period=meta["period"], frequency=frequency, regions=regions,
                        landmask_file=landmask_file, reference_root=reference_root,
                        reference=references.get(variable), max_members=max_members,
                        force_compute=force_compute, log=log,
                    )
                except Exception as error:  # keep going; report at the end
                    failures.append((season, variable, experiment, f"{type(error).__name__}: {error}"))
                    log(f"[FAIL] {season} {variable} {experiment}: {type(error).__name__}: {error}")
                    continue
                results.setdefault(season, {}).setdefault(variable, {})[experiment] = ds
    return results, registries, failures


def shock_summary_table(results, *, shock_days=2.0, settled_days=(10.0, 30.0), tolerance=0.2):
    """``shock_summary`` for every season and variable, with a leading season column."""
    tables = [
        shock_summary(datasets, shock_days=shock_days, settled_days=settled_days,
                      tolerance=tolerance).assign(season=season)
        for season, by_variable in results.items()
        for datasets in by_variable.values()
    ]
    if not tables:
        return pd.DataFrame()
    table = pd.concat(tables, ignore_index=True)
    return table[["season"] + [c for c in table.columns if c != "season"]]


def _ylabel(ds, name):
    return f"{ds.attrs['variable']} ({ds[name].attrs.get('units', '')})"


def plot_ensemble_timeseries(datasets, region, *, labels, colors, zoom_days=10,
                             member_display="members", title=""):
    """Ensemble mean and reference over the full forecast and the first ``zoom_days``.

    ``member_display``: "members" draws every member, "quartiles" shades the
    25th-75th percentile envelope, "none" shows only ensemble means.
    """
    if member_display not in MEMBER_DISPLAYS:
        raise ValueError(f"member_display must be one of {MEMBER_DISPLAYS}")
    last_day = max(float(ds.lead.max()) for ds in datasets.values()) / 24.0
    zoomed = zoom_days is not None and zoom_days < last_day
    if zoomed:
        fig, axes = plt.subplots(1, 2, figsize=(14, 4.5), gridspec_kw={"width_ratios": (3, 2)})
    else:  # the zoom would repeat the full forecast
        fig, axis = plt.subplots(figsize=(10, 4.5))
        axes = [axis]
    reference = None
    for experiment, ds in datasets.items():
        sel = ds.sel(region=region)
        days = sel.lead / 24.0
        color = colors[experiment]
        if member_display == "quartiles":
            low, high = sel.member_mean.quantile([0.25, 0.75], "member")
        for ax in axes:
            if member_display == "members":
                for member in sel.member.values:
                    ax.plot(days, sel.member_mean.sel(member=member), color=color, lw=0.5, alpha=0.3)
            elif member_display == "quartiles":
                ax.fill_between(days, low, high, color=color, alpha=0.2, lw=0)
            ax.plot(days, sel.member_mean.mean("member"), color=color, lw=2,
                    label=labels[experiment])
        if reference is None and "reference_mean" in sel:
            reference, reference_name = sel.reference_mean, ds.attrs.get("reference") or "reference"
    for ax in axes:
        if reference is not None:
            ax.plot(reference.lead / 24.0, reference, color="black", lw=1.5, ls="--", label=reference_name)
        ax.set_xlabel("Lead (days)")
        ax.grid(alpha=0.3)
    first = next(iter(datasets.values()))
    axes[0].set_ylabel(_ylabel(first, "member_mean"))
    if zoomed:
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
    reference_name = None
    for experiment, ds in datasets.items():
        sel = ds.sel(region=region)
        days = sel.lead / 24.0
        color, label = colors[experiment], labels[experiment]
        ax_tend.plot(days, sel.member_tendency.mean("member"), color=color, lw=1.5, label=label)
        mean = sel.member_mean.mean("member")
        if "reference_mean" in sel:
            ax_bias.plot(days, mean - sel.reference_mean, color=color, lw=1.5)
            ax_err.plot(days, sel.ensemble_mean_rmse, color=color, lw=1.5)
            if reference_name is None:
                reference_name = ds.attrs.get("reference") or "reference"
                if bool(sel.reference_tendency.notnull().any()):
                    ax_tend.plot(days, sel.reference_tendency, color="black", ls="--", lw=1.2,
                                 label=reference_name)
        else:
            ax_bias.plot(days, mean - mean.isel(lead=0), color=color, lw=1.5)
        ax_err.plot(days, sel.ensemble_spread, color=color, lw=1.5, ls=":")
        ax_range.plot(days, sel.member_max.max("member"), color=color, lw=1.2)
        ax_range.plot(days, sel.member_min.min("member"), color=color, lw=1.2, ls="--")
    first = next(iter(datasets.values()))
    units = first.member_mean.attrs.get("units", "")
    ax_tend.set_ylabel(f"|d{first.attrs['variable']}/dt| ({units} day-1)")
    ax_tend.set_title("Regional mean absolute tendency")
    has_reference = reference_name is not None
    ax_bias.set_title(f"Ensemble-mean bias vs {reference_name}" if has_reference
                      else "Ensemble-mean change since lead 0")
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
        if metric in {"tendency_ratio", "tendency_ratio_ref", "spread_ratio"}:
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


def results_colors(results):
    """Fixed color per experiment across every season and variable."""
    return experiment_colors(e for by_variable in results.values()
                             for datasets in by_variable.values() for e in datasets)


def plot_ensemble_ts_figures(results, registries, *, regions, figure_dir, tag, plots=PLOT_KINDS,
                             colors=None, zoom_days=10, member_display="members", dpi=200,
                             show=True, log=print):
    """Draw and save the requested figure kinds per season/variable/region; return the paths."""
    colors = colors or results_colors(results)
    figure_dir = Path(figure_dir)
    figure_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for season, by_variable in results.items():
        registry = registries[season]
        for variable, datasets in by_variable.items():
            labels = {name: registry[name]["label"] for name in datasets}
            for region in regions:
                for kind in plots:
                    if kind == "timeseries":
                        fig = plot_ensemble_timeseries(
                            datasets, region, labels=labels, colors=colors, zoom_days=zoom_days,
                            member_display=member_display,
                            title=f"{variable} {region} ({season} fc, {tag})")
                    else:
                        fig = plot_shock_diagnostics(
                            datasets, region, labels=labels, colors=colors, zoom_days=zoom_days,
                            title=f"{variable} {region}: initial shock ({season} fc, {tag})")
                    output = figure_dir / f"{variable}_{region}_{season}_fc_{tag}_{_FIGURE_SUFFIX[kind]}.png"
                    fig.savefig(output, dpi=dpi, bbox_inches="tight")
                    written.append(output)
                    log(f"[FIGURE] {output}")
                    if not show:
                        plt.close(fig)
    return written


def plot_shock_summary_figures(summary, *, figure_dir, tag, colors, metrics=None, dpi=200,
                               show=True, log=print):
    """One grouped-bar summary figure per season and region; return the paths."""
    metrics = list(metrics or SUMMARY_METRICS)
    figure_dir = Path(figure_dir)
    figure_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for (season, region), table in summary.groupby(["season", "region"], sort=False):
        fig = plot_shock_summary(table, metrics, colors=colors,
                                 title=f"Initial-shock summary: {region} ({season} fc, {tag})")
        output = figure_dir / f"shock_summary_{region}_{season}_fc_{tag}.png"
        fig.savefig(output, dpi=dpi, bbox_inches="tight")
        written.append(output)
        log(f"[FIGURE] {output}")
        if not show:
            plt.close(fig)
    return written
