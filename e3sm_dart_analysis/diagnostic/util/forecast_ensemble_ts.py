"""Hindcast ensemble time series and initial-shock diagnostics.

Each hindcast is one continuous forecast from the analysis time. Model records
carry ``time_bnds``: the first record has zero width and holds the initial
state; later records are interval means (``time: mean``) or end-of-interval
values (``time: point``). Every record is placed at its lead time (interval
midpoint for means, interval end for point values) and, where available,
verified against a reference averaged over the same interval: ERA5 6-hourly
samples, or daily means (GPCP or GPM IMERG precipitation, NOAA interpolated OLR).
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
from configs.experiment_labels import experiment_label

FREQUENCIES = ("6hourly", "daily", "daily_mean", "daily_00z", "h1")  # ts/<subdirectory>
REGIONS = ("global", "land", "ocean") + tuple(r for r in BOX_REGIONS if r != "global")
# Reference products on the model's 180x360 grid. "samples": 6-hourly values at
# 00/06/12/18Z; "daily": calendar-day means (any stamp within the day).
REFERENCE_PRODUCTS = {
    "ERA5": {"kind": "samples", "file": "ERA5/6hourly/ERA5.6hourly.en00.{variable}.{year}01-{year}12.nc"},
    "GPM": {"kind": "daily", "file": "IMERG/daily/{variable}.daily.{year}.nc"},
    "GPCP": {"kind": "daily", "file": "GPCP/daily/{variable}.daily.{year}.nc"},
    "NOAA OLR": {"kind": "daily", "file": "NOAA_IOLR/daily/{variable}.daily.{year}.nc"},
}
# Products available per variable; the first is the default.
REFERENCES = {
    "PSL": ("ERA5",),
    "TS": ("ERA5",),
    "TREFHT": ("ERA5",),
    "PRECT": ("ERA5", "GPM", "GPCP"),
    "FLUT": ("NOAA OLR",),
}
REFERENCE_VARIABLES = tuple(REFERENCES)
MEMBER_DISPLAYS = ("members", "quartiles", "none")
PLOT_KINDS = ("tendency", "timeseries")
MODE_DEFAULTS = {"references": {}, "member_display": "members", "zoom_days": 10, "plots": ("tendency",),
                 "core": (), "kind": "shock"}
MODE_OPTIONS = ("frequency", "variables", *MODE_DEFAULTS)
# Bump when cached content changes meaning; older caches are rebuilt.
CACHE_VERSION = 4
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
    "reference_tendency": ("reference regional mean of |dX/dt| (units per day); missing when the reference "
                           "is coarser than the model output; from raw sample differences when model means "
                           "are no longer than the sample spacing"),
    "ensemble_mean_rmse": "gridpoint RMSE of the ensemble mean against the reference",
}
SUMMARY_METRICS = {
    "tendency_ratio": "mean |dX/dt| over the shock window / over the settled window",
    "tendency_ratio_ref": "mean |dX/dt| over the shock window / the reference's over the same leads",
    "adjustment_days": "last lead before the settled window where the 1-day running |dX/dt| exceeds (1 + tolerance) x the settled level",
    "initial_bias": "ensemble-mean minus reference regional mean at the first lead with a reference (lead 0 for ERA5)",
    "drift": "bias change from that lead to the last day of the shock window",
    "spread_growth_ratio": "ensemble spread at the end of the shock window / at lead 0 (not spread/RMSE)",
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
    elif variable in _PRECIPITATION and units in {"mm/day", "mm d-1"}:
        units = "mm day-1"
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
            else:
                # Members may differ only by trailing records (e.g. one member
                # post-processed with the final 6-h interval kept); keep the
                # leads all members share.
                n = min(len(lead), len(reference_lead))
                if not np.array_equal(lead[:n], reference_lead[:n]):
                    raise ValueError(f"{member} lead times differ from EN01")
                if n < len(reference_lead):
                    print(f"[WARN] {member}: {len(reference_lead) - n} trailing record(s) missing; "
                          f"using the {n} leads every member shares")
                reference_lead, interval = reference_lead[:n], (interval[0][:n], interval[1][:n])
            field = _normalize(ds[variable].isel(time=keep)).astype("float32")
            field = field.drop_vars("time").rename(time="lead").assign_coords(lead=lead)
            fields.append(field.expand_dims(member=[member]))
    except Exception:
        for ds in opened:
            ds.close()
        raise
    fields = [f.isel(lead=slice(0, len(reference_lead))) for f in fields]
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


def _load_reference(reference_root, name, variable, start, end):
    spec = REFERENCE_PRODUCTS[name]
    files = [Path(reference_root) / spec["file"].format(variable=variable, year=year)
             for year in range(start.min().year, end.max().year + 1)]
    missing = [str(p) for p in files if not p.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing {name} files: " + ", ".join(missing))
    with xr.open_mfdataset(files, combine="by_coords") as ds:
        data = _normalize(ds[variable]).sel(time=slice(start.min().floor("D"), end.max())).load()
    return _units(data, variable)


def read_reference_sample_tendency(reference_root, variable, field, product=None):
    """|dX/dt| per day between the samples bounding each model interval, on the model grid.

    For model means no longer than the sample spacing, where ``read_reference``
    cannot give a matching tendency. The value at a lead spans that record's
    interval (centred on its midpoint); the model tendency at the same lead
    spans the previous and current midpoints, i.e. it is centred half an
    interval earlier. Zero-width records and products that are not sampled
    return None.
    """
    name = reference_product(variable, product)
    if REFERENCE_PRODUCTS[name]["kind"] != "samples":
        return None
    start = pd.DatetimeIndex(field.interval_start.values)
    end = pd.DatetimeIndex(field.interval_end.values)
    data = _load_reference(reference_root, name, variable, start, end)
    width = np.asarray((end - start) / pd.Timedelta(days=1), dtype=float)
    valid = width > 0
    at_start = data.sel(time=start[valid].values).values
    at_end = data.sel(time=end[valid].values).values
    values = np.full((len(start),) + data.shape[1:], np.nan, dtype="float32")
    values[valid] = np.abs(at_end - at_start) / width[valid][:, None, None]
    tendency = xr.DataArray(values, dims=("lead", "lat", "lon"),
                            coords={"lead": field.lead, "lat": field.lat, "lon": field.lon})
    tendency.attrs.update(units=data.attrs["units"], name=name)
    return tendency


def read_reference(reference_root, variable, field, product=None):
    """Reference averaged over each model interval, on the model grid.

    ``samples`` (ERA5 6-hourly): trapezoid mean of the samples in the interval;
    zero-width and instantaneous records take the sample at their time. For
    6-hourly model means this is the mean of the two bounding samples, an
    approximation of the interval mean. ``daily`` (GPM, GPCP, NOAA OLR): overlap-weighted
    mean of the calendar days the interval covers; zero-width records have no
    reference. Gridpoints the product does not cover stay missing.
    """
    name = reference_product(variable, product)
    spec = REFERENCE_PRODUCTS[name]
    start = pd.DatetimeIndex(field.interval_start.values)
    end = pd.DatetimeIndex(field.interval_end.values)
    data = _load_reference(reference_root, name, variable, start, end)
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
                           resolution_hours=_RESOLUTION_HOURS[spec["kind"]], kind=spec["kind"])
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


def derive_ensemble_ts(field, masks, reference=None, sample_tendency=None):
    """Regional member, spread, tendency, and (optional) reference diagnostics.

    ``sample_tendency`` (``read_reference_sample_tendency``) replaces the
    reference tendency where interval means of the samples cannot match the model.
    """
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
        # Sampled references cannot match model means no longer than the sample
        # spacing: the mean of two bounding samples turns consecutive differences
        # into halved 2-step differences, damping the diurnal cycle in |dX/dt|.
        resolution = reference.attrs["resolution_hours"]
        resolved = resolution <= step * 1.01
        if reference.attrs["kind"] == "samples" and not point:
            resolved &= step > resolution * 1.01
        instantaneous = not resolved and sample_tendency is not None and reference.attrs["kind"] == "samples"
        if instantaneous:
            tendency = _area_mean(sample_tendency, weights)
        out["reference_tendency"] = tendency if resolved or instantaneous else xr.full_like(tendency, np.nan)
        out["ensemble_mean_rmse"] = np.sqrt(_area_mean(error ** 2, weights))
    out = out.compute().assign_coords(intervals)
    for name in out.data_vars:
        units = field.attrs["units"]
        out[name].attrs["units"] = f"{units} day-1" if name.endswith("tendency") else units
        out[name].attrs["long_name"] = DIAGNOSTICS[name]
    if reference is not None and instantaneous:
        out["reference_tendency"].attrs["method"] = "instantaneous"
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
        point = "time: point" in str(field.attrs.get("cell_methods", ""))
        sample_tendency = (read_reference_sample_tendency(reference_root, variable, field, reference)
                           if reference and not point else None)
        ds = derive_ensemble_ts(field, masks, reference_field, sample_tendency)
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
            row = {"experiment": experiment, "label": experiment_label(experiment),
                   "variable": ds.attrs["variable"], "region": str(region),
                   "units": sel.member_mean.attrs.get("units", ""),
                   "tendency_ratio": float(tendency.where(shock).mean()) / settled_level}
            smooth = _running_daily_mean(tendency)
            above = smooth.where((days > 0) & (days < settled_days[0])) > (1 + tolerance) * settled_level
            exceed = days.values[above.values]
            row["adjustment_days"] = float(exceed.max()) if exceed.size else 0.0
            spread = sel.ensemble_spread
            row["spread_growth_ratio"] = float(spread.where(shock_end).mean()) / float(spread.isel(lead=0))
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
                             member_display="members", title="", fontz=16):
    """Ensemble mean and reference over the full forecast and the first ``zoom_days``.

    ``member_display``: "members" draws every member, "quartiles" shades the
    25th-75th percentile envelope, "none" shows only ensemble means. ``fontz`` is
    the base font size (titles 0.9x, axis labels 0.8x, ticks 0.75x, legend 0.7x).
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
        ax.set_xlabel("Lead (days)", fontsize=fontz * 0.8)
        ax.tick_params(labelsize=fontz * 0.75)
        ax.grid(alpha=0.3)
    first = next(iter(datasets.values()))
    axes[0].set_ylabel(_ylabel(first, "member_mean"), fontsize=fontz * 0.8)
    if zoomed:
        axes[1].set_xlim(0, zoom_days)
        axes[1].set_title(f"First {zoom_days:g} days", fontsize=fontz * 0.9)
    axes[0].legend(fontsize=fontz * 0.7, ncol=2)
    fig.suptitle(title, fontsize=fontz)
    fig.tight_layout()
    return fig


def _is_core(experiment, core):
    return not core or experiment.startswith(tuple(core))


_UNIT_LABELS = {"mm day-1": r"mm day$^{-1}$", "W m-2": r"W m$^{-2}$"}


def plot_shock_tendency(datasets_by_variable, region, *, labels, colors, core=(), settled_days=(10.0, 30.0),
                        zoom_days=10, title="", fontz=16):
    """Rows = variables; |dX/dt| with the reference (left) and |dX/dt| / settled level (right).

    Panels are labelled "(a) |dX/dt| (PSL)", "(b) Normalized |dX/dt| (PSL)", ...;
    ``core`` experiment-name prefixes are drawn solid, other experiments are
    dashed sensitivity curves. ``fontz`` is the base font size.
    """
    variables = list(datasets_by_variable)
    fig, axes = plt.subplots(len(variables), 2, figsize=(15, 3.6 * len(variables)), sharex=True, squeeze=False)
    handles, references = {}, set()
    for row, variable in enumerate(variables):
        ax_abs, ax_norm = axes[row]
        reference_drawn = False
        for experiment, ds in datasets_by_variable[variable].items():
            sel = ds.sel(region=region)
            days = sel.lead / 24.0
            tendency = sel.member_tendency.mean("member")
            settled = float(tendency.where((days >= settled_days[0]) & (days <= settled_days[1])).mean())
            style = dict(color=colors[experiment], lw=2.2, ls="-") if _is_core(experiment, core) \
                else dict(color=colors[experiment], lw=1.5, ls="--")
            line, = ax_abs.plot(days, tendency, **style)
            ax_norm.plot(days, tendency / settled, **style)
            handles[labels[experiment]] = line
            if not reference_drawn and "reference_tendency" in sel and bool(sel.reference_tendency.notnull().any()):
                name = ds.attrs.get("reference") or "reference"
                reference = sel.reference_tendency
                ls = ":"
                if reference.attrs.get("method") == "instantaneous":
                    name, ls = f"{name} (instantaneous)", "-."
                ref_line, = ax_abs.plot(days, reference, color="black", lw=1.5, ls=ls)
                ref_settled = float(reference.where((days >= settled_days[0]) & (days <= settled_days[1])).mean())
                ax_norm.plot(days, reference / ref_settled, color="black", lw=1.5, ls=ls)
                handles[name] = ref_line
                references.add(name)
                reference_drawn = True
        units = next(iter(datasets_by_variable[variable].values())).member_mean.attrs.get("units", "")
        ax_abs.set_title(f"({chr(97 + 2 * row)}) |dX/dt| ({variable})", loc="left", fontsize=fontz * 0.9, color="#2b2b29")
        ax_norm.set_title(f"({chr(98 + 2 * row)}) Normalized |dX/dt| ({variable})", loc="left",
                          fontsize=fontz * 0.9, color="#2b2b29")
        ax_abs.set_ylabel(f"{_UNIT_LABELS.get(units, units)} day$^{{-1}}$", fontsize=fontz * 0.8)
        ax_norm.axhline(1.0, color="#6f6e69", lw=1.0)
        for ax in (ax_abs, ax_norm):
            ax.set_xlim(0, zoom_days)
            ax.grid(alpha=0.3)
            ax.tick_params(labelsize=fontz * 0.75)
    for ax in axes[-1]:
        ax.set_xlabel("Lead (days)", fontsize=fontz * 0.8)
    order = [k for k in handles if k not in references] + [k for k in handles if k in references]
    fit_legend(fig, lambda size: fig.legend([handles[k] for k in order], order, loc="lower center",
                                            ncol=int(np.ceil(len(order) / 2)), frameon=False, fontsize=size,
                                            handlelength=2.5, columnspacing=1.8), fontz * 0.85)
    if title:
        fig.suptitle(title, fontsize=fontz)
    fig.tight_layout(rect=(0, 0.075, 1, 0.98 if title else 1))
    return fig


def plot_shock_summary(summary, metrics, *, colors, title="", fontz=16):
    """Grouped bars: one panel per metric, variables on the x axis, one bar per experiment."""
    metrics = [m for m in metrics if m in summary and summary[m].notna().any()]
    variables = list(dict.fromkeys(summary["variable"]))
    experiments = list(dict.fromkeys(summary["experiment"]))
    width = 0.8 / len(experiments)
    fig, axes = plt.subplots(1, len(metrics), figsize=(4.4 * len(metrics), 4.8), squeeze=False)
    x = np.arange(len(variables))
    handles = {}
    for index_m, (ax, metric) in enumerate(zip(axes[0], metrics)):
        for index, experiment in enumerate(experiments):
            rows = summary[summary["experiment"] == experiment].set_index("variable")
            values = [rows[metric].get(v, np.nan) for v in variables]
            label = rows["label"].iloc[0] if len(rows) else experiment
            bars = ax.bar(x + (index - (len(experiments) - 1) / 2) * width, values, width,
                          color=colors[experiment], label=label)
            handles[label] = bars
        if metric in {"tendency_ratio", "tendency_ratio_ref", "spread_growth_ratio"}:
            ax.axhline(1, color="gray", lw=0.8, ls="--")
        else:
            ax.axhline(0, color="gray", lw=0.8)
        ax.set_xticks(x, variables)
        ax.set_title(f"({chr(97 + index_m)}) {metric}", loc="left", fontsize=fontz * 0.9, color="#2b2b29")
        ax.tick_params(labelsize=fontz * 0.75)
        ax.grid(axis="y", alpha=0.3)
    fit_legend(fig, lambda size: fig.legend(list(handles.values()), list(handles), loc="lower center",
                                            ncol=int(np.ceil(len(handles) / 2)), frameon=False, fontsize=size),
               fontz * 0.85)
    if title:
        fig.suptitle(title, fontsize=fontz)
    fig.tight_layout(rect=(0, 0.14, 1, 0.95 if title else 1))
    return fig


def fit_legend(fig, make, fontsize, *, margin=0.98, steps=6):
    """Figure legend ``make(fontsize)``, its font shrunk until the box fits within the figure width.

    ``make`` draws the legend at a font size and returns it; a legend that fits
    is kept as drawn. Returns ``(legend, fontsize)``; spacing set in font units
    (handle length, padding) shrinks with the font.
    """
    legend = make(fontsize)
    renderer = fig.canvas.get_renderer()
    for _ in range(steps):
        width, limit = legend.get_window_extent(renderer).width, margin * fig.bbox.width
        if width <= limit:
            break
        legend.remove()
        fontsize *= limit / width
        legend = make(fontsize)
    return legend, fontsize


def results_colors(results):
    """Fixed color per experiment across every season and variable."""
    return experiment_colors(e for by_variable in results.values()
                             for datasets in by_variable.values() for e in datasets)


def plot_ensemble_ts_figures(results, registries, *, regions, figure_dir, tag, plots=("tendency",),
                             colors=None, core=(), settled_days=(10.0, 30.0), zoom_days=10,
                             member_display="members", dpi=200, show=True, log=print, fontz=16, title=True,
                             fmt="pdf"):
    """Draw and save the requested figure kinds as ``fmt`` files (pdf, png, ...); return the paths.

    ``tendency``: one figure per season/region with every variable;
    ``timeseries``: one figure per season/variable/region.
    """
    colors = colors or results_colors(results)
    figure_dir = Path(figure_dir)
    figure_dir.mkdir(parents=True, exist_ok=True)
    written = []

    def save(fig, name):
        output = figure_dir / name
        fig.savefig(output, dpi=dpi, bbox_inches="tight")
        written.append(output)
        log(f"[FIGURE] {output}")
        if not show:
            plt.close(fig)

    for season, by_variable in results.items():
        registry = registries[season]
        labels = {name: registry[name]["label"] for datasets in by_variable.values() for name in datasets}
        for region in regions:
            if "tendency" in plots:
                save(plot_shock_tendency(by_variable, region, labels=labels, colors=colors, core=core,
                                         settled_days=settled_days, zoom_days=zoom_days, fontz=fontz,
                                         title=f"Initial-shock tendency: {region} ({season} fc)" if title else ""),
                     f"fig_shock_tendency_{region}_{season}_fc_{tag}.{fmt}")
            if "timeseries" in plots:
                for variable, datasets in by_variable.items():
                    save(plot_ensemble_timeseries(datasets, region, labels=labels, colors=colors,
                                                  zoom_days=zoom_days, member_display=member_display,
                                                  fontz=fontz, title=f"{variable} {region} ({season} fc, {tag})"),
                         f"fig_{variable}_{region}_{season}_fc_{tag}_ensemble_ts.{fmt}")
    return written


def plot_shock_summary_figures(summary, *, figure_dir, tag, colors, metrics=None, dpi=200,
                               show=True, log=print, fontz=16, title=True, fmt="pdf"):
    """One grouped-bar summary figure per season and region, saved as ``fmt``; return the paths."""
    metrics = list(metrics or SUMMARY_METRICS)
    figure_dir = Path(figure_dir)
    figure_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for (season, region), table in summary.groupby(["season", "region"], sort=False):
        fig = plot_shock_summary(table, metrics, colors=colors, fontz=fontz,
                                 title=f"Initial-shock summary: {region} ({season} fc)" if title else "")
        output = figure_dir / f"fig_shock_summary_{region}_{season}_fc_{tag}.{fmt}"
        fig.savefig(output, dpi=dpi, bbox_inches="tight")
        written.append(output)
        log(f"[FIGURE] {output}")
        if not show:
            plt.close(fig)
    return written
