"""Biweekly-window scorecard of hindcast errors: RMSE and CRPS skill, anomaly pattern correlation.

For each season, variable, and experiment the daily ensemble mean is averaged
over the lead windows (week1-2 ... week7-8, day 1 = initial date) and compared
with the daily reference averaged over the same days. Days where the reference
is missing at a gridpoint (IMERG, ESA CCI gaps) are left out of both window
means there. Scores over a region (cos-lat weighted):

* ``rmse``: RMSE of the window-mean ensemble mean against the reference;
* ``rmse_clim``: RMSE of the climatology forecast (the reference climatology of
  those days) against the reference;
* ``msss``: mean-square skill score ``1 - (rmse / rmse_clim)**2``: 1 perfect,
  0 no better than climatology, negative worse than climatology;
* ``bias``: area-mean model minus reference;
* ``crps``: area-mean fair CRPS of the ensemble members' window means, which
  corrects for ensemble size so 10- and 40-member ensembles compare fairly;
* ``crps_clim``: the same for the climatological ensemble, the window means of
  each climatology year;
* ``crpss``: CRPS skill score ``1 - crps / crps_clim``;
* ``pcc``: centered pattern correlation of window-mean anomalies, model and
  reference both relative to the same reference climatology of those days.

With ``remove_mean_bias`` the region's area-mean bias is removed from both the
model and the climatology forecast before their RMSE (centered RMSE), so a
uniform offset such as a global surface-pressure loss does not dominate. CRPS
is always scored on the full fields.

Daily values match in time between model and reference:

* daily means (``sampling="daily_mean"``): model ``ts/daily_mean`` 00Z-24Z means
  (ELM ``ts/daily`` for land) against daily-mean references;
* PSL, TREFHT, TS, U10 and the surface fluxes (``"6hourly_mean"``): also
  ``ts/daily_mean``, which builds them from the 6-hourly output (6-h means
  weighted over the day; instantaneous PSL as the mean of the 00, 06, 12, and
  18Z states, the same construction as the ERA5 daily means they are verified
  against). Runs without ``daily_mean`` get the same values from ``ts/6hourly``
  via ``util.model_daily``; the label is kept so existing caches stay valid;
* U, V, T, Q, Z on pressure levels (``"00z"``): only 00 UTC states exist
  (``ts/daily_00z``), so day d is the state at the end of the day (00 UTC on
  d + 1; the lead-0 state is never used) and ERA5 is sampled at the same
  instants, for both the reference and the climatology.

The climatology comes from one shared period (``clim_years``) for every
variable: daily references are averaged per calendar day and smoothed with a
31-day running mean; references kept only as monthly means (ERA5 TREFHT, TS,
LHFLX, SHFLX) are averaged per month and interpolated to days. Land variables
are scored over land only. Pressure-level variables are scored only where the
level is above ground: the ERA5 surface-pressure climatology of the season's
months must exceed the level (850 hPa is below ground over Antarctica, Tibet,
the Andes, and Greenland, where model and ERA5 fields are extrapolations).
"""
from __future__ import annotations

import multiprocessing
import os
import re
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.collections import PolyCollection
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap

from configs.s2s_config import LEAD_WINDOWS
from util.bias_metrics import experiment_registry
from util.forecast_ensemble_ts import region_masks
from util.lnd_sm_daily_skill import top_layer_weights

WINDOWS = LEAD_WINDOWS  # {"week1-2": (1, 14), ...}; day 1 = initial date
CACHE_VERSION = 3  # 2: per-member and per-climatology-year window means; 3: time-matched sampling
# Recorded in the 6-hourly variables' model caches; changing it rebuilds only those.
SIX_HOURLY_DAILY = "00Z-00Z calendar days; instantaneous: mean of 00/06/12/18Z (as ERA5 daily)"
_G = 9.80665
_ERA5_00Z = "ERA5/6hourly/ERA5.6hourly.en00.{name}.{year}01-{year}12.nc"
_ERA5_DAILY = "ERA5/daily/{name}_{year}01_{year}12.nc"
_ERA5_MONTHLY = "ERA5/monthly/ERA5_analysis_monthly_{year}.nc"
_IMERG = "IMERG/daily/{name}.daily.{year}.nc"
_NOAA_OLR = "NOAA_IOLR/daily/{name}.daily.{year}.nc"
_ESA_CCI = "ESA_CCI/H20SOI/daily/{name}.daily.{year}.nc"

# component: atm or lnd (lnd is scored over land). model: (realm, variable).
# reference: (label, file template, variable). climatology: "daily" from the
# reference itself, or ("monthly", ERA5 monthly variable). sampling (default
# "daily_mean"): see the module docstring.
VARIABLES = {
    "PSL": dict(component="atm", model=("atm", "PSL"), reference=("ERA5", _ERA5_DAILY, "PSL"), climatology="daily", units="hPa", sampling="6hourly_mean"),
    **{f"{name}{level}": dict(component="atm", model=("atm", f"{name}{level}"),
                              reference=("ERA5 00Z", _ERA5_00Z, f"{name}{level}"), climatology="daily", units=units,
                              sampling="00z")
       for name, units in (("U", "m s-1"), ("V", "m s-1"), ("T", "K"), ("Q", "g kg-1"), ("Z", "m"))
       for level in (200, 500, 700, 850, 925)},
    "TREFHT": dict(component="atm", model=("atm", "TREFHT"), reference=("ERA5", _ERA5_DAILY, "TREFHT"), climatology=("monthly", "TREFHT"), units="K", sampling="6hourly_mean"),
    "PRECT": dict(component="atm", model=("atm", "PRECT"), reference=("GPM IMERG", _IMERG, "PRECT"), climatology="daily", units="mm day-1"),
    "FLUT": dict(component="atm", model=("atm", "FLUT"), reference=("NOAA OLR", _NOAA_OLR, "FLUT"), climatology="daily", units="W m-2"),
    # Surface: 00Z-00Z daily means of the 6-hourly output vs ERA5 daily means. EAM U10 is
    # the 10-m wind speed; ERA5 WSPD10 is the daily mean of the 6-hourly sqrt(u10^2 + v10^2).
    # ERA5 daily TS and fluxes exist only for 2011-2012, so their climatology is monthly.
    # ERA5 daily FLNS is corrupt (values 0-21 W m-2), so net long-wave is not scored.
    "TS": dict(component="atm", model=("atm", "TS"), reference=("ERA5", _ERA5_DAILY, "TS"), climatology=("monthly", "TS"), units="K", sampling="6hourly_mean"),
    "WSPD10": dict(component="atm", model=("atm", "U10"), reference=("ERA5", _ERA5_DAILY, "WSPD10"), climatology="daily", units="m s-1", sampling="6hourly_mean"),
    **{name: dict(component="atm", model=("atm", name), reference=("ERA5", _ERA5_DAILY, name), climatology=("monthly", name),
                  units="W m-2", sampling="6hourly_mean")
       for name in ("LHFLX", "SHFLX", "FSNS", "FSDS", "FLDS")},
    "TSA": dict(component="lnd", model=("lnd", "TSA"), reference=("ERA5", _ERA5_DAILY, "TREFHT"), climatology=("monthly", "TREFHT"), units="K"),
    "TG": dict(component="lnd", model=("lnd", "TG"), reference=("ERA5", _ERA5_DAILY, "TS"), climatology=("monthly", "TS"), units="K"),
    "SOILM_5CM": dict(component="lnd", model=("lnd", "H2OSOI"), reference=("ESA CCI", _ESA_CCI, "H2OSOI"), climatology="daily", units="m3 m-3"),
    "EFLX_LH_TOT": dict(component="lnd", model=("lnd", "EFLX_LH_TOT"), reference=("ERA5", _ERA5_DAILY, "LHFLX"), climatology=("monthly", "LHFLX"), units="W m-2"),
    "FSH": dict(component="lnd", model=("lnd", "FSH"), reference=("ERA5", _ERA5_DAILY, "SHFLX"), climatology=("monthly", "SHFLX"), units="W m-2"),
}
# Years each reference source covers; the shared climatology period must fit all used.
SOURCE_YEARS = {"ERA5 daily": (1979, 2017), "ERA5 00Z": (1979, 2017), "ERA5 monthly": (2001, 2019),
                "GPM IMERG": (2001, 2020), "NOAA OLR": (2001, 2020), "ESA CCI": (2001, 2020)}


def _climatology_source(variable):
    spec = VARIABLES[variable]
    if spec["climatology"] != "daily":
        return "ERA5 monthly"
    label = spec["reference"][0]
    return "ERA5 daily" if label == "ERA5" else label


def climatology_years(variables, clim_years, exclude_years=()):
    """Years used for the shared climatology, after checking every source covers them."""
    y0, y1 = clim_years
    for variable in variables:
        source = _climatology_source(variable)
        s0, s1 = SOURCE_YEARS[source]
        if y0 < s0 or y1 > s1:
            raise ValueError(f"{variable}: {source} covers {s0}-{s1}, outside clim_years {clim_years}")
    years = [year for year in range(y0, y1 + 1) if year not in set(exclude_years)]
    if not years:
        raise ValueError("No climatology years left after exclusions")
    return years


def window_days(init_date):
    """{window: DatetimeIndex of its days}; day 1 is ``init_date``."""
    start = pd.Timestamp(init_date)
    return {label: pd.date_range(start + pd.Timedelta(days=d0 - 1), start + pd.Timedelta(days=d1 - 1))
            for label, (d0, d1) in WINDOWS.items()}


_PRESSURE_LEVEL = re.compile(r"[UVTQZ]\d{3}")


def _convert(field, variable):
    """Model and reference fields in the scorecard units."""
    units = str(field.attrs.get("units", ""))
    if variable == "PSL" and units == "Pa":
        field = field * 0.01
    elif _PRESSURE_LEVEL.fullmatch(variable) and variable[0] == "Z" and units.startswith("m**2"):
        field = field / _G
    elif _PRESSURE_LEVEL.fullmatch(variable) and variable[0] == "Q":
        field = field * 1000.0
    elif variable == "PRECT" and units in {"m/s", "m s-1"}:
        field = field * 86400 * 1000
    return field.astype("float32")


def _normalize(field):
    rename = {k: v for k, v in (("latitude", "lat"), ("longitude", "lon")) if k in field.dims}
    field = field.rename(rename) if rename else field
    extra = [d for d in field.dims if d not in ("time", "lat", "lon") and field.sizes[d] == 1]
    if extra:
        field = field.squeeze(extra, drop=True)
    return field.sortby("lat")


def _daily(field):
    return field.assign_coords(time=field.time.dt.floor("D"))


def _sampling(variable):
    return VARIABLES[variable].get("sampling", "daily_mean")


def _on_days(field, days, variable, *, reindex=False):
    """``field`` on ``days`` (time labelled by day): the daily value, or for "00z"
    variables the snapshot at 00 UTC ending each day."""
    if _sampling(variable) == "00z":
        instants = days + pd.Timedelta(days=1)
        field = field.reindex(time=instants) if reindex else field.sel(time=instants)
        return field.assign_coords(time=days)
    field = _daily(field)
    return field.reindex(time=days) if reindex else field.sel(time=days)


def _file_years(days, variable):
    shift = pd.Timedelta(days=1) if _sampling(variable) == "00z" else pd.Timedelta(0)
    return sorted({(day + shift).year for day in days})


def _top_soil(field, land_support_dir, depth_cm=5.0):
    with xr.open_dataset(Path(land_support_dir) / "dzsoi_elm.nc") as ds:
        weights = top_layer_weights(ds["DZSOI"].load(), depth_cm)
    weights = weights.isel(levgrnd=slice(0, int(weights.notnull().any(("lat", "lon")).sum())))
    field = field.isel(levgrnd=slice(0, weights.sizes["levgrnd"]))
    return (field * weights.assign_coords(levgrnd=field.levgrnd)).sum("levgrnd", min_count=1)


def read_reference_days(reference_root, variable, days):
    """Daily reference on ``days`` (time, lat, lon), in scorecard units."""
    label, template, name = VARIABLES[variable]["reference"]
    files = [Path(reference_root) / template.format(name=name, year=year) for year in _file_years(days, variable)]
    missing = [str(p) for p in files if not p.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing {label} files: " + ", ".join(missing))
    with xr.open_mfdataset(files, combine="by_coords") as ds:
        field = _on_days(_normalize(ds[name]), days, variable).load()
    return _convert(field, variable)


def read_model_days(config, variable, days, *, land_support_dir, max_members=None, log=print):
    """Daily member fields on ``days`` (member, time, lat, lon), each record labelled by the day it covers.

    Files are found by ``util.model_daily`` (``ts/daily_mean``, else 6-hourly
    output averaged to 00Z-00Z days, else ``daily_00z`` 00Z states; the older
    ``ts/daily`` layout and the ELM daily series are read the same way).
    """
    from util.model_daily import daily_source, open_model_daily

    realm, name = VARIABLES[variable]["model"]
    expected = "00z" if _sampling(variable) == "00z" else "daily_mean"
    count = config["nens"] if max_members is None else min(max_members, config["nens"])
    years = _file_years(days, variable)
    fields, kinds = [], set()
    for number in range(1, count + 1):
        member = f"EN{number:02d}"
        ts_root = Path(config["path"]) / member / "archive/post" / realm / "180x360_aave/ts"
        kind, _, sampling = daily_source(ts_root, name, member, years)
        if sampling != expected:
            raise ValueError(f"{variable}: {member} {kind} gives {sampling!r} values; the scorecard verifies "
                             f"{expected!r} against {VARIABLES[variable]['reference'][0]}")
        kinds.add(kind)
        field = open_model_daily(ts_root, name, member, years)
        field = _normalize(_on_days(field, days, variable, reindex=True))
        if name == "H2OSOI":
            field = _top_soil(field, land_support_dir)
        fields.append(_convert(field.load(), variable))
    log(f"  {variable}: {count} members ({', '.join(sorted(kinds))})")
    return xr.concat(fields, pd.Index(np.arange(1, count + 1), name="member"))


def _climatology_daily(reference_root, variable, days, years):
    """Calendar-day reference climatology on ``days``, 31-day running mean, and
    each year's unsmoothed values on the same calendar days (year, time, lat, lon)."""
    label, template, name = VARIABLES[variable]["reference"]
    pad = pd.Timedelta(days=15)
    offsets = (days - days[0]).days
    stack = []
    for year in years:
        start = days[0].replace(year=year)
        span = pd.date_range(start - pad, start + pd.Timedelta(days=int(offsets[-1])) + pad)
        files = [Path(reference_root) / template.format(name=name, year=y) for y in _file_years(span, variable)]
        files = [f for f in files if f.is_file()]  # padding may reach before the record
        with xr.open_mfdataset(files, combine="by_coords") as ds:
            field = _on_days(_normalize(ds[name]), span, variable, reindex=True).load()
        stack.append(_convert(field, variable).assign_coords(time=np.arange(len(span))))
    stack = xr.concat(stack, pd.Index(years, name="year"))
    clim = stack.mean("year", skipna=True).rolling(time=31, center=True, min_periods=16).mean()
    inside = slice(15, 15 + len(offsets))
    return clim.isel(time=inside).assign_coords(time=days), stack.isel(time=inside).assign_coords(time=days)


def _climatology_monthly(reference_root, variable, days, years):
    """ERA5 monthly climatology linearly interpolated from mid-month to ``days``,
    and each year's monthly means interpolated the same way."""
    name = VARIABLES[variable]["climatology"][1]
    months = []
    for year in years:
        with xr.open_dataset(Path(reference_root) / _ERA5_MONTHLY.format(year=year)) as ds:
            months.append(_normalize(ds[name]).load().assign_coords(time=np.arange(1, 13)))
    stack = _convert(xr.concat(months, pd.Index(years, name="year")), variable)
    # Cyclic monthly series at mid-month of the forecast year, then to days.
    year = days[0].year
    mid = [pd.Timestamp(year, m, 15) for m in range(1, 13)]
    cyclic = xr.concat([stack.isel(time=[11]), stack, stack.isel(time=[0])], "time")
    cyclic = cyclic.assign_coords(time=[pd.Timestamp(year - 1, 12, 15), *mid, pd.Timestamp(year + 1, 1, 15)])
    stack = cyclic.interp(time=days)
    return stack.mean("year", skipna=True), stack


def above_ground(reference_root, level, days, years):
    """True where the ERA5 surface-pressure climatology of the months of ``days``
    (``years``) is above ``level`` hPa."""
    months = sorted({day.month for day in days})
    stack = []
    for year in years:
        with xr.open_dataset(Path(reference_root) / _ERA5_MONTHLY.format(year=year)) as ds:
            ps = _normalize(ds["PS"])
            stack.append(ps.isel(time=[m - 1 for m in months]).mean("time").load())
    return (xr.concat(stack, "year").mean("year") * 0.01 > level).rename("above_ground")


def read_climatology_days(reference_root, variable, days, years):
    """``(climatology, yearly)``: the climatology on ``days`` and each year's values on them."""
    if VARIABLES[variable]["climatology"] == "daily":
        return _climatology_daily(reference_root, variable, days, years)
    return _climatology_monthly(reference_root, variable, days, years)


def _window_means(daily, valid, windows, min_days):
    means = []
    for days in windows.values():
        sel = daily.sel(time=days).where(valid.sel(time=days))
        count = valid.sel(time=days).sum("time")
        means.append(sel.mean("time", skipna=True).where(count >= min_days))
    return xr.concat(means, pd.Index(list(windows), name="window"))


def scorecard_cache_path(cache_dir, season, variable, experiment):
    return Path(cache_dir) / season / f"{experiment}_{variable}_window_means.nc"


def reference_cache_path(cache_dir, season, variable):
    return Path(cache_dir) / season / f"reference_{variable}_window_means.nc"


def above_ground_cache_path(cache_dir, season, level):
    return Path(cache_dir) / season / f"above_ground_{level}hPa.nc"


def _read_cache(path, settings):
    """The cache at ``path``, loaded, if its recorded settings match ``settings``; else None.

    Each file is opened once. Experiment caches carry the settings on their
    "members" variable, the others on the file; an unreadable file counts as
    missing (rebuilt).
    """
    if not Path(path).is_file():
        return None
    try:
        cached = xr.load_dataset(path)
    except (OSError, ValueError):
        return None
    attrs = {**cached["members"].attrs, **cached.attrs} if "members" in cached else cached.attrs
    return cached if all(str(attrs.get(k)) == str(v) for k, v in settings.items()) else None


def _read_caches(requests, workers):
    """{path: Dataset or None} for ``requests`` {path: settings}, read by ``workers`` processes.

    Reading many cache files is latency-bound on the shared file system, so
    parallel reads are about 6x faster than one at a time with 8 workers. They
    must be processes: the netCDF-C library is not thread-safe (concurrent opens
    from threads crash it).
    """
    paths = list(requests)
    if workers <= 1:
        return {path: _read_cache(path, requests[path]) for path in paths}
    with ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context("fork")) as pool:
        return dict(zip(paths, pool.map(_read_cache, paths, [requests[p] for p in paths])))


def _write(ds, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".nc.tmp")
    ds.to_netcdf(temporary)
    os.replace(temporary, path)


def build_window_means(cases, variables, *, data_root, reference_root, land_support_dir, cache_dir,
                       clim_years=(2001, 2017), clim_exclude_years=(2012,), min_days=10,
                       max_members=None, force_compute=False, workers=8, log=print):
    """Reuse or build window means of the references, climatology, and experiments.

    Returns ``(window_means, registries, failures)`` with
    ``window_means[season][variable] -> Dataset`` holding ``reference``,
    ``climatology`` (window, lat, lon), ``climatology_members`` (year, window,
    lat, lon), ``model`` (experiment, window, lat, lon), the ensemble mean,
    ``members`` (experiment, member, window, lat, lon), NaN past an experiment's
    size, and for pressure-level variables ``above_ground`` (lat, lon).
    Valid caches of a season are read up front, ``workers`` files at a time;
    the above-ground masks are cached too.
    """
    years = climatology_years(variables, clim_years, clim_exclude_years)
    window_means, registries, failures = {}, {}, []
    for season, experiments in cases.items():
        registry = registries[season] = experiment_registry(data_root, season, "fc")
        unknown = set(experiments) - set(registry)
        if unknown:
            raise ValueError(f"{season}: {sorted(unknown)} unavailable; choose from {list(registry)}")
        periods = {registry[e]["period"] for e in experiments}
        if len(periods) != 1:
            raise ValueError(f"{season}: experiments have different periods {periods}")
        period = periods.pop()
        init = pd.Timestamp(f"{period[:4]}-{period[4:6]}-01")
        windows = window_days(init)
        days = pd.DatetimeIndex(sorted({d for w in windows.values() for d in w}))
        settings = dict(cache_version=CACHE_VERSION, init_date=str(init.date()), min_days=min_days,
                        clim_years=",".join(map(str, years)))
        surface_settings = dict(cache_version=CACHE_VERSION, clim_years=",".join(map(str, years)),
                                months=",".join(str(m) for m in sorted({d.month for d in days})))

        def model_settings(variable, experiment):
            members = registry[experiment]["nens"] if max_members is None else min(max_members, registry[experiment]["nens"])
            recorded = {**settings, "ensemble_members": members}
            if _sampling(variable) == "6hourly_mean":
                recorded["daily_construction"] = SIX_HOURLY_DAILY
            return recorded

        # Read every cache the season needs at once (in parallel); None marks one to build.
        requests = {}
        if not force_compute:
            for variable in variables:
                requests[reference_cache_path(cache_dir, season, variable)] = settings
                for experiment in experiments:
                    requests[scorecard_cache_path(cache_dir, season, variable, experiment)] = \
                        model_settings(variable, experiment)
            for level in sorted({int(v[1:]) for v in variables if _PRESSURE_LEVEL.fullmatch(v)}):
                requests[above_ground_cache_path(cache_dir, season, level)] = surface_settings
        cached = _read_caches(requests, workers) if requests else {}
        surfaces = {}  # level -> above-ground mask, shared by the season's variables
        for variable in variables:
            try:
                ref_path = reference_cache_path(cache_dir, season, variable)
                daily_reference = None
                reference = cached.get(ref_path)
                if reference is not None:
                    log(f"[REUSE] {season} {variable} reference")
                else:
                    daily_reference = read_reference_days(reference_root, variable, days)
                    valid = daily_reference.notnull()
                    clim, clim_years_days = read_climatology_days(reference_root, variable, days, years)
                    reference = xr.Dataset({
                        "reference": _window_means(daily_reference, valid, windows, min_days),
                        "climatology": _window_means(clim, valid, windows, min_days),
                        "climatology_members": _window_means(clim_years_days, valid, windows, min_days),
                        "valid_days": valid.astype("int8"),
                    }, attrs={**settings, "reference_label": VARIABLES[variable]["reference"][0]})
                    _write(reference, ref_path)
                    log(f"[WRITE] {season} {variable} reference")
                valid = reference["valid_days"].astype(bool)
                models = []
                for experiment in experiments:
                    path = scorecard_cache_path(cache_dir, season, variable, experiment)
                    if cached.get(path) is not None:
                        members = cached[path]["members"]
                        log(f"[REUSE] {season} {variable} {experiment}")
                    else:
                        daily = read_model_days(registry[experiment], variable, days,
                                                land_support_dir=land_support_dir,
                                                max_members=max_members, log=log)
                        members = _window_means(daily, valid, windows, min_days).rename("members")
                        members.attrs.update(model_settings(variable, experiment))
                        _write(members.to_dataset(), path)
                        log(f"[WRITE] {season} {variable} {experiment}")
                    models.append(members.expand_dims(experiment=[experiment]))
                # Experiments differ in size; missing members are NaN.
                members = xr.concat(models, "experiment", join="outer")
                members.attrs = {}
                ds = reference[["reference", "climatology", "climatology_members"]].assign(
                    model=members.mean("member", skipna=True), members=members)
                if _PRESSURE_LEVEL.fullmatch(variable):
                    level = int(variable[1:])
                    if level not in surfaces:
                        path = above_ground_cache_path(cache_dir, season, level)
                        if cached.get(path) is not None:
                            surfaces[level] = cached[path]["above_ground"].astype(bool)
                        else:
                            surfaces[level] = above_ground(reference_root, level, days, years)
                            mask = surfaces[level].astype("int8").to_dataset()
                            mask.attrs.update(surface_settings)
                            _write(mask, path)
                            log(f"[WRITE] {season} above-ground mask {level} hPa")
                    ds["above_ground"] = surfaces[level].reindex_like(ds.reference, method="nearest", tolerance=0.01)
                window_means.setdefault(season, {})[variable] = ds
            except Exception as error:  # keep going; report at the end
                failures.append((season, variable, f"{type(error).__name__}: {error}"))
                log(f"[FAIL] {season} {variable}: {type(error).__name__}: {error}")
    return window_means, registries, failures


def _weighted_pcc(x, y, w):
    valid = x.notnull() & y.notnull() & (w > 0)
    w = w.where(valid, 0.0)
    total = w.sum(("lat", "lon"))
    mx = (x.fillna(0) * w).sum(("lat", "lon")) / total
    my = (y.fillna(0) * w).sum(("lat", "lon")) / total
    dx, dy = (x - mx).fillna(0), (y - my).fillna(0)
    cov = (dx * dy * w).sum(("lat", "lon"))
    return cov / np.sqrt((dx ** 2 * w).sum(("lat", "lon")) * (dy ** 2 * w).sum(("lat", "lon")))


def _weighted_error(x, y, w, remove_mean_bias=False):
    """Area-mean bias and RMSE of ``x - y``; the RMSE is centered if ``remove_mean_bias``."""
    valid = x.notnull() & y.notnull() & (w > 0)
    w = w.where(valid, 0.0)
    total = w.sum(("lat", "lon"))
    error = (x - y).fillna(0)
    bias = (error * w).sum(("lat", "lon")) / total
    if remove_mean_bias:
        error = error - bias
    return bias, np.sqrt((error ** 2 * w).sum(("lat", "lon")) / total)


def _fair_crps(members, obs, dim):
    """Fair CRPS of an ensemble along ``dim`` (NaN members ignored; NaN with fewer than 2).

    CRPS = mean|x_i - y| - sum_ij |x_i - x_j| / (2 n (n - 1)), with the pair sum
    from the sorted members x_(k): sum_ij |x_i - x_j| = 2 sum_k (2k - n - 1) x_(k).
    """
    def crps(x, y):
        n = np.isfinite(x).sum(-1)
        ordered = np.sort(x, axis=-1)  # NaN last
        k = np.arange(1, x.shape[-1] + 1)
        weight = np.where(k <= n[..., None], 2 * k - n[..., None] - 1, 0)
        pairs = (weight * np.nan_to_num(ordered)).sum(-1)
        with np.errstate(invalid="ignore", divide="ignore"):
            error = np.nanmean(np.abs(x - y[..., None]), axis=-1)
            return np.where(n >= 2, error - pairs / (n * (n - 1)), np.nan)

    return xr.apply_ufunc(crps, members, obs, input_core_dims=[[dim], []], dask="forbidden")


def _weighted_mean_pair(x, y, w):
    """Area means of ``x`` and ``y`` over the points where both are valid."""
    valid = x.notnull() & y.notnull() & (w > 0)
    w = w.where(valid, 0.0)
    total = w.sum(("lat", "lon"))
    return (x.fillna(0) * w).sum(("lat", "lon")) / total, (y.fillna(0) * w).sum(("lat", "lon")) / total


def window_scores(window_means, *, regions, landmask_file, remove_mean_bias=False):
    """Long table of ``rmse``, ``rmse_clim``, ``msss``, ``bias``, ``crps``,
    ``crps_clim``, ``crpss``, and ``pcc`` per season/variable/experiment/region/window."""
    rows = []
    masks = None
    for season, by_variable in window_means.items():
        for variable, ds in by_variable.items():
            if masks is None:
                masks = region_masks(ds.lat, ds.lon, regions, landmask_file)
                land = region_masks(ds.lat, ds.lon, ("land",), landmask_file).isel(region=0, drop=True)
            component = VARIABLES[variable]["component"]
            cos = np.cos(np.deg2rad(ds.lat)).clip(min=0)
            obs_anomaly = ds.reference - ds.climatology
            model_anomaly = ds.model - ds.climatology
            crps_model = _fair_crps(ds.members, ds.reference, "member")
            crps_reference = _fair_crps(ds.climatology_members, ds.reference, "year")
            for region in regions:
                mask = masks.sel(region=region, drop=True)
                if component == "lnd":
                    mask = mask & land
                if "above_ground" in ds:
                    mask = mask & ds.above_ground
                weights = (cos * mask).broadcast_like(mask).astype(float)
                bias, rmse = _weighted_error(ds.model, ds.reference, weights, remove_mean_bias)
                _, rmse_clim = _weighted_error(ds.climatology, ds.reference, weights, remove_mean_bias)
                msss = 1.0 - (rmse / rmse_clim) ** 2
                crps, crps_clim = _weighted_mean_pair(crps_model, crps_reference, weights)
                crpss = 1.0 - crps / crps_clim
                pcc = _weighted_pcc(model_anomaly, obs_anomaly, weights)
                for experiment in ds.experiment.values:
                    for window in ds.window.values:
                        at = dict(experiment=experiment, window=window)
                        rows.append(dict(season=season, component=component, variable=variable,
                                         units=VARIABLES[variable]["units"], experiment=str(experiment),
                                         region=region, window=str(window),
                                         rmse=float(rmse.sel(at)),
                                         rmse_clim=float(rmse_clim.sel(window=window)),
                                         msss=float(msss.sel(at)), bias=float(bias.sel(at)),
                                         crps=float(crps.sel(at)), crps_clim=float(crps_clim.sel(at)),
                                         crpss=float(crpss.sel(at)),
                                         pcc=float(pcc.sel(at)),
                                         mean_bias_removed=bool(remove_mean_bias)))
    return pd.DataFrame(rows)


# Colors: sequential blue for PCC; diverging red (worse than climatology) - gray - blue (better).
_PCC_CMAP = ListedColormap(["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
                            "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281"])
_SKILL_CMAP = LinearSegmentedColormap.from_list(
    "msss", ["#a32a2a", "#e34948", "#f2aaa6", "#f0efec", "#9ec5f4", "#3987e5", "#104281"])
_MISSING = "#d9d8d4"
_TRIANGLE = {  # cell corners with y growing downward: TL, TR, BR, BL, center
    0: lambda x, y: [(x, y), (x + 1, y), (x + 0.5, y + 0.5)],          # top
    1: lambda x, y: [(x + 1, y), (x + 1, y + 1), (x + 0.5, y + 0.5)],  # right
    2: lambda x, y: [(x + 1, y + 1), (x, y + 1), (x + 0.5, y + 0.5)],  # bottom
    3: lambda x, y: [(x, y + 1), (x, y), (x + 0.5, y + 0.5)],          # left
}


def _group_rows(groups, available, gap):
    """Row positions (top edges) for the variables of each group, with ``gap`` rows between groups."""
    rows, spans, y = [], [], 0.0
    for name, variables in groups.items():
        present = [v for v in variables if v in available]
        if not present:
            continue
        if rows:
            y += gap
        spans.append((name, y, y + len(present)))
        rows.extend((v, y + i) for i, v in enumerate(present))
        y += len(present)
    return rows, spans, y


def plot_window_scorecard_paper(scores, *, region, cases, groups, labels=None, season_names=None, title="",
                                width=7.2, row_height=0.18, fontsize=7, key_size=0.60):
    """Every season and variable in one figure: PCC, MSSS, and CRPSS blocks side by side.

    ``cases`` {season: experiments}; inside each block the seasons sit side by
    side (``season_names`` {season: name} heads them, e.g. "January"; None:
    the season key), one column per experiment. Rows are the variables of
    ``groups`` ({"Surface": (...), "Upper air": (...), "Land": (...)}),
    bracketed and separated by a gap; inside a group, rows of a new
    pressure-level variable (U, V, T, Q, Z) start below a thin line. Experiment
    names (``labels`` {experiment: label}, e.g. "CTRL10-S0"; None: the
    experiment key) sit under the panels, season names above them. Each cell holds the
    four lead-window triangles (key at the bottom left). ``width`` is
    the figure width in inches (7.2: double column), ``row_height`` inches
    per variable row, ``key_size`` the side of the lead-week key in inches;
    margins and offsets scale with ``fontsize`` (sized at 7 pt), the left
    margin with the longest variable name. Save as PDF.
    """
    table = scores[scores["region"] == region]
    rows, spans, n_rows = _group_rows(groups, set(table["variable"]), gap=0.5)
    present = set(zip(table["season"], table["experiment"]))
    columns = [(season, e) for season, experiments in cases.items() for e in experiments if (season, e) in present]
    seasons = list(dict.fromkeys(season for season, _ in columns))
    season_gap = 0.35  # cell widths between seasons inside a block
    xs, x = [], 0.0
    for k, (season, _) in enumerate(columns):
        if k and season != columns[k - 1][0]:
            x += season_gap
        xs.append(x)
        x += 1
    span = x
    lookup = table.set_index(["variable", "season", "experiment", "window"])
    pcc_cmap = _PCC_CMAP.with_extremes(under="#f0efec", bad=_MISSING)
    pcc_norm = BoundaryNorm(np.linspace(0, 1, 13), _PCC_CMAP.N)
    skill_cmap = _SKILL_CMAP.with_extremes(bad=_MISSING)
    skill_norm = BoundaryNorm(np.arange(-1.1, 1.11, 0.2), 256, extend="min")
    centered = bool(table["mean_bias_removed"].any()) if "mean_bias_removed" in table else False
    metrics = [("pcc", "Anomaly PCC", pcc_cmap, pcc_norm),
               ("msss", "MSSS" + (" (mean bias removed)" if centered else ""), skill_cmap, skill_norm),
               ("crpss", "CRPSS", skill_cmap, skill_norm)]
    style = {"font.size": fontsize, "axes.titlesize": fontsize + 1, "xtick.labelsize": fontsize - 1,
             "ytick.labelsize": fontsize, "pdf.fonttype": 42, "ps.fonttype": 42}
    with plt.rc_context(style):
        scale = fontsize / 7  # margins below are sized for 7 pt text
        top = 0.62 * scale  # inches for headings and seasons
        names = max(len((labels or {}).get(e, e)) for _, e in columns)  # rotated experiment names under the panels
        bottom = max(0.62 * scale + 0.62 * names * (fontsize - 1) / 72,  # experiment names and colorbars
                     key_size + 0.62 * scale)  # the window key (with its labels) in the corner left of them
        height = top + bottom + row_height * n_rows
        fig = plt.figure(figsize=(width, height))
        longest = max(len(v) for v, _ in rows)  # variable names left of the rows, group names left of them
        left = (0.45 * scale + 0.6 * longest * fontsize / 72) / width
        right, gap = 0.08 / width, 0.12 * scale / width
        block = (1 - left - right - 2 * gap) / 3
        y0, y1 = bottom / height, 1 - top / height
        collections = []
        for col, (metric, name, cmap, norm) in enumerate(metrics):
            ax = fig.add_axes([left + col * (block + gap), y0, block, y1 - y0])
            polygons, values = [], []
            for variable, y in rows:
                for (season, experiment), x0 in zip(columns, xs):
                    for index, window in enumerate(WINDOWS):
                        polygons.append(_TRIANGLE[index](x0, y))
                        key = (variable, season, experiment, window)
                        values.append(lookup[metric].get(key, np.nan) if key in lookup.index else np.nan)
            collection = PolyCollection(polygons, array=np.ma.masked_invalid(np.asarray(values, float)),
                                        cmap=cmap, norm=norm, edgecolors="#fcfcfb", linewidths=0.4)
            ax.add_collection(collection)
            collections.append(collection)
            ax.set_xlim(0, span)
            ax.set_ylim(n_rows, 0)
            ax.set_xticks(np.asarray(xs) + 0.5, [(labels or {}).get(e, e) for _, e in columns], rotation=90)
            ax.set_yticks([y + 0.5 for _, y in rows], [v for v, _ in rows] if col == 0 else [])
            ax.tick_params(length=0, pad=2)
            for spine in ax.spines.values():
                spine.set_visible(False)
            # Season names above the experiment labels, with a rule over their columns.
            trans = ax.get_xaxis_transform()
            label_top = 1 + 0.08 * scale / (y1 - y0) / height  # just above the panel, in axes units
            for season in seasons:
                spots = [x0 for (s, _), x0 in zip(columns, xs) if s == season]
                ax.plot([spots[0] + 0.08, spots[-1] + 0.92], [label_top] * 2, color="#2b2b29", lw=0.5,
                        transform=trans, clip_on=False)
                ax.text((spots[0] + spots[-1] + 1) / 2, label_top + 0.02 * scale / (y1 - y0) / height,
                        (season_names or {}).get(season, season), ha="center", va="bottom", transform=trans)
            ax.set_title(f"({chr(97 + col)}) {name}", fontweight="bold", pad=0.30 * scale * 72)
            for (v_prev, _), (variable, y) in zip(rows, rows[1:]):
                same_group = any(a <= y < b and a <= y - 1 < b for _, a, b in spans)
                if same_group and _PRESSURE_LEVEL.fullmatch(variable) and variable[0] != v_prev[0]:
                    for season in seasons:
                        spots = [x0 for (s, _), x0 in zip(columns, xs) if s == season]
                        ax.plot([spots[0], spots[-1] + 1], [y, y], color="#2b2b29", lw=0.6, clip_on=False)
        # Group brackets left of the variable names.
        axes0 = fig.axes[0]
        trans = axes0.get_yaxis_transform()
        for name, a, b in spans:
            x = -(left * width - 0.2 * scale) / (block * width)  # near the figure edge, in block units
            axes0.plot([x + 0.04, x, x, x + 0.04], [a + 0.1, a + 0.1, b - 0.1, b - 0.1], color="#2b2b29",
                       lw=0.6, transform=trans, clip_on=False)
            axes0.text(x - 0.03, (a + b) / 2, name, rotation=90, ha="right", va="center", transform=trans,
                       fontweight="bold")
        # Colorbars under the blocks.
        cb_h, cb_y = 0.09 * scale / height, 0.36 * scale / height
        cax = fig.add_axes([left + 0.1 * block, cb_y, 0.8 * block, cb_h])
        fig.colorbar(collections[0], cax=cax, orientation="horizontal", extend="min",
                     ticks=np.linspace(0, 1, 5 if fontsize <= 8 else 3)).set_label("PCC (below 0: gray)", labelpad=2)
        cax = fig.add_axes([left + block + gap + 0.1 * block, cb_y, 1.8 * block + gap, cb_h])
        fig.colorbar(collections[1], cax=cax, orientation="horizontal", extend="min",
                     ticks=np.arange(-1.0, 1.01, 0.5)).set_label(
            "MSSS, CRPSS (blue: beats climatology)", labelpad=2)
        for cax in fig.axes[-2:]:
            cax.tick_params(length=2, pad=1.5)
        # Window key in the empty corner below the variable names, its name above it.
        key_x = 0.55 * left * width - key_size / 2  # inches from the left edge (center in the margin)
        key_y = bottom / 2 - 0.10 * scale - key_size / 2  # inches from the bottom edge
        kax = fig.add_axes([key_x / width, key_y / height, key_size / width, key_size / height])
        for index in range(len(WINDOWS)):
            kax.add_patch(plt.Polygon(_TRIANGLE[index](0, 0), facecolor="#f0efec", edgecolor="#6f6e69", lw=0.6))
        places = {0: (0.5, -0.06, "center", "bottom"), 1: (1.06, 0.5, "left", "center"),
                  2: (0.5, 1.06, "center", "top"), 3: (-0.06, 0.5, "right", "center")}
        for index, window in enumerate(WINDOWS):
            x, y, ha, va = places[index]
            kax.text(x, y, window.replace("week", ""), ha=ha, va=va, fontsize=fontsize - 0.5)
        kax.text(0.5, -(0.06 + 0.17 * scale) / key_size, "Lead weeks", ha="center", va="bottom", fontsize=fontsize,
                 fontweight="bold")
        kax.set_xlim(0, 1)
        kax.set_ylim(1, 0)
        kax.axis("off")
        if title:
            fig.suptitle(title, fontsize=fontsize + 2, y=1 - 0.05 / height, va="top")
    return fig
