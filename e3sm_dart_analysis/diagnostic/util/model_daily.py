"""Daily model fields for hindcast evaluation, with the matching reference sampling.

Each member's ``.../180x360_aave/ts`` directory offers a variable in some of:

* ``daily_mean``: 00Z-24Z daily means built from the 6-hourly output (interval
  means averaged over the day; instantaneous PS/PSL as the mean of the 00, 06,
  12, and 18Z states, the same construction as the ERA5 daily files), stamped
  at 12Z with time_bnds [day, day + 1];
* ``daily_max`` / ``daily_min``: 00Z-24Z daily extremes (PRECTMX, TREFHTMX;
  TREFHTMN), stamped like ``daily_mean``;
* ``6hourly``: the 6-hourly series those means come from;
* ``daily_00z``: once-daily 00Z instantaneous states (pressure-level U/V/T/Q/Z/
  OMEGA, UBOT/VBOT/TBOT/ZBOT, CAPE/CIN, and 00Z values of PSL, TREFHT, ...),
  stamped at their instant; the first record is the lead-0 initial state;
* ``h1`` / ``daily`` (older layout, still used by the DA runs and ELM):
  the EAM h1 daily series, daily means of a few fields and 00Z states of the
  rest, each stamped at the end of its interval.

Rule used by every fcst workflow (``daily_source``):

1. If a daily mean can be had (``daily_mean``, the ``daily_max``/``daily_min``
   extremes, else built from ``6hourly``, else a daily-mean ``h1``/``daily``
   file), use it and compare with a daily-mean reference.
2. Otherwise use the 00Z state (``daily_00z``, else a ``time: point``
   ``h1``/``daily`` file), labelled by its instant, and compare with the
   reference sampled at the same instant (ERA5 6-hourly at 00Z), never with a
   daily mean. The lead-0 initial state is dropped.

Daily means are labelled by the 00Z start of the day they cover;
``DataArray.attrs["sampling"]`` is ``"daily_mean"`` or ``"00z"``.

Workflows that must combine fields as daily values (e.g. the RMM index) use
``open_model_daily_estimate``: daily means where they exist, else the mean of
the two 00Z states bracketing the day (``"00z_pair_mean"``).
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import xarray as xr

from util.time_labels import label_by_interval_start

SAMPLINGS = ("daily_mean", "00z")
_HOUR = np.timedelta64(3600, "s")
_DAY = np.timedelta64(86400, "s")


_SOURCES = ("daily_mean", "daily_max", "daily_min", "6hourly", "daily_00z", "h1", "daily")
_DAILY_AGGREGATES = ("daily_mean", "daily_max", "daily_min", "6hourly")


def _source_files(ts_root, kind, variable, member, years):
    files = [Path(ts_root) / kind / f"{variable}.{member}.{year}.nc" for year in years]
    return files if all(f.is_file() for f in files) else None


def daily_source(ts_root, variable, member, years):
    """``(kind, files, sampling)`` for ``variable`` under a member's ``ts_root``.

    ``kind`` is the ts subdirectory used (``_SOURCES`` order: daily means first,
    00Z states last); ``sampling`` is ``"daily_mean"`` or ``"00z"``.
    """
    found = {}
    for kind in _SOURCES:
        files = _source_files(ts_root, kind, variable, member, years)
        if files is None:
            continue
        if kind in _DAILY_AGGREGATES:
            return kind, files, "daily_mean"
        if kind == "daily_00z":
            found.setdefault("00z", (kind, files))
            continue
        with xr.open_dataset(files[0]) as ds:
            point = "time: point" in str(ds[variable].attrs.get("cell_methods", ""))
            name = ds.time.attrs.get("bounds", "time_bnds")
            later = slice(1, None)  # record 0 may be the zero-width lead-0 state
            if not point and name in ds:
                bounds = ds[name].values[later]
                point = bool(len(bounds)) and bool((bounds[:, 1] == bounds[:, 0]).all())  # "mean" of one step
            hours = set(ds.time.dt.hour.values[later].tolist())
        if not point:
            return kind, files, "daily_mean"
        if hours != {0}:
            raise ValueError(f"{variable}: {kind} holds instantaneous values at {sorted(hours)} UTC, not 00Z "
                             f"(e.g. the 09Z first-step snapshots of the DA runs); no matching reference sampling")
        found.setdefault("00z", (kind, files))
    if "00z" in found:
        kind, files = found["00z"]
        return kind, files, "00z"
    raise FileNotFoundError(f"No {'/'.join(_SOURCES)} {variable} files for {member} {list(years)} under {ts_root}")


def model_daily_files(ts_root, variable, member, years):
    """The files ``open_model_daily`` reads for ``variable`` (see ``daily_source``)."""
    return daily_source(ts_root, variable, member, years)[1]


def daily_from_6hourly(ds, variable):
    """00Z-00Z daily means of a 6-hourly series; only complete days are kept."""
    field = ds[variable]
    point = "time: point" in str(field.attrs.get("cell_methods", ""))
    if point:
        times = ds.time.values.astype("datetime64[s]")
        days = np.unique(times.astype("datetime64[D]")).astype("datetime64[s]")
        rows, kept = [], []
        for day in days:
            hits = [np.flatnonzero(times == day + h * _HOUR) for h in (0, 6, 12, 18)]
            if all(len(h) == 1 for h in hits):
                weights = np.zeros(len(times))
                weights[[int(h[0]) for h in hits]] = 0.25
                rows.append(weights)
                kept.append(day)
    else:
        field = label_by_interval_start(ds, field)
        name = ds.time.attrs.get("bounds", "time_bnds")
        bounds = ds[name].values.astype("datetime64[s]")
        width = (bounds[:, 1] - bounds[:, 0]) / _HOUR
        keep = width >= 0.5 * np.median(width[width > 0])
        start, width = bounds[keep, 0], width[keep]
        days = np.unique(start.astype("datetime64[D]")).astype("datetime64[s]")
        rows, kept = [], []
        for day in days:
            inside = (start >= day) & (start + width * _HOUR <= day + _DAY)
            if abs(width[inside].sum() - 24.0) > 1e-6:
                continue
            rows.append(np.where(inside, width, 0.0) / 24.0)
            kept.append(day)
    if not kept:
        raise ValueError(f"{variable}: no complete day in the 6-hourly series")
    weights = xr.DataArray(np.array(rows), dims=("day", "time"))
    daily = xr.dot(weights, field.fillna(0), dim="time").where(
        xr.dot(weights, field.notnull().astype(float), dim="time") > 0.999)
    daily = daily.rename(day="time").assign_coords(time=np.array(kept).astype("datetime64[ns]"))
    daily.attrs = dict(field.attrs, sampling="daily_mean")
    daily.name = variable
    return daily.transpose("time", ...)


def open_model_daily(ts_root, variable, member, years):
    """Daily model field of one member: daily means (rule 1) or 00Z states (rule 2)."""
    kind, files, sampling = daily_source(ts_root, variable, member, years)
    ds = xr.open_mfdataset(files, combine="by_coords")
    if kind == "6hourly":
        return daily_from_6hourly(ds, variable)
    field = ds[variable]
    if sampling == "00z":
        field = field.isel(time=slice(1, None))  # drop the lead-0 initial state
        field.attrs = dict(field.attrs, sampling="00z")
        return field
    # Daily means: h1/daily stamp them at the day end, daily_mean/max/min at 12Z; label by the day.
    field = label_by_interval_start(ds, field)
    field = field.assign_coords(time=field.time.dt.floor("D"))
    field.attrs = dict(field.attrs, sampling="daily_mean")
    return field



def open_model_daily_estimate(ts_root, variable, member, years):
    """Daily values of one member for workflows that combine fields day by day.

    Daily means (rule 1) as ``open_model_daily``. A field held only as 00Z
    states (e.g. U850, U200) is estimated for day d as the mean of the 00Z
    states on d and d + 1, centred on the same 00Z-24Z window; day 1 uses the
    lead-0 initial state, and a day missing either state is dropped.
    ``attrs["sampling"]`` is ``"daily_mean"`` or ``"00z_pair_mean"``.
    """
    kind, files, sampling = daily_source(ts_root, variable, member, years)
    if sampling == "daily_mean":
        return open_model_daily(ts_root, variable, member, years)
    ds = xr.open_mfdataset(files, combine="by_coords")
    field = ds[variable]
    field = field.sel(time=field.time.dt.hour == 0)  # 00Z states, lead 0 included
    days = field.time.dt.floor("D").values
    states = field.assign_coords(time=days)
    following = states.shift(time=-1)
    gap = np.append(np.diff(days) == _DAY.astype("timedelta64[ns]"), False)  # next state is 24 h later
    daily = (0.5 * (states + following)).isel(time=np.flatnonzero(gap))
    daily.attrs = dict(field.attrs, sampling="00z_pair_mean",
                       cell_methods="time: mean (of the 00Z states at the start and end of the day)")
    daily.name = variable
    return daily


def reference_at_sampling(daily_mean, six_hourly, sampling):
    """Reference matching a model field's ``sampling``.

    ``daily_mean``: the reference daily means, labelled by day (00Z).
    ``00z``: the 6-hourly reference at 00Z instants. Pass whichever product
    the sampling needs; the other may be None.
    """
    if sampling not in SAMPLINGS:
        raise ValueError(f"sampling must be one of {SAMPLINGS}")
    if sampling == "daily_mean":
        return daily_mean.assign_coords(time=daily_mean.time.dt.floor("D"))
    return six_hourly.sel(time=six_hourly.time.dt.hour == 0)


# --------------------------------------------------------------------------
# File-list interface (workflows that collect ts/<kind> file paths, e.g. S2S)
# --------------------------------------------------------------------------
def _parse_member_files(files, variable):
    """(ts_root, member, years) from ``.../ts/<kind>/<variable>.<member>.<year>.nc`` paths."""
    paths = [Path(f) for f in files]
    if not paths:
        raise ValueError(f"No files for {variable}")
    members, years = set(), []
    for path in paths:
        stem = path.name[len(variable) + 1:-3] if path.name.startswith(variable + ".") else ""
        member, _, year = stem.partition(".")
        if not member or not year.isdigit():
            raise ValueError(f"Not a <variable>.<member>.<year>.nc file: {path}")
        members.add(member)
        years.append(int(year))
    if len(members) != 1:
        raise ValueError(f"Files mix members {sorted(members)}")
    return paths[0].parent.parent, members.pop(), sorted(years)


def model_sampling(files, variable):
    """``"daily_mean"`` or ``"00z"``: how ``open_model_daily_files`` will sample ``variable``."""
    ts_root, member, years = _parse_member_files(files, variable)
    return daily_source(ts_root, variable, member, years)[2]


def open_model_daily_files(files, variable):
    """``open_model_daily`` for a list of one member's ts/<kind> files (any kind; the rule picks the source)."""
    ts_root, member, years = _parse_member_files(files, variable)
    return open_model_daily(ts_root, variable, member, years)


_ERA5_DAILY_NAME = "{variable}_{year}01_{year}12.nc"
_ERA5_6HOURLY_NAME = "ERA5.6hourly.en00.{variable}.{year}01-{year}12.nc"


def reference_files(files, variable, sampling):
    """Reference files for ``sampling``: the daily files as given, or for ``"00z"``
    the ERA5 6-hourly files of the same variable and years (select 00Z after opening).

    Only ERA5 provides instantaneous 6-hourly fields; other products cannot
    verify ``"00z"`` model values.
    """
    if sampling == "daily_mean":
        return list(files)
    out = []
    for f in map(Path, files):
        year = next((y for y in range(1900, 2101) if f.name == _ERA5_DAILY_NAME.format(variable=variable, year=y)), None)
        if f.parent.name != "daily" or year is None:
            raise ValueError(f"{variable}: model values are 00Z states, but {f} is not an ERA5 daily file; "
                             "verify them against an instantaneous (6-hourly) reference at 00Z")
        six_hourly = f.parent.parent / "6hourly" / _ERA5_6HOURLY_NAME.format(variable=variable, year=year)
        if not six_hourly.is_file():
            raise FileNotFoundError(f"{variable}: no ERA5 6-hourly file for 00Z sampling: {six_hourly}")
        out.append(six_hourly)
    return out


def select_00z(field):
    """The 00Z records of a 6-hourly reference."""
    return field.sel(time=field.time.dt.hour == 0)


# --------------------------------------------------------------------------
# Pressure levels below ground
# --------------------------------------------------------------------------
_LEVEL_VARIABLE = re.compile(r"(U|V|T|Q|Z|OMEGA)(\d{3,4})")
ABOVE_GROUND_YEARS = range(2001, 2018)  # same ERA5 PS climatology as the window scorecard
_ERA5_MONTHLY_NAME = "ERA5_analysis_monthly_{year}.nc"
_above_ground_cache = {}


def pressure_level(variable):
    """Pressure level (hPa) of a level variable such as ``U850``, else None."""
    match = _LEVEL_VARIABLE.fullmatch(str(variable))
    return int(match.group(2)) if match else None


def above_ground(era5_monthly_dir, level, months, years=ABOVE_GROUND_YEARS):
    """(lat, lon) True where the ERA5 surface-pressure climatology of ``months`` exceeds ``level`` hPa.

    Elsewhere (e.g. 850 hPa over Tibet, the Andes, Greenland, Antarctica) model
    and ERA5 level values are extrapolations below the surface.
    """
    key = (str(era5_monthly_dir), level, tuple(sorted(months)), tuple(years))
    if key not in _above_ground_cache:
        stack = []
        for year in years:
            with xr.open_dataset(Path(era5_monthly_dir) / _ERA5_MONTHLY_NAME.format(year=year)) as ds:
                ps = ds["PS"].isel(time=[m - 1 for m in sorted(months)]).mean("time").load()
            stack.append(ps)
        ps = xr.concat(stack, "year").mean("year")
        ps = ps.rename({k: v for k, v in (("latitude", "lat"), ("longitude", "lon")) if k in ps.dims})
        ps = ps.assign_coords(lon=ps.lon % 360).sortby("lat").sortby("lon")
        _above_ground_cache[key] = (ps * 0.01 > level).rename("above_ground")
    return _above_ground_cache[key]


def mask_below_ground(field, variable, reference_files_used, *, lat="lat", lon="lon"):
    """``field`` with gridpoints where ``variable``'s level is below ground set to NaN.

    Non-level variables are returned unchanged. The ERA5 monthly directory is
    found next to the ERA5 daily/6hourly reference files.
    """
    level = pressure_level(variable)
    if level is None:
        return field
    era5_root = Path(list(reference_files_used)[0]).parent.parent
    months = sorted(set(int(m) for m in field.time.dt.month.values))
    mask = above_ground(era5_root / "monthly", level, months)
    mask = mask.rename({"lat": lat, "lon": lon}) if (lat, lon) != ("lat", "lon") else mask
    target = field[lon] % 360
    mask = mask.reindex({lat: field[lat].values, lon: target.values}, method="nearest", tolerance=0.01)
    mask = mask.assign_coords({lon: field[lon].values})
    return field.where(mask.fillna(False).astype(bool))
