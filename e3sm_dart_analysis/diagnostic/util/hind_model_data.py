"""Model fields of the 2012-2013 AMIP hindcasts (``hind_*`` workflows).

All readers open lazily through ``util.dask_helpers`` and return fields on a
common time convention so model and reference can be aligned by label:

* daily means are labelled by the 00Z start of the day they cover (zppy
  stamps them at the interval end, Jan 1 -> Jan 2 00Z; see
  ``util.time_labels``), on ``datetime64`` (the noleap model calendar maps
  onto real dates; a reference Feb 29 simply has no model partner);
* monthly means are labelled ``YYYY-MM-16`` on ``cftime.DatetimeNoLeap``
  (``monthly_label``), the convention used for every monthly reference.
"""
from __future__ import annotations

import re
from pathlib import Path

import cftime
import numpy as np
import xarray as xr

from configs.hind_experiment_config import (
    DEFAULT_INPUT_DATA_ROOT,
    END_DATE,
    START_DATE,
    output_dir,
)
from util.dask_helpers import open_mfdataset_lazy
from util.time_labels import label_by_interval_start

_CFTIME = xr.coders.CFDatetimeCoder(use_cftime=True)


# --------------------------------------------------------------------------
# Time helpers shared with the reference readers
# --------------------------------------------------------------------------
def to_datetime64(values) -> np.ndarray:
    """cftime (any calendar whose dates exist in the real calendar) or datetime64 -> datetime64[ns]."""
    values = np.asarray(values)
    if np.issubdtype(values.dtype, np.datetime64):
        return values.astype("datetime64[ns]")
    return np.array(
        [np.datetime64(f"{t.year:04d}-{t.month:02d}-{t.day:02d}T{t.hour:02d}:{t.minute:02d}:{t.second:02d}")
         for t in values.ravel()],
        dtype="datetime64[ns]",
    ).reshape(values.shape)


def monthly_label(da: xr.DataArray, time: str = "time") -> xr.DataArray:
    """Relabel a monthly series as ``YYYY-MM-16`` on ``cftime.DatetimeNoLeap``.

    Model months are stamped at their end, references at the middle; labelling
    by year and month after labelling by interval start gives both the same axis.
    """
    vals = da[time].values
    if np.issubdtype(np.asarray(vals).dtype, np.datetime64):
        years = vals.astype("datetime64[Y]").astype(int) + 1970
        months = vals.astype("datetime64[M]").astype(int) % 12 + 1
    else:
        years = np.array([t.year for t in vals])
        months = np.array([t.month for t in vals])
    labels = np.array([cftime.DatetimeNoLeap(int(y), int(m), 16) for y, m in zip(years, months)])
    out = da.assign_coords({time: labels})
    if np.unique(out[time].values).size != out.sizes[time]:
        raise ValueError(f"{da.name}: duplicate months after monthly labelling")
    return out


def daily_label(da: xr.DataArray, time: str = "time") -> xr.DataArray:
    """Floor a daily series to its 00Z day label on datetime64 and drop duplicate days."""
    days = to_datetime64(da[time].values).astype("datetime64[D]").astype("datetime64[ns]")
    da = da.assign_coords({time: days}).sortby(time)
    _, keep = np.unique(da[time].values, return_index=True)
    return da.isel({time: np.sort(keep)}) if keep.size != da.sizes[time] else da


def standardize_latlon(da: xr.DataArray) -> xr.DataArray:
    """``lat``/``lon`` names, longitude in [0, 360), both ascending; singleton levels dropped."""
    rename = {k: v for k, v in (("latitude", "lat"), ("longitude", "lon")) if k in da.dims}
    if rename:
        da = da.rename(rename)
    for lev in ("level", "lev", "plev"):
        if lev in da.dims and da.sizes[lev] == 1:
            da = da.isel({lev: 0}, drop=True)
    if "lon" in da.dims and float(da.lon.min()) < 0:
        da = da.assign_coords(lon=da.lon % 360)
    for dim in ("lat", "lon"):
        if dim in da.dims:
            da = da.sortby(dim)
    return da


def pressure_in_pa(da: xr.DataArray, dim: str = "plev") -> xr.DataArray:
    """Rename ``lev``/``level`` to ``plev`` and express it in Pa."""
    for cand in ("lev", "level"):
        if cand in da.dims and dim not in da.dims:
            da = da.rename({cand: dim})
    if dim not in da.coords:
        return da
    units = (da[dim].attrs.get("units") or "").strip().lower()
    vals = da[dim].values.astype(np.float64)
    if units in ("hpa", "mb", "mbar", "millibar", "millibars") or (not units and np.nanmax(vals) < 2000):
        vals = vals * 100.0
    elif units not in ("pa", "pascal", "pascals", ""):
        raise ValueError(f"{da.name}: unrecognized {dim} units {units!r}")
    da = da.assign_coords({dim: vals})
    da[dim].attrs["units"] = "Pa"
    return da


def select_plev(da: xr.DataArray, plev_pa, dim: str = "plev", tol: float = 1.0) -> xr.DataArray:
    """Select pressure levels (Pa) exactly, else nearest within ``tol`` Pa."""
    target = np.asarray(plev_pa, dtype=np.float64)
    out = da.sel({dim: target}, method="nearest", tolerance=tol)
    return out.assign_coords({dim: target})


def _daily_mean_model_calendar(field: xr.DataArray, per_day: int = 4) -> xr.DataArray:
    """00Z-00Z means of interval means labelled by interval start, in the model calendar.

    Averaging before any conversion to real dates keeps noleap days intact
    (converted, the last Feb 28 interval of 2012 would span 30 h). Days with
    fewer than ``per_day`` records are dropped.
    """
    day = field.time.dt.floor("D")
    grouped = field.assign_coords(day=("time", day.values)).groupby("day")
    daily = grouped.mean("time")
    count = grouped.count("time").isel({d: 0 for d in field.dims if d not in ("time",)})
    daily = daily.sel(day=count.day.values[count.values >= per_day]).rename(day="time")
    daily.attrs = dict(field.attrs, cell_methods="time: mean (00Z-00Z day)")
    return daily


# --------------------------------------------------------------------------
# Model readers
# --------------------------------------------------------------------------
def _files(directory: Path, variable: str) -> list[Path]:
    pattern = re.compile(rf"^{re.escape(variable)}_\d{{6}}_\d{{6}}\.nc$")
    files = sorted(p for p in Path(directory).glob(f"{variable}_*.nc") if pattern.match(p.name))
    if not files:
        raise FileNotFoundError(f"No {variable}_YYYYMM_YYYYMM.nc files under {directory}")
    return files


def open_model_daily(
    experiment: str,
    variable: str,
    *,
    data_root=DEFAULT_INPUT_DATA_ROOT,
    start: str = START_DATE,
    end: str = END_DATE,
    chunks=None,
) -> xr.DataArray:
    """Daily-mean model field (time, lat, lon) labelled by the day it covers.

    Read from ``ts/daily/1yr``; a variable held only as 6-hourly means
    (e.g. PSL) is averaged to complete 00Z-00Z days from ``ts/6hourly/1yr``.
    """
    try:
        files = _files(output_dir(experiment, "daily", data_root), variable)
        six_hourly = False
    except FileNotFoundError:
        files = _files(output_dir(experiment, "6hourly", data_root), variable)
        six_hourly = True
    ds = open_mfdataset_lazy(files, chunks=chunks or {"time": 92}, decode_times=_CFTIME,
                             data_vars="minimal", coords="minimal", compat="override")
    field = label_by_interval_start(ds, ds[variable])
    if six_hourly:
        field = _daily_mean_model_calendar(field)
    field = daily_label(standardize_latlon(field)).sel(time=slice(start, end))
    field.attrs = dict(field.attrs, experiment=experiment, sampling="daily_mean")
    field.name = variable
    return field


def open_model_monthly(
    experiment: str,
    variable: str,
    *,
    data_root=DEFAULT_INPUT_DATA_ROOT,
    start: str = START_DATE,
    end: str = END_DATE,
    cmip: bool = True,
    chunks=None,
) -> xr.DataArray:
    """Monthly model field labelled ``YYYY-MM-16`` (noleap), pressure in Pa.

    ``cmip=True`` reads the CMIP-named ``cmip_ts/monthly/<var>_Amon_*`` files
    (``ua``, ``ta``, ``zg``, ``psl`` on 19 levels); ``cmip=False`` reads the
    native names from ``ts/monthly/1yr`` and its ``<YYYY>_<YYYY>`` pressure-level
    subdirectories (``U``, ``T``, ``Z3``).
    """
    if cmip:
        directory = output_dir(experiment, "cmip_monthly", data_root)
        files = sorted(directory.glob(f"{variable}_A*mon_*.nc"))
        if not files:
            raise FileNotFoundError(f"No {variable}_Amon_*.nc under {directory}")
    else:
        directory = output_dir(experiment, "monthly", data_root)
        files = sorted(directory.glob(f"*_*/{variable}_*.nc")) or _files(directory, variable)
    ds = open_mfdataset_lazy(files, chunks=chunks or {"time": 12}, decode_times=_CFTIME,
                             data_vars="minimal", coords="minimal", compat="override")
    field = label_by_interval_start(ds, ds[variable])
    field = pressure_in_pa(standardize_latlon(monthly_label(field)))
    field = field.sel(time=slice(cftime.DatetimeNoLeap(*map(int, start.split("-")[:2]), 1),
                                 cftime.DatetimeNoLeap(*map(int, end.split("-")[:2]), 28)))
    field.attrs = dict(field.attrs, experiment=experiment)
    field.name = variable
    return field
