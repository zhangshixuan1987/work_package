"""Reference fields for the 2012-2013 AMIP hindcasts (``hind_*`` workflows).

Observations are read from one root, ``OBSERVATION_ROOT`` (a notebook
parameter), by default the reference-data archive
``/global/cfs/cdirs/e3sm/zhan391/data``:

* ``ERA5/1x1/6hourly/ERA5.6hourly.en00.<VAR>.<YYYY>01-<YYYY>12.nc`` (1979-2019;
  U850, U200, PSL, ...), averaged to 00Z-00Z days;
* ``obs_for_e3sm_diags/time-series/<product>/<var>_<YYYYMM>_<YYYYMM>.nc``
  (the E3SM Diags time series: NOAA OLR and IMERG daily; ERA5, ERA-Interim and
  NOAA-20C monthly).

``REFERENCES`` lists, per product and variable, candidate paths relative to
the root in search order; the first that exists is read (a warning reports
when it does not cover the requested dates). The later candidates follow the
layout of the per-machine copies (``reference`` on Perlmutter,
``acme_init/Observations`` on Compy), so the same code runs there. Daily fields
are labelled by the 00Z start of their day (``datetime64``), monthly fields
``YYYY-MM-16`` (``cftime.DatetimeNoLeap``), as in ``util.hind_model_data``.
"""
from __future__ import annotations

import re
from pathlib import Path

import cftime
import numpy as np
import xarray as xr

from util.hind_model_data import daily_label, monthly_label, pressure_in_pa, standardize_latlon

DEFAULT_OBSERVATION_ROOT = Path("/global/cfs/cdirs/e3sm/zhan391/data")
_TS = "obs_for_e3sm_diags/time-series"

# (product, variable) -> [(glob relative to OBSERVATION_ROOT, variable in file), ...] in search order
REFERENCES = {
    ("NOAA_OLR", "FLUT"): [(f"{_TS}/NOAA-OLR_Daily/FLUT_*.nc", "FLUT"),
                           ("NOAA_IOLR/daily/FLUT.daily.*.nc", "FLUT")],
    ("IMERG", "PRECT"): [(f"{_TS}/IMERG_Daily/PRECT_*.nc", "PRECT"),
                         ("NASA_GPM/daily/PRECT.daily.*.nc", "PRECT"),
                         ("IMERG/daily/PRECT.daily.*.nc", "PRECT")],
    ("NOAA20C", "psl"): [(f"{_TS}/NOAA-20C/psl_*.nc", "psl")],
}
# ERA5 daily from the 6-hourly 1-degree archive (one source for all variables)
for _var in ("U850", "U200", "V850", "V200", "PSL"):
    REFERENCES[("ERA5_DAILY", _var)] = [(f"ERA5/1x1/6hourly/ERA5.6hourly.en00.{_var}.*.nc", _var),
                                        (f"ERA5/daily/{_var}_*.nc", _var),
                                        (f"ERA5/6hourly/ERA5.6hourly.en00.{_var}.*.nc", _var)]
for _var in ("ua", "ta", "zg", "va", "psl"):
    REFERENCES[("ERA5", _var)] = [(f"{_TS}/ERA5/{_var}_*.nc", _var)]
    REFERENCES[("ERAI", _var)] = [(f"{_TS}/ERA-Interim/{_var}_*.nc", _var)]

# Monthly reanalyses usable as QBO/SSW references
REANALYSES = ("ERA5", "ERAI")

_RANGE = re.compile(r"(\d{4})\d{2}[-_](\d{4})\d{2}")
_YEAR = re.compile(r"\.(\d{4})\.nc$")


def _overlapping(files: list[Path], start: str | None, end: str | None) -> list[Path]:
    """Files whose YYYYMM-YYYYMM (or .YYYY.) name range overlaps [start - 1 yr, end]."""
    if start is None or end is None:
        return files
    y0, y1 = int(start[:4]) - 1, int(end[:4])
    keep = []
    for f in files:
        m = _RANGE.search(f.name)
        a, b = (int(m.group(1)), int(m.group(2))) if m else ((int(_YEAR.search(f.name).group(1)),) * 2
                                                            if _YEAR.search(f.name) else (y0, y1))
        if b >= y0 and a <= y1:
            keep.append(f)
    return keep or files


def _open(files, chunks):
    """Open lazily with dask chunks set at open time.

    The archive files are unchunked (contiguous) NetCDF, so a dimension left
    unchunked becomes one task that reads the whole multi-GB variable; chunking
    at open lets each task read only the slices it needs.
    """
    return xr.open_mfdataset(files, combine="by_coords", parallel=len(files) > 1, chunks=chunks,
                             data_vars="minimal", coords="minimal", compat="override")


def reference_files(product: str, variable: str, *, observation_root=DEFAULT_OBSERVATION_ROOT,
                    start: str | None = None, end: str | None = None) -> tuple[list[Path], str]:
    """Files (those overlapping start..end) and in-file variable name of the first candidate that exists."""
    key = (product, variable)
    if key not in REFERENCES:
        raise KeyError(f"No reference registered for {key}; known: {sorted(REFERENCES)}")
    root = Path(observation_root)
    tried = []
    for pattern, name in REFERENCES[key]:
        files = sorted(root.glob(pattern))
        if files:
            return _overlapping(files, start, end), name
        tried.append(str(root / pattern))
    raise FileNotFoundError(f"No files for {product} {variable}; tried:\n  " + "\n  ".join(tried))


def open_reference_daily(product: str, variable: str, *, start: str, end: str,
                         observation_root=DEFAULT_OBSERVATION_ROOT, lat_range=None,
                         chunks=None) -> xr.DataArray:
    """Daily reference (time, lat, lon) labelled by the 00Z start of each day.

    Sub-daily files (ERA5 6-hourly) are averaged to complete 00Z-00Z days.
    ``lat_range`` subsets before anything is read.
    """
    files, name = reference_files(product, variable, observation_root=observation_root, start=start, end=end)
    ds = _open(files, chunks or {"time": 124})
    field = standardize_latlon(ds[name])
    if lat_range is not None:
        field = field.sel(lat=slice(*lat_range))
    field = field.sel(time=slice(np.datetime64(start), np.datetime64(end) + np.timedelta64(1, "D")))
    steps = np.diff(field.time.values[:8].astype("datetime64[h]").astype(np.int64))
    if steps.size and np.median(steps) < 24:
        per_day = int(round(24 / np.median(steps)))
        # Complete days from the time axis alone (no data are read), then a blockwise mean of each
        # day's ``per_day`` consecutive steps: without flox a dask groupby makes one task per day.
        day = field.time.dt.floor("D").values
        days, steps_per_day = np.unique(day, return_counts=True)
        complete = days[steps_per_day == per_day]
        field = field.isel(time=np.flatnonzero(np.isin(day, complete)))
        field = (field.chunk(time=per_day * 31).coarsen(time=per_day).mean()
                 .assign_coords(time=complete.astype("datetime64[ns]")))
    field = daily_label(field).sel(time=slice(start, end))
    if field.sizes["time"] == 0 or field.time.values[0] > np.datetime64(start) \
            or field.time.values[-1] < np.datetime64(end):
        span = f"{str(field.time.values[0])[:10]}..{str(field.time.values[-1])[:10]}" if field.sizes["time"] else "none"
        print(f"[WARN] {product} {variable}: {files[0].parent} covers {span} of {start}..{end}")
    field.attrs = dict(field.attrs, reference=product, sampling="daily_mean")
    field.name = variable
    return field


def open_reference_monthly(product: str, variable: str, *, start: str, end: str,
                           observation_root=DEFAULT_OBSERVATION_ROOT, chunks=None) -> xr.DataArray:
    """Monthly reference labelled ``YYYY-MM-16`` (noleap), pressure in Pa.

    Opened lazily; select levels and latitudes before computing (the ERA5
    pressure-level files are 0.25 degree, ~75 GB each).
    """
    files, name = reference_files(product, variable, observation_root=observation_root, start=start, end=end)
    ds = _open(files, chunks or {"time": 1, "plev": 1})
    field = ds[name]
    s = np.datetime64(start[:7]) if np.issubdtype(field.time.dtype, np.datetime64) else None
    if s is not None:
        field = field.sel(time=slice(s, np.datetime64(end[:7]) + np.timedelta64(1, "M")))
    field = pressure_in_pa(standardize_latlon(monthly_label(field)))
    field = field.sel(time=slice(cftime.DatetimeNoLeap(*map(int, start.split("-")[:2]), 1),
                                 cftime.DatetimeNoLeap(*map(int, end.split("-")[:2]), 28)))
    field.attrs = dict(field.attrs, reference=product)
    field.name = variable
    return field
