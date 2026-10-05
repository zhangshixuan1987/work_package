"""Reusable data and statistical helpers for MOV pathway notebooks."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import xarray as xr


def normalize_coords(obj):
    rename = {}
    for old, new in (("latitude", "lat"), ("longitude", "lon"), ("Latitude", "lat"), ("Longitude", "lon"), ("plev", "lev")):
        if old in obj.coords or old in obj.dims:
            rename[old] = new
    if rename:
        obj = obj.rename(rename)
    if "lon" in obj.coords:
        obj = obj.assign_coords(lon=((obj.lon + 180) % 360) - 180).sortby("lon")
    return obj


def select_level(da, level_hpa):
    names = [name for name in ("lev", "level", "plev") if name in da.coords or name in da.dims]
    if not names:
        return da
    name = names[0]
    target = level_hpa * 100.0 if float(da[name].max()) > 2000 else level_hpa
    return da.sel({name: target}, method="nearest")


def case_monthly_dirs(case, *, datasets, data_root: Path, search_subdirs: Sequence[Path]):
    root = Path(data_root) / datasets[case]["name"]
    return [root / Path(subdir) for subdir in search_subdirs]


def find_case_files(case, var_names, *, datasets, data_root, search_subdirs):
    directories = case_monthly_dirs(case, datasets=datasets, data_root=data_root, search_subdirs=search_subdirs)
    files = [path for base in directories if base.exists() for var in var_names for path in sorted(base.glob(f"{var}_*.nc"))]
    if not files:
        raise FileNotFoundError(f"No files found for {case} with variables {var_names}. Searched: {', '.join(map(str, directories))}")
    return files


def open_case_var(case, var_names, *, datasets, data_root, search_subdirs, period, level_hpa=None):
    files = find_case_files(case, var_names, datasets=datasets, data_root=data_root, search_subdirs=search_subdirs)
    ds = normalize_coords(xr.open_mfdataset(files, combine="by_coords", use_cftime=True))
    var = next((name for name in var_names if name in ds.data_vars), list(ds.data_vars)[0])
    da = normalize_coords(ds[var])
    if level_hpa is not None:
        da = select_level(da, level_hpa)
    return da.sel(time=slice(f"{period[0]}-01-01", f"{period[1]}-12-31"))


def standardize_units(da, kind):
    units = str(da.attrs.get("units", "")).lower()
    if kind == "precip" and (("kg" in units and "s-1" in units) or units in {"m/s", "m s-1"}):
        da = da * 86400.0
        da.attrs["units"] = "mm day-1"
    return da


def seasonal_year_mean(da, *, period, season="DJF"):
    da = da.sortby("time")
    if season.upper() == "DJF":
        month, year = da.time.dt.month, da.time.dt.year
        keep = month.isin([12, 1, 2])
        sub = da.where(keep, drop=True)
        winter_year = xr.where(month == 12, year + 1, year).where(keep, drop=True)
        return sub.assign_coords(season_year=("time", winter_year.data)).groupby("season_year").mean("time", skipna=True).rename({"season_year": "year"}).sel(year=slice(period[0] + 1, period[1]))
    sub = da.where(da.time.dt.season == season.upper(), drop=True)
    return sub.groupby("time.year").mean("time", skipna=True).sel(year=slice(*period))


def subset_box(da, bounds):
    da = normalize_coords(da)
    lat0, lat1 = bounds["lat"]
    lat_slice = slice(lat0, lat1) if da.lat[0] < da.lat[-1] else slice(lat1, lat0)
    return da.sel(lat=lat_slice, lon=slice(*bounds["lon"]))


def area_mean(da, bounds=None):
    sub = subset_box(da, bounds) if bounds else da
    return sub.weighted(np.cos(np.deg2rad(sub.lat))).mean(("lat", "lon"), skipna=True)


def regress_map(y, x, stats_module=None):
    y, x = xr.align(y, x, join="inner")
    x = (x - x.mean("year")) / x.std("year")
    slope = ((y - y.mean("year")) * x).mean("year")
    corr = xr.corr(y, x, dim="year")
    if stats_module is None:
        return slope, corr, xr.full_like(slope, np.nan)
    n = y.year.size
    tval = corr * np.sqrt((n - 2) / xr.where(1 - corr**2 <= 0, np.nan, 1 - corr**2))
    pval = xr.apply_ufunc(lambda value: 2 * stats_module.t.sf(np.abs(value), n - 2), tval)
    return slope, corr, pval


def highpass_transient(da, window=8):
    return da - da.rolling(time=window, center=True, min_periods=window // 2).mean()
