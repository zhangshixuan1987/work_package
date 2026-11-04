"""Daily member fields and matching references for the F20TR initialization hindcasts.

Layout and experiments: ``configs/fcst_init_experiment_config.py``. Every field is
indexed by forecast lead day ``lead`` = 1..``max_lead``; lead day 1 is the start date
(00Z-24Z of the initial day), and ``time`` (lead) is the verification day. Fields
of one configuration and season stack to (init, member, lead, lat, lon).

Sampling, as in ``util.forecast_window_scorecard`` (``VARIABLES[...]["sampling"]``):

* ``daily_mean``: 00Z-24Z means. From the 6-hourly stream where it exists
  (``util.model_daily.daily_from_6hourly``: interval means averaged over the day,
  instantaneous PSL as the mean of the 00/06/12/18Z states), else the daily means
  of the daily stream (FLUT, PRECT; ELM fields), which are stamped at the day end
  and relabelled by their start (``util.time_labels.label_by_interval_start``).
  References: ERA5 daily means of the 6-hourly 1-degree archive, IMERG, NOAA OLR,
  ESA CCI.
* ``00z``: the pressure-level fields are 00Z states in the daily stream; lead day
  d is represented by the state at 00Z ending it (24 h x d after the start), and
  verified with ERA5 at the same instant.

Units follow the scorecard: Z in m, Q in g kg-1, PRECT in mm day-1, PSL in hPa,
soil moisture of the top 5 cm in m3 m-3.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from configs.fcst_init_experiment_config import (DEFAULT_INPUT_DATA_ROOT, MEMBERS, init_dates, ts_file)

DEFAULT_OBSERVATION_ROOT = Path("/global/cfs/cdirs/e3sm/zhan391/data")
DEFAULT_LAND_SUPPORT_DIR = DEFAULT_OBSERVATION_ROOT / "lnd_sea_mask"
MAX_LEAD = 56          # 8 weeks; the runs cover about 62 days
_G = 9.80665
_DAY = np.timedelta64(1, "D")

_ERA5 = "ERA5/1x1/6hourly/ERA5.6hourly.en00.{name}.{year}01-{year}12.nc"
_IMERG = "obs_for_e3sm_diags/time-series/IMERG_Daily/PRECT_*.nc"
_OLR = "obs_for_e3sm_diags/time-series/NOAA-OLR_Daily/FLUT_*.nc"
_ESA = "ESA_CCI_H2OSOI/1x1/daily/H2OSOI.daily.{year}.nc"


@dataclass(frozen=True)
class VarSpec:
    component: str      # atm | lnd (lnd variables are scored over land)
    model: str          # model variable
    stream: str         # 6hourly | daily (post-processed stream holding it)
    sampling: str       # daily_mean | 00z
    reference: str      # ERA5 | IMERG | NOAA_OLR | ESA_CCI
    ref_name: str       # variable in the reference files
    units: str


VARIABLES: dict[str, VarSpec] = {
    "PSL": VarSpec("atm", "PSL", "6hourly", "daily_mean", "ERA5", "PSL", "hPa"),
    "TREFHT": VarSpec("atm", "TREFHT", "6hourly", "daily_mean", "ERA5", "TREFHT", "K"),
    "TS": VarSpec("atm", "TS", "6hourly", "daily_mean", "ERA5", "TS", "K"),
    "PRECT": VarSpec("atm", "PRECT", "daily", "daily_mean", "IMERG", "PRECT", "mm day-1"),
    "FLUT": VarSpec("atm", "FLUT", "daily", "daily_mean", "NOAA_OLR", "FLUT", "W m-2"),
    **{v: VarSpec("atm", v, "6hourly", "daily_mean", "ERA5", v, "W m-2")
       for v in ("LHFLX", "SHFLX", "FSNS", "FSDS", "FLDS")},
    **{f"{n}{p}": VarSpec("atm", f"{n}{p}", "daily", "00z", "ERA5", f"{n}{p}", u)
       for n, u in (("U", "m s-1"), ("V", "m s-1"), ("T", "K"), ("Q", "g kg-1"), ("Z", "m"))
       for p in (200, 500, 700, 850, 925)},
    "TSA": VarSpec("lnd", "TSA", "daily", "daily_mean", "ERA5", "TREFHT", "K"),
    "SOILM_5CM": VarSpec("lnd", "H2OSOI", "daily", "daily_mean", "ESA_CCI", "H2OSOI", "m3 m-3"),
    "SOILWATER_10CM": VarSpec("lnd", "SOILWATER_10CM", "daily", "daily_mean", "", "", "kg m-2"),
    "TWS": VarSpec("lnd", "TWS", "daily", "daily_mean", "", "", "mm"),
    "H2OSNO": VarSpec("lnd", "H2OSNO", "daily", "daily_mean", "", "", "mm"),
}


def spec(variable: str) -> VarSpec:
    if variable not in VARIABLES:
        raise KeyError(f"Unknown variable {variable!r}; choose from {list(VARIABLES)}")
    return VARIABLES[variable]


def _convert(field: xr.DataArray, variable: str) -> xr.DataArray:
    """Model and reference fields in catalog units."""
    units = str(field.attrs.get("units", ""))
    name = variable.rstrip("0123456789")
    if variable == "PSL" and units == "Pa":
        field = field * 0.01
    elif name == "Z" and units.startswith("m**2"):
        field = field / _G
    elif name == "Q" and variable != "Q" and units.lower().startswith("kg"):
        field = field * 1000.0
    elif variable == "PRECT" and units in {"m/s", "m s-1"}:
        field = field * 86400 * 1000
    field = field.astype("float32")
    field.attrs = dict(units=spec(variable).units)
    return field


def _normalize(field: xr.DataArray) -> xr.DataArray:
    rename = {k: v for k, v in (("latitude", "lat"), ("longitude", "lon")) if k in field.dims}
    field = field.rename(rename) if rename else field
    extra = [d for d in field.dims if d not in ("time", "lat", "lon", "levgrnd") and field.sizes[d] == 1]
    if extra:
        field = field.squeeze(extra, drop=True)
    if float(field.lon.min()) < 0:
        field = field.assign_coords(lon=field.lon % 360)
    return field.sortby("lat").sortby("lon")


def lead_days(init, max_lead: int = MAX_LEAD) -> pd.DatetimeIndex:
    """Verification day of each lead day 1..max_lead (lead 1 = the start date)."""
    return pd.date_range(pd.Timestamp(np.datetime64(init, "D")), periods=max_lead, freq="D")


def _on_leads(field: xr.DataArray, init, max_lead: int, sampling: str) -> xr.DataArray:
    """(lead, ...) field of the verification days; days not present become NaN."""
    days = lead_days(init, max_lead)
    if sampling == "00z":
        instants = days + pd.Timedelta(days=1)  # 00Z state ending each day
        field = field.reindex(time=instants.values)
    else:
        field = field.assign_coords(time=field.time.dt.floor("D")).reindex(time=days.values)
    return (field.assign_coords(time=days.values).rename(time="lead")
            .assign_coords(lead=np.arange(1, max_lead + 1), time=("lead", days.values)))


_dz_cache: dict[str, xr.DataArray] = {}


def top_soil_weights(land_support_dir=DEFAULT_LAND_SUPPORT_DIR, depth_cm: float = 5.0) -> xr.DataArray:
    """Thickness weights of the ELM layers within the top ``depth_cm`` (ESA CCI sensing depth)."""
    from util.lnd_sm_daily_skill import top_layer_weights

    key = f"{land_support_dir}:{depth_cm}"
    if key not in _dz_cache:
        with xr.open_dataset(Path(land_support_dir) / "dzsoi_elm.nc") as ds:
            w = top_layer_weights(ds["DZSOI"].load(), depth_cm)
        _dz_cache[key] = w.isel(levgrnd=slice(0, int(w.notnull().any(("lat", "lon")).sum())))
    return _dz_cache[key]


def open_member_daily(config: str, init, member: str, variable: str, *, max_lead: int = MAX_LEAD,
                      data_root=DEFAULT_INPUT_DATA_ROOT, land_support_dir=DEFAULT_LAND_SUPPORT_DIR,
                      pair_mean: bool = False) -> xr.DataArray:
    """(lead, lat, lon) daily field of one member in catalog units (lazy).

    ``pair_mean``: for a ``00z`` variable, estimate each day's mean as the mean of the 00Z
    states at its start and end (lead day 1 uses the initial state), for workflows that
    combine fields as daily values (RMM; as ``util.model_daily.open_model_daily_estimate``).
    """
    from util.model_daily import daily_from_6hourly
    from util.time_labels import label_by_interval_start

    s = spec(variable)
    sampling = s.sampling
    path = ts_file(config, init, member, s.model, s.component, s.stream, data_root)
    ds = xr.open_dataset(path, chunks={"time": -1})
    if s.stream == "6hourly":
        field = daily_from_6hourly(ds, s.model)
    elif s.sampling == "00z":
        field = ds[s.model]
        field = field.sel(time=field.time.dt.hour == 0)
        if pair_mean:
            states = field.assign_coords(time=field.time.dt.floor("D"))
            field = 0.5 * (states + states.shift(time=-1))
            field.attrs = dict(ds[s.model].attrs)
            sampling = "daily_mean"
    else:
        field = label_by_interval_start(ds, ds[s.model])
    field.attrs["units"] = ds[s.model].attrs.get("units", "")
    field = _normalize(field)
    if variable == "SOILM_5CM":
        w = top_soil_weights(land_support_dir)
        field = field.isel(levgrnd=slice(0, w.sizes["levgrnd"]))
        units = field.attrs.get("units", "")
        field = (field * w.assign_coords(levgrnd=field.levgrnd).values[None]).sum("levgrnd", min_count=1)
        field.attrs["units"] = units
    out = _convert(_on_leads(field, init, max_lead, sampling), variable).rename(variable)
    out.attrs["sampling"] = "00z_pair_mean" if sampling != s.sampling else s.sampling
    return out


def open_ensemble(config: str, init, variable: str, *, members=None, max_lead: int = MAX_LEAD,
                  data_root=DEFAULT_INPUT_DATA_ROOT, land_support_dir=DEFAULT_LAND_SUPPORT_DIR,
                  pair_mean: bool = False) -> xr.DataArray:
    """(member, lead, lat, lon) of one configuration and start date."""
    members = MEMBERS[config] if members is None else members
    fields = [open_member_daily(config, init, m, variable, max_lead=max_lead, data_root=data_root,
                                land_support_dir=land_support_dir, pair_mean=pair_mean) for m in members]
    return xr.concat(fields, pd.Index(list(members), name="member"), coords="minimal", compat="override")


def open_campaign(config: str, season: str, variable: str, *, members=None, max_lead: int = MAX_LEAD,
                  data_root=DEFAULT_INPUT_DATA_ROOT, land_support_dir=DEFAULT_LAND_SUPPORT_DIR) -> xr.DataArray:
    """(init, member, lead, lat, lon) of every start date of ``season``; ``time`` is (init, lead)."""
    inits = init_dates([season])[season]
    fields = [open_ensemble(config, i, variable, members=members, max_lead=max_lead, data_root=data_root,
                            land_support_dir=land_support_dir).drop_vars("time") for i in inits]
    out = xr.concat(fields, pd.Index(np.array(inits, dtype="datetime64[ns]"), name="init"))
    return out.assign_coords(time=valid_time(out.init.values, max_lead))


def valid_time(inits, max_lead: int = MAX_LEAD) -> xr.DataArray:
    """(init, lead) verification day."""
    inits = np.asarray(inits, dtype="datetime64[ns]")
    t = inits[:, None] + (np.arange(max_lead) * _DAY).astype("timedelta64[ns]")[None]
    return xr.DataArray(t, dims=("init", "lead"),
                        coords={"init": inits, "lead": np.arange(1, max_lead + 1)})


# --------------------------------------------------------------------------
# References
# --------------------------------------------------------------------------
def _reference_series(variable: str, start, end, observation_root) -> xr.DataArray:
    """Daily reference (time labelled by its day; 00Z instants for ``00z``) covering start..end."""
    s = spec(variable)
    root = Path(observation_root)
    y0, y1 = pd.Timestamp(start).year, (pd.Timestamp(end) + pd.Timedelta(days=1)).year
    if s.reference == "ERA5":
        files = [root / _ERA5.format(name=s.ref_name, year=y) for y in range(y0, y1 + 1)]
    elif s.reference == "ESA_CCI":
        files = [root / _ESA.format(year=y) for y in range(y0, y1 + 1)]
    elif s.reference in ("IMERG", "NOAA_OLR"):
        files = sorted(root.glob(_IMERG if s.reference == "IMERG" else _OLR))
    else:
        raise ValueError(f"{variable} has no reference")
    missing = [str(f) for f in files if not f.is_file()]
    if missing or not files:
        raise FileNotFoundError(f"{variable}: missing {s.reference} files {missing or files}")
    ds = xr.open_mfdataset(files, combine="by_coords", chunks={"time": 124}, data_vars="minimal",
                           coords="minimal", compat="override")
    field = _normalize(ds[s.ref_name])
    field.attrs["units"] = ds[s.ref_name].attrs.get("units", "")
    t0, t1 = np.datetime64(pd.Timestamp(start)), np.datetime64(pd.Timestamp(end)) + 2 * _DAY
    field = field.sel(time=slice(t0, t1))
    if s.reference == "ERA5":
        if s.sampling == "00z":
            field = field.sel(time=field.time.dt.hour == 0)
        else:  # complete 00Z-24Z days of the 6-hourly analyses
            day = field.time.dt.floor("D").values
            days, n = np.unique(day, return_counts=True)
            keep = np.isin(day, days[n == 4])
            field = field.isel(time=np.flatnonzero(keep))
            field = field.coarsen(time=4).mean().assign_coords(time=days[n == 4])
            field.attrs["units"] = ds[s.ref_name].attrs.get("units", "")
    else:
        field = field.assign_coords(time=field.time.dt.floor("D"))
    return field


def open_reference(variable: str, inits, *, max_lead: int = MAX_LEAD,
                   observation_root=DEFAULT_OBSERVATION_ROOT) -> xr.DataArray:
    """(init, lead, lat, lon) reference on every start date's verification days (lazy)."""
    inits = np.asarray(inits, dtype="datetime64[D]")
    s = spec(variable)
    series = _reference_series(variable, inits.min(), inits.max() + max_lead * _DAY, observation_root)
    out = [_on_leads(series, i, max_lead, s.sampling).drop_vars("time") for i in inits]
    out = xr.concat(out, pd.Index(inits.astype("datetime64[ns]"), name="init"))
    return _convert(out.assign_coords(time=valid_time(inits, max_lead)), variable).rename(variable)


def reference_climatology(variable: str, inits, years, *, max_lead: int = MAX_LEAD, smooth_days: int = 31,
                          observation_root=DEFAULT_OBSERVATION_ROOT) -> xr.DataArray:
    """(init, lead, lat, lon) calendar-day climatology of the reference over ``years``.

    For each year the reference is read on the same calendar days as the verification
    days (padded by ``smooth_days // 2``), averaged over the years, smoothed with a
    centred ``smooth_days`` running mean, and put back on each start date's leads.
    Feb 29 is skipped in non-leap years (its value is the Feb 28/Mar 1 neighbours' mean).
    """
    inits = np.asarray(inits, dtype="datetime64[D]")
    s = spec(variable)
    half = smooth_days // 2
    first, last = pd.Timestamp(inits.min()) - pd.Timedelta(days=half), \
        pd.Timestamp(inits.max() + max_lead * _DAY) + pd.Timedelta(days=half)
    ref_year = first.year
    target = pd.date_range(first, last, freq="D")
    source = "daily"
    if s.reference == "ERA5" and not _era5_6hourly_covers(s.ref_name, years, first, last, observation_root):
        clim = _era5_monthly_climatology(variable, years, target, observation_root)
        source = "ERA5 monthly"
    else:
        stack = []
        for y in years:
            shift = y - ref_year
            a, b = first + pd.DateOffset(years=shift), last + pd.DateOffset(years=shift)
            series = _reference_series(variable, a, b, observation_root)
            days = pd.DatetimeIndex([d + pd.DateOffset(years=shift) for d in target])
            if s.sampling == "00z":
                series = series.assign_coords(time=series.time.dt.floor("D") - _DAY)
            series = series.reindex(time=days.values).assign_coords(time=target.values)
            stack.append(series.load())
        clim = xr.concat(stack, "year").mean("year")
        if smooth_days > 1:
            clim = clim.rolling(time=smooth_days, center=True, min_periods=half + 1).mean()
        clim.attrs["units"] = stack[0].attrs.get("units", "")
    shifted = clim.assign_coords(time=clim.time + _DAY) if s.sampling == "00z" else clim
    out = [_on_leads(shifted, i, max_lead, s.sampling).drop_vars("time") for i in inits]
    out = xr.concat(out, pd.Index(inits.astype("datetime64[ns]"), name="init"))
    out = _convert(out.assign_coords(time=valid_time(inits, max_lead)), variable).rename(variable)
    out.attrs.update(climatology_years=f"{min(years)}-{max(years)}", smooth_days=smooth_days,
                     climatology_source=source)
    return out


_ERA5_MONTHLY = "ERA5/1x1/monthly/ERA5.en00.{name}.197901-201812.nc"


def _era5_6hourly_covers(name, years, first, last, observation_root) -> bool:
    root = Path(observation_root)
    span = range(first.year - first.year, last.year - first.year + 1)
    return all((root / _ERA5.format(name=name, year=y + k)).is_file() for y in years for k in span)


def _era5_monthly_climatology(variable: str, years, target: pd.DatetimeIndex, observation_root) -> xr.DataArray:
    """Daily climatology on ``target`` interpolated linearly between the mid-month means of
    the ERA5 monthly series (for TREFHT, TS and fluxes, whose 6-hourly files cover only
    2011-2012; the same choice as ``util.forecast_window_scorecard``)."""
    s = spec(variable)
    path = Path(observation_root) / _ERA5_MONTHLY.format(name=s.ref_name)
    if not path.is_file():
        raise FileNotFoundError(f"{variable}: no ERA5 6-hourly files for {min(years)}-{max(years)} and no {path}")
    with xr.open_dataset(path) as ds:
        field = _normalize(ds[s.ref_name])
        units = ds[s.ref_name].attrs.get("units", "")
        field = field.sel(time=field.time.dt.year.isin(list(years))).load()
    monthly = field.groupby("time.month").mean("time")
    # mid-month anchors spanning the target period (one month of padding on each side)
    anchors = pd.date_range(target[0] - pd.offsets.MonthBegin(2), target[-1] + pd.offsets.MonthBegin(2), freq="MS") \
        + pd.Timedelta(days=14)
    values = monthly.sel(month=xr.DataArray(anchors.month, dims="time")).assign_coords(time=anchors.values)
    clim = values.interp(time=target.values)
    clim.attrs["units"] = units
    return clim


# --------------------------------------------------------------------------
# Region masks
# --------------------------------------------------------------------------
def region_masks(regions, lat, lon, *, land_support_dir=DEFAULT_LAND_SUPPORT_DIR) -> xr.DataArray:
    """(region, lat, lon) boolean masks of ``configs.s2s_config``-style regions (``mask``: land/ocean)."""
    with xr.open_dataset(Path(land_support_dir) / "LANDFRAC_E3SM_180x360.nc") as ds:
        landfrac = _normalize(ds["LANDFRAC"].squeeze(drop=True)).load()
    landfrac = landfrac.reindex(lat=lat, lon=lon, method="nearest", tolerance=0.01)
    out = []
    for r in regions:
        lo0, lo1 = r["lon"]
        m = (lat >= r["lat"][0]) & (lat <= r["lat"][1])
        lon_ok = ((lon >= lo0) & (lon <= lo1)) if lo0 <= lo1 else ((lon >= lo0) | (lon <= lo1))
        mask = xr.DataArray(m[:, None] & lon_ok[None, :], dims=("lat", "lon"), coords={"lat": lat, "lon": lon})
        if r.get("mask") == "land":
            mask = mask & (landfrac > 0.5)
        elif r.get("mask") == "ocean":
            mask = mask & (landfrac < 0.5)
        out.append(mask)
    return xr.concat(out, pd.Index([r["name"] for r in regions], name="region"))
