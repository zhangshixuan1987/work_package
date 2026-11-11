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
  References: the ERA5 1-degree daily archive (``ERA5/1x1/daily``, daily means of the
  6-hourly archive: exact for fluxes, trapezoid rule over 00Z..00Z for states), GPCP V3.3 daily,
  NOAA OLR, ESA CCI.
* ``00z``: the pressure-level fields are 00Z states in the daily stream; lead day
  d is represented by the state at 00Z ending it (24 h x d after the start), and
  verified with ERA5 at the same instant (``ERA5/1x1/6hourly``).

Units follow the scorecard: Z in m, Q in g kg-1, PRECT in mm day-1, PSL in hPa, TMQ in
kg m-2, soil moisture of the top 5 cm in m3 m-3.

Climatologies (``reference_climatology``, ``model_climatology``) are built the same way
for every variable and for both sides: calendar-month means over the climatology years,
interpolated linearly between mid-month values to the (init, lead) verification days.
The reference uses the monthly means of its daily series (ERA5: the daily archive, so
the climatology comes from the same product as the verifying days), the model the
monthly output of the free-running E3SMv3 AMIP runs (``month_of_records`` reads the
month each monthly record averages). Both are all-hours means,
also for the ``00z`` fields, so each 00Z anomaly carries its own system's 00Z departure
from the daily mean (ERA5: RMS about 6 m Z500, 0.5 K T850, 0.5 m s-1 U850/U200).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from configs.fcst_init_experiment_config import (DEFAULT_INPUT_DATA_ROOT, init_dates, members as config_members,
                                                  season_of, ts_file)

# Locations are required arguments set by the notebooks (they depend on the machine):
# observation_root (the reference-data archive), land_support_dir (its lnd_sea_mask folder,
# land fraction and ELM layer thicknesses) and amip_root (its E3SMv3 folder). Paths below
# are relative to observation_root.
MAX_LEAD = 56          # 8 weeks; the runs cover about 62 days
_G = 9.80665
_DAY = np.timedelta64(1, "D")

_ERA5_6H = "ERA5/1x1/6hourly/ERA5.6hourly.en00.{name}.{year}01-{year}12.nc"   # 00Z states
_ERA5_DAILY = "ERA5/1x1/daily/ERA5.daily.en00.{name}.{year}01-{year}12.nc"     # daily means
_IMERG = "NASA_GPM/1x1/daily/PRECT.daily.*.nc"   # IMERG, 60S-60N, 2001-2020, m s-1 (= e3sm_diags IMERG_Daily)
# GPCP V3.3 daily, global, 1998-2025, mm/day: built 2026-11-05 by GPCPv3.3/gpcp_v33_daily_1x1.py
# (2 x 2 area-weighted box mean of the 0.5-degree GES DISC files). Not obs_for_e3sm_diags GPCPDAY_3.3,
# whose 1-degree file has its latitudes reversed, nor GPCPv3.3/1x1/daily/GPCP_daily_* (V3.1, 2007-2017).
_GPCP = "GPCPv3.3/1x1/daily/PRECT.daily.*.nc"
_OLR = "NOAA_CDR/1x1/daily/FLUT.daily.*.nc"     # NOAA OLR CDR v02r00, 1979-2025, W m-2
_ESA = "ESA_CCI_H2OSOI/1x1/daily/H2OSOI.daily.{year}.nc"


@dataclass(frozen=True)
class VarSpec:
    component: str      # atm | lnd (lnd variables are scored over land)
    model: str          # model variable
    stream: str         # 6hourly | daily (post-processed stream holding it)
    sampling: str       # daily_mean | 00z
    reference: str      # ERA5 | GPCP | IMERG | NOAA_OLR | ESA_CCI
    ref_name: str       # variable in the reference files
    units: str


VARIABLES: dict[str, VarSpec] = {
    "PSL": VarSpec("atm", "PSL", "6hourly", "daily_mean", "ERA5", "PSL", "hPa"),
    "TREFHT": VarSpec("atm", "TREFHT", "6hourly", "daily_mean", "ERA5", "TREFHT", "K"),
    "TS": VarSpec("atm", "TS", "6hourly", "daily_mean", "ERA5", "TS", "K"),
    "PRECT": VarSpec("atm", "PRECT", "daily", "daily_mean", "GPCP", "PRECT", "mm day-1"),   # GPCP V3.3 (global)
    "FLUT": VarSpec("atm", "FLUT", "daily", "daily_mean", "NOAA_OLR", "FLUT", "W m-2"),
    "TMQ": VarSpec("atm", "TMQ", "6hourly", "daily_mean", "ERA5", "TMQ", "kg m-2"),
    # vertically integrated water vapour flux (eastward, northward): ERA5 VIWVE/VIWVN, same units and sign
    "TUQ": VarSpec("atm", "TUQ", "6hourly", "daily_mean", "ERA5", "VIWVE", "kg m-1 s-1"),
    "TVQ": VarSpec("atm", "TVQ", "6hourly", "daily_mean", "ERA5", "VIWVN", "kg m-1 s-1"),
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
    elif variable == "PRECT" and units in {"kg m-2 s-1", "kg/m2/s"}:
        field = field * 86400
    elif variable == "TMQ" and units.startswith("g m"):  # ERA5 monthly product
        field = field * 0.001
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


def top_soil_weights(land_support_dir, depth_cm: float = 5.0) -> xr.DataArray:
    """Thickness weights of the ELM layers within the top ``depth_cm`` (ESA CCI sensing depth)."""
    from util.lnd_sm_daily_skill import top_layer_weights

    key = f"{land_support_dir}:{depth_cm}"
    if key not in _dz_cache:
        with xr.open_dataset(Path(land_support_dir) / "dzsoi_elm.nc") as ds:
            w = top_layer_weights(ds["DZSOI"].load(), depth_cm)
        _dz_cache[key] = w.isel(levgrnd=slice(0, int(w.notnull().any(("lat", "lon")).sum())))
    return _dz_cache[key]


def _subset_box(ds: xr.Dataset, box: dict | None) -> xr.Dataset:
    """``ds`` cut to the lat/lon ``box`` (``lat``, ``lon`` ranges in degrees north / east 0-360; a
    ``lon`` range with lon0 > lon1 wraps through 0) by index, whatever the file's coordinate names,
    longitude convention or ordering; ``box=None`` returns ``ds`` unchanged."""
    if box is None:
        return ds
    lat = "lat" if "lat" in ds.dims else "latitude"
    lon = "lon" if "lon" in ds.dims else "longitude"
    la, lo = ds[lat].values, ds[lon].values % 360
    lo0, lo1 = box["lon"]
    lon_ok = ((lo >= lo0) & (lo <= lo1)) if lo0 <= lo1 else ((lo >= lo0) | (lo <= lo1))
    lat_ok = (la >= box["lat"][0]) & (la <= box["lat"][1])
    return ds.isel({lat: np.flatnonzero(lat_ok), lon: np.flatnonzero(lon_ok)})


def open_member_daily(config: str, init, member: str, variable: str, *, max_lead: int = MAX_LEAD,
                      data_root=DEFAULT_INPUT_DATA_ROOT, land_support_dir,
                      pair_mean: bool = False, box: dict | None = None) -> xr.DataArray:
    """(lead, lat, lon) daily field of one member in catalog units (lazy).

    ``box``: cut the file to this lat/lon box (as ``heatwave_case.to_box``) right after opening, so
    the 6-hourly-to-daily averaging and the reads touch only the box (default: the whole grid).

    ``pair_mean``: for a ``00z`` variable, estimate each day's mean as the mean of the 00Z
    states at its start and end (lead day 1 uses the initial state), for workflows that
    combine fields as daily values (RMM; as ``util.model_daily.open_model_daily_estimate``).
    """
    from util.model_daily import daily_from_6hourly
    from util.time_labels import label_by_interval_start

    s = spec(variable)
    sampling = s.sampling
    path = ts_file(config, init, member, s.model, s.component, s.stream, data_root)
    ds = _subset_box(xr.open_dataset(path, chunks={"time": -1}), box)
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
        if box is not None:
            w = _normalize(w).reindex(lat=field.lat, lon=field.lon, method="nearest", tolerance=0.01)
        field = field.isel(levgrnd=slice(0, w.sizes["levgrnd"]))
        units = field.attrs.get("units", "")
        field = (field * w.assign_coords(levgrnd=field.levgrnd).values[None]).sum("levgrnd", min_count=1)
        field.attrs["units"] = units
    out = _convert(_on_leads(field, init, max_lead, sampling), variable).rename(variable)
    out.attrs["sampling"] = "00z_pair_mean" if sampling != s.sampling else s.sampling
    return out


def open_ensemble(config: str, init, variable: str, *, members=None, max_lead: int = MAX_LEAD,
                  data_root=DEFAULT_INPUT_DATA_ROOT, land_support_dir,
                  pair_mean: bool = False, box: dict | None = None) -> xr.DataArray:
    """(member, lead, lat, lon) of one configuration and start date."""
    members = config_members(config, season_of(init)) if members is None else members
    fields = [open_member_daily(config, init, m, variable, max_lead=max_lead, data_root=data_root,
                                land_support_dir=land_support_dir, pair_mean=pair_mean, box=box) for m in members]
    return xr.concat(fields, pd.Index(list(members), name="member"), coords="minimal", compat="override")


def open_campaign(config: str, season: str, variable: str, *, members=None, max_lead: int = MAX_LEAD,
                  data_root=DEFAULT_INPUT_DATA_ROOT, land_support_dir, box: dict | None = None) -> xr.DataArray:
    """(init, member, lead, lat, lon) of every start date of ``season``; ``time`` is (init, lead)."""
    inits = init_dates([season])[season]
    fields = [open_ensemble(config, i, variable, members=members, max_lead=max_lead, data_root=data_root,
                            land_support_dir=land_support_dir, box=box).drop_vars("time") for i in inits]
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
def _reference_series(variable: str, start, end, observation_root, sampling: str | None = None) -> xr.DataArray:
    """Daily reference (time labelled by its day; 00Z instants for ``00z``) covering start..end.

    ``sampling`` overrides the variable's sampling (``"daily_mean"`` gives the ERA5 daily
    means of a ``00z`` field). ERA5 daily means come from ``ERA5/1x1/daily`` (stamped at
    00Z of their day, written only for complete days), 00Z states from ``ERA5/1x1/6hourly``.
    ERA5 records whose field is identically zero are treated as missing (a guard against
    the zero records of the pre-2026-11-05 archive).
    """
    s = spec(variable)
    sampling = s.sampling if sampling is None else sampling
    root = Path(observation_root)
    # the 00Z state ending the last day lies in the next year's file
    y0 = pd.Timestamp(start).year
    y1 = (pd.Timestamp(end) + pd.Timedelta(days=1)).year if sampling == "00z" else pd.Timestamp(end).year
    if s.reference == "ERA5":
        pattern = _ERA5_6H if sampling == "00z" else _ERA5_DAILY
        files = [root / pattern.format(name=s.ref_name, year=y) for y in range(y0, y1 + 1)]
    elif s.reference == "ESA_CCI":
        files = [root / _ESA.format(year=y) for y in range(y0, y1 + 1)]
    elif s.reference in ("GPCP", "IMERG", "NOAA_OLR"):
        files = sorted(root.glob({"GPCP": _GPCP, "IMERG": _IMERG, "NOAA_OLR": _OLR}[s.reference]))
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
        field = field.where((field != 0).any(("lat", "lon")))
        if sampling == "00z":
            field = field.sel(time=field.time.dt.hour == 0)
        else:
            field = field.assign_coords(time=field.time.dt.floor("D"))
        field.attrs["units"] = ds[s.ref_name].attrs.get("units", "")
    else:
        field = field.assign_coords(time=field.time.dt.floor("D"))
    return field


def open_reference(variable: str, inits, *, max_lead: int = MAX_LEAD,
                   observation_root) -> xr.DataArray:
    """(init, lead, lat, lon) reference on every start date's verification days (lazy)."""
    inits = np.asarray(inits, dtype="datetime64[D]")
    s = spec(variable)
    series = _reference_series(variable, inits.min(), inits.max() + max_lead * _DAY, observation_root)
    out = [_on_leads(series, i, max_lead, s.sampling).drop_vars("time") for i in inits]
    out = xr.concat(out, pd.Index(inits.astype("datetime64[ns]"), name="init"))
    return _convert(out.assign_coords(time=valid_time(inits, max_lead)), variable).rename(variable)


def reference_climatology(variable: str, inits, years, *, max_lead: int = MAX_LEAD, smooth_days: int = 31,
                          observation_root, source: str = "auto") -> xr.DataArray:
    """(init, lead, lat, lon) climatology of the reference over ``years`` on the verification days.

    ``source``:

    * ``"daily"``: for each year the reference is read on the same calendar days as the
      verification days (padded by ``smooth_days // 2``), averaged over the years,
      smoothed with a centred ``smooth_days`` running mean, and put back on each start
      date's leads (Feb 29 is skipped in non-leap years). Exact for the sampling of the
      verified field, e.g. the 00Z states of the pressure-level fields.
    * ``"monthly"``: calendar-month means over ``years`` interpolated linearly between
      mid-month values to the verification days (``_monthly_on_leads``), built the same
      way as the model climatology (``model_climatology``), from the monthly means of the
      daily series (``_reference_monthly``; ERA5 daily means also for ``00z`` fields).
    * ``"auto"`` (default): ``"monthly"`` for every variable, the treatment of
      ``model_climatology``; ``"daily"`` is kept for comparison.
    """
    s = spec(variable)
    if source == "auto":
        source = "monthly"
    if source == "monthly":
        monthly, note = _reference_monthly(variable, years, observation_root)
        out = _monthly_on_leads(monthly, variable, inits, max_lead)
        out.attrs.update(climatology_years=f"{min(years)}-{max(years)}", climatology_source=note,
                         climatology_method="monthly")
        return out
    if source != "daily":
        raise ValueError("source must be 'auto', 'daily' or 'monthly'")
    inits = np.asarray(inits, dtype="datetime64[D]")
    half = smooth_days // 2
    first, last = pd.Timestamp(inits.min()) - pd.Timedelta(days=half), \
        pd.Timestamp(inits.max() + max_lead * _DAY) + pd.Timedelta(days=half)
    ref_year = first.year
    target = pd.date_range(first, last, freq="D")
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
                     climatology_method="daily", climatology_source=f"{s.reference} daily, {smooth_days}-day running mean")
    return out


# --------------------------------------------------------------------------
# Monthly climatologies (reference and AMIP model)
# --------------------------------------------------------------------------


def month_of_records(ds: xr.Dataset) -> pd.DatetimeIndex:
    """Month (as its first day) that each record of a monthly-mean file is labelled with.

    The rule of ``esp_lab.data_access_obs.transform_to_mid_month``: the lower time bound
    when the file has bounds, otherwise the month of the stamp itself; a stamp is never
    reinterpreted (a day-1 stamp is a month start in CDS/C3S products but the end of the
    interval in E3SM history files). Labels that may still be wrong are checked against
    the data with ``verify_month_labels``.
    """
    def months(values):
        out = []
        for v in values:  # numpy datetime64 or cftime (e.g. the no-leap AMIP files)
            v = pd.Timestamp(v) if isinstance(v, np.datetime64) else v
            out.append(pd.Timestamp(v.year, v.month, 1))
        return pd.DatetimeIndex(out)

    for name in ("time_bnds", "time_bounds", "bounds_time", "time_bound"):
        if name in ds.variables and "time" in ds[name].dims:
            return months(ds[name].transpose("time", ...).values[:, 0])
    t = ds["time"].values
    return months(pd.DatetimeIndex(t) if np.issubdtype(t.dtype, np.datetime64) else t)


def _band_means(field: xr.DataArray, mask) -> np.ndarray:
    """(time, band) cos-lat weighted means in 10-degree latitude bands (NaN where a band is empty)."""
    lat = field.lat.values
    w = np.cos(np.deg2rad(lat))
    v = field.transpose("time", "lat", "lon").values.astype(np.float64)
    ok = np.isfinite(v) if mask is None else np.isfinite(v) & mask[None]
    out = []
    for lo in range(-90, 90, 10):
        sel = (lat >= lo) & (lat < lo + 10)
        ww = w[sel][None, :, None] * ok[:, sel]
        tot = ww.sum((1, 2))
        with np.errstate(invalid="ignore", divide="ignore"):
            out.append(np.where(tot > 0, (np.where(ok[:, sel], v[:, sel], 0) * ww).sum((1, 2)) / tot, np.nan))
    return np.stack(out, 1)


def detect_month_shift(field: xr.DataArray, anchor: xr.DataArray, label: str = "", *,
                       candidates=(-1, 0, 1), max_nrmse: float = 0.3, min_ratio: float = 2.0,
                       scale_free: bool = False) -> tuple[int, dict]:
    """Shift (in months) to add to the labels of ``field`` so that it matches ``anchor``.

    Both are monthly (time, lat, lon) series of the same quantity in the same units, with
    ``time`` at month starts; ``anchor`` has trusted labels (independent data, e.g. the
    monthly means of a daily archive, or a file with time bounds). The comparison uses
    10-degree latitude-band means (any grids; on a common grid only the points valid in
    both), with each band's mean over the compared months removed, so a constant offset
    between products does not matter. For each candidate shift: the RMS difference
    normalised by the anchor's variability (nRMSE). The best shift must have nRMSE <
    ``max_nrmse`` and beat every other by ``min_ratio``, else ``ValueError``. For the ERA5
    monthly files the right shift gives nRMSE 0.002-0.03, a wrong one more than 0.5.
    ``scale_free``: compare a related quantity of another sign or scale (e.g. a cloud
    radiative effect against another product's cloud forcing): both series are
    standardised and the anchor's sign is matched, so nRMSE = sqrt(2 (1 - |r|)).
    """
    same_grid = (field.sizes.get("lat") == anchor.sizes.get("lat") and field.sizes.get("lon") == anchor.sizes.get("lon")
                 and np.allclose(field.lat, anchor.lat) and np.allclose(field.lon, anchor.lon))
    def valid(x):  # points valid in every non-empty month
        filled = x.notnull().any(("lat", "lon"))
        return x.isel(time=filled.values).notnull().all("time")

    mask = (valid(field) & valid(anchor)).values if same_grid else None
    fb, ab = _band_means(field, mask), _band_means(anchor, mask)
    fi = field.time.dt.year.values * 12 + field.time.dt.month.values - 1
    ai = anchor.time.dt.year.values * 12 + anchor.time.dt.month.values - 1
    scores = {}
    fi, fb = fi[np.isfinite(fb).any(1)], fb[np.isfinite(fb).any(1)]   # drop empty months
    ai, ab = ai[np.isfinite(ab).any(1)], ab[np.isfinite(ab).any(1)]
    for k in candidates:
        common = np.intersect1d(fi + k, ai)
        if common.size < 12:
            continue
        f = fb[np.searchsorted(fi + k, common)]
        a = ab[np.searchsorted(ai, common)]
        good = np.isfinite(f).all(0) & np.isfinite(a).all(0)
        f, a = f[:, good] - f[:, good].mean(0), a[:, good] - a[:, good].mean(0)
        if scale_free:
            f, a = f / np.sqrt((f ** 2).mean()), a / np.sqrt((a ** 2).mean())
            f = f * np.sign((f * a).sum())
        score = float(np.sqrt(((f - a) ** 2).mean() / (a ** 2).mean())) if good.any() else np.nan
        if np.isfinite(score):
            scores[k] = score
    if not scores:
        raise ValueError(f"{label}: fewer than 12 months overlap the anchor")
    order = sorted(scores, key=scores.get)
    best = order[0]
    if scores[best] >= max_nrmse or any(scores[k] < min_ratio * scores[best] for k in order[1:]):
        raise ValueError(f"{label}: cannot identify the month of each record (nRMSE by shift {scores})")
    return best, scores


def _shift_months(field: xr.DataArray, k: int) -> xr.DataArray:
    if not k:
        return field
    return field.assign_coords(time=(pd.DatetimeIndex(field.time.values) + pd.DateOffset(months=k)).values)


def verify_month_labels(field: xr.DataArray, anchor: xr.DataArray, label: str) -> tuple[xr.DataArray, str]:
    """``field`` with its month labels checked against ``anchor`` (``detect_month_shift``) and
    corrected if needed (with a warning); also returns a note for the provenance attributes."""
    import warnings

    k, scores = detect_month_shift(field, anchor, label)
    note = f"month labels checked (shift {k:+d}, nRMSE {scores[k]:.3f})"
    if k:
        warnings.warn(f"{label}: records are labelled {abs(k)} month(s) {'early' if k > 0 else 'late'}; "
                      f"relabelled by {k:+d} (nRMSE by shift {', '.join(f'{s:+d}: {v:.3f}' for s, v in scores.items())})",
                      stacklevel=2)
    return _shift_months(field, k), note


def _check_months(field: xr.DataArray, years, label: str) -> None:
    have = set(zip(field.time.dt.year.values.tolist(), field.time.dt.month.values.tolist()))
    missing = [(y, m) for y in years for m in range(1, 13) if (y, m) not in have]
    if missing:
        raise FileNotFoundError(f"{label}: {len(missing)} months of {min(years)}-{max(years)} missing, e.g. {missing[:3]}")


def _calendar_months(field: xr.DataArray) -> xr.DataArray:
    out = field.groupby("time.month").mean("time")
    out.attrs = dict(field.attrs)
    return out


def _reference_monthly(variable: str, years, observation_root) -> tuple[xr.DataArray, str]:
    """(month, lat, lon) calendar-month climatology of the reference over ``years``, catalog units:
    the monthly means of its daily series (ERA5: ``ERA5/1x1/daily``, also for ``00z`` fields, so
    the climatology is an all-hours mean like the AMIP one), every month of ``years`` required."""
    s = spec(variable)
    parts = []
    for y in years:
        series = _reference_series(variable, f"{y}-01-01", f"{y}-12-31", observation_root, sampling="daily_mean")
        series = _convert(series.sel(time=series.time.dt.year == y), variable).load()
        months = series.resample(time="MS").mean()
        _check_months(months, [y], f"{variable} {s.reference} daily")
        parts.append(_calendar_months(months))
    monthly = xr.concat(parts, "year").mean("year")
    monthly.attrs = dict(parts[0].attrs)
    src = "ERA5 1x1 daily" if s.reference == "ERA5" else s.reference
    return monthly, f"{src}: monthly means of the daily series"


def _monthly_to_days(monthly: xr.DataArray, instants: pd.DatetimeIndex) -> xr.DataArray:
    """Linear interpolation of a (month, ...) climatology between mid-month values (day 15)."""
    anchors = pd.date_range(instants[0].normalize() - pd.offsets.MonthBegin(2),
                            instants[-1].normalize() + pd.offsets.MonthBegin(2), freq="MS") + pd.Timedelta(days=14)
    values = monthly.sel(month=xr.DataArray(anchors.month, dims="time")).assign_coords(time=anchors.values)
    return values.drop_vars("month").interp(time=instants.values)


def _monthly_on_leads(monthly: xr.DataArray, variable: str, inits, max_lead: int) -> xr.DataArray:
    """(init, lead, lat, lon) of a calendar-month climatology: the value at the midday of each
    verification day (``daily_mean``) or at the 00Z instant ending it (``00z``)."""
    s = spec(variable)
    inits = np.asarray(inits, dtype="datetime64[D]")
    days = pd.date_range(pd.Timestamp(inits.min()), pd.Timestamp(inits.max() + max_lead * _DAY), freq="D")
    if s.sampling == "00z":
        clim = _monthly_to_days(monthly, days).assign_coords(time=days.values)  # 00Z instants, labelled as such
    else:
        clim = _monthly_to_days(monthly, days + pd.Timedelta(hours=12)).assign_coords(time=days.values)
    clim.attrs = dict(monthly.attrs)
    out = [_on_leads(clim, i, max_lead, s.sampling).drop_vars("time") for i in inits]
    out = xr.concat(out, pd.Index(inits.astype("datetime64[ns]"), name="init"))
    return _convert(out.assign_coords(time=valid_time(inits, max_lead)), variable).rename(variable)


AMIP_CASES = ("v3.LR.amip_0101", "v3.LR.amip_0151")
_AMIP_CMIP = "{case}/post/atm/180x360_aave/cmip_ts/monthly/{name}_Amon_*.nc"
_AMIP_LND = "{case}/post/lnd/180x360_aave/ts/monthly/5yr/{name}_*.nc"
# catalog variable -> CMIP monthly name (single-level) of the E3SMv3 DECK AMIP runs
_AMIP_ATM = {"TREFHT": "tas", "TS": "ts", "PSL": "psl", "PRECT": "pr", "FLUT": "rlut", "LHFLX": "hfls",
             "SHFLX": "hfss", "FLDS": "rlds", "FSDS": "rsds", "TMQ": "prw"}
_AMIP_PLEV = {"Z": "zg", "T": "ta", "U": "ua", "V": "va"}  # on plev19 (no specific humidity)
_DENH2O, _DENICE = 1000.0, 917.0


def _year_files(pattern: Path, years) -> list[str]:
    """Files of ``pattern`` whose ``_YYYYMM[-_]YYYYMM`` span overlaps ``years``."""
    import re

    out = []
    for f in sorted(pattern.parent.glob(pattern.name)):
        m = re.search(r"_(\d{4})\d{2}[-_](\d{4})\d{2}\.nc$", f.name)
        if m and int(m.group(1)) <= max(years) and int(m.group(2)) >= min(years):
            out.append(str(f))
    return out


def _open_monthly(files: list[str], name: str, years) -> xr.DataArray:
    """Monthly series of ``name`` in ``years`` padded by a year each side (for relabelling)."""
    span = list(range(min(years) - 1, max(years) + 2))
    parts = []
    for f in files:
        with xr.open_dataset(f) as ds:
            field = _normalize(ds[name]).assign_coords(time=month_of_records(ds).values)
            field = field.sel(time=field.time.dt.year.isin(span)).load()
            field.attrs["units"] = ds[name].attrs.get("units", "")
            parts.append(field)
    out = xr.concat(parts, "time")
    return out.isel(time=np.unique(out.time.values, return_index=True)[1])


_LND_SHIFT: dict[tuple[str, str], tuple[int, str]] = {}


def _amip_land_shift(case: str, years, amip_root) -> tuple[int, str]:
    """Month shift of the ELM monthly ``ts`` files of ``case`` relative to its ``cmip_ts``
    (time bounds): ``TSA`` against ``tas`` over land (``detect_month_shift``); once per case."""
    key = (str(amip_root), case)
    if key not in _LND_SHIFT:
        root = Path(amip_root)
        test = [max(years) - 4, max(years)]
        tsa = _open_monthly(_year_files(root / _AMIP_LND.format(case=case, name="TSA"), test), "TSA", test)
        tas = _open_monthly(_year_files(root / _AMIP_CMIP.format(case=case, name="tas"), test), "tas", test)
        k, scores = detect_month_shift(tsa, tas, f"{case} ELM TSA vs cmip_ts tas")
        if k:
            import warnings
            warnings.warn(f"{case}: ELM monthly files labelled {k:+d} month(s) off cmip_ts; relabelled", stacklevel=2)
        _LND_SHIFT[key] = (k, f"ELM labels checked against cmip_ts tas (shift {k:+d}, nRMSE {scores[k]:.3f})")
    return _LND_SHIFT[key]


def _amip_monthly_series(variable: str, case: str, years, amip_root, land_support_dir) -> xr.DataArray:
    """(time, lat, lon) monthly means of one AMIP case in ``years``, catalog units.

    ``cmip_ts`` files carry time bounds (their labels are taken as they are); the ELM
    ``ts`` labels are checked against them once per case (``_amip_land_shift``)."""
    s = spec(variable)
    root = Path(amip_root)
    name = variable.rstrip("0123456789")
    level = variable[len(name):]
    if variable in _AMIP_ATM:
        cmip = _AMIP_ATM[variable]
        field = _open_monthly(_year_files(root / _AMIP_CMIP.format(case=case, name=cmip), years), cmip, years)
    elif s.component == "atm" and name in _AMIP_PLEV and level:
        cmip = _AMIP_PLEV[name]
        field = _open_monthly(_year_files(root / _AMIP_CMIP.format(case=case, name=cmip), years), cmip, years)
        field = field.sel(plev=float(level) * 100.0, method="nearest", tolerance=1.0).drop_vars("plev")
    elif variable == "SOILM_5CM":
        w = top_soil_weights(land_support_dir)
        nlev = w.sizes["levgrnd"]
        with xr.open_dataset(Path(land_support_dir) / "dzsoi_elm.nc") as ds:
            dz = _normalize(ds["DZSOI"]).isel(levgrnd=slice(0, nlev)).load()
        liq, ice = (_open_monthly(_year_files(root / _AMIP_LND.format(case=case, name=n), years), n, years)
                    .isel(levgrnd=slice(0, nlev)) for n in ("SOILLIQ", "SOILICE"))
        vol = liq / (dz.values[None] * _DENH2O) + ice / (dz.values[None] * _DENICE)  # H2OSOI, m3 m-3
        field = (vol * w.values[None]).sum("levgrnd", min_count=1)
        field.attrs["units"] = "m3 m-3"
    elif s.component == "lnd" and s.model:
        field = _open_monthly(_year_files(root / _AMIP_LND.format(case=case, name=s.model), years), s.model, years)
    else:
        raise KeyError(f"{variable}: no AMIP monthly counterpart")
    if not field.sizes.get("time"):
        raise FileNotFoundError(f"{variable}: no {case} monthly files for {min(years)}-{max(years)}")
    note = "cmip_ts (time bounds)"
    if s.component == "lnd":
        k, note = _amip_land_shift(case, years, amip_root)
        field = _shift_months(field, k)
    field = field.sel(time=field.time.dt.year.isin(list(years)))
    _check_months(field, years, f"{variable} {case}")
    return _convert(field, variable).assign_attrs(month_labels=note)


def model_climatology(variable: str, inits, years, *, max_lead: int = MAX_LEAD, amip_root,
                      cases=AMIP_CASES, land_support_dir) -> xr.DataArray:
    """(init, lead, lat, lon) climatology of the free-running E3SMv3 AMIP runs (``cases``) over
    ``years``, on the verification days.

    Calendar-month means of each case, averaged over the cases and interpolated as the
    reference monthly climatology (``_monthly_on_leads``). The AMIP runs keep the
    hindcasts' grid and boundary forcing (observed SST and sea ice), so forecast
    anomalies from it contain no observed information and remove the model's mean bias
    (not the lead-dependent initialization drift). For ``00z`` fields the monthly mean
    covers all hours (the AMIP runs have no 00Z T/U/Z output), as does the reference
    climatology, so the forecast and observed 00Z anomalies each carry their own
    system's 00Z departure from the daily mean.
    """
    per_case = [_calendar_months(_amip_monthly_series(variable, c, years, amip_root, land_support_dir))
                for c in cases]
    monthly = xr.concat(per_case, "case").mean("case").assign_attrs(per_case[0].attrs)
    out = _monthly_on_leads(monthly, variable, inits, max_lead)
    out.attrs.update(climatology_years=f"{min(years)}-{max(years)}", climatology_method="monthly",
                     climatology_source=f"E3SMv3 AMIP monthly ({', '.join(cases)}; "
                                        f"{per_case[0].attrs.get('month_labels', '')})")
    return out


# --------------------------------------------------------------------------
# Region selection and masks
# --------------------------------------------------------------------------
def _region_box(r: dict) -> dict:
    return dict(lat=[float(x) for x in r["lat"]], lon=[float(x) for x in r["lon"]], mask=r.get("mask"))


def region_catalog(custom: dict | None = None) -> dict[str, dict]:
    """Named regions: ``configs.s2s_config.TCC_REGIONS`` plus ``custom`` boxes.

    ``custom`` = {name: dict(lat=(lat0, lat1), lon=(lon0, lon1), mask=None | "land" | "ocean")}, longitudes
    0-360 (lon0 > lon1 wraps through 0E, e.g. lon=(350, 40) for Europe).
    """
    from configs.s2s_config import TCC_REGIONS

    catalog = {r["name"]: dict(r) for r in TCC_REGIONS}
    for name, r in (custom or {}).items():
        if name in catalog:
            raise ValueError(f"custom region {name!r} would replace the standard region; choose another name")
        lat, lon = tuple(r["lat"]), tuple(r["lon"])
        if not (-90 <= lat[0] < lat[1] <= 90) or not all(0 <= x <= 360 for x in lon) \
                or r.get("mask") not in (None, "land", "ocean"):
            raise ValueError(f"custom region {name!r}: need -90 <= lat0 < lat1 <= 90, lon in 0-360, "
                             f"mask None/'land'/'ocean'; got {r}")
        catalog[name] = dict(name=name, lat=lat, lon=lon, **({"mask": r["mask"]} if r.get("mask") else {}))
    return catalog


def select_regions(names, custom: dict | None = None) -> list[dict]:
    """The regions ``names`` (in that order) from :func:`region_catalog`; unknown names raise, listing the choices."""
    catalog = region_catalog(custom)
    names = [names] if isinstance(names, str) else list(names)
    unknown = [n for n in names if n not in catalog]
    if unknown:
        raise ValueError(f"unknown regions {unknown}; available: {list(catalog)} (or add them to CUSTOM_REGIONS)")
    if len(set(names)) != len(names):
        raise ValueError(f"duplicate regions in {names}")
    return [catalog[n] for n in names]


def _deg(x: float, pos: str, neg: str) -> str:
    return f"{abs(x):g}°{pos if x >= 0 else neg}" if x else "0°"


def region_label(r: dict) -> str:
    """'CONUS (25°N-50°N, 125°W-66°W)'; 'LAND (global, land)'."""
    lat, lon = r["lat"], r["lon"]
    whole_lon = (lon[1] - lon[0]) % 360 == 0
    if tuple(lat) == (-90, 90) and whole_lon:
        box = "global"
    else:
        lo = [x - 360 if x > 180 else x for x in lon]
        lon_txt = "all longitudes" if whole_lon else f"{_deg(lo[0], 'E', 'W')}-{_deg(lo[1], 'E', 'W')}"
        box = f"{_deg(lat[0], 'N', 'S')}-{_deg(lat[1], 'N', 'S')}, {lon_txt}"
    return f"{r['name']} ({box}{', ' + r['mask'] if r.get('mask') else ''})"


def regions_table(regions) -> pd.DataFrame:
    return pd.DataFrame([dict(region=r["name"], lat=tuple(r["lat"]), lon=tuple(r["lon"]), mask=r.get("mask") or "",
                              label=region_label(r)) for r in regions]).set_index("region")


def region_attrs(regions) -> dict:
    """NetCDF attribute recording the region boxes a regional file was computed for."""
    import json

    return {"region_boxes": json.dumps({r["name"]: _region_box(r) for r in regions})}


def cache_has_regions(path, regions) -> bool:
    """True when the regional file ``path`` holds every region of ``regions`` with the same box.

    Files written before ``region_boxes`` was recorded are accepted for the standard regions only.
    """
    import json

    from configs.s2s_config import TCC_REGIONS

    path = Path(path)
    if not path.exists():
        return False
    with xr.open_dataset(path) as ds:
        have = set(ds["region"].values.tolist()) if "region" in ds.coords else set()
        boxes = json.loads(ds.attrs.get("region_boxes", "{}"))
    standard = {r["name"]: _region_box(r) for r in TCC_REGIONS}
    for r in regions:
        recorded = boxes.get(r["name"], standard.get(r["name"]))
        if r["name"] not in have or recorded != _region_box(r):
            print(f"[STALE] {path.name}: region {r['name']} not computed with this box, recomputing")
            return False
    return True


def region_masks(regions, lat, lon, *, land_support_dir) -> xr.DataArray:
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
