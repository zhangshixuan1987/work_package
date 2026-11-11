"""Initial-state checks, daily regional drift and lead-window maps of the ERA5 initial-condition tests.

Model output: the 1-degree time series written by ``post_process/regen_f20tr_ts.sh`` under
``<member dir>/post/{atm,lnd}/180x360_aave/ts/{daily,6hourly}/0yr`` (the same files the hind_lead_eval
notebooks read). Daily values per lead day (lead 1 = the start date), as in ``diagnostic/util/fcst_init_data.py``:
00Z-00Z means of the 6-hourly output for the surface, flux and column fields, and the 00Z state ending each
day for the pressure-level fields (verified at the same instant). Reference: ERA5 1-degree (daily means from
``ERA5/1x1/daily``, 00Z states from ``ERA5/1x1/6hourly``) under the notebook's ``OBSERVATION_ROOT``.

The time handling (dropping the zero-width initial records, labelling means by their interval start) is
reused, unchanged, from ``diagnostic/util/model_daily.py`` and ``diagnostic/util/time_labels.py``; the notebook
puts ``WORK_DIR/diagnostic`` on ``sys.path``.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

GRAVITY = 9.80616
DU_KG_M2 = 2.1415e-5            # one Dobson unit of ozone in kg m-2
M_AIR = 28.966
STREAM_DIRS = {("atm", "daily"): "atm/180x360_aave/ts/daily/0yr",
               ("atm", "6hourly"): "atm/180x360_aave/ts/6hourly/0yr",
               ("lnd", "daily"): "lnd/180x360_aave/ts/daily/0yr"}
_ERA5 = {"daily_mean": "ERA5/1x1/daily/ERA5.daily.en00.{name}.{year}01-{year}12.nc",
         "00z": "ERA5/1x1/6hourly/ERA5.6hourly.en00.{name}.{year}01-{year}12.nc"}


# ---------------------------------------------------------------------------------------------------------
# 1. Atmosphere initial files
# ---------------------------------------------------------------------------------------------------------
# Fields a complete E3SMv3 initial file holds and the ERA5 file lacks; grouped by what they do
KEY_FIELDS = {
    "radiative gases": ("H2OLNZ", "O3", "N2OLNZ", "CH4LNZ"),
    "aerosols": ("num_a1", "so4_a1", "bc_a1", "dst_a1", "pom_a1"),
    "hydrometeor numbers": ("NUMLIQ", "NUMICE"),
    "sea ice / surface": ("ICEFRAC", "TSICE", "TS1", "SNOWHICE"),
}


def ic_grid(full_file) -> xr.Dataset:
    """lat/lon/area of the dynamics (``ncol_d``) and physics (``ncol``) columns of a complete initial file.

    The bare ERA5 files hold the dynamics columns as ``ncol`` without coordinates, in the ``ncol_d`` order
    (``ensemble_perturbations/util/complete_era5_ic.py`` copies them unchanged)."""
    with xr.open_dataset(full_file, decode_times=False) as ds:
        return ds[["lat_d", "lon_d", "area_d", "lat", "lon", "area", "hyai", "hybi", "lev"]].load()


def _dyn(ds: xr.Dataset, grid: xr.Dataset) -> tuple[str, xr.DataArray, xr.DataArray]:
    """Name of the dynamics-column dimension of ``ds``, with its latitude and area."""
    if "ncol_d" in ds.dims:
        return "ncol_d", grid["lat_d"], grid["area_d"]
    if ds.sizes.get("ncol") == grid.sizes["ncol_d"]:
        return "ncol", grid["lat_d"].rename(ncol_d="ncol"), grid["area_d"].rename(ncol_d="ncol")
    raise ValueError("unrecognized column layout")


def _dp(ds: xr.Dataset, grid: xr.Dataset, dim: str) -> xr.DataArray:
    ps = ds["PS"].isel(time=0)
    pi = grid["hyai"] * 1.0e5 + grid["hybi"] * ps
    return pi.diff("ilev").rename(ilev="lev").assign_coords(lev=grid["lev"].values)


def _column(ds, name, dp) -> xr.DataArray | None:
    """Column integral (kg m-2) of a mass mixing ratio; None when the field is absent."""
    if name not in ds:
        return None
    return (ds[name].isel(time=0).assign_coords(lev=dp.lev.values) * dp).sum("lev") / GRAVITY


def _sum_columns(ds, names, dp):
    cols = [c for c in (_column(ds, n, dp) for n in names) if c is not None]
    return sum(cols) if cols else None


def _gmean(x, area):
    return float((x * area).sum() / area.sum()) if x is not None else np.nan


def _level(ds, grid, p_hpa):
    return int(np.abs(grid["lev"].values - p_hpa).argmin())


def ic_summary(path, grid: xr.Dataset) -> tuple[dict, xr.Dataset]:
    """Global (area-weighted) summary of one atmosphere initial file and its zonal means (2-degree bands)."""
    with xr.open_dataset(path, decode_times=False) as ds:
        dim, lat, area = _dyn(ds, grid)
        dp = _dp(ds, grid, dim)
        state = [v for v in ds.data_vars if "time" in ds[v].dims and ({"ncol", "ncol_d"} & set(ds[v].dims))
                 and ds[v].dtype.kind == "f"]
        zero = [v for v in state if not bool((ds[v] != 0).any())]
        absent = {g: [v for v in names if v not in ds or v in zero] for g, names in KEY_FIELDS.items()}
        q = _column(ds, "Q", dp)
        h2o = _column(ds, "H2OLNZ", dp)
        o3 = _column(ds, "O3", dp)
        cols = {
            "TOZ (DU)": None if o3 is None else o3 / DU_KG_M2,
            "column Q (kg m-2)": q,
            "column H2OLNZ (kg m-2)": h2o,
            "sulfate burden (mg m-2)": _scale(_sum_columns(ds, [v for v in ds if v.startswith("so4_a")], dp), 1e6),
            "BC burden (mg m-2)": _scale(_sum_columns(ds, [v for v in ds if v.startswith("bc_a")], dp), 1e6),
            "dust burden (mg m-2)": _scale(_sum_columns(ds, [v for v in ds if v.startswith("dst_a")], dp), 1e6),
            "column NUMLIQ (1e10 m-2)": _scale(_column(ds, "NUMLIQ", dp), 1e-10),
        }
        row = {"file": Path(path).name, "state fields": len(state), "all-zero fields": len(zero)}
        row.update({f"missing {g}": ", ".join(v) or "-" for g, v in absent.items()})
        row.update({k: _gmean(v, area) for k, v in cols.items()})
        for name, mw in (("N2OLNZ", 44.013), ("CH4LNZ", 16.043)):
            c = _column(ds, name, dp)
            row[f"{name[:3]} column mean (ppbv)"] = np.nan if c is None else \
                _gmean(c / (dp.sum("lev") / GRAVITY), area) * M_AIR / mw * 1e9
        for p in (850, 50):
            k = _level(ds, grid, p)
            row[f"T level {grid['lev'].values[k]:.0f} hPa (K)"] = _gmean(ds["T"].isel(time=0, lev=k), area)
        zonal = {k: _zonal(v, lat, area) for k, v in cols.items() if v is not None}
        zonal["T near 50 hPa (K)"] = _zonal(ds["T"].isel(time=0, lev=_level(ds, grid, 50)), lat, area)
    out = xr.Dataset(zonal)
    return row, out


def _scale(x, f):
    return None if x is None else x * f


def _zonal(x: xr.DataArray, lat: xr.DataArray, area: xr.DataArray, width: float = 2.0) -> xr.DataArray:
    """Area-weighted means of a column field in latitude bands of ``width`` degrees."""
    edges = np.arange(-90, 90 + width, width)
    idx = np.clip(np.digitize(np.asarray(lat), edges) - 1, 0, len(edges) - 2)
    w = np.asarray(area, dtype=float)
    num = np.bincount(idx, weights=np.asarray(x, dtype=float) * w, minlength=len(edges) - 1)
    den = np.bincount(idx, weights=w, minlength=len(edges) - 1)
    with np.errstate(invalid="ignore"):
        return xr.DataArray(num / den, dims="lat", coords={"lat": 0.5 * (edges[:-1] + edges[1:])})


def era5_field_check(path, bare_path, grid, fields=("T", "U", "V", "Q", "PS", "CLDLIQ", "CLDICE")) -> dict:
    """Max |difference| of the ERA5 fields between a completed file and the bare ERA5 file (0 = kept)."""
    out = {}
    with xr.open_dataset(path, decode_times=False) as a, xr.open_dataset(bare_path, decode_times=False) as b:
        da, _, _ = _dyn(a, grid)
        db, _, _ = _dyn(b, grid)
        for v in fields:
            if v in a and v in b:
                x = a[v].isel(time=0).values
                y = b[v].isel(time=0).values
                out[v] = float(np.nanmax(np.abs(x - y))) if x.shape == y.shape else np.nan
            else:
                out[v] = np.nan
    return out


# ---------------------------------------------------------------------------------------------------------
# 2. Daily fields on lead days
# ---------------------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class VarSpec:
    component: str
    stream: str        # 6hourly | daily
    sampling: str      # daily_mean | 00z
    era5: str          # ERA5 archive name ("" = no reference)
    units: str
    long_name: str


VARIABLES: dict[str, VarSpec] = {
    "TREFHT": VarSpec("atm", "6hourly", "daily_mean", "TREFHT", "K", "2-m temperature"),
    "TS": VarSpec("atm", "6hourly", "daily_mean", "TS", "K", "surface temperature"),
    "PSL": VarSpec("atm", "6hourly", "daily_mean", "PSL", "hPa", "sea-level pressure"),
    "TMQ": VarSpec("atm", "6hourly", "daily_mean", "TMQ", "kg m-2", "precipitable water"),
    "PRECT": VarSpec("atm", "6hourly", "daily_mean", "PRECT", "mm day-1", "precipitation"),
    "TOZ": VarSpec("atm", "6hourly", "daily_mean", "TCO", "DU", "total column ozone"),
    "FLDS": VarSpec("atm", "6hourly", "daily_mean", "FLDS", "W m-2", "surface downward longwave"),
    "FLDSC": VarSpec("atm", "6hourly", "daily_mean", "FLDSC", "W m-2", "clear-sky surface downward longwave"),
    "FSDSC": VarSpec("atm", "6hourly", "daily_mean", "FSDSC", "W m-2", "clear-sky surface downward shortwave"),
    "FSNTOAC": VarSpec("atm", "6hourly", "daily_mean", "FSNTOAC", "W m-2", "clear-sky net TOA shortwave"),
    "FLUTC": VarSpec("atm", "6hourly", "daily_mean", "FLUTC", "W m-2", "clear-sky outgoing longwave"),
    "SCO": VarSpec("atm", "daily", "daily_mean", "", "DU", "stratospheric column ozone"),
    "TCO": VarSpec("atm", "daily", "daily_mean", "", "DU", "tropospheric column ozone"),
    "FSDS": VarSpec("atm", "6hourly", "daily_mean", "FSDS", "W m-2", "surface downward shortwave"),
    "FLUT": VarSpec("atm", "6hourly", "daily_mean", "FLUT", "W m-2", "outgoing longwave"),
    "LHFLX": VarSpec("atm", "6hourly", "daily_mean", "LHFLX", "W m-2", "latent heat flux"),
    "SHFLX": VarSpec("atm", "6hourly", "daily_mean", "SHFLX", "W m-2", "sensible heat flux"),
    "TGCLDLWP": VarSpec("atm", "6hourly", "daily_mean", "", "kg m-2", "cloud liquid water path"),
    "CLDTOT": VarSpec("atm", "6hourly", "daily_mean", "CLDTOT", "fraction", "total cloud fraction"),
    "T850": VarSpec("atm", "daily", "00z", "T850", "K", "850-hPa temperature (00Z)"),
    "T500": VarSpec("atm", "daily", "00z", "T500", "K", "500-hPa temperature (00Z)"),
    "T100": VarSpec("atm", "daily", "00z", "T100", "K", "100-hPa temperature (00Z)"),
    "T050": VarSpec("atm", "daily", "00z", "", "K", "50-hPa temperature (00Z)"),
    "T010": VarSpec("atm", "daily", "00z", "", "K", "10-hPa temperature (00Z)"),
    "Z500": VarSpec("atm", "daily", "00z", "Z500", "m", "500-hPa height (00Z)"),
    "U200": VarSpec("atm", "daily", "00z", "U200", "m s-1", "200-hPa zonal wind (00Z)"),
    "TSA": VarSpec("lnd", "daily", "daily_mean", "TREFHT", "K", "land 2-m temperature (ELM)"),
}


def spec(variable: str) -> VarSpec:
    if variable not in VARIABLES:
        raise KeyError(f"unknown variable {variable!r}; choose from {list(VARIABLES)}")
    return VARIABLES[variable]


def _convert(field: xr.DataArray, variable: str, source: str) -> xr.DataArray:
    units = str(field.attrs.get("units", ""))
    if variable == "PSL" and units == "Pa":
        field = field * 0.01
    elif variable == "Z500" and units.startswith("m**2"):
        field = field / 9.80665
    elif variable == "PRECT" and units in ("m/s", "m s-1"):
        field = field * 86400.0e3
    elif variable == "PRECT" and source == "ERA5" and units == "m":        # daily accumulation
        field = field * 1.0e3
    elif variable == "TOZ" and units.startswith("kg"):
        field = field / DU_KG_M2
    field = field.astype("float32")
    field.attrs = dict(units=spec(variable).units, long_name=spec(variable).long_name)
    return field


def _normalize(field: xr.DataArray) -> xr.DataArray:
    field = field.rename({k: v for k, v in (("latitude", "lat"), ("longitude", "lon")) if k in field.dims})
    if float(field.lon.min()) < 0:
        field = field.assign_coords(lon=field.lon % 360)
    return field.sortby("lat").sortby("lon")


def lead_days(init, max_lead: int) -> pd.DatetimeIndex:
    return pd.date_range(pd.Timestamp(np.datetime64(init, "D")), periods=max_lead, freq="D")


def _on_leads(field: xr.DataArray, init, max_lead: int, sampling: str) -> xr.DataArray:
    """(lead, lat, lon) of lead days 1..max_lead; days not present are NaN."""
    days = lead_days(init, max_lead)
    if sampling == "00z":
        field = field.reindex(time=(days + pd.Timedelta(days=1)).values)   # 00Z state ending each day
    else:
        field = field.assign_coords(time=field.time.dt.floor("D")).reindex(time=days.values)
    return (field.assign_coords(time=days.values).rename(time="lead")
            .assign_coords(lead=np.arange(1, max_lead + 1), date=("lead", days.values)))


def ts_file(member_path, variable: str) -> Path:
    s = spec(variable)
    d = Path(member_path) / "post" / STREAM_DIRS[s.component, s.stream]
    hits = sorted(d.glob(f"{variable}_[0-9][0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][0-9][0-9].nc"))
    if len(hits) != 1:
        raise FileNotFoundError(f"{'no' if not hits else 'several'} {variable} files in {d}")
    return hits[0]


def stream_ok(member_path, stream: str) -> bool:
    """True when ``regen_f20tr_ts.sh`` finished the member's stream (``atm_daily``, ``atm_6hourly``, ``lnd_daily``)."""
    f = Path(member_path) / "post" / "scripts" / f"regen_ts_{stream}.status"
    return f.is_file() and f.read_text().startswith("OK")


def open_member(member_path, init, variable: str, max_lead: int) -> xr.DataArray:
    """(lead, lat, lon) daily field of one member in catalog units."""
    from util.model_daily import daily_from_6hourly
    from util.time_labels import label_by_interval_start

    s = spec(variable)
    ds = xr.open_dataset(ts_file(member_path, variable), chunks={"time": -1})
    if s.stream == "6hourly":
        field = daily_from_6hourly(ds, variable)
    elif s.sampling == "00z":
        field = ds[variable]
        field = field.sel(time=field.time.dt.hour == 0)
    else:
        field = label_by_interval_start(ds, ds[variable])
    field.attrs["units"] = ds[variable].attrs.get("units", "")
    out = _on_leads(_normalize(field), init, max_lead, s.sampling)
    return _convert(out, variable, "model").rename(variable)


def open_era5(variable: str, init, max_lead: int, observation_root) -> xr.DataArray | None:
    """(lead, lat, lon) ERA5 on the verification days; None when the variable has no reference."""
    s = spec(variable)
    if not s.era5:
        return None
    days = lead_days(init, max_lead)
    years = range(days[0].year, (days[-1] + pd.Timedelta(days=1)).year + 1)
    files = [Path(observation_root) / _ERA5[s.sampling].format(name=s.era5, year=y) for y in years]
    missing = [str(f) for f in files if not f.is_file()]
    if missing:
        raise FileNotFoundError(f"{variable}: missing ERA5 files {missing}")
    ds = xr.open_mfdataset(files, combine="by_coords", chunks={"time": 124})
    field = _normalize(ds[s.era5])
    field.attrs["units"] = ds[s.era5].attrs.get("units", "")
    if s.sampling == "00z":
        field = field.sel(time=field.time.dt.hour == 0)
    units = field.attrs["units"]
    out = _on_leads(field, init, max_lead, s.sampling)
    out.attrs["units"] = units
    return _convert(out, variable, "ERA5").rename(variable)


# ---------------------------------------------------------------------------------------------------------
# 3. Regions and regional means
# ---------------------------------------------------------------------------------------------------------
REGIONS = {
    "GLOBAL": dict(lat=(-90, 90), lon=(0, 360), mask=None),
    "LAND": dict(lat=(-90, 90), lon=(0, 360), mask="land"),
    "OCEAN": dict(lat=(-90, 90), lon=(0, 360), mask="ocean"),
    "NH_LAND": dict(lat=(30, 90), lon=(0, 360), mask="land"),
    "TROPICS": dict(lat=(-20, 20), lon=(0, 360), mask=None),
    "NH_POLAR": dict(lat=(60, 90), lon=(0, 360), mask=None),
    "SH_POLAR": dict(lat=(-90, -60), lon=(0, 360), mask=None),
}


def region_masks(names, lat, lon, land_support_dir) -> xr.DataArray:
    """(region, lat, lon) masks; land/ocean from ``LANDFRAC_E3SM_180x360.nc`` (land fraction > 0.5)."""
    with xr.open_dataset(Path(land_support_dir) / "LANDFRAC_E3SM_180x360.nc") as ds:
        lf = _normalize(ds["LANDFRAC"].squeeze(drop=True)).load()
    lf = lf.reindex(lat=lat, lon=lon, method="nearest", tolerance=0.01)
    out = []
    for n in names:
        r = REGIONS[n]
        lo0, lo1 = r["lon"]
        lat_ok = (lat >= r["lat"][0]) & (lat <= r["lat"][1])
        lon_ok = ((lon >= lo0) & (lon <= lo1)) if lo0 <= lo1 else ((lon >= lo0) | (lon <= lo1))
        m = lat_ok & lon_ok
        if r["mask"] == "land":
            m = m & (lf > 0.5)
        elif r["mask"] == "ocean":
            m = m & (lf < 0.5)
        out.append(m.transpose("lat", "lon"))
    return xr.concat(out, pd.Index(list(names), name="region"))


def regional_mean(field: xr.DataArray, masks: xr.DataArray) -> xr.DataArray:
    """cos-lat weighted mean over each region's valid points; dims (region, ...)."""
    w = np.cos(np.deg2rad(field.lat)) * masks
    valid = field.notnull()
    num = (field.fillna(0) * w).sum(("lat", "lon"))
    den = (valid * w).sum(("lat", "lon"))
    return (num / den.where(den > 0)).transpose("region", ...)


# ---------------------------------------------------------------------------------------------------------
# 4. Lead windows
# ---------------------------------------------------------------------------------------------------------
def window_means(field: xr.DataArray, windows: dict) -> xr.DataArray:
    """Mean over the lead days of each window {label: (first, last)} (dim ``window``)."""
    out = [field.sel(lead=slice(a, b)).mean("lead") for a, b in windows.values()]
    out = [o.drop_vars([c for c in ("date",) if c in o.coords]) for o in out]
    return xr.concat(out, pd.Index(list(windows), name="window"))


def noise_scale(members: xr.DataArray) -> xr.DataArray:
    """Standard deviation of the difference of two independent members, sqrt(2) x member spread.

    Used as the yardstick for a single-member difference A - B: beyond the first days, differences
    smaller than this are indistinguishable from the chaotic divergence of two members."""
    return np.sqrt(2.0) * members.std("member", ddof=1)


def member_vs_ensemble(member: xr.DataArray, ensemble: xr.DataArray, windows: dict) -> pd.DataFrame:
    """Score one run against the members of another ensemble, by region and lead window.

    ``member`` is (region, lead), ``ensemble`` (member, region, lead): regional bias against ERA5 (or the
    raw value without a reference). Returns per region and window the run's value, the ensemble mean and
    spread, their difference, ``z`` = difference / (spread x sqrt(1 + 1/N)), the spread of the difference
    between a new member and the mean of N members, and ``rank``, the number of members below the run
    (0 or N: outside every member, probability 1/(N+1) on each side for an exchangeable member)."""
    m = window_means(member, windows)
    e = window_means(ensemble, windows)
    n = ensemble.sizes["member"]
    mean, sd = e.mean("member"), e.std("member", ddof=1)
    rank = (e < m).sum("member")
    rows = []
    for r in m.region.values:
        for w in windows:
            x, mu, s = (float(a.sel(region=r, window=w)) for a in (m, mean, sd))
            rows.append(dict(region=r, window=w, run=x, ensemble_mean=mu, ensemble_spread=s, difference=x - mu,
                             z=(x - mu) / (s * np.sqrt(1 + 1 / n)) if s > 0 else np.nan,
                             rank=int(rank.sel(region=r, window=w)), members=n))
    return pd.DataFrame(rows)


__all__ = ["KEY_FIELDS", "VARIABLES", "REGIONS", "ic_grid", "ic_summary", "era5_field_check", "spec",
           "open_member", "open_era5", "stream_ok", "region_masks", "regional_mean", "window_means",
           "noise_scale", "member_vs_ensemble", "lead_days"]
