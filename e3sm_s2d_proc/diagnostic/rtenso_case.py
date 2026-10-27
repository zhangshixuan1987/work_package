"""Shared helpers for the RTENSO2026 2015/16 El Nino case-study notebooks.

The RTENSO2026 2015 hindcasts (``clim`` and ``obs`` land/sea-ice/river ICs,
Feb/May/Aug/Nov starts, 10 members, 24 months) use the same ERA5 atmosphere
and CFSR ocean ICs as the Reanalysis (BruteForce) set, so Reanalysis provides
the lead-time climatology that removes model drift:

* May and Nov starts use the Reanalysis 1981-2010 climatology for the same
  start month and lead.
* Feb and Aug starts have no Reanalysis climatology (only 1980 and 1985), so
  each lead uses the May/Nov climatology for the same target calendar month at
  the nearest lead; when a May and a Nov lead are equally near (3 months either
  side) the two are averaged.

Paths follow ``esp_lab.env_paths``; caches are written with
``atomic_to_netcdf`` and reused only when their provenance attributes match
(``derivation_mode`` = "auto" | "rebuild" | "require", as in ESP-Lab).
"""

from __future__ import annotations

import glob
from datetime import timedelta
from pathlib import Path

import cftime
import numpy as np
import xarray as xr

from esp_lab import env_paths
from esp_lab import data_access_obs as dao
from esp_lab.diagnostics.sst_index import BASIC_REGIONS as SST_INDEX_REGIONS
from esp_lab.utils.filename_utils import figure_filename
from esp_lab.utils.netcdf_utils import atomic_to_netcdf, load_netcdf
from esp_lab.utils.regrid_utils import make_regridder

CACHE_VERSION = 1
OBS_CACHE_VERSION = 2   # 2: TS (HadISST2 SST) masked over land, as ESP-Lab's sst_land_mask
GRID = "180x360_aave"
GRID_TAG = "1x1deg"
NENS = 10
NLEAD = 24
CASE_YEAR = 2015
CLIM_YEARS = (1981, 2010)
OBS_YEARS = (2015, 2017)
MODEL_PREFIX = "WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5"

# Display name -> case settings. "source" is the s2d_diag folder/file prefix.
CASES = {
    "RTENSO2026-clim": {
        "case_prefix": f"{MODEL_PREFIX}_RTENSO2026_clim",
        "source": "RTENSO2026_clim",
        "init_months": (2, 5, 8, 11),
        "color": "tab:blue",
    },
    "RTENSO2026-obs": {
        "case_prefix": f"{MODEL_PREFIX}_RTENSO2026_obs",
        "source": "RTENSO2026_obs",
        "init_months": (2, 5, 8, 11),
        "color": "tab:orange",
    },
    "Reanalysis-2015": {
        "case_prefix": f"{MODEL_PREFIX}_BruteForce",
        "source": "Reanalysis",
        "init_months": (5, 11),
        "color": "tab:green",
    },
}
REFERENCE = {"case_prefix": f"{MODEL_PREFIX}_BruteForce", "source": "Reanalysis", "init_months": (5, 11)}

# Field -> realm, file variable, obs product/variable, unit conversion to display units.
FIELDS = {
    "TS": {"realm": "atm", "obs": ("HadISST2", "sst"), "scale": 1.0, "offset": -273.15, "units": "degC", "long_name": "Surface temperature (SST over ocean)"},
    "TREFHT": {"realm": "atm", "obs": ("ERA5", "TREFHT"), "scale": 1.0, "offset": -273.15, "units": "degC", "long_name": "2-m air temperature"},
    "PRECT": {"realm": "atm", "obs": ("GPCP_v2.3", "PRECT"), "scale": 86400.0 * 1000.0, "offset": 0.0, "units": "mm/day", "long_name": "Precipitation"},
    "PSL": {"realm": "atm", "obs": ("ERA5", "PSL"), "scale": 0.01, "offset": 0.0, "units": "hPa", "long_name": "Sea-level pressure"},
    "U850": {"realm": "atm", "var": "U", "plev": 85000.0, "obs": ("ERA5", "ua"), "scale": 1.0, "offset": 0.0, "units": "m/s", "long_name": "850-hPa zonal wind"},
    "SOILWATER_10CM": {"realm": "lnd", "obs": None, "scale": 1.0, "offset": 0.0, "units": "kg/m2", "long_name": "Top 10-cm soil water"},
    "H2OSNO": {"realm": "lnd", "obs": None, "scale": 1.0, "offset": 0.0, "units": "mm", "long_name": "Snow water equivalent"},
    "TWS": {"realm": "lnd", "obs": None, "scale": 1.0, "offset": 0.0, "units": "mm", "long_name": "Terrestrial water storage"},
    "TSA": {"realm": "lnd", "obs": None, "scale": 1.0, "offset": -273.15, "units": "degC", "long_name": "Land 2-m air temperature"},
    "SICE": {"realm": "ice", "var": "timeMonthly_avg_iceAreaCell", "obs": ("HadISST2", "sic"), "scale": 1.0, "offset": 0.0, "units": "fraction", "long_name": "Sea-ice concentration"},
}

# Verification regions [lon0, lon1, lat0, lat1] (0-360 longitudes).
REGIONS = {
    "Global": [0.0, 360.0, -60.0, 60.0],
    "Tropics": [0.0, 360.0, -20.0, 20.0],
    "TropPacific": [120.0, 280.0, -20.0, 20.0],
    "NHextratropics": [0.0, 360.0, 20.0, 90.0],
    "SHextratropics": [0.0, 360.0, -90.0, -20.0],
    "NorthAmerica": [190.0, 310.0, 15.0, 75.0],
}
NINO_REGIONS = {k: SST_INDEX_REGIONS[k]["lonlat"] for k in ("Nino3.4", "Nino3", "Nino4", "Nino12")}


# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
def diag_dir(source, *parts):
    return env_paths.s2d_diag_root().joinpath(source, *parts)


def figure_dir():
    path = env_paths.figure_root().parent / "e3sm-s2d_rtenso2026_case"
    path.mkdir(parents=True, exist_ok=True)
    return path


def save_figure(fig, *parts, dpi=150):
    path = figure_dir() / figure_filename(*parts)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    return path


def case_dir(case_prefix, year, init_month):
    return env_paths.raw_model_root() / f"{case_prefix}_{year}{init_month:02d}0100"


def member_file(case_prefix, year, init_month, member, field):
    info = FIELDS[field]
    var = info.get("var", field)
    tsdir = case_dir(case_prefix, year, init_month) / f"EN{member:02d}" / "post" / info["realm"] / GRID / "ts" / "monthly" / "2yr"
    matches = sorted(glob.glob(str(tsdir / f"{var}_{year}{init_month:02d}_*.nc")))
    return Path(matches[0]) if matches else None


# -----------------------------------------------------------------------------
# Model data
# -----------------------------------------------------------------------------
def _target_times(year, init_month, nlead=NLEAD):
    times = []
    for lead in range(nlead):
        months = init_month - 1 + lead
        times.append(cftime.DatetimeNoLeap(year + months // 12, months % 12 + 1, 15))
    return times


def _to_display_units(da, field):
    info = FIELDS[field]
    out = (da * info["scale"] + info["offset"]).astype("float32")
    out.attrs = {"units": info["units"], "long_name": info["long_name"]}
    return out


def open_member(case_prefix, year, init_month, member, field, nlead=NLEAD):
    """One member's monthly series as (lead, lat, lon) in display units, or None."""
    path = member_file(case_prefix, year, init_month, member, field)
    if path is None:
        return None
    info = FIELDS[field]
    var = info.get("var", field)
    with xr.open_dataset(path, decode_times=True) as ds:
        # Month of each record: EAM and ELM stamp the month end and carry bounds
        # (time_bnds / time_bounds), so use the bounds midpoint; MPAS-SI has no
        # bounds and stamps near the month start, so shift by two weeks.
        bnds_name = ds["time"].attrs.get("bounds")
        if bnds_name in ds:
            b = ds[bnds_name].values
            months = [b0 + (b1 - b0) / 2 for b0, b1 in b]
        else:
            months = [t + timedelta(days=14) for t in ds["time"].values]
        if len(months) < nlead or (months[0].year, months[0].month) != (year, init_month):
            raise ValueError(f"{path} does not start at {year}-{init_month:02d} with {nlead} months")
        da = ds[var].isel(time=slice(0, nlead))
        if "plev" in da.dims:
            da = da.sel(plev=info["plev"])
        da = da.load()
    da = da.rename(time="lead").assign_coords(lead=np.arange(nlead))
    return _to_display_units(da.drop_vars([c for c in da.coords if c not in ("lead", "lat", "lon")]), field)


def open_hindcast(case_prefix, year, init_month, field, members=range(NENS), nlead=NLEAD):
    """Ensemble as (member, lead, lat, lon) with a target-month ``time`` coordinate."""
    arrays, found = [], []
    for member in members:
        da = open_member(case_prefix, year, init_month, member, field, nlead)
        if da is not None:
            arrays.append(da)
            found.append(member)
    if not arrays:
        return None
    out = xr.concat(arrays, dim=xr.DataArray(found, dims="member", name="member"))
    out = out.assign_coords(time=("lead", _target_times(year, init_month, nlead)))
    out.attrs.update(case_prefix=case_prefix, init_year=year, init_month=init_month)
    return out


def land_fraction():
    """Static EAM land fraction on the 1x1 grid, from the first Reanalysis TS file."""
    path = member_file(REFERENCE["case_prefix"], CLIM_YEARS[0], REFERENCE["init_months"][0], 0, "TS")
    with xr.open_dataset(path) as ds:
        return ds["LANDFRAC"].isel(time=0, drop=True).load()


# -----------------------------------------------------------------------------
# Reanalysis lead-time climatology
# -----------------------------------------------------------------------------
def _cache_ok(path, expected):
    if not path.is_file():
        return False
    with xr.open_dataset(path) as ds:
        return all(str(ds.attrs.get(k)) == str(v) for k, v in expected.items())


def leadtime_climatology(field, init_month, years=CLIM_YEARS, mode="auto"):
    """Reanalysis ensemble-mean climatology (lead, lat, lon) for a May or Nov start."""
    if init_month not in REFERENCE["init_months"]:
        raise ValueError(f"Reanalysis climatology exists only for init months {REFERENCE['init_months']}")
    realm = FIELDS[field]["realm"]
    path = diag_dir(REFERENCE["source"], "leadtime_clim", realm, field,
                    f"{REFERENCE['source']}_init{init_month:02d}_{field}_leadclim_y{years[0]}-{years[1]}_m{NENS}_l{NLEAD}_{GRID_TAG}.nc")
    expected = {"cache_version": CACHE_VERSION, "case_prefix": REFERENCE["case_prefix"],
                "clim_start_year": years[0], "clim_end_year": years[1], "nlead": NLEAD, "units": FIELDS[field]["units"]}
    if mode != "rebuild" and _cache_ok(path, expected):
        return load_netcdf(path)[field]
    if mode == "require":
        raise FileNotFoundError(f"missing or stale lead-time climatology cache {path}")

    total, count, used, skipped = None, None, [], []
    for year in range(years[0], years[1] + 1):
        ens = open_hindcast(REFERENCE["case_prefix"], year, init_month, field)
        if ens is None:
            skipped.append(year)
            continue
        ens_mean = ens.mean("member")
        valid = ens_mean.notnull()
        total = ens_mean.fillna(0.0) if total is None else total + ens_mean.fillna(0.0)
        count = valid.astype("int32") if count is None else count + valid
        used.append(year)
        print(f"  {field} init{init_month:02d}: {year} ({ens.sizes['member']} members)")
    if not used:
        raise FileNotFoundError(f"no Reanalysis {field} hindcasts for init{init_month:02d} in {years}")
    clim = (total / count).where(count > 0).drop_vars("time", errors="ignore")
    clim.attrs = {"units": FIELDS[field]["units"], "long_name": FIELDS[field]["long_name"]}
    ds = clim.to_dataset(name=field)
    ds.attrs.update(expected, years_used=",".join(map(str, used)), years_skipped=",".join(map(str, skipped)) or "none")
    atomic_to_netcdf(ds, path)
    return clim


def nearest_reference_leads(init_month, lead, ref_months=REFERENCE["init_months"], nlead=NLEAD):
    """(ref_init, ref_lead) pairs with the same target month as (init_month, lead) and the nearest lead."""
    candidates = []
    for ref in ref_months:
        offset = (init_month - ref) % 12
        for ref_lead in (lead + offset, lead + offset - 12):
            if 0 <= ref_lead < nlead:
                candidates.append((abs(ref_lead - lead), ref, ref_lead))
    best = min(c[0] for c in candidates)
    return [(ref, ref_lead) for dist, ref, ref_lead in candidates if dist == best]


def _target_months(init_month, nlead=NLEAD):
    return [(init_month - 1 + lead) % 12 + 1 for lead in range(nlead)]


def obs_monthly_climatology(field, clim_years=CLIM_YEARS, mode="auto"):
    da = obs_monthly(field, mode=mode)
    return da.sel(time=slice(f"{clim_years[0]}-01", f"{clim_years[1]}-12")).groupby("time.month").mean("time")


def matched_climatology(field, init_month, method="lead_drift", mode="auto"):
    """Lead-time climatology for any start month.

    May/Nov use the Reanalysis climatology directly. Other starts use
    ``method``:

    * "lead_drift" (default): observed 1981-2010 climatology of each target
      month plus the model drift at the same lead, averaged over the May and
      Nov starts (drift = Reanalysis lead climatology minus the observed
      climatology of its target month). Keeps the fast early-lead drift.
    * "target_month": the May/Nov climatology for the same target month at
      the nearest lead (averaged when equally near). Simple, but at short
      leads it borrows a climatology that has already drifted.

    Fields without an observational product always use "target_month".
    """
    if init_month in REFERENCE["init_months"]:
        clim = leadtime_climatology(field, init_month, mode=mode)
        clim.attrs["climatology_method"] = f"Reanalysis init{init_month:02d} {CLIM_YEARS[0]}-{CLIM_YEARS[1]}"
        return clim
    refs = {m: leadtime_climatology(field, m, mode=mode) for m in REFERENCE["init_months"]}
    if method == "lead_drift" and FIELDS[field]["obs"] is not None:
        obs_clim = obs_monthly_climatology(field, mode=mode)
        drifts = [refs[m] - obs_clim.sel(month=xr.DataArray(_target_months(m), dims="lead")).drop_vars("month")
                  for m in REFERENCE["init_months"]]
        drift = sum(drifts) / len(drifts)
        clim = obs_clim.sel(month=xr.DataArray(_target_months(init_month), dims="lead")).drop_vars("month") + drift
        clim = clim.assign_coords(lead=np.arange(NLEAD))
        clim.attrs = {"units": FIELDS[field]["units"],
                      "climatology_method": (f"{FIELDS[field]['obs'][0]} {CLIM_YEARS[0]}-{CLIM_YEARS[1]} target-month "
                                             "climatology + Reanalysis init05/init11 mean drift at the same lead")}
        return clim
    if method not in ("lead_drift", "target_month"):
        raise ValueError(f"unknown climatology method {method!r}")
    leads, mapping = [], []
    for lead in range(NLEAD):
        pairs = nearest_reference_leads(init_month, lead)
        leads.append(sum(refs[ref].sel(lead=ref_lead, drop=True) for ref, ref_lead in pairs) / len(pairs))
        mapping.append("+".join(f"init{ref:02d}L{ref_lead}" for ref, ref_lead in pairs))
    clim = xr.concat(leads, dim="lead").assign_coords(lead=np.arange(NLEAD))
    clim.attrs = {"units": FIELDS[field]["units"],
                  "climatology_method": "target-month matched Reanalysis init05/init11, nearest lead",
                  "lead_map": ";".join(mapping)}
    return clim


def hindcast_anomaly(case_prefix, init_month, field, year=CASE_YEAR, method="lead_drift", mode="auto"):
    """Drift-corrected anomalies (member, lead, lat, lon), or None when the case is missing."""
    ens = open_hindcast(case_prefix, year, init_month, field)
    if ens is None:
        return None
    clim = matched_climatology(field, init_month, method=method, mode=mode)
    anom = ens - clim
    anom.attrs = {"units": FIELDS[field]["units"], "case": case_prefix, "init_month": init_month,
                  "climatology_method": clim.attrs.get("climatology_method", "")}
    return anom


# -----------------------------------------------------------------------------
# Observations
# -----------------------------------------------------------------------------
def model_grid():
    lat = np.arange(-89.5, 90.0, 1.0)
    lon = np.arange(0.5, 360.0, 1.0)
    return xr.Dataset(coords={"lat": ("lat", lat), "lon": ("lon", lon)})


def obs_monthly(field, years=(CLIM_YEARS[0], OBS_YEARS[1]), mode="auto"):
    """Observed monthly values on the model 1x1 grid (time, lat, lon) in display units."""
    product, obs_field = FIELDS[field]["obs"]
    path = diag_dir("observations", "case_study", FIELDS[field]["realm"], field,
                    f"{product}_{field}_monthly_y{years[0]}-{years[1]}_{GRID_TAG}.nc")
    expected = {"cache_version": OBS_CACHE_VERSION, "product": product, "start_year": years[0],
                "end_year": years[1], "units": FIELDS[field]["units"]}
    if mode != "rebuild" and _cache_ok(path, expected):
        return load_netcdf(path)[field]
    if mode == "require":
        raise FileNotFoundError(f"missing or stale observation cache {path}")

    # Chunk in time so 0.25-degree ERA5 is regridded a year at a time on shared login nodes.
    ds = dao.get_monthly_data(obs_dir=str(env_paths.obs_root()), product=product, field=obs_field,
                              start_year=str(years[0]), end_year=str(years[1]), chunks={"time": 12})
    da = ds[list(ds.data_vars)[0]]
    if "plev" in da.dims:
        plev = da["plev"]
        target = FIELDS[field]["plev"] if float(plev.max()) > 2000 else FIELDS[field]["plev"] / 100.0
        da = da.sel(plev=target)
    da = _obs_to_display_units(da, field)
    if da.sizes["lat"] != 180 or da.sizes["lon"] != 360:
        src = da.to_dataset(name=field)
        da = make_regridder(src, model_grid(), method="conservative")(src)[field]
    da = da.drop_vars([c for c in da.coords if c not in ("time", "lat", "lon")])
    if field in ("TS", "SICE"):  # HadISST2 sst/sic are filled over land; keep ocean points (EAM LANDFRAC < 0.5)
        da = da.where(land_fraction() < 0.5)
    da.attrs = {"units": FIELDS[field]["units"], "long_name": FIELDS[field]["long_name"], "product": product}
    out = da.to_dataset(name=field)
    out.attrs.update(expected)
    atomic_to_netcdf(out, path)
    return da


def obs_anomaly(field, years=OBS_YEARS, clim_years=CLIM_YEARS, mode="auto"):
    """Observed monthly anomalies from the 1981-2010 monthly climatology, on the model grid."""
    product = FIELDS[field]["obs"][0]
    path = diag_dir("observations", "case_study", FIELDS[field]["realm"], field,
                    f"{product}_{field}_anomaly_y{years[0]}-{years[1]}_clim_{clim_years[0]}_{clim_years[1]}_{GRID_TAG}.nc")
    expected = {"cache_version": OBS_CACHE_VERSION, "product": product, "clim_start_year": clim_years[0],
                "clim_end_year": clim_years[1], "units": FIELDS[field]["units"]}
    if mode != "rebuild" and _cache_ok(path, expected):
        return load_netcdf(path)[field]
    if mode == "require":
        raise FileNotFoundError(f"missing or stale observation cache {path}")

    da = obs_monthly(field, years=(min(clim_years[0], years[0]), max(clim_years[1], years[1])), mode=mode)
    clim = da.sel(time=slice(f"{clim_years[0]}-01", f"{clim_years[1]}-12")).groupby("time.month").mean("time")
    anom = (da.sel(time=slice(f"{years[0]}-01", f"{years[1]}-12")).groupby("time.month") - clim).drop_vars("month")
    anom.attrs = {"units": FIELDS[field]["units"], "long_name": FIELDS[field]["long_name"], "product": product}
    out = anom.to_dataset(name=field)
    out.attrs.update(expected)
    atomic_to_netcdf(out, path)
    return anom


def _obs_to_display_units(da, field):
    units = str(da.attrs.get("units", "")).lower()
    if field in ("TS", "TREFHT") and units in ("k", "kelvin"):
        da = da - 273.15
    elif field == "PSL" and units in ("pa",):
        da = da / 100.0
    elif field == "SICE" and (units in ("%", "percent") or float(da.max()) > 1.5):
        da = da / 100.0  # HadISST2 sic is in percent although its units say "Fraction"
    da.attrs["units"] = FIELDS[field]["units"]
    return da


def match_obs(model, obs):
    """Observed values at each model lead's target month (dims follow ``model.time``)."""
    target = model["time"]
    picked = obs.sel(time=[cftime.DatetimeNoLeap(t.year, t.month, 15) for t in target.values], method="nearest")
    return picked.rename(time="lead").assign_coords(lead=model["lead"].values, time=("lead", target.values))


# -----------------------------------------------------------------------------
# Statistics
# -----------------------------------------------------------------------------
def area_weights(lat):
    return np.cos(np.deg2rad(lat))


def region_subset(da, region):
    lon0, lon1, lat0, lat1 = REGIONS[region] if isinstance(region, str) else region
    sub = da.sel(lat=slice(lat0, lat1))
    if lon0 <= lon1:
        return sub.sel(lon=slice(lon0, lon1))
    return xr.concat([sub.sel(lon=slice(lon0, 360.0)), sub.sel(lon=slice(0.0, lon1))], dim="lon")


def box_mean(da, lonlat):
    sub = region_subset(da, lonlat)
    return sub.weighted(area_weights(sub["lat"]).fillna(0)).mean(("lat", "lon"))


def pattern_acc(fcst, obs, region="Global"):
    """Area-weighted, centered spatial anomaly correlation over a region."""
    f, o = region_subset(fcst, region), region_subset(obs, region)
    valid = f.notnull() & o.notnull()
    w = (area_weights(f["lat"]) * xr.ones_like(f)).where(valid)
    f, o = f.where(valid), o.where(valid)
    fa = f - (f * w).sum(("lat", "lon")) / w.sum(("lat", "lon"))
    oa = o - (o * w).sum(("lat", "lon")) / w.sum(("lat", "lon"))
    cov = (w * fa * oa).sum(("lat", "lon"))
    return cov / np.sqrt((w * fa**2).sum(("lat", "lon")) * (w * oa**2).sum(("lat", "lon")))


def pattern_rmse(fcst, obs, region="Global"):
    f, o = region_subset(fcst, region), region_subset(obs, region)
    valid = f.notnull() & o.notnull()
    w = (area_weights(f["lat"]) * xr.ones_like(f)).where(valid)
    return np.sqrt((w * (f - o) ** 2).sum(("lat", "lon")) / w.sum(("lat", "lon")))


def seasonal_leads(da, nseason=8):
    """Non-overlapping 3-month means over leads (season 1 = leads 0-2), keeping the centre-month time."""
    seasons = []
    for s in range(nseason):
        chunk = da.isel(lead=slice(3 * s, 3 * s + 3))
        if chunk.sizes["lead"] < 3:
            break
        mean = chunk.mean("lead")
        if "time" in chunk.coords:
            mean = mean.assign_coords(time=chunk["time"].values[1])
        seasons.append(mean)
    return xr.concat(seasons, dim="season").assign_coords(season=np.arange(1, len(seasons) + 1))


def target_mean(da, months):
    """Mean over the leads whose target month is in ``months`` [(year, month), ...]; None if not all covered."""
    want = {(y, m) for y, m in months}
    keep = [i for i, t in enumerate(da["time"].values) if (t.year, t.month) in want]
    if len(keep) != len(want):
        return None
    return da.isel(lead=keep).mean("lead")


def obs_target_mean(obs, months):
    times = [cftime.DatetimeNoLeap(y, m, 15) for y, m in months]
    return obs.sel(time=times, method="nearest").mean("time")


def hemispheric_area(frac):
    """Sea-ice area (10^6 km^2) per hemisphere from a concentration field on the 1x1 grid."""
    radius_km = 6371.0
    cell = (radius_km**2) * np.deg2rad(1.0) ** 2 * np.cos(np.deg2rad(frac["lat"]))
    area = (frac.fillna(0.0) * cell).sum("lon")
    return xr.concat([area.sel(lat=slice(0, 90)).sum("lat"), area.sel(lat=slice(-90, 0)).sum("lat")],
                     dim=xr.DataArray(["NH", "SH"], dims="hemisphere", name="hemisphere")) / 1e6


def target_season_label(time):
    """Label a 3-month season by its centre month, e.g. Jan 2016 -> 'DJF 2015/16'."""
    months = "JFMAMJJASOND"
    m = time.month - 1
    name = months[(m - 1) % 12] + months[m] + months[(m + 1) % 12]
    if m in (0, 11):
        first = time.year - 1 if m == 0 else time.year
        return f"{name} {first}/{str(first + 1)[2:]}"
    return f"{name} {time.year}"


def welch_t(a, b, dim="member"):
    """Welch t statistic and two-sided p value along ``dim`` (vectorised, NaN-aware)."""
    from scipy import stats as sstats
    na, nb = a.count(dim), b.count(dim)
    va, vb = a.var(dim, ddof=1) / na, b.var(dim, ddof=1) / nb
    t = (a.mean(dim) - b.mean(dim)) / np.sqrt(va + vb)
    dof = (va + vb) ** 2 / (va**2 / (na - 1) + vb**2 / (nb - 1))
    p = xr.apply_ufunc(lambda tt, df: 2.0 * sstats.t.sf(np.abs(tt), df), t, dof)
    return t, p


__all__ = [
    "CASES", "FIELDS", "REGIONS", "NINO_REGIONS", "REFERENCE", "CASE_YEAR", "CLIM_YEARS", "NENS", "NLEAD",
    "open_hindcast", "leadtime_climatology", "matched_climatology", "nearest_reference_leads",
    "hindcast_anomaly", "obs_monthly", "obs_anomaly", "match_obs", "box_mean", "pattern_acc", "pattern_rmse",
    "seasonal_leads", "target_mean", "obs_target_mean", "hemispheric_area", "target_season_label",
    "welch_t", "region_subset", "area_weights",
    "figure_dir", "save_figure", "diag_dir",
]
