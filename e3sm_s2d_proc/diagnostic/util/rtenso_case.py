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
import time
from datetime import timedelta
from functools import lru_cache
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
# RTENSO2026 post-processing has its own root; other sets stay under env_paths.raw_model_root().
RTENSO2026_POST_ROOT = Path("/global/cfs/cdirs/e3sm/S2S2D/RTENSO2026/post_process")

# Display name -> case settings. "source" is the s2d_diag folder/file prefix.
CASES = {
    "Reanalysis-CLIM": {
        "case_prefix": f"{MODEL_PREFIX}_RTENSO2026_clim",
        "source": "RTENSO2026_clim",
        "init_months": (2, 5, 8, 11),
        "color": "tab:blue",
    },
    "Reanalysis-OBS": {
        "case_prefix": f"{MODEL_PREFIX}_RTENSO2026_obs",
        "source": "RTENSO2026_obs",
        "init_months": (2, 5, 8, 11),
        "color": "tab:orange",
    },
    # Placeholder: obs ICs with NSIDCv6 sea ice; set init_months once the runs are post-processed.
    "Reanalysis-NSIDC": {
        "case_prefix": f"{MODEL_PREFIX}_RTENSO2026_obs_ice_nsidc",
        "source": "RTENSO2026_obs_ice_nsidc",
        "init_months": (),
        "color": "tab:red",
    },
    "Reanalysis-FOSIRL": {
        "case_prefix": f"{MODEL_PREFIX}_BruteForce",
        "source": "Reanalysis",
        "init_months": (2, 5, 8, 11),
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
    "U200": {"realm": "atm", "var": "U", "plev": 20000.0, "obs": ("ERA5", "ua"), "scale": 1.0, "offset": 0.0, "units": "m/s", "long_name": "200-hPa zonal wind"},
    "V200": {"realm": "atm", "var": "V", "plev": 20000.0, "obs": ("ERA5", "va"), "scale": 1.0, "offset": 0.0, "units": "m/s", "long_name": "200-hPa meridional wind"},
    "Z500": {"realm": "atm", "var": "Z3", "plev": 50000.0, "obs": ("ERA5", "zg"), "scale": 1.0, "offset": 0.0, "units": "m", "long_name": "500-hPa geopotential height"},
    "Z10": {"realm": "atm", "var": "Z3", "plev": 1000.0, "obs": ("ERA5", "zg"), "scale": 1.0, "offset": 0.0, "units": "m", "long_name": "10-hPa geopotential height"},
    "T925": {"realm": "atm", "var": "T", "plev": 92500.0, "obs": ("ERA5", "ta"), "scale": 1.0, "offset": 0.0, "units": "K", "long_name": "925-hPa temperature"},
    "FLUT": {"realm": "atm", "obs": ("ERA5", "rlut"), "scale": 1.0, "offset": 0.0, "units": "W/m2", "long_name": "Outgoing longwave radiation"},
    "SOILWATER_10CM": {"realm": "lnd", "obs": None, "scale": 1.0, "offset": 0.0, "units": "kg/m2", "long_name": "Top 10-cm soil water"},
    "H2OSNO": {"realm": "lnd", "obs": None, "scale": 1.0, "offset": 0.0, "units": "mm", "long_name": "Snow water equivalent"},
    "TWS": {"realm": "lnd", "obs": None, "scale": 1.0, "offset": 0.0, "units": "mm", "long_name": "Terrestrial water storage"},
    "TSA": {"realm": "lnd", "obs": None, "scale": 1.0, "offset": -273.15, "units": "degC", "long_name": "Land 2-m air temperature"},
    "SICE": {"realm": "ice", "var": "timeMonthly_avg_iceAreaCell", "obs": ("HadISST2", "sic"), "scale": 1.0, "offset": 0.0, "units": "fraction", "long_name": "Sea-ice concentration"},
    # MPAS-O writes the 0-700 m heat content per native cell (J); "per_native_cell" divides it by the
    # native cell area remapped to the 1x1 grid (mpas_cell_area) to give J m-2, as EN4.
    "OHC700": {"realm": "ocn", "var": "timeMonthly_avg_oceanHeatContentSfcTo700m", "per_native_cell": True,
               "obs": ("EN4.2.2", "OHC700"), "scale": 1e-9, "offset": 0.0, "units": "GJ/m2",
               "long_name": "0-700 m ocean heat content"},
}

# Last year of each e3sm_diags obs time series. Past it, obs_product() switches to the
# archive products below (one file per year on the 1x1 grid, updated to the latest month).
OBS_LAST_YEAR = {"HadISST2": 2022, "ERA5": 2019, "GPCP_v2.3": 2017, "EN4.2.2": 0}   # EN4: archive only
OBS_ARCHIVE_ROOT = Path("/global/cfs/cdirs/e3sm/zhan391/data")
OBS_ARCHIVE = {
    "TS": ("OISSTv2.1", "NOAA_OISST/1x1/monthly/SST.monthly.*.nc", "SST"),       # 1981-09 on, degC
    "PRECT": ("GPCP_v2.3", "GPCPv2.3/1x1/monthly/PRECT.monthly.*.nc", "PRECT"),   # mm/day
    "SICE": ("NSIDCv6", "NOAA_NSIDC/1x1/monthly/ICEFRAC.monthly.*.nc", "ICEFRAC"),  # fraction
    "FLUT": ("NOAA_CDR_OLR", "NOAA_CDR/1x1/monthly/FLUT.monthly.*.nc", "FLUT"),     # 1979 on, W m-2
    "OHC700": ("EN4.2.2", "ENv4.2.2/1x1/monthly/OHC700.monthly.*.nc", "OHC700"),   # 1979 on, J m-2
}
# Alternative archive products, selected with obs_anomaly(field, product=...).
OBS_ALT_PRODUCTS = {
    "TS": {
        "HadISST1.1": ("HadSST/1x1/monthly/TS.monthly.*.nc", "TS"),              # 1870 on, K
        "OISSTv2.1": ("NOAA_OISST/1x1/monthly/SST.monthly.*.nc", "SST"),         # 1981-09 on, degC
        "ERSSTv5": ("ERSSTv5/1x1/monthly/sst.monthly.*.nc", "sst"),              # 1854 on, degC
        "ESA_CCI_SST": ("ESA_CCI_SST/1x1/monthly/sst.monthly.*.nc", "sst"),      # 1980 on, K
    },
    "FLUT": {
        "NOAA_CDR_OLR": ("NOAA_CDR/1x1/monthly/FLUT.monthly.*.nc", "FLUT"),      # NOAA CDR OLR, 1979 on, W m-2
    },
    "TREFHT": {
        "ERA5_Land": ("ERA5_LND/1x1/monthly/t2m_*.nc", "t2m"),                    # land only, 1950 on, K
    },
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


# Figures and tables of the case-study notebooks: OUTPUT_ROOT/<notebook>/{figures,data}/.
# Reusable caches (observation anomalies, lead-time climatologies) stay under diag_dir().
OUTPUT_ROOT = Path("/global/cfs/cdirs/e3sm/www/zhan391/e3sm-rtenso2026")


def output_dir(notebook, kind, *parts):
    """OUTPUT_ROOT/<notebook>/<kind>/<parts...>, created if missing (kind: "figures" or "data")."""
    path = OUTPUT_ROOT.joinpath(notebook, kind, *parts)
    path.mkdir(parents=True, exist_ok=True)
    return path


def figure_dir(notebook):
    return output_dir(notebook, "figures")


def save_figure(fig, notebook, *parts, dpi=150, ext="png"):
    path = figure_dir(notebook) / figure_filename(*parts, ext=ext)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    return path


def post_root(case_prefix):
    return RTENSO2026_POST_ROOT if "_RTENSO2026_" in case_prefix else env_paths.raw_model_root()


def case_dir(case_prefix, year, init_month):
    return post_root(case_prefix) / f"{case_prefix}_{year}{init_month:02d}0100"


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


@lru_cache(maxsize=4)
def mpas_cell_area(map_file):
    """Native MPAS cell area (m2) averaged onto the 1x1 grid with the file's own remap weights.

    An extensive per-cell field (e.g. heat content in J) remapped with "aave" is the area-weighted
    mean of native-cell totals, so dividing it by this mean native area gives a per-m2 value.
    """
    with xr.open_dataset(map_file) as m:
        row, col, S = m["row"].values - 1, m["col"].values - 1, m["S"].values
        area_a, frac_b = m["area_a"].values, m["frac_b"].values
        nlat, nlon = (int(n) for n in m["dst_grid_dims"].values[::-1])
    total = np.bincount(row, weights=S * area_a[col], minlength=frac_b.size)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(frac_b > 0, total / frac_b, np.nan) * EARTH_RADIUS**2
    grid = model_grid()
    return xr.DataArray(mean.reshape(nlat, nlon), coords={"lat": grid["lat"], "lon": grid["lon"]}, dims=("lat", "lon"))


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
        if info.get("per_native_cell"):
            da = da / mpas_cell_area(ds.attrs["map_file"]).values
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


def obs_monthly_climatology(field, clim_years=CLIM_YEARS, mode="auto", product=None):
    da = obs_monthly(field, mode=mode, product=product)
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


def obs_product(field, years=OBS_YEARS, product=None):
    """(product, file variable, archive glob or None) covering ``years`` for ``field``.

    ``product`` picks an entry of ``OBS_ALT_PRODUCTS``; None (or the default product) keeps
    the e3sm_diags default with its archive fallback.
    """
    default, obs_field = FIELDS[field]["obs"]
    if product is not None and product != default:
        if product not in OBS_ALT_PRODUCTS.get(field, {}):
            raise ValueError(f"unknown {field} product {product}; "
                             f"choose from {[default] + sorted(OBS_ALT_PRODUCTS.get(field, {}))}")
        pattern, obs_field = OBS_ALT_PRODUCTS[field][product]
        return product, obs_field, str(OBS_ARCHIVE_ROOT / pattern)
    product = default
    if years[1] <= OBS_LAST_YEAR.get(product, years[1]):
        return product, obs_field, None
    if field not in OBS_ARCHIVE:
        raise ValueError(f"no {field} observations after {OBS_LAST_YEAR[product]} ({product}); "
                         f"archive products cover {sorted(OBS_ARCHIVE)}")
    product, pattern, obs_field = OBS_ARCHIVE[field]
    return product, obs_field, str(OBS_ARCHIVE_ROOT / pattern)


def _source_stamp(pattern):
    """Latest modification time (s) of an archive product's files; 0 for e3sm_diags series.

    Part of the observation cache keys, so the caches rebuild when a product is extended
    (e.g. a new month appended to the current year's file).
    """
    files = glob.glob(pattern) if pattern else []
    return int(max(Path(f).stat().st_mtime for f in files)) if files else 0


def _archive_preproc(ds, **kwargs):
    """ESP-Lab monthly preprocessing that places each month by its time stamp.

    ESP-Lab takes the month from the lower time bound, but the GPCP v2.3 archive's bounds run
    from the 17th of the previous month, which would shift every month back by one. The archive
    time stamps (month start for OISST/GPCP, mid-month for NSIDC/ERSST/ESA CCI) name the right month.
    """
    if "valid_time" in ds.dims:   # ERA5(-Land) CDS files
        ds = ds.rename(valid_time="time")
    ds = ds.drop_vars([v for v in ("time_bnds", "time_bounds", "expver", "number") if v in ds.variables])
    ds["time"].attrs.pop("bounds", None)
    return dao.preprocessor_monthly(ds, **kwargs)


def obs_monthly(field, years=(CLIM_YEARS[0], OBS_YEARS[1]), mode="auto", product=None):
    """Observed monthly values on the model 1x1 grid (time, lat, lon) in display units."""
    product, obs_field, pattern = obs_product(field, years, product)
    path = diag_dir("observations", "case_study", FIELDS[field]["realm"], field,
                    f"{product}_{field}_monthly_y{years[0]}-{years[1]}_{GRID_TAG}.nc")
    expected = {"cache_version": OBS_CACHE_VERSION, "product": product, "start_year": years[0],
                "end_year": years[1], "units": FIELDS[field]["units"], "source_mtime": _source_stamp(pattern)}
    if mode != "rebuild" and _cache_ok(path, expected):
        return load_netcdf(path)[field]
    if mode == "require":
        raise FileNotFoundError(f"missing or stale observation cache {path}")

    # Chunk in time so 0.25-degree ERA5 is regridded a year at a time on shared login nodes.
    if pattern is None:
        ds = dao.get_monthly_data(obs_dir=str(env_paths.obs_root()), product=product, field=obs_field,
                                  start_year=str(years[0]), end_year=str(years[1]), chunks={"time": 12})
    else:
        ds = dao.get_monthly_data_from_pattern(pattern, field=obs_field, start_year=str(years[0]),
                                               end_year=str(years[1]), chunks={"time": 12},
                                               preproc=_archive_preproc)
    da = ds[list(ds.data_vars)[0]]
    if "plev" in da.dims:
        plev = da["plev"]
        target = FIELDS[field]["plev"] if float(plev.max()) > 2000 else FIELDS[field]["plev"] / 100.0
        da = da.sel(plev=target)
    da = _obs_to_display_units(da, field)
    if da.sizes["lat"] != 180 or da.sizes["lon"] != 360:
        src = da.to_dataset(name=field)
        da = make_regridder(src, model_grid(), method="conservative")(src)[field]
    # Load so a fresh build matches a cache read: dask chunks (one year) carried into the lead-drift
    # climatology are shorter than the ONI rolling window
    da = da.drop_vars([c for c in da.coords if c not in ("time", "lat", "lon")]).load()
    if field in ("TS", "SICE"):  # SST/sea ice are filled or missing over land; keep ocean points (EAM LANDFRAC < 0.5)
        da = da.where(land_fraction() < 0.5)
    da.attrs = {"units": FIELDS[field]["units"], "long_name": FIELDS[field]["long_name"], "product": product}
    out = da.to_dataset(name=field)
    out.attrs.update(expected)
    atomic_to_netcdf(out, path)
    return da


def obs_anomaly(field, years=OBS_YEARS, clim_years=CLIM_YEARS, mode="auto", product=None):
    """Observed monthly anomalies from the 1981-2010 monthly climatology, on the model grid.

    Past the e3sm_diags coverage the archive product supplies both the anomaly years and the
    climatology (OISST starts in 1981-09, so its Jan-Aug climatology is 1982-2010). ``product``
    selects an ``OBS_ALT_PRODUCTS`` entry instead (e.g. "ERSSTv5" or "ESA_CCI_SST" for TS).
    """
    span = (min(clim_years[0], years[0]), max(clim_years[1], years[1]))
    product, _, pattern = obs_product(field, span, product)
    path = diag_dir("observations", "case_study", FIELDS[field]["realm"], field,
                    f"{product}_{field}_anomaly_y{years[0]}-{years[1]}_clim_{clim_years[0]}_{clim_years[1]}_{GRID_TAG}.nc")
    expected = {"cache_version": OBS_CACHE_VERSION, "product": product, "clim_start_year": clim_years[0],
                "clim_end_year": clim_years[1], "units": FIELDS[field]["units"], "source_mtime": _source_stamp(pattern)}
    if mode != "rebuild" and _cache_ok(path, expected):
        return load_netcdf(path)[field]
    if mode == "require":
        raise FileNotFoundError(f"missing or stale observation cache {path}")

    da = obs_monthly(field, years=span, mode=mode, product=product)
    clim = da.sel(time=slice(f"{clim_years[0]}-01", f"{clim_years[1]}-12")).groupby("time.month").mean("time")
    # Load so a fresh build matches a cache read: dask chunks shorter than a rolling window (ONI) fail
    anom = (da.sel(time=slice(f"{years[0]}-01", f"{years[1]}-12")).groupby("time.month") - clim).drop_vars("month").load()
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
    elif field == "OHC700" and units in ("j m-2", "j/m2", "j m**-2"):
        da = da * FIELDS[field]["scale"]
    elif field == "SICE" and (units in ("%", "percent") or float(da.max()) > 1.5):
        da = da / 100.0  # HadISST2 sic is in percent although its units say "Fraction"
    da.attrs["units"] = FIELDS[field]["units"]
    return da


def match_obs(model, obs):
    """Observed values at each model lead's target month (dims follow ``model.time``).

    Target months past the end of the observations are NaN.
    """
    target = model["time"]
    picked = obs.reindex(time=[cftime.DatetimeNoLeap(t.year, t.month, 15) for t in target.values])
    return picked.rename(time="lead").assign_coords(lead=model["lead"].values, time=("lead", target.values))


# -----------------------------------------------------------------------------
# Cache builds (0_build_caches.ipynb)
# -----------------------------------------------------------------------------
CACHE_FIELDS = ["TS", "PRECT", "TREFHT", "PSL", "U850"]


def build_caches(fields=CACHE_FIELDS, mode="auto", verbose=True):
    """Build the Reanalysis lead-time climatologies (every ``REFERENCE`` start month) and the
    observation anomalies of ``fields`` ahead of the notebooks.

    ``mode``: "auto" reuses valid caches, "rebuild" recomputes them. Returns one row per cache
    built or checked: {"field", "cache", "seconds"}.
    """
    if mode not in ("auto", "rebuild"):
        raise ValueError(f"mode must be 'auto' or 'rebuild', not {mode!r}")
    rows = []

    def timed(field, cache, build):
        t0 = time.time()
        build()
        rows.append({"field": field, "cache": cache, "seconds": round(time.time() - t0, 1)})
        if verbose:
            print(f"{field} {cache} ready ({rows[-1]['seconds']:.0f} s)", flush=True)

    for field in fields:
        for init_month in REFERENCE["init_months"]:
            timed(field, f"init{init_month:02d} climatology",
                  lambda: leadtime_climatology(field, init_month, mode=mode))
        if FIELDS[field]["obs"] is not None:
            timed(field, "observation anomalies", lambda: obs_anomaly(field, mode=mode))
    return rows


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


def nino_indices(anom, names):
    """Nino box means of an SST anomaly, plus "ONI": the 3-month centred running mean
    of Nino3.4 along lead (hindcasts) or time (observations), as in ESP-Lab sst_index."""
    out = {name: box_mean(anom, NINO_REGIONS[name]) for name in names if name != "ONI"}
    if "ONI" in names:
        nino34 = out["Nino3.4"] if "Nino3.4" in out else box_mean(anom, NINO_REGIONS["Nino3.4"])
        dim = "lead" if "lead" in nino34.dims else "time"
        out["ONI"] = nino34.rolling({dim: 3}, center=True, min_periods=1).mean().assign_attrs(
            long_name="Oceanic Nino Index (3-month running mean of Nino3.4 anomalies)")
    return xr.Dataset({name: out[name] for name in names})


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


def seasonal_leads(da, nseason=8, complete=False):
    """Non-overlapping 3-month means over leads (season 1 = leads 0-2), keeping the centre-month time.

    ``complete=True`` leaves a season NaN unless each of its months has data somewhere (e.g.
    observations that end inside the season).
    """
    seasons = []
    for s in range(nseason):
        chunk = da.isel(lead=slice(3 * s, 3 * s + 3))
        if chunk.sizes["lead"] < 3:
            break
        mean = chunk.mean("lead")
        if complete and not bool(chunk.notnull().any([d for d in chunk.dims if d != "lead"]).all()):
            mean = mean.where(False)
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
    """Mean over ``months`` [(year, month), ...]; NaN unless every month is observed."""
    times = [cftime.DatetimeNoLeap(y, m, 15) for y, m in months]
    return obs.reindex(time=times).mean("time", skipna=False)


EARTH_RADIUS = 6.371e6


def velocity_potential(u, v):
    """Velocity potential (m2/s) of the wind (u, v) on a regular global lat-lon grid.

    Solves lap(chi) = div(u, v) for each zonal wavenumber: Fourier in longitude and a
    finite-volume Laplacian in latitude with no flux through the poles. chi has zero global
    mean. Linear, so the potential of anomalous winds is the velocity potential anomaly.
    """
    u, v = u.transpose(..., "lat", "lon"), v.transpose(..., "lat", "lon")
    lat = np.deg2rad(u["lat"].values)
    nlat, nlon = lat.size, u.sizes["lon"]
    dlat, a = lat[1] - lat[0], EARTH_RADIUS
    coslat = np.cos(lat)
    faces = np.cos(np.concatenate([[lat[0] - dlat / 2], (lat[:-1] + lat[1:]) / 2, [lat[-1] + dlat / 2]]))
    faces[[0, -1]] = 0.0   # no flux through the poles

    m = np.arange(nlon // 2 + 1)
    uh = np.fft.rfft(u.values, axis=-1)
    vh = np.fft.rfft(v.values * coslat[:, None], axis=-1)
    vface = np.zeros(vh.shape[:-2] + (nlat + 1, vh.shape[-1]), dtype=complex)   # v cos(lat) on the faces
    vface[..., 1:-1, :] = 0.5 * (vh[..., 1:, :] + vh[..., :-1, :])
    div = (1j * m * uh + (vface[..., 1:, :] - vface[..., :-1, :]) / dlat) / (a * coslat[:, None])

    lead = div.shape[:-2]
    rhs = div.reshape(-1, nlat, m.size)
    chi = np.empty_like(rhs)
    upper, lower = faces[1:-1] / (a**2 * dlat**2 * coslat[:-1]), faces[1:-1] / (a**2 * dlat**2 * coslat[1:])
    for k in m:
        lap = np.zeros((nlat, nlat))
        idx = np.arange(nlat)
        lap[idx, idx] = -(faces[:-1] + faces[1:]) / (a**2 * dlat**2 * coslat) - k**2 / (a * coslat) ** 2
        lap[idx[:-1], idx[1:]] = upper
        lap[idx[1:], idx[:-1]] = lower
        b = rhs[:, :, k].T
        chi[:, :, k] = (np.linalg.lstsq(lap, b, rcond=None)[0] if k == 0 else np.linalg.solve(lap, b)).T
    out = np.fft.irfft(chi.reshape(lead + (nlat, m.size)), n=nlon, axis=-1)
    out = xr.DataArray(out, coords=u.coords, dims=u.dims)
    w = area_weights(u["lat"])
    out = out - out.weighted(w).mean(("lat", "lon"))
    return out.assign_attrs(units="m2/s", long_name="Velocity potential")


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
    "hindcast_anomaly", "obs_monthly", "obs_anomaly", "match_obs", "box_mean", "nino_indices", "pattern_acc", "pattern_rmse",
    "seasonal_leads", "target_mean", "obs_target_mean", "hemispheric_area", "target_season_label",
    "welch_t", "region_subset", "area_weights",
    "figure_dir", "save_figure", "diag_dir",
]
