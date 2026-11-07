"""Stratosphere and polar-vortex helpers for ``9_ssw_polar_vortex.ipynb``.

Every product is a zonal mean on 1-degree latitude bands (centres -19.5..89.5, kept from the
QBO band to the pole) and pressure levels (hPa):

* **Model, daily**: the native ``eam.h1`` daily means (not post-processed) of
  U/V/T at 22 levels (1000-1 hPa) and Z at 16 levels (1000-10 hPa); each native column goes
  to the 1-degree band holding its centre, weighted by its area. ``VT`` is the eddy heat flux
  [v*T*] of the daily means within each band (no sub-daily transients). Cached per case and
  start under ``S2D_DIAG_ROOT/<source>/strat/``.
* **Model, monthly**: the post-processed 1x1 monthly U, T, Z3 (19 levels, 1000-1 hPa). The
  Reanalysis (BruteForce) 1981-2010 May/Nov starts give the lead-time climatology that removes
  drift; Feb/Aug starts use the "lead_drift" method of ``rtenso_case`` (observed target-month
  climatology plus the May/Nov mean drift at the same lead).
* **ERA5**: ``ERA5/1x1/daily_zonal`` (ERA5/scripts/daily_zonal; zonal means of the 00Z and
  12Z analyses, 37 levels).

SSW central dates follow Charlton and Polvani (2007) as applied by Butler et al. (2017): the
daily zonal-mean U at 60N, 10 hPa turns easterly between 1 Nov and 31 Mar; events are
separated by at least ``sep_days`` consecutive westerly days; and the wind must return to
westerly for ``recover_days`` consecutive days before 30 Apr (else it is the final warming).
The final warming date (Butler et al. 2019) is the spring onset of easterlies after which the
wind does not stay westerly for ``recover_days`` consecutive days before 31 Jul.

Data roots (raw simulations, the ERA5 archive) are arguments, set in the notebook.
"""

from __future__ import annotations

import glob
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import cftime
import netCDF4
import numpy as np
import pandas as pd
import xarray as xr
from scipy import sparse

from esp_lab.utils.netcdf_utils import atomic_to_netcdf, load_netcdf

from util import rtenso_case as rc

CACHE_VERSION = 1
LAT = np.arange(-89.5, 90.0, 1.0)
LAT_MIN = -19.5                      # zonal means kept from 19.5S (QBO band) to the pole
KEEP = LAT >= LAT_MIN
DAILY_LEVELS = (1000, 975, 950, 925, 900, 850, 800, 700, 600, 500, 400, 300, 250, 200, 150, 100, 50, 25, 10, 5, 2, 1)
Z_LEVELS = (1000, 975, 950, 925, 900, 850, 800, 700, 600, 500, 400, 300, 200, 100, 50, 10)
DAILY_FIELDS = ("U", "T", "Z", "VT")
MONTHLY_VARS = {"U": "U", "T": "T", "Z": "Z3"}
ERA5_FIELDS = ("U", "T", "Z", "VT")

# Index bands (centres of the 1-degree bands inside each range)
U60 = (59.0, 61.0)                   # SSW wind at 60N: the 59.5N and 60.5N bands
PCAP = (60.0, 90.0)                  # polar cap
EQ = (-5.0, 5.0)                     # QBO
VT_BAND = (45.0, 75.0)               # 100-hPa eddy heat flux

UNITS = {"U": "m/s", "T": "K", "Z": "m", "VT": "K m/s", "PSL": "hPa"}
LONG_NAMES = {"U": "zonal-mean zonal wind", "T": "zonal-mean temperature", "Z": "zonal-mean geopotential height",
              "VT": "eddy heat flux [v*T*]", "PSL": "zonal-mean sea-level pressure"}


# -----------------------------------------------------------------------------
# Generic helpers
# -----------------------------------------------------------------------------
def lat_mean(da, band):
    """cos(lat)-weighted mean over the 1-degree bands whose centres lie in ``band``."""
    sub = da.sel(lat=slice(*band))
    return sub.weighted(np.cos(np.deg2rad(sub["lat"]))).mean("lat")


def interp_logp(da, plev):
    """Linear interpolation in log(pressure) to ``plev`` (hPa), dimension ``plev``."""
    src = da.assign_coords(logp=("plev", np.log(da["plev"].values))).swap_dims(plev="logp").drop_vars("plev")
    out = src.interp(logp=np.log(np.asarray(plev, dtype=float)))
    return out.assign_coords(plev=("logp", np.asarray(plev, dtype=float))).swap_dims(logp="plev").drop_vars("logp")


def to_datetime64(times):
    """cftime (any calendar) or datetime64 -> datetime64[ns] (no-leap dates are valid dates)."""
    times = np.atleast_1d(times)
    if np.issubdtype(np.asarray(times).dtype, np.datetime64):
        return np.asarray(times, dtype="datetime64[ns]")
    return np.array([np.datetime64(f"{t.year:04d}-{t.month:02d}-{t.day:02d}") for t in times], dtype="datetime64[ns]")


def noleap_doy(times):
    """Day of a 365-day year (0-364); 29 Feb takes the index of 28 Feb."""
    t = pd.DatetimeIndex(to_datetime64(times))
    cum = np.cumsum([0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30])
    day = np.where((t.month == 2) & (t.day == 29), 28, t.day)
    return cum[t.month - 1] + day - 1


def _cache_ok(path, expected):
    if not Path(path).is_file():
        return False
    with xr.open_dataset(path) as ds:
        return all(str(ds.attrs.get(k)) == str(v) for k, v in expected.items())


def _filled(x):
    x = np.ma.filled(x.astype(np.float64), np.nan) if np.ma.isMaskedArray(x) else np.asarray(x, dtype=np.float64)
    x[np.abs(x) > 1e30] = np.nan
    return x


# -----------------------------------------------------------------------------
# Model: daily zonal means from the native eam.h1 files
# -----------------------------------------------------------------------------
def h1_files(simulation_root, case_prefix, year, init_month, member):
    hist = Path(simulation_root) / f"{case_prefix}_{year}{init_month:02d}0100" / f"EN{member:02d}" / "archive" / "atm" / "hist"
    return sorted(hist.glob("*.eam.h1.*.nc"))


def band_matrix(lat, area):
    """Sparse (n_kept_bands, ncol) area weights averaging native columns into 1-degree bands."""
    idx = np.clip(np.floor(np.asarray(lat) + 90.0).astype(int), 0, 179)
    rows = -np.ones(180, dtype=int)
    rows[KEEP] = np.arange(KEEP.sum())
    r = rows[idx]
    ok = r >= 0
    W = sparse.lil_matrix(sparse.csr_matrix((np.asarray(area)[ok], (r[ok], np.flatnonzero(ok))),
                                            shape=(int(KEEP.sum()), len(lat))))
    total = np.asarray(W.sum(1)).ravel()
    full = np.flatnonzero(total > 0)
    for i in np.flatnonzero(total == 0):      # e.g. 89.5N: no ne30pg2 column centre in 89-90N
        W[i] = W[full[np.argmin(np.abs(full - i))]]
    W = W.tocsr()
    return sparse.diags(1.0 / np.asarray(W.sum(1)).ravel()) @ W


def _read_member_daily(args):
    """One member's daily zonal means: (dates, {field: (time, plev, lat)}, psl (time, lat))."""
    files, label = args
    W, dates, out, psl = None, [], {f: [] for f in DAILY_FIELDS}, []
    for path in files:
        with netCDF4.Dataset(path) as nc:
            if W is None:
                W = band_matrix(nc["lat"][:], nc["area"][:])
            bnds = np.asarray(nc["time_bnds"][:])
            keep = np.flatnonzero(bnds[:, 1] > bnds[:, 0])          # drop the zero-length initial record
            if keep.size == 0:
                continue
            tv = nc["time"]
            dates.extend(to_datetime64(cftime.num2date(bnds[keep, 0], tv.units, getattr(tv, "calendar", "noleap"))))

            def native(name):
                return _filled(nc[name][keep])                      # (time, ncol)

            def zm(x):
                return (W @ x.T).T                                   # (time, band)

            nt, nb = keep.size, int(KEEP.sum())
            blk = {f: np.full((nt, len(DAILY_LEVELS), nb), np.nan) for f in DAILY_FIELDS}
            for k, p in enumerate(DAILY_LEVELS):
                u, v, t = native(f"U{p:03d}"), native(f"V{p:03d}"), native(f"T{p:03d}")
                blk["U"][:, k], blk["T"][:, k] = zm(u), zm(t)
                blk["VT"][:, k] = zm(v * t) - zm(v) * zm(t)
                if p in Z_LEVELS:
                    blk["Z"][:, k] = zm(native(f"Z{p:03d}"))
            for f in DAILY_FIELDS:
                out[f].append(blk[f].astype(np.float32))
            psl.append((zm(native("PSL")) / 100.0).astype(np.float32))
    if not dates:
        raise FileNotFoundError(f"{label}: no daily records")
    dates = np.asarray(dates, dtype="datetime64[ns]")
    order = np.argsort(dates, kind="stable")
    dates, first = np.unique(dates[order], return_index=True)
    sel = order[first]
    fields = {f: np.concatenate(out[f])[sel] for f in DAILY_FIELDS}
    gaps = np.diff(dates).astype("timedelta64[D]").astype(int)
    if (gaps != 1).any():
        raise ValueError(f"{label}: daily records are not consecutive")
    return dates, fields, np.concatenate(psl)[sel]


def model_daily_zonal(simulation_root, case_prefix, source, year, init_month, members=range(rc.NENS),
                      mode="auto", workers=8, verbose=True):
    """Daily zonal means of one case and start: Dataset U, T, Z, VT (member, time, plev, lat), PSL (member, time, lat).

    Reads the native ``eam.h1`` files under ``simulation_root`` (one process per member) and
    caches the result; None when no member has files.
    """
    path = rc.diag_dir(source, "strat", "daily_zonal",
                       f"{source}_{year}{init_month:02d}_daily_zonal_lat{LAT_MIN:g}-90_v{CACHE_VERSION}.nc")
    expected = {"cache_version": CACHE_VERSION, "case_prefix": case_prefix, "init_year": year, "init_month": init_month}
    if mode != "rebuild" and _cache_ok(path, expected):
        return load_netcdf(path)
    if mode == "require":
        raise FileNotFoundError(f"missing or stale daily zonal-mean cache {path}")
    jobs, found = [], []
    for m in members:
        files = h1_files(simulation_root, case_prefix, year, init_month, m)
        if files:
            jobs.append((files, f"{source} {year}{init_month:02d} EN{m:02d}"))
            found.append(m)
    if not jobs:
        return None
    if verbose:
        print(f"[{source} {year}-{init_month:02d}] reading {len(jobs)} members ({len(jobs[0][0])} h1 files each)", flush=True)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(_read_member_daily, jobs))
    n = min(r[0].size for r in results)                              # common length
    dates = results[0][0][:n]
    lat = LAT[KEEP]
    ds = xr.Dataset(coords={"member": found, "time": dates, "plev": np.asarray(DAILY_LEVELS, float), "lat": lat})
    for f in DAILY_FIELDS:
        ds[f] = (("member", "time", "plev", "lat"), np.stack([r[1][f][:n] for r in results]))
        ds[f].attrs = {"units": UNITS[f], "long_name": LONG_NAMES[f]}
    ds["PSL"] = (("member", "time", "lat"), np.stack([r[2][:n] for r in results]))
    ds["PSL"].attrs = {"units": UNITS["PSL"], "long_name": LONG_NAMES["PSL"]}
    ds["plev"].attrs = {"units": "hPa"}
    ds.attrs.update(expected, source=f"{simulation_root} eam.h1 daily means, native columns averaged into 1-degree bands")
    atomic_to_netcdf(ds, path)
    return ds


def pooled_daily_zonal(simulation_root, pool, year, init_month, mode="auto", workers=8, verbose=True):
    """Members of the ``pool`` cases {name: case_prefix} for one start, as one ``member`` dimension
    labelled "<case>:<member>". A case repeating an already pooled case member for member (the
    1 Feb 2026 clim and obs runs) is counted once."""
    parts = []
    for name, prefix in pool.items():
        source = prefix.split(f"{rc.MODEL_PREFIX}_")[-1]
        ds = model_daily_zonal(simulation_root, prefix, source, year, init_month, mode=mode, workers=workers, verbose=verbose)
        if ds is None:
            if verbose:
                print(f"[{name} {year}-{init_month:02d}] no daily output")
            continue
        u = ds["U"].sel(plev=10.0)
        if any(u.sizes == p["U"].sel(plev=10.0).sizes and np.array_equal(u.values, p["U"].sel(plev=10.0).values, equal_nan=True)
               for p in parts):
            if verbose:
                print(f"[{name} {year}-{init_month:02d}] repeats an already pooled case member for member; not pooled again")
            continue
        parts.append(ds.assign_coords(member=[f"{name}:{m}" for m in ds["member"].values]))
    if not parts:
        return None
    n = min(p.sizes["time"] for p in parts)
    return xr.concat([p.isel(time=slice(0, n)) for p in parts], dim="member")


# -----------------------------------------------------------------------------
# Model: monthly zonal means from the post-processed 1x1 files
# -----------------------------------------------------------------------------
def _read_member_monthly(args):
    case_prefix, year, init_month, member, nlead = args
    tsdir = rc.case_dir(case_prefix, year, init_month) / f"EN{member:02d}" / "post" / "atm" / rc.GRID / "ts" / "monthly" / "2yr"
    out, plev = {}, None
    for name, var in MONTHLY_VARS.items():
        files = sorted(glob.glob(str(tsdir / f"{var}_{year}{init_month:02d}_*.nc")))
        if not files:
            return None
        with netCDF4.Dataset(files[0]) as nc:
            lat = np.asarray(nc["lat"][:])
            if not np.allclose(lat, LAT):
                raise ValueError(f"{files[0]}: unexpected latitudes")
            tv = nc["time"]
            b = np.asarray(nc[tv.bounds][:])
            first = cftime.num2date(0.5 * (b[0, 0] + b[0, 1]), tv.units, getattr(tv, "calendar", "noleap"))
            if b.shape[0] < nlead or (first.year, first.month) != (year, init_month):
                raise ValueError(f"{files[0]} does not start at {year}-{init_month:02d} with {nlead} months")
            j = np.flatnonzero(KEEP)
            x = _filled(nc[var][:nlead, :, j[0]:j[-1] + 1, :])
            with np.errstate(all="ignore"):
                out[name] = np.nanmean(x, axis=-1).astype(np.float32)   # (lead, plev, lat)
            p = np.asarray(nc["plev"][:], dtype=float) / 100.0
            if plev is not None and not np.allclose(p, plev):
                raise ValueError("pressure levels differ between U, T, Z3")
            plev = p
    return plev, out


def monthly_zonal(case_prefix, source, years, init_month, members=range(rc.NENS), nlead=rc.NLEAD,
                  mode="auto", workers=8, verbose=True):
    """Monthly zonal means U, T, Z (year, member, lead, plev, lat) of a case for ``years`` (cached).
    Missing members or years are NaN; None when nothing is found."""
    years = list(years)
    path = rc.diag_dir(source, "strat", "monthly_zonal",
                       f"{source}_init{init_month:02d}_monthly_zonal_y{years[0]}-{years[-1]}_m{len(list(members))}_l{nlead}_v{CACHE_VERSION}.nc")
    expected = {"cache_version": CACHE_VERSION, "case_prefix": case_prefix, "init_month": init_month,
                "years": ",".join(map(str, years)), "nlead": nlead}
    if mode != "rebuild" and _cache_ok(path, expected):
        return load_netcdf(path)
    if mode == "require":
        raise FileNotFoundError(f"missing or stale monthly zonal-mean cache {path}")
    members = list(members)
    jobs = [(case_prefix, y, init_month, m, nlead) for y in years for m in members]
    if verbose:
        print(f"[{source} init{init_month:02d}] reading {len(jobs)} member-years", flush=True)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(_read_member_monthly, jobs))
    got = [r for r in results if r is not None]
    if not got:
        return None
    plev = got[0][0]
    shape = (len(years), len(members), nlead, plev.size, int(KEEP.sum()))
    data = {f: np.full(shape, np.nan, dtype=np.float32) for f in MONTHLY_VARS}
    for (_, y, _, m, _), r in zip(jobs, results):
        if r is not None:
            for f in MONTHLY_VARS:
                data[f][years.index(y), members.index(m)] = r[1][f]
    ds = xr.Dataset({f: (("year", "member", "lead", "plev", "lat"), data[f]) for f in MONTHLY_VARS},
                    coords={"year": years, "member": members, "lead": np.arange(nlead), "plev": plev, "lat": LAT[KEEP]})
    for f in MONTHLY_VARS:
        ds[f].attrs = {"units": UNITS[f], "long_name": LONG_NAMES[f]}
    ds["plev"].attrs = {"units": "hPa"}
    n_found = int(np.isfinite(data["U"][..., 0, -1, 0]).sum())
    ds.attrs.update(expected, member_years_found=n_found)
    atomic_to_netcdf(ds, path)
    return ds


def target_months(init_month, nlead=rc.NLEAD):
    return [(init_month - 1 + lead) % 12 + 1 for lead in range(nlead)]


def lead_climatology(init_month, ref, obs_month_clim, method="lead_drift"):
    """Monthly lead-time climatology (lead, plev, lat) for any start month.

    ``ref``: {init_month: Reanalysis monthly_zonal Dataset} for the starts that have one (May,
    Nov). Other starts: ``obs_month_clim`` (month, plev, lat; ERA5 on the model levels) of each
    target month plus the mean Reanalysis drift at the same lead ("lead_drift"), or the
    Reanalysis value for the same target month at the nearest lead ("target_month").
    """
    clim = {m: ds.mean(("year", "member")) for m, ds in ref.items()}
    if init_month in clim:
        return clim[init_month].assign_attrs(climatology_method=f"Reanalysis init{init_month:02d} {rc.CLIM_YEARS[0]}-{rc.CLIM_YEARS[1]}")
    if method == "lead_drift":
        drift = sum(clim[m] - obs_month_clim.sel(month=xr.DataArray(target_months(m), dims="lead")).drop_vars("month")
                    for m in clim) / len(clim)
        out = obs_month_clim.sel(month=xr.DataArray(target_months(init_month), dims="lead")).drop_vars("month") + drift
        return out.assign_coords(lead=np.arange(out.sizes["lead"])).assign_attrs(
            climatology_method="ERA5 target-month climatology + Reanalysis init05/init11 mean drift at the same lead")
    leads = []
    for lead in range(rc.NLEAD):
        pairs = rc.nearest_reference_leads(init_month, lead, ref_months=tuple(clim))
        leads.append(sum(clim[r].sel(lead=rl, drop=True) for r, rl in pairs) / len(pairs))
    return xr.concat(leads, dim="lead").assign_coords(lead=np.arange(rc.NLEAD)).assign_attrs(
        climatology_method="target-month matched Reanalysis init05/init11, nearest lead")


def daily_from_lead(clim, year, init_month, dates):
    """A monthly lead climatology (lead, ...) interpolated linearly in time to daily ``dates``
    (anchored at the 15th of each target month; held constant beyond the first/last anchor)."""
    anchors = []
    for lead in range(clim.sizes["lead"]):
        m = init_month - 1 + lead
        anchors.append(np.datetime64(f"{year + m // 12:04d}-{m % 12 + 1:02d}-15"))
    c = clim.assign_coords(time=("lead", np.asarray(anchors, dtype="datetime64[ns]"))).swap_dims(lead="time").drop_vars("lead")
    out = c.interp(time=to_datetime64(dates))
    return out.ffill("time").bfill("time")


# -----------------------------------------------------------------------------
# ERA5 daily zonal means (ERA5/1x1/daily_zonal)
# -----------------------------------------------------------------------------
def era5_daily_zonal(era5_zonal_dir, fields=ERA5_FIELDS, years=None, plev=None):
    """ERA5 daily zonal means (time, plev, lat) from ``era5_zonal_dir`` (ERA5/1x1/daily_zonal),
    latitudes LAT_MIN..89.5, optionally restricted to ``years`` (first, last) and levels ``plev``."""
    out = {}
    for f in fields:
        files = sorted(Path(era5_zonal_dir).glob(f"ERA5.daily_zonal.en00.{f}.*.nc"))
        if years is not None:
            files = [p for p in files if years[0] <= int(p.name.split(".")[-2][:4]) <= years[1]]
        if not files:
            raise FileNotFoundError(f"no ERA5 daily zonal {f} files in {era5_zonal_dir}")
        da = xr.open_mfdataset(files, combine="by_coords", data_vars="minimal", coords="minimal", compat="override")[f]
        da = da.sel(lat=slice(LAT_MIN - 0.01, 90.0))
        if plev is not None:
            da = da.sel(plev=list(plev))
        out[f] = da.astype("float32").load()
        out[f].attrs = {"units": UNITS[f], "long_name": LONG_NAMES[f], "product": "ERA5 daily zonal mean (00Z, 12Z)"}
    return xr.Dataset(out)


def daily_climatology(da, years, nharm=4):
    """Smoothed daily climatology and standard deviation (doy 0-364, ...) over ``years`` (first, last),
    keeping the annual mean and ``nharm`` harmonics of the raw calendar-day statistics."""
    sub = da.sel(time=slice(f"{years[0]}-01-01", f"{years[1]}-12-31"))
    doy = xr.DataArray(noleap_doy(sub["time"].values), dims="time", name="doy")
    mean = sub.groupby(doy).mean("time")
    anom = sub.groupby(doy) - mean
    std = np.sqrt((anom**2).groupby(doy).mean("time"))

    def smooth(x):
        spec = np.fft.rfft(x.values, axis=x.get_axis_num("doy"))
        idx = [slice(None)] * spec.ndim
        idx[x.get_axis_num("doy")] = slice(nharm + 1, None)
        spec[tuple(idx)] = 0.0
        return x.copy(data=np.fft.irfft(spec, n=x.sizes["doy"], axis=x.get_axis_num("doy")))

    return smooth(mean), smooth(std)


def daily_anomaly(da, clim):
    """Anomaly of a daily series from a (doy, ...) climatology."""
    doy = xr.DataArray(noleap_doy(da["time"].values), dims="time", coords={"time": da["time"]})
    return da - clim.sel(doy=doy).drop_vars("doy")


def monthly_means(da):
    """Calendar-month means of a daily series, stamped at the 15th."""
    m = da.resample(time="MS").mean()
    return m.assign_coords(time=m["time"].values + np.timedelta64(14, "D"))


def obs_month_climatology(da, years):
    """(month, ...) climatology of calendar-month means over ``years``."""
    m = monthly_means(da).sel(time=slice(f"{years[0]}-01", f"{years[1]}-12"))
    return m.groupby("time.month").mean("time")


# -----------------------------------------------------------------------------
# SSW and final warming
# -----------------------------------------------------------------------------
def _runs(mask):
    """(start, length) of the runs of True in a boolean array."""
    m = np.concatenate([[False], np.asarray(mask, bool), [False]])
    d = np.diff(m.astype(int))
    starts, ends = np.flatnonzero(d == 1), np.flatnonzero(d == -1)
    return starts, ends - starts


def winter_of(date):
    """Winter labelled by the year of its November (Jan 2027 -> 2026)."""
    d = pd.Timestamp(date)
    return d.year if d.month >= 7 else d.year - 1


def covered_winters(times, season=((11, 1), (4, 30))):
    """Winters (November year) whose 1 Nov - 30 Apr is fully inside ``times``."""
    t = pd.DatetimeIndex(to_datetime64(times))
    out = []
    for w in range(t[0].year - 1, t[-1].year + 1):
        start = pd.Timestamp(w, *season[0])
        end = pd.Timestamp(w + 1, *season[1])
        if t[0] <= start and t[-1] >= end:
            out.append(w)
    return out


def ssw_central_dates(u, times, sep_days=20, recover_days=10, season=((11, 1), (3, 31)), recover_by=(4, 30)):
    """SSW central dates of a daily U(60N, 10 hPa) series (Charlton and Polvani 2007).

    Only winters fully covered from 1 Nov to 30 Apr are searched (``covered_winters``).
    """
    t = pd.DatetimeIndex(to_datetime64(times))
    u = np.asarray(u, dtype=float)
    west = u > 0
    events = []
    for w in covered_winters(t):
        start, stop = pd.Timestamp(w, *season[0]), pd.Timestamp(w + 1, *season[1])
        recover_end = t.get_loc(pd.Timestamp(w + 1, *recover_by))
        last = None
        for i in range(t.get_loc(start), t.get_loc(stop) + 1):
            if not (u[i] < 0 and i > 0 and u[i - 1] >= 0):
                continue
            # separation: a run of >= sep_days westerly days since the previous event
            if last is not None:
                _, lengths = _runs(west[last:i])
                if not (lengths >= sep_days).any():
                    continue
            # not a final warming: >= recover_days consecutive westerly days, ending by 30 Apr
            s, lengths = _runs(west[i:recover_end + 1])
            if not (lengths >= recover_days).any():
                continue
            events.append(t[i])
            last = i
    return events


def final_warming_date(u, times, spring_year, recover_days=10, search=((1, 1), (6, 30)), until=(7, 31)):
    """Final warming of the spring of ``spring_year`` (Butler et al. 2019): the first day in
    ``search`` on which U(60N, 10 hPa) turns easterly (or the first day of the series, if it starts
    easterly) after which it never stays westerly for ``recover_days`` consecutive days before
    ``until`` (or the end of the series). NaT when there is none, or when the series ends before
    30 Apr or starts after the search window. ``until`` ends before the autumn westerlies
    (E3SM's return in early August, ERA5's in late August)."""
    t = pd.DatetimeIndex(to_datetime64(times))
    u = np.asarray(u, dtype=float)
    lo, hi = pd.Timestamp(spring_year, *search[0]), pd.Timestamp(spring_year, *search[1])
    end = min(pd.Timestamp(spring_year, *until), t[-1])
    if t[-1] < pd.Timestamp(spring_year, 4, 30) or t[0] > hi:
        return pd.NaT
    i0, i1, ie = t.searchsorted(lo), t.searchsorted(min(hi, end), side="right") - 1, t.get_loc(end)
    west = u > 0
    for i in range(i0, i1 + 1):
        if not (u[i] < 0 and (i == 0 or u[i - 1] >= 0)):
            continue
        _, lengths = _runs(west[i:ie + 1])
        if not (lengths >= recover_days).any():
            return t[i]
    return pd.NaT


def ssw_table(u, times, label="", **kw):
    """DataFrame of SSW central dates of one series with its covered winters."""
    rows = [{"series": label, "winter": winter_of(d), "central_date": d} for d in ssw_central_dates(u, times, **kw)]
    return pd.DataFrame(rows, columns=["series", "winter", "central_date"])


def lag_composite(da, dates, lags):
    """(event, lag, ...) windows of a daily series around ``dates``; days outside the series are NaN."""
    t = pd.DatetimeIndex(to_datetime64(da["time"].values))
    out = []
    lags = np.asarray(lags)
    for d in dates:
        i = t.get_loc(pd.Timestamp(d))
        idx = i + lags
        ok = (idx >= 0) & (idx < t.size)
        win = da.isel(time=np.clip(idx, 0, t.size - 1)).where(xr.DataArray(ok, dims="time"))
        out.append(win.drop_vars("time").rename(time="lag").assign_coords(lag=lags))
    if not out:
        return None
    return xr.concat(out, dim="event")
