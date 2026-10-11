"""Scorecard metrics for EAM-DART analysis-period evaluation.

Every experiment is scored on the same 64 assimilation cycles per season
(S0: analyses 16 Dec 2011 00Z-31 Dec 18Z; S1: 16 May 2012 00Z-31 May 18Z):

* group A, upper-air T/U/V/Q: DART ensemble mean/sd (remapped to pressure
  levels) at the analysis times vs ERA5 6-hourly; ``analysis`` uses the DART
  posterior, ``prior`` the 6-h forecast DART started from;
* group B, surface state: ``ts/6hourly_da`` members at the analysis times vs
  ERA5 6-hourly (PSL also has a ``prior``: the 6-h forecast valid at that time);
* groups C/D, water cycle, fluxes and radiation: ``ts/6hourly`` forecast means
  aggregated to daily means vs daily references (GPCP, IMERG, NOAA OLR, or
  ERA5 averaged over its 00/06/12/18Z samples).

For each verification time and region the cosine-latitude weighted bias, MSE,
ensemble variance, and Gaussian CRPS are stored, on a mask shared by every
experiment and state of the season, so differences between experiments are
not caused by different sampling. Window scores and paired DART-vs-baseline
changes with block-bootstrap confidence intervals are derived from those
per-time scores.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy.special import erf

from util.bias_metrics import experiment_registry

N_CYCLES = 64
WINDOW_START = {"S0": "2011-12-16 00:00", "S1": "2012-05-16 00:00"}
REGIONS = {
    "Global": (-90, 90),
    "NH extratropics": (20, 90),
    "Tropics": (-20, 20),
    "SH extratropics": (-90, -20),
}
# analysis: posterior / 6hourly_da; prior: 6-h forecast valid at the analysis
# time; forecast: daily means of the 6-h forecast means (process variables).
STATES = ("analysis", "prior", "forecast")


@dataclass(frozen=True)
class Row:
    """One heat-map row: what to score and against which reference."""

    group: str
    name: str            # row label, e.g. "T850", "PRECT (GPCP)"
    variable: str        # model variable
    kind: str            # "dart" | "ts_6h" | "ts_daily"
    reference: str       # "era5_6h" | "era5_daily4" | "gpcp" | "imerg" | "noaa_olr"
    level: int | None = None        # hPa, DART rows
    reference_variable: str | None = None
    units: str = ""
    has_prior: bool = False


def default_rows():
    rows = [Row("A Upper-air state", f"{v}{p}", v, "dart", "era5_6h", level=p,
                reference_variable=f"{v}{p}", units={"T": "K", "U": "m s-1", "V": "m s-1", "Q": "g kg-1"}[v],
                has_prior=True)
            for v in ("T", "U", "V", "Q") for p in (850, 500, 200)]
    rows += [
        Row("B Surface state", "PSL", "PSL", "ts_6h", "era5_6h", reference_variable="PSL", units="hPa", has_prior=True),
        Row("B Surface state", "TS", "TS", "ts_6h", "era5_6h", reference_variable="TS", units="K"),
        Row("B Surface state", "TREFHT", "TREFHT", "ts_6h", "era5_6h", reference_variable="TREFHT", units="K"),
        Row("B Surface state", "U10", "U10", "ts_6h", "era5_6h", reference_variable="U10", units="m s-1"),
        Row("C Water cycle & fluxes", "PRECT (GPCP)", "PRECT", "ts_daily", "gpcp", units="mm day-1"),
        Row("C Water cycle & fluxes", "PRECT (IMERG)", "PRECT", "ts_daily", "imerg", units="mm day-1"),
        Row("C Water cycle & fluxes", "PRECT (ERA5)", "PRECT", "ts_daily", "era5_daily4", reference_variable="PRECT", units="mm day-1"),
        Row("C Water cycle & fluxes", "LHFLX", "LHFLX", "ts_daily", "era5_daily4", reference_variable="LHFLX", units="W m-2"),
        Row("C Water cycle & fluxes", "SHFLX", "SHFLX", "ts_daily", "era5_daily4", reference_variable="SHFLX", units="W m-2"),
        Row("D Radiation", "FLUT (NOAA)", "FLUT", "ts_daily", "noaa_olr", units="W m-2"),
        Row("D Radiation", "FLDS", "FLDS", "ts_daily", "era5_daily4", reference_variable="FLDS", units="W m-2"),
        Row("D Radiation", "FSDS", "FSDS", "ts_daily", "era5_daily4", reference_variable="FSDS", units="W m-2"),
    ]
    return rows


def land_rows():
    """Land-surface state and land-atmosphere exchange, scored over land only.

    No soil-moisture or snow series exist for every experiment (CTRL10/DART10
    ELM history output lacks them), so the land state is verified through
    what ERA5 6-hourly provides: skin and 2 m temperature (analysis) and the
    surface energy and water fluxes (daily means of the 6-h forecasts). FLNS
    is omitted: the ERA5 6-hourly FLNS/FLNSC files are mostly zeros.
    """
    era5 = dict(reference="era5_daily4", units="W m-2")
    return [
        Row("A Land surface state", "TS", "TS", "ts_6h", "era5_6h", reference_variable="TS", units="K"),
        Row("A Land surface state", "TREFHT", "TREFHT", "ts_6h", "era5_6h", reference_variable="TREFHT", units="K"),
        Row("B Surface fluxes", "LHFLX", "LHFLX", "ts_daily", reference_variable="LHFLX", **era5),
        Row("B Surface fluxes", "SHFLX", "SHFLX", "ts_daily", reference_variable="SHFLX", **era5),
        Row("B Surface fluxes", "FSNS", "FSNS", "ts_daily", reference_variable="FSNS", **era5),
        Row("C Forcing", "PRECT (GPCP)", "PRECT", "ts_daily", "gpcp", units="mm day-1"),
        Row("C Forcing", "PRECT (IMERG)", "PRECT", "ts_daily", "imerg", units="mm day-1"),
        Row("C Forcing", "PRECT (ERA5)", "PRECT", "ts_daily", "era5_daily4", reference_variable="PRECT", units="mm day-1"),
        Row("C Forcing", "FSDS", "FSDS", "ts_daily", reference_variable="FSDS", **era5),
        Row("C Forcing", "FLDS", "FLDS", "ts_daily", reference_variable="FLDS", **era5),
    ]


def land_mask_like(field, land_fraction, *, threshold=0.5, south=-60.0):
    """Boolean land mask on ``field``'s grid: land fraction >= ``threshold``, north of ``south``.

    ``south`` drops Antarctica so the ice sheet does not dominate SH land.
    """
    land = _standard_grid(land_fraction)
    if float(field.lon.min()) < 0 <= float(land.lon.min()):
        land = land.assign_coords(lon=((land.lon + 180) % 360) - 180).sortby("lon")
    land = land.interp(lat=field.lat, lon=field.lon, method="nearest")
    return (land >= threshold) & (field.lat >= south)


# ---------------------------------------------------------------- time axes
def analysis_times(season):
    return pd.date_range(WINDOW_START[season], periods=N_CYCLES, freq="6h")


def window_days(season):
    return pd.date_range(WINDOW_START[season], periods=N_CYCLES // 4, freq="D")


def _stamp(time):
    time = pd.Timestamp(time)
    return f"{time:%Y-%m-%d}-{time.hour * 3600 + time.minute * 60:05d}"


# ------------------------------------------------------------ unit handling
_SCALE = {"PSL": 0.01, "PRECT": 86400e3, "Q": 1e3}  # Pa->hPa, m/s->mm/day, kg/kg->g/kg


def _to_row_units(field, variable, units_hint=""):
    """Convert a model or reference field to the row units."""
    units = str(units_hint or field.attrs.get("units", "")).lower()
    if variable == "PRECT":
        if "mm" in units:
            return field
        if "kg" in units:                       # kg m-2 s-1 == mm/s
            return field * 86400.0
        return field * 86400e3                  # m/s
    return field * _SCALE.get(variable, 1.0)


def _standard_grid(field):
    field = field.rename({k: v for k, v in {"latitude": "lat", "longitude": "lon"}.items() if k in field.dims})
    if "level" in field.dims:
        field = field.squeeze("level", drop=True)
    return field.reset_coords(drop=True)


# ---------------------------------------------------------------- loaders
def load_dart(remap_root, experiment, row, times, *, state="analysis"):
    """Ensemble mean and sd of a remapped DART summary at ``times``."""
    product = "output" if state == "analysis" else "forecast"
    means, sds = [], []
    for time in times:
        stamp = _stamp(time)
        for kind, store in (("mean", means), ("sd", sds)):
            path = Path(remap_root) / experiment / f"{experiment}.{product}_{kind}.{stamp}.nc"
            with xr.open_dataset(path) as ds:
                field = ds[row.variable].sel(plev=row.level * 100.0).squeeze(drop=True).load()
            store.append(_standard_grid(field))
    mean = xr.concat(means, "time").assign_coords(time=times)
    sd = xr.concat(sds, "time").assign_coords(time=times)
    return _to_row_units(mean, row.variable), _to_row_units(sd, row.variable)


def _member_files(config, source, variable, year):
    return [config["path"] / f"EN{m:02d}" / "archive/post/atm/180x360_aave/ts" / source / f"{variable}.EN{m:02d}.{year}.nc"
            for m in range(1, config["nens"] + 1)]


def _open_members(config, source, variable, year, times=None):
    """Member fields stacked along ``ens`` (optionally only at ``times``), plus cycle starts."""
    files = _member_files(config, source, variable, year)
    missing = [f for f in files if not f.is_file()]
    if missing:
        raise FileNotFoundError(f"{len(missing)} member files missing, e.g. {missing[0]}")
    fields = []
    for path in files:
        with xr.open_dataset(path) as ds:
            field = ds[variable] if times is None else ds[variable].sel(time=times)
            fields.append(field.load())
            if path == files[0]:
                bounds = ds["time_bnds"] if times is None else ds["time_bnds"].sel(time=times)
                cycle_start = bounds.isel({bounds.dims[-1]: 0}).load()
    return xr.concat(fields, "ens"), cycle_start


def load_ts_6h(config, row, times, *, state="analysis"):
    """Members at analysis times: ``6hourly_da`` (analysis) or the 6-h forecast valid then (prior)."""
    year = times[0].year
    source = "6hourly_da" if state == "analysis" else "6hourly"
    field, _ = _open_members(config, source, row.variable, year, times)
    field = _to_row_units(field, row.variable)
    return _standard_grid(field.mean("ens")), _standard_grid(field.std("ens", ddof=1))


def load_ts_daily(config, row, days):
    """Daily means of the 6-h forecast means, grouped by the cycle-start day."""
    year = days[0].year
    field, start = _open_members(config, "6hourly", row.variable, year)
    day = start.dt.floor("D")
    keep = day.isin(days.values)
    field = field.isel(time=np.flatnonzero(keep.values))
    field = field.assign_coords(day=("time", day.values[keep.values]))
    counts = field.day.to_series().value_counts()
    if (counts != 4).any():
        raise ValueError(f"{row.variable}: expected 4 cycles per day, got {counts[counts != 4].to_dict()}")
    daily = field.groupby("day").mean("time").rename(day="time")
    daily = _to_row_units(daily, row.variable)
    return _standard_grid(daily.mean("ens")), _standard_grid(daily.std("ens", ddof=1))


def _era5_6h(reference_root, variable, times):
    root = Path(reference_root) / "ERA5/6hourly"
    years = sorted({t.year for t in times})
    fields = []
    for year in years:
        path = root / f"ERA5.6hourly.en00.{variable}.{year}01-{year}12.nc"
        with xr.open_dataset(path) as ds:
            fields.append(_standard_grid(ds[variable]).sel(time=[t for t in times if t.year == year]).load())
    field = xr.concat(fields, "time")
    return field, field.attrs.get("units", "")


def load_reference(reference_root, row, times_or_days):
    """Reference field on the row's verification times (6-h) or days."""
    root = Path(reference_root)
    if row.reference == "era5_6h":
        if row.variable == "U10":   # EAM U10 is the 10 m wind speed
            u, _ = _era5_6h(root, "U10", times_or_days)
            v, _ = _era5_6h(root, "V10", times_or_days)
            return np.hypot(u, v)
        field, units = _era5_6h(root, row.reference_variable, times_or_days)
        return _to_row_units(field, row.variable, units)
    if row.reference == "era5_daily4":
        samples = pd.DatetimeIndex([d + pd.Timedelta(hours=h) for d in times_or_days for h in (0, 6, 12, 18)])
        field, units = _era5_6h(root, row.reference_variable, samples)
        daily = field.groupby(field.time.dt.floor("D")).mean("time").rename(floor="time")
        return _to_row_units(daily.assign_coords(time=times_or_days.values), row.variable, units)
    product = {"gpcp": ("GPCP/daily", "PRECT.daily.{y}.nc", "PRECT"),
               "imerg": ("IMERG/daily", "PRECT.daily.{y}.nc", "PRECT"),
               "noaa_olr": ("NOAA_IOLR/daily", "FLUT.daily.{y}.nc", "FLUT")}[row.reference]
    year = times_or_days[0].year
    with xr.open_dataset(root / product[0] / product[1].format(y=year)) as ds:
        field = _standard_grid(ds[product[2]])
        field = field.assign_coords(time=field.time.dt.floor("D")).sel(time=times_or_days.values).load()
    return _to_row_units(field, row.variable, field.attrs.get("units", ""))


# ----------------------------------------------------------------- scoring
def gaussian_crps(mean, sd, obs):
    """CRPS of N(mean, sd^2) for ``obs``; reduces to |obs - mean| when sd == 0."""
    sd = np.maximum(sd, 1e-12)
    z = (obs - mean) / sd
    pdf = np.exp(-0.5 * z ** 2) / np.sqrt(2 * np.pi)
    cdf = 0.5 * (1 + erf(z / np.sqrt(2)))
    return sd * (z * (2 * cdf - 1) + 2 * pdf - 1 / np.sqrt(np.pi))


def per_time_scores(mean, sd, obs, mask):
    """Weighted regional bias/MSE/variance/CRPS at every verification time."""
    weights = np.cos(np.deg2rad(obs.lat)).broadcast_like(obs).where(mask)
    error = mean - obs
    terms = {"bias": error, "mse": error ** 2, "var": sd ** 2,
             "crps": xr.apply_ufunc(gaussian_crps, mean, sd, obs)}
    out = {}
    for region, (south, north) in REGIONS.items():
        band = (obs.lat >= south) & (obs.lat <= north)
        w = weights.where(band)
        total = w.sum(("lat", "lon"))
        out[region] = xr.Dataset({k: (v * w).sum(("lat", "lon")) / total for k, v in terms.items()})
    return xr.concat(list(out.values()), pd.Index(list(out), name="region"))


def score_row(row, season, experiments, *, data_root, remap_root, reference_root, states=STATES,
              land_fraction=None, log=print):
    """Per-time scores for every experiment/state of one row and season.

    With ``land_fraction`` (a 2-D land-fraction field) only land cells are
    scored; see :func:`land_mask_like`.

    Returns a Dataset with dims (experiment, state, time, region) and the
    ensemble size per experiment; states without data are NaN.
    """
    registry = experiment_registry(data_root, season, "da")
    times = analysis_times(season) if row.kind in ("dart", "ts_6h") else window_days(season)
    if row.kind == "ts_daily":
        states = ["forecast"]
    else:
        states = [s for s in states if s == "analysis" or (s == "prior" and row.has_prior)]
    reference = load_reference(reference_root, row, times)
    fields = {}
    for experiment in experiments:
        for state in states:
            if row.kind == "dart":
                fields[experiment, state] = load_dart(remap_root, experiment, row, times, state=state)
            elif row.kind == "ts_6h":
                fields[experiment, state] = load_ts_6h(registry[experiment], row, times, state=state)
            else:
                fields[experiment, state] = load_ts_daily(registry[experiment], row, times)
            log(f"  loaded {season} {row.name} {experiment} {state}")
    # one mask for every experiment and state: compare like with like
    mask = reference.notnull()
    if land_fraction is not None:
        mask = mask & land_mask_like(reference, land_fraction)
    for mean, sd in fields.values():
        mask = mask & mean.notnull() & sd.notnull()
    scores = []
    for experiment in experiments:
        per_state = [per_time_scores(*fields[experiment, s], reference, mask) for s in states]
        scores.append(xr.concat(per_state, pd.Index(states, name="state")))
    out = xr.concat(scores, pd.Index(list(experiments), name="experiment"))
    out = out.reindex(state=list(STATES))
    out["ensemble_members"] = ("experiment", [registry[e]["nens"] for e in experiments])
    out["valid_fraction"] = mask.mean(("lat", "lon")).mean("time")
    out.attrs.update(domain="land" if land_fraction is not None else "all", row=row.name, group=row.group, variable=row.variable, units=row.units,
                     reference=row.reference, season=season, kind=row.kind,
                     verification="6-hourly" if row.kind != "ts_daily" else "daily",
                     n_times=len(times))
    return out


def score_cache_path(cache_dir, season, row):
    safe = row.name.replace(" ", "_").replace("(", "").replace(")", "")
    return Path(cache_dir) / f"scorecard_{season}_{safe}.nc"


def compute_scorecard(rows, seasons, *, data_root, remap_root, reference_root, cache_dir,
                      force_compute=False, land_fraction=None, log=print):
    """Score every row/season (cached); return ``{(season, row_name): Dataset}`` and failures."""
    results, failures = {}, []
    for season, experiments in seasons.items():
        for row in rows:
            path = score_cache_path(cache_dir, season, row)
            if path.is_file() and not force_compute:
                cached = xr.load_dataset(path)
                if list(cached.experiment.values) == list(experiments):
                    results[season, row.name] = cached
                    log(f"[REUSE] {path.name}")
                    continue
            try:
                ds = score_row(row, season, experiments, data_root=data_root, remap_root=remap_root,
                               reference_root=reference_root, land_fraction=land_fraction, log=log)
            except Exception as error:
                failures.append((season, row.name, f"{type(error).__name__}: {error}"))
                log(f"[FAIL] {season} {row.name}: {type(error).__name__}: {error}")
                continue
            path.parent.mkdir(parents=True, exist_ok=True)
            tmp = path.with_suffix(".nc.tmp")
            ds.to_netcdf(tmp)
            os.replace(tmp, path)
            results[season, row.name] = ds
            log(f"[WRITE] {path.name}")
    return results, failures


# ------------------------------------------------------- window aggregation
def primary_state(ds):
    """State shown in the with/without-DA scorecard for this row."""
    return "forecast" if ds.attrs.get("kind") == "ts_daily" else "analysis"


def window_scores(ds):
    """Window-mean bias, RMSE, spread, spread-skill ratio, and CRPS."""
    n = ds["ensemble_members"]
    rmse = np.sqrt(ds["mse"].mean("time"))
    spread = np.sqrt(ds["var"].mean("time"))
    return xr.Dataset({
        "bias": ds["bias"].mean("time"),
        "rmse": rmse,
        "spread": spread,
        "ssr": np.sqrt((n + 1) / n) * spread / rmse,
        "crps": ds["crps"].mean("time"),
    })


def _block_indices(n_times, block, rng):
    starts = rng.integers(0, n_times - block + 1, size=int(np.ceil(n_times / block)))
    return np.concatenate([np.arange(s, s + block) for s in starts])[:n_times]


def relative_change(test, base, *, score="rmse", n_boot=1000, block=4, seed=0):
    """Percent change of ``score`` (rmse or crps) from ``base`` to ``test``.

    ``test`` and ``base`` are per-time score DataArrays over the same times
    (mse or crps). The 95% interval comes from a paired moving-block bootstrap
    over time (default 4 cycles = 1 day) to respect autocorrelation.
    """
    def change(t, b):
        # Times missing in the shared mask are missing for both, so nanmean
        # keeps the comparison paired.
        if score == "rmse":
            return 100 * (np.sqrt(np.nanmean(t, -1)) / np.sqrt(np.nanmean(b, -1)) - 1)
        return 100 * (np.nanmean(t, -1) / np.nanmean(b, -1) - 1)

    t = np.asarray(test); b = np.asarray(base)
    value = change(t, b)
    rng = np.random.default_rng(seed)
    n = t.shape[-1]
    samples = np.stack([change(t[..., idx], b[..., idx]) for idx in
                        (_block_indices(n, block, rng) for _ in range(n_boot))])
    lower, upper = np.nanpercentile(samples, [2.5, 97.5], axis=0)
    return value, lower, upper


def change_table(results, *, comparisons, score="rmse", n_boot=1000, block_cycles=4):
    """Tidy table of % changes for (season, test, base, test_state, base_state) comparisons."""
    rows = []
    for (season, row_name), ds in results.items():
        block = block_cycles if ds.attrs["verification"] == "6-hourly" else max(1, block_cycles // 4)
        field = "mse" if score == "rmse" else "crps"
        for comp in comparisons:
            c_season, test, base, test_state, base_state = comp
            if c_season != season or test not in ds.experiment or base not in ds.experiment:
                continue
            t = ds[field].sel(experiment=test, state=test_state)
            b = ds[field].sel(experiment=base, state=base_state)
            if bool(t.isnull().all()) or bool(b.isnull().all()):
                continue
            for region in ds.region.values:
                value, lower, upper = relative_change(t.sel(region=region).values, b.sel(region=region).values,
                                                      score=score, n_boot=n_boot, block=block)
                rows.append({"season": season, "group": ds.attrs["group"], "row": row_name, "region": region,
                             "test": test, "test_state": test_state, "base": base, "base_state": base_state,
                             "score": score, "change_pct": float(value), "ci_lower": float(lower),
                             "ci_upper": float(upper), "significant": bool(lower > 0 or upper < 0)})
    return pd.DataFrame(rows)


# ------------------------------------------------- normalised RMSE (seasons)
def reference_std(row, season, *, reference_root):
    """Regional std of the reference about its window-mean at each grid point.

    sqrt of the cos-latitude weighted mean over the window and region of
    (obs - time-mean obs)^2: the natural variability the window's errors are
    measured against. Uses the reference's own valid mask.
    """
    times = analysis_times(season) if row.kind in ("dart", "ts_6h") else window_days(season)
    obs = load_reference(reference_root, row, times)
    anomaly2 = (obs - obs.mean("time")) ** 2
    weights = np.cos(np.deg2rad(obs.lat)).broadcast_like(obs).where(obs.notnull())
    out = {}
    for region, (south, north) in REGIONS.items():
        w = weights.where((obs.lat >= south) & (obs.lat <= north))
        out[region] = float(np.sqrt((anomaly2 * w).sum() / w.sum()))
    return pd.Series(out, name="reference_std").rename_axis("region")


def compute_reference_std(rows, seasons, *, reference_root, cache_path, force_compute=False, log=print):
    """``reference_std`` for every row and season, cached as one CSV."""
    cache_path = Path(cache_path)
    if cache_path.is_file() and not force_compute:
        log(f"[REUSE] {cache_path.name}")
        return pd.read_csv(cache_path)
    frames = []
    for season in seasons:
        for row in rows:
            try:
                frame = reference_std(row, season, reference_root=reference_root).reset_index()
            except Exception as error:
                log(f"[FAIL] {season} {row.name}: {type(error).__name__}: {error}")
                continue
            frame.insert(0, "row", row.name); frame.insert(0, "season", season)
            frames.append(frame)
    table = pd.concat(frames, ignore_index=True)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(cache_path, index=False)
    log(f"[WRITE] {cache_path.name}")
    return table


def normalized_rmse_table(results, reference_stds, experiments):
    """Window RMSE / reference std of ``experiments[season]`` in its primary state.

    ``experiments`` maps season -> experiment. Columns match
    :func:`change_table` where they overlap, with the value in ``nrmse``.
    """
    stds = reference_stds.set_index(["season", "row", "region"])["reference_std"]
    rows = []
    for (season, row_name), ds in results.items():
        if season not in experiments or experiments[season] not in ds.experiment:
            continue
        state = primary_state(ds)
        rmse = np.sqrt(ds["mse"].sel(experiment=experiments[season], state=state).mean("time"))
        for region in ds.region.values:
            if (season, row_name, region) not in stds.index:
                continue
            ref = float(stds[season, row_name, region])
            value = float(rmse.sel(region=region))
            rows.append({"season": season, "group": ds.attrs["group"], "row": row_name, "region": region,
                         "test": experiments[season], "test_state": state, "rmse": value,
                         "reference_std": ref, "nrmse": value / ref})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ plotting
def _latitude(value):
    return f"{abs(value)}\u00b0{'N' if value > 0 else 'S' if value < 0 else ''}"


def plot_change_heatmap(table, *, variable_order, experiment_order, experiment_label, region_order=tuple(REGIONS),
                        limit=40, title="", cbar_label="RMSE change (%)", fontsize=8,
                        value="change_pct", levels=None, cmap="RdBu_r", fmt="{:.0f}", region_bounds=None):
    """Scorecard heat map: experiments (within regions) on y, variables on x.

    ``table`` has one ``value`` per (row, region, experiment), e.g. a
    :func:`change_table` or :func:`normalized_rmse_table` subset;
    ``experiment_label(record) -> str`` names its experiment. Only variables in
    ``variable_order`` are drawn. By default (% changes) the colours are
    discrete levels symmetric about 0 out to about +-``limit``, blue = error
    reduced; pass ``levels`` and ``cmap`` for a sequential metric. If the table
    has a ``significant`` column, cells whose 95% bootstrap interval includes 0
    are hatched with grey text. ``region_bounds`` overrides the latitude range
    shown in the panel titles (default :data:`REGIONS`).
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm
    from matplotlib.ticker import MaxNLocator
    from matplotlib.patches import Rectangle

    table = table.assign(experiment=table.apply(experiment_label, axis=1))
    columns = [v for v in variable_order if v in set(table["row"])]
    rows = [(region, e) for region in region_order for e in experiment_order]
    values = np.full((len(rows), len(columns)), np.nan)
    significant = np.zeros_like(values, dtype=bool)
    for record in table.itertuples():
        if record.row in columns and (record.region, record.experiment) in rows:
            i, j = rows.index((record.region, record.experiment)), columns.index(record.row)
            values[i, j] = getattr(record, value)
            significant[i, j] = getattr(record, "significant", True)

    # One panel per region, stacked; fixed sizes in inches: experiment labels
    # left, panel titles above each panel, variable labels below the last
    # panel, colour bar right. Variable groups are separated by white lines
    # only; name them in the caption.
    n_exp = len(experiment_order)
    # ``fontsize`` is the tick-label size; all other text and every length
    # in inches scale with it (8 pt = reference layout).
    scale = fontsize / 8
    panel_h, gap, top = 0.32 * n_exp * scale, 0.3 * scale, (0.55 if title else 0.3) * scale
    bottom = (0.075 * max(len(c) for c in columns) + 0.25) * scale   # rotated labels
    left, right = 1.0 * scale, 1.2 * scale
    fig_w = 0.42 * scale * len(columns) + left + right
    fig_h = top + len(region_order) * panel_h + (len(region_order) - 1) * gap + bottom
    fig = plt.figure(figsize=(fig_w, fig_h))
    if levels is None:   # discrete levels on a round step (~16 bins) symmetric about 0
        levels = MaxNLocator(nbins=16, steps=[1, 2, 2.5, 5, 10], symmetric=True).tick_values(-limit, limit)
        extend = "both"
    else:
        levels = np.asarray(levels)
        extend = {(True, True): "both", (True, False): "min", (False, True): "max", (False, False): "neither"}[
            bool(np.nanmin(values) < levels[0]), bool(np.nanmax(values) > levels[-1])]
    cmap = plt.get_cmap(cmap)
    norm = BoundaryNorm(levels, cmap.N, extend=extend)
    group_of = table.drop_duplicates("row").set_index("row")["group"]
    groups = list(dict.fromkeys(group_of[c] for c in columns))
    axes = []
    for k, region in enumerate(region_order):
        y0 = fig_h - top - (k + 1) * panel_h - k * gap
        ax = fig.add_axes([left / fig_w, y0 / fig_h, 1 - (left + right) / fig_w, panel_h / fig_h])
        axes.append(ax)
        block = values[k * n_exp:(k + 1) * n_exp]
        sig = significant[k * n_exp:(k + 1) * n_exp]
        image = ax.imshow(block, cmap=cmap, norm=norm, aspect="auto")
        for i in range(n_exp):
            for j in range(len(columns)):
                if np.isnan(block[i, j]):
                    continue
                if not sig[i, j]:
                    ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, hatch="///",
                                           edgecolor="#9a9a9a", linewidth=0))
                r, g, b, _ = cmap(norm(block[i, j]))
                text_color = "#1a1a1a" if 0.299 * r + 0.587 * g + 0.114 * b > 0.5 else "white"
                text = fmt.format(block[i, j])
                text = fmt.format(0.0) if float(text) == 0 else text   # no "-0"
                ax.text(j, i, text, ha="center", va="center", fontsize=0.875 * fontsize,
                        color=text_color if sig[i, j] else "#6b6b6b")
        for g in groups[1:]:                           # group separators
            first = min(j for j, c in enumerate(columns) if group_of[c] == g)
            ax.axvline(first - 0.5, color="white", linewidth=3 * scale)
        ax.set_xticks(range(len(columns)))
        ax.set_xticklabels(columns if k == len(region_order) - 1 else [], rotation=90, fontsize=fontsize)
        ax.set_yticks(range(n_exp))
        ax.set_yticklabels(experiment_order, fontsize=fontsize)
        south, north = (region_bounds or REGIONS)[region]
        label = f"({'abcdefgh'[k]}) {region}"
        if (south, north) != (-90, 90):
            label += f" ({_latitude(south)}–{_latitude(north)})"
        ax.set_title(label, loc="left", fontsize=1.125 * fontsize, pad=3 * scale)
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
    first, last = axes[0].get_position(), axes[-1].get_position()
    cax = fig.add_axes([first.x1 + 0.2 * scale / fig_w, last.y0, 0.14 * scale / fig_w, first.y1 - last.y0])
    cbar = fig.colorbar(image, cax=cax, orientation="vertical", extend=extend, ticks=levels[::2],
                        drawedges=True)
    cbar.dividers.set_color("white")
    cbar.dividers.set_linewidth(0.6 * scale)
    cbar.outline.set_visible(False)
    cbar.set_label(cbar_label, fontsize=fontsize)
    cbar.ax.tick_params(labelsize=0.875 * fontsize)
    if title:
        fig.suptitle(title, fontsize=1.25 * fontsize, x=first.x0 + first.width / 2, y=1 - 0.15 * scale / fig_h)
    return fig
