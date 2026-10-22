"""Daily ensemble performance of the hindcasts: RMSE, spread/error, bias, fair CRPS.

Workflow: 6-hourly member fields -> daily member fields -> ensemble diagnostics.
Daily fields are formed on 24-hour forecast-lead bins (day k covers lead hours
(24(k-1), 24k]), not calendar days, so starts at any hour are handled:

* interval means (``time: mean``, e.g. TREFHT, PRECT rate, FLUT) average the
  intervals inside the bin, weighted by width; the zero-width lead-0 snapshot
  is never used;
* instantaneous values (``time: point``, e.g. PSL) form a time-weighted
  (trapezoid) mean of the samples at the bin edges and inside, which includes
  the initial state at lead 0 for day 1.

References are averaged over the same bins (``util.forecast_ensemble_ts.read_reference``).
For each day, over the gridpoints a region and the reference share (same mask
and cos-latitude weights for every metric):

* ``rmse``: RMSE of the ensemble mean;
* ``bias``: mean of ensemble mean minus reference;
* ``spread``: RMS ensemble spread, sqrt(<var_members>), unbiased member variance;
* ``spread_corr``: spread * sqrt((N + 1) / N), the finite-ensemble correction, so
  spread_corr / rmse has expectation 1 for a reliable N-member ensemble;
* ``crps_fair``: fair CRPS (ensemble-size adjusted), area mean;
* ``member_mean`` (region, member, day) and ``reference_mean``: regional means
  of each member and of the reference, for the field time series.

Ensembles larger than ``subset_size`` are also scored on random subsets of that
size (``*_subset``), to separate initialization effects from averaging more members.
"""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from util.bias_metrics import experiment_registry
from util.forecast_ensemble_ts import (
    REFERENCES, fit_legend, read_member_ensemble, read_reference, reference_product, region_masks,
)
from util.model_evaluation import experiment_colors

CACHE_VERSION = 2  # 2: regional means of each member and of the reference (field time series)
PRODUCT = "daily_skill"
SPREAD_DEFINITION = "rms of unbiased member variance; corrected by sqrt((N+1)/N)"
AGGREGATION = "24-h forecast-lead bins from 6-hourly member fields"
# Expected temporal semantics of the 6-hourly model output (checked, not assumed).
EXPECTED_SEMANTICS = {
    "PSL": {"cell_method": "time: point", "rate": False},
    "TREFHT": {"cell_method": "time: mean", "rate": False},
    "TS": {"cell_method": "time: mean", "rate": False},
    "PRECT": {"cell_method": "time: mean", "rate": True},
    "FLUT": {"cell_method": "time: mean", "rate": False},
}
METRICS = ("rmse", "bias", "spread", "spread_corr", "crps_fair")


DAILY_MODE_DEFAULTS = {"references": {}, "max_day": 15, "n_subsets": 100, "subset_size": 10, "seed": 42,
                       "baselines": ("CTRL10", "CAPT10")}
DAILY_MODE_OPTIONS = ("kind", "variables", *DAILY_MODE_DEFAULTS)


def validate_daily_mode(mode, modes):
    """Settings of a daily-performance mode with defaults filled in, after checking every option."""
    if mode not in modes:
        raise ValueError(f"Unknown mode {mode!r}; choose from {list(modes)}")
    config = {**DAILY_MODE_DEFAULTS, **modes[mode]}
    unknown = set(config) - set(DAILY_MODE_OPTIONS)
    if unknown:
        raise ValueError(f"{mode}: unknown options {sorted(unknown)}; choose from {list(DAILY_MODE_OPTIONS)}")
    if "variables" not in config:
        raise ValueError(f"Mode {mode!r} lacks variables")
    for variable in config["variables"]:
        reference_product(variable, config["references"].get(variable))
        if variable not in REFERENCES:
            raise ValueError(f"{variable} has no reference; choose from {list(REFERENCES)}")
    if config["subset_size"] < 2:
        raise ValueError("subset_size must be at least 2")
    return config


# --------------------------------------------------------------------------
# Temporal semantics and daily aggregation
# --------------------------------------------------------------------------
def check_temporal_semantics(config, variables, period, *, frequency="6hourly"):
    """Table of each variable's stored cell method, units, and record layout (first member).

    Raises when a variable's cell method or rate units differ from ``EXPECTED_SEMANTICS``.
    """
    rows = []
    for variable in variables:
        field, opened = read_member_ensemble(config, variable, period, frequency=frequency, max_members=1)
        try:
            start = pd.DatetimeIndex(field.interval_start.values)
            end = pd.DatetimeIndex(field.interval_end.values)
            width = (end - start) / pd.Timedelta(hours=1)
            method = str(field.attrs.get("cell_methods", ""))
            units = str(opened[0][variable].attrs.get("units", ""))
        finally:
            for ds in opened:
                ds.close()
        expected = EXPECTED_SEMANTICS.get(variable)
        point = "time: point" in method
        row = dict(variable=variable, cell_methods=method, units=units,
                   first_record_width_h=float(width[0]), record_width_h=float(np.median(width[1:])),
                   aggregation="trapezoid of instantaneous values" if point else "width-weighted mean of interval means")
        if expected is not None:
            if expected["cell_method"] not in method:
                raise ValueError(f"{variable}: cell_methods {method!r}, expected {expected['cell_method']!r}")
            is_rate = "/s" in units or "s-1" in units or "s**-1" in units
            if expected["rate"] != is_rate:
                raise ValueError(f"{variable}: units {units!r} {'are not' if expected['rate'] else 'are'} a rate")
            row["check"] = "ok"
        else:
            row["check"] = "not in EXPECTED_SEMANTICS"
        rows.append(row)
    return pd.DataFrame(rows)


def lead_day_weights(field):
    """(day, lead) weights that turn ``field`` records into complete 24-h lead bins.

    Returns ``(weights, starts, ends)``; incomplete bins (e.g. the last day of
    the run) are left out.
    """
    start = pd.DatetimeIndex(field.interval_start.values)
    end = pd.DatetimeIndex(field.interval_end.values)
    init = start[0]
    s = np.asarray((start - init) / pd.Timedelta(hours=1), dtype=float)
    e = np.asarray((end - init) / pd.Timedelta(hours=1), dtype=float)
    point = "time: point" in str(field.attrs.get("cell_methods", ""))
    rows, starts, ends = [], [], []
    for day in range(1, int(e.max() // 24) + 1):
        a, b = 24.0 * (day - 1), 24.0 * day
        weights = np.zeros(len(s))
        if point:
            inside = np.where((e >= a - 1e-6) & (e <= b + 1e-6))[0]
            times = e[inside]
            if inside.size < 2 or abs(times[0] - a) > 1e-6 or abs(times[-1] - b) > 1e-6:
                continue
            gaps = np.diff(times)
            weights[inside] = np.concatenate([[0.0], gaps / 2]) + np.concatenate([gaps / 2, [0.0]])
        else:
            inside = (e - s > 0) & (s >= a - 1e-6) & (e <= b + 1e-6)
            weights[inside] = (e - s)[inside]
            if abs(weights.sum() - 24.0) > 1e-6:
                continue
        rows.append(weights / weights.sum())
        starts.append(init + pd.Timedelta(hours=a))
        ends.append(init + pd.Timedelta(hours=b))
    if not rows:
        raise ValueError("No complete 24-h lead bin")
    days = np.arange(1, len(rows) + 1)
    weights = xr.DataArray(np.array(rows), dims=("day", "lead"), coords={"day": days, "lead": field.lead})
    return weights, pd.DatetimeIndex(starts), pd.DatetimeIndex(ends)


def daily_member_fields(field):
    """Daily (member, day, lat, lon) fields and their bin edges."""
    weights, starts, ends = lead_day_weights(field)
    daily = xr.dot(weights, field.drop_vars(["interval_start", "interval_end"]), dim="lead")
    daily = daily.transpose("member", "day", "lat", "lon")
    daily.attrs = dict(field.attrs)
    return daily, starts, ends


def daily_reference(reference_root, variable, product, starts, ends, lat, lon):
    """Reference averaged over the same 24-h lead bins, (day, lat, lon)."""
    days = np.arange(1, len(starts) + 1)
    template = xr.Dataset(coords=dict(lead=days, lat=lat, lon=lon,
                                      interval_start=("lead", starts.values),
                                      interval_end=("lead", ends.values)),
                          attrs={"cell_methods": "time: mean"})
    reference = read_reference(reference_root, variable, template, product)
    return reference.rename(lead="day").drop_vars(["interval_start", "interval_end"], errors="ignore")


def _stored_daily_kind(config):
    """ts subdirectory holding the model's own daily series: daily_mean (current layout), else h1/daily."""
    root = Path(config["path"]) / "EN01" / "archive/post/atm/180x360_aave/ts"
    return next((k for k in ("daily_mean", "h1", "daily") if (root / k).is_dir()), "daily")


def validate_daily_aggregation(config, variable, period, *, days=3):
    """Compare lead-bin daily means with the model's stored daily file (member 1, first ``days``).

    The stored file is ``ts/daily_mean`` (00Z-24Z means built from the 6-hourly
    output; instantaneous PSL as the mean of the 00/06/12/18Z states) or, in the
    older layout, the EAM h1 series (``ts/h1`` or ``ts/daily``), whose state
    variables are 00Z snapshots and are reported but not compared.
    """
    kind = _stored_daily_kind(config)
    stored, opened = read_member_ensemble(config, variable, period, frequency=kind, max_members=1)
    try:
        method = str(stored.attrs.get("cell_methods", ""))
        result = {"variable": variable, "stored_daily_file": f"ts/{kind}", "daily_file_cell_methods": method}
        if "time: point" in method:
            result["comparison"] = "daily file holds snapshots; not a daily mean"
            return result
        stored = stored.isel(member=0)
        keep = (stored.interval_end > stored.interval_start).values
        stored = stored.isel(lead=keep)
        stored_start = pd.DatetimeIndex(stored.interval_start.values)
    finally:
        for ds in opened:
            ds.close()
    field, opened = read_member_ensemble(config, variable, period, frequency="6hourly", max_members=1)
    try:
        daily, starts, _ = daily_member_fields(field)
        ours = daily.isel(member=0, day=slice(0, days)).values
    finally:
        for ds in opened:
            ds.close()
    stored, opened = read_member_ensemble(config, variable, period, frequency=kind, max_members=1)
    try:
        stored = stored.isel(member=0, lead=np.flatnonzero(keep)[[int(np.flatnonzero(stored_start == t)[0]) for t in starts[:days]]]).values
    finally:
        for ds in opened:
            ds.close()
    difference = np.abs(ours - stored)
    result.update(comparison="compared", max_abs_difference=float(np.nanmax(difference)),
                  mean_abs_difference=float(np.nanmean(difference)),
                  mean_abs_value=float(np.nanmean(np.abs(ours))))
    return result


# --------------------------------------------------------------------------
# Metrics
# --------------------------------------------------------------------------
def _area_mean(values, weights):
    """Weighted area mean over (lat, lon) of numpy arrays (..., lat, lon)."""
    total = weights.sum(axis=(-2, -1))
    return np.where(total > 0, (np.nan_to_num(values) * weights).sum(axis=(-2, -1)) / np.where(total > 0, total, 1), np.nan)


def _fair_crps(members, observed):
    """Fair CRPS per day/gridpoint; ``members`` (N, ...) and ``observed`` (...)."""
    n = members.shape[0]
    term1 = np.mean(np.abs(members - observed[None]), axis=0)
    ordered = np.sort(members, axis=0)
    rank = (2 * np.arange(1, n + 1) - n - 1).reshape((n,) + (1,) * (members.ndim - 1))
    return term1 - (rank * ordered).sum(axis=0) / (n * (n - 1))


def daily_metrics(daily, reference, masks, *, n_subsets=100, subset_size=10, seed=42):
    """Regional daily metrics (region, day); subset metrics (region, subset, day) when N > subset_size."""
    members = np.asarray(daily.values, dtype=np.float64)  # (N, day, lat, lon)
    observed = np.asarray(reference.transpose("day", "lat", "lon").values, dtype=np.float64)
    n = members.shape[0]
    cos = np.cos(np.deg2rad(daily.lat.values)).clip(min=0)[:, None]
    valid = np.isfinite(observed) & np.isfinite(members).all(axis=0)
    regions = masks.region.values.tolist()
    region_weights = [np.where(valid & masks.sel(region=r).values[None], cos, 0.0) for r in regions]

    mean = members.mean(axis=0)
    error = mean - observed
    variance = members.var(axis=0, ddof=1) if n > 1 else np.full_like(mean, np.nan)
    crps = _fair_crps(members, observed) if n > 1 else np.abs(error)
    out = {}
    for name, field in (("rmse", error ** 2), ("bias", error), ("spread", variance), ("crps_fair", crps)):
        values = np.stack([_area_mean(field, w) for w in region_weights])
        out[name] = np.sqrt(values) if name in ("rmse", "spread") else values
    out["spread_corr"] = out["spread"] * np.sqrt((n + 1) / n)
    out["spread_ratio"] = out["spread_corr"] / out["rmse"]
    out["spread_ratio_raw"] = out["spread"] / out["rmse"]
    coords = {"region": regions, "day": daily.day.values}
    ds = xr.Dataset({k: (("region", "day"), v) for k, v in out.items()}, coords=coords)
    ds["member_mean"] = (("region", "member", "day"),
                         np.stack([_area_mean(members, w[None]) for w in region_weights]))
    ds["reference_mean"] = (("region", "day"), np.stack([_area_mean(observed, w) for w in region_weights]))

    if n > subset_size and n_subsets > 0:
        rng = np.random.default_rng(seed)
        rmse, ratio = [], []
        for _ in range(n_subsets):
            pick = members[rng.choice(n, subset_size, replace=False)]
            sub_error = pick.mean(axis=0) - observed
            sub_var = pick.var(axis=0, ddof=1)
            r = np.sqrt(np.stack([_area_mean(sub_error ** 2, w) for w in region_weights]))
            s = np.sqrt(np.stack([_area_mean(sub_var, w) for w in region_weights])) * np.sqrt((subset_size + 1) / subset_size)
            rmse.append(r)
            ratio.append(s / r)
        ds["rmse_subset"] = (("subset", "region", "day"), np.array(rmse))
        ds["spread_ratio_subset"] = (("subset", "region", "day"), np.array(ratio))
        ds = ds.transpose("region", "subset", "member", "day")
    ds["n_points"] = (("region", "day"), np.stack([(w > 0).sum(axis=(-2, -1)) for w in region_weights]))
    return ds


# --------------------------------------------------------------------------
# Caches and driver
# --------------------------------------------------------------------------
def daily_skill_cache_path(cache_dir, experiment, variable, period):
    return Path(cache_dir) / PRODUCT / f"{experiment}_{variable}_{period}_daily_skill.nc"


def _settings(*, reference, members, n_subsets, subset_size, seed):
    return dict(product=PRODUCT, cache_version=CACHE_VERSION, aggregation=AGGREGATION,
                spread_definition=SPREAD_DEFINITION, reference=reference, ensemble_members=members,
                n_subsets=n_subsets, subset_size=subset_size, seed=seed)


def _cache_status(path, settings, regions):
    if not Path(path).is_file():
        return "missing"
    with xr.open_dataset(path) as cached:
        for key, value in settings.items():
            if str(cached.attrs.get(key)) != str(value):
                return f"{key}: cached {cached.attrs.get(key)}, requested {value}"
        missing = set(regions) - set(cached.region.values.tolist())
        if missing:
            return f"regions {sorted(missing)} not cached"
    return ""


def read_or_process_daily_skill(cases, variables, *, data_root, reference_root, cache_dir, regions,
                                landmask_file, references=None, max_members=None, n_subsets=100,
                                subset_size=10, seed=42, force_compute=False, log=print):
    """Reuse or build daily metrics per season/variable/experiment.

    Returns ``(results, registries, failures)`` with
    ``results[season][variable][experiment] -> Dataset``.
    """
    references = references or {}
    results, registries, failures = {}, {}, []
    masks = None
    for season, experiments in cases.items():
        registry = registries[season] = experiment_registry(data_root, season, "fc")
        unknown = set(experiments) - set(registry)
        if unknown:
            raise ValueError(f"{season}: {sorted(unknown)} unavailable; choose from {list(registry)}")
        for variable in variables:
            product = reference_product(variable, references.get(variable))
            if product is None:
                raise ValueError(f"{variable} has no reference; choose from {list(REFERENCES)}")
            for experiment in experiments:
                meta = registry[experiment]
                members = meta["nens"] if max_members is None else min(max_members, meta["nens"])
                settings = _settings(reference=product, members=members, n_subsets=n_subsets,
                                     subset_size=subset_size, seed=seed)
                path = daily_skill_cache_path(cache_dir, experiment, variable, meta["period"])
                try:
                    reason = "forced" if force_compute else _cache_status(path, settings, regions)
                    if not reason:
                        ds = xr.load_dataset(path).sel(region=list(regions))
                        log(f"[REUSE] {path.name}")
                    else:
                        field, opened = read_member_ensemble(meta, variable, meta["period"],
                                                             frequency="6hourly", max_members=max_members)
                        try:
                            daily, starts, ends = daily_member_fields(field)
                            daily = daily.load()
                        finally:
                            for opened_ds in opened:
                                opened_ds.close()
                        reference = daily_reference(reference_root, variable, product, starts, ends,
                                                    daily.lat, daily.lon)
                        if masks is None:
                            masks = region_masks(daily.lat, daily.lon, regions, landmask_file)
                        ds = daily_metrics(daily, reference, masks, n_subsets=n_subsets,
                                           subset_size=subset_size, seed=seed)
                        ds = ds.assign_coords(day_start=("day", starts.values), day_end=("day", ends.values))
                        ds.attrs.update(settings, experiment=experiment, label=meta["label"],
                                        variable=variable, period=meta["period"],
                                        units=str(daily.attrs.get("units", "")),
                                        cell_methods=str(daily.attrs.get("cell_methods", "")))
                        path.parent.mkdir(parents=True, exist_ok=True)
                        temporary = path.with_suffix(".nc.tmp")
                        ds.to_netcdf(temporary)
                        os.replace(temporary, path)
                        log(f"[WRITE] {path.name} ({reason})")
                except Exception as error:  # keep going; report at the end
                    failures.append((season, variable, experiment, f"{type(error).__name__}: {error}"))
                    log(f"[FAIL] {season} {variable} {experiment}: {type(error).__name__}: {error}")
                    continue
                results.setdefault(season, {}).setdefault(variable, {})[experiment] = ds
    return results, registries, failures


# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------
SUMMARY_WINDOWS = {"days1-3": (1, 3), "days4-10": (4, 10)}


def _window_mean(values, days, window):
    d0, d1 = window
    return float(values.where((days >= d0) & (days <= d1)).mean())


def daily_skill_summary(results, *, windows=SUMMARY_WINDOWS, ratio_days=(1, 10),
                        baselines=("CTRL10", "CAPT10")):
    """One row per season/variable/experiment/region with window means and % RMSE changes.

    ``baselines`` are experiment-name prefixes; ``rmse_<window>_vs_<baseline>_pct``
    compares with the season's matching experiment. Subset columns use the
    ``subset_size``-member random subsets (mean and 5-95% range).
    """
    rows = []
    for season, by_variable in results.items():
        for variable, datasets in by_variable.items():
            for experiment, ds in datasets.items():
                for region in ds.region.values:
                    sel = ds.sel(region=region)
                    row = dict(season=season, variable=variable, units=ds.attrs.get("units", ""),
                               experiment=experiment, members=int(ds.attrs.get("ensemble_members", 0)),
                               region=str(region), reference=ds.attrs.get("reference", ""))
                    for name, window in windows.items():
                        for metric in ("rmse", "bias", "crps_fair"):
                            row[f"{metric}_{name}"] = _window_mean(sel[metric], sel.day, window)
                        if "rmse_subset" in sel:
                            per_subset = sel.rmse_subset.where((sel.day >= window[0]) & (sel.day <= window[1])).mean("day")
                            row[f"rmse_{name}_subset_mean"] = float(per_subset.mean())
                            row[f"rmse_{name}_subset_p05"] = float(per_subset.quantile(0.05))
                            row[f"rmse_{name}_subset_p95"] = float(per_subset.quantile(0.95))
                    label = f"days{ratio_days[0]}-{ratio_days[1]}"
                    row[f"spread_ratio_{label}"] = _window_mean(sel.spread_ratio, sel.day, ratio_days)
                    row[f"spread_ratio_raw_{label}"] = _window_mean(sel.spread_ratio_raw, sel.day, ratio_days)
                    rows.append(row)
    table = pd.DataFrame(rows)
    if table.empty:
        return table
    for prefix in baselines:
        for name in windows:
            column = f"rmse_{name}"
            reference = {}
            for (season, variable, region), group in table.groupby(["season", "variable", "region"]):
                match = group[group["experiment"].str.startswith(prefix)]
                if len(match) == 1:
                    reference[(season, variable, region)] = float(match[column].iloc[0])
            base = table.apply(lambda r: reference.get((r["season"], r["variable"], r["region"]), np.nan), axis=1)
            table[f"{column}_vs_{prefix}_pct"] = 100.0 * (table[column] / base - 1.0)
            if f"{column}_subset_mean" in table:
                table[f"{column}_subset_vs_{prefix}_pct"] = 100.0 * (table[f"{column}_subset_mean"] / base - 1.0)
    return table


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------
PANEL_VARIABLES = ("PSL", "TREFHT", "PRECT", "FLUT")


def results_colors(results):
    return experiment_colors(e for by_variable in results.values()
                             for datasets in by_variable.values() for e in datasets)


_UNIT_LABELS = {"mm day-1": r"mm day$^{-1}$", "mm/day": r"mm day$^{-1}$", "W/m2": r"W m$^{-2}$",
                "W m-2": r"W m$^{-2}$", "m/s": r"m s$^{-1}$", "m s-1": r"m s$^{-1}$"}


def _days_since_label(datasets_by_variable):
    """"Days since YYYY-MM-DD" from the start of day 1 of the first dataset ("Forecast day" without one)."""
    ds = next(iter(next(iter(datasets_by_variable.values())).values()))
    if "day_start" not in ds.coords:
        return "Forecast day"
    return f"Days since {pd.Timestamp(ds.day_start.sel(day=ds.day.min()).values):%Y-%m-%d}"


def _initial_month(datasets_by_variable):
    """Month name of the start of day 1 ("January"), or "" without ``day_start``."""
    ds = next(iter(next(iter(datasets_by_variable.values())).values()))
    if "day_start" not in ds.coords:
        return ""
    return f"{pd.Timestamp(ds.day_start.sel(day=ds.day.min()).values):%B}"


def plot_daily_performance(datasets_by_variable, region, *, labels, colors, variables=PANEL_VARIABLES,
                           max_day=15, title="", fontz=16, line_width=2.2, line_styles=None, marker="o",
                           markers=None, marker_size=4, subset_line_width=1.6, subset_line_style="--",
                           subset_alpha=0.18, reference_line=None, grid_alpha=0.3, panel_size=(7.5, 3.6),
                           legend_frame=True, legend_ncol=None, legend_edgecolor="#b5b4ae", xlabel=None):
    """Fixed rows x 2 figure: daily ensemble-mean RMSE (left) and spread_corr/RMSE (right).

    Panels are labelled "(a) Ensemble-Mean RMSE (PSL)", "(b) Spread/RMSE (PSL)",
    ...; references, the spread correction, and the subsets belong in the
    caption. For ensembles larger than the subset size, the RMSE panels also
    show the mean RMSE of the random equal-size subsets (dashed) and their
    5-95% range (shaded), so every experiment can be compared at the same
    ensemble size. ``fontz`` is the base font size; every label scales from it.

    Line style: ``line_width``, ``marker`` (None: no markers), and ``marker_size``
    apply to every experiment; ``line_styles`` and ``markers`` ({experiment:
    value}) override them per experiment. Subsets: ``subset_line_width``,
    ``subset_line_style``, ``subset_alpha`` (shading; 0 hides it).
    ``reference_line``: keyword arguments of the spread/RMSE = 1 line.
    ``panel_size``: (width, height) of one panel in inches.
    ``xlabel``: None gives "Days since <initial date>" (the start of day 1;
    day k is the 24 h ending k days after it); any string replaces it.
    Legend (below the panels): ``legend_frame`` draws a rounded box with edge
    ``legend_edgecolor``; ``legend_ncol`` columns (None: one row). The subset
    entry shows the shading behind the subset line, as in the panels.
    """
    variables = [v for v in variables if v in datasets_by_variable]
    line_styles, markers = line_styles or {}, markers or {}
    reference_line = {"color": "#6f6e69", "lw": 1.2, "ls": "--", **(reference_line or {})}
    fig, axes = plt.subplots(len(variables), 2, figsize=(2 * panel_size[0], panel_size[1] * len(variables)),
                             sharex=True, squeeze=False)
    handles = {}
    subset_size = None
    for row, variable in enumerate(variables):
        datasets = datasets_by_variable[variable]
        ax_rmse, ax_ratio = axes[row]
        units = ""
        for experiment, ds in datasets.items():
            sel = ds.sel(region=region).where(ds.day <= max_day, drop=True)
            style = {"color": colors[experiment], "lw": line_width, "ls": line_styles.get(experiment, "-"),
                     "marker": markers.get(experiment, marker), "ms": marker_size}
            line, = ax_rmse.plot(sel.day, sel.rmse, **style)
            ax_ratio.plot(sel.day, sel.spread_ratio, **style)
            if "rmse_subset" in sel:
                low, high = sel.rmse_subset.quantile([0.05, 0.95], "subset")
                if subset_alpha:
                    ax_rmse.fill_between(sel.day, low, high, color=colors[experiment], alpha=subset_alpha, lw=0)
                ax_rmse.plot(sel.day, sel.rmse_subset.mean("subset"), color=colors[experiment],
                             lw=subset_line_width, ls=subset_line_style)
                subset_size = ds.attrs.get("subset_size", subset_size)
            handles[experiment] = line
            units = ds.attrs.get("units", "")
        ax_ratio.axhline(1.0, **reference_line)
        ax_rmse.set_title(f"({chr(97 + 2 * row)}) Ensemble-Mean RMSE ({variable})", loc="left",
                          fontsize=fontz * 0.9, color="#2b2b29")
        ax_ratio.set_title(f"({chr(98 + 2 * row)}) Spread/RMSE ({variable})", loc="left",
                           fontsize=fontz * 0.9, color="#2b2b29")
        ax_rmse.set_ylabel(_UNIT_LABELS.get(units, units), fontsize=fontz * 0.8)
        for ax in (ax_rmse, ax_ratio):
            ax.grid(alpha=grid_alpha)
            ax.set_xlim(0.5, max_day + 0.5)
            ax.tick_params(labelsize=fontz * 0.75)
    for ax in axes[-1]:
        ax.set_xlabel(xlabel or _days_since_label(datasets_by_variable), fontsize=fontz * 0.8)
    legend_handles = [handles[e] for e in handles]
    legend_labels = [labels[e] for e in handles]
    if subset_size is not None:
        subset_line = plt.Line2D([], [], color="#6f6e69", lw=subset_line_width, ls=subset_line_style)
        legend_handles.append((plt.Rectangle((0, 0), 1, 1, color="#9a9a9a", alpha=min(1, 2 * subset_alpha), lw=0),
                               subset_line) if subset_alpha else subset_line)
        legend_labels.append(f"{subset_size}-member subsets" + (" (5-95%)" if subset_alpha else ""))
    _figure_legend(fig, legend_handles, legend_labels, height=panel_size[1] * len(variables), fontz=fontz,
                   title=title, ncol=legend_ncol, frame=legend_frame, edgecolor=legend_edgecolor)
    return fig


DAILY_METRIC_NAMES = {"rmse": "Ensemble-Mean RMSE", "spread_ratio": "Spread/RMSE"}


def plot_daily_metric(results, region, *, metric, labels, colors, variables=PANEL_VARIABLES, seasons=None,
                      max_day=15, title="", fontz=16, line_width=2.2, line_styles=None, marker="o",
                      markers=None, marker_size=4, subset_line_width=1.6, subset_line_style="--",
                      subset_alpha=0.18, reference_line=None, grid_alpha=0.3, panel_size=(7.5, 3.6),
                      legend_frame=True, legend_ncol=None, legend_edgecolor="#b5b4ae", xlabel=None,
                      season_names=None, legend_by="experiment", legend_fontz=None):
    """One daily metric, one row per variable and one column per season (S0 left, S1 right).

    ``results[season][variable][experiment]`` as returned by
    ``read_or_process_daily_skill``; ``metric`` is "rmse" (ensemble-mean RMSE,
    with the equal-size subsets of the larger ensembles as in
    ``plot_daily_performance``) or "spread_ratio" (spread_corr/RMSE, with the
    ``reference_line`` at 1). Panels are labelled "(a) PSL (January)", row by
    row; the metric is on the y axis. ``season_names``: {season: name} for the
    panel labels; None names each season by the month of its initial date.
    ``legend_by``: "experiment" lists every experiment (``CTRL10-S0``,
    ``CTRL10-S1``, ...), one legend row per season with each experiment under
    its type; "type" lists each type once (``CTRL10``). Experiments run for one
    season only share a column (DART10-S0 above DART20-S0). The subset entry
    comes last, in a column of its own after the last experiment. ``legend_fontz``:
    legend font size (None: 0.85 x ``fontz``); the box scales with it.
    ``xlabel``: None gives "Days since <initial date>" under each season column.
    Other style arguments as in ``plot_daily_performance``.
    """
    if metric not in DAILY_METRIC_NAMES:
        raise ValueError(f"metric must be one of {list(DAILY_METRIC_NAMES)}")
    if legend_by not in ("experiment", "type"):
        raise ValueError('legend_by must be "experiment" or "type"')
    seasons = [s for s in (seasons or results) if s in results]
    variables = [v for v in variables if any(v in results[s] for s in seasons)]
    line_styles, markers = line_styles or {}, markers or {}
    reference_line = {"color": "#6f6e69", "lw": 1.2, "ls": "--", **(reference_line or {})}
    nrows, ncols = len(variables), len(seasons)
    names = {season: (season_names or {}).get(season) or _initial_month(results[season]) for season in seasons}
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * panel_size[0], nrows * panel_size[1]),
                             sharex="col", squeeze=False)
    handles = {}
    subset_size = None
    for col, season in enumerate(seasons):
        for row, variable in enumerate(variables):
            ax = axes[row, col]
            units = ""
            for experiment, ds in results[season].get(variable, {}).items():
                sel = ds.sel(region=region).where(ds.day <= max_day, drop=True)
                line, = ax.plot(sel.day, sel[metric], color=colors[experiment], lw=line_width,
                                ls=line_styles.get(experiment, "-"), marker=markers.get(experiment, marker),
                                ms=marker_size)
                if metric == "rmse" and "rmse_subset" in sel:
                    low, high = sel.rmse_subset.quantile([0.05, 0.95], "subset")
                    if subset_alpha:
                        ax.fill_between(sel.day, low, high, color=colors[experiment], alpha=subset_alpha, lw=0)
                    ax.plot(sel.day, sel.rmse_subset.mean("subset"), color=colors[experiment],
                            lw=subset_line_width, ls=subset_line_style)
                    subset_size = ds.attrs.get("subset_size", subset_size)
                name = labels[experiment] if legend_by == "experiment" else labels[experiment].rsplit("-", 1)[0]
                handles.setdefault(season if legend_by == "experiment" else "", {}).setdefault(name, line)
                units = ds.attrs.get("units", "")
            if metric == "spread_ratio":
                ax.axhline(1.0, **reference_line)
                ax.set_ylabel("Spread/RMSE", fontsize=fontz * 0.8)
            else:
                ax.set_ylabel(f"RMSE ({_UNIT_LABELS.get(units, units)})" if units else "RMSE",
                              fontsize=fontz * 0.8)
            ax.set_title(f"({chr(97 + row * ncols + col)}) {variable} ({names[season]})", loc="left",
                         fontsize=fontz * 0.9, color="#2b2b29")
            ax.grid(alpha=grid_alpha)
            ax.set_xlim(0.5, max_day + 0.5)
            ax.tick_params(labelsize=fontz * 0.75)
        axes[-1, col].set_xlabel(xlabel or _days_since_label(results[season]), fontsize=fontz * 0.8)
    extra = []
    if subset_size is not None:
        subset_line = plt.Line2D([], [], color="#6f6e69", lw=subset_line_width, ls=subset_line_style)
        extra.append((f"{subset_size}-member subsets" + (" (5-95%)" if subset_alpha else ""),
                      (plt.Rectangle((0, 0), 1, 1, color="#9a9a9a", alpha=min(1, 2 * subset_alpha), lw=0),
                       subset_line) if subset_alpha else subset_line))
    _season_legend(fig, handles, extra, height=panel_size[1] * nrows, fontz=fontz, title=title,
                   ncol=legend_ncol, frame=legend_frame, edgecolor=legend_edgecolor, legend_fontz=legend_fontz,
                   extra_last=True, pack=True)
    return fig


def _figure_legend(fig, handles, labels, *, height, fontz, title, ncol=None, frame=True, edgecolor="#b5b4ae",
                   legend_fontz=None):
    """Boxed legend below the panels (one row unless ``ncol``; filled row by row), title, and layout.

    ``legend_fontz``: legend font size (None: 0.85 x ``fontz``); the box and its spacing scale with it.
    A legend wider than the figure is shrunk to fit (``fit_legend``).
    """
    ncol = ncol or len(handles)
    rows = int(np.ceil(len(handles) / ncol))
    # Fill row by row (matplotlib fills legend columns first).
    order = [i for c in range(ncol) for i in range(c, len(handles), ncol)]
    legend, legend_fontz = fit_legend(fig, lambda size: fig.legend(
        [handles[i] for i in order], [labels[i] for i in order],
        loc="lower center", bbox_to_anchor=(0.5, 0.005), ncol=ncol, fontsize=size,
        frameon=frame, fancybox=True, framealpha=1, facecolor="white",
        edgecolor=edgecolor, borderpad=0.6, handlelength=2.6, handletextpad=0.6,
        columnspacing=1.6), legend_fontz or fontz * 0.85)
    legend.get_frame().set_linewidth(0.8)
    bottom = (0.45 + 0.32 * rows) * legend_fontz / 0.85 / 16 / height  # legend band, fraction of the figure height
    if title:
        fig.suptitle(title, fontsize=fontz)
    fig.tight_layout(rect=(0, bottom, 1, 0.98 if title else 1))


def _season_legend(fig, handles, extra, *, height, fontz, title, ncol=None, frame=True, edgecolor="#b5b4ae",
                   legend_fontz=None, extra_last=False, pack=False):
    """Legend with one row per season, each experiment in the column of its type (CTRL10-S1 under CTRL10-S0).

    ``handles``: {group: {name: line}} (one group per season, or a single group
    of types); ``extra``: [(label, handle)] entries (observations, member range,
    subsets) placed in the empty slots of the last row, or in a row of their own
    when they do not fit; ``extra_last=True`` places them after the last
    experiment instead, stacked in end columns of their own (a lone entry on the
    last row), so the experiment columns keep their own widths. ``pack=True`` stacks the
    types missing from some seasons into shared columns, in order (DART10-S0
    above DART20-S0 when both are S0 only). ``ncol`` None: one column per
    experiment type, or per packed column (plus the extra columns).
    """
    types = list(dict.fromkeys(name.rsplit("-", 1)[0] for group in handles.values() for name in group))
    blank = (" ", plt.Line2D([], [], alpha=0))  # a non-empty label, so the slot is kept
    rows = []
    for group in handles.values():
        by_type = {name.rsplit("-", 1)[0]: (name, line) for name, line in group.items()}
        rows.append([by_type.get(kind, blank) for kind in types])
    if pack and len(rows) > 1:  # stack the partly empty columns, keeping the full ones in place
        columns = [list(column) for column in zip(*rows)]
        partial = [k for k, column in enumerate(columns) if any(entry is blank for entry in column)]
        if partial:
            entries = [entry for k in partial for entry in columns[k] if entry is not blank]
            stacked = [entries[i:i + len(rows)] for i in range(0, len(entries), len(rows))]
            stacked[-1] += [blank] * (len(rows) - len(stacked[-1]))
            columns = ([c for k, c in enumerate(columns) if k not in partial and k < partial[0]] + stacked
                       + [c for k, c in enumerate(columns) if k not in partial and k > partial[0]])
            rows = [list(row) for row in zip(*columns)]
    free = [k for k, entry in enumerate(rows[-1]) if entry is blank] if rows else []
    if extra and extra_last and rows:  # extra columns at the end, stacked; a lone entry on the last row
        depth = len(rows)
        columns = [list(extra[i:i + depth]) for i in range(0, len(extra), depth)]
        columns[-1] = [blank] * (depth - len(columns[-1])) + columns[-1]
        rows = [row + [column[r] for column in columns] for r, row in enumerate(rows)]
    elif extra and len(extra) <= len(free):  # fill the empty slots of the last season row (no extra row)
        for k, entry in zip(free, extra):
            rows[-1][k] = entry
    elif extra:
        rows.append(extra)
    ncol = ncol or max(len(r) for r in rows)
    entries = []
    for row in rows:  # pad each row to whole legend rows so every group starts a new row
        entries += row + [blank] * (-len(row) % ncol)
    _figure_legend(fig, [h for _, h in entries], [label for label, _ in entries], height=height,
                   fontz=fontz, title=title, ncol=ncol, frame=frame, edgecolor=edgecolor, legend_fontz=legend_fontz)


def plot_daily_timeseries(results, region, *, labels, colors, variables=("TREFHT", "PRECT", "FLUT"),
                          seasons=None, max_day=15, title="", fontz=16, line_width=2.2, line_styles=None,
                          marker="o", markers=None, marker_size=4, member_display="range", member_alpha=0.15,
                          reference_line=None, reference_name=None, grid_alpha=0.3, panel_size=(5.5, 3.6),
                          legend_frame=True, legend_ncol=None, legend_edgecolor="#b5b4ae", xlabel=None,
                          legend_by="experiment", season_names=None, show_reference=True, legend_fontz=None):
    """Regional-mean daily field: ensemble mean of each experiment and the observations.

    ``results[season][variable][experiment]`` as returned by
    ``read_or_process_daily_skill``; one row per variable and one column per
    season (``seasons``, default all), on the same 24-h forecast-lead days,
    gridpoints, and references as ``plot_daily_performance``.
    ``legend_by``: "experiment" lists every experiment (``CTRL10-S0``,
    ``CTRL10-S1``, ...), one legend row per season; "type" lists each type once
    (``CTRL10``), for figures where a type looks the same in every season.
    ``member_display``: "range" shades the member min-max, "members" draws every
    member thinly, "none" shows only ensemble means (``member_alpha`` sets the
    shading/line opacity). ``reference_line``: keyword arguments of the
    observation line (black by default; ``show_reference=False`` omits it and its
    legend entry); ``reference_name`` its legend label
    (None: "Observations"; the products belong in the caption).
    ``xlabel``: None gives "Days since <initial date>" under each season column.
    ``season_names``: {season: name} for the panel labels, "(a) TREFHT (January)";
    None names each season by the month of its initial date.
    Experiments run for one season only share a legend column (DART10-S0 above
    DART20-S0); observations and member range come last, stacked in a column
    of their own (``legend_ncol`` None: one column per experiment column).
    ``legend_fontz``: legend font size (None: 0.85 x ``fontz``); the box scales with it.
    Other style arguments as in ``plot_daily_performance``.
    """
    if member_display not in ("range", "members", "none"):
        raise ValueError('member_display must be "range", "members", or "none"')
    if legend_by not in ("experiment", "type"):
        raise ValueError('legend_by must be "experiment" or "type"')
    seasons = [s for s in (seasons or results) if s in results]
    variables = [v for v in variables if any(v in results[s] for s in seasons)]
    line_styles, markers = line_styles or {}, markers or {}
    reference_line = {"color": "black", "lw": line_width, "ls": "-", **(reference_line or {})}
    nrows, ncols = len(variables), len(seasons)
    names = {season: (season_names or {}).get(season) or _initial_month(results[season]) for season in seasons}
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * panel_size[0], nrows * panel_size[1]),
                             sharex="col", squeeze=False)
    handles, products, observed = {}, {}, None
    for col, season in enumerate(seasons):
        for row, variable in enumerate(variables):
            ax = axes[row, col]
            datasets = results[season].get(variable, {})
            units, reference = "", None
            for experiment, ds in datasets.items():
                if "member_mean" not in ds:
                    raise KeyError(f"{experiment} {variable}: cache has no member_mean; rebuild the daily caches")
                sel = ds.sel(region=region).where(ds.day <= max_day, drop=True)
                color = colors[experiment]
                if member_display == "range":
                    ax.fill_between(sel.day, sel.member_mean.min("member"), sel.member_mean.max("member"),
                                    color=color, alpha=member_alpha, lw=0)
                elif member_display == "members":
                    for member in sel.member.values:
                        ax.plot(sel.day, sel.member_mean.sel(member=member), color=color, lw=0.6,
                                alpha=member_alpha)
                line, = ax.plot(sel.day, sel.member_mean.mean("member"), color=color, lw=line_width,
                                ls=line_styles.get(experiment, "-"), marker=markers.get(experiment, marker),
                                ms=marker_size)
                name = labels[experiment] if legend_by == "experiment" else labels[experiment].rsplit("-", 1)[0]
                handles.setdefault(season if legend_by == "experiment" else "", {}).setdefault(name, line)
                if reference is None:
                    reference = sel.reference_mean
                    products[variable] = ds.attrs.get("reference", "")
                units = ds.attrs.get("units", "")
            if reference is not None and show_reference:
                observed, = ax.plot(reference.day, reference, zorder=5, **reference_line)
            ax.set_title(f"({chr(97 + row * ncols + col)}) {variable} ({names[season]})", loc="left",
                         fontsize=fontz * 0.9, color="#2b2b29")
            ax.set_ylabel(_UNIT_LABELS.get(units, units), fontsize=fontz * 0.8)
            ax.grid(alpha=grid_alpha)
            ax.set_xlim(0.5, max_day + 0.5)
            ax.tick_params(labelsize=fontz * 0.75)
        axes[-1, col].set_xlabel(xlabel or _days_since_label(results[season]), fontsize=fontz * 0.8)
    extra = []
    if observed is not None:
        extra.append((reference_name or "Observations", observed))
    if member_display == "range":
        extra.append(("Member range",
                      plt.Rectangle((0, 0), 1, 1, color="#9a9a9a", alpha=min(1, 2 * member_alpha), lw=0)))
    _season_legend(fig, handles, extra, height=panel_size[1] * nrows, fontz=fontz, title=title,
                   ncol=legend_ncol, frame=legend_frame, edgecolor=legend_edgecolor, legend_fontz=legend_fontz,
                   extra_last=True, pack=True)
    return fig
