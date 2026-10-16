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
* ``crps_fair``: fair CRPS (ensemble-size adjusted), area mean.

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
    REFERENCES, read_member_ensemble, read_reference, reference_product, region_masks,
)
from util.model_evaluation import experiment_colors

CACHE_VERSION = 1
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
        ds = ds.transpose("region", "subset", "day")
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


def plot_daily_performance(datasets_by_variable, region, *, labels, colors, variables=PANEL_VARIABLES,
                           max_day=15, title="", fontz=16):
    """Fixed rows x 2 figure: daily ensemble-mean RMSE (left) and spread_corr/RMSE (right).

    Panels are labelled "(a) Ensemble-Mean RMSE (PSL)", "(b) Spread/RMSE (PSL)",
    ...; references, the spread correction, and the subsets belong in the
    caption. For ensembles larger than the subset size, the RMSE panels also
    show the mean RMSE of the random equal-size subsets (dashed) and their
    5-95% range (shaded), so every experiment can be compared at the same
    ensemble size. ``fontz`` is the base font size; every label scales from it.
    """
    variables = [v for v in variables if v in datasets_by_variable]
    fig, axes = plt.subplots(len(variables), 2, figsize=(15, 3.6 * len(variables)), sharex=True, squeeze=False)
    handles = {}
    subset_size = None
    for row, variable in enumerate(variables):
        datasets = datasets_by_variable[variable]
        ax_rmse, ax_ratio = axes[row]
        units = ""
        for experiment, ds in datasets.items():
            sel = ds.sel(region=region).where(ds.day <= max_day, drop=True)
            line, = ax_rmse.plot(sel.day, sel.rmse, color=colors[experiment], lw=2.2, marker="o", ms=4)
            ax_ratio.plot(sel.day, sel.spread_ratio, color=colors[experiment], lw=2.2, marker="o", ms=4)
            if "rmse_subset" in sel:
                low, high = sel.rmse_subset.quantile([0.05, 0.95], "subset")
                ax_rmse.fill_between(sel.day, low, high, color=colors[experiment], alpha=0.18, lw=0)
                ax_rmse.plot(sel.day, sel.rmse_subset.mean("subset"), color=colors[experiment], lw=1.6, ls="--")
                subset_size = ds.attrs.get("subset_size", subset_size)
            handles[experiment] = line
            units = ds.attrs.get("units", "")
        ax_ratio.axhline(1.0, color="#6f6e69", lw=1.2, ls="--")
        ax_rmse.set_title(f"({chr(97 + 2 * row)}) Ensemble-Mean RMSE ({variable})", loc="left",
                          fontsize=fontz * 0.9, color="#2b2b29")
        ax_ratio.set_title(f"({chr(98 + 2 * row)}) Spread/RMSE ({variable})", loc="left",
                           fontsize=fontz * 0.9, color="#2b2b29")
        ax_rmse.set_ylabel(_UNIT_LABELS.get(units, units), fontsize=fontz * 0.8)
        for ax in (ax_rmse, ax_ratio):
            ax.grid(alpha=0.3)
            ax.set_xlim(0.5, max_day + 0.5)
            ax.tick_params(labelsize=fontz * 0.75)
    for ax in axes[-1]:
        ax.set_xlabel("Forecast day", fontsize=fontz * 0.8)
    legend_handles = [handles[e] for e in handles]
    legend_labels = [labels[e] for e in handles]
    if subset_size is not None:
        legend_handles.append(plt.Line2D([], [], color="#6f6e69", lw=1.6, ls="--"))
        legend_labels.append(f"{subset_size}-member subsets (5-95% shaded)")
    fig.legend(legend_handles, legend_labels, loc="lower center", ncol=int(np.ceil(len(legend_handles) / 2)),
               frameon=False, fontsize=fontz * 0.85, handlelength=2.5, columnspacing=1.8)
    if title:
        fig.suptitle(title, fontsize=fontz)
    fig.tight_layout(rect=(0, 0.075, 1, 0.98 if title else 1))
    return fig
