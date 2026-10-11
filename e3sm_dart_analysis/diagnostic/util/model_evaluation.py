"""Extended ensemble model-evaluation metrics built on the maintained bias readers."""
from __future__ import annotations

import os
from pathlib import Path
import numpy as np
import xarray as xr
import xskillscore as xs
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

import pandas as pd

from configs.input_paths import OBSERVATION_ROOT, product_dir

from util.bias_metrics import (
    ATM_VARIABLES, LND_VARIABLES, experiment_registry, read_model_ensemble,
    read_reference, derive_metric_dataset, default_model_source, MODEL_SOURCES,
    _REFERENCE_PATTERNS, _period_years,
)

MAP_METRICS = (
    "mean", "reference_mean", "bias", "bias_stddev", "rmse", "spread",
    "crps", "coverage",
)
SCALAR_METRICS = (
    "model_global_mean", "reference_global_mean", "model_rms", "reference_rms",
    "rmse_global", "pattern_correlation", "bias_global", "spread_global",
    "crps_global", "coverage_global", "ssr", "bias_squared", "spread_squared",
    "mse", "sem",
)
OPTIONAL_METRICS = ("rank_histogram",)
METRIC_DESCRIPTIONS = {
    "bias": "ensemble-mean model minus reference",
    "rmse": "member RMSE against the reference",
    "spread": "ensemble standard deviation",
    "crps": "continuous ranked probability score",
    "coverage": "fraction of times the reference lies inside the ensemble envelope",
    "pattern_correlation": "cosine-latitude-weighted spatial correlation of climatologies",
    "ssr": "global spread divided by global RMSE",
    "rank_histogram": "reference rank among ensemble members over all valid samples",
}


def _weights(field):
    return np.cos(np.deg2rad(field.lat)).clip(min=0).broadcast_like(field)


def _weighted_mean(field):
    w = _weights(field)
    return field.weighted(w).mean(("lat", "lon"), skipna=True)


def derive_evaluation_dataset(model, reference, *, experiment, variable, component,
                              run_segment, season, period, n_bootstrap=0,
                              random_seed=42, include_rank_histogram=True,
                              model_source=None):
    """Return cached map and scalar ensemble verification diagnostics."""
    model, reference = xr.align(model, reference, join="inner")
    if model.sizes.get("time", 0) == 0:
        raise ValueError("Model and reference have no overlapping times")
    base = derive_metric_dataset(
        model, reference, experiment=experiment, variable=variable,
        component=component, run_segment=run_segment, season=season,
        period=period, n_bootstrap=n_bootstrap, random_seed=random_seed,
        model_source=model_source,
    )
    ensemble_mean = model.mean("ens", skipna=True)
    model_clim = ensemble_mean.mean("time", skipna=True)
    reference_clim = reference.mean("time", skipna=True)
    difference = ensemble_mean - reference
    ensemble_spread = model.std("ens", skipna=True)
    coverage = ((reference >= model.min("ens", skipna=True)) &
                (reference <= model.max("ens", skipna=True))).astype("float32")
    base["mean_map"] = model_clim
    base["reference_mean_map"] = reference_clim
    base["bias_stddev_map"] = difference.std("time", skipna=True)
    base["coverage_map"] = coverage.mean("time", skipna=True)
    w = _weights(model_clim)
    model_global = model_clim.weighted(w).mean(("lat", "lon"), skipna=True)
    reference_global = reference_clim.weighted(w).mean(("lat", "lon"), skipna=True)
    model_rms = np.sqrt(((model_clim-model_global)**2).weighted(w).mean(("lat","lon"),skipna=True))
    reference_rms = np.sqrt(((reference_clim-reference_global)**2).weighted(w).mean(("lat","lon"),skipna=True))
    rmse_global = xs.rmse(reference_clim, model_clim, dim=["lat","lon"], weights=w, skipna=True)
    pcor = xs.pearson_r(reference_clim-reference_global, model_clim-model_global,
                        dim=["lat","lon"], weights=w, skipna=True)
    bias_global = _weighted_mean(base["bias_map"])
    spread_global = _weighted_mean(base["spread_map"])
    crps_global = _weighted_mean(base["crps_map"])
    coverage_global = _weighted_mean(base["coverage_map"])
    scalars = {
        "model_global_mean": model_global, "reference_global_mean": reference_global,
        "model_rms": model_rms, "reference_rms": reference_rms,
        "rmse_global": rmse_global, "pattern_correlation": pcor,
        "bias_global": bias_global, "spread_global": spread_global,
        "crps_global": crps_global, "coverage_global": coverage_global,
        "ssr": spread_global/rmse_global, "bias_squared": bias_global**2,
        "spread_squared": spread_global**2, "mse": rmse_global**2,
        "sem": spread_global/np.sqrt(model.sizes["ens"]),
    }
    for name,value in scalars.items(): base[name]=value
    if include_rank_histogram:
        ranks=(model < reference).sum("ens").where(reference.notnull())
        values=np.asarray(ranks.values).ravel(); values=values[np.isfinite(values)].astype(int)
        counts=np.bincount(values,minlength=model.sizes["ens"]+1).astype("float64")
        if counts.sum(): counts/=counts.sum()
        base["rank_histogram"] = xr.DataArray(counts,dims="rank",coords={"rank":np.arange(counts.size)})
    base.attrs["available_map_metrics"] = ",".join(MAP_METRICS)
    base.attrs["available_scalar_metrics"] = ",".join(SCALAR_METRICS)
    return base


def evaluation_cache_path(cache_dir, experiment, variable, period):
    return Path(cache_dir)/f"{experiment}_{variable}_{period}_model_evaluation.nc"


def evaluation_cache_is_current(path, *, n_bootstrap=0, model_source="clim"):
    """True when ``path`` exists and was built with these settings."""
    path=Path(path)
    if not path.is_file():
        return False
    with xr.open_dataset(path) as cached:
        return (int(cached.attrs.get("bootstrap_samples",0)) == int(n_bootstrap)
                and cached.attrs.get("model_source","clim") == model_source)


def cache_evaluation(path, *, model=None, reference=None, force_compute=False, **metadata):
    path=Path(path)
    if path.is_file() and not force_compute:
        cached=xr.load_dataset(path)
        requested=int(metadata.get("n_bootstrap",0))
        actual=int(cached.attrs.get("bootstrap_samples",0))
        if requested != actual:
            raise ValueError(f"Cached bootstrap_samples={actual}, requested={requested}; set force_compute=True")
        return cached
    if model is None or reference is None:
        raise ValueError("model and reference are required for a missing cache")
    ds=derive_evaluation_dataset(model,reference,**metadata)
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+".tmp")
    ds.to_netcdf(temporary); os.replace(temporary,path)
    return ds


def plot_evaluation_maps(datasets, metric, variable, *, label_map=None, title="", units=None):
    key=f"{metric}_map"
    absent=[name for name,ds in datasets.items() if key not in ds]
    if absent: raise KeyError(f"{key} missing for {absent}; map metrics are {MAP_METRICS}")
    fields={name:ds[key] for name,ds in datasets.items()}
    finite=np.concatenate([x.values[np.isfinite(x.values)] for x in fields.values()])
    if finite.size==0: raise ValueError(f"No finite values for {variable} {metric}")
    diverging=metric in {"bias","bias_stddev"}
    if diverging:
        limit=float(np.nanpercentile(np.abs(finite),98)); vmin,vmax,cmap=-limit,limit,"RdBu_r"
    else:
        vmin,vmax,cmap=0.0,float(np.nanpercentile(finite,98)),"YlGnBu"
    n=len(fields); cols=2 if n==4 else min(3,n); rows=int(np.ceil(n/cols)); projection=ccrs.PlateCarree()
    fig,axes=plt.subplots(rows,cols,figsize=(5*cols,3.0*rows),subplot_kw={"projection":projection},squeeze=False)
    image=None
    for i,(name,field) in enumerate(fields.items()):
        ax=axes.flat[i]; image=ax.pcolormesh(field.lon,field.lat,field,shading="auto",transform=projection,cmap=cmap,vmin=vmin,vmax=vmax)
        if metric=="bias" and {"bias_ci_lower","bias_ci_upper"} <= set(datasets[name]):
            sig=xr.where((datasets[name].bias_ci_lower>0)|(datasets[name].bias_ci_upper<0),1.0,np.nan)
            ax.contourf(sig.lon,sig.lat,sig,levels=[0.5,1.5],colors="none",hatches=[".."],transform=projection)
        ax.coastlines(linewidth=.6); ax.set_global(); ax.set_title((label_map or {}).get(name,name))
    for i in range(n,rows*cols): axes.flat[i].axis("off")
    if image is not None: fig.colorbar(image,ax=axes.ravel().tolist(),shrink=.75,label=f"{variable} {metric} ({units or metric_units(variable, metric)})")
    if title: fig.suptitle(title)
    return fig


# ---------------------------------------------------------------------------
# Split map / scalar workflow
#
# Both kinds are derived from the same ensemble-verification dataset, so the
# global scalars are always consistent with the maps, but each kind is cached
# and reused on its own (``*_map_metrics.nc`` / ``*_scalar_metrics.nc``).
# ---------------------------------------------------------------------------

METRIC_KINDS = ("map", "scalar")
UNITLESS_METRICS = {"coverage", "coverage_global", "pattern_correlation", "ssr", "rank_histogram"}
SQUARED_METRICS = {"bias_squared", "spread_squared", "mse"}
CATEGORICAL_COLORS = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")


def variable_config(variable):
    """Return the ATM or LND configuration for ``variable``."""
    config = ATM_VARIABLES.get(variable) or LND_VARIABLES.get(variable)
    if config is None:
        raise KeyError(f"Unknown variable {variable!r}")
    return config


def metric_units(variable, metric):
    """Units of ``metric`` for ``variable`` (e.g. 'degC', 'degC^2', '1')."""
    name = metric.removesuffix("_map")
    if name in UNITLESS_METRICS:
        return "1"
    units = variable_config(variable)["units"]
    return f"({units})^2" if name in SQUARED_METRICS else units


def metric_cache_path(cache_dir, kind, experiment, variable, period):
    if kind not in METRIC_KINDS:
        raise ValueError(f"kind must be one of {METRIC_KINDS}")
    return Path(cache_dir) / f"{experiment}_{variable}_{period}_{kind}_metrics.nc"


def select_metric_kind(dataset, kind):
    """Keep the map (``*_map`` and bootstrap CI) or scalar variables of ``dataset``."""
    if kind == "map":
        keep = [name for name in dataset.data_vars if name.endswith("_map") or "_ci_" in name]
    else:
        keep = [name for name in dataset.data_vars if name in SCALAR_METRICS or name in OPTIONAL_METRICS]
    selected = dataset[keep]
    selected.attrs = dict(dataset.attrs, metric_kind=kind)
    return selected


def model_input_files(config, component, variable, period, *, source, max_members=None):
    """Model files the readers open for one experiment (for cache freshness)."""
    count = config["nens"] if max_members is None else min(max_members, config["nens"])
    *_, years = _period_years(period)
    names = [variable] + (["PRECC", "PRECL"] if variable == "PRECT" else [])
    files = []
    for number in range(1, count + 1):
        member = f"EN{number:02d}"
        post = config["path"] / member / "archive/post" / component / "180x360_aave"
        if source == "clim":
            files += sorted((post / "clim").glob("*.nc"))
        else:
            files += [post / "ts" / source / f"{name}.{member}.{year}.nc"
                      for year in years for name in names]
    return [path for path in files if path.is_file()]


def reference_input_files(reference_root, variable, period):
    product = variable_config(variable)["reference"]
    *_, years = _period_years(period)
    return [product_dir(product, reference_root) / "monthly" / _REFERENCE_PATTERNS[product].format(year=year)
            for year in years]


def _cache_settings(*, model_source, n_bootstrap, ensemble_members, lat_range, lon_range,
                    include_rank_histogram):
    return {
        "model_source": model_source,
        "bootstrap_samples": int(n_bootstrap),
        "ensemble_members": int(ensemble_members),
        "lat_range": f"{lat_range[0]},{lat_range[1]}",
        "lon_range": f"{lon_range[0]},{lon_range[1]}",
        "rank_histogram": int(bool(include_rank_histogram)),
    }


def metric_cache_status(path, settings, inputs=()):
    """Return '' when ``path`` can be reused, otherwise why it must be rebuilt."""
    path = Path(path)
    if not path.is_file():
        return "missing"
    newest = max((Path(p).stat().st_mtime for p in inputs), default=0.0)
    if newest > path.stat().st_mtime:
        return "inputs changed"
    with xr.open_dataset(path) as cached:
        for key, value in settings.items():
            if str(cached.attrs.get(key)) != str(value):
                return f"{key}: cached {cached.attrs.get(key)}, requested {value}"
    return ""


def validate_cases(cases, variables, *, data_root, component, run_segment, model_source):
    """Check cases/variables up front and return ``{season: registry}``."""
    known = ATM_VARIABLES if component == "atm" else LND_VARIABLES
    unknown = [v for v in variables if v not in known]
    if unknown:
        raise ValueError(f"{unknown} unavailable for {component}; choose from {list(known)}")
    invalid = {v: s for v, s in (model_source or {}).items() if s not in MODEL_SOURCES}
    if invalid:
        raise ValueError(f"model sources {invalid} invalid; choose from {MODEL_SOURCES}")
    registries = {}
    for season, experiments in cases.items():
        registry = experiment_registry(data_root, season, run_segment)
        missing = [e for e in experiments if e not in registry]
        if missing:
            raise ValueError(f"{season}: {missing} unavailable; choose from {list(registry)}")
        registries[season] = registry
    return registries


def evaluate_cases(kind, cases, variables, *, data_root, cache_dir, component="atm",
                   run_segment="da", model_source=None, n_bootstrap=0, random_seed=42,
                   include_rank_histogram=True, max_members=None, lat_range=(-90, 90),
                   lon_range=(-180, 180), force_compute=False, reference_root=None, log=print):
    """Compute or reuse ``kind`` ('map' or 'scalar') metrics for every case.

    ``cases`` maps a season to its experiments. ``model_source=None`` uses each
    variable's default input (see ``default_model_source``); a string applies
    to every variable. A cache is rebuilt when it is
    missing, older than any model/reference input, or was built with different
    settings (model source, bootstrap, ensemble size, region, rank histogram).
    One failing experiment/variable is logged and skipped, not fatal.

    Returns ``(results, registries, failures)`` with
    ``results[season][variable][experiment] -> Dataset``.
    """
    if kind not in METRIC_KINDS:
        raise ValueError(f"kind must be one of {METRIC_KINDS}")
    sources = {v: model_source or default_model_source(component, run_segment, v) for v in variables}
    reference_root = Path(reference_root or OBSERVATION_ROOT)
    registries = validate_cases(cases, variables, data_root=data_root, component=component,
                                run_segment=run_segment, model_source=sources)
    results, failures, references = {}, [], {}
    for season, experiments in cases.items():
        registry = registries[season]
        for variable in variables:
            source = sources[variable]
            for experiment in experiments:
                config = registry[experiment]
                period = config["period"]
                target = metric_cache_path(cache_dir, kind, experiment, variable, period)
                members = config["nens"] if max_members is None else min(max_members, config["nens"])
                settings = _cache_settings(
                    model_source=source, n_bootstrap=n_bootstrap, ensemble_members=members,
                    lat_range=lat_range, lon_range=lon_range,
                    include_rank_histogram=include_rank_histogram and kind == "scalar")
                try:
                    inputs = (model_input_files(config, component, variable, period,
                                                source=source, max_members=max_members)
                              + reference_input_files(reference_root, variable, period))
                    reason = "forced" if force_compute else metric_cache_status(target, settings, inputs)
                    if not reason:
                        ds = xr.load_dataset(target)
                        log(f"[REUSE] {target.name}")
                    else:
                        model, _ = read_model_ensemble(
                            data_root, experiment, season, run_segment, component, variable,
                            period=period, lat_range=lat_range, lon_range=lon_range,
                            max_members=max_members, source=source)
                        key = (variable, period)
                        if key not in references:
                            references[key] = read_reference(
                                reference_root, component, variable, period,
                                lat_range=lat_range, lon_range=lon_range)
                        full = derive_evaluation_dataset(
                            model, references[key], experiment=experiment, variable=variable,
                            component=component, run_segment=run_segment, season=season,
                            period=period, n_bootstrap=n_bootstrap, random_seed=random_seed,
                            include_rank_histogram=settings["rank_histogram"],
                            model_source=source)
                        ds = select_metric_kind(full, kind)
                        ds.attrs.update(settings, label=config["label"])
                        target.parent.mkdir(parents=True, exist_ok=True)
                        temporary = target.with_suffix(".nc.tmp")
                        ds.to_netcdf(temporary)
                        os.replace(temporary, target)
                        log(f"[WRITE] {target.name} ({reason})")
                except Exception as error:  # keep going; report at the end
                    failures.append((season, variable, experiment, f"{type(error).__name__}: {error}"))
                    log(f"[FAIL] {season} {variable} {experiment}: {type(error).__name__}: {error}")
                    continue
                results.setdefault(season, {}).setdefault(variable, {})[experiment] = ds
    return results, registries, failures


def scalar_summary(results, registries):
    """One row per season/variable/experiment with every cached scalar metric."""
    rows = []
    for season, by_variable in results.items():
        for variable, by_experiment in by_variable.items():
            for experiment, ds in by_experiment.items():
                config = registries[season][experiment]
                row = {"season": season, "variable": variable, "units": variable_config(variable)["units"],
                       "experiment": experiment, "label": config["label"], "period": config["period"],
                       "ensemble_members": int(ds.attrs.get("ensemble_members", config["nens"]))}
                row.update({m: float(ds[m]) for m in SCALAR_METRICS if m in ds})
                rows.append(row)
    return pd.DataFrame(rows)


def experiment_colors(experiments):
    """Fixed categorical color per experiment, in first-seen order (never cycled)."""
    order = list(dict.fromkeys(experiments))
    if len(order) > len(CATEGORICAL_COLORS):
        raise ValueError(f"At most {len(CATEGORICAL_COLORS)} experiments can be colored; facet instead")
    return dict(zip(order, CATEGORICAL_COLORS))


def plot_scalar_metrics(summary, metrics, *, variables=None, title="", colors=None):
    """Small multiples: one row per variable, one column per scalar metric.

    Each panel has one bar per season/experiment; bars are labeled on the axis
    (color only repeats that identity). Reference lines: 0 for bias, 1 for
    spread-skill ratio, the nominal ensemble coverage is not drawn.
    """
    variables = list(variables or dict.fromkeys(summary["variable"]))
    missing = [m for m in metrics if m not in summary]
    if missing:
        raise KeyError(f"{missing} not in summary; scalar metrics are {SCALAR_METRICS}")
    colors = colors or experiment_colors(summary["experiment"])
    fig, axes = plt.subplots(len(variables), len(metrics), squeeze=False,
                             figsize=(2.7 * len(metrics) + 0.6, 2.4 * len(variables) + 0.8))
    for row, variable in enumerate(variables):
        table = summary[summary["variable"] == variable]
        names = [f"{r.experiment}" for r in table.itertuples()]
        x = np.arange(len(table))
        for col, metric in enumerate(metrics):
            ax = axes[row, col]
            ax.bar(x, table[metric], width=0.7, color=[colors[e] for e in table["experiment"]],
                   edgecolor="white", linewidth=1.0, zorder=2)
            if metric == "bias_global" or table[metric].min() < 0:
                ax.axhline(0, color="#555555", linewidth=0.8, zorder=1)
            if metric == "ssr":
                ax.axhline(1, color="#555555", linewidth=0.8, linestyle="--", zorder=1)
            ax.grid(axis="y", color="#e5e5e5", linewidth=0.6, zorder=0)
            ax.spines[["top", "right"]].set_visible(False)
            ax.set_xticks(x)
            ax.set_xticklabels(names if row == len(variables) - 1 else [], rotation=45, ha="right", fontsize=8)
            if row == 0:
                ax.set_title(metric, fontsize=10)
            units = metric_units(variable, metric)
            ax.set_ylabel(f"{variable} ({units})" if col == 0 else f"({units})", fontsize=8)
            ax.tick_params(axis="y", labelsize=8)
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_rank_histograms(datasets, variable, *, label_map=None, title="", colors=None):
    """Rank histogram per experiment; dashed line = flat (well-dispersed) frequency."""
    names = [n for n, ds in datasets.items() if "rank_histogram" in ds]
    if not names:
        raise KeyError("No rank_histogram in these datasets; set include_rank_histogram=True")
    colors = colors or experiment_colors(datasets)
    cols = min(len(names), 3 if len(names) != 4 else 2)
    rows = int(np.ceil(len(names) / cols))
    fig, axes = plt.subplots(rows, cols, squeeze=False, figsize=(3.6 * cols, 2.6 * rows))
    for ax, name in zip(axes.flat, names):
        counts = datasets[name]["rank_histogram"]
        ax.bar(counts["rank"], counts, width=0.85, color=colors[name], zorder=2)
        ax.axhline(1 / counts.size, color="#555555", linewidth=0.8, linestyle="--", zorder=3)
        ax.set_title((label_map or {}).get(name, name), fontsize=10)
        ax.set_xlabel(f"rank of reference among {counts.size - 1} members", fontsize=8)
        ax.set_ylabel("frequency", fontsize=8)
        ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
        ax.grid(axis="y", color="#e5e5e5", linewidth=0.6, zorder=0)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=8)
    for ax in axes.flat[len(names):]:
        ax.axis("off")
    fig.suptitle(title or f"{variable}: rank histograms")
    fig.tight_layout()
    return fig
