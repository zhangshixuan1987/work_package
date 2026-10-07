"""Shared analysis- and forecast-period ensemble bias diagnostics."""

from __future__ import annotations

import os
import re
from pathlib import Path

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import xskillscore as xs

from util.array_stats import nanquantile


ATM_VARIABLES = {
    "FLUT": {"reference": "CERES-OAFlux", "candidates": ["FLUT", "OLR"], "units": "W m-2"},
    "PRECT": {"reference": "GPCP", "candidates": ["PRECT", "pr", "precip"], "units": "mm day-1"},
    "TMQ": {"reference": "ERA5", "candidates": ["TMQ", "tcwv"], "units": "kg m-2"},
    "TS": {"reference": "ERA5", "candidates": ["TS", "skt", "skin_temperature"], "units": "degC"},
    "TREFHT": {"reference": "ERA5", "candidates": ["TREFHT", "t2m"], "units": "degC"},
    "PSL": {"reference": "ERA5", "candidates": ["PSL", "msl"], "units": "hPa"},
}

LND_VARIABLES = {
    "TSA": {"reference": "MODIS_LST", "candidates": ["T2M", "TSA"], "units": "degC"},
    "H2OSOI": {"reference": "ESA_CCI", "candidates": ["H2OSOI"], "units": "m3 m-3", "reducer": "top_5cm"},
    "SOILWATER_10CM": {"reference": "CPC_SOM", "candidates": ["SOILWATER_10CM"], "units": "mm"},
}

METRICS = {
    "bias": "ensemble-mean model minus reference, averaged over time",
    "rmse": "root-mean-square member error over ensemble and time",
    "spread": "ensemble standard deviation averaged over time",
    "crps": "continuous ranked probability score averaged over time",
}

_REFERENCE_PATTERNS = {
    "CERES-OAFlux": "CERES-OAFlux_{year}.nc",
    "GPCP": "PRECT.monthly.{year}.nc",
    "ERA5": "ERA5_analysis_monthly_{year}.nc",
    "MODIS_LST": "MOD11C3.monthly.{year}.nc",
    "ESA_CCI": "H2OSOI.monthly.{year}.nc",
    "CPC_SOM": "SOILWATER_10CM.monthly.{year}.nc",
}

# S0 is the boreal-winter cycle; S1 is the boreal-summer cycle.
_EXPERIMENTS = {
    "S0": {
        "CTRL10-S0": {"nens": 10, "label": "CTRL (ENS=10)", "da": ("CTRLEN10", "201112-201112"), "fc": ("CTRLEN10_15day", "201201-201202")},
        "CAPT10-S0": {"nens": 10, "label": "CAPT (ENS=10)", "da": None, "fc": ("CAPTEN10_15day", "201201-201202")},
        "DART10-S0": {"nens": 10, "label": "EAM-DART (ENS=10)", "da": ("DARTEN10", "201112-201112"), "fc": None},
        "DART20-S0": {"nens": 20, "label": "EAM-DART (ENS=20)", "da": ("DARTEN20", "201112-201112"), "fc": ("DARTEN20_15day", "201201-201202")},
        "DART40-S0": {"nens": 40, "label": "EAM-DART (ENS=40)", "da": ("DARTEN40", "201112-201112"), "fc": ("DARTEN40_15day", "201201-201202")},
    },
    "S1": {
        "CTRL10-S1": {"nens": 10, "label": "CTRL (ENS=10)", "da": None, "fc": ("CTRLEN10s1_15day", "201206-201207")},
        "CAPT10-S1": {"nens": 10, "label": "CAPT (ENS=10)", "da": None, "fc": ("CAPTEN10S1_15day", "201206-201207")},
        "DART40-S1": {"nens": 40, "label": "EAM-DART (ENS=40)", "da": ("DARTEN40S1", "201205-201205"), "fc": ("DARTEN40S1_15day", "201206-201207")},
    },
}


def experiment_registry(data_root, season, run_segment, *, resolution="ne30pg2_r05_IcoswISC30E3r5", machine="compy"):
    """Return experiments available for one season and analysis/forecast segment."""
    if season not in _EXPERIMENTS:
        raise ValueError(f"season must be one of {sorted(_EXPERIMENTS)}")
    if run_segment not in {"da", "fc"}:
        raise ValueError("run_segment must be 'da' (analysis) or 'fc' (forecast)")
    registry = {}
    for name, metadata in _EXPERIMENTS[season].items():
        run = metadata[run_segment]
        if run is None:
            continue
        run_name, period = run
        run_id = f"{run_name}_F20TR_{resolution}_{machine}"
        registry[name] = {
            "label": metadata["label"],
            "nens": metadata["nens"],
            "period": period,
            "run_id": run_id,
            "path": Path(data_root) / run_id,
        }
    return registry


def available_experiments(season, run_segment):
    return list(experiment_registry(".", season, run_segment))


def _period_years(period):
    match = re.fullmatch(r"(\d{4})(\d{2})-(\d{4})(\d{2})", period)
    if not match:
        raise ValueError(f"Invalid period {period!r}; expected YYYYMM-YYYYMM")
    y0, m0, y1, m1 = map(int, match.groups())
    return y0, m0, y1, m1, range(y0, y1 + 1)


def _month_keys(time):
    return (time.dt.year * 100 + time.dt.month).astype("int32")


def _select_period(data, period):
    y0, m0, y1, m1, _ = _period_years(period)
    keys = _month_keys(data.time)
    selected = data.where((keys >= y0 * 100 + m0) & (keys <= y1 * 100 + m1), drop=True)
    return selected.assign_coords(time=_month_keys(selected.time).values)


def _normalize_grid(data):
    rename = {}
    if "latitude" in data.dims:
        rename["latitude"] = "lat"
    if "longitude" in data.dims:
        rename["longitude"] = "lon"
    if rename:
        data = data.rename(rename)
    if "lon" in data.coords and float(data.lon.max()) > 180.0:
        data = data.assign_coords(lon=((data.lon + 180.0) % 360.0) - 180.0).sortby("lon")
    return data


def _subset_region(data, lat_range, lon_range):
    data = _normalize_grid(data)
    lat_slice = slice(*lat_range) if float(data.lat[0]) <= float(data.lat[-1]) else slice(*lat_range[::-1])
    return data.sel(lat=lat_slice, lon=slice(*lon_range))


def _convert_units(data, variable, target_units, *, is_reference=False, reference=None):
    data = data.astype("float32")
    source_units = str(data.attrs.get("units", "")).lower()
    if target_units == "degC":
        if "k" in source_units or (data.size and float(data.mean(skipna=True).compute()) > 100.0):
            data = data - 273.15
    elif target_units == "hPa":
        if "pa" in source_units and "hpa" not in source_units:
            data = data * 0.01
    elif variable == "PRECT" and target_units == "mm day-1":
        if not (is_reference and reference == "GPCP") and ("s-1" in source_units or "s**-1" in source_units or source_units in {"m/s", "m s-1"}):
            data = data * 86400.0 * 1000.0
    data.attrs["units"] = target_units
    return data


def _reduce_land_variable(dataset, variable, reducer):
    if reducer != "top_5cm":
        return dataset[variable]
    field = dataset[variable]
    level_name = next((name for name in ("levgrnd", "levsoi", "level") if name in field.dims), None)
    if level_name is None:
        return field
    if "DZSOI" in dataset:
        thickness = dataset["DZSOI"]
        cumulative = thickness.cumsum(level_name)
        mask = cumulative <= 0.05
        if not bool(mask.any()):
            mask = xr.zeros_like(thickness, dtype=bool)
            mask[{level_name: 0}] = True
        weights = thickness.where(mask, 0.0)
        return (field * weights).sum(level_name) / weights.sum(level_name)
    return field.isel({level_name: 0})


def _open_files(files):
    if not files:
        raise FileNotFoundError("No input files")
    if len(files) == 1:
        return xr.open_dataset(files[0], chunks={})
    return xr.open_mfdataset(files, combine="by_coords", parallel=True, chunks={})


def read_model_ensemble(data_root, experiment, season, run_segment, component, variable, period=None, lat_range=(-90, 90), lon_range=(-180, 180), max_members=None):
    """Read monthly climatology files beneath each EN## member directory."""
    registry = experiment_registry(data_root, season, run_segment)
    if experiment not in registry:
        raise ValueError(f"{experiment!r} is unavailable; choose from {list(registry)}")
    config = registry[experiment]
    period = period or config["period"]
    variable_config = (ATM_VARIABLES if component == "atm" else LND_VARIABLES)[variable]
    count = config["nens"] if max_members is None else min(max_members, config["nens"])
    members = []
    opened = []
    try:
        for member_number in range(1, count + 1):
            members.append(_lazy_member_field(
                config, member_number, component, variable, variable_config,
                period, lat_range, lon_range, opened,
            ))
        # One dask graph over all members/months so files are read in parallel.
        ensemble = xr.concat(members, dim="ens").compute()
    finally:
        for dataset in opened:
            dataset.close()
    return ensemble, config


def _lazy_member_field(config, member_number, component, variable, variable_config, period, lat_range, lon_range, opened):
    """Build one member's monthly field lazily; open datasets are appended to ``opened``."""
    member = f"EN{member_number:02d}"
    clim = config["path"] / member / "archive/post" / component / "180x360_aave/clim"
    y0, m0, y1, m1, _ = _period_years(period)
    lower, upper = y0 * 100 + m0, y1 * 100 + m1
    candidates = sorted(clim.glob("*.nc"))
    files = []
    for candidate in candidates:
        match = re.search(r"(\d{4})-(\d{2})\.nc$", candidate.name)
        if match and lower <= int(match.group(1) + match.group(2)) <= upper:
            files.append((candidate, int(match.group(1) + match.group(2))))
    if not files:
        raise FileNotFoundError(f"No monthly files for {period} in {clim}")

    monthly_fields = []
    for monthly_file, month_key in files:
        dataset = xr.open_dataset(monthly_file, chunks={})
        opened.append(dataset)
        if variable not in dataset:
            if variable == "PRECT" and {"PRECC", "PRECL"} <= set(dataset):
                field = dataset["PRECC"] + dataset["PRECL"]
            else:
                raise KeyError(f"{variable} not found in {monthly_file}")
        else:
            field = _reduce_land_variable(
                dataset, variable, variable_config.get("reducer")
            )
        if "time" in field.dims:
            field = field.mean("time", skipna=True)
        field = _subset_region(field, lat_range, lon_range)
        field = _convert_units(field, variable, variable_config["units"])
        monthly_fields.append(field.reset_coords(drop=True).expand_dims(time=[month_key]))
    return xr.concat(monthly_fields, dim="time").expand_dims(ens=[member])


def read_reference(reference_root, component, variable, period, lat_range=(-90, 90), lon_range=(-180, 180)):
    """Read the configured monthly reference product for a variable."""
    config = (ATM_VARIABLES if component == "atm" else LND_VARIABLES)[variable]
    product = config["reference"]
    _, _, _, _, years = _period_years(period)
    pattern = _REFERENCE_PATTERNS[product]
    files = [Path(reference_root) / product / "monthly" / pattern.format(year=year) for year in years]
    missing = [path for path in files if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing reference files: " + ", ".join(map(str, missing)))
    dataset = _open_files(files)
    try:
        name = next((candidate for candidate in config["candidates"] if candidate in dataset), None)
        if name is None:
            raise KeyError(f"None of {config['candidates']} found in {files[0]}")
        field = _select_period(dataset[name], period)
        field = _subset_region(field, lat_range, lon_range)
        field = _convert_units(field, variable, config["units"], is_reference=True, reference=product).load()
        field = field.reset_coords(drop=True)
    finally:
        dataset.close()
    return field


def derive_metric_dataset(model, reference, *, experiment, variable, component, run_segment, season, period, n_bootstrap=0, random_seed=42):
    """Compute bias, RMSE, spread, CRPS, and optional ensemble-bootstrap CIs."""
    model, reference = xr.align(model, reference, join="inner")
    if model.sizes.get("time", 0) == 0:
        raise ValueError("Model and reference have no overlapping months")

    def maps(sample):
        mean = sample.mean("ens", skipna=True)
        error = sample - reference
        bias = (mean - reference).mean("time", skipna=True)
        rmse = np.sqrt((error ** 2).mean(("ens", "time"), skipna=True))
        spread = sample.std("ens", skipna=True).mean("time", skipna=True)
        crps = xs.crps_ensemble(reference, sample, member_dim="ens", dim=[]).mean("time", skipna=True)
        return {"bias": bias, "rmse": rmse, "spread": spread, "crps": crps}

    results = maps(model)
    dataset = xr.Dataset({f"{name}_map": field for name, field in results.items()})
    if n_bootstrap:
        rng = np.random.default_rng(random_seed)
        samples = {name: [] for name in METRICS}
        for _ in range(n_bootstrap):
            indices = rng.integers(0, model.sizes["ens"], size=model.sizes["ens"])
            boot = maps(model.isel(ens=indices))
            for name, field in boot.items():
                samples[name].append(field)
        for name, fields in samples.items():
            stack = xr.concat(fields, dim="bootstrap")
            dataset[f"{name}_ci_lower"] = nanquantile(stack, 0.025, "bootstrap").reset_coords(drop=True)
            dataset[f"{name}_ci_upper"] = nanquantile(stack, 0.975, "bootstrap").reset_coords(drop=True)

    dataset.attrs.update({
        "experiment": experiment,
        "variable": variable,
        "component": component,
        "run_segment": run_segment,
        "season": season,
        "period": period,
        "reference": (ATM_VARIABLES if component == "atm" else LND_VARIABLES)[variable]["reference"],
        "ensemble_members": model.sizes["ens"],
        "bootstrap_samples": n_bootstrap,
    })
    return dataset


def cache_path(cache_dir, experiment, variable, period):
    return Path(cache_dir) / f"{experiment}_{variable}_{period}_ensemble_bias_metrics.nc"


def cache_bias_metrics(path, *, model=None, reference=None, force_compute=False, **metadata):
    """Reuse a cached dataset unless force_compute is true."""
    path = Path(path)
    if path.is_file() and not force_compute:
        return xr.load_dataset(path)
    if model is None or reference is None:
        raise ValueError("model and reference are required to derive a missing cache")
    dataset = derive_metric_dataset(model, reference, **metadata)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    dataset.to_netcdf(temporary)
    os.replace(temporary, path)
    return dataset


def plot_metric_maps(datasets, metric, variable, label_map=None, title=""):
    """Plot one cached metric across experiments with a shared color scale."""
    key = f"{metric}_map"
    fields = {name: dataset[key] for name, dataset in datasets.items()}
    values = np.concatenate([field.values[np.isfinite(field.values)] for field in fields.values()])
    if values.size == 0:
        raise ValueError(f"No finite values for {variable} {metric}")
    if metric == "bias":
        limit = float(np.nanpercentile(np.abs(values), 98))
        vmin, vmax, cmap = -limit, limit, "RdBu_r"
    else:
        vmin, vmax, cmap = 0.0, float(np.nanpercentile(values, 98)), "YlGnBu"
    count = len(fields)
    columns = min(3, count)
    rows = int(np.ceil(count / columns))
    projection = ccrs.PlateCarree()
    figure, axes = plt.subplots(rows, columns, figsize=(5 * columns, 3.8 * rows), subplot_kw={"projection": projection}, squeeze=False)
    image = None
    for index, (experiment, field) in enumerate(fields.items()):
        axis = axes.flat[index]
        image = axis.pcolormesh(field.lon, field.lat, field, shading="auto", transform=projection, cmap=cmap, vmin=vmin, vmax=vmax)
        if metric == "bias":
            lower_key, upper_key = "bias_ci_lower", "bias_ci_upper"
            dataset = datasets[experiment]
            if lower_key in dataset and upper_key in dataset:
                significant = xr.where((dataset[lower_key] > 0) | (dataset[upper_key] < 0), 1.0, np.nan)
                axis.contourf(significant.lon, significant.lat, significant, levels=[0.5, 1.5], colors="none", hatches=[".."], transform=projection)
        axis.coastlines(linewidth=0.6)
        axis.set_global()
        axis.set_title(label_map.get(experiment, experiment) if label_map else experiment)
    for index in range(count, rows * columns):
        axes.flat[index].axis("off")
    if image is not None:
        figure.colorbar(image, ax=axes.ravel().tolist(), shrink=0.75, label=f"{variable} {metric}")
    if title:
        figure.suptitle(title)
    return figure


def _validated_registry(data_root, season, run_segment, component, experiments, variables):
    variable_registry = ATM_VARIABLES if component == "atm" else LND_VARIABLES
    unknown_variables = sorted(set(variables) - set(variable_registry))
    if unknown_variables:
        raise ValueError(f"Unsupported {component} variables: {unknown_variables}")
    registry = experiment_registry(data_root, season, run_segment)
    unknown_experiments = sorted(set(experiments) - set(registry))
    if unknown_experiments:
        raise ValueError(
            f"Experiments unavailable for {season}/{run_segment}: "
            f"{unknown_experiments}; choose from {list(registry)}"
        )
    return registry


def run_bias_metrics(
    *, data_root, reference_root, cache_dir, season, run_segment, component,
    experiments, variables, lat_range=(-90, 90), lon_range=(-180, 180),
    n_bootstrap=0, random_seed=42, max_members=None, force_compute=False,
):
    """Compute and cache metrics for every experiment/variable not yet cached."""
    registry = _validated_registry(data_root, season, run_segment, component, experiments, variables)
    Path(cache_dir).mkdir(parents=True, exist_ok=True)
    for variable in variables:
        references = {}
        try:
            for experiment in experiments:
                period = registry[experiment]["period"]
                output = cache_path(cache_dir, experiment, variable, period)
                if output.is_file() and not force_compute:
                    print(f"Reusing {output.name}")
                    continue
                if period not in references:
                    references[period] = read_reference(
                        reference_root, component, variable, period, lat_range, lon_range,
                    )
                model, _ = read_model_ensemble(
                    data_root, experiment, season, run_segment, component, variable,
                    period=period, lat_range=lat_range, lon_range=lon_range,
                    max_members=max_members,
                )
                try:
                    cache_bias_metrics(
                        output, model=model, reference=references[period],
                        force_compute=force_compute, experiment=experiment,
                        variable=variable, component=component, run_segment=run_segment,
                        season=season, period=period, n_bootstrap=n_bootstrap,
                        random_seed=random_seed,
                    ).close()
                finally:
                    model.close()
        finally:
            for reference in references.values():
                reference.close()


def plot_cached_bias_metrics(
    *, data_root, cache_dir, figure_dir, season, run_segment, component,
    experiments, variables, metrics, dpi=200,
):
    """Plot cached metric maps; return the list of written figure paths."""
    unknown_metrics = sorted(set(metrics) - set(METRICS))
    if unknown_metrics:
        raise ValueError(f"Unknown metrics: {unknown_metrics}")
    registry = _validated_registry(data_root, season, run_segment, component, experiments, variables)
    figure_dir = Path(figure_dir)
    figure_dir.mkdir(parents=True, exist_ok=True)
    label_map = {name: registry[name]["label"] for name in experiments}
    written = []
    for variable in variables:
        paths = {
            experiment: cache_path(cache_dir, experiment, variable, registry[experiment]["period"])
            for experiment in experiments
        }
        missing = [str(path) for path in paths.values() if not path.is_file()]
        if missing:
            raise FileNotFoundError("Missing cached bias metrics:\n" + "\n".join(missing))
        cached = {experiment: xr.load_dataset(path) for experiment, path in paths.items()}
        try:
            for metric in metrics:
                figure = plot_metric_maps(
                    cached, metric, variable, label_map=label_map,
                    title=f"{variable} {metric.upper()} | {season} {run_segment.upper()}",
                )
                output = figure_dir / f"{variable}_{metric}_{season}_{run_segment}.png"
                figure.savefig(output, dpi=dpi, bbox_inches="tight")
                plt.close(figure)
                written.append(output)
        finally:
            for dataset in cached.values():
                dataset.close()
    return written
