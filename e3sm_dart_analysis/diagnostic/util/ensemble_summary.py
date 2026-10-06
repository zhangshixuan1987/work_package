"""Cached spatial range summaries for local, Zarr, or Intake-ESM ensembles."""
from __future__ import annotations

import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


SUMMARY_STATISTICS = ("spatial_min", "spatial_mean", "spatial_max")
CACHE_FIELDS = (
    "spatial_min",
    "spatial_mean",
    "spatial_max",
    "ensemble_min",
    "ensemble_mean",
    "ensemble_max",
)


def open_ensemble_source(
    source,
    *,
    source_type="zarr",
    variable=None,
    dataset_key=None,
    storage_options=None,
    open_kwargs=None,
    member_dim="member_id",
):
    """Open a NetCDF file/glob, a Zarr store, or one Intake-ESM catalog entry."""
    storage_options = storage_options or {}
    open_kwargs = open_kwargs or {}
    if source_type == "zarr":
        return xr.open_zarr(source, storage_options=storage_options, **open_kwargs)
    if source_type == "member_netcdf":
        import glob
        matches = sorted(glob.glob(str(source)))
        if not matches:
            raise FileNotFoundError(source)
        labels = []
        for path in matches:
            match = re.search(r"EN\d+", Path(path).name)
            if match is None:
                raise ValueError(f"Cannot parse EN## member from {path}")
            labels.append(match.group(0))
        dataset = xr.open_mfdataset(
            matches,
            combine="nested",
            concat_dim=member_dim,
            **open_kwargs,
        )
        return dataset.assign_coords({member_dim: labels})
    if source_type == "netcdf":
        if any(char in str(source) for char in "*?[]"):
            import glob
            matches = sorted(glob.glob(str(source)))
            if not matches:
                raise FileNotFoundError(source)
            return xr.open_mfdataset(matches, combine="by_coords", **open_kwargs)
        return xr.open_dataset(source, **open_kwargs)
    if source_type == "intake_esm":
        if variable is None:
            raise ValueError("variable is required for an Intake-ESM source")
        try:
            import intake
        except ImportError as error:
            raise ImportError("intake-esm is required only for source_type=intake_esm") from error
        catalog = intake.open_esm_datastore(source)
        datasets = catalog.search(variable=variable).to_dataset_dict(
            xarray_open_kwargs={"consolidated": True, **open_kwargs},
            storage_options=storage_options,
        )
        if not datasets:
            raise FileNotFoundError(f"No Intake-ESM entry found for {variable}")
        key = dataset_key or sorted(datasets)[0]
        if key not in datasets:
            raise KeyError(f"{key!r} not found; choose one of {sorted(datasets)}")
        return datasets[key]
    raise ValueError("source_type must be member_netcdf, netcdf, zarr, or intake_esm")


def _select_level(field, level_dim=None, level=None):
    if level_dim is None:
        level_dim = next((name for name in ("lev", "plev", "level") if name in field.dims), None)
    if level_dim is None:
        return field
    if level is None:
        return field.isel({level_dim: -1}, drop=True)
    if isinstance(level, (int, np.integer)):
        return field.isel({level_dim: int(level)}, drop=True)
    return field.sel({level_dim: level}, method="nearest").squeeze(drop=True)


def compute_ensemble_summary(
    dataset,
    variable,
    *,
    member_dim="member_id",
    time_dim="time",
    spatial_dims=None,
    level_dim=None,
    level=None,
):
    """Reduce every ensemble member to spatial min/mean/max time series."""
    if variable not in dataset:
        raise KeyError(f"{variable!r} not found; choose one of {list(dataset.data_vars)}")
    field = _select_level(dataset[variable], level_dim=level_dim, level=level)
    if member_dim not in field.dims or time_dim not in field.dims:
        raise ValueError(f"{variable} must contain {member_dim!r} and {time_dim!r}; dims={field.dims}")
    if spatial_dims is None:
        spatial_dims = tuple(dim for dim in field.dims if dim not in {member_dim, time_dim})
    else:
        spatial_dims = tuple(spatial_dims)
    if not spatial_dims:
        raise ValueError("At least one spatial dimension is required")
    missing = [dim for dim in spatial_dims if dim not in field.dims]
    if missing:
        raise ValueError(f"Missing spatial dimensions: {missing}")

    reduced = xr.Dataset(
        {
            "spatial_min": field.min(spatial_dims, skipna=True),
            "spatial_mean": field.mean(spatial_dims, skipna=True),
            "spatial_max": field.max(spatial_dims, skipna=True),
        }
    )
    stacked = reduced.to_array("statistic")
    reduced["ensemble_min"] = stacked.min(member_dim, skipna=True)
    reduced["ensemble_mean"] = stacked.mean(member_dim, skipna=True)
    reduced["ensemble_max"] = stacked.max(member_dim, skipna=True)
    units = str(field.attrs.get("units", ""))
    for name in SUMMARY_STATISTICS:
        reduced[name].attrs["units"] = units
    reduced.attrs.update(
        variable=variable,
        member_dim=member_dim,
        time_dim=time_dim,
        spatial_dims=",".join(spatial_dims),
        level_dim=level_dim or "",
        level="surface" if level is None else str(level),
        units=units,
    )
    return reduced


def cache_ensemble_summary(
    path,
    *,
    dataset=None,
    force_compute=False,
    variable=None,
    member_dim="member_id",
    time_dim="time",
    spatial_dims=None,
    level_dim=None,
    level=None,
):
    path = Path(path)
    if path.is_file() and not force_compute:
        return xr.load_dataset(path)
    if dataset is None or variable is None:
        raise ValueError("dataset and variable are required for a missing cache")
    result = compute_ensemble_summary(
        dataset,
        variable,
        member_dim=member_dim,
        time_dim=time_dim,
        spatial_dims=spatial_dims,
        level_dim=level_dim,
        level=level,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    result.to_netcdf(tmp)
    os.replace(tmp, path)
    return result


def plot_ensemble_summary(dataset, *, title=None):
    """Plot member spatial minima, means, and maxima with ensemble envelopes."""
    member_dim = dataset.attrs["member_dim"]
    time_dim = dataset.attrs["time_dim"]
    time = dataset[time_dim]
    units = dataset.attrs.get("units", "")
    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    labels = {
        "spatial_max": "Ensemble-member spatial maxima",
        "spatial_mean": "Ensemble-member spatial means",
        "spatial_min": "Ensemble-member spatial minima",
    }
    stat_values = list(dataset.statistic.values.astype(str))
    for axis, name in zip(axes, ("spatial_max", "spatial_mean", "spatial_min")):
        field = dataset[name]
        for member in field[member_dim].values:
            axis.plot(time, field.sel({member_dim: member}), color="tab:red", alpha=0.15, linewidth=0.7)
        index = stat_values.index(name)
        lower = dataset.ensemble_min.isel(statistic=index)
        upper = dataset.ensemble_max.isel(statistic=index)
        mean = dataset.ensemble_mean.isel(statistic=index)
        axis.fill_between(time, lower, upper, color="0.85", zorder=0)
        axis.plot(time, mean, color="black", linewidth=1.5, label="ensemble mean")
        axis.set_ylabel(units)
        axis.set_title(labels[name])
        axis.grid(alpha=0.25)
    axes[0].legend()
    axes[-1].set_xlabel("Time")
    fig.suptitle(title or dataset.attrs["variable"])
    fig.tight_layout()
    return fig


def plot_ensemble_summary_comparison(
    datasets,
    *,
    labels=None,
    colors=None,
    title=None,
    show_members=True,
):
    """Compare spatial-range summaries from two or more experiments."""
    if len(datasets) < 2:
        raise ValueError("At least two ensemble summaries are required for comparison")

    labels = labels or {}
    colors = colors or {}
    default_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    statistics = ("spatial_max", "spatial_mean", "spatial_min")
    panel_titles = {
        "spatial_max": "Ensemble-member spatial maxima",
        "spatial_mean": "Ensemble-member spatial means",
        "spatial_min": "Ensemble-member spatial minima",
    }

    variables = {
        str(dataset.attrs.get("variable", "")) for dataset in datasets.values()
    }
    if len(variables) != 1:
        raise ValueError(
            f"All summaries must use one variable; found {sorted(variables)}"
        )
    units = {str(dataset.attrs.get("units", "")) for dataset in datasets.values()}
    if len(units) != 1:
        raise ValueError(f"All summaries must use one unit; found {sorted(units)}")

    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    for experiment_index, (experiment, dataset) in enumerate(datasets.items()):
        member_dim = dataset.attrs["member_dim"]
        time_dim = dataset.attrs["time_dim"]
        time = dataset[time_dim]
        statistic_values = list(dataset.statistic.values.astype(str))
        label = labels.get(experiment, experiment)
        color = colors.get(
            experiment,
            default_colors[experiment_index % len(default_colors)],
        )

        for axis, statistic in zip(axes, statistics):
            field = dataset[statistic]
            if show_members:
                for member in field[member_dim].values:
                    axis.plot(
                        time,
                        field.sel({member_dim: member}),
                        color=color,
                        alpha=0.08,
                        linewidth=0.6,
                    )
            statistic_index = statistic_values.index(statistic)
            lower = dataset.ensemble_min.isel(statistic=statistic_index)
            upper = dataset.ensemble_max.isel(statistic=statistic_index)
            center = dataset.ensemble_mean.isel(statistic=statistic_index)
            axis.fill_between(time, lower, upper, color=color, alpha=0.14)
            axis.plot(time, center, color=color, linewidth=1.8, label=label)

    unit = next(iter(units))
    for axis, statistic in zip(axes, statistics):
        axis.set_ylabel(unit)
        axis.set_title(panel_titles[statistic])
        axis.grid(alpha=0.25)
    axes[0].legend()
    axes[-1].set_xlabel("Time")
    variable = next(iter(variables))
    fig.suptitle(title or f"{variable} ensemble spatial-range comparison")
    fig.tight_layout()
    return fig
