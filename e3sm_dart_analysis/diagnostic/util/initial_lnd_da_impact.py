"""Diagnose how atmospheric DA changes the initialized land state."""

from __future__ import annotations

from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
import xarray as xr


def _required(dataset: xr.Dataset, names: tuple[str, ...]) -> None:
    missing = [name for name in names if name not in dataset]
    if missing:
        raise KeyError(f"Dataset is missing required variables: {', '.join(missing)}")


def _weighted_mean(field: xr.DataArray, weights: xr.DataArray) -> float:
    valid = np.isfinite(field)
    denominator = weights.where(valid).sum(skipna=True)
    if float(denominator) <= 0.0:
        return np.nan
    return float((field * weights).where(valid).sum(skipna=True) / denominator)


def derive_land_da_impact(
    dataset: xr.Dataset,
    *,
    baseline: str = "CTRL10-S0",
    experiments: tuple[str, ...] = ("DART10-S0", "DART20-S0", "DART40-S0"),
    time_index: int = 0,
) -> xr.Dataset:
    """Derive land-state improvements relative to a no-DA baseline.

    Positive ``absolute_error_reduction`` means that the experiment is closer
    to the observed land state than the baseline. All summary statistics use
    the common finite mask for each experiment/baseline pair.
    """
    _required(dataset, ("obs", "model_mean"))
    if "exp" not in dataset.dims:
        raise ValueError("Expected an 'exp' dimension in the input dataset.")

    available = {str(value) for value in dataset.exp.values}
    requested = {baseline, *experiments}
    missing = sorted(requested - available)
    if missing:
        raise KeyError(
            f"Experiments missing from cache: {missing}; available: {sorted(available)}"
        )

    selected = dataset.isel(time=time_index, drop=True) if "time" in dataset.dims else dataset
    observation = selected["obs"].sel(exp=baseline, drop=True)
    baseline_mean = selected["model_mean"].sel(exp=baseline, drop=True)
    baseline_error = baseline_mean - observation
    latitude = observation["lat"]
    latitude_weights = np.cos(np.deg2rad(latitude))

    outputs = []
    for experiment in experiments:
        experiment_observation = selected["obs"].sel(exp=experiment, drop=True)
        experiment_mean = selected["model_mean"].sel(exp=experiment, drop=True)
        baseline_aligned, experiment_aligned, observation_aligned, experiment_obs_aligned = xr.align(
            baseline_mean,
            experiment_mean,
            observation,
            experiment_observation,
            join="inner",
        )
        # Observation is duplicated along exp in the source cache. Requiring a
        # shared mask prevents missing-data differences from becoming a DA signal.
        valid = (
            np.isfinite(baseline_aligned)
            & np.isfinite(experiment_aligned)
            & np.isfinite(observation_aligned)
            & np.isfinite(experiment_obs_aligned)
        )
        baseline_error = (baseline_aligned - observation_aligned).where(valid)
        experiment_error = (experiment_aligned - observation_aligned).where(valid)
        benefit = np.abs(baseline_error) - np.abs(experiment_error)

        zonal_baseline_rmse = np.sqrt((baseline_error ** 2).mean("lon", skipna=True))
        zonal_experiment_rmse = np.sqrt((experiment_error ** 2).mean("lon", skipna=True))
        zonal_reduction = zonal_baseline_rmse - zonal_experiment_rmse

        weights = latitude_weights.broadcast_like(baseline_error)
        baseline_rmse = np.sqrt(_weighted_mean(baseline_error ** 2, weights))
        experiment_rmse = np.sqrt(_weighted_mean(experiment_error ** 2, weights))
        skill = 100.0 * (baseline_rmse - experiment_rmse) / baseline_rmse
        improved_fraction = 100.0 * _weighted_mean((benefit > 0).astype(float), weights)

        output = xr.Dataset(
            {
                "absolute_error_reduction": benefit,
                "zonal_rmse_reduction": zonal_reduction,
                "baseline_rmse": xr.DataArray(baseline_rmse),
                "experiment_rmse": xr.DataArray(experiment_rmse),
                "rmse_skill": xr.DataArray(skill),
                "improved_land_fraction": xr.DataArray(improved_fraction),
            }
        ).expand_dims(experiment=[experiment])
        outputs.append(output)

    result = xr.concat(outputs, dim="experiment")
    units = selected.attrs.get("units", selected["model_mean"].attrs.get("units", ""))
    result["absolute_error_reduction"].attrs.update(
        long_name="reduction in absolute land-state error",
        units=units,
        interpretation="positive values indicate improvement relative to baseline",
    )
    result["zonal_rmse_reduction"].attrs.update(
        long_name="zonal RMSE reduction relative to baseline", units=units
    )
    result["baseline_rmse"].attrs["units"] = units
    result["experiment_rmse"].attrs["units"] = units
    result["rmse_skill"].attrs["units"] = "%"
    result["improved_land_fraction"].attrs["units"] = "%"
    result.attrs.update(
        baseline=baseline,
        target_date=str(selected.attrs.get("target_date", "")),
        variable=str(selected.attrs.get("out_var", selected.attrs.get("model_var", ""))),
        depth_cm=float(selected.attrs.get("depth_cm", np.nan)),
        units=units,
        caveat=(
            "Comparisons with a different ensemble size combine atmospheric-DA "
            "and ensemble-size effects."
        ),
    )
    return result


def derive_land_da_impact_from_evaluations(
    datasets: dict[str, xr.Dataset],
    *,
    baseline: str = "CTRL10-S0",
    experiments: tuple[str, ...] = ("DART10-S0", "DART20-S0", "DART40-S0"),
) -> xr.Dataset:
    """Derive impact from per-experiment model-evaluation caches.

    The comparison uses each experiment's ensemble-mean climatology and the
    corresponding reference climatology. The same-size DART10/CTRL10 pair is
    therefore the clean atmospheric-DA attribution; DART20/40 additionally
    show sensitivity to ensemble size.
    """
    requested = (baseline, *experiments)
    missing = [name for name in requested if name not in datasets]
    if missing:
        raise KeyError(f"Evaluation datasets missing for: {missing}")
    for name in requested:
        _required(datasets[name], ("mean_map", "reference_mean_map"))

    means = xr.concat(
        [datasets[name]["mean_map"] for name in requested],
        dim=pd.Index(requested, name="exp"),
        join="inner",
    )
    observations = xr.concat(
        [datasets[name]["reference_mean_map"] for name in requested],
        dim=pd.Index(requested, name="exp"),
        join="inner",
    )
    first = datasets[baseline]
    combined = xr.Dataset({"model_mean": means, "obs": observations})
    combined.attrs.update(
        target_date=f"{first.attrs.get('period', '')} mean",
        out_var=first.attrs.get("variable", ""),
        units=first["mean_map"].attrs.get("units", ""),
    )
    return derive_land_da_impact(
        combined, baseline=baseline, experiments=experiments
    )


def impact_summary(impact: xr.Dataset) -> pd.DataFrame:
    """Return one compact row per experiment for display or CSV output."""
    rows = []
    for experiment in impact.experiment.values:
        selected = impact.sel(experiment=experiment)
        rows.append(
            {
                "experiment": str(experiment),
                "baseline": impact.attrs["baseline"],
                "baseline_rmse": float(selected["baseline_rmse"]),
                "experiment_rmse": float(selected["experiment_rmse"]),
                "rmse_skill_percent": float(selected["rmse_skill"]),
                "improved_land_fraction_percent": float(
                    selected["improved_land_fraction"]
                ),
            }
        )
    return pd.DataFrame(rows)


def plot_land_da_impact(
    impact: xr.Dataset,
    *,
    primary_experiment: str | None = None,
    experiment_labels: dict[str, str] | None = None,
    colors: dict[str, str] | None = None,
    percentile: float = 98.0,
    figsize: tuple[float, float] = (16.0, 5.2),
):
    """Plot where, by latitude, and by how much atmospheric DA helps land state."""
    experiments = [str(value) for value in impact.experiment.values]
    if not experiments:
        raise ValueError("At least one experiment is required.")
    primary_experiment = primary_experiment or experiments[-1]
    if primary_experiment not in experiments:
        raise KeyError(f"Primary experiment {primary_experiment!r} is unavailable.")
    labels = experiment_labels or {}
    palette = colors or {}

    figure = plt.figure(figsize=figsize, layout="constrained")
    grid = figure.add_gridspec(1, 3, width_ratios=(1.55, 0.75, 0.90))
    map_axis = figure.add_subplot(grid[0, 0], projection=ccrs.Robinson())
    zonal_axis = figure.add_subplot(grid[0, 1])
    summary_axis = figure.add_subplot(grid[0, 2])

    field = impact["absolute_error_reduction"].sel(experiment=primary_experiment)
    finite = np.abs(np.asarray(field.values)[np.isfinite(field.values)])
    limit = float(np.percentile(finite, percentile)) if finite.size else 1.0
    if not np.isfinite(limit) or limit <= 0.0:
        limit = 1.0
    image = map_axis.pcolormesh(
        field.lon,
        field.lat,
        field,
        transform=ccrs.PlateCarree(),
        cmap="RdBu",
        norm=TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit),
        shading="auto",
        rasterized=True,
    )
    map_axis.set_global()
    map_axis.add_feature(cfeature.LAND, facecolor="none", edgecolor="none")
    map_axis.coastlines(linewidth=0.55)
    colorbar = figure.colorbar(image, ax=map_axis, orientation="horizontal", pad=0.04, shrink=0.88)
    units = impact.attrs.get("units", "")
    colorbar.set_label(f"Absolute-error reduction ({units}); blue = improvement")
    map_axis.set_title(f"a  Where DA helps: {labels.get(primary_experiment, primary_experiment)}", loc="left")

    for index, experiment in enumerate(experiments):
        line = impact["zonal_rmse_reduction"].sel(experiment=experiment)
        zonal_axis.plot(
            line,
            line.lat,
            linewidth=2.0,
            label=labels.get(experiment, experiment),
            color=palette.get(experiment, f"C{index}"),
        )
    zonal_axis.axvline(0.0, color="0.25", linewidth=0.8)
    zonal_axis.axhline(0.0, color="0.75", linewidth=0.6)
    zonal_axis.set_ylim(-90, 90)
    zonal_axis.set_yticks((-90, -60, -30, 0, 30, 60, 90))
    zonal_axis.grid(alpha=0.22)
    zonal_axis.set_xlabel(f"Zonal RMSE reduction ({units})")
    zonal_axis.set_ylabel("Latitude")
    zonal_axis.set_title("b  Latitudinal response", loc="left")
    zonal_axis.legend(frameon=False, fontsize=9)

    x = np.arange(len(experiments), dtype=float)
    width = 0.36
    skill = impact["rmse_skill"].values
    fraction = impact["improved_land_fraction"].values
    bars_skill = summary_axis.bar(
        x - width / 2,
        skill,
        width,
        label="Global RMSE reduction",
        color=[palette.get(name, f"C{i}") for i, name in enumerate(experiments)],
    )
    bars_fraction = summary_axis.bar(
        x + width / 2,
        fraction,
        width,
        label="Land area improved",
        color="0.72",
        edgecolor="0.35",
    )
    summary_axis.axhline(0.0, color="0.25", linewidth=0.8)
    summary_axis.set_xticks(x, [labels.get(name, name) for name in experiments], rotation=18, ha="right")
    summary_axis.set_ylabel("Percent")
    summary_axis.set_title("c  Global benefit", loc="left")
    summary_axis.grid(axis="y", alpha=0.22)
    summary_axis.legend(frameon=False, fontsize=9)
    summary_axis.bar_label(bars_skill, fmt="%.1f%%", padding=2, fontsize=9)
    summary_axis.bar_label(bars_fraction, fmt="%.0f%%", padding=2, fontsize=9)

    depth = impact.attrs.get("depth_cm", np.nan)
    depth_label = f"top {depth:g} cm " if np.isfinite(depth) else ""
    date = impact.attrs.get("target_date", "")
    baseline = labels.get(impact.attrs["baseline"], impact.attrs["baseline"])
    figure.suptitle(
        f"Atmospheric DA impact on {depth_label}soil-moisture initialization ({date})\n"
        f"Benefit relative to {baseline}; positive values indicate smaller error",
        fontsize=14,
    )
    figure.text(
        0.995,
        0.005,
        impact.attrs.get("caveat", ""),
        ha="right",
        va="bottom",
        fontsize=8,
        color="0.35",
    )
    return figure


def save_impact_products(
    impact: xr.Dataset,
    *,
    data_path: str | Path,
    csv_path: str | Path,
) -> pd.DataFrame:
    """Save the derived fields and scalar summary."""
    data_path = Path(data_path)
    csv_path = Path(csv_path)
    data_path.parent.mkdir(parents=True, exist_ok=True)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    impact.to_netcdf(data_path)
    summary = impact_summary(impact)
    summary.to_csv(csv_path, index=False)
    return summary
