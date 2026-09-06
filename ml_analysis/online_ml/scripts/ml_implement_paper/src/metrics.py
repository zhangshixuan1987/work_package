"""Weighted error and stability metrics for online ML experiments."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import xarray as xr


def _dims_present(array: xr.DataArray, candidates: Sequence[str]) -> list[str]:
    return [dimension for dimension in candidates if dimension in array.dims]


def weighted_mean(
    field: xr.DataArray,
    weights: xr.DataArray,
    dims: Sequence[str] | None = None,
) -> xr.DataArray:
    """NaN-safe weighted mean over the requested horizontal dimensions."""
    dims = list(dims or weights.dims)
    valid_weights = weights.broadcast_like(field).where(field.notnull())
    numerator = (field * valid_weights).sum(dims, skipna=True)
    denominator = valid_weights.sum(dims, skipna=True)
    return numerator / denominator


def weighted_rmse(
    experiment: xr.DataArray,
    reference: xr.DataArray,
    weights: xr.DataArray,
    dims: Sequence[str] | None = None,
) -> xr.DataArray:
    """Area-weighted root-mean-square error."""
    return np.sqrt(weighted_mean((experiment - reference) ** 2, weights, dims))


def pressure_weights(levels: xr.DataArray) -> xr.DataArray:
    """Approximate pressure-layer thickness weights from pressure midpoints."""
    values = np.asarray(levels, dtype=float)
    order = np.argsort(values)
    sorted_values = values[order]
    edges = np.empty(len(values) + 1)
    edges[1:-1] = 0.5 * (sorted_values[:-1] + sorted_values[1:])
    edges[0] = max(0.0, sorted_values[0] - 0.5 * (sorted_values[1] - sorted_values[0]))
    edges[-1] = sorted_values[-1] + 0.5 * (sorted_values[-1] - sorted_values[-2])
    thickness = np.diff(edges)
    restored = np.empty_like(thickness)
    restored[order] = thickness
    return xr.DataArray(restored, dims=levels.dims, coords=levels.coords, name="dp")


def pressure_area_rmse(
    experiment: xr.DataArray,
    reference: xr.DataArray,
    area: xr.DataArray,
    pressure_dim: str = "plev",
) -> xr.DataArray:
    """Area- and pressure-layer-weighted atmospheric RMSE."""
    dp = pressure_weights(experiment[pressure_dim])
    combined = area * dp
    return weighted_rmse(experiment, reference, combined, list(combined.dims))


def vertical_rmse_profile(
    experiment: xr.DataArray,
    reference: xr.DataArray,
    area: xr.DataArray,
) -> xr.DataArray:
    """RMSE at every pressure level over time and horizontal space."""
    squared = (experiment - reference) ** 2
    horizontal = list(area.dims)
    mean_square = weighted_mean(squared, area, horizontal)
    if "time" in mean_square.dims:
        mean_square = mean_square.mean("time", skipna=True)
    return np.sqrt(mean_square)


def spatial_rmse(experiment: xr.DataArray, reference: xr.DataArray) -> xr.DataArray:
    """Time RMSE at each horizontal grid cell."""
    return np.sqrt(((experiment - reference) ** 2).mean("time", skipna=True))


def temporal_pearson_correlation(
    experiment: xr.DataArray,
    reference: xr.DataArray,
) -> xr.DataArray:
    """Pearson correlation over time at each horizontal grid cell."""
    experiment_anomaly = experiment - experiment.mean("time", skipna=True)
    reference_anomaly = reference - reference.mean("time", skipna=True)
    covariance = (experiment_anomaly * reference_anomaly).mean("time", skipna=True)
    experiment_std = np.sqrt((experiment_anomaly**2).mean("time", skipna=True))
    reference_std = np.sqrt((reference_anomaly**2).mean("time", skipna=True))
    denominator = experiment_std * reference_std
    return covariance / denominator.where(denominator > 0)


def rmse_ratio(ml_rmse: xr.DataArray, control_rmse: xr.DataArray) -> xr.DataArray:
    return ml_rmse / control_rmse.where(control_rmse > 0)


def percentage_improvement(
    ml_rmse: xr.DataArray,
    control_rmse: xr.DataArray,
    minimum_control: float = 0.0,
) -> xr.DataArray:
    valid_control = control_rmse.where(control_rmse > minimum_control)
    return 100.0 * (1.0 - ml_rmse / valid_control)


def tendency_rms(
    tendency: xr.DataArray,
    area: xr.DataArray,
    pressure_dim: str = "plev",
) -> tuple[xr.DataArray, xr.DataArray | None]:
    """Return global and pressure-dependent RMS of an ML tendency."""
    profile = np.sqrt(weighted_mean(tendency**2, area, list(area.dims)))
    profile = profile.resample(time="1D").mean() if "time" in profile.dims else profile
    if pressure_dim not in profile.dims:
        return profile, None
    dp = pressure_weights(profile[pressure_dim])
    global_rms = np.sqrt((profile**2 * dp).sum(pressure_dim) / dp.sum())
    return global_rms, profile


def regional_weighted_statistics(
    field: xr.DataArray,
    area: xr.DataArray,
    latitude: xr.DataArray,
    land_fraction: xr.DataArray | None = None,
) -> xr.DataArray:
    """Area means for global, latitude-band, land, and ocean regions."""
    masks: dict[str, xr.DataArray] = {
        "Global": xr.ones_like(latitude, dtype=bool),
        "Tropics": abs(latitude) <= 30,
        "NH_extratropics": latitude > 30,
        "SH_extratropics": latitude < -30,
    }
    if land_fraction is not None:
        masks["Land"] = land_fraction >= 0.5
        masks["Ocean"] = land_fraction < 0.5
    values = [weighted_mean(field.where(mask), area.where(mask), list(area.dims)) for mask in masks.values()]
    return xr.concat(values, dim=xr.IndexVariable("region", list(masks)))


def improvement_timeseries(
    rmse: xr.DataArray,
    control_name: str = "CTRL",
    ml_names: Sequence[str] = ("UNET", "UNETXTR"),
) -> xr.DataArray:
    return xr.concat(
        [percentage_improvement(rmse.sel(experiment=name), rmse.sel(experiment=control_name)) for name in ml_names],
        dim=xr.IndexVariable("experiment", list(ml_names)),
    )


def global_rmse_diagnostic(
    experiments: Mapping[str, xr.Dataset],
    reference: xr.Dataset,
    variables: Sequence[str],
    area: xr.DataArray,
) -> xr.Dataset:
    """Build Figure 1 daily global RMSE and improvement diagnostics."""
    by_variable = []
    for variable in variables:
        by_experiment = [
            pressure_area_rmse(dataset[variable], reference[variable], area)
            for dataset in experiments.values()
        ]
        by_variable.append(
            xr.concat(by_experiment, dim=xr.IndexVariable("experiment", list(experiments)))
        )
    rmse = xr.concat(by_variable, dim=xr.IndexVariable("variable", list(variables)))
    rmse.name = "rmse"
    improvement = improvement_timeseries(rmse)
    improvement.name = "improvement_percent"
    return xr.Dataset({"rmse": rmse, "improvement_percent": improvement})


def vertical_rmse_diagnostic(
    experiments: Mapping[str, xr.Dataset],
    reference: xr.Dataset,
    variables: Sequence[str],
    area: xr.DataArray,
) -> xr.Dataset:
    """Build Figure 2 pressure-level RMSE and ML/control ratios."""
    profiles = []
    for variable in variables:
        values = [
            vertical_rmse_profile(dataset[variable], reference[variable], area)
            for dataset in experiments.values()
        ]
        profiles.append(xr.concat(values, dim=xr.IndexVariable("experiment", list(experiments))))
    rmse = xr.concat(profiles, dim=xr.IndexVariable("variable", list(variables)))
    ratios = xr.concat(
        [rmse_ratio(rmse.sel(experiment=name), rmse.sel(experiment="CTRL")) for name in ("UNET", "UNETXTR")],
        dim=xr.IndexVariable("experiment", ["UNET", "UNETXTR"]),
    )
    return xr.Dataset({"rmse": rmse, "rmse_ratio": ratios})


def spatial_reduction_diagnostic(
    experiments: Mapping[str, xr.Dataset],
    reference: xr.Dataset,
    selections: Mapping[str, tuple[str, float]],
    area: xr.DataArray,
    latitude: xr.DataArray,
    land_fraction: xr.DataArray | None = None,
    minimum_control_rmse: float = 1.0e-12,
) -> xr.Dataset:
    """Build Figure 3 spatial absolute/percentage reductions and regional means."""
    reductions = []
    percentages = []
    summaries = []
    for label, (variable, level) in selections.items():
        reference_field = reference[variable].sel(plev=level, method="nearest")
        control = spatial_rmse(experiments["CTRL"][variable].sel(plev=level, method="nearest"), reference_field)
        variable_reductions = []
        variable_percentages = []
        variable_summaries = []
        for name in ("UNET", "UNETXTR"):
            ml = spatial_rmse(experiments[name][variable].sel(plev=level, method="nearest"), reference_field)
            reduction = control - ml
            percentage = 100.0 * reduction / control.where(control > minimum_control_rmse)
            variable_reductions.append(reduction)
            variable_percentages.append(percentage)
            variable_summaries.append(
                regional_weighted_statistics(reduction, area, latitude, land_fraction)
            )
        experiment_axis = xr.IndexVariable("experiment", ["UNET", "UNETXTR"])
        reductions.append(xr.concat(variable_reductions, dim=experiment_axis))
        percentages.append(xr.concat(variable_percentages, dim=experiment_axis))
        summaries.append(xr.concat(variable_summaries, dim=experiment_axis))
    variable_axis = xr.IndexVariable("variable", list(selections))
    return xr.Dataset({
        "rmse_reduction": xr.concat(reductions, dim=variable_axis),
        "improvement_percent": xr.concat(percentages, dim=variable_axis),
        "regional_mean_reduction": xr.concat(summaries, dim=variable_axis),
    })


def stability_diagnostic(
    experiments: Mapping[str, xr.Dataset],
    area: xr.DataArray,
    tendencies: Sequence[str] = ("Nudge_U", "Nudge_V", "Nudge_T", "Nudge_Q"),
) -> xr.Dataset:
    """Build Figure 4 state-drift and tendency diagnostics."""
    result = xr.Dataset()
    experiment_axis = xr.IndexVariable("experiment", list(experiments))
    for output_name, source_name in (
        ("PS_anomaly", "PS"),
        ("TMQ_anomaly", "TMQ"),
        ("TOA_imbalance", None),
    ):
        series = []
        for dataset in experiments.values():
            field = dataset["FSNT"] - dataset["FLNT"] if source_name is None else dataset[source_name]
            mean = weighted_mean(field, area, list(area.dims))
            if output_name.endswith("anomaly"):
                mean = mean - mean.isel(time=0)
            series.append(mean)
        result[output_name] = xr.concat(series, dim=experiment_axis)

    global_values = []
    profile_values = []
    tendency_experiments = []
    for experiment_name, dataset in experiments.items():
        available = [name for name in tendencies if name in dataset]
        if not available:
            continue
        tendency_experiments.append(experiment_name)
        globals_for_experiment = []
        profiles_for_experiment = []
        for name in tendencies:
            if name not in dataset:
                template = globals_for_experiment[0] if globals_for_experiment else weighted_mean(dataset["PS"], area)
                globals_for_experiment.append(xr.full_like(template, np.nan))
                continue
            global_rms, profile = tendency_rms(dataset[name], area)
            globals_for_experiment.append(global_rms)
            if profile is not None:
                profiles_for_experiment.append(profile.expand_dims(tendency=[name]))
        global_values.append(xr.concat(globals_for_experiment, dim=xr.IndexVariable("tendency", list(tendencies))))
        if profiles_for_experiment:
            profile_values.append(xr.concat(profiles_for_experiment, dim="tendency"))
    if global_values:
        rms = xr.concat(global_values, dim=xr.IndexVariable("experiment", tendency_experiments))
        result["tendency_rms"] = rms
        result["normalized_tendency_rms"] = rms / rms.mean("time", skipna=True)
    if profile_values:
        result["tendency_rms_profile"] = xr.concat(
            profile_values, dim=xr.IndexVariable("experiment", tendency_experiments)
        )
    return result
