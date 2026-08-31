"""Preprocessing utilities for E3SM and reference atmospheric datasets."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import xarray as xr


def subset_time(dataset: xr.Dataset, start: str, end: str) -> xr.Dataset:
    """Select the inclusive analysis period."""
    return dataset.sel(time=slice(start, end))


def match_common_times(datasets: Mapping[str, xr.Dataset]) -> dict[str, xr.Dataset]:
    """Restrict every dataset to exactly the timestamps shared by all datasets."""
    common = None
    for dataset in datasets.values():
        values = dataset.indexes["time"]
        common = values if common is None else common.intersection(values)
    if common is None or len(common) == 0:
        raise ValueError("Datasets have no common timestamps")
    return {name: dataset.sel(time=common) for name, dataset in datasets.items()}


def daily_mean(dataset: xr.Dataset) -> xr.Dataset:
    """Convert subdaily data to daily means."""
    return dataset.resample(time="1D").mean(keep_attrs=True)


def grid_cell_weights(dataset: xr.Dataset, area_name: str = "area") -> xr.DataArray:
    """Return normalized E3SM cell area, falling back to cosine latitude."""
    if area_name in dataset:
        weights = dataset[area_name].astype("float64")
    elif "lat" in dataset.coords:
        weights = np.cos(np.deg2rad(dataset["lat"])).clip(min=0.0)
    else:
        raise KeyError(f"Neither {area_name!r} nor latitude is available for weighting")
    weights = weights.where(np.isfinite(weights) & (weights > 0))
    return weights / weights.sum(skipna=True)


def _pressure_to_hpa(pressure: xr.DataArray) -> xr.DataArray:
    units = pressure.attrs.get("units", "").lower()
    if units in {"pa", "pascal", "pascals"} or float(pressure.max()) > 2_000:
        pressure = pressure / 100.0
    pressure.attrs["units"] = "hPa"
    return pressure


def hybrid_pressure(
    dataset: xr.Dataset,
    *,
    ps_name: str = "PS",
    hyam_name: str = "hyam",
    hybm_name: str = "hybm",
    p0_name: str = "P0",
) -> xr.DataArray:
    """Calculate E3SM hybrid-level mid-point pressure in hPa."""
    p0 = dataset[p0_name] if p0_name in dataset else 100_000.0
    pressure = dataset[hyam_name] * p0 + dataset[hybm_name] * dataset[ps_name]
    pressure.name = "pressure"
    pressure.attrs["units"] = "Pa"
    return _pressure_to_hpa(pressure)


def _interp_profile(values: np.ndarray, pressure: np.ndarray, targets: np.ndarray) -> np.ndarray:
    valid = np.isfinite(values) & np.isfinite(pressure)
    if valid.sum() < 2:
        return np.full(targets.shape, np.nan)
    xp = pressure[valid]
    fp = values[valid]
    order = np.argsort(xp)
    return np.interp(targets, xp[order], fp[order], left=np.nan, right=np.nan)


def interpolate_to_pressure(
    field: xr.DataArray,
    target_levels: Sequence[float],
    *,
    pressure: xr.DataArray | None = None,
    level_dim: str = "lev",
) -> xr.DataArray:
    """Interpolate a field to common pressure levels in hPa."""
    targets = xr.DataArray(np.asarray(target_levels, dtype=float), dims="plev", name="plev")
    targets.attrs["units"] = "hPa"
    if pressure is None:
        if level_dim not in field.coords:
            raise KeyError(f"{field.name} has no {level_dim!r} pressure coordinate")
        coordinate = _pressure_to_hpa(field[level_dim])
        result = field.assign_coords({level_dim: coordinate}).interp({level_dim: targets})
        return result.rename({level_dim: "plev"}) if level_dim in result.dims else result
    return xr.apply_ufunc(
        _interp_profile,
        field,
        _pressure_to_hpa(pressure),
        targets,
        input_core_dims=[[level_dim], [level_dim], ["plev"]],
        output_core_dims=[["plev"]],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[float],
        dask_gufunc_kwargs={"output_sizes": {"plev": len(target_levels)}},
    ).assign_coords(plev=targets)


def extract_level(field: xr.DataArray, pressure_hpa: float) -> xr.DataArray:
    """Extract the nearest available common pressure level."""
    return field.sel(plev=pressure_hpa, method="nearest")


def interpolate_dataset_fields(
    dataset: xr.Dataset,
    variables: Sequence[str],
    target_levels: Sequence[float],
    *,
    level_dim: str = "lev",
    pressure_name: str | None = None,
    hybrid_names: Mapping[str, str] | None = None,
) -> xr.Dataset:
    """Interpolate configured 3-D fields while retaining all other variables."""
    pressure = dataset[pressure_name] if pressure_name else None
    if pressure is None and hybrid_names and all(
        hybrid_names[key] in dataset for key in ("ps", "hyam", "hybm")
    ):
        pressure = hybrid_pressure(
            dataset,
            ps_name=hybrid_names["ps"],
            hyam_name=hybrid_names["hyam"],
            hybm_name=hybrid_names["hybm"],
            p0_name=hybrid_names.get("p0", "P0"),
        )
    result = dataset.drop_vars([name for name in variables if name in dataset])
    for name in variables:
        if name not in dataset:
            raise KeyError(f"Required atmospheric variable {name!r} is missing")
        result[name] = interpolate_to_pressure(
            dataset[name], target_levels, pressure=pressure, level_dim=level_dim
        )
    return result


def basic_qc(dataset: xr.Dataset, variables: Sequence[str]) -> dict[str, int]:
    """Count NaN and infinite values and reject all-missing variables."""
    report: dict[str, int] = {}
    for name in variables:
        if name not in dataset:
            raise KeyError(f"Required variable {name!r} is missing")
        values = dataset[name]
        finite = np.isfinite(values)
        if not bool(finite.any().compute() if hasattr(finite.data, "compute") else finite.any()):
            raise ValueError(f"Variable {name!r} contains no finite values")
        invalid = (~finite).sum()
        report[name] = int(invalid.compute() if hasattr(invalid.data, "compute") else invalid)
    return report
