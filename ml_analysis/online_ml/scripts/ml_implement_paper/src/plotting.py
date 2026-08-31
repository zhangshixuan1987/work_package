"""Publication-style plotting for the four online ML diagnostics."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


COLORS = {"CTRL": "#222222", "UNET": "#2878B5", "UNETXTR": "#D95319"}


def set_publication_style() -> None:
    plt.rcParams.update({
        "figure.dpi": 120,
        "savefig.dpi": 300,
        "font.size": 10,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def save_figure(figure: plt.Figure, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, bbox_inches="tight", facecolor="white")
    return path


def plot_rmse_timeseries(dataset: xr.Dataset, variables: Sequence[str], path: str | Path):
    set_publication_style()
    figure, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    for axis, variable in zip(axes.flat, variables):
        for experiment in dataset.experiment.values:
            values = dataset.rmse.sel(variable=variable, experiment=experiment)
            axis.plot(values.time, values, label=str(experiment), color=COLORS.get(str(experiment)))
        axis.set_title(variable)
        axis.set_ylabel(f"RMSE [{dataset.rmse.sel(variable=variable).attrs.get('units', '')}]")
    axes[-1, 0].set_xlabel("Date")
    axes[-1, 1].set_xlabel("Date")
    axes[0, 0].legend(frameon=False)
    figure.suptitle("Daily global atmospheric RMSE")
    return save_figure(figure, path)


def plot_vertical_rmse_ratio(dataset: xr.Dataset, variables: Sequence[str], path: str | Path):
    set_publication_style()
    figure, axes = plt.subplots(2, 2, figsize=(9, 9), sharex=True, sharey=True)
    for axis, variable in zip(axes.flat, variables):
        for experiment in dataset.experiment.values:
            values = dataset.rmse_ratio.sel(variable=variable, experiment=experiment)
            axis.plot(values, values.plev, label=str(experiment), color=COLORS.get(str(experiment)))
        axis.axvline(1.0, color="0.35", linestyle="--", linewidth=1)
        axis.set_title(variable)
        axis.set_xlabel("RMSE / CTRL RMSE")
        axis.set_yscale("log")
        axis.invert_yaxis()
        axis.set_ylabel("Pressure [hPa]")
    axes[0, 0].legend(frameon=False)
    figure.suptitle("Three-month vertical RMSE ratio")
    return save_figure(figure, path)


def plot_spatial_error_reduction(dataset: xr.Dataset, variables: Sequence[str], path: str | Path):
    set_publication_style()
    try:
        import cartopy.crs as ccrs
        projection = ccrs.Robinson()
        transform = ccrs.PlateCarree()
        subplot_kw = {"projection": projection}
    except ImportError:
        transform = None
        subplot_kw = {}
    figure, axes = plt.subplots(2, 4, figsize=(16, 7), subplot_kw=subplot_kw, squeeze=False)
    values = dataset.rmse_reduction
    limit = float(np.nanmax(abs(values.values)))
    image = None
    for row, experiment in enumerate(dataset.experiment.values):
        for column, variable in enumerate(variables):
            axis = axes[row, column]
            field = values.sel(experiment=experiment, variable=variable)
            kwargs = {"ax": axis, "cmap": "RdBu_r", "vmin": -limit, "vmax": limit, "add_colorbar": False}
            if transform is not None and {"lat", "lon"}.issubset(field.coords):
                kwargs["transform"] = transform
                axis.set_global()
                axis.coastlines(linewidth=0.5)
            image = field.plot(**kwargs)
            axis.set_title(f"{experiment}: {variable}")
    figure.colorbar(image, ax=axes, orientation="horizontal", fraction=0.05, pad=0.08, label="CTRL RMSE − ML RMSE")
    figure.suptitle("Spatial error reduction (positive is improvement)")
    return save_figure(figure, path)


def plot_stability(dataset: xr.Dataset, path: str | Path):
    set_publication_style()
    figure, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    state_fields = [("PS_anomaly", "Surface pressure anomaly"), ("TMQ_anomaly", "TMQ anomaly"), ("TOA_imbalance", "FSNT − FLNT")]
    for axis, (name, title) in zip(axes.flat[:3], state_fields):
        for experiment in dataset.experiment.values:
            if name in dataset:
                values = dataset[name].sel(experiment=experiment)
                axis.plot(values.time, values, label=str(experiment), color=COLORS.get(str(experiment)))
        axis.set_title(title)
    if "normalized_tendency_rms" in dataset:
        for experiment in dataset.experiment.values:
            for tendency in dataset.tendency.values:
                values = dataset.normalized_tendency_rms.sel(experiment=experiment, tendency=tendency)
                axes[1, 1].plot(values.time, values, label=f"{experiment} {tendency}", alpha=0.8)
    axes[1, 1].set_title("Normalized ML tendency RMS")
    axes[0, 0].legend(frameon=False)
    axes[1, 1].legend(frameon=False, fontsize=7, ncol=2)
    figure.suptitle("Physical stability and ML forcing")
    return save_figure(figure, path)

