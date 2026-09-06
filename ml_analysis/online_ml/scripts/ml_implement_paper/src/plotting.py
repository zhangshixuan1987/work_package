"""Publication-style plotting for the four online ML diagnostics."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence
import warnings

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import numpy as np
import xarray as xr


COLORS = {
    "CTRL": "#222222",
    "UNET-IMT": "#2878B5",
    "UNETXTR-IMT": "#D95319",
    "UNETXTR-IMT-A15": "#9467BD",
    "UNETXTR-IMT-C05": "#2CA02C",
    "UNETXTR-IMT-C05A15": "#8C564B",
    "UNETXTR-LCZ-C05A15": "#17BECF",
    "UNETXTR-LCZ-C05A15-PTAP100": "#E377C2",
    "REF": "#7F7F7F",
}


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


def _matplotlib_time(coordinate: xr.DataArray) -> np.ndarray:
    """Convert a cftime coordinate to values Matplotlib can plot."""
    index = coordinate.to_index()
    if hasattr(index, "to_datetimeindex"):
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message="Converting a CFTimeIndex", category=RuntimeWarning
            )
            try:
                index = index.to_datetimeindex(time_unit="s")
            except TypeError:  # Compatibility with older xarray versions.
                index = index.to_datetimeindex()
    return np.asarray(index)


def save_figure(figure: plt.Figure, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, bbox_inches="tight", facecolor="white")
    return path


def plot_rmse_timeseries(
    dataset: xr.Dataset,
    variables: Sequence[str],
    path: str | Path,
    *,
    figsize: tuple[float, float] = (11, 7),
    layout: tuple[int, int] = (2, 2),
    sharex: bool = True,
    colors: Mapping[str, str] | None = None,
    title: str = "Daily global atmospheric RMSE",
    xlabel: str = "Date",
    legend_loc: str = "best",
    legend_frame: bool = False,
    linewidth: float = 1.5,
    font_size: float = 10,
):
    set_publication_style()
    plt.rcParams.update({"font.size": font_size})
    nrows, ncols = layout
    if nrows * ncols < len(variables):
        raise ValueError(f"layout {layout} has fewer panels than {len(variables)} variables")
    plot_colors = COLORS if colors is None else colors
    figure, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=sharex, squeeze=False)
    for axis, variable in zip(axes.flat, variables):
        for experiment in dataset.experiment.values:
            values = dataset.rmse.sel(variable=variable, experiment=experiment)
            name = str(experiment)
            axis.plot(
                _matplotlib_time(values.time),
                values,
                label=name,
                color=plot_colors.get(name),
                linewidth=linewidth,
            )
        axis.set_title(variable)
        axis.set_ylabel(f"RMSE [{dataset.rmse.sel(variable=variable).attrs.get('units', '')}]")
    for axis in axes.flat[len(variables):]:
        axis.set_visible(False)
    for axis in axes[-1, :]:
        if axis.get_visible():
            axis.set_xlabel(xlabel)
    axes[0, 0].legend(loc=legend_loc, frameon=legend_frame)
    figure.suptitle(title)
    return save_figure(figure, path)


def plot_vertical_rmse_ratio(
    dataset: xr.Dataset,
    variables: Sequence[str],
    path: str | Path,
    *,
    figsize: tuple[float, float] = (9, 9),
    layout: tuple[int, int] = (2, 2),
    colors: Mapping[str, str] | None = None,
    title: str = "Three-month vertical RMSE ratio",
    xlabel: str = "RMSE / CTRL RMSE",
    ylabel: str = "Pressure [hPa]",
    reference_value: float = 1.0,
    pressure_scale: str = "log",
    legend_loc: str = "best",
    legend_frame: bool = False,
    linewidth: float = 1.5,
    font_size: float = 10,
):
    set_publication_style()
    plt.rcParams.update({"font.size": font_size})
    nrows, ncols = layout
    if nrows * ncols < len(variables):
        raise ValueError(f"layout {layout} has fewer panels than {len(variables)} variables")
    plot_colors = COLORS if colors is None else colors
    figure, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True, squeeze=False)
    for axis, variable in zip(axes.flat, variables):
        for experiment in dataset.experiment.values:
            values = dataset.rmse_ratio.sel(variable=variable, experiment=experiment)
            name = str(experiment)
            axis.plot(values, values.plev, label=name, color=plot_colors.get(name), linewidth=linewidth)
        axis.axvline(reference_value, color="0.35", linestyle="--", linewidth=1)
        axis.set_title(variable)
        axis.set_xlabel(xlabel)
        axis.set_yscale(pressure_scale)
        axis.invert_yaxis()
        axis.set_ylabel(ylabel)
    for axis in axes.flat[len(variables):]:
        axis.set_visible(False)
    axes[0, 0].legend(loc=legend_loc, frameon=legend_frame)
    figure.suptitle(title)
    return save_figure(figure, path)


def plot_spatial_error_reduction(
    dataset: xr.Dataset,
    variables: Sequence[str],
    path: str | Path,
    *,
    figsize: tuple[float, float] = (16, 7),
    cmap: str = "RdBu_r",
    color_limit: float | None = None,
    color_levels: Sequence[float] | None = None,
    data_variable: str = "rmse_reduction",
    font_size: float = 10,
    map_projection: str = "robinson",
    longitude_ticks: Sequence[float] = (-180, -120, -60, 0, 60, 120, 180),
    latitude_ticks: Sequence[float] = (-90, -60, -30, 0, 30, 60, 90),
    colorbar_label: str = "CTRL RMSE − ML RMSE",
    title: str = "Spatial error reduction (positive is improvement)",
    coastline_linewidth: float = 0.5,
    annotate_statistics: bool = True,
):
    set_publication_style()
    plt.rcParams.update({"font.size": font_size})
    try:
        import cartopy.crs as ccrs
        if map_projection.lower() in {"platecarree", "latlon", "lat-lon"}:
            projection = ccrs.PlateCarree()
        elif map_projection.lower() == "robinson":
            projection = ccrs.Robinson()
        else:
            raise ValueError(f"Unsupported map_projection: {map_projection!r}")
        transform = ccrs.PlateCarree()
        subplot_kw = {"projection": projection}
    except ImportError:
        transform = None
        subplot_kw = {}
    figure, axes = plt.subplots(len(dataset.experiment), len(variables), figsize=figsize, subplot_kw=subplot_kw, squeeze=False)
    if data_variable not in dataset:
        raise KeyError(f"Dataset does not contain {data_variable!r}")
    values = dataset[data_variable]
    levels = None
    norm = None
    plot_cmap = cmap
    if color_levels is not None:
        levels = np.asarray(color_levels, dtype=float)
        if levels.ndim != 1 or len(levels) < 2 or np.any(np.diff(levels) <= 0):
            raise ValueError("color_levels must be a strictly increasing one-dimensional sequence")
        plot_cmap = plt.colormaps.get_cmap(cmap).resampled(len(levels) - 1)
        norm = BoundaryNorm(levels, plot_cmap.N)
    else:
        limit = float(np.nanmax(abs(values.values))) if color_limit is None else color_limit
    image = None
    for row, experiment in enumerate(dataset.experiment.values):
        for column, variable in enumerate(variables):
            axis = axes[row, column]
            field = values.sel(experiment=experiment, variable=variable)
            kwargs = {"ax": axis, "cmap": plot_cmap, "add_colorbar": False}
            if norm is None:
                kwargs.update({"vmin": -limit, "vmax": limit})
            else:
                kwargs["norm"] = norm
            if transform is not None and {"lat", "lon"}.issubset(field.coords):
                kwargs["transform"] = transform
                axis.set_global()
                axis.coastlines(linewidth=coastline_linewidth)
                if map_projection.lower() in {"platecarree", "latlon", "lat-lon"}:
                    gridlines = axis.gridlines(
                        crs=transform,
                        draw_labels=True,
                        xlocs=list(longitude_ticks),
                        ylocs=list(latitude_ticks),
                        linewidth=0.4,
                        color="0.4",
                        alpha=0.5,
                        linestyle=":",
                    )
                    gridlines.top_labels = False
                    gridlines.right_labels = False
                    gridlines.bottom_labels = row == len(dataset.experiment) - 1
                    gridlines.left_labels = column == 0
                    gridlines.xlabel_style = {"size": font_size}
                    gridlines.ylabel_style = {"size": font_size}
            image = field.plot(**kwargs)
            axis.set_title(f"{experiment}: {variable}")
            if annotate_statistics:
                statistic_names = {
                    "global_rmse": "RMSE",
                    "global_pcc": "PCC",
                    "improved_grid_fraction": "Improve",
                    "degraded_grid_fraction": "Degrade",
                }
                selected = {}
                for name, label in statistic_names.items():
                    if name in dataset:
                        value = float(
                            dataset[name].sel(experiment=experiment, variable=variable).squeeze()
                        )
                        if np.isfinite(value):
                            selected[label] = value
                if selected:
                    first_line = "  ".join(
                        f"{label}={selected[label]:.3g}" for label in ("RMSE", "PCC")
                        if label in selected
                    )
                    second_line = "  ".join(
                        f"{label}={selected[label]:.1f}%" for label in ("Improve", "Degrade")
                        if label in selected
                    )
                    annotation = "\n".join(line for line in (first_line, second_line) if line)
                    axis.text(
                        0.01, 0.02, annotation, transform=axis.transAxes,
                        ha="left", va="bottom", fontsize=0.72 * font_size,
                        bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none", "pad": 1.5},
                        zorder=10,
                    )
    colorbar_kwargs = {
        "ax": axes, "orientation": "horizontal", "fraction": 0.05,
        "pad": 0.08, "label": colorbar_label,
    }
    if levels is not None:
        colorbar_kwargs.update({"boundaries": levels, "ticks": levels, "extend": "both"})
    figure.colorbar(image, **colorbar_kwargs)
    figure.suptitle(title)
    return save_figure(figure, path)


def plot_stability(
    dataset: xr.Dataset,
    path: str | Path,
    *,
    figsize: tuple[float, float] = (11, 7),
    colors: Mapping[str, str] | None = None,
    title: str = "Physical stability and ML forcing",
    state_fields: Sequence[tuple[str, str]] | None = None,
    tendency_title: str = "Normalized ML tendency RMS",
    legend_frame: bool = False,
    tendency_legend_fontsize: float = 7,
    tendency_legend_columns: int = 2,
    linewidth: float = 1.5,
    font_size: float = 10,
):
    set_publication_style()
    plt.rcParams.update({"font.size": font_size})
    plot_colors = COLORS if colors is None else colors
    figure, axes = plt.subplots(2, 2, figsize=figsize, sharex=True)
    if state_fields is None:
        state_fields = [("PS_anomaly", "Surface pressure anomaly"), ("TMQ_anomaly", "TMQ anomaly"), ("TOA_imbalance", "FSNT − FLNT")]
    for axis, (field_name, panel_title) in zip(axes.flat[:3], state_fields):
        for experiment in dataset.experiment.values:
            if field_name in dataset:
                values = dataset[field_name].sel(experiment=experiment)
                experiment_name = str(experiment)
                axis.plot(
                    _matplotlib_time(values.time),
                    values,
                    label=experiment_name,
                    color=plot_colors.get(experiment_name),
                    linewidth=linewidth,
                )
        axis.set_title(panel_title)
    if "normalized_tendency_rms" in dataset:
        for experiment in dataset.experiment.values:
            for tendency in dataset.tendency.values:
                values = dataset.normalized_tendency_rms.sel(experiment=experiment, tendency=tendency)
                axes[1, 1].plot(
                    _matplotlib_time(values.time),
                    values,
                    label=f"{experiment} {tendency}",
                    alpha=0.8,
                )
        axes[1, 1].set_title(tendency_title)
        axes[1, 1].legend(frameon=legend_frame, fontsize=tendency_legend_fontsize, ncol=tendency_legend_columns)
    else:
        axes[1, 1].set_visible(False)
    axes[0, 0].legend(frameon=legend_frame)
    figure.suptitle(title)
    return save_figure(figure, path)


def plot_monthly_ml_forcing(
    dataset: xr.Dataset,
    path: str | Path,
    *,
    figsize: tuple[float, float] = (11, 7),
    colors: Mapping[str, str] | None = None,
    title: str = "Monthly normalized ML-forcing RMS",
    legend_frame: bool = False,
    linewidth: float = 1.5,
    font_size: float = 10,
):
    """Plot one monthly normalized forcing-RMS panel per tendency."""
    set_publication_style()
    plt.rcParams.update({"font.size": font_size})
    plot_colors = COLORS if colors is None else colors
    tendencies = list(dataset.tendency.values)
    figure, axes = plt.subplots(2, 2, figsize=figsize, sharex=True, squeeze=False)
    for axis, tendency in zip(axes.flat, tendencies):
        for experiment in dataset.experiment.values:
            values = dataset.normalized_tendency_rms.sel(
                experiment=experiment, tendency=tendency
            )
            experiment_name = str(experiment)
            axis.plot(
                _matplotlib_time(values.time),
                values,
                label=experiment_name,
                color=plot_colors.get(experiment_name),
                linewidth=linewidth,
            )
        axis.axhline(1.0, color="0.35", linestyle="--", linewidth=1)
        axis.set_title(str(tendency))
        axis.set_ylabel("RMS / time-mean RMS")
    for axis in axes.flat[len(tendencies):]:
        axis.set_visible(False)
    axes[0, 0].legend(frameon=legend_frame)
    figure.suptitle(title)
    return save_figure(figure, path)


def plot_monthly_vertical_bias(
    dataset: xr.Dataset,
    path: str | Path,
    *,
    variables: Sequence[str] = ('U', 'V', 'T', 'Q'),
    month: int | None = None,
    figsize: tuple[float, float] | None = None,
    color_limits: Mapping[str, float] | None = None,
    cmap: str = 'RdBu_r',
    pressure_scale: str = 'log',
    font_size: float = 12,
):
    """Latitude-pressure biases with row labels and dedicated column colorbars.

    month is YYYYMM; None gives a day-weighted mean over cached months.
    Each variable uses one symmetric color scale shared across experiments.
    """
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
    from .monthly_bias import cross_section_statistics

    set_publication_style()
    plt.rcParams.update({'font.size': font_size})
    required_fields = [field for variable in variables
                       for field in (variable, f'{variable}_model', f'{variable}_reference')]
    missing = set(required_fields) - set(dataset.data_vars)
    if missing:
        raise ValueError('Recompute Figure 6 caches to add paired fields for RMSE/PCC: '
                         + ', '.join(sorted(missing)))
    if month is None:
        fields = dataset[required_fields].weighted(dataset.days_in_month).mean('month')
        label = f'{int(dataset.month.min())}–{int(dataset.month.max())} day-weighted mean'
    else:
        fields = dataset[required_fields].sel(month=month)
        label = str(month)
    nrows, ncols = dataset.sizes['experiment'], len(variables)
    figure = plt.figure(
        figsize=figsize or (2.4 + 4 * ncols, 2.2 * nrows + 1.6),
        constrained_layout=True,
    )
    # A separate label column avoids repeating long case names in every panel.
    # Fixed-height colorbar axes prevent the large gap created by spanning rows.
    grid = figure.add_gridspec(
        nrows + 1, ncols + 1,
        width_ratios=[0.58] + [1] * ncols,
        height_ratios=[1] * nrows + [0.09],
        hspace=0.04, wspace=0.04,
    )
    axes = np.empty((nrows, ncols), dtype=object)
    for row, experiment in enumerate(dataset.experiment.values):
        label_axis = figure.add_subplot(grid[row, 0])
        label_axis.set_axis_off()
        case_label = str(experiment).replace('-IMT', '\nIMT').replace('-LCZ', '\nLCZ')
        case_label = case_label.replace('-PTAP100', '\nPTAP100')
        label_axis.text(0.98, 0.5, case_label, ha='right', va='center',
                        fontsize=font_size, fontweight='semibold', linespacing=1.4)
        for col in range(ncols):
            shared = axes[0, 0] if row or col else None
            axes[row, col] = figure.add_subplot(grid[row, col + 1], sharex=shared, sharey=shared)

    for col, variable in enumerate(variables):
        values = fields[variable].load()
        units = dataset[variable].attrs.get('units', '')
        display_factor = 1.0
        if variable == 'Q' and units.lower().replace(' ', '') in {'kg/kg', 'kgkg-1', 'kgkg^-1', 'kgkg**-1'}:
            display_factor = 1000.0
            values = values * display_factor
            units = 'g kg⁻¹'
        finite = np.asarray(values)[np.isfinite(values)]
        if not finite.size:
            raise ValueError(f'{variable}: no finite bias values to plot')
        limit = (color_limits or {}).get(variable, max(float(np.max(np.abs(finite))), 1e-12))
        if not np.isfinite(limit) or limit <= 0:
            raise ValueError(f'{variable}: color limit must be positive and finite')
        for row, experiment in enumerate(dataset.experiment.values):
            axis = axes[row, col]
            field = values.sel(experiment=experiment).transpose('plev', 'lat').sortby('plev').sortby('lat')
            panel = axis.contourf(field.lat, field.plev, field,
                                  levels=np.linspace(-limit, limit, 21), cmap=cmap, extend='both')
            model_section = fields[f'{variable}_model'].sel(experiment=experiment) * display_factor
            reference_section = fields[f'{variable}_reference'].sel(experiment=experiment) * display_factor
            rmse, pcc = cross_section_statistics(model_section, reference_section)
            pcc_label = f'{pcc:.2f}' if np.isfinite(pcc) else 'N/A'
            axis.text(0.02, 0.025, f'RMSE={rmse:.3g} {units}  PCC={pcc_label}',
                      transform=axis.transAxes, ha='left', va='bottom',
                      fontsize=0.78 * font_size,
                      bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.85, 'pad': 1.5},
                      zorder=10)
            axis.set_yscale(pressure_scale)
            bottom, top = float(field.plev.max()), float(field.plev.min())
            axis.set_ylim(bottom, top)
            tick_candidates = (1, 3, 10, 30, 100, 300, 1000) if top < 100 else (100, 200, 300, 500, 700, 1000)
            pressure_ticks = [p for p in tick_candidates if top <= p <= bottom]
            axis.yaxis.set_major_locator(FixedLocator(pressure_ticks))
            axis.yaxis.set_major_formatter(FuncFormatter(lambda value, pos: f'{value:g}'))
            axis.yaxis.set_minor_locator(NullLocator())
            axis.xaxis.set_major_locator(FixedLocator([-90, -60, -30, 0, 30, 60, 90]))
            axis.xaxis.set_major_formatter(FuncFormatter(
                lambda value, pos: '0°' if value == 0 else f'{abs(value):g}°{"S" if value < 0 else "N"}'
            ))
            axis.tick_params(labelleft=col == 0, labelbottom=row == nrows - 1,
                             labelsize=font_size - 1, length=3)
            axis.grid(True, color='0.4', linewidth=0.4, alpha=0.2)
            if row == 0:
                axis.set_title(variable, fontsize=font_size + 2, fontweight='semibold', pad=9)
            if col == 0:
                axis.set_ylabel('Pressure [hPa]', fontsize=font_size - 1)
            if row == nrows - 1:
                axis.set_xlabel('Latitude', fontsize=font_size)
        colorbar_axis = figure.add_subplot(grid[-1, col + 1])
        colorbar = figure.colorbar(
            panel, cax=colorbar_axis, orientation='horizontal',
            ticks=np.linspace(-limit, limit, 5), format='%.2g',
        )
        colorbar.set_label(f'{variable} bias [{units}]', fontsize=font_size)
        colorbar.ax.tick_params(labelsize=font_size - 1, length=3)
    reference_label = dataset.attrs.get('reference_label', 'REF')
    figure.suptitle(f'Monthly zonal-mean bias (experiment − {reference_label})\n{label}',
                   fontsize=font_size + 2)
    return save_figure(figure, path)
