"""2D spatial mapping diagnostics and publication-quality plotting.

Provides unified mapping workflows for 2D fields ('T500', 'U200', 'U850', 'TMQ', 'PRECT'),
spatial error reductions (Figure 3 style), and spatial climatology / bias maps
following the layout and visual aesthetics of Figure 6 (closed axes, left experiment
label column, top variable headers, dedicated per-column colorbars, and metric annotations).
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import numpy as np
import xarray as xr

from .plotting import save_figure, set_publication_style
from . import spatial_cache

DEFAULT_VARIABLES: tuple[str, ...] = ('T500', 'U200', 'U850', 'TMQ', 'PRECT')

VARIABLE_UNITS: dict[str, str] = {
    'T500': 'K',
    'U200': 'm s⁻¹',
    'U850': 'm s⁻¹',
    'TMQ': 'kg m⁻²',
    'PRECT': 'mm day⁻¹',
}

VARIABLE_LONG_NAMES: dict[str, str] = {
    'T500': 'Temperature at 500 hPa',
    'U200': 'Zonal wind at 200 hPa',
    'U850': 'Zonal wind at 850 hPa',
    'TMQ': 'Total Precipitable Water',
    'PRECT': 'Total Precipitation Rate',
}

VARIABLE_SCALES: dict[str, float] = {
    'PRECT': 86400.0 * 1000.0,  # converts m s⁻¹ -> mm day⁻¹
}

DEFAULT_BIAS_LIMITS: dict[str, float] = {
    'T500': 3.0,
    'U200': 5.0,
    'U850': 3.0,
    'TMQ': 5.0,
    'PRECT': 5.0,
}

DEFAULT_COLOR_LEVELS: tuple[float, ...] = (-100, -50, -25, -10, -5, 0, 5, 10, 25, 50, 100)

DEFAULT_SEASONS: tuple[str, ...] = ('ANN', 'DJF', 'MAM', 'JJA', 'SON')


def format_case_label(experiment: str) -> str:
    """Format experiment names with newlines for clean multi-line display on label axes."""
    label = str(experiment)
    for token in ('-IMT', '-LCZ', '-PTAP100'):
        label = label.replace(token, f'\n{token[1:]}')
    return label


def spatial_field_statistics(
    model: xr.DataArray | np.ndarray,
    reference: xr.DataArray | np.ndarray,
    lat: xr.DataArray | np.ndarray | None = None,
    area: xr.DataArray | np.ndarray | None = None,
) -> tuple[float, float]:
    """Calculate area-weighted RMSE and Pearson pattern correlation (PCC) of two 2D fields.

    Parameters
    ----------
    model : xr.DataArray or np.ndarray
        Model 2D horizontal field (lat, lon).
    reference : xr.DataArray or np.ndarray
        Reference 2D horizontal field (lat, lon).
    lat : array-like, optional
        Latitude coordinates used to compute cos(lat) area weights if area is None.
    area : array-like, optional
        Precomputed grid cell area weights.

    Returns
    -------
    tuple[float, float]
        (rmse, pcc)
    """
    if isinstance(model, xr.DataArray) and isinstance(reference, xr.DataArray):
        model, reference = xr.align(model, reference, join='exact')
        if lat is None and 'lat' in model.coords:
            lat = model.lat.values
    x = np.asarray(model, dtype=float)
    y = np.asarray(reference, dtype=float)

    if area is not None:
        weights = np.asarray(area, dtype=float)
    elif lat is not None:
        lat_rad = np.deg2rad(np.asarray(lat, dtype=float))
        weights = np.cos(lat_rad)
        if weights.ndim == 1 and x.ndim == 2:
            weights = weights[:, None] * np.ones((1, x.shape[1]))
    else:
        weights = np.ones_like(x)

    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(weights) & (weights > 0)
    if not valid.any():
        return float('nan'), float('nan')

    x_v, y_v, w_v = x[valid], y[valid], weights[valid]
    w_norm = w_v / np.sum(w_v)

    rmse = float(np.sqrt(np.sum(w_norm * (x_v - y_v) ** 2)))

    x_mean = np.sum(w_norm * x_v)
    y_mean = np.sum(w_norm * y_v)
    dx = x_v - x_mean
    dy = y_v - y_mean
    denominator = np.sqrt(np.sum(w_norm * dx ** 2) * np.sum(w_norm * dy ** 2))
    if x_v.size > 1 and np.ptp(x_v) > 0 and np.ptp(y_v) > 0 and denominator > 0:
        pcc = float(np.clip(np.sum(w_norm * dx * dy) / denominator, -1.0, 1.0))
    else:
        pcc = float('nan')

    return rmse, pcc


def extract_seasonal_means(
    da: xr.DataArray,
    seasons: Sequence[str] = DEFAULT_SEASONS,
) -> xr.DataArray:
    """Extract annual and meteorological seasonal means from a 3-hourly DataArray."""
    if 'time' not in da.dims:
        if 'season' in da.dims:
            return da
        return da.expand_dims(season=['ANN'])

    parts = []
    season_names = []
    if 'ANN' in seasons:
        parts.append(da.mean('time').compute().expand_dims(season=['ANN']))
        season_names.append('ANN')

    met_seasons = [s for s in seasons if s in {'DJF', 'MAM', 'JJA', 'SON'}]
    if met_seasons:
        grouped = da.groupby('time.season').mean('time').compute()
        avail = [s for s in met_seasons if s in grouped.season.values]
        if avail:
            parts.append(grouped.sel(season=avail))
            season_names.extend(avail)

    combined = xr.concat(parts, dim='season')
    combined.attrs.update(da.attrs)
    return combined


def compute_spatial_bias(
    experiment: xr.Dataset,
    reference: xr.Dataset,
    variables: Sequence[str] = DEFAULT_VARIABLES,
    *,
    seasons: Sequence[str] = DEFAULT_SEASONS,
    scales: Mapping[str, float] | None = None,
    units: Mapping[str, str] | None = None,
) -> xr.Dataset:
    """Compute 2D spatial climatology / mean bias across seasons (experiment minus reference).

    Preserves model mean, reference mean, and bias field across requested seasons,
    applying standard display scaling (e.g. m/s to mm/day for PRECT).
    """
    scale_map = dict(VARIABLE_SCALES)
    if scales is not None:
        scale_map.update(scales)
    unit_map = dict(VARIABLE_UNITS)
    if units is not None:
        unit_map.update(units)

    result = xr.Dataset()
    for var in variables:
        if var not in experiment or var not in reference:
            continue
        exp_da = experiment[var]
        ref_da = reference[var]

        # Extract annual and seasonal means if time dimension is present
        if 'time' in exp_da.dims:
            exp_mean = extract_seasonal_means(exp_da, seasons)
        elif 'season' in exp_da.dims:
            exp_mean = exp_da
        else:
            exp_mean = exp_da.expand_dims(season=['ANN'])

        if 'time' in ref_da.dims:
            ref_mean = extract_seasonal_means(ref_da, seasons)
        elif 'season' in ref_da.dims:
            ref_mean = ref_da
        else:
            ref_mean = ref_da.expand_dims(season=['ANN'])

        exp_mean, ref_mean = xr.align(exp_mean, ref_mean, join='exact')

        factor = scale_map.get(var, 1.0)
        target_unit = unit_map.get(var, exp_da.attrs.get('units', ''))

        exp_scaled = exp_mean * factor
        ref_scaled = ref_mean * factor
        bias_field = exp_scaled - ref_scaled

        exp_scaled.attrs.update(units=target_unit, long_name=f'{var} model mean')
        ref_scaled.attrs.update(units=target_unit, long_name=f'{var} reference mean')
        bias_field.attrs.update(units=target_unit, long_name=f'{var} mean bias (model − ref)')

        result[f'{var}_model'] = exp_scaled
        result[f'{var}_reference'] = ref_scaled
        result[var] = bias_field

        # Record metrics for annotation per season
        if 'season' in bias_field.dims:
            rmse_vals, pcc_vals = [], []
            for s in bias_field.season.values:
                r_val, p_val = spatial_field_statistics(
                    exp_scaled.sel(season=s), ref_scaled.sel(season=s), lat=bias_field.lat
                )
                rmse_vals.append(r_val)
                pcc_vals.append(p_val)
            result[f'{var}_rmse'] = xr.DataArray(rmse_vals, coords={'season': bias_field.season}, dims=['season'])
            result[f'{var}_pcc'] = xr.DataArray(pcc_vals, coords={'season': bias_field.season}, dims=['season'])
        else:
            rmse, pcc = spatial_field_statistics(exp_scaled, ref_scaled, lat=bias_field.lat)
            result[f'{var}_rmse'] = xr.DataArray(rmse)
            result[f'{var}_pcc'] = xr.DataArray(pcc)

    result.attrs['bias_definition'] = 'experiment minus reference 2D time mean'
    return result


def spatial_bias_signature(
    case: str,
    case_dirs: Mapping[str, Path | str],
    post_subdir: Path | str,
    period: str,
    variables: Sequence[str] = DEFAULT_VARIABLES,
    ref_case: str = 'REF',
    seasons: Sequence[str] = DEFAULT_SEASONS,
) -> str:
    """Unique signature for 2D spatial bias checkpoint."""
    import json
    return json.dumps({
        'version': 2,
        'case': case,
        'period': period,
        'variables': list(variables),
        'seasons': list(seasons),
        'model_sources': {v: str(Path(case_dirs[case]) / post_subdir / f'{v}_{period}.nc') for v in variables},
        'ref_sources': {v: str(Path(case_dirs[ref_case]) / post_subdir / f'{v}_{period}.nc') for v in variables},
    }, sort_keys=True)


def spatial_bias_cache_valid(
    path: str | Path,
    variables: Sequence[str] = DEFAULT_VARIABLES,
    expected_signature: str | None = None,
) -> bool:
    """Validate Figure 7 NetCDF checkpoint file."""
    path = Path(path)
    if not path.exists():
        return False
    try:
        with xr.open_dataset(path) as ds:
            for v in variables:
                if v not in ds.data_vars:
                    return False
            if 'season' not in ds.dims:
                return False
            if expected_signature and ds.attrs.get('cache_signature') != expected_signature:
                return False
            return True
    except Exception:
        return False


def compute_case_spatial_bias(
    case: str,
    case_dirs: Mapping[str, Path | str],
    post_subdir: Path | str,
    period: str,
    variables: Sequence[str] = DEFAULT_VARIABLES,
    *,
    ref_case: str = 'REF',
    ref_dataset: xr.Dataset | None = None,
    seasons: Sequence[str] = DEFAULT_SEASONS,
    scales: Mapping[str, float] | None = None,
    units: Mapping[str, str] | None = None,
) -> xr.Dataset:
    """Compute annual and seasonal 2D bias for one experiment vs reference."""
    scale_map = dict(VARIABLE_SCALES)
    if scales:
        scale_map.update(scales)
    unit_map = dict(VARIABLE_UNITS)
    if units:
        unit_map.update(units)

    case_dir = Path(case_dirs[case])
    post_sub = Path(post_subdir)

    model_ds = xr.Dataset()
    for v in variables:
        p = case_dir / post_sub / f'{v}_{period}.nc'
        with xr.open_dataset(p) as ds:
            model_ds[v] = extract_seasonal_means(ds[v], seasons)

    if ref_dataset is None:
        ref_dir = Path(case_dirs[ref_case])
        ref_ds = xr.Dataset()
        for v in variables:
            p = ref_dir / post_sub / f'{v}_{period}.nc'
            with xr.open_dataset(p) as ds:
                ref_ds[v] = extract_seasonal_means(ds[v], seasons)
    else:
        ref_ds = ref_dataset

    bias_ds = compute_spatial_bias(model_ds, ref_ds, variables, seasons=seasons, scales=scale_map, units=unit_map)
    bias_ds = bias_ds.expand_dims(experiment=[case])
    bias_ds.attrs['experiment'] = case
    bias_ds.attrs['period'] = period
    bias_ds.attrs['reference_case'] = str(case_dirs[ref_case])
    return bias_ds


def load_spatial_products(
    pair_paths: Mapping[tuple[str, str], str | Path],
    cases: Sequence[str],
    variables: Sequence[str] = DEFAULT_VARIABLES,
) -> xr.Dataset:
    """Assemble Figure 3 multi-case diagnostic dataset from checkpoint files."""
    return spatial_cache.load_products(pair_paths, cases, variables)


def plot_spatial_mapping(
    dataset: xr.Dataset,
    target: str | Path | Sequence[str],
    destination: str | Path | None = None,
    *,
    variables: Sequence[str] | None = None,
    path: str | Path | None = None,
    experiments: Sequence[str] | None = None,
    data_variable: str | None = None,
    season: str | None = None,
    color_limits: Mapping[str, float] | float | None = None,
    color_levels: Sequence[float] | None = None,
    cmap: str = 'RdBu_r',
    map_projection: str = 'platecarree',
    closed_axis: bool = True,
    font_size: float = 12,
    figsize: tuple[float, float] | None = None,
    shared_colorbar: bool | None = None,
    colorbar_label: str | None = None,
    coastline_linewidth: float = 0.5,
    annotate_statistics: bool = True,
    display_scales: Mapping[str, float] | None = None,
    unit_overrides: Mapping[str, str] | None = None,
    longitude_ticks: Sequence[float] = (-180, -120, -60, 0, 60, 120, 180),
    latitude_ticks: Sequence[float] = (-90, -60, -30, 0, 30, 60, 90),
    title: str | None = None,
    case_label_formatter: Callable[[str], str] | None = None,
) -> Path:
    """Publication-quality 2D horizontal mapping plot following Figure 6's aesthetics.

    Features:
    - Left label column for experiment names (with clean multiline breaks).
    - Variable titles along the top row.
    - Closed-axis frames (all 4 spines visible, black 0.8 pt outline).
    - Dedicated per-column horizontal colorbars with specific units, or a single
      shared colorbar across all columns (e.g. for relative improvement %).
    - Corner metric annotations for RMSE, PCC, Improve %, Degrade %.
    - PlateCarree or Robinson map projections with clean tick degree formatting.

    Parameters
    ----------
    dataset : xr.Dataset
        Dataset containing spatial fields. Supported structures:
        - Figure 3 format: dims ('experiment', 'variable', 'lat', 'lon') with
          data_vars such as 'improvement_percent', 'rmse_reduction', etc.
        - Multi-case 2D fields: dims ('experiment', 'lat', 'lon') with variables
          named after the requested variables (e.g. 'T500', 'PRECT', etc.).
    path : str or Path
        Destination image path.
    variables : Sequence[str], default DEFAULT_VARIABLES
        Variables to plot across columns.
    experiments : Sequence[str], optional
        List/order of experiments to plot. If None, inferred from dataset.
    data_variable : str, optional
        Data variable to plot if dataset has Figure 3 structure
        (e.g., 'improvement_percent' or 'rmse_reduction').
    color_limits : Mapping[str, float] or float, optional
        Symmetric color limit [-limit, limit] per variable or single float.
    color_levels : Sequence[float], optional
        Explicit discrete contour / boundary levels (e.g. for improvement_percent).
    cmap : str, default 'RdBu_r'
        Colormap name.
    map_projection : str, default 'platecarree'
        Map projection ('platecarree' or 'robinson').
    closed_axis : bool, default True
        Whether to enclose map subplots with complete closed spines.
    font_size : float, default 12
        Base typography size.
    figsize : tuple[float, float], optional
        Figure dimensions (width, height). Auto-scaled if None.
    shared_colorbar : bool, optional
        If True, renders one unified colorbar at the bottom. If False, renders
        dedicated column colorbars per variable. If None, auto-detected:
        True if color_levels is set or data_variable is 'improvement_percent',
        False if plotting variables with differing physical units.
    colorbar_label : str, optional
        Custom label for shared colorbar.
    coastline_linewidth : float, default 0.5
        Linewidth for coastlines.
    annotate_statistics : bool, default True
        Whether to display RMSE, PCC, and grid fractions in lower-left corner bbox.
    display_scales : Mapping[str, float], optional
        Multiplier overrides per variable.
    unit_overrides : Mapping[str, str], optional
        Unit string overrides per variable.
    longitude_ticks : Sequence[float]
        Longitude ticks for PlateCarree gridlines.
    latitude_ticks : Sequence[float]
        Latitude ticks for PlateCarree gridlines.
    title : str, optional
        Overall figure suptitle.
    case_label_formatter : Callable[[str], str], optional
        Formatter function for experiment labels on the left axis.

    Returns
    -------
    Path
        Path to the saved figure.
    """
    import cartopy.crs as ccrs

    set_publication_style(closed_spines=closed_axis)
    plt.rcParams.update({'font.size': font_size})

    proj_lower = map_projection.lower().strip()
    if proj_lower in {'platecarree', 'latlon', 'lat-lon'}:
        projection = ccrs.PlateCarree()
    elif proj_lower == 'robinson':
        projection = ccrs.Robinson()
    else:
        raise ValueError(f"Unsupported map_projection: {map_projection!r}")
    data_transform = ccrs.PlateCarree()

    # Resolve arguments between:
    # plot_spatial_mapping(dataset, path, variables=...) and
    # plot_spatial_mapping(dataset, variables, path, ...)
    if isinstance(target, (str, Path)):
        dest_path = Path(target)
        var_list = list(variables if variables is not None else DEFAULT_VARIABLES)
    elif isinstance(target, (Sequence, tuple, list)):
        var_list = list(target)
        if destination is not None:
            dest_path = Path(destination)
        elif path is not None:
            dest_path = Path(path)
        else:
            raise ValueError("A destination path must be provided")
    else:
        dest_path = Path(path if path is not None else destination)
        var_list = list(variables if variables is not None else DEFAULT_VARIABLES)

    # Handle seasonal selection if dataset has a season dimension
    sel_season = season
    if 'season' in dataset.dims or 'season' in dataset.coords:
        if sel_season is None:
            if 'ANN' in dataset.season.values:
                sel_season = 'ANN'
            else:
                sel_season = str(dataset.season.values[0])
        dataset = dataset.sel(season=sel_season)

    # Determine experiments list
    if experiments is not None:
        case_names = list(experiments)
    elif 'experiment' in dataset.coords or 'experiment' in dataset.dims:
        case_names = list(dataset.experiment.values)
    else:
        case_names = ['Model']

    nrows = len(case_names)
    ncols = len(var_list)

    if nrows == 0 or ncols == 0:
        raise ValueError("Must have at least one experiment and one variable to plot")

    # Determine colorbar sharing
    if shared_colorbar is None:
        if color_levels is not None or data_variable == 'improvement_percent':
            shared_colorbar = True
        else:
            shared_colorbar = False

    # Default figure size
    if figsize is None:
        figsize = (2.2 + 3.8 * ncols, 2.1 * nrows + 1.6)

    figure = plt.figure(figsize=figsize, constrained_layout=True)

    # GridSpec: left label column + ncols data columns; nrows rows + spacer row + 1 colorbar row
    grid = figure.add_gridspec(
        nrows + 2, ncols + 1,
        width_ratios=[0.58] + [1.0] * ncols,
        height_ratios=[1.0] * nrows + [0.06, 0.08],
        hspace=0.04, wspace=0.04,
    )

    formatter = case_label_formatter or format_case_label
    axes = np.empty((nrows, ncols), dtype=object)

    # Create subplots and left labels
    for row, experiment in enumerate(case_names):
        label_axis = figure.add_subplot(grid[row, 0])
        label_axis.set_axis_off()
        case_label = formatter(experiment)
        label_axis.text(
            0.98, 0.5, case_label,
            ha='right', va='center',
            fontsize=font_size, fontweight='semibold', linespacing=1.4,
        )
        for col in range(ncols):
            axes[row, col] = figure.add_subplot(grid[row, col + 1], projection=projection)

    scale_map = dict(VARIABLE_SCALES)
    if display_scales:
        scale_map.update(display_scales)
    unit_map = dict(VARIABLE_UNITS)
    if unit_overrides:
        unit_map.update(unit_overrides)

    # Shared norm / levels if applicable
    shared_norm = None
    shared_cmap = cmap
    if shared_colorbar and color_levels is not None:
        clevels = np.asarray(color_levels, dtype=float)
        shared_cmap = plt.colormaps.get_cmap(cmap).resampled(len(clevels) - 1)
        shared_norm = BoundaryNorm(clevels, shared_cmap.N)

    last_panel = None
    col_panels = [None] * ncols
    col_limits = [None] * ncols
    col_units = [None] * ncols

    for col, var in enumerate(var_list):
        unit = unit_map.get(var, '')
        factor = scale_map.get(var, 1.0)

        # Determine color limits for this column if not using shared levels
        var_limit = None
        if not shared_norm:
            if isinstance(color_limits, Mapping):
                var_limit = color_limits.get(var)
            elif isinstance(color_limits, (int, float)):
                var_limit = float(color_limits)
            if var_limit is None:
                var_limit = DEFAULT_BIAS_LIMITS.get(var, 5.0)

        col_limits[col] = var_limit

        for row, experiment in enumerate(case_names):
            axis = axes[row, col]

            # Extract 2D field
            if data_variable is not None and data_variable in dataset:
                # Figure 3 format
                sub_ds = dataset[data_variable]
                if 'experiment' in sub_ds.dims and 'variable' in sub_ds.dims:
                    field = sub_ds.sel(experiment=experiment, variable=var).squeeze()
                elif 'experiment' in sub_ds.dims:
                    field = sub_ds.sel(experiment=experiment).squeeze()
                else:
                    field = sub_ds.squeeze()
            elif var in dataset:
                sub_da = dataset[var]
                if 'experiment' in sub_da.dims:
                    field = sub_da.sel(experiment=experiment).squeeze()
                else:
                    field = sub_da.squeeze()
            else:
                raise KeyError(f"Cannot find variable {var!r} or {data_variable!r} in dataset")

            # Apply scale factor only if input units are raw m/s or user explicitly overrides
            raw_units = field.attrs.get('units', '').strip().lower()
            needs_scaling = (
                (raw_units in {'m/s', 'ms-1', 'ms^-1', 'ms**-1', 'm s**-1', 'm s^-1', 'm s-1'})
                or (display_scales is not None and var in display_scales)
            )
            if needs_scaling and factor != 1.0:
                field = field * factor

            if not unit and 'units' in field.attrs:
                unit = field.attrs['units']
            col_units[col] = unit

            # Map boundaries & coastlines
            axis.set_global()
            axis.coastlines(linewidth=coastline_linewidth, color='0.2')

            # Closed spines
            if closed_axis:
                for spine in axis.spines.values():
                    spine.set_visible(True)
                    spine.set_color('black')
                    spine.set_linewidth(0.8)

            # Gridlines
            if proj_lower == 'platecarree':
                gl = axis.gridlines(
                    crs=data_transform,
                    draw_labels=True,
                    xlocs=list(longitude_ticks),
                    ylocs=list(latitude_ticks),
                    linewidth=0.35,
                    color='0.4',
                    alpha=0.3,
                    linestyle=':',
                )
                gl.top_labels = False
                gl.right_labels = False
                gl.bottom_labels = (row == nrows - 1)
                gl.left_labels = (col == 0)
                gl.xlabel_style = {'size': font_size - 1}
                gl.ylabel_style = {'size': font_size - 1}

            # Plotting the 2D field
            plot_kwargs: dict[str, Any] = {
                'transform': data_transform,
                'cmap': shared_cmap,
            }
            if shared_norm is not None:
                plot_kwargs['norm'] = shared_norm
            elif var_limit is not None:
                plot_kwargs['vmin'] = -var_limit
                plot_kwargs['vmax'] = var_limit

            panel = axis.pcolormesh(field.lon, field.lat, field.values, **plot_kwargs)
            last_panel = panel
            col_panels[col] = panel

            # Top row title
            if row == 0:
                axis.set_title(var, fontsize=font_size + 2, fontweight='semibold', pad=9)

            # Annotation box (RMSE, PCC, Improve%, Degrade%)
            if annotate_statistics:
                # Check if Figure 3 metrics exist in dataset
                stat_names = {
                    'global_rmse': ('RMSE', '{:.3g}'),
                    'global_pcc': ('PCC', '{:.2f}'),
                    'improved_grid_fraction': ('Improve', '{:.1f}%'),
                    'degraded_grid_fraction': ('Degrade', '{:.1f}%'),
                }
                found_stats = {}
                for key, (short_name, fmt) in stat_names.items():
                    if key in dataset:
                        sub_stat = dataset[key]
                        try:
                            if 'experiment' in sub_stat.dims and 'variable' in sub_stat.dims:
                                val = float(sub_stat.sel(experiment=experiment, variable=var).squeeze())
                            elif 'experiment' in sub_stat.dims:
                                val = float(sub_stat.sel(experiment=experiment).squeeze())
                            else:
                                val = float(sub_stat.squeeze())
                            if np.isfinite(val):
                                found_stats[short_name] = fmt.format(val)
                        except (KeyError, ValueError):
                            pass

                # If precomputed RMSE & PCC exist (e.g. from compute_spatial_bias), use them
                if not found_stats and f'{var}_rmse' in dataset and f'{var}_pcc' in dataset:
                    try:
                        r_sub = dataset[f'{var}_rmse']
                        p_sub = dataset[f'{var}_pcc']
                        if 'experiment' in r_sub.dims:
                            r_val = float(r_sub.sel(experiment=experiment).squeeze())
                            p_val = float(p_sub.sel(experiment=experiment).squeeze())
                        else:
                            r_val = float(r_sub.squeeze())
                            p_val = float(p_sub.squeeze())
                        if np.isfinite(r_val):
                            found_stats['RMSE'] = f'{r_val:.3g}'
                        if np.isfinite(p_val):
                            found_stats['PCC'] = f'{p_val:.2f}'
                    except Exception:
                        pass
                # If model/reference paired fields exist in dataset, compute RMSE & PCC
                elif not found_stats and f'{var}_model' in dataset and f'{var}_reference' in dataset:
                    try:
                        m_field = dataset[f'{var}_model'].sel(experiment=experiment).squeeze()
                        r_field = dataset[f'{var}_reference'].sel(experiment=experiment).squeeze()
                        r_val, p_val = spatial_field_statistics(m_field, r_field)
                        if np.isfinite(r_val):
                            found_stats['RMSE'] = f'{r_val:.3g}'
                        if np.isfinite(p_val):
                            found_stats['PCC'] = f'{p_val:.2f}'
                    except Exception:
                        pass

                if found_stats:
                    first_line = '  '.join(
                        f'{name}={found_stats[name]}' for name in ('RMSE', 'PCC')
                        if name in found_stats
                    )
                    second_line = '  '.join(
                        f'{name}={found_stats[name]}' for name in ('Improve', 'Degrade')
                        if name in found_stats
                    )
                    annotation = '\n'.join(line for line in (first_line, second_line) if line)
                    axis.text(
                        0.02, 0.035, annotation,
                        transform=axis.transAxes,
                        ha='left', va='bottom',
                        fontsize=0.74 * font_size,
                        bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.88, 'pad': 1.8},
                        zorder=10,
                    )

    # Colorbar handling
    if shared_colorbar:
        colorbar_axis = figure.add_subplot(grid[-1, 1:])
        cb_kwargs: dict[str, Any] = {
            'cax': colorbar_axis,
            'orientation': 'horizontal',
        }
        if color_levels is not None:
            cb_kwargs.update({
                'boundaries': list(color_levels),
                'ticks': list(color_levels),
                'extend': 'both',
            })
        cb = figure.colorbar(last_panel, **cb_kwargs)
        if closed_axis:
            cb.outline.set_visible(True)
            cb.outline.set_linewidth(0.8)
            cb.outline.set_edgecolor('black')
        label_text = colorbar_label or (
            'RMSE improvement relative to CTRL [%]' if data_variable == 'improvement_percent'
            else 'Value'
        )
        cb.set_label(label_text, fontsize=font_size)
        cb.ax.tick_params(labelsize=font_size - 1, length=3)
    else:
        # Dedicated column colorbars for each variable
        for col, var in enumerate(var_list):
            panel = col_panels[col]
            limit = col_limits[col]
            unit = col_units[col]
            colorbar_axis = figure.add_subplot(grid[-1, col + 1])
            cb = figure.colorbar(
                panel, cax=colorbar_axis, orientation='horizontal',
                ticks=np.linspace(-limit, limit, 5) if limit else None,
                format='%.2g',
            )
            if closed_axis:
                cb.outline.set_visible(True)
                cb.outline.set_linewidth(0.8)
                cb.outline.set_edgecolor('black')
            unit_str = f' [{unit}]' if unit else ''
            cb_label = f'{var}{unit_str}'
            cb.set_label(cb_label, fontsize=font_size)
            cb.ax.tick_params(labelsize=font_size - 1, length=3)

    if title:
        figure.suptitle(title, fontsize=font_size + 2)
    elif sel_season:
        figure.suptitle(f'{sel_season} Mean 2D Bias Relative to Reference', fontsize=font_size + 2)

    return save_figure(figure, dest_path)

