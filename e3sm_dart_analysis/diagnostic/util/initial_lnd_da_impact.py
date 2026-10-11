"""Diagnose how atmospheric DA changes the initialized land state.

Input is the processed soil-moisture snapshot written by
``analysis_lnd_init/02_soil_moisture_distribution.ipynb`` (``model_mean``, ``obs``, ``spread``
with an ``experiment`` dimension). Land is not assimilated, so any difference
between a DART experiment and the no-DA control comes from the atmospheric
forcing the land saw during the assimilation window.

Three aspects are derived relative to the baseline:

* ensemble spread: ratio maps, zonal profiles, and land-mean change;
* ensemble mean: the DA-induced shift (experiment minus baseline);
* error: RMSE against the observation on a mask shared by all experiments.
"""

from __future__ import annotations

from pathlib import Path

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import numpy as np
import pandas as pd
import xarray as xr

REQUIRED = ("model_mean", "obs", "spread")


def _weighted_mean(field: xr.DataArray) -> float:
    """cos-latitude weighted mean over finite cells."""
    weights = np.cos(np.deg2rad(field["lat"])).broadcast_like(field).where(np.isfinite(field))
    return float((field * weights).sum() / weights.sum())


def derive_land_da_impact(
    processed: xr.Dataset,
    *,
    baseline: str = "CTRL10-S0",
    experiments: tuple[str, ...] = ("DART10-S0", "DART20-S0", "DART40-S0"),
    min_spread: float = 0.05,
) -> xr.Dataset:
    """Spread, mean-shift, and error diagnostics of ``experiments`` vs ``baseline``.

    The spread ratio is left undefined where the baseline spread is below
    ``min_spread`` (field units), e.g. ice sheets and deserts, where a ratio of
    two near-zero spreads is noise.
    """
    missing = [name for name in REQUIRED if name not in processed]
    if missing:
        raise KeyError(
            f"Processed cache lacks {missing}; rerun analysis_lnd_init/02_soil_moisture_distribution "
            "with FORCE_COMPUTE = True."
        )
    names = (baseline, *experiments)
    unknown = sorted(set(names) - {str(e) for e in processed.experiment.values})
    if unknown:
        raise KeyError(f"Experiments missing from the processed cache: {unknown}")
    data = processed.sel(experiment=list(names))
    obs = data["obs"].isel(experiment=0, drop=True)
    # One observation mask for every experiment: compare like with like.
    common = obs.notnull() & data["model_mean"].notnull().all("experiment")
    error = (data["model_mean"] - obs).where(common)

    spread = data["spread"]
    land_spread = xr.concat([xr.DataArray(_weighted_mean(spread.sel(experiment=e))) for e in names],
                            pd.Index(names, name="experiment"))
    rmse = xr.concat([xr.DataArray(np.sqrt(_weighted_mean(error.sel(experiment=e) ** 2))) for e in names],
                     pd.Index(names, name="experiment"))
    base_abs = np.abs(error.sel(experiment=baseline))
    improved = xr.concat(
        [xr.DataArray(100 * _weighted_mean((np.abs(error.sel(experiment=e)) < base_abs).astype(float)
                                           .where(common)))
         for e in experiments],
        pd.Index(experiments, name="comparison"))

    def compare(field):
        return field.sel(experiment=list(experiments)).rename(experiment="comparison")

    units = processed["model_mean"].attrs.get("units", "kg/m2")
    impact = xr.Dataset({
        "spread": spread,
        "zonal_spread": spread.mean("lon", skipna=True),
        "land_spread": land_spread,
        "rmse": rmse,
        "spread_ratio": compare(spread) / spread.sel(experiment=baseline, drop=True).where(
            lambda base_spread: base_spread >= min_spread),
        "mean_shift": compare(data["model_mean"]) - data["model_mean"].sel(experiment=baseline, drop=True),
        "spread_change_pct": 100 * (compare(land_spread) / land_spread.sel(experiment=baseline) - 1),
        "rmse_change_pct": 100 * (compare(rmse) / rmse.sel(experiment=baseline) - 1),
        "improved_fraction_pct": improved,
    })
    for name in ("spread", "zonal_spread", "land_spread", "rmse", "mean_shift"):
        impact[name].attrs["units"] = units
    impact.attrs.update(
        baseline=baseline,
        observation=str(processed.attrs.get("observation", "")),
        target_date=str(processed.attrs.get("target_date", "")),
        depth_cm=float(processed.attrs.get("depth_cm", np.nan)),
        units=units,
        n_common_cells=int(common.sum()),
        min_spread=min_spread,
    )
    return impact


def impact_summary(impact: xr.Dataset) -> pd.DataFrame:
    """One row per comparison experiment."""
    base = impact.attrs["baseline"]
    rows = []
    for e in impact.comparison.values:
        rows.append({
            "experiment": str(e),
            "baseline": base,
            "land_spread": float(impact["land_spread"].sel(experiment=e)),
            "baseline_land_spread": float(impact["land_spread"].sel(experiment=base)),
            "spread_change_pct": float(impact["spread_change_pct"].sel(comparison=e)),
            "rmse": float(impact["rmse"].sel(experiment=e)),
            "baseline_rmse": float(impact["rmse"].sel(experiment=base)),
            "rmse_change_pct": float(impact["rmse_change_pct"].sel(comparison=e)),
            "improved_fraction_pct": float(impact["improved_fraction_pct"].sel(comparison=e)),
            "mean_abs_shift": _weighted_mean(np.abs(impact["mean_shift"].sel(comparison=e))),
        })
    return pd.DataFrame(rows)


def _map_panel(fig, spec, field, levels, cmap, label, title, fontsize, extend="both"):
    ax = fig.add_subplot(spec, projection=ccrs.Robinson())
    cmap = plt.get_cmap(cmap)
    norm = BoundaryNorm(levels, cmap.N, extend=extend)
    image = ax.pcolormesh(field.lon, field.lat, field, transform=ccrs.PlateCarree(),
                          cmap=cmap, norm=norm, shading="auto", rasterized=True)
    ax.set_global()
    ax.coastlines(linewidth=0.5)
    ax.set_title(title, loc="left", fontsize=1.125 * fontsize)
    cbar = fig.colorbar(image, ax=ax, orientation="horizontal", pad=0.04, shrink=0.85,
                        extend=extend, ticks=levels, drawedges=True)
    cbar.dividers.set_color("white")
    cbar.outline.set_visible(False)
    cbar.set_label(label, fontsize=fontsize)
    cbar.ax.tick_params(labelsize=0.875 * fontsize)
    return ax


def plot_land_da_impact(
    impact: xr.Dataset,
    *,
    primary_experiment: str = "DART10-S0",
    colors: dict[str, str] | None = None,
    ratio_levels=(0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.4, 2.0, 3.3, 5.0, 10.0),   # ~symmetric in log about 1
    shift_levels=(-3, -2, -1, -0.5, -0.2, 0.2, 0.5, 1, 2, 3),
    fontsize: float = 12,
):
    """2x2 figure: spread ratio map, mean-shift map, zonal spread, land-mean bars."""
    colors = colors or {}
    base = impact.attrs["baseline"]
    comparisons = [str(e) for e in impact.comparison.values]
    experiments = [base, *comparisons]
    color = {name: colors.get(name, f"C{i}") for i, name in enumerate(experiments)}
    units = impact.attrs.get("units", "kg/m2").replace("kg/m2", "kg m$^{-2}$")
    depth = impact.attrs.get("depth_cm", np.nan)
    layer = f"top {depth:g} cm " if np.isfinite(depth) else ""

    scale = fontsize / 12
    fig = plt.figure(figsize=(15 * scale, 10 * scale), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=(1.15, 1))

    _map_panel(fig, grid[0, 0], impact["spread_ratio"].sel(comparison=primary_experiment),
               list(ratio_levels), "RdBu_r", "Spread ratio",
               f"(a) Ensemble spread, {primary_experiment} / {base}", fontsize)
    _map_panel(fig, grid[0, 1], impact["mean_shift"].sel(comparison=primary_experiment),
               list(shift_levels), "BrBG", f"Mean shift ({units})",
               f"(b) Ensemble mean, {primary_experiment} − {base}", fontsize)

    zonal = fig.add_subplot(grid[1, 0])
    for name in experiments:
        line = impact["zonal_spread"].sel(experiment=name)
        zonal.plot(line.lat, line, color=color[name], linewidth=2 * scale,
                   linestyle="--" if name == base else "-", label=name)
    zonal.set_xlim(-60, 90)
    zonal.set_xticks(range(-60, 91, 30))
    zonal.set_xticklabels([f"{abs(v)}°{'N' if v > 0 else 'S' if v < 0 else ''}"
                           for v in range(-60, 91, 30)])
    zonal.set_ylim(bottom=0)
    zonal.set_ylabel(f"Zonal-mean spread ({units})", fontsize=fontsize)
    zonal.tick_params(labelsize=0.875 * fontsize)
    zonal.grid(alpha=0.25)
    zonal.legend(frameon=False, fontsize=0.875 * fontsize)
    zonal.set_title(f"(c) {layer}soil-moisture spread by latitude", loc="left", fontsize=1.125 * fontsize)

    bars = fig.add_subplot(grid[1, 1])
    observation = impact.attrs.get("observation", "obs").replace("_", " ")
    metrics = [("spread_change_pct", "Land-mean spread"), ("rmse_change_pct", f"RMSE vs {observation}")]
    x = np.arange(len(metrics))
    width = 0.8 / len(comparisons)
    for k, name in enumerate(comparisons):
        values = [float(impact[m].sel(comparison=name)) for m, _ in metrics]
        offset = (k - (len(comparisons) - 1) / 2) * width
        drawn = bars.bar(x + offset, values, width, color=color[name], label=name)
        bars.bar_label(drawn, fmt="%+.0f", padding=2, fontsize=0.8 * fontsize)
    bars.axhline(0, color="0.25", linewidth=0.8)
    bars.set_xticks(x, [label for _, label in metrics], fontsize=fontsize)
    bars.set_ylabel(f"Change vs {base} (%)", fontsize=fontsize)
    bars.tick_params(axis="y", labelsize=0.875 * fontsize)
    bars.grid(axis="y", alpha=0.25)
    bars.legend(frameon=False, fontsize=0.875 * fontsize, loc="lower right")
    bars.set_title("(d) Land-mean change", loc="left", fontsize=1.125 * fontsize)
    return fig


def save_impact_products(impact: xr.Dataset, *, data_path: str | Path, csv_path: str | Path) -> pd.DataFrame:
    """Save the derived fields and the scalar summary."""
    data_path, csv_path = Path(data_path), Path(csv_path)
    data_path.parent.mkdir(parents=True, exist_ok=True)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    impact.to_netcdf(data_path)
    summary = impact_summary(impact)
    summary.to_csv(csv_path, index=False)
    return summary
