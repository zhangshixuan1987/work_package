"""Daily top-layer soil-moisture skill of an ensemble against ESA CCI.

Over a verification window (the 16-day DA scorecard window) each land cell
gets, from the ensemble-mean daily series and the ESA CCI daily retrievals on
the days both are valid:

* ``r``: anomaly correlation, i.e. the correlation of departures from each
  series' own window mean (day-to-day wetting and drying);
* ``ubrmse``: unbiased RMSE, the RMSE after removing the window-mean
  difference (the standard soil-moisture DA metric; insensitive to the
  model-satellite climatological offset);
* ``bias``: window-mean model minus ESA CCI, for context;
* ``n``: number of valid days.

Model "daily" values are the ELM daily time series (one value per day at
09 UTC); soil moisture is the thickness-weighted mean of the layers within the
top ``depth_cm`` (two layers, 4.5 cm, for 5 cm), in m3 m-3 like ESA CCI.
"""

from __future__ import annotations

from pathlib import Path

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import numpy as np
import pandas as pd
import xarray as xr

from configs.input_paths import product_dir

REGIONS = {
    "Global": (-60, 90),
    "Tropics": (-20, 20),
    "NH extratropics": (20, 90),
    "SH extratropics": (-60, -20),
}
METRICS = ("r", "ubrmse", "bias")


def top_layer_weights(dz: xr.DataArray, depth_cm: float = 5.0) -> xr.DataArray:
    """Normalised thickness weights of the layers whose bottom lies within ``depth_cm``."""
    top = dz.where(dz.cumsum("levgrnd") <= depth_cm / 100.0 + 1e-6)
    return top / top.sum("levgrnd")


def load_ensemble_top_sm(config, days, weights, *, log=print) -> tuple[xr.DataArray, xr.DataArray]:
    """Ensemble mean and spread of daily top-layer soil moisture on ``days``.

    ``config`` has ``path`` (case directory) and ``nens``; set
    ``member_dirs=False`` for a single-instance case whose archive sits directly
    under ``path`` (e.g. the free-running CTRL runs). The spread is NaN for one
    member.
    """
    year = days[0].year
    members = []
    for m in range(1, config["nens"] + 1):
        case = Path(config["path"]) / (f"EN{m:02d}" if config.get("member_dirs", True) else "")
        path = case / "archive/post/lnd/180x360_aave/ts/daily" / f"H2OSOI.EN{m:02d}.{year}.nc"
        with xr.open_dataset(path) as ds:
            field = ds["H2OSOI"].isel(levgrnd=slice(0, weights.sizes["levgrnd"]))
            field = field.assign_coords(time=field.time.dt.floor("D")).sel(time=days.values)
            members.append((field * weights).sum("levgrnd", min_count=1).load())
    log(f"  loaded {config['nens']} members")
    ensemble = xr.concat(members, "ens")
    return ensemble.mean("ens"), ensemble.std("ens", ddof=1)


def load_esa_cci(reference_root, days) -> xr.DataArray:
    path = product_dir("ESA_CCI", reference_root) / "daily" / f"H2OSOI.daily.{days[0].year}.nc"
    with xr.open_dataset(path) as ds:
        obs = ds["H2OSOI"]
        return obs.assign_coords(time=obs.time.dt.floor("D")).sel(time=days.values).load()


def skill_maps(model: xr.DataArray, obs: xr.DataArray, *, min_days: int = 8) -> xr.Dataset:
    """Per-cell anomaly correlation, ubRMSE, bias over the days both are valid."""
    model, obs = xr.align(model, obs, join="inner")
    valid = model.notnull() & obs.notnull()
    m, o = model.where(valid), obs.where(valid)
    n = valid.sum("time")
    ma, oa = m - m.mean("time"), o - o.mean("time")
    m_std, o_std = np.sqrt((ma ** 2).mean("time")), np.sqrt((oa ** 2).mean("time"))
    ok = (n >= min_days) & (o_std > 1e-4) & (m_std > 1e-5)   # correlation undefined for flat series
    out = xr.Dataset({
        "r": (ma * oa).mean("time") / (m_std * o_std),
        "ubrmse": np.sqrt(((ma - oa) ** 2).mean("time")),
        "bias": (m - o).mean("time"),
        "obs_std": o_std,
        "model_std": m_std,
    }).where(ok)
    out["n"] = n.where(ok)
    for name in ("ubrmse", "bias", "obs_std", "model_std"):
        out[name].attrs["units"] = "m3 m-3"
    return out


def regional_summary(maps: xr.Dataset, regions=REGIONS) -> pd.DataFrame:
    """cos-latitude weighted regional means of the per-cell metrics."""
    rows = []
    weights = np.cos(np.deg2rad(maps.lat))
    for region, (south, north) in regions.items():
        band = (maps.lat >= south) & (maps.lat <= north)
        row = {"region": region, "cells": int(maps["r"].where(band).notnull().sum())}
        for name in (*METRICS, "obs_std", "model_std"):
            field = maps[name].where(band)
            w = weights.broadcast_like(field).where(field.notnull())
            row[name] = float((field * w).sum() / w.sum())
        row["positive_r_pct"] = float(100 * ((maps["r"] > 0).where(maps["r"].notnull()) * weights)
                                      .where(band).sum() / weights.broadcast_like(maps["r"])
                                      .where(maps["r"].notnull() & band).sum())
        rows.append(row)
    return pd.DataFrame(rows)


MAP_EXTENT = (-180, 180, -60, 80)


def _map(ax, field, levels, cmap, title, fontsize, extend, land=None):
    """Land-only map on a rectangular grid; land without data is shaded grey."""
    if land is not None:   # grey underlay marks land with no verification data
        ax.pcolormesh(land.lon, land.lat, land.where(land > 0.5), transform=ccrs.PlateCarree(),
                      cmap=plt.matplotlib.colors.ListedColormap(["0.85"]), shading="auto", rasterized=True)
    cmap = plt.get_cmap(cmap)
    image = ax.pcolormesh(field.lon, field.lat, field, transform=ccrs.PlateCarree(), cmap=cmap,
                          norm=BoundaryNorm(levels, cmap.N, extend=extend), shading="auto", rasterized=True)
    ax.set_extent(MAP_EXTENT, crs=ccrs.PlateCarree())
    ax.coastlines(linewidth=0.4, color="0.3")
    ax.spines["geo"].set_linewidth(0.6)
    ax.set_title(title, loc="left", fontsize=1.125 * fontsize, pad=4)
    return image


def plot_sm_skill(maps: dict[str, xr.Dataset], summary: pd.DataFrame, *, labels: dict[str, str] | None = None,
                  colors: dict[str, str] | None = None, fontsize: float = 12, map_experiments=None,
                  land: xr.DataArray | None = None,
                  r_levels=(-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8),
                  ubrmse_levels=(0, 0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06)):
    """Maps of r and ubRMSE for ``map_experiments`` (columns) and regional-mean bars.

    ``land`` (land fraction) shades land cells without data grey. On the
    ubRMSE bars a black tick marks the observed day-to-day variability (ESA
    CCI std about its window mean): ubRMSE near it means little skill.
    """
    labels, colors = labels or {}, colors or {}
    names = list(map_experiments or maps)
    bar_names = list(maps)
    scale = fontsize / 12
    map_w = 7.0 * scale
    map_h = map_w * (MAP_EXTENT[3] - MAP_EXTENT[2]) / (MAP_EXTENT[1] - MAP_EXTENT[0])
    fig = plt.figure(figsize=(len(names) * map_w + 1.3 * scale, 2 * map_h + 4.6 * scale), layout="constrained")
    fig.get_layout_engine().set(h_pad=0.06 * scale, w_pad=0.06 * scale, hspace=0.02, wspace=0.02)
    grid = fig.add_gridspec(3, len(names) + 1, width_ratios=(*[1] * len(names), 0.035),
                            height_ratios=(map_h, map_h, 3.4 * scale))
    letters = iter("abcdefgh")
    rows = [("r", r_levels, "RdBu_r", "both", "Anomaly correlation"),
            ("ubrmse", ubrmse_levels, "YlOrRd", "max", "ubRMSE (m$^3$ m$^{-3}$)")]
    for i, (metric, levels, cmap, extend, label) in enumerate(rows):
        for j, name in enumerate(names):
            ax = fig.add_subplot(grid[i, j], projection=ccrs.PlateCarree())
            image = _map(ax, maps[name][metric], list(levels), cmap,
                         f"({next(letters)}) {label.split(' (')[0]}, {labels.get(name, name)}", fontsize, extend, land)
        cbar = fig.colorbar(image, cax=fig.add_subplot(grid[i, -1]), extend=extend, ticks=list(levels),
                            drawedges=True)
        cbar.dividers.set_color("white")
        cbar.outline.set_visible(False)
        cbar.set_label(label, fontsize=fontsize)
        cbar.ax.tick_params(labelsize=0.875 * fontsize)

    regions = list(dict.fromkeys(summary["region"]))
    x = np.arange(len(regions))
    width = 0.8 / len(bar_names)
    sub = grid[2, :].subgridspec(1, 2, wspace=0.08)
    for k, metric in enumerate(("r", "ubrmse")):
        ax = fig.add_subplot(sub[0, k])
        for j, name in enumerate(bar_names):
            rows_ = summary[summary.experiment == name].set_index("region").loc[regions]
            xs = x + (j - (len(bar_names) - 1) / 2) * width
            bars = ax.bar(xs, rows_[metric], width * 0.92, color=colors.get(name, f"C{j}"),
                          label=labels.get(name, name))
            if metric == "r":
                ax.bar_label(bars, fmt="%.2f", padding=2, fontsize=0.75 * fontsize)
            else:   # mid-bar, clear of the variability tick near its top
                ax.bar_label(bars, fmt="%.3f", label_type="center", color="white", fontsize=0.75 * fontsize)
            if metric == "ubrmse":
                ax.hlines(rows_["obs_std"], xs - width * 0.46, xs + width * 0.46, color="black", linewidth=2 * scale,
                          label="ESA CCI day-to-day std" if j == 0 else None)
        ax.set_xticks(x, regions, fontsize=0.875 * fontsize)
        ax.tick_params(axis="y", labelsize=0.875 * fontsize)
        ax.set_ylabel("Anomaly correlation" if metric == "r" else "ubRMSE (m$^3$ m$^{-3}$)", fontsize=fontsize)
        ax.set_ylim(0, None)
        ax.margins(y=0.15)
        ax.grid(axis="y", alpha=0.25)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.set_title(f"({next(letters)}) Land mean, {'anomaly correlation' if metric == 'r' else 'ubRMSE'}",
                     loc="left", fontsize=1.125 * fontsize)
        ax.legend(frameon=False, fontsize=0.8 * fontsize, loc="upper left")
    return fig
