"""Figures of the ``hind_lead_eval`` workflow (styles from ``configs.fcst_init_experiment_config``)."""
from __future__ import annotations

import string

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from configs.fcst_init_experiment_config import STYLE_MAP
from util.fcst_init_skill import WEEK_LABELS

_LETTERS = string.ascii_lowercase


def _style(key):
    return STYLE_MAP.get(key, dict(color="0.3"))


def plot_lead_curves(scores: dict, variables, *, region: str = "GLOBAL", metrics=("rmse", "spread_ratio"),
                     max_lead: int = 56, season: str = "") -> plt.Figure:
    """Rows: variables; columns: metrics. ``scores[(config, variable)]`` = ``lead_scores`` dataset.

    RMSE panels also show each configuration's spread (dashed).
    """
    configs = list(dict.fromkeys(c for c, _ in scores))
    variables = [v for v in variables if any((c, v) in scores for c in configs)]
    fig, axes = plt.subplots(len(variables), len(metrics), figsize=(5.2 * len(metrics), 3.0 * len(variables)),
                             squeeze=False, constrained_layout=True)
    k = 0
    for i, v in enumerate(variables):
        for j, m in enumerate(metrics):
            ax = axes[i, j]
            for c in configs:
                if (c, v) not in scores:
                    continue
                ds = scores[c, v].sel(region=region).sel(lead=slice(1, max_lead))
                ax.plot(ds.lead, ds[m], label=c, **_style(c))
                if m == "rmse":
                    ax.plot(ds.lead, ds["spread"], **{**_style(c), "linestyle": "--", "linewidth": 1.2})
            if m == "spread_ratio":
                ax.axhline(1, color="k", linewidth=0.6)
            units = next(iter(scores[c, v].attrs.get("units", "") for c in configs if (c, v) in scores), "")
            ax.set_title(f"({_LETTERS[k]}) {v} {m}" + (f" [{units}]" if m in ("rmse", "bias", "crps_fair") and units else ""),
                         loc="left")
            ax.set_xlabel("Lead day")
            k += 1
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle(f"{season} daily ensemble performance, {region} (dashed: spread)")
    return fig


def plot_weekly_regional(table: pd.DataFrame, variables, *, region: str = "GLOBAL",
                         metrics=("pattern_acc", "msss"), season: str = "") -> plt.Figure:
    """Weekly regional scores from a long table (config, season, variable, region, week, metric columns)."""
    t = table[(table.region == region) & (table.season == season)] if season else table[table.region == region]
    variables = [v for v in variables if v in set(t.variable)]
    fig, axes = plt.subplots(len(variables), len(metrics), figsize=(4.8 * len(metrics), 2.8 * len(variables)),
                             squeeze=False, constrained_layout=True)
    k = 0
    for i, v in enumerate(variables):
        for j, m in enumerate(metrics):
            ax = axes[i, j]
            for c, g in t[t.variable == v].groupby("config", sort=False):
                ax.plot(g.week, g[m], marker="o", label=c, **_style(c))
            ax.axhline(0, color="k", linewidth=0.6)
            ax.set_xticks(range(1, 9))
            ax.set_title(f"({_LETTERS[k]}) {v} {m}", loc="left")
            ax.set_xlabel("Lead week")
            k += 1
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle(f"{season or 'All starts'}: weekly anomaly skill, {region}")
    return fig


def _map_axes(nrows, ncols, width=4.2, height=2.3):
    import cartopy.crs as ccrs

    fig, axes = plt.subplots(nrows, ncols, figsize=(width * ncols, height * nrows + 0.8), squeeze=False,
                             subplot_kw={"projection": ccrs.PlateCarree(central_longitude=180)}, constrained_layout=True)
    for ax in axes.flat:
        ax.coastlines(linewidth=0.5, color="0.25")
        ax.set_global()
    return fig, axes


def _shade(ax, m: xr.DataArray, levels, cmap):
    """Discrete-level shading of a (lat, lon) map. pcolormesh, not contourf: cartopy takes minutes
    to project the contour polygons of a noisy 1-degree field."""
    import cartopy.crs as ccrs
    from matplotlib.colors import BoundaryNorm

    norm = BoundaryNorm(levels, plt.get_cmap(cmap).N, extend="both")
    return ax.pcolormesh(m.lon, m.lat, m, cmap=cmap, norm=norm, shading="nearest", rasterized=True,
                         transform=ccrs.PlateCarree())


def plot_week_maps(maps: xr.DataArray, *, title: str, weeks=(1, 2, 3, 4, 5, 6, 7, 8), levels=None,
                   cmap="RdBu_r", label="", stipple: xr.DataArray | None = None) -> plt.Figure:
    """2 x 4 panels of a (week, lat, lon) map (ESP-Lab s2s_skill layout); ``stipple``: bool mask."""
    import cartopy.crs as ccrs

    weeks = [w for w in weeks if w in maps.week.values]
    fig, axes = _map_axes(2, 4)
    levels = np.linspace(-1, 1, 21) if levels is None else levels
    cf = None
    for ax, w in zip(axes.flat, weeks):
        cf = _shade(ax, maps.sel(week=w), levels, cmap)
        if stipple is not None:
            s = stipple.sel(week=w).isel(lat=slice(None, None, 4), lon=slice(None, None, 4))
            lo, la = np.meshgrid(s.lon.values, s.lat.values)
            hit = np.asarray(s.values, dtype=bool)
            ax.plot(lo[hit], la[hit], ".", color="k", markersize=1.0, alpha=0.4, transform=ccrs.PlateCarree())
        ax.set_title(WEEK_LABELS[w], fontsize=10)
    for ax in list(axes.flat)[len(weeks):]:
        ax.set_visible(False)
    if cf is not None:
        fig.colorbar(cf, ax=axes, orientation="horizontal", shrink=0.5, pad=0.03, label=label)
    fig.suptitle(title)
    return fig


def plot_config_week_maps(maps: dict[str, xr.DataArray], *, week: int, title: str, levels=None, cmap="RdBu_r",
                          label="") -> plt.Figure:
    """One lead week, one panel per configuration."""
    keys = list(maps)
    ncols = min(3, len(keys))
    nrows = int(np.ceil(len(keys) / ncols))
    fig, axes = _map_axes(nrows, ncols)
    levels = np.linspace(-1, 1, 21) if levels is None else levels
    cf = None
    for k, (ax, c) in enumerate(zip(axes.flat, keys)):
        cf = _shade(ax, maps[c].sel(week=week), levels, cmap)
        ax.set_title(f"({_LETTERS[k]}) {c}", loc="left", fontsize=10)
    for ax in list(axes.flat)[len(keys):]:
        ax.set_visible(False)
    if cf is not None:
        fig.colorbar(cf, ax=axes, orientation="horizontal", shrink=0.5, pad=0.03, label=label)
    fig.suptitle(title)
    return fig


def plot_scorecard(change: pd.DataFrame, *, title: str, higher_is_better: bool, fmt="{:+.2f}", vmax=None,
                   units="") -> plt.Figure:
    """Heat table: rows configurations, columns (variable, week); blue = better than the baseline."""
    vals = change.values.astype(float)
    vmax = vmax or np.nanmax(np.abs(vals)) or 1.0
    cmap = plt.get_cmap("RdBu" if higher_is_better else "RdBu_r")
    fig, ax = plt.subplots(figsize=(0.55 * vals.shape[1] + 2.5, 0.5 * vals.shape[0] + 1.8), constrained_layout=True)
    im = ax.imshow(vals, cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto")
    for (i, j), x in np.ndenumerate(vals):
        if np.isfinite(x):
            ax.text(j, i, fmt.format(x), ha="center", va="center", fontsize=7)
    ax.set_yticks(range(len(change.index)), change.index)
    cols = [f"{a}\nW{b}" if isinstance(a, str) else str(a) for a, b in change.columns] \
        if isinstance(change.columns, pd.MultiIndex) else list(map(str, change.columns))
    ax.set_xticks(range(len(cols)), cols, fontsize=7)
    fig.colorbar(im, ax=ax, shrink=0.8, label=units)
    ax.set_title(title)
    return fig


def plot_rmm_scores(scores: dict[str, xr.Dataset], *, max_lead: int = 56, season: str = "",
                    cor_ref: float = 0.5) -> plt.Figure:
    """Bivariate COR and RMSE (ensemble mean) and spread by lead, mean over the start dates."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)
    for c, ds in scores.items():
        d = ds.sel(lead=slice(1, max_lead))
        axes[0].plot(d.lead, d.cor_ensmean, label=c, **_style(c))
        axes[1].plot(d.lead, d.rmse_ensmean, label=c, **_style(c))
        axes[1].plot(d.lead, d.spread, **{**_style(c), "linestyle": "--", "linewidth": 1.2})
        axes[2].plot(d.lead, d.amp_error_ensmean, label=c, **_style(c))
    axes[0].axhline(cor_ref, color="k", linewidth=0.6)
    axes[1].axhline(np.sqrt(2), color="k", linewidth=0.6)
    axes[2].axhline(0, color="k", linewidth=0.6)
    for ax, t in zip(axes, ("(a) Bivariate COR", "(b) Bivariate RMSE (dashed: spread)", "(c) Amplitude error")):
        ax.set_title(t, loc="left")
        ax.set_xlabel("Lead day")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle(f"{season} RMM skill of the ensemble mean, pooled over start dates")
    return fig
