"""Figures of the ``hind_lead_eval`` workflow (styles from ``configs.fcst_init_experiment_config``).

Layout of the line figures: one figure per metric and region, variables in rows, the start seasons in
columns (Dec 2011 and May 2012 starts side by side, y-axis shared along each row), one legend below the
panels. Every figure is drawn with ``readable_rc()``: all font sizes scale with one base size, ``FONT_SIZE``
(``set_font_size(fontz)`` in a notebook's figure settings).
"""
from __future__ import annotations

import functools
import string

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib import transforms
from matplotlib.lines import Line2D

from configs.fcst_init_experiment_config import STYLE_MAP, config_label
from util.fcst_init_skill import WEEK_LABELS

_LETTERS = string.ascii_lowercase

FONT_SIZE = 18.0   # base font size (pt); every label scales with it
# size of each kind of text relative to FONT_SIZE
FONT_SCALES = dict(label=1.0, title=1.0, tick=0.9, legend=1.0, suptitle=1.15, cell=0.75)


def set_font_size(fontz: float, **scales) -> None:
    """Base font size of every figure drawn afterwards, and optionally the scale factors of
    ``FONT_SCALES`` (label, title, tick, legend, suptitle, cell = scorecard numbers and columns)."""
    global FONT_SIZE
    unknown = set(scales) - set(FONT_SCALES)
    if unknown:
        raise ValueError(f"unknown font scales {sorted(unknown)}; choose from {list(FONT_SCALES)}")
    FONT_SIZE = float(fontz)
    FONT_SCALES.update({k: float(v) for k, v in scales.items()})


def fs(kind: str = "label") -> float:
    """Font size of one kind of text: ``FONT_SCALES[kind]`` x ``FONT_SIZE``."""
    return round(FONT_SCALES[kind] * FONT_SIZE, 1)


LEGEND_BOX = dict(frameon=True, fancybox=False, edgecolor="k", framealpha=1.0)   # boxed legends


def fig_size(width: float, height: float) -> tuple[float, float]:
    """A figure size designed for 12-pt text, grown with ``FONT_SIZE`` so larger text still fits."""
    k = FONT_SIZE / 12.0
    return (width * k, height * k)


def readable_rc() -> dict:
    """Matplotlib settings of the workflow's figures, font sizes from ``FONT_SIZE``."""
    return {
        "font.size": fs("label"), "axes.labelsize": fs("label"), "axes.titlesize": fs("title"),
        "xtick.labelsize": fs("tick"), "ytick.labelsize": fs("tick"), "legend.fontsize": fs("legend"),
        "figure.titlesize": fs("suptitle"), "lines.linewidth": 2.0, "lines.markersize": 5,
        "axes.grid": True, "grid.alpha": 0.3, "grid.linewidth": 0.6,
        "axes.spines.top": True, "axes.spines.right": True,      # closed box around each panel
        "figure.constrained_layout.use": True,                   # labels never overlap at large fonts
    }
SEASON_TITLES = {"DJF": "Dec 2011 starts (DJF)", "JJA": "May 2012 starts (JJA)", "ALL": "All starts"}
SEASON_PANEL = {"DJF": "Dec init.", "JJA": "May init.", "ALL": "All inits"}   # right-hand panel titles
METRIC_LABELS = {
    "rmse": "RMSE", "bias": "Bias", "spread": "Spread", "spread_ratio": "Spread / RMSE", "crps_fair": "CRPS",
    "pattern_acc": "Pattern ACC", "msss": "MSSS", "rmse_anom": "Anomaly RMSE",
    "memory": "SM memory r(week 1, week w)",
}
_UNITLESS = {"spread_ratio", "pattern_acc", "msss", "memory"}


def readable(func):
    """Draw the figure of ``func`` with ``readable_rc()``."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        with plt.rc_context(readable_rc()):
            return func(*args, **kwargs)
    return wrapper


def _style(key):
    return STYLE_MAP.get(key, dict(color="0.3"))


def season_title(season: str) -> str:
    return SEASON_TITLES.get(season, season)


def metric_label(metric: str, units: str = "") -> str:
    if metric.startswith("coupling_"):
        return f"SM coupling r(SM, {metric.removeprefix('coupling_')})"
    label = METRIC_LABELS.get(metric, metric)
    return f"{label} [{units}]" if units and metric not in _UNITLESS else label


def _season_grid(rows, seasons, *, panel=(6.4, 3.2)):   # panel height:width = 1:2
    """Rows x seasons axes, y shared along each row; each panel titled "(a) <row>" on the left and
    "Dec init." / "May init." on the right."""
    fig, axes = plt.subplots(len(rows), len(seasons), figsize=(panel[0] * len(seasons), panel[1] * len(rows) + 0.9),
                             squeeze=False, sharex=True, sharey="row", layout="constrained")
    _panel_titles(axes, [[r] * len(seasons) for r in rows], [list(seasons)] * len(rows))
    return fig, axes


def _panel_titles(axes, left, seasons):
    """``(a) left[i][j]`` on the left of panel (i, j), its season on the right; letters run along the rows."""
    k = 0
    for i in range(axes.shape[0]):
        for j in range(axes.shape[1]):
            axes[i, j].set_title(f"({_LETTERS[k]}) {left[i][j]}", loc="left")
            axes[i, j].set_title(SEASON_PANEL.get(seasons[i][j], seasons[i][j]), loc="right")
            k += 1


def _legend(fig, configs, extra: list[tuple[Line2D, str]] = ()):
    handles = [Line2D([], [], **_style(c)) for c in configs] + [h for h, _ in extra]
    labels = [config_label(c) for c in configs] + [t for _, t in extra]
    fig.legend(handles, labels, loc="outside lower center", ncols=min(len(labels), 6), **LEGEND_BOX)


def _lead_days(ax, max_lead):
    ax.set_xlim(1, max_lead)
    ax.set_xticks(range(7, max_lead + 1, 7))      # weekly ticks


def _weeks(ax, weeks):
    ax.set_xticks(list(weeks))


def _suptitle(fig, text, region_label):
    fig.suptitle(f"{text}, {region_label}" if region_label else text)


# --------------------------------------------------------------------------------------------------
# Lead-day and lead-week curves
# --------------------------------------------------------------------------------------------------

@readable
def plot_lead_metric(scores: dict, variables, *, metric: str, seasons, region: str, region_label: str = "",
                     max_lead: int = 56, show_spread: bool = False) -> plt.Figure:
    """One metric by lead day: rows ``variables``, columns ``seasons``.

    ``scores[(season, config, variable)]`` is a ``lead_scores`` dataset (``region`` dimension).
    ``show_spread`` adds each configuration's spread (dashed) to RMSE panels; the spread/RMSE figure
    shows the same information as a ratio.
    """
    configs = list(dict.fromkeys(c for _, c, _ in scores))
    seasons = [s for s in seasons if any(k[0] == s for k in scores)]
    variables = [v for v in variables if any(k[2] == v for k in scores)]
    fig, axes = _season_grid(variables, seasons)
    for i, v in enumerate(variables):
        units = ""
        for j, s in enumerate(seasons):
            ax = axes[i, j]
            for c in configs:
                if (s, c, v) not in scores:
                    continue
                ds = scores[s, c, v]
                units = units or ds.attrs.get("units", "")
                d = ds.sel(region=region).sel(lead=slice(1, max_lead))
                ax.plot(d.lead, d[metric], **_style(c))
                if show_spread and metric == "rmse":
                    ax.plot(d.lead, d["spread"], **{**_style(c), "linestyle": "--", "linewidth": 1.2})
            if metric == "spread_ratio":
                ax.axhline(1, color="k", linewidth=0.8)
            elif metric == "bias":
                ax.axhline(0, color="k", linewidth=0.8)
            _lead_days(ax, max_lead)
        axes[i, 0].set_ylabel(metric_label(metric, units))
    for ax in axes[-1]:
        ax.set_xlabel("Lead day")
    extra = [(Line2D([], [], color="0.3", linestyle="--", linewidth=1.2), "spread")] \
        if show_spread and metric == "rmse" else []
    _legend(fig, configs, extra)
    _suptitle(fig, f"Daily ensemble {metric_label(metric)} by lead day", region_label or region)
    return fig


def _weekly_panels(fig, axes, frame: pd.DataFrame, variables, seasons, value: str, label: str, zero=True):
    for i, v in enumerate(variables):
        for j, s in enumerate(seasons):
            ax = axes[i, j]
            t = frame[(frame.variable == v) & (frame.season == s)]
            for c, g in t.groupby("config", sort=False):
                g = g.sort_values("week")
                ax.plot(g.week, g[value], marker="o", **_style(c))
            if zero:
                ax.axhline(0, color="k", linewidth=0.8)
            _weeks(ax, range(1, 9))
        axes[i, 0].set_ylabel(label)
    for ax in axes[-1]:
        ax.set_xlabel("Lead week")


@readable
def plot_weekly_metric(table: pd.DataFrame, variables, *, metric: str, seasons, region: str,
                       region_label: str = "", title: str = "") -> plt.Figure:
    """One weekly score (``table``: config, season, variable, region, week, metric columns): rows
    ``variables``, columns ``seasons``."""
    t = table[table.region == region]
    seasons = [s for s in seasons if s in set(t.season)]
    variables = [v for v in variables if v in set(t.variable)]
    fig, axes = _season_grid(variables, seasons)
    _weekly_panels(fig, axes, t, variables, seasons, metric, metric_label(metric))
    _legend(fig, list(dict.fromkeys(t.config)))
    _suptitle(fig, title or f"Weekly anomaly {metric_label(metric)}", region_label or region)
    return fig


@readable
def plot_bias_impact(table: pd.DataFrame, variables, *, metric: str, seasons, region: str, region_label: str = "",
                     biases=("none", "amip")) -> plt.Figure:
    """Weekly ``metric`` under ``biases[1]`` minus ``biases[0]`` (long table with a ``bias`` column).

    With ("none", "amip"): positive = skill gained by removing the AMIP mean bias (the cost of the model
    bias, once the forecast has drifted to the model climate); negative in the first weeks, when the
    forecast is still near the observed state.
    """
    t = table[table.region == region]
    wide = t.pivot_table(index=["season", "variable", "config", "week"], columns="bias", values=metric, sort=False)
    if not set(biases) <= set(wide.columns):
        raise ValueError(f"bias treatments {biases} not all in the table ({list(wide.columns)})")
    diff = (wide[biases[1]] - wide[biases[0]]).rename("delta").reset_index()
    seasons = [s for s in seasons if s in set(diff.season)]
    variables = [v for v in variables if v in set(diff.variable)]
    fig, axes = _season_grid(variables, seasons)
    _weekly_panels(fig, axes, diff, variables, seasons, "delta", f"Δ{metric_label(metric)}")
    _legend(fig, list(dict.fromkeys(diff.config)))
    _suptitle(fig, f"Effect of removing the model mean bias: {metric_label(metric)} ({biases[1]} − {biases[0]})",
              region_label or region)
    return fig


@readable
def plot_regional_weekly(frame: pd.DataFrame, metrics, *, seasons, region: str, region_label: str = "",
                         title: str = "") -> plt.Figure:
    """Weekly regional values of several metrics (``frame``: config, season, region, metric, week, value):
    rows ``metrics``, columns ``seasons``."""
    t = frame[frame.region == region].rename(columns={"metric": "variable"})
    seasons = [s for s in seasons if s in set(t.season)]
    metrics = [m for m in metrics if m in set(t.variable)]
    fig, axes = _season_grid([metric_label(m) for m in metrics], seasons)
    _weekly_panels(fig, axes, t, metrics, seasons, "value", "Correlation")
    _legend(fig, list(dict.fromkeys(t.config)))
    _suptitle(fig, title, region_label or region)
    return fig


LEAD_SKILL_LABELS = {"acc": "ACC", "nrmse": "nRMSE", "acc_mean": "ACC", "nrmse_region": "nRMSE"}
LEAD_SKILL_REFS = {"acc": 0.0, "nrmse": 1.0, "acc_mean": 0.0, "nrmse_region": 1.0}   # reference lines


@readable
def plot_lead_skill(table: pd.DataFrame, *, variable: str, region: str, seasons, metrics=("acc", "nrmse"),
                    sig: str | None = None, title: str = "", region_label: str = "", ylims: dict | None = None,
                    panel=(6.4, 3.2)) -> plt.Figure:
    """Skill by lead week in the layout of ESP-Lab ``2a_regional_acc_skill_ts`` / ``3a_sst_skill_ts``: rows
    ``seasons`` (Dec init., May init., All inits), columns ``metrics`` (ACC left, nRMSE right), one line per
    configuration. ``table``: config, season, variable, region, week and the metric columns; ``sig``: a bool
    column, markers filled where true and open otherwise (3a: ACC significant). ``ylims``: {metric: (lo, hi)}."""
    t = table[(table.variable == variable) & (table.region == region)]
    seasons = [x for x in seasons if x in set(t.season)]
    configs = list(dict.fromkeys(t.config))
    k = FONT_SIZE / 12.0
    fig, axes = plt.subplots(len(seasons), len(metrics), figsize=(panel[0] * k * len(metrics), panel[1] * k * len(seasons) + k),
                             squeeze=False, sharex=True, sharey="col", layout="constrained")
    n = 0
    for i, season in enumerate(seasons):
        for j, m in enumerate(metrics):
            ax = axes[i, j]
            for c in configs:
                d = t[(t.season == season) & (t.config == c)].sort_values("week")
                if d.empty:
                    continue
                st = _style(c)
                ax.plot(d.week, d[m], **{**st, "marker": None})
                filled = d[sig].astype(bool).values if sig else np.ones(len(d), bool)
                col = st.get("color", "0.3")
                ax.plot(d.week[filled], d[m][filled], "o", color=col, markersize=7, linestyle="none")
                ax.plot(d.week[~filled], d[m][~filled], "o", markerfacecolor="w", markeredgecolor=col,
                        markersize=7, linestyle="none")
            if m in LEAD_SKILL_REFS:
                ax.axhline(LEAD_SKILL_REFS[m], color="k", linewidth=0.8)
            if ylims and m in ylims:
                ax.set_ylim(*ylims[m])
            ax.set_xticks(sorted(set(t.week)))
            ax.set_ylabel(LEAD_SKILL_LABELS.get(m, metric_label(m)))
            if i == len(seasons) - 1:
                ax.set_xlabel("Lead time (weeks)")
            ax.set_title(f"({_LETTERS[n]}) {LEAD_SKILL_LABELS.get(m, metric_label(m))}", loc="left")
            ax.set_title(SEASON_PANEL.get(season, season), loc="right")
            n += 1
    extra = []
    if sig:
        extra = [(Line2D([], [], marker="o", color="0.3", linestyle="none"), "ACC significant (95%)"),
                 (Line2D([], [], marker="o", markerfacecolor="w", markeredgecolor="0.3", linestyle="none"), "not significant")]
    _legend(fig, configs, extra)
    _suptitle(fig, title or f"{variable} skill", region_label or region)
    return fig


@readable
def plot_index_series(series: pd.DataFrame, *, variable: str, region: str, windows, seasons, units: str = "",
                      title: str = "", region_label: str = "", panel=(8.0, 4.0)) -> plt.Figure:
    """Regional index by start date at fixed leads (ESP-Lab ``3a_sst_skill_ts`` DJF series): one panel per
    lead window, the starts of ``seasons`` along x (one segment per season), the observed index in black and
    the ensemble mean of each configuration; one box per season lists each configuration's ACC and nRMSE over
    that season's starts. ``series``: config, season, variable, region, init, week, model, obs."""
    t = series[(series.variable == variable) & (series.region == region) & series.season.isin(seasons)]
    configs = list(dict.fromkeys(t.config))
    k = FONT_SIZE / 12.0
    fig, axes = plt.subplots(1, len(windows), figsize=(panel[0] * k * len(windows), panel[1] * k + 1.2 * k),
                             squeeze=False, sharey=True, layout="constrained")
    starts = (t[["season", "init"]].drop_duplicates()
              .assign(order=lambda d: d.season.map({x: n for n, x in enumerate(seasons)})).sort_values(["order", "init"]))
    xpos, x, prev = {}, 0, None
    for season, init in zip(starts.season, starts.init):
        x += 1.0 if prev in (None, season) else 1.8   # gap between seasons
        xpos[season, init], prev = x, season
    nbox = 0
    for j, (ax, w) in enumerate(zip(axes[0], windows)):
        sel = t[t.week.isin(list(w))]
        mean = sel.groupby(["config", "season", "init"], as_index=False)[["model", "obs"]].mean()
        obs = mean.drop_duplicates(["season", "init"])
        for season in seasons:
            o = obs[obs.season == season].sort_values("init")
            ax.plot([xpos[season, i] for i in o.init], o.obs, "-o", color="k", markersize=5, linewidth=1.5, zorder=5)
        boxes = {season: [] for season in seasons}
        for c in configs:
            d = mean[mean.config == c]
            for season in seasons:
                e = d[d.season == season].sort_values("init")
                ax.plot([xpos[season, i] for i in e.init], e.model, **{**_style(c), "marker": "o", "markersize": 5})
                if len(e) > 2:   # skill over the starts of this season (one start month, as ESP-Lab 3a)
                    r = np.corrcoef(e.model, e.obs)[0, 1]
                    nr = np.sqrt(((e.model - e.obs) ** 2).mean() / (e.obs ** 2).mean())
                    boxes[season].append(f"{config_label(c)}: ACC {r:+.2f}, nRMSE {nr:.2f}")
        ax.axhline(0, color="0.5", linewidth=0.8)
        ax.set_xticks(list(xpos.values()), [pd.Timestamp(i).strftime("%d %b") for _, i in xpos], rotation=45, ha="right")
        ax.set_title(f"({_LETTERS[j]}) {window_title(w)} lead", loc="left")
        if j == 0:
            ax.set_ylabel(f"{variable} anomaly" + (f" [{units}]" if units else ""))
        blend = transforms.blended_transform_factory(ax.transData, ax.transAxes)
        for season, lines in boxes.items():
            if lines:
                x0 = min(v for (s_, _), v in xpos.items() if s_ == season) - 0.3
                ax.text(x0, 0.98, "\n".join(lines), transform=blend, ha="left", va="top", fontsize=fs("cell"),
                        bbox=dict(boxstyle="square,pad=0.3", fc="w", ec="0.6", alpha=0.9))
                nbox = max(nbox, len(lines))
        ax.set_xlim(min(xpos.values()) - 0.4, max(xpos.values()) + 0.4)
    if nbox:   # headroom above the data for the skill box
        lo, hi = axes[0, 0].get_ylim()
        frac = min(0.6, (nbox * 1.3 * fs("cell") / 72.0 + 0.15) / (panel[1] * k))
        axes[0, 0].set_ylim(lo, lo + (hi - lo) / (1.0 - frac))
    _legend(fig, configs, [(Line2D([], [], color="k", marker="o", markersize=5), "Observed")])
    _suptitle(fig, title or f"{variable} regional anomaly by start date", region_label or region)
    return fig


def window_title(weeks) -> str:
    weeks = tuple(weeks)
    return f"Week {weeks[0]}" if len(weeks) == 1 else f"Weeks {weeks[0]}-{weeks[-1]}"


@readable
def plot_rmm_scores(scores: dict, *, seasons, max_lead: int = 56, cor_ref: float = 0.5) -> plt.Figure:
    """RMM skill of the ensemble mean: rows ``seasons``; columns bivariate COR, RMSE (dashed: spread),
    amplitude error. ``scores[(config, season)]``."""
    configs = list(dict.fromkeys(c for c, _ in scores))
    seasons = [s for s in seasons if any(k[1] == s for k in scores)]
    fig, axes = plt.subplots(len(seasons), 3, figsize=(18, 3.0 * len(seasons) + 0.9), squeeze=False,
                             sharex=True, sharey="col", layout="constrained")
    columns = ("Bivariate COR", "Bivariate RMSE (dashed: spread)", "Amplitude error")
    for i, s in enumerate(seasons):
        for c in configs:
            if (c, s) not in scores:
                continue
            d = scores[c, s].sel(lead=slice(1, max_lead))
            axes[i, 0].plot(d.lead, d.cor_ensmean, **_style(c))
            axes[i, 1].plot(d.lead, d.rmse_ensmean, **_style(c))
            axes[i, 1].plot(d.lead, d.spread, **{**_style(c), "linestyle": "--", "linewidth": 1.2})
            axes[i, 2].plot(d.lead, d.amp_error_ensmean, **_style(c))
        for ax, y in zip(axes[i], (cor_ref, np.sqrt(2), 0)):
            ax.axhline(y, color="k", linewidth=0.8)
            _lead_days(ax, max_lead)
    _panel_titles(axes, [list(columns)] * len(seasons), [[s] * 3 for s in seasons])
    for j in range(3):
        axes[-1, j].set_xlabel("Lead day")
    _legend(fig, configs)
    fig.suptitle("RMM skill of the ensemble mean, pooled over the start dates")
    return fig


# --------------------------------------------------------------------------------------------------
# Maps
# --------------------------------------------------------------------------------------------------

def _map_axes(nrows, ncols, width=None, height=None, extent=None):
    """Map panels; ``extent`` (lon0, lon1, lat0, lat1) in degrees, None for the globe; the default panel
    size grows with ``FONT_SIZE`` (4.8 x 2.7 inches at 12 pt)."""
    import cartopy.crs as ccrs

    width = width or 0.4 * FONT_SIZE
    height = height or 0.225 * FONT_SIZE
    fig, axes = plt.subplots(nrows, ncols, figsize=(width * ncols, height * nrows + 1.0), squeeze=False,
                             subplot_kw={"projection": ccrs.PlateCarree(central_longitude=180)}, layout="constrained")
    for ax in axes.flat:
        ax.coastlines(linewidth=0.6, color="0.25")
        if extent is None:
            ax.set_global()
        else:
            ax.set_extent(extent, crs=ccrs.PlateCarree())
    return fig, axes


def _shade(ax, m: xr.DataArray, levels, cmap):
    """Discrete-level shading of a (lat, lon) map. pcolormesh, not contourf: cartopy takes minutes
    to project the contour polygons of a noisy 1-degree field."""
    import cartopy.crs as ccrs
    from matplotlib.colors import BoundaryNorm

    norm = BoundaryNorm(levels, plt.get_cmap(cmap).N, extend="both")
    return ax.pcolormesh(m.lon, m.lat, m, cmap=cmap, norm=norm, shading="nearest", rasterized=True,
                         transform=ccrs.PlateCarree())


@readable
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
        ax.set_title(WEEK_LABELS[w])
    for ax in list(axes.flat)[len(weeks):]:
        ax.set_visible(False)
    if cf is not None:
        fig.colorbar(cf, ax=axes, orientation="horizontal", shrink=0.5, pad=0.03, label=label)
    fig.suptitle(title)
    return fig


@readable
def plot_config_week_maps(maps: dict[str, xr.DataArray], *, week: int, title: str, levels=None, cmap="RdBu_r",
                          label="", extent=None, ncols: int = 3, panel_size=(None, None)) -> plt.Figure:
    """One lead week, one panel per configuration; ``extent`` (lon0, lon1, lat0, lat1) zooms the maps
    (None: the globe), ``panel_size`` (width, height) in inches."""
    keys = list(maps)
    ncols = min(ncols, len(keys))
    nrows = int(np.ceil(len(keys) / ncols))
    fig, axes = _map_axes(nrows, ncols, *panel_size, extent=extent)
    levels = np.linspace(-1, 1, 21) if levels is None else levels
    cf = None
    for k, (ax, c) in enumerate(zip(axes.flat, keys)):
        cf = _shade(ax, maps[c].sel(week=week), levels, cmap)
        ax.set_title(f"({_LETTERS[k]}) {config_label(c)}", loc="left")
    for ax in list(axes.flat)[len(keys):]:
        ax.set_visible(False)
    if cf is not None:
        fig.colorbar(cf, ax=axes, orientation="horizontal", shrink=0.5, pad=0.03, label=label)
    fig.suptitle(title)
    return fig


@readable
def plot_map_grid(panels: dict[str, dict[str, xr.DataArray]], *, title: str, levels=None, cmap="RdBu_r", label="",
                  extent=None, panel_size=(None, None)) -> plt.Figure:
    """Rows x columns of (lat, lon) maps on one colour scale: ``panels[row][column]``; rows labelled on the
    left, columns titled on the top row. ``extent`` (lon0, lon1, lat0, lat1) zooms the maps (None: the globe)."""
    rows = list(panels)
    cols = list(dict.fromkeys(c for r in rows for c in panels[r]))
    fig, axes = _map_axes(len(rows), len(cols), *panel_size, extent=extent)
    levels = np.linspace(-1, 1, 21) if levels is None else levels
    cf, k = None, 0
    for i, r in enumerate(rows):
        for j, c in enumerate(cols):
            ax = axes[i, j]
            if c not in panels[r]:
                ax.set_visible(False)
                continue
            cf = _shade(ax, panels[r][c], levels, cmap)
            ax.set_title(f"({_LETTERS[k]}) {r}: {config_label(c)}", loc="left")
            k += 1
    if cf is not None:
        fig.colorbar(cf, ax=axes, orientation="horizontal", shrink=0.5, pad=0.03, label=label)
    fig.suptitle(title)
    return fig


@readable
def plot_window_maps(panels: dict, *, title: str, levels, cmap="RdBu_r", label="", extent=None,
                     panel_size=(None, None), stipple: dict | None = None, note=None, extend="both") -> plt.Figure:
    """Lead-time skill maps in the layout of ESP-Lab ``1a_*_leadtime_acc_skill_map`` / ``1b_*_rmse_skill_map``.

    ``panels[season][config]``: (window, lat, lon) maps. One section per season (e.g. Dec and May starts),
    stacked and titled "Dec init." / "May init."; rows = lead windows (labelled on the left), columns =
    configurations (titled on each section's top row); one colour scale. ``stipple`` has the structure of
    ``panels`` (bool maps); ``note(season, config, window)`` returns the text of each panel's lower-right box
    (e.g. an area fraction) or None. ``extent`` (lon0, lon1, lat0, lat1) zooms the maps (None: the globe);
    ``panel_size`` (width, height) in inches.
    """
    import cartopy.crs as ccrs
    from cartopy.mpl.ticker import LatitudeFormatter
    from matplotlib.colors import BoundaryNorm

    seasons = list(panels)
    cols = list(dict.fromkeys(c for s in seasons for c in panels[s]))
    windows = list(dict.fromkeys(w for s in seasons for c in panels[s] for w in panels[s][c].window.values))
    lon0, lon1, lat0, lat1 = extent if extent is not None else (0, 360, -90, 90)
    aspect = (lat1 - lat0) / ((lon1 - lon0) % 360 or 360)
    width = panel_size[0] or 0.25 * FONT_SIZE
    height = panel_size[1] or width * aspect + 0.03 * FONT_SIZE
    sec_h = height * len(windows) + 0.08 * FONT_SIZE
    fig = plt.figure(figsize=(width * len(cols) + 0.1 * FONT_SIZE, sec_h * len(seasons) + 0.1 * FONT_SIZE),
                     layout="constrained")
    subs = fig.subfigures(len(seasons) + 1, 1, height_ratios=[sec_h] * len(seasons) + [0.07 * FONT_SIZE],
                          hspace=0.02)
    data_crs = ccrs.PlateCarree()
    proj = ccrs.PlateCarree(central_longitude=180 if extent is None else 0.5 * (lon0 + lon1 + (360 if lon1 < lon0 else 0)))
    norm = BoundaryNorm(levels, plt.get_cmap(cmap).N, extend=extend)
    cf = None
    for sf, s in zip(subs, seasons):
        k = 0   # panel letters restart in each section
        sf.suptitle(SEASON_PANEL.get(s, s), fontsize=fs("title"))
        axes = sf.subplots(len(windows), len(cols), squeeze=False, subplot_kw={"projection": proj})
        for i, win in enumerate(windows):
            for j, c in enumerate(cols):
                ax = axes[i, j]
                if c not in panels[s] or win not in panels[s][c].window.values:
                    ax.set_visible(False)
                    continue
                m = panels[s][c].sel(window=win)
                cf = ax.pcolormesh(m.lon, m.lat, m, cmap=cmap, norm=norm, shading="nearest", rasterized=True,
                                   transform=data_crs)
                ax.coastlines(linewidth=0.6, color="0.25")
                if extent is None:
                    ax.set_global()
                else:
                    ax.set_extent(extent, crs=data_crs)
                if stipple is not None and c in stipple.get(s, {}):
                    st = stipple[s][c].sel(window=win)
                    step = 4 if extent is None else 1
                    st = st.isel(lat=slice(None, None, step), lon=slice(None, None, step))
                    lo, la = np.meshgrid(st.lon.values, st.lat.values)
                    hit = np.asarray(st.fillna(False).values, dtype=bool)
                    ax.plot(lo[hit], la[hit], ".", color="k", markersize=1.2, alpha=0.5, transform=data_crs)
                if extent is None:   # outer-edge lat/lon labels (ESP-Lab layout)
                    if j == 0:
                        ax.set_yticks([-60, 0, 60], crs=data_crs)
                        ax.yaxis.set_major_formatter(LatitudeFormatter())
                    if i == len(windows) - 1:
                        ax.set_xticks([-120, 0, 120], crs=proj)   # map centred on 180E
                        ax.set_xticklabels(["60°E", "180°", "60°W"])
                    ax.tick_params(labelsize=fs("cell"), length=2)
                if i == 0:
                    ax.set_title(config_label(c), fontsize=fs("title"))
                if j == 0:
                    ax.text(-0.2 if extent is None else -0.05, 0.5, str(win), rotation=90, ha="right", va="center",
                            transform=ax.transAxes, fontsize=fs("label"))
                ax.text(0.015, 0.97, f"({_LETTERS[k % 26]})", ha="left", va="top", transform=ax.transAxes,
                        fontsize=fs("cell"), bbox=dict(boxstyle="square,pad=0.15", fc="w", ec="none", alpha=0.8))
                k += 1
                txt = note(s, c, win) if note is not None else None
                if txt:
                    ax.text(0.985, 0.03, txt, ha="right", va="bottom", transform=ax.transAxes, fontsize=fs("cell"),
                            bbox=dict(boxstyle="square,pad=0.15", fc="w", ec="none", alpha=0.8))
                ax.grid(False)
    if cf is not None:
        cax = subs[-1].add_axes([0.25, 0.55, 0.5, 0.3])
        subs[-1].colorbar(cf, cax=cax, orientation="horizontal", label=label)
    fig.suptitle(title)
    return fig


# --------------------------------------------------------------------------------------------------
# Scorecards
# --------------------------------------------------------------------------------------------------

@readable
def plot_scorecard(changes, *, title: str, higher_is_better: bool, fmt="{:+.2f}", vmax=None,
                   units="") -> plt.Figure:
    """Heat tables, rows configurations, columns (variable, week); blue = better than the baseline.

    ``changes``: one DataFrame, or {season: DataFrame} drawn one above the other on one colour scale.
    """
    changes = changes if isinstance(changes, dict) else {"": changes}
    allv = np.concatenate([c.values.astype(float).ravel() for c in changes.values()])
    vmax = vmax or np.nanmax(np.abs(allv)) or 1.0
    cmap = plt.get_cmap("RdBu" if higher_is_better else "RdBu_r")
    first = next(iter(changes.values()))
    cell = FONT_SIZE / 12.0                              # cell size grows with the font
    width = 0.75 * cell * first.shape[1] + 3.0 * cell
    heights = [0.55 * cell * c.shape[0] + 1.2 * cell for c in changes.values()]
    fig, axes = plt.subplots(len(changes), 1, figsize=(width, sum(heights) + 1.0), squeeze=False,
                             gridspec_kw={"height_ratios": heights}, layout="constrained")
    im = None
    for ax, (season, change) in zip(axes[:, 0], changes.items()):
        vals = change.values.astype(float)
        im = ax.imshow(vals, cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto")
        for (i, j), x in np.ndenumerate(vals):
            if np.isfinite(x):
                ax.text(j, i, fmt.format(x), ha="center", va="center", fontsize=fs("cell"),
                        color="w" if abs(x) > 0.6 * vmax else "k")
        ax.set_yticks(range(len(change.index)), [config_label(c) for c in change.index])
        cols = [f"{a}\nW{b}" if isinstance(a, str) else str(a) for a, b in change.columns] \
            if isinstance(change.columns, pd.MultiIndex) else list(map(str, change.columns))
        ax.set_xticks(range(len(cols)), cols, fontsize=fs("cell"))
        ax.grid(False)
        if season:
            ax.set_title(season_title(season), loc="left")
    fig.colorbar(im, ax=axes[:, 0], shrink=0.8, label=units)
    fig.suptitle(title)
    return fig
