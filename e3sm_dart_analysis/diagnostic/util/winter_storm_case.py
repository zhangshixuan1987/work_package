"""Case study of the January 2012 Pacific Northwest storm in the F20TR hindcasts (``hind_lead_eval/07``).

The Dec 2011 starts are verified on a fixed observed window (16-20 Jan 2012): regional daily
temperature and precipitation of every member against ERA5 and GPCP, member event metrics, event
maps, and the large-scale storm environment (sea-level pressure, 500-hPa height, integrated vapour
transport) on chosen days. The fields, boxes, regional means and calendar handling are those of the
heat-wave case (``util.heatwave_case``).

Temperature is used as a daily mean in deg C. ``bias="none"`` keeps the model value; ``"amip"``
subtracts the model's mean bias (AMIP minus observed climatology of the day), so freezing-related
counts are not set by the model's climate offset.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from util.heatwave_case import _date_axis, runs  # noqa: F401  (re-exported helpers)

KELVIN = 273.15


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
def pattern_scores(model: xr.DataArray, obs: xr.DataArray, weights: xr.DataArray) -> tuple[float, float]:
    """Centred pattern correlation and RMSE of two (lat, lon) maps over the weighted points."""
    w = weights.values
    m, o = model.values, obs.values
    ok = (w > 0) & np.isfinite(m) & np.isfinite(o)
    if ok.sum() < 3:
        return np.nan, np.nan
    w, m, o = w[ok] / w[ok].sum(), m[ok], o[ok]
    dm, do = m - (w * m).sum(), o - (w * o).sum()
    r = (w * dm * do).sum() / np.sqrt((w * dm ** 2).sum() * (w * do ** 2).sum())
    return float(r), float(np.sqrt((w * (m - o) ** 2).sum()))


def storm_metrics(t_c: np.ndarray, p: np.ndarray, times, window, obs_peak, *, tolerance_days: int = 3,
                  freeze_c: float = 0.0, near_freeze_c: float = 1.0) -> dict:
    """Metrics of one regional series: daily-mean temperature ``t_c`` (deg C) and precipitation ``p``
    (mm/day) on calendar days ``times``. Cold-wet coincidence: precipitation on days with daily-mean
    temperature at or below ``freeze_c`` / ``near_freeze_c`` (not a diagnosed snowfall)."""
    t = pd.DatetimeIndex(times)
    w0, w1 = pd.Timestamp(window[0]), pd.Timestamp(window[1])
    inw = (t >= w0) & (t <= w1)
    if not inw.any() or np.isnan(p[inw]).all():
        return {}
    tol = (t >= w0 - pd.Timedelta(days=tolerance_days)) & (t <= w1 + pd.Timedelta(days=tolerance_days))
    pk = int(np.argmax(np.where(tol & np.isfinite(p), p, -np.inf)))
    tw, pw = t_c[inw], p[inw]
    return dict(event_precip_mm=float(np.nansum(pw)), max_daily_precip_mm=float(np.nanmax(pw)),
                peak_precip_date_error_days=int((t[pk] - pd.Timestamp(obs_peak)) / pd.Timedelta(days=1)),
                event_mean_t_c=float(np.nanmean(tw)), event_min_t_c=float(np.nanmin(tw)),
                precip_freezing_days_mm=float(np.nansum(np.where(tw <= freeze_c, pw, 0.0))),
                precip_near_freezing_days_mm=float(np.nansum(np.where(tw <= near_freeze_c, pw, 0.0))),
                window_days_covered=int(np.isfinite(pw).sum()))


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def plot_region_series(panels: dict, obs: dict, window, styles: dict, flagged=(), *, ylabels: dict,
                       title: str, spread=(0.1, 0.9), xlim=None):
    """Rows: regions; columns: variables. ``panels[region][var][config]`` is (member, time) with
    members of all starts stacked; ``obs[region][var]`` the observed series (None: model only)."""
    regions = list(panels)
    variables = list(panels[regions[0]])
    fig, axes = plt.subplots(len(regions), len(variables), figsize=(5.4 * len(variables), 2.9 * len(regions)),
                             sharex=True, squeeze=False)
    for i, reg in enumerate(regions):
        for j, v in enumerate(variables):
            ax = axes[i, j]
            ax.axvspan(pd.Timestamp(window[0]), pd.Timestamp(window[1]) + pd.Timedelta(hours=23), color="0.9", zorder=0)
            for c, da in panels[reg][v].items():
                st = dict(styles.get(c, {}))
                if c in flagged:
                    st["linestyle"] = "--"
                tt = pd.DatetimeIndex(da.time.values)
                lo, hi = (da.quantile(q, "member") for q in spread)
                ax.fill_between(tt, lo, hi, color=st.get("color"), alpha=0.12, lw=0, step="mid")
                ax.plot(tt, da.mean("member"), drawstyle="steps-mid", label=c + (" (IC defect)" if c in flagged else ""), **st)
            o = obs[reg].get(v)
            if o is not None:
                ax.plot(pd.DatetimeIndex(o.time.values), o, drawstyle="steps-mid", label="OBS", **styles.get("OBS", {}))
            if v == "TREFHT":
                ax.axhline(0, color="tab:blue", lw=0.8, ls=":")
            ax.set_ylabel(ylabels[v], fontsize=8)
            if j == 0:
                ax.text(0.01, 0.95, reg, transform=ax.transAxes, va="top", fontsize=8, weight="bold")
            if xlim is not None:
                ax.set_xlim(pd.Timestamp(xlim[0]), pd.Timestamp(xlim[1]))
            _date_axis(ax)
    axes[0, 0].legend(fontsize=7, ncol=2, frameon=False, loc="lower left")
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    return fig


def plot_event_maps(obs: xr.DataArray, model: dict, regions: list[dict], flagged=(), *, title: str, cmap: str,
                    vmin: float, vmax: float, units: str, extent):
    """Event map of the observations and each configuration (ensemble and start mean)."""
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    names, fields = ["OBS", *model], [obs, *model.values()]
    ncol = 3
    nrow = int(np.ceil(len(names) / ncol))
    proj = ccrs.LambertConformal(central_longitude=-120, standard_parallels=(40, 50))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.1 * ncol, 3.6 * nrow), subplot_kw=dict(projection=proj), squeeze=False)
    for ax, n, f in zip(axes.flat, names, fields):
        im = ax.pcolormesh(f.lon, f.lat, f, transform=ccrs.PlateCarree(), cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
        ax.coastlines(lw=0.6)
        ax.add_feature(cfeature.STATES, lw=0.3, edgecolor="0.4")
        for r in regions:
            (a0, a1), (o0, o1) = r["lat"], r["lon"]
            ax.plot([o0, o1, o1, o0, o0], [a0, a0, a1, a1, a0], color="k", lw=1, transform=ccrs.PlateCarree())
        ax.set_extent(extent, ccrs.PlateCarree())
        ax.set_title(n + (" (IC defect)" if n in flagged else ""), fontsize=9)
    for ax in list(axes.flat)[len(names):]:
        ax.set_visible(False)
    fig.colorbar(im, ax=axes.ravel().tolist(), orientation="horizontal", shrink=0.5, pad=0.04, label=units)
    fig.suptitle(title, fontsize=10)
    return fig


def plot_environment(days, obs: dict, model: dict, flagged=(), *, title: str, extent, z_levels=None,
                     psl_levels=None, ivt_max: float = 600.0, ivt_step: int = 4):
    """Rows: days; columns: OBS and configurations. Shading: IVT magnitude; vectors: IVT; black
    contours: 500-hPa height; white contours: sea-level pressure. ``obs[day]`` and ``model[config][day]``
    are datasets with ``Z500``, ``PSL``, ``TUQ``, ``TVQ``."""
    import cartopy.crs as ccrs

    names = ["OBS", *model]
    z_levels = np.arange(4800, 6000, 60) if z_levels is None else z_levels
    psl_levels = np.arange(960, 1050, 8) if psl_levels is None else psl_levels
    proj = ccrs.LambertConformal(central_longitude=-135, standard_parallels=(35, 55))
    fig, axes = plt.subplots(len(days), len(names), figsize=(3.4 * len(names), 2.9 * len(days)),
                             subplot_kw=dict(projection=proj), squeeze=False)
    for i, d in enumerate(days):
        for j, n in enumerate(names):
            ax = axes[i, j]
            ds = obs[d] if n == "OBS" else model[n][d]
            ivt = np.hypot(ds.TUQ, ds.TVQ)
            im = ax.pcolormesh(ds.lon, ds.lat, ivt, transform=ccrs.PlateCarree(), cmap="YlGnBu", vmin=0, vmax=ivt_max,
                               shading="auto")
            s = slice(None, None, ivt_step)
            ax.quiver(ds.lon.values[s], ds.lat.values[s], ds.TUQ.values[s, s], ds.TVQ.values[s, s],
                      transform=ccrs.PlateCarree(), scale=12000, width=0.003, color="0.25")
            ax.contour(ds.lon, ds.lat, ds.Z500, levels=z_levels, colors="k", linewidths=0.6, transform=ccrs.PlateCarree())
            cs = ax.contour(ds.lon, ds.lat, ds.PSL, levels=psl_levels, colors="w", linewidths=0.8,
                            transform=ccrs.PlateCarree())
            ax.clabel(cs, fontsize=5, fmt="%d")
            ax.coastlines(lw=0.6, color="0.3")
            ax.set_extent(extent, ccrs.PlateCarree())
            if i == 0:
                ax.set_title(n + (" (IC defect)" if n in flagged else ""), fontsize=8)
            if j == 0:
                ax.text(-0.05, 0.5, str(pd.Timestamp(d).date()), transform=ax.transAxes, rotation=90, va="center",
                        ha="right", fontsize=8)
    fig.colorbar(im, ax=axes.ravel().tolist(), orientation="horizontal", shrink=0.4, pad=0.02,
                 label="integrated vapour transport (kg m-1 s-1)")
    fig.suptitle(title, fontsize=10)
    return fig


__all__ = ["KELVIN", "pattern_scores", "storm_metrics", "plot_region_series", "plot_event_maps", "plot_environment"]
