"""Case study of the summer 2012 central U.S. heat wave in the F20TR hindcasts (``hind_lead_eval/06``).

Steps: the simulation coverage of the observed event, the event defined from observations
(calendar-day percentile of daily-mean 2-m temperature), regional daily series of every member
against the observations, and heat metrics of every member on the observed event window.

Fields are kept on a map box (e.g. CONUS) as (init, member, lead, lat, lon) and cached once per
configuration and variable; regional means are cos(lat)-weighted over the land points of the
region. Calendar dates of a start's leads are ``init + lead - 1`` (lead 1 = the start date).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from util.fcst_init_data import (_reference_series, _convert, model_climatology, open_campaign,
                                 reference_climatology)

_DAY = np.timedelta64(1, "D")


# ---------------------------------------------------------------------------
# Boxes, masks and regional means
# ---------------------------------------------------------------------------
def to_box(da: xr.DataArray, box: dict) -> xr.DataArray:
    """``da`` on the lat/lon box ``box`` (``lat``, ``lon`` ranges in degrees north / east 0-360)."""
    return da.sel(lat=slice(*box["lat"]), lon=slice(*box["lon"]))


def region_weights(mask: xr.DataArray) -> xr.DataArray:
    """cos(lat) weights on the region's points, zero elsewhere."""
    return (np.cos(np.deg2rad(mask.lat)) * mask).astype("float64")


def region_mean(da: xr.DataArray, weights: xr.DataArray) -> xr.DataArray:
    """Weighted mean over lat/lon with ``weights`` (``region_weights``) on ``da``'s grid."""
    w = weights.reindex_like(da.isel({d: 0 for d in da.dims if d not in ("lat", "lon")}, drop=True))
    return (da * w).sum(("lat", "lon")) / w.sum()


def calendar(init, lead) -> np.ndarray:
    """Calendar day of each lead (lead 1 = the start date)."""
    return np.datetime64(pd.Timestamp(init), "D") + (np.asarray(lead) - 1) * _DAY


def coverage_table(inits, n_days: int, window: tuple[str, str]) -> pd.DataFrame:
    """Start, last simulated day and the lead days of the event window for each start."""
    w0, w1 = (np.datetime64(pd.Timestamp(x), "D") for x in window)
    rows = []
    for i in inits:
        i = np.datetime64(pd.Timestamp(i), "D")
        last = i + (n_days - 1) * _DAY
        rows.append(dict(init=str(i), last_day=str(last), lead_at_onset=int((w0 - i) / _DAY) + 1,
                         lead_at_end=int((w1 - i) / _DAY) + 1,
                         window_covered="full" if last >= w1 else ("onset only" if last >= w0 else "no")))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Model and observed fields (cached)
# ---------------------------------------------------------------------------
def model_box(config: str, season: str, variable: str, members, box: dict, *, n_days: int, data_root,
              land_support_dir, cache: Path | None = None, force: bool = False) -> xr.DataArray:
    """(init, member, lead, lat, lon) daily field of ``members`` on ``box``; cached in ``cache``."""
    key = ",".join(members)
    if cache is not None and cache.exists() and not force:
        da = xr.load_dataarray(cache)
        if da.attrs.get("members") == key and da.attrs.get("box") == str(box) and da.sizes["lead"] == n_days:
            return da
    field = open_campaign(config, season, variable, members=list(members), max_lead=n_days, data_root=data_root,
                          land_support_dir=land_support_dir)
    out = to_box(field, box).drop_vars("time").load().transpose("init", "member", "lead", "lat", "lon")
    out = out.astype("float32").rename(variable).assign_attrs(config=config, season=season, members=key, box=str(box),
                                                              units=field.attrs.get("units", ""))
    if cache is not None:
        out.to_netcdf(cache)
    return out


def observed_box(variable: str, start: str, end: str, box: dict, *, observation_root) -> xr.DataArray:
    """(time, lat, lon) daily observations (the variable's reference product) on ``box``."""
    days = pd.date_range(start, end, freq="D")
    series = _reference_series(variable, days[0], days[-1], observation_root)
    series = series.reindex(time=days.values)
    out = _convert(to_box(series, box).load(), variable).rename(variable)
    return out


def climatologies(variable: str, start: str, end: str, years, box: dict, *, observation_root, amip_root,
                  amip_cases, land_support_dir) -> xr.Dataset:
    """Observed and AMIP climatologies (``clim_obs``, ``clim_amip``) on the calendar days start..end."""
    days = pd.date_range(start, end, freq="D")
    inits = [np.datetime64(days[0], "D")]
    kw = dict(max_lead=len(days))
    obs = reference_climatology(variable, inits, years, observation_root=observation_root, **kw)
    mod = model_climatology(variable, inits, years, amip_root=amip_root, cases=amip_cases,
                            land_support_dir=land_support_dir, **kw)
    out = xr.Dataset({n: to_box(d.isel(init=0, drop=True), box).drop_vars("time", errors="ignore")
                      .rename(lead="time").assign_coords(time=days.values).load()
                      for n, d in (("clim_obs", obs), ("clim_amip", mod))})
    out.attrs.update(clim_obs=obs.attrs.get("climatology_source", ""), clim_amip=mod.attrs.get("climatology_source", ""))
    return out


def percentile_threshold(variable: str, start: str, end: str, years, box: dict, weights: xr.DataArray, *,
                         q: float = 0.9, half_window: int = 7, observation_root) -> xr.Dataset:
    """Calendar-day ``q`` quantile over ``years`` and +-``half_window`` days, per grid point and of the
    regional mean, on the calendar days start..end of the event year (``thr_grid``, ``thr_region``)."""
    days = pd.date_range(start, end, freq="D")
    ref_year = days[0].year
    grid, region = [], []
    for y in years:
        a = days[0] - pd.Timedelta(days=half_window) + pd.DateOffset(years=y - ref_year)
        b = days[-1] + pd.Timedelta(days=half_window) + pd.DateOffset(years=y - ref_year)
        s = observed_box(variable, str(a.date()), str(b.date()), box, observation_root=observation_root)
        grid.append(s)
        region.append(region_mean(s, weights))
    def stack(series_list):
        windows = []
        for s in series_list:
            v = s.values
            n = len(days)
            windows.append(np.stack([v[k:k + n] for k in range(2 * half_window + 1)]))  # (window, day, ...)
        return np.concatenate(windows)                                                 # (sample, day, ...)
    g = np.nanquantile(stack(grid), q, axis=0)
    r = np.nanquantile(stack(region), q, axis=0)
    lat, lon = grid[0].lat, grid[0].lon
    return xr.Dataset(dict(thr_grid=(("time", "lat", "lon"), g), thr_region=(("time",), r)),
                      coords=dict(time=days.values, lat=lat, lon=lon),
                      attrs=dict(quantile=q, half_window=half_window, years=f"{min(years)}-{max(years)}"))


# ---------------------------------------------------------------------------
# Events and metrics
# ---------------------------------------------------------------------------
def runs(flag: np.ndarray) -> list[tuple[int, int]]:
    """(first, last) index of each run of True values."""
    out, k, n = [], 0, len(flag)
    while k < n:
        if flag[k]:
            j = k
            while j + 1 < n and flag[j + 1]:
                j += 1
            out.append((k, j))
            k = j + 1
        else:
            k += 1
    return out


def observed_events(series: xr.DataArray, threshold: xr.DataArray, anomaly: xr.DataArray,
                    min_days: int = 3) -> pd.DataFrame:
    """Runs of at least ``min_days`` days with ``series`` above ``threshold`` (regional means)."""
    flag = (series > threshold).values
    t = pd.DatetimeIndex(series.time.values)
    rows = []
    for a, b in runs(flag):
        if b - a + 1 < min_days:
            continue
        seg = anomaly.isel(time=slice(a, b + 1))
        k = int(np.argmax(seg.values))
        rows.append(dict(onset=t[a].date(), end=t[b].date(), days=b - a + 1, peak_date=t[a + k].date(),
                         peak_anomaly=float(seg.max()), mean_anomaly=float(seg.mean()),
                         accumulated_exceedance=float((series - threshold).isel(time=slice(a, b + 1)).sum())))
    return pd.DataFrame(rows)


def heat_metrics(anom: np.ndarray, times, thr_anom: np.ndarray, area: np.ndarray, window, obs_peak,
                 tolerance_days: int = 7) -> dict:
    """Metrics of one regional anomaly series ``anom`` (same calendar days ``times`` as ``thr_anom``,
    the observed threshold minus climatology, and ``area``, the fraction of region points above their
    threshold) on the event ``window``; peak timing searched within +-``tolerance_days`` of it."""
    t = pd.DatetimeIndex(times)
    w0, w1 = pd.Timestamp(window[0]), pd.Timestamp(window[1])
    inw = (t >= w0) & (t <= w1)
    tol = (t >= w0 - pd.Timedelta(days=tolerance_days)) & (t <= w1 + pd.Timedelta(days=tolerance_days))
    if not inw.any() or np.isnan(anom[inw]).all():
        return {}
    a, th = anom[inw], thr_anom[inw]
    exceed = a > th
    longest = max((b - c + 1 for c, b in runs(exceed)), default=0)
    at = np.where(tol & np.isfinite(anom), anom, -np.inf)
    k = int(np.argmax(at))
    return dict(event_mean_anomaly=float(np.nanmean(a)), peak_anomaly=float(anom[k]),
                peak_date_error_days=int((t[k] - pd.Timestamp(obs_peak)) / pd.Timedelta(days=1)),
                exceedance_days=int(exceed.sum()), longest_run_days=int(longest),
                accumulated_exceedance=float(np.nansum(np.clip(a - th, 0, None))),
                area_fraction_mean=float(np.nanmean(area[inw])), window_days_covered=int(np.isfinite(a).sum()))


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def _date_axis(ax):
    ax.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=mdates.MO, interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))


def plot_coverage(cov: pd.DataFrame, window, obs_series: xr.DataArray, threshold: xr.DataArray, *, title: str):
    """Observed regional anomaly with its threshold and the window, and each start's simulated span."""
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(10, 5.2), sharex=True, gridspec_kw=dict(height_ratios=(3, 1.3)))
    t = pd.DatetimeIndex(obs_series.time.values)
    ax.plot(t, obs_series, color="k", lw=2, label="ERA5 anomaly (regional mean)")
    ax.plot(t, threshold, color="firebrick", lw=1.2, ls="--", label="90th percentile minus climatology")
    ax.fill_between(t, threshold, obs_series, where=obs_series > threshold, color="firebrick", alpha=0.3, lw=0)
    for a in (ax, bx):
        a.axvspan(pd.Timestamp(window[0]), pd.Timestamp(window[1]), color="0.85", zorder=0)
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_ylabel("TREFHT anomaly (K)")
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    ax.set_title(title, fontsize=10)
    for k, r in cov.iterrows():
        bx.plot([pd.Timestamp(r.init), pd.Timestamp(r.last_day)], [k, k], lw=6, color="tab:purple", solid_capstyle="butt")
        bx.text(pd.Timestamp(r.last_day), k, f"  lead {r.lead_at_onset}-{r.lead_at_end}", va="center", fontsize=8)
    bx.set_yticks(range(len(cov)), [f"start {r.init[5:]}" for _, r in cov.iterrows()], fontsize=8)
    bx.set_ylim(-0.7, len(cov) - 0.3)
    _date_axis(bx)
    fig.tight_layout()
    return fig


def plot_series(panels: dict, obs: dict, window, styles: dict, flagged=(), *, ylabels: dict, title: str,
                spread=(0.1, 0.9)):
    """Rows: starts; columns: variables. ``panels[var][init][config]`` is (member, time) on calendar days;
    ``obs[var]`` the observed series (time), or ``{init: series}`` when it depends on the start
    (e.g. precipitation accumulated from the start)."""
    variables = list(panels)
    inits = list(panels[variables[0]])
    fig, axes = plt.subplots(len(inits), len(variables), figsize=(6.2 * len(variables), 2.6 * len(inits)),
                             sharex=True, squeeze=False)
    for i, init in enumerate(inits):
        for j, v in enumerate(variables):
            ax = axes[i, j]
            ax.axvspan(pd.Timestamp(window[0]), pd.Timestamp(window[1]), color="0.9", zorder=0)
            for c, da in panels[v][init].items():
                st = dict(styles.get(c, {}))
                if c in flagged:
                    st["linestyle"] = "--"
                t = pd.DatetimeIndex(da.time.values)
                lo, hi = (da.quantile(q, "member") for q in spread)
                ax.fill_between(t, lo, hi, color=st.get("color"), alpha=0.12, lw=0)
                ax.plot(t, da.mean("member"), label=c + (" (IC defect)" if c in flagged else ""), **st)
            o = obs[v][init] if isinstance(obs[v], dict) else obs[v].sel(time=slice(pd.Timestamp(init), None))
            ax.plot(pd.DatetimeIndex(o.time.values), o, **styles.get("OBS", dict(color="k")), label="OBS")
            ax.set_ylabel(ylabels[v], fontsize=8)
            if j == 0:
                ax.text(0.01, 0.95, f"start {str(init)[:10]}", transform=ax.transAxes, va="top", fontsize=8)
            _date_axis(ax)
    axes[0, 0].legend(fontsize=7, ncol=3, frameon=False, loc="lower left")
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    return fig


def plot_metrics(table: pd.DataFrame, obs_metrics: dict, metrics: list[str], configs: list[str], styles: dict,
                 flagged=(), *, title: str):
    """Member values (points, one marker per start) and the observed value (black line) of each metric."""
    inits = sorted(table.init.unique())
    markers = dict(zip(inits, "o^sDvP"))
    ncol = min(3, len(metrics))
    nrow = int(np.ceil(len(metrics) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.3 * ncol, 3.2 * nrow), squeeze=False)
    for ax, m in zip(axes.flat, metrics):
        for x, c in enumerate(configs):
            col = styles.get(c, {}).get("color")
            for k, init in enumerate(inits):
                v = table[(table.config == c) & (table.init == init)][m].values
                jitter = (k - (len(inits) - 1) / 2) * 0.12
                ax.scatter(np.full(v.size, x + jitter), v, s=10, color=col, marker=markers[init], alpha=0.55,
                           label=str(init)[:10] if x == 0 else None)
            ens = table[table.config == c][m]
            ax.hlines(ens.median(), x - 0.32, x + 0.32, color=col, lw=2.5)
        if m in obs_metrics:
            ax.axhline(obs_metrics[m], color="k", lw=1.6)
        ax.set_xticks(range(len(configs)), [c + ("*" if c in flagged else "") for c in configs], rotation=30, fontsize=8)
        ax.set_title(m.replace("_", " "), fontsize=9)
    for ax in list(axes.flat)[len(metrics):]:
        ax.set_visible(False)
    axes[0, 0].legend(fontsize=7, frameon=False, title="start", title_fontsize=7)
    fig.suptitle(title + ("  (* IC defect: missing radiative tracers)" if flagged else ""), fontsize=10)
    fig.tight_layout()
    return fig


def plot_maps(obs_map: xr.DataArray, model_maps: dict, region: dict, flagged=(), *, title: str, vmax: float = 6.0,
              units: str = "K"):
    """Event-window mean anomaly: observations and each configuration's ensemble mean."""
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    names = ["OBS", *model_maps]
    fields = [obs_map, *model_maps.values()]
    ncol = 3
    nrow = int(np.ceil(len(names) / ncol))
    proj = ccrs.LambertConformal(central_longitude=-97, standard_parallels=(33, 45))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 3.2 * nrow), subplot_kw=dict(projection=proj), squeeze=False)
    for ax, n, f in zip(axes.flat, names, fields):
        im = ax.pcolormesh(f.lon, f.lat, f, transform=ccrs.PlateCarree(), cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                           shading="auto")
        ax.coastlines(lw=0.5)
        ax.add_feature(cfeature.STATES, lw=0.3, edgecolor="0.5")
        (a0, a1), (o0, o1) = region["lat"], region["lon"]
        ax.plot([o0, o1, o1, o0, o0], [a0, a0, a1, a1, a0], color="k", lw=1, transform=ccrs.PlateCarree())
        ax.set_extent([float(f.lon.min()), float(f.lon.max()), float(f.lat.min()), float(f.lat.max())], ccrs.PlateCarree())
        ax.set_title(n + (" (IC defect)" if n in flagged else ""), fontsize=9)
    for ax in list(axes.flat)[len(names):]:
        ax.set_visible(False)
    fig.colorbar(im, ax=axes.ravel().tolist(), orientation="horizontal", shrink=0.5, pad=0.04, label=f"anomaly ({units})")
    fig.suptitle(title, fontsize=10)
    return fig


__all__ = ["to_box", "region_weights", "region_mean", "calendar", "coverage_table", "model_box", "observed_box",
           "climatologies", "percentile_threshold", "runs", "observed_events", "heat_metrics", "plot_coverage",
           "plot_series", "plot_metrics", "plot_maps"]
