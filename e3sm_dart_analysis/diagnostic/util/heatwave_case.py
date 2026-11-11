"""Case study of the summer 2012 central U.S. heat wave in the F20TR hindcasts (``hind_lead_eval/06``).

Steps: the simulation coverage of the observed event, the event defined from observations
(calendar-day percentile of daily-mean 2-m temperature), regional daily series of every member
against the observations, heat metrics of every member on the observed event window
(event-centred), ensemble scores (CRPS) on standard S2S lead windows (lead-centred), event
probabilities, and unpaired ensemble contrasts between configurations.

Fields are kept on a map box (e.g. CONUS) as (init, member, lead, lat, lon) and cached once per
configuration and variable; regional means are cos(lat)-weighted over the land points of the
region. Calendar dates of a start's leads are ``init + lead - 1`` (lead 1 = the start date).
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from configs.fcst_init_experiment_config import (case_attrs, case_current, case_names, config_label, init_dates, relabel,
                                                set_case_names)
from util.fcst_init_data import (_reference_series, _convert, model_climatology,
                                 open_member_daily, reference_climatology)
from util.fcst_init_plotting import LEGEND_BOX, fig_size, fs, readable

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
def _load_member_box(args):
    """One (start, member) daily field on the box, loaded: (lead, lat, lon) (a process-pool task)."""
    config, init, member, variable, n_days, box, data_root, land_support_dir, names = args
    set_case_names(names)   # the caller's case names, whatever the worker start method
    da = open_member_daily(config, init, member, variable, max_lead=n_days, data_root=data_root,
                           land_support_dir=land_support_dir, box=box)
    return to_box(da, box).drop_vars("time", errors="ignore").load()


def _parallel_map(func, jobs, workers):
    """``[func(j) for j in jobs]`` over ``workers`` processes; if a worker dies (e.g. a login-node
    memory or process limit), the unfinished jobs are rerun serially in this process."""
    jobs = list(jobs)
    if workers <= 1 or len(jobs) <= 1:
        return [func(j) for j in jobs]
    results, done = [None] * len(jobs), [False] * len(jobs)
    try:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for i, r in enumerate(pool.map(func, jobs)):
                results[i], done[i] = r, True
    except (BrokenProcessPool, OSError) as err:     # a worker killed, or the pool could not start
        todo = [i for i, d in enumerate(done) if not d]
        print(f"[WARN] reader pool failed ({type(err).__name__}); reading the remaining {len(todo)} files serially", flush=True)
        for i in todo:
            results[i] = func(jobs[i])
    return results


def model_box(config: str, season: str, variable: str, members, box: dict, *, n_days: int, data_root,
              land_support_dir, cache: Path | None = None, force: bool = False, workers: int = 4) -> xr.DataArray:
    """(init, member, lead, lat, lon) daily field of ``members`` on ``box``; cached in ``cache``.

    Each (start, member) file is cut to the box before the daily averaging (``open_member_daily(box=)``)
    and the files are read by ``workers`` processes (``workers=1``: serially in this process). On a login
    node keep ``workers`` small (memory is limited per user); on a compute node it can be much larger.
    """
    key = ",".join(members)
    if cache is not None and cache.exists() and not force:
        da = xr.load_dataarray(cache)
        if (da.attrs.get("members") == key and da.attrs.get("box") == str(box) and da.sizes["lead"] == n_days
                and case_current(da.attrs, config)):
            return da
    inits = list(init_dates([season])[season])
    members = list(members)
    names = {config: case_names()[config]}
    jobs = [(config, i, m, variable, n_days, box, data_root, land_support_dir, names) for i in inits for m in members]
    fields = _parallel_map(_load_member_box, jobs, workers)
    units = fields[0].attrs.get("units", "")
    per_init = [xr.concat(fields[k * len(members):(k + 1) * len(members)], pd.Index(members, name="member"),
                          coords="minimal", compat="override") for k in range(len(inits))]
    out = xr.concat(per_init, pd.Index(np.array(inits, dtype="datetime64[ns]"), name="init"), coords="minimal",
                    compat="override").transpose("init", "member", "lead", "lat", "lon")
    out = out.astype("float32").rename(variable).assign_attrs(config=config, season=season, members=key, box=str(box),
                                                              units=units, **case_attrs(config))
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
    threshold) on the event ``window``.

    Primary (exact dates, the observed window): event-mean anomaly, peak anomaly and its date error,
    exceedance days, longest run, accumulated exceedance (K day), mean area fraction. Secondary
    (timing-tolerant, the window widened by +-``tolerance_days``): ``peak_anomaly_tol`` and its date
    error, ``longest_run_tol_days`` (the longest run of hot days the model produces near the event).
    The peak is that of the anomaly series, searched by the same rules for observations and models.
    """
    t = pd.DatetimeIndex(times)
    w0, w1 = pd.Timestamp(window[0]), pd.Timestamp(window[1])
    inw = (t >= w0) & (t <= w1)
    tol = (t >= w0 - pd.Timedelta(days=tolerance_days)) & (t <= w1 + pd.Timedelta(days=tolerance_days))
    if not inw.any() or np.isnan(anom[inw]).all():
        return {}
    a, th = anom[inw], thr_anom[inw]
    exceed = a > th
    longest = max((b - c + 1 for c, b in runs(exceed)), default=0)
    longest_tol = max((b - c + 1 for c, b in runs(tol & (anom > thr_anom))), default=0)

    def peak(sel):
        k = int(np.argmax(np.where(sel & np.isfinite(anom), anom, -np.inf)))
        return float(anom[k]), int((t[k] - pd.Timestamp(obs_peak)) / pd.Timedelta(days=1))

    pk, pk_err = peak(inw)
    pt, pt_err = peak(tol)
    return dict(event_mean_anomaly=float(np.nanmean(a)), peak_anomaly=pk, peak_date_error_days=pk_err,
                exceedance_days=int(exceed.sum()), longest_run_days=int(longest),
                accumulated_exceedance=float(np.nansum(np.clip(a - th, 0, None))),
                area_fraction_mean=float(np.nanmean(area[inw])), window_days_covered=int(np.isfinite(a).sum()),
                peak_anomaly_tol=pt, peak_date_error_tol_days=pt_err, longest_run_tol_days=int(longest_tol))


# ---------------------------------------------------------------------------
# Ensemble scores: CRPS, S2S lead windows, unpaired contrasts
# ---------------------------------------------------------------------------
def crps_ensemble(members, obs: float) -> float:
    """CRPS of the ensemble ``members`` (1-D; NaN dropped) against ``obs``:
    mean|F_m - O| - 1/(2 M^2) sum_m sum_n |F_m - F_n| (the empirical-distribution form)."""
    f = np.asarray(members, dtype="float64")
    f = f[np.isfinite(f)]
    if f.size == 0 or not np.isfinite(obs):
        return np.nan
    return float(np.mean(np.abs(f - obs)) - 0.5 * np.mean(np.abs(f[:, None] - f[None, :])))


def lead_window_scores(fc: xr.DataArray, obs: xr.DataArray, init, windows: dict) -> list[dict]:
    """Scores of one start's regional series ``fc`` (member, time; calendar days) against ``obs`` (time)
    on lead windows ``windows`` ({name: (first lead, last lead)}, lead 1 = the start date): the
    window means of observations and members, ensemble-mean error, member spread (s.d.) and CRPS.
    A window not fully simulated is skipped."""
    i = np.datetime64(pd.Timestamp(init), "D")
    lead = ((fc.time.values.astype("datetime64[D]") - i) / _DAY).astype(int) + 1
    rows = []
    for name, (l0, l1) in windows.items():
        sel = (lead >= l0) & (lead <= l1)
        if sel.sum() < l1 - l0 + 1:
            continue
        days = fc.time.values[sel]
        f = fc.sel(time=days).mean("time").values
        o = float(obs.sel(time=days).mean())
        rows.append(dict(window=name, lead_first=l0, lead_last=l1, date_first=str(days[0])[:10],
                         date_last=str(days[-1])[:10], obs=o, ens_mean=float(np.mean(f)),
                         ens_mean_error=float(np.mean(f) - o), spread=float(np.std(f, ddof=1)),
                         crps=crps_ensemble(f, o)))
    return rows


def bootstrap_contrast(a, b, stat=np.mean, *, n_boot: int = 2000, ci=(0.05, 0.95), seed: int = 0) -> dict:
    """``stat(a) - stat(b)`` of two unpaired ensembles (1-D member values) and its bootstrap interval,
    resampling the members of each ensemble independently (no member pairing is assumed)."""
    a, b = np.asarray(a, dtype="float64"), np.asarray(b, dtype="float64")
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    rng = np.random.default_rng(seed)
    boot = np.array([stat(rng.choice(a, a.size)) - stat(rng.choice(b, b.size)) for _ in range(n_boot)])
    lo, hi = np.quantile(boot, ci)
    return dict(diff=float(stat(a) - stat(b)), lo=float(lo), hi=float(hi))


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def _date_axis(ax):
    ax.xaxis.set_major_locator(mdates.WeekdayLocator(byweekday=mdates.MO, interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))


@readable
def plot_coverage(cov: pd.DataFrame, window, obs_series: xr.DataArray, threshold: xr.DataArray, *, title: str):
    """Observed regional anomaly with its threshold and the window, and each start's simulated span."""
    fig, (ax, bx) = plt.subplots(2, 1, figsize=fig_size(10, 5.2), sharex=True, gridspec_kw=dict(height_ratios=(3, 1.3)))
    t = pd.DatetimeIndex(obs_series.time.values)
    ax.plot(t, obs_series, color="k", lw=2, label="ERA5 anomaly (regional mean)")
    ax.plot(t, threshold, color="firebrick", lw=1.2, ls="--", label="90th percentile minus climatology")
    ax.fill_between(t, threshold, obs_series, where=obs_series > threshold, color="firebrick", alpha=0.3, lw=0)
    for a in (ax, bx):
        a.axvspan(pd.Timestamp(window[0]), pd.Timestamp(window[1]), color="0.85", zorder=0)
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_ylabel("TREFHT anomaly (K)")
    ax.legend(loc="upper left", **LEGEND_BOX)
    ax.set_title(title)
    for k, r in cov.iterrows():
        bx.plot([pd.Timestamp(r.init), pd.Timestamp(r.last_day)], [k, k], lw=6, color="tab:purple", solid_capstyle="butt")
        bx.text(pd.Timestamp(r.last_day), k, f"  lead {r.lead_at_onset}-{r.lead_at_end}", va="center")
    bx.set_yticks(range(len(cov)), [f"start {r.init[5:]}" for _, r in cov.iterrows()])
    bx.set_ylim(-0.7, len(cov) - 0.3)
    _date_axis(bx)
    return fig


@readable
def plot_series(panels: dict, obs: dict, window, styles: dict, flagged=(), *, ylabels: dict, title: str,
                spread=(0.1, 0.9)):
    """Rows: starts; columns: variables. ``panels[var][init][config]`` is (member, time) on calendar days;
    ``obs[var]`` the observed series (time), or ``{init: series}`` when it depends on the start
    (e.g. precipitation accumulated from the start)."""
    variables = list(panels)
    inits = list(panels[variables[0]])
    fig, axes = plt.subplots(len(inits), len(variables), figsize=fig_size(6.2 * len(variables), 2.6 * len(inits)),
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
                ax.plot(t, da.mean("member"), label=config_label(c) + (" (IC defect)" if c in flagged else ""), **st)
            o = obs[v][init] if isinstance(obs[v], dict) else obs[v].sel(time=slice(pd.Timestamp(init), None))
            ax.plot(pd.DatetimeIndex(o.time.values), o, **styles.get("OBS", dict(color="k")), label="OBS")
            ax.set_ylabel(ylabels[v])
            if j == 0:
                ax.text(0.01, 0.95, f"start {str(init)[:10]}", transform=ax.transAxes, va="top")
            _date_axis(ax)
    axes[0, 0].legend(ncol=3, **LEGEND_BOX, loc="lower left")
    fig.suptitle(title)
    return fig


def _metric_title(m: str, units: dict | None) -> str:
    u = (units or {}).get(m, "")
    return m.replace("_", " ") + (f" ({u})" if u else "")


@readable
def plot_by_start(table: pd.DataFrame, obs_metrics: dict, metrics: list[str], configs: list[str], styles: dict, *,
                  title: str, units: dict | None = None, spread=(0.1, 0.9)):
    """Start-resolved metrics: for each start (x) and configuration, the member median (point) and the
    ``spread`` member quantiles (bar); the observed value as a black line."""
    inits = sorted(table.init.unique())
    ncol = min(3, len(metrics))
    nrow = int(np.ceil(len(metrics) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=fig_size(4.3 * ncol, 3.0 * nrow), squeeze=False)
    for ax, m in zip(axes.flat, metrics):
        for k, c in enumerate(configs):
            col = styles.get(c, {}).get("color")
            off = (k - (len(configs) - 1) / 2) * 0.15
            for x, init in enumerate(inits):
                v = table[(table.config == c) & (table.init == init)][m].dropna()
                if v.empty:
                    continue
                lo, hi = v.quantile(spread[0]), v.quantile(spread[1])
                ax.plot([x + off] * 2, [lo, hi], color=col, lw=1.2)
                ax.plot(x + off, v.median(), "o", color=col, ms=4, label=config_label(c) if x == 0 else None)
        if m in obs_metrics:
            ax.axhline(obs_metrics[m], color="k", lw=1.4, label="OBS")
        ax.set_xticks(range(len(inits)), [f"start {str(i)[5:10]}" for i in inits])
        ax.set_title(_metric_title(m, units))
    for ax in list(axes.flat)[len(metrics):]:
        ax.set_visible(False)
    axes[0, 0].legend(**LEGEND_BOX)
    fig.suptitle(title)
    return fig


@readable
def plot_lead_scores(scores: pd.DataFrame, configs: list[str], styles: dict, columns: dict, *, title: str,
                     event_leads: dict | None = None):
    """Rows: scores (``columns``: {column: y label}); columns: starts; x: lead windows. Windows that
    overlap the observed event's leads of a start (``event_leads[init] = (first, last)``) are shaded."""
    inits = sorted(scores.init.unique())
    wins = list(dict.fromkeys(scores.sort_values("lead_first").window))
    fig, axes = plt.subplots(len(columns), len(inits), figsize=fig_size(3.4 * len(inits), 2.5 * len(columns)),
                             sharex=True, sharey="row", squeeze=False)
    for j, init in enumerate(inits):
        s = scores[scores.init == init]
        for i, (col, ylab) in enumerate(columns.items()):
            ax = axes[i, j]
            if event_leads and init in event_leads:
                e0, e1 = event_leads[init]
                for x, w in enumerate(wins):
                    r = s[s.window == w]
                    if not r.empty and r.lead_first.iloc[0] <= e1 and r.lead_last.iloc[0] >= e0:
                        ax.axvspan(x - 0.5, x + 0.5, color="0.9", zorder=0)
            for c in configs:
                r = s[s.config == c].set_index("window").reindex(wins)
                ax.plot(range(len(wins)), r[col], marker="o", ms=4, label=config_label(c), **styles.get(c, {}))
            if col == "ens_mean_error":
                ax.axhline(0, color="0.5", lw=0.8)
            if j == 0:
                ax.set_ylabel(ylab)
            if i == 0:
                ax.set_title(f"start {str(init)[:10]}")
            ax.set_xticks(range(len(wins)), wins, rotation=30)
    axes[0, 0].legend(**LEGEND_BOX)
    fig.suptitle(title)
    return fig


@readable
def plot_event_probability(prob: pd.DataFrame, configs: list[str], styles: dict, *, title: str):
    """Grouped bars: the member fraction with the event (``p_event``) per start and configuration;
    the Brier score of each bar is written above it."""
    inits = sorted(prob.init.unique())
    fig, ax = plt.subplots(figsize=fig_size(2.2 * len(inits) + 2, 3.4))
    width = 0.8 / len(configs)
    for k, c in enumerate(configs):
        r = prob[prob.config == c].set_index("init").reindex(inits)
        x = np.arange(len(inits)) + (k - (len(configs) - 1) / 2) * width
        ax.bar(x, r.p_event, width * 0.9, color=styles.get(c, {}).get("color"), label=config_label(c))
        for xi, (p, bs) in zip(x, zip(r.p_event, r.brier)):
            if np.isfinite(p):
                ax.text(xi, p + 0.02, f"{bs:.2f}", ha="center", fontsize=fs("cell"), rotation=90)
    ax.set_xticks(range(len(inits)), [f"start {str(i)[:10]}" for i in inits])
    ax.set_ylim(0, 1.15)
    ax.set_ylabel("member fraction with the event")
    ax.legend(**LEGEND_BOX, ncol=len(configs), loc="upper left")
    ax.set_title(title)
    return fig


@readable
def plot_contrasts(table: pd.DataFrame, metrics: list[str], *, title: str, units: dict | None = None):
    """One panel per metric: each contrast (y) and start (marker), the difference and its bootstrap
    interval; zero marked."""
    contrasts = list(dict.fromkeys(table.contrast))
    inits = sorted(table.init.unique())
    markers = dict(zip(inits, "o^sDvP"))
    fig, axes = plt.subplots(1, len(metrics), figsize=fig_size(4.0 * len(metrics), 0.55 * len(contrasts) + 1.6),
                             sharey=True, squeeze=False)
    for ax, m in zip(axes[0], metrics):
        for y, con in enumerate(contrasts):
            for k, init in enumerate(inits):
                r = table[(table.contrast == con) & (table.init == init) & (table.metric == m)]
                if r.empty:
                    continue
                yy = y + (k - (len(inits) - 1) / 2) * 0.15
                ax.plot([r.lo.iloc[0], r.hi.iloc[0]], [yy, yy], color="0.4", lw=1)
                ax.plot(r["diff"].iloc[0], yy, markers[init], color="tab:blue", ms=4,
                        label=f"start {str(init)[:10]}" if y == 0 else None)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_title(_metric_title(m, units))
        ax.set_yticks(range(len(contrasts)), [relabel(str(c)) for c in contrasts])
        ax.invert_yaxis()
    axes[0, 0].legend(**LEGEND_BOX)
    fig.suptitle(title)
    return fig


@readable
def plot_metrics(table: pd.DataFrame, obs_metrics: dict, metrics: list[str], configs: list[str], styles: dict,
                 flagged=(), *, title: str, units: dict | None = None):
    """Member values (points, one marker per start) and the observed value (black line) of each metric."""
    inits = sorted(table.init.unique())
    markers = dict(zip(inits, "o^sDvP"))
    ncol = min(3, len(metrics))
    nrow = int(np.ceil(len(metrics) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=fig_size(4.3 * ncol, 3.2 * nrow), squeeze=False)
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
        ax.set_xticks(range(len(configs)), [config_label(c) + ("*" if c in flagged else "") for c in configs], rotation=30)
        ax.set_title(_metric_title(m, units))
    for ax in list(axes.flat)[len(metrics):]:
        ax.set_visible(False)
    axes[0, 0].legend(**LEGEND_BOX, title="start")
    fig.suptitle(title + ("  (* IC defect: missing radiative tracers)" if flagged else ""))
    return fig


@readable
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
    fig, axes = plt.subplots(nrow, ncol, figsize=fig_size(4.6 * ncol, 3.2 * nrow), subplot_kw=dict(projection=proj), squeeze=False)
    for ax, n, f in zip(axes.flat, names, fields):
        im = ax.pcolormesh(f.lon, f.lat, f, transform=ccrs.PlateCarree(), cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                           shading="auto")
        ax.coastlines(lw=0.5)
        ax.add_feature(cfeature.STATES, lw=0.3, edgecolor="0.5")
        (a0, a1), (o0, o1) = region["lat"], region["lon"]
        ax.plot([o0, o1, o1, o0, o0], [a0, a0, a1, a1, a0], color="k", lw=1, transform=ccrs.PlateCarree())
        ax.set_extent([float(f.lon.min()), float(f.lon.max()), float(f.lat.min()), float(f.lat.max())], ccrs.PlateCarree())
        ax.set_title(config_label(n) + (" (IC defect)" if n in flagged else ""))
    for ax in list(axes.flat)[len(names):]:
        ax.set_visible(False)
    fig.colorbar(im, ax=axes.ravel().tolist(), orientation="horizontal", shrink=0.5, pad=0.04, label=f"anomaly ({units})")
    fig.suptitle(title)
    return fig


@readable
def plot_maps_by_start(obs_map: xr.DataArray, maps: dict, region: dict, *, title: str, vmax: float = 6.0,
                       units: str = "K"):
    """Rows: starts; columns: OBS (first row only) and each configuration's ensemble-mean map
    (``maps[init][config]``)."""
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    inits = list(maps)
    configs = list(maps[inits[0]])
    ncol = len(configs) + 1
    proj = ccrs.LambertConformal(central_longitude=-97, standard_parallels=(33, 45))
    fig, axes = plt.subplots(len(inits), ncol, figsize=fig_size(3.6 * ncol, 2.5 * len(inits)),
                             subplot_kw=dict(projection=proj), squeeze=False)
    im = None
    for i, init in enumerate(inits):
        for j, (n, f) in enumerate([("OBS", obs_map if i == 0 else None), *maps[init].items()]):
            ax = axes[i, j]
            if f is None:
                ax.set_visible(False)
                continue
            im = ax.pcolormesh(f.lon, f.lat, f, transform=ccrs.PlateCarree(), cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                               shading="auto")
            ax.coastlines(lw=0.5)
            ax.add_feature(cfeature.STATES, lw=0.3, edgecolor="0.5")
            (a0, a1), (o0, o1) = region["lat"], region["lon"]
            ax.plot([o0, o1, o1, o0, o0], [a0, a0, a1, a1, a0], color="k", lw=1, transform=ccrs.PlateCarree())
            ax.set_extent([float(f.lon.min()), float(f.lon.max()), float(f.lat.min()), float(f.lat.max())],
                          ccrs.PlateCarree())
            ax.set_title(n if n == "OBS" else f"{config_label(n)}, start {str(init)[5:10]}")
    fig.colorbar(im, ax=axes.ravel().tolist(), orientation="horizontal", shrink=0.5, pad=0.03, label=f"anomaly ({units})")
    fig.suptitle(title)
    return fig


__all__ = ["to_box", "region_weights", "region_mean", "calendar", "coverage_table", "model_box", "observed_box",
           "climatologies", "percentile_threshold", "runs", "observed_events", "heat_metrics", "crps_ensemble",
           "lead_window_scores", "bootstrap_contrast", "plot_coverage", "plot_series", "plot_metrics", "plot_by_start",
           "plot_lead_scores", "plot_event_probability", "plot_contrasts", "plot_maps", "plot_maps_by_start"]
