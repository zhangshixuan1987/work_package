"""Subseasonal skill of tropical winds in the AMIP hindcasts.

Daily U850/U200 (model vs ERA5) over ``lat_range`` are reduced to anomalies:

* ``harmonics`` (default): departure from each record's mean plus first
  ``n_harmonics`` annual harmonics over the analysis period (as for the MJO
  Hovmöller); valid from lead day 0;
* ``bandpass_20_100``: 20-100 day Lanczos band-pass. The filter needs
  ``numtaps // 2`` days on each side, and the model record starts at the
  initial date, so only leads >= ``numtaps // 2`` (100 for 201 taps) are
  valid; the reference is filtered over its whole record;
* ``rm60``: departure from a centred running mean (``min_periods`` 1/3 of
  the window, so the record ends are kept).
 Per day: cos-lat weighted spatial anomaly correlation
(ACC) and RMSE; a ``smooth_days`` running mean for display; window means of
ACC with a moving-block bootstrap CI (within-window sampling uncertainty of
the single realization). Windows start at lead 0, where the initial state
matters most.

Time scales (``skill_timescales``): day-0 ACC and the e-folding lead, the
first (interpolated) lead at which the ACC falls below ACC(0)/e.
``crossover_lead``: first lead from which the mean ACC of one atmosphere
initial-state group stays above another's (e.g. DART above ERA5), a measure of
initial shock: a state closer to the reference at day 0 that loses skill faster
is overtaken. ``plot_wind_skill_groups`` draws the ACC grouped by atmosphere
initial state and the DART - ERA5 difference with its factorial 95% interval
(``util.hind_factorial``).

The verifying analysis (ERA5) is also the ERA5-initialized runs' initial
state, so their day-0 ACC is favoured; read the first days with that in mind.
"""
from __future__ import annotations

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from configs.hind_experiment_config import STYLE_MAP
from util.hind_factorial import GROUP_COLORS, additive_fit, groups_by, plot_group_curves
from util.mjo_hovmoller_coupling import remove_seasonal_cycle


@dataclass(frozen=True)
class WindSkillConfig:
    lat_range: tuple[float, float] = (-20.0, 20.0)
    lon_range: tuple[float, float] | None = None
    lead_days_max: int = 120
    smooth_days: int = 7
    anom_method: str = "harmonics"   # "harmonics" | "bandpass_20_100" | "rm60"
    n_harmonics: int = 3
    bp_low: float = 20.0
    bp_high: float = 100.0
    bp_numtaps: int = 201
    rm_window: int = 60
    lead_windows: tuple[tuple[int, int], ...] = ((0, 4), (5, 14), (15, 28), (29, 42))
    block_len: int = 7
    n_boot: int = 2000
    seed: int = 12345


def lanczos_bandpass_weights(low: float, high: float, numtaps: int) -> np.ndarray:
    if numtaps % 2 == 0:
        raise ValueError("numtaps must be odd")
    n = np.arange(numtaps) - (numtaps - 1) / 2.0
    h = 2 / low * np.sinc(2 / low * n) - 2 / high * np.sinc(2 / high * n)
    w = h * np.sinc(n / ((numtaps - 1) / 2.0))
    return w - w.mean()


def _filter(da: xr.DataArray, cfg: WindSkillConfig) -> xr.DataArray:
    if cfg.anom_method == "harmonics":
        return remove_seasonal_cycle(da, "harmonics", cfg.n_harmonics)
    if cfg.anom_method == "rm60":
        return da - da.rolling(time=cfg.rm_window, center=True, min_periods=max(5, cfg.rm_window // 3)).mean()
    if cfg.anom_method != "bandpass_20_100":
        raise ValueError(f"Unknown anom_method {cfg.anom_method!r}")
    w = lanczos_bandpass_weights(cfg.bp_low, cfg.bp_high, cfg.bp_numtaps)
    half = len(w) // 2

    def conv(x):
        y = np.convolve(x, w, mode="same")
        y[:half] = np.nan
        y[-half:] = np.nan
        return y
    return xr.apply_ufunc(conv, da, input_core_dims=[["time"]], output_core_dims=[["time"]], vectorize=True)


def _fill_daily(da: xr.DataArray) -> xr.DataArray:
    days = np.arange(da.time.values[0].astype("datetime64[D]"), da.time.values[-1].astype("datetime64[D]") + 1)
    return da.reindex(time=days.astype("datetime64[ns]")).interpolate_na("time")


def domain_subset(da: xr.DataArray, cfg: WindSkillConfig) -> xr.DataArray:
    """Lazy ``lat_range``/``lon_range`` subset; load it (e.g. with ``compute_all``) before ``daily_skill``."""
    da = da.sel(lat=slice(*cfg.lat_range))
    if cfg.lon_range is not None:
        da = da.sel(lon=slice(*cfg.lon_range))
    return da


def daily_skill(model: xr.DataArray, ref: xr.DataArray, init: str, cfg: WindSkillConfig) -> xr.Dataset:
    """Daily ACC and RMSE of one variable vs lead day (day of ``init`` = lead 0)."""
    m = _filter(_fill_daily(domain_subset(model, cfg).compute()), cfg)
    r = _filter(_fill_daily(domain_subset(ref, cfg).compute()), cfg)
    t0 = np.datetime64(init, "D")
    sl = slice(t0, t0 + np.timedelta64(cfg.lead_days_max, "D"))
    m, r = xr.align(m.sel(time=sl), r.sel(time=sl), join="inner")
    w = np.cos(np.deg2rad(m.lat))
    w = w / w.mean()
    wm = lambda x: (x * w).mean(("lat", "lon"))
    a, b = m - wm(m), r - wm(r)
    acc = wm(a * b) / np.sqrt(wm(a * a) * wm(b * b))
    rmse = np.sqrt(wm((m - r) ** 2))
    lead = ((m.time.values.astype("datetime64[D]") - t0) / np.timedelta64(1, "D")).astype(int)
    out = xr.Dataset({"ACC": acc, "RMSE": rmse}).assign_coords(lead=("time", lead)).swap_dims(time="lead")
    sm = lambda x: x.rolling(lead=cfg.smooth_days, center=True, min_periods=max(3, cfg.smooth_days // 2)).mean()
    return out.assign(ACC_smooth=sm(out.ACC), RMSE_smooth=sm(out.RMSE))


def block_bootstrap_ci(x: np.ndarray, block_len: int, n: int, seed: int) -> tuple[float, float]:
    x = x[np.isfinite(x)]
    if x.size < 2:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    starts = np.arange(max(1, x.size - block_len + 1))
    nblk = int(np.ceil(x.size / block_len))
    means = [x[np.concatenate([np.arange(s, min(s + block_len, x.size)) for s in rng.choice(starts, nblk)])[:x.size]]
             .mean() for _ in range(n)]
    return tuple(np.percentile(means, [2.5, 97.5]))


def window_table(skill: xr.Dataset, cfg: WindSkillConfig) -> pd.DataFrame:
    rows = []
    for a, b in cfg.lead_windows:
        seg = skill.ACC.sel(lead=slice(a, b)).values
        lo, hi = block_bootstrap_ci(seg, cfg.block_len, cfg.n_boot, cfg.seed)
        rows.append(dict(window=f"{a:02d}-{b:02d}", ACC=float(np.nanmean(seg)), ci_low=lo, ci_high=hi,
                         n_days=int(np.isfinite(seg).sum())))
    return pd.DataFrame(rows)


def efold_lead(acc: np.ndarray) -> float:
    """First lead (linearly interpolated, in days) at which ``acc`` falls below acc[0] / e."""
    acc = np.asarray(acc, dtype=float)
    if not np.isfinite(acc[0]) or acc[0] <= 0:
        return np.nan
    thr = acc[0] / np.e
    below = np.flatnonzero(acc < thr)
    if below.size == 0:
        return np.nan
    k = int(below[0])
    return float(k - 1 + (acc[k - 1] - thr) / (acc[k - 1] - acc[k]))


def skill_timescales(skills: dict[str, xr.Dataset]) -> pd.DataFrame:
    """Day-0 ACC and e-folding lead of every experiment (index experiment)."""
    return pd.DataFrame({e: dict(ACC_day0=float(s.ACC.sel(lead=0)), efold_days=efold_lead(s.ACC.values))
                         for e, s in skills.items()}).T


def crossover_lead(skills: dict[str, xr.Dataset], a: str = "DART", b: str = "ERA5", persist: int = 3) -> float:
    """First lead from which the mean ACC of atm group ``a`` exceeds group ``b`` for ``persist`` days."""
    g = groups_by(skills, "atm")
    if not g[a] or not g[b]:
        return np.nan
    d = np.mean([skills[e].ACC.values for e in g[a]], 0) - np.mean([skills[e].ACC.values for e in g[b]], 0)
    lead = next(iter(skills.values())).lead.values
    for k in range(len(d) - persist + 1):
        if np.all(d[k:k + persist] > 0):
            return float(lead[k])
    return np.nan


def plot_wind_skill_groups(skills: dict[str, xr.Dataset], variable: str, *, max_lead: int = 30,
                           a: str = "DART", b: str = "ERA5") -> plt.Figure:
    """(a) daily ACC by atmosphere initial state (mean, member range); (b) group ``a`` - ``b``
    as the factorial effect with its 95% interval (noise from the interaction residual)."""
    sk = {e: s.sel(lead=slice(0, max_lead)) for e, s in skills.items()}
    lead = next(iter(sk.values())).lead.values
    curves = {e: s.ACC.values for e, s in sk.items()}
    fit = additive_fit(curves)
    k = fit["contrasts"].index(f"atm {a}−{b}")
    est, half = fit["effect"][k], fit["t_crit"] * fit["se"][k]
    cross = crossover_lead(skills, a, b)
    fig, axes = plt.subplots(1, 2, figsize=(15, 4.6))
    plot_group_curves(axes[0], lead, curves)
    axes[0].axhline(0, color="k", linewidth=0.5)
    axes[0].set_ylabel("ACC")
    axes[0].set_title(f"(a) {variable} daily ACC by atmosphere initial state", loc="left")
    axes[0].legend(fontsize=9, frameon=False)
    axes[1].plot(lead, est, color="k", linewidth=2, label=f"atm {a} − {b} (additive fit)")
    axes[1].fill_between(lead, est - half, est + half, color="0.6", alpha=0.4, linewidth=0, label="95% CI")
    axes[1].axhline(0, color="k", linewidth=0.5)
    if np.isfinite(cross):
        for ax in axes:
            ax.axvline(cross, color=GROUP_COLORS[a], linestyle="--", linewidth=1)
        axes[1].text(cross + 0.3, 0.95, f"{a} ahead from day {cross:g}", transform=axes[1].get_xaxis_transform(),
                     va="top", fontsize=9, color=GROUP_COLORS[a])
    axes[1].set_ylabel(r"$\Delta$ACC")
    axes[1].set_title(f"(b) Effect of the atmosphere initial state, {a} − {b}", loc="left")
    axes[1].legend(fontsize=9, frameon=False, loc="lower right")
    for ax in axes:
        ax.set_xlabel("Lead day")
        ax.set_xlim(0, max_lead)
        ax.grid(True, alpha=0.3)
    fig.text(0.01, 0.005, "Verified against ERA5, which is also the ERA5 runs' initial state (day-0 ACC favours them).",
             fontsize=8)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    return fig


def plot_wind_skill(skills: dict[str, xr.Dataset], tables: pd.DataFrame, variable: str, *, style_map=None) -> plt.Figure:
    """Smoothed daily ACC and RMSE vs lead, and window-mean ACC with bootstrap CI."""
    style_map = style_map or STYLE_MAP
    fig, axes = plt.subplots(1, 3, figsize=(18, 4.6))
    for name, ds in skills.items():
        st = style_map.get(name, {})
        axes[0].plot(ds.lead, ds.ACC_smooth, label=name, **st)
        axes[1].plot(ds.lead, ds.RMSE_smooth, **st)
    exps = list(skills)
    windows = list(dict.fromkeys(tables.window))
    x = np.arange(len(windows))
    width = 0.8 / len(exps)
    for j, e in enumerate(exps):
        sub = tables[tables.experiment == e].set_index("window").loc[windows]
        axes[2].bar(x - 0.4 + width * (j + 0.5), sub.ACC, width, color=style_map.get(e, {}).get("color"),
                    yerr=[sub.ACC - sub.ci_low, sub.ci_high - sub.ACC], capsize=2, label=e)
    axes[2].set_xticks(x)
    axes[2].set_xticklabels([f"days {w}" for w in windows])
    for ax, t in zip(axes, ("(a) ACC (smoothed)", "(b) RMSE (m s$^{-1}$, smoothed)", "(c) Window-mean ACC, 95% block CI")):
        ax.set_title(f"{t}: {variable}", loc="left")
        ax.grid(True, alpha=0.3)
    axes[0].axhline(0, color="k", linewidth=0.5)
    axes[0].set_xlabel("Lead day")
    axes[1].set_xlabel("Lead day")
    axes[0].legend(ncol=2, fontsize=8, frameon=False)
    fig.tight_layout()
    return fig
