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
the single realization).
"""
from __future__ import annotations

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from configs.hind_experiment_config import STYLE_MAP
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
    lead_windows: tuple[tuple[int, int], ...] = ((5, 14), (15, 28), (29, 42))
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
