"""Radiative/chemical tracer initialization of the F20TR hindcasts (``hind_lead_eval/08``).

HYB_BASE starts from the bare ERA5 eam.i (T, U, V, Q, PS, CLDLIQ, CLDICE): ozone, the linoz gases the radiation
reads (``rad_climate``: H2OLNZ, O3, N2OLNZ, CH4LNZ), aerosols and chemistry start near zero. ERA5CLIMFILL is the
same ERA5 state completed with the CLIM_BASE fields. Their difference is the effect of the tracer initialization;
it is traced from the tracers (ozone columns) through the clear-sky radiation to the stratospheric, tropospheric and
surface temperature and circulation, and set against the initialization signals it would contaminate.

Every member is reduced as it is read, one start at a time (``reduce_start``): regional means (region, lead),
zonal means (lead, lat), lead-window means (window, lat, lon) and, where a reference exists, the ensemble-mean
error against it by region. The statistics of a contrast a - b (``contrast_stats``) treat the members of each
ensemble as exchangeable and unpaired, and the starts as independent:

* difference of the ensemble means, averaged over the starts;
* standard error: sqrt(sum over starts of (var_a / n_a + var_b / n_b)) / n_starts, from the member variances;
* member spread: sqrt of the mean of the two member variances (the internal variability of one member);
* ``z`` = difference / standard error.

A difference is "outside the spread" when |difference| > k x spread; ``recovery`` gives the last lead day on which
it is (the defect has recovered from the next day on) and the last day on which it is significant (|z| > z_crit).
"""
from __future__ import annotations

import string
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.colors import BoundaryNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from configs.fcst_init_experiment_config import STYLE_MAP, config_label, relabel
from util.fcst_init_plotting import LEGEND_BOX, SEASON_PANEL, _lead_days, _season_grid, fig_size, fs, readable

_LETTERS = string.ascii_lowercase
GLOBE = dict(lat=(-90.0, 90.0), lon=(0.0, 360.0))
CONTRAST_COLORS = ("k", "tab:orange", "tab:blue", "tab:green", "tab:purple", "tab:brown")
CHAIN_COLORS = ("tab:purple", "tab:red", "tab:cyan", "tab:green", "tab:orange", "tab:brown", "tab:pink", "tab:olive")


# ---------------------------------------------------------------------------
# Reductions of one start (all members)
# ---------------------------------------------------------------------------
def _area_mean(x: np.ndarray, w: np.ndarray) -> np.ndarray:
    """cos-lat weighted mean over the last two axes with weights ``w`` (lat, lon), skipping NaN."""
    ok = np.isfinite(x)
    num = np.sum(np.where(ok, x, 0.0) * w, axis=(-2, -1))
    den = np.sum(np.where(ok, w, 0.0), axis=(-2, -1))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num / den, np.nan)


def region_weights(masks: xr.DataArray) -> np.ndarray:
    """(region, lat, lon) cos-lat weights inside each region mask, zero outside."""
    return np.cos(np.deg2rad(masks.lat.values)).clip(min=0)[None, :, None] * masks.values.astype(float)


def reduce_start(members: xr.DataArray, masks: xr.DataArray, windows: dict,
                 reference: xr.DataArray | None = None) -> xr.Dataset:
    """Reductions of one start's (member, lead, lat, lon) field: ``reg`` (member, lead, region) regional means,
    ``zonal`` (member, lead, lat) zonal means, ``win_mean`` / ``win_var`` (window, lat, lon) ensemble mean and
    member variance of the lead-window means, and with a (lead, lat, lon) ``reference`` the regional ``err``
    (lead, region; ensemble mean minus reference), ``sqerr`` (its square) and ``var`` (member variance), averaged
    over the points where the reference and every member are valid."""
    members = members.transpose("member", "lead", "lat", "lon")
    f = np.asarray(members.values, dtype=np.float64)
    w = region_weights(masks)
    leads = members.lead.values
    reg = np.stack([_area_mean(f, wr) for wr in w], axis=-1)                     # (member, lead, region)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)   # all-NaN columns
        zonal = np.nanmean(f, axis=-1)                                           # (member, lead, lat)
        win = np.stack([np.nanmean(f[:, (leads >= a) & (leads <= b)], axis=1) for a, b in windows.values()], 1)
    coords = dict(member=members.member.values, lead=leads, region=masks.region.values, lat=members.lat.values,
                  lon=members.lon.values, window=list(windows))
    ds = xr.Dataset({"reg": (("member", "lead", "region"), reg.astype("float32")),
                     "zonal": (("member", "lead", "lat"), zonal.astype("float32")),
                     "win_mean": (("window", "lat", "lon"), win.mean(0).astype("float32")),
                     "win_var": (("window", "lat", "lon"), win.var(0, ddof=1).astype("float32"))},
                    coords=coords)
    if reference is not None:
        o = np.asarray(reference.transpose("lead", "lat", "lon").values, dtype=np.float64)
        valid = np.isfinite(o) & np.isfinite(f).all(axis=0)
        err = np.where(valid, f.mean(0) - o, np.nan)
        var = np.where(valid, f.var(0, ddof=1), np.nan)
        for name, x in (("err", err), ("sqerr", err ** 2), ("var", var)):
            ds[name] = (("lead", "region"), np.stack([_area_mean(x, wr) for wr in w], axis=-1).astype("float32"))
    return ds


def reduce_reference(reference: xr.DataArray, masks: xr.DataArray) -> xr.Dataset:
    """``reg`` (init, lead, region) and ``zonal`` (init, lead, lat) of an (init, lead, lat, lon) reference."""
    o = np.asarray(reference.transpose("init", "lead", "lat", "lon").values, dtype=np.float64)
    reg = np.stack([_area_mean(o, wr) for wr in region_weights(masks)], axis=-1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)   # all-NaN columns
        zonal = np.nanmean(o, axis=-1)
    return xr.Dataset({"reg": (("init", "lead", "region"), reg.astype("float32")),
                       "zonal": (("init", "lead", "lat"), zonal.astype("float32"))},
                      coords=dict(init=reference.init.values, lead=reference.lead.values,
                                  region=masks.region.values, lat=reference.lat.values))


# ---------------------------------------------------------------------------
# Contrasts, recovery and skill
# ---------------------------------------------------------------------------
def _member_stats(field: xr.DataArray) -> tuple[xr.DataArray, xr.DataArray, int]:
    return field.mean("member"), field.var("member", ddof=1), field.sizes["member"]


def contrast_stats(a: xr.DataArray, b: xr.DataArray) -> xr.Dataset:
    """``delta``, ``se``, ``spread``, ``z`` of a - b from (init, member, ...) fields of two ensembles started
    on the same dates (members unpaired; starts independent)."""
    ma, va, na = _member_stats(a)
    mb, vb, nb = _member_stats(b)
    return _stats(ma - mb, va / na + vb / nb, 0.5 * (va + vb))


def contrast_stats_windows(a: xr.Dataset, b: xr.Dataset) -> xr.Dataset:
    """As ``contrast_stats`` for the lead-window maps of two summaries (``win_mean``, ``win_var``)."""
    na, nb = a.sizes["member"], b.sizes["member"]
    return _stats(a.win_mean - b.win_mean, a.win_var / na + b.win_var / nb, 0.5 * (a.win_var + b.win_var))


def _stats(diff: xr.DataArray, var_of_diff: xr.DataArray, member_var: xr.DataArray) -> xr.Dataset:
    n = diff.sizes["init"]
    delta = diff.mean("init")
    se = np.sqrt(var_of_diff.sum("init")) / n
    spread = np.sqrt(member_var.mean("init"))
    return xr.Dataset({"delta": delta, "se": se, "spread": spread, "z": delta / se.where(se > 0)})


def last_lead(flag: np.ndarray, leads: np.ndarray) -> int:
    """Last lead day on which ``flag`` is True (0 when never)."""
    hit = np.flatnonzero(np.asarray(flag, dtype=bool))
    return int(leads[hit[-1]]) if hit.size else 0


def recovery(stats: xr.Dataset, k: float, z_crit: float) -> dict:
    """Day-1 and peak difference, and the last lead days on which |delta| > k x spread and |z| > z_crit
    (``stats``: (lead,) ``contrast_stats`` of one region)."""
    leads = stats.lead.values
    d = stats.delta.values
    peak = int(np.nanargmax(np.abs(d))) if np.isfinite(d).any() else 0
    return dict(day1=float(d[0]), peak=float(d[peak]), peak_lead=int(leads[peak]),
                spread_at_peak=float(stats.spread.values[peak]),
                outside_spread_until=last_lead(np.abs(d) > k * stats.spread.values, leads),
                significant_until=last_lead(np.abs(stats.z.values) > z_crit, leads))


def pooled_skill(summary: xr.Dataset) -> xr.Dataset:
    """(region, lead) ensemble-mean bias, RMSE, spread (unbiased member variance x (N+1)/N) and spread/RMSE,
    pooled over the starts of a summary with ``err``, ``sqerr``, ``var``."""
    n = summary.sizes["member"]
    ds = xr.Dataset({"bias": summary["err"].mean("init"), "rmse": np.sqrt(summary["sqerr"].mean("init")),
                     "spread": np.sqrt(summary["var"].mean("init") * (n + 1) / n)})
    ds["spread_ratio"] = ds.spread / ds.rmse
    return ds.transpose("region", "lead")


def window_change(skill: pd.DataFrame, control: str, windows: dict, *, region: str, metric: str = "rmse",
                  variables=None) -> dict[str, pd.DataFrame]:
    """{season: (config, "<variable>\\n<window>") table} of the % change of the lead-window mean ``metric``
    against ``control`` (configurations other than the control as rows)."""
    t = skill[skill.region == region]
    variables = list(dict.fromkeys(t.variable)) if variables is None else [v for v in variables if v in set(t.variable)]
    out = {}
    for s, ts in t.groupby("season", sort=False):
        cols = {}
        for v in variables:
            tv = ts[ts.variable == v]
            for wname, (a, b) in windows.items():
                m = tv[(tv.lead >= a) & (tv.lead <= b)].groupby("config", sort=False)[metric].mean()
                if control in m:
                    cols[f"{v}\n{wname}"] = 100.0 * (m.drop(control) - m[control]) / m[control]
        out[s] = pd.DataFrame(cols)
    return out


def symmetric_levels(values, n: int = 10, q: float = 98.0) -> np.ndarray:
    """2n+1 levels symmetric about 0 up to the ``q`` percentile of |values| (rounded to 2 significant digits)."""
    v = np.abs(np.asarray(values, dtype=float))
    top = np.nanpercentile(v, q) if np.isfinite(v).any() else 1.0
    top = float(f"{top:.2g}") if top > 0 else 1.0
    return np.linspace(-top, top, 2 * n + 1)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def contrast_label(pair) -> str:
    return relabel(f"{pair[0]} - {pair[1]}")


@readable
def plot_regional_series(series: dict, obs: dict, variables, *, seasons, region_label: str, units: dict,
                         quantiles=(0.1, 0.9), max_lead: int = 56) -> plt.Figure:
    """Regional means by lead day: rows ``variables``, columns ``seasons``. ``series[(season, config, variable)]``:
    (init, member, lead); each configuration's mean over starts and members (line) and the ``quantiles`` of its
    member-starts (shaded); ``obs[(season, variable)]``: (init, lead) reference, black."""
    configs = list(dict.fromkeys(c for _, c, _ in series))
    fig, axes = _season_grid(variables, seasons)
    for i, v in enumerate(variables):
        for j, s in enumerate(seasons):
            ax = axes[i, j]
            for c in configs:
                if (s, c, v) not in series:
                    continue
                d = series[s, c, v].stack(sample=("init", "member"))
                st = STYLE_MAP.get(c, dict(color="0.3"))
                ax.fill_between(d.lead, d.quantile(quantiles[0], "sample"), d.quantile(quantiles[1], "sample"),
                                color=st["color"], alpha=0.15, linewidth=0)
                ax.plot(d.lead, d.mean("sample"), **st)
            if (s, v) in obs:
                o = obs[s, v].mean("init")
                ax.plot(o.lead, o, **STYLE_MAP["ERA5"])
            _lead_days(ax, max_lead)
        axes[i, 0].set_ylabel(f"{v} [{units.get(v, '')}]")
    for ax in axes[-1]:
        ax.set_xlabel("Lead day")
    handles = [Line2D([], [], **STYLE_MAP.get(c, dict(color="0.3"))) for c in configs]
    labels = [config_label(c) for c in configs]
    if any(k[1] in variables for k in obs):
        handles.append(Line2D([], [], **STYLE_MAP["ERA5"]))
        labels.append("ERA5")
    handles.append(Patch(color="0.5", alpha=0.3))
    labels.append(f"{quantiles[0]:.0%}-{quantiles[1]:.0%} of member-starts")
    fig.legend(handles, labels, loc="outside lower center", ncols=min(len(labels), 6), **LEGEND_BOX)
    fig.suptitle(f"Regional mean by lead day (mean over starts and members), {region_label}")
    return fig


@readable
def plot_contrast_curves(table: pd.DataFrame, variables, *, contrasts: dict, seasons, region: str, region_label: str,
                         units: dict, k: float, z_crit: float, max_lead: int = 56) -> plt.Figure:
    """Contrasts a - b by lead day (rows ``variables``, columns ``seasons``): the difference with its +-z_crit
    standard-error band, and the +-k x member spread of the first contrast (grey), its recovery day dotted."""
    t = table[table.region == region]
    names = [n for n in contrasts if n in set(t.contrast)]
    fig, axes = _season_grid(variables, seasons)
    for i, v in enumerate(variables):
        for j, s in enumerate(seasons):
            ax = axes[i, j]
            tv = t[(t.variable == v) & (t.season == s)]
            for n, color in zip(names, CONTRAST_COLORS):
                d = tv[tv.contrast == n].sort_values("lead")
                if d.empty:
                    continue
                if n == names[0]:
                    ax.fill_between(d.lead, -k * d.spread, k * d.spread, color="0.85", linewidth=0, zorder=0)
                    out = d.lead[np.abs(d.delta) > k * d.spread]
                    if len(out) and out.max() < max_lead:
                        ax.axvline(out.max() + 1, color=color, linestyle=":", linewidth=1.4)
                ax.fill_between(d.lead, d.delta - z_crit * d.se, d.delta + z_crit * d.se, color=color, alpha=0.25,
                                linewidth=0)
                ax.plot(d.lead, d.delta, color=color, linewidth=2.2 if n == names[0] else 1.6)
            ax.axhline(0, color="k", linewidth=0.8)
            _lead_days(ax, max_lead)
        axes[i, 0].set_ylabel(f"Δ{v} [{units.get(v, '')}]")
    for ax in axes[-1]:
        ax.set_xlabel("Lead day")
    handles = [Line2D([], [], color=c, linewidth=2) for _, c in zip(names, CONTRAST_COLORS)]
    labels = [contrast_label(contrasts[n]) for n in names]
    handles += [Patch(color="0.85"), Line2D([], [], color=CONTRAST_COLORS[0], linestyle=":")]
    labels += [f"±{k:g} x member spread", "recovered (inside the spread from then on)"]
    fig.legend(handles, labels, loc="outside lower center", ncols=min(len(labels), 3), **LEGEND_BOX)
    fig.suptitle(f"Ensemble-mean differences by lead day (bands: ±{z_crit:g} standard errors), {region_label}")
    return fig


@readable
def plot_recovery_chain(table: pd.DataFrame, chain, *, contrast: str, pair, seasons, k: float,
                        max_lead: int = 56) -> plt.Figure:
    """The difference of ``contrast`` normalized by its peak |difference| for each (variable, region) of
    ``chain``, one panel per season; dotted: the day each recovers within k x member spread."""
    fig, axes = plt.subplots(1, len(seasons), figsize=fig_size(6.5 * len(seasons), 4.6), squeeze=False,
                             sharey=True, layout="constrained")
    t = table[table.contrast == contrast]
    for j, s in enumerate(seasons):
        ax = axes[0, j]
        for (v, r), color in zip(chain, CHAIN_COLORS):
            d = t[(t.season == s) & (t.variable == v) & (t.region == r)].sort_values("lead")
            if d.empty:
                continue
            peak = np.nanmax(np.abs(d.delta))
            ax.plot(d.lead, d.delta / peak, color=color, linewidth=2)
            out = d.lead[np.abs(d.delta) > k * d.spread]
            if len(out) and out.max() < max_lead:
                ax.axvline(out.max() + 1, color=color, linestyle=":", linewidth=1.4)
        ax.axhline(0, color="k", linewidth=0.8)
        _lead_days(ax, max_lead)
        ax.set_xlabel("Lead day")
        ax.set_title(f"({_LETTERS[j]})", loc="left")
        ax.set_title(SEASON_PANEL.get(s, s), loc="right")
    axes[0, 0].set_ylabel("Δ / max |Δ|")
    handles = [Line2D([], [], color=c, linewidth=2) for _, c in zip(chain, CHAIN_COLORS)]
    fig.legend(handles, [f"{v} ({r})" for v, r in chain], loc="outside lower center", ncols=min(len(chain), 3),
               **LEGEND_BOX)
    fig.suptitle(f"From tracers to the surface: {contrast_label(pair)}, normalized\n"
                 f"(dotted: inside ±{k:g} x member spread from then on)")
    return fig


@readable
def plot_zonal_lead(panels: dict, variables, *, seasons, title: str, units: dict, k: float,
                    max_lead: int = 56) -> plt.Figure:
    """Zonal-mean difference, latitude x lead day: rows ``variables``, columns ``seasons``;
    ``panels[(season, variable)]``: (lead, lat) ``contrast_stats``. Contour: |Δ| = k x member spread."""
    fig, axes = plt.subplots(len(variables), len(seasons), figsize=fig_size(5.2 * len(seasons), 2.7 * len(variables)),
                             squeeze=False, sharex=True, sharey=True, layout="constrained")
    n = 0
    for i, v in enumerate(variables):
        levels = symmetric_levels(np.concatenate([panels[s, v].delta.values.ravel() for s in seasons
                                                  if (s, v) in panels]))
        cf = None
        for j, s in enumerate(seasons):
            ax = axes[i, j]
            if (s, v) not in panels:
                ax.set_visible(False)
                continue
            d = panels[s, v].transpose("lat", "lead")
            cf = ax.pcolormesh(d.lead, d.lat, d.delta, cmap="RdBu_r",
                               norm=BoundaryNorm(levels, plt.get_cmap("RdBu_r").N, extend="both"), shading="nearest",
                               rasterized=True)
            ratio = (np.abs(d.delta) / (k * d.spread)).where(d.spread > 0)
            if np.isfinite(ratio).any():
                ax.contour(d.lead, d.lat, ratio, levels=[1.0], colors="k", linewidths=0.8)
            ax.set_title(f"({_LETTERS[n % 26]}) Δ{v}", loc="left")
            ax.set_title(SEASON_PANEL.get(s, s), loc="right")
            ax.set_xlim(1, max_lead)
            ax.set_xticks(range(7, max_lead + 1, 7))
            ax.set_yticks([-60, -30, 0, 30, 60])
            n += 1
        if cf is not None:
            fig.colorbar(cf, ax=axes[i, :].tolist(), shrink=0.9, pad=0.01, label=units.get(v, ""))
        axes[i, 0].set_ylabel("Latitude")
    for ax in axes[-1]:
        ax.set_xlabel("Lead day")
    fig.suptitle(f"{title}\n(contour: |Δ| = {k:g} x member spread)")
    return fig


__all__ = ["GLOBE", "reduce_start", "reduce_reference", "contrast_stats", "contrast_stats_windows", "recovery",
           "pooled_skill", "window_change", "symmetric_levels", "contrast_label", "plot_regional_series",
           "plot_contrast_curves", "plot_recovery_chain", "plot_zonal_lead"]
