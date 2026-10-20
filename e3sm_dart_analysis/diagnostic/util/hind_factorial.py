"""Atmosphere x land initial-state attribution of the AMIP hindcast scores.

The experiments (``configs.hind_experiment_config.INITIAL_STATES``) fill 8 of
the 3 x 3 cells {piC, DART, ERA5} atm x {piC, DART, ERA5} lnd, one run each.
Any score (a number per experiment, or a curve such as ACC per lead day) is
fitted with the additive model

    score = mu + a[atm] + l[lnd] + residual,     a[piC] = l[piC] = 0,

by least squares: a two-way ANOVA without replication. With one run per cell
the residual (interaction plus chaotic run-to-run noise, 8 - 5 = 3 degrees of
freedom) is the only noise estimate, so it sets the standard error of every
effect; an effect is flagged when its 95% interval (Student t, df = 3,
t = 3.18) excludes zero. Differences between single runs are not shown as
significant unless they exceed this floor.

Effects (contrasts of the fitted parameters):

* ``atm DART−piC``, ``atm ERA5−piC``, ``atm DART−ERA5``: atmosphere initial state,
  averaged over the land states;
* ``lnd DART−piC``, ``lnd ERA5−piC``, ``lnd ERA5−DART``: land initial state,
  averaged over the atmosphere states.

``plot_factorial_scorecard`` draws, per metric, the 3 x 3 grid of scores and
the effects with their 95% intervals; ``plot_group_curves`` draws curves
grouped by atmosphere initial state (mean and range of the members).

``model_range_check``: for statistics of the free-running model (coupling,
balance, climate), where the initial state should not matter, the 8 runs act
as a pseudo-ensemble; the observed value is compared with their range.
"""
from __future__ import annotations

import re
from typing import Mapping

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from configs.hind_experiment_config import INIT_LEVELS, INITIAL_STATES

PARAMS = ("mu", "atm_DART", "atm_ERA5", "lnd_DART", "lnd_ERA5")
CONTRASTS = {
    "atm DART−piC": {"atm_DART": 1},
    "atm ERA5−piC": {"atm_ERA5": 1},
    "atm DART−ERA5": {"atm_DART": 1, "atm_ERA5": -1},
    "lnd DART−piC": {"lnd_DART": 1},
    "lnd ERA5−piC": {"lnd_ERA5": 1},
    "lnd ERA5−DART": {"lnd_ERA5": 1, "lnd_DART": -1},
}
GROUP_COLORS = {"piC": "0.45", "DART": "tab:red", "ERA5": "tab:blue"}


def _design(experiments) -> np.ndarray:
    rows = []
    for e in experiments:
        atm, lnd = INITIAL_STATES[e]
        rows.append([1.0, atm == "DART", atm == "ERA5", lnd == "DART", lnd == "ERA5"])
    return np.asarray(rows, dtype=float)


def _contrast_matrix() -> np.ndarray:
    C = np.zeros((len(CONTRASTS), len(PARAMS)))
    for i, c in enumerate(CONTRASTS.values()):
        for p, w in c.items():
            C[i, PARAMS.index(p)] = w
    return C


def _fit_column(X: np.ndarray, y: np.ndarray, C: np.ndarray):
    ok = np.isfinite(y)
    X, y = X[ok], y[ok]
    nan = np.full(C.shape[0], np.nan)
    if X.shape[0] == 0:
        return nan, nan, np.nan, 0
    used = X.any(axis=0)                       # drop the levels no experiment has
    Xu = X[:, used]
    rank = np.linalg.matrix_rank(Xu)
    df = X.shape[0] - rank
    beta = np.full(X.shape[1], np.nan)
    beta[used] = np.linalg.lstsq(Xu, y, rcond=None)[0]
    est = C @ np.nan_to_num(beta)
    # a contrast is estimable when it lies in the row space of the design (e.g. not atm DART−piC
    # when no piC-atm run has a value)
    Cu = C[:, used]
    estimable = np.isclose(Cu @ np.linalg.pinv(Xu) @ Xu, Cu, atol=1e-8).all(axis=1)
    est[~estimable | (C[:, ~used] != 0).any(axis=1)] = np.nan
    if df <= 0:
        return est, nan, np.nan, df
    resid = y - Xu @ beta[used]
    sd = float(np.sqrt(np.sum(resid ** 2) / df))
    cov = np.zeros((X.shape[1], X.shape[1]))
    cov[np.ix_(used, used)] = np.linalg.pinv(Xu.T @ Xu)
    se = sd * np.sqrt(np.einsum("ij,jk,ik->i", C, cov, C))
    se[np.isnan(est)] = np.nan
    return est, se, sd, df


def additive_fit(values: Mapping[str, np.ndarray | float]) -> dict:
    """Fit ``values`` (experiment -> scalar or array of equal shape) with the additive model.

    Returns ``effect``/``se`` (contrast, *shape), ``resid_sd`` and ``df`` (*shape),
    ``t_crit`` (*shape) and the contrast names. NaN values drop that experiment
    from that element's fit.
    """
    exps = [e for e in values if e in INITIAL_STATES]
    Y = np.stack([np.asarray(values[e], dtype=float) for e in exps])
    shape = Y.shape[1:]
    Y = Y.reshape(len(exps), -1)
    X, C = _design(exps), _contrast_matrix()
    out = [_fit_column(X, Y[:, j], C) for j in range(Y.shape[1])]
    est = np.array([o[0] for o in out]).T.reshape((len(CONTRASTS), *shape))
    se = np.array([o[1] for o in out]).T.reshape((len(CONTRASTS), *shape))
    sd = np.array([o[2] for o in out]).reshape(shape)
    df = np.array([o[3] for o in out]).reshape(shape)
    t = np.where(df > 0, stats.t.ppf(0.975, np.maximum(df, 1)), np.nan)
    return dict(contrasts=list(CONTRASTS), effect=est, se=se, resid_sd=sd, df=df, t_crit=t)


def factorial_effects(scores: pd.DataFrame) -> pd.DataFrame:
    """Effects of every metric (``scores``: index experiment, one column per metric).

    Columns: metric, effect, estimate, se, ci95 (half-width), significant,
    resid_sd (noise floor of one run), df.
    """
    rows = []
    for metric in scores.columns:
        fit = additive_fit(scores[metric].to_dict())
        for i, name in enumerate(fit["contrasts"]):
            est, se = float(fit["effect"][i]), float(fit["se"][i])
            half = float(fit["t_crit"]) * se
            rows.append(dict(metric=metric, effect=name, estimate=est, se=se, ci95=half,
                             significant=bool(np.isfinite(half) and abs(est) > half),
                             resid_sd=float(fit["resid_sd"]), df=int(fit["df"])))
    return pd.DataFrame(rows)


def factorial_grid(score: Mapping[str, float]) -> pd.DataFrame:
    """Scores laid out on the atm (rows) x lnd (columns) grid; NaN where no run."""
    grid = pd.DataFrame(np.nan, index=pd.Index(INIT_LEVELS, name="atm"), columns=pd.Index(INIT_LEVELS, name="lnd"))
    for e, v in score.items():
        if e in INITIAL_STATES:
            grid.loc[INITIAL_STATES[e]] = v
    return grid


def pair_noise(values: Mapping[str, np.ndarray | float]) -> tuple[np.ndarray, np.ndarray]:
    """95% half-width for the difference of two single runs, t * sqrt(2) * resid_sd, and resid_sd."""
    fit = additive_fit(values)
    return fit["t_crit"] * np.sqrt(2.0) * fit["resid_sd"], fit["resid_sd"]


def groups_by(experiments, factor: str = "atm") -> dict[str, list[str]]:
    k = 0 if factor == "atm" else 1
    return {lev: [e for e in experiments if e in INITIAL_STATES and INITIAL_STATES[e][k] == lev]
            for lev in INIT_LEVELS}


def plot_group_curves(ax, x, curves: Mapping[str, np.ndarray], *, factor: str = "atm", alpha: float = 0.18,
                      linewidth: float = 2.2, label_suffix: str = "") -> None:
    """Mean (line) and member range (band) of the curves grouped by ``factor`` initial state."""
    for lev, members in groups_by(curves, factor).items():
        if not members:
            continue
        Y = np.stack([np.asarray(curves[m], dtype=float) for m in members])
        c = GROUP_COLORS[lev]
        ax.plot(x, np.nanmean(Y, 0), color=c, linewidth=linewidth,
                label=f"{factor} {lev} (n={len(members)}){label_suffix}")
        if len(members) > 1:
            ax.fill_between(x, np.nanmin(Y, 0), np.nanmax(Y, 0), color=c, alpha=alpha, linewidth=0)


def plot_factorial_scorecard(scores: pd.DataFrame, effects: pd.DataFrame | None = None, *,
                             labels: Mapping[str, str] | None = None, fmt: str = "{:.2f}",
                             cmap: str = "viridis", higher_is_better: Mapping[str, bool] | None = None,
                             title: str = "", font_scale: float = 11) -> plt.Figure:
    """One row per metric: (left) the atm x lnd grid of scores, (right) effects ± 95% CI.

    Effects whose interval excludes zero are coloured (green: towards better,
    when ``higher_is_better`` gives the metric's sense; otherwise black); the rest grey.
    """
    effects = factorial_effects(scores) if effects is None else effects
    labels = labels or {}
    higher_is_better = higher_is_better or {}
    metrics = list(scores.columns)
    fig, axes = plt.subplots(len(metrics), 2, figsize=(11.5, 2.9 * len(metrics)), squeeze=False,
                             gridspec_kw=dict(width_ratios=(1, 1.5)))
    for i, m in enumerate(metrics):
        ax, ax2 = axes[i]
        grid = factorial_grid(scores[m].to_dict())
        vals = grid.values
        ax.imshow(vals, cmap=cmap if higher_is_better.get(m, True) else f"{cmap}_r", aspect="auto")
        for (r, c), v in np.ndenumerate(vals):
            ax.text(c, r, "—" if np.isnan(v) else fmt.format(v), ha="center", va="center",
                    fontsize=font_scale, color="k" if np.isnan(v) else "w",
                    bbox=None if np.isnan(v) else dict(boxstyle="round,pad=0.15", facecolor="0.2", alpha=0.35,
                                                       linewidth=0))
        ax.set_xticks(range(3), [f"lnd {lev}" for lev in INIT_LEVELS], fontsize=font_scale * 0.9)
        ax.set_yticks(range(3), [f"atm {lev}" for lev in INIT_LEVELS], fontsize=font_scale * 0.9)
        ax.set_title(f"({chr(97 + 2 * i)}) {labels.get(m, m)}", loc="left", fontsize=font_scale)
        sub = effects[effects.metric == m].set_index("effect").loc[list(CONTRASTS)]
        y = np.arange(len(sub))[::-1]
        sense = higher_is_better.get(m)
        colors = []
        for est, sig in zip(sub.estimate, sub.significant):
            if not sig:
                colors.append("0.7")
            elif sense is None:
                colors.append("k")
            else:
                colors.append("tab:green" if (est > 0) == sense else "tab:orange")
        ax2.barh(y, sub.estimate, xerr=sub.ci95, color=colors, capsize=3, error_kw=dict(ecolor="0.3"))
        ax2.axvline(0, color="k", linewidth=0.6)
        ax2.axhline(2.5, color="0.6", linewidth=0.6, linestyle=":")
        ax2.set_yticks(y, list(sub.index), fontsize=font_scale * 0.9)
        rs = sub.resid_sd.iloc[0]
        ax2.set_title(f"({chr(98 + 2 * i)}) effects ± 95% CI (noise σ = {fmt.format(rs)}, df = {sub.df.iloc[0]})",
                      loc="left", fontsize=font_scale)
        ax2.grid(True, axis="x", alpha=0.3)
    if title:
        fig.suptitle(title, fontsize=font_scale * 1.15)
    fig.text(0.01, 0.003, "Additive fit score = μ + atm + lnd; noise = interaction residual of the 8 single runs. "
             "Coloured: interval excludes 0 (green = better, orange = worse).", fontsize=font_scale * 0.75)
    fig.tight_layout(rect=(0, 0.02, 1, 0.98 if title else 1))
    return fig


def model_range_check(table: pd.DataFrame, refs, experiments) -> pd.DataFrame:
    """Per metric (column of ``table``, index source): reference values, model mean/std/min/max over
    ``experiments`` and whether each reference lies inside the model range (a pseudo-ensemble check)."""
    rows = []
    for m in table.columns:
        x = table.loc[list(experiments), m].astype(float)
        row = dict(metric=m, model_mean=x.mean(), model_std=x.std(ddof=1), model_min=x.min(), model_max=x.max())
        for r in refs:
            v = float(table.loc[r, m])
            row[r] = v
            row[f"{r}_in_range"] = bool(x.min() <= v <= x.max())
        rows.append(row)
    return pd.DataFrame(rows).set_index("metric")


def relabel_figure(fig: plt.Figure, labels: Mapping[str, str]) -> plt.Figure:
    """Replace experiment keys (whole words) by their display labels in every text of ``fig``:
    titles, legends, annotations and fixed tick labels."""
    from matplotlib.ticker import FixedFormatter, FuncFormatter
    labels = {k: v for k, v in labels.items() if k != v}
    if not labels:
        return fig
    pattern = re.compile(r"\b(" + "|".join(map(re.escape, sorted(labels, key=len, reverse=True))) + r")\b")
    sub = lambda text: pattern.sub(lambda m: labels[m.group(1)], text)
    for t in fig.findobj(lambda o: hasattr(o, "get_text") and hasattr(o, "set_text")):
        t.set_text(sub(t.get_text()))
    for ax in fig.axes:  # tick labels are redrawn from the formatter (FuncFormatter for fixed labels since mpl 3.9)
        for axis in (ax.xaxis, ax.yaxis):
            fmt = axis.get_major_formatter()
            if isinstance(fmt, FixedFormatter):
                fmt.seq = [sub(str(x)) for x in fmt.seq]
            elif isinstance(fmt, FuncFormatter) and not getattr(fmt.func, "_relabelled", False):
                func = lambda x, pos=None, f=fmt.func: sub(str(f(x, pos)))
                func._relabelled = True
                fmt.func = func
    return fig
