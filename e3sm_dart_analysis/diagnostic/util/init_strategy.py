"""Initialization-strategy skill of the F20TR hindcasts (``hind_lead_eval/05a-05d``).

Scores of several ensembles (experiments) against one verification, and of strategy contrasts
``a - b`` between them, with a paired bootstrap that draws the same starts (within each season) and
the same member indices for every experiment:

* **pattern ACC**: centred, area-weighted spatial correlation over the domain of the ensemble-mean
  anomaly, for each start and lead, averaged over the starts;
* **temporal ACC**: local correlation across the starts at every grid cell of the ensemble-mean
  anomaly, area-weighted over the domain (area-averaged local temporal ACC, not the ACC of the
  domain-mean series);
* **CRPSS**: ``1 - CRPS / CRPS_clim`` of the area-weighted mean fair CRPS, summed over the starts,
  against a common climatological ensemble (observed anomalies of other years on the same calendar
  days);
* **held-out-member ACC**: pattern ACC of the mean of the other members against each member of the
  same experiment in turn (within-ensemble agreement).

Inputs are on the domain's points: forecast anomalies ``(start, member, step, point)``, observed
anomalies ``(start, step, point)``, climatological ensemble ``(start, year, step, point)`` and area
weights ``(point,)``. Pattern ACC and area-mean CRPS are linear or quadratic in the member counts of a
bootstrap draw, so they are computed from per-start member matrices and every draw costs only a few
small products; the temporal ACC is computed from the drawn ensemble means in batches of draws.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

QUANTILES = (2.5, 97.5)


# ---------------------------------------------------------------------------
# Draws
# ---------------------------------------------------------------------------
def season_strata(start_seasons) -> list[np.ndarray]:
    """Indices of the starts of each season (starts are resampled within their season)."""
    s = np.asarray(start_seasons)
    return [np.flatnonzero(s == x) for x in dict.fromkeys(s)]


def draws(start_seasons, n_members: int, n_boot: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Paired bootstrap draws: start counts ``(n_boot, start)`` (resampled within seasons) and member
    counts ``(n_boot, member)``, the same for every experiment."""
    rng = np.random.default_rng(seed)
    strata = season_strata(start_seasons)
    n_start = sum(len(s) for s in strata)
    k = np.zeros((n_boot, n_start))
    c = np.zeros((n_boot, n_members))
    for b in range(n_boot):
        ii = np.concatenate([rng.choice(s, s.size) for s in strata])
        mm = rng.integers(0, n_members, n_members)
        k[b] = np.bincount(ii, minlength=n_start)
        c[b] = np.bincount(mm, minlength=n_members)
    return k, c


def interval(samples: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lo, hi = np.nanpercentile(samples, QUANTILES, axis=0)
    return lo, hi


# ---------------------------------------------------------------------------
# Per-experiment statistics
# ---------------------------------------------------------------------------
def _valid_weights(x: np.ndarray, o: np.ndarray, w: np.ndarray) -> np.ndarray:
    """(start, step, point) weights: ``w`` where the observation and every member are finite, normalized."""
    ok = np.isfinite(o) & np.isfinite(x).all(1)
    ww = np.where(ok, w, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return ww / ww.sum(-1, keepdims=True)


def pattern_stats(x: np.ndarray, o: np.ndarray, w: np.ndarray) -> dict:
    """Weighted centred cross products for the pattern ACC of any member-count weighting:
    ``xo`` (start, step, member), ``xx`` (start, step, member, member), ``oo`` (start, step)."""
    ww = _valid_weights(x, o, w)
    x = np.nan_to_num(x.astype("float64"))
    o = np.nan_to_num(o.astype("float64"))
    xc = x - (x * ww[:, None]).sum(-1, keepdims=True)
    oc = o - (o * ww).sum(-1, keepdims=True)
    xs = np.moveaxis(xc * np.sqrt(ww)[:, None], 1, 2)                     # (start, step, member, point)
    os_ = oc * np.sqrt(ww)
    return dict(xo=np.einsum("slmp,slp->slm", xs, os_), xx=xs @ np.swapaxes(xs, -1, -2), oo=(os_ ** 2).sum(-1))


def pattern_acc_from_counts(st: dict, c: np.ndarray) -> np.ndarray:
    """(draw, start, step) pattern ACC of the ensemble mean weighted by member counts ``c`` (draw, member)."""
    num = np.einsum("slm,bm->bsl", st["xo"], c)
    var = np.einsum("bm,slmn,bn->bsl", c, st["xx"], c)
    with np.errstate(invalid="ignore", divide="ignore"):
        return num / np.sqrt(var * st["oo"][None])


def held_out_acc(st: dict) -> np.ndarray:
    """(member, start, step) pattern ACC of the mean of the other members against each member."""
    xx = st["xx"]
    m = xx.shape[-1]
    out = []
    for k in range(m):
        c = np.ones(m)
        c[k] = 0.0
        num = np.einsum("slm,m->sl", xx[..., k], c)
        var = np.einsum("m,slmn,n->sl", c, xx, c)
        with np.errstate(invalid="ignore", divide="ignore"):
            out.append(num / np.sqrt(var * xx[..., k, k]))
    return np.stack(out)


def crps_stats(x: np.ndarray, o: np.ndarray, w: np.ndarray) -> dict:
    """Area-weighted absolute differences for the fair CRPS of any member-count weighting:
    ``a`` (start, step, member) = sum_p w |x_m - o|, ``d`` (start, step, member, member) = sum_p w |x_m - x_n|."""
    ww = _valid_weights(x, o, w)
    x = np.nan_to_num(x.astype("float32"))
    o = np.nan_to_num(o.astype("float32"))
    a = np.einsum("smlp,slp->slm", np.abs(x - o[:, None]), ww)
    m = x.shape[1]
    d = np.zeros(a.shape + (m,))
    for k in range(m):
        d[..., k] = np.einsum("smlp,slp->slm", np.abs(x - x[:, k:k + 1]), ww)
    return dict(a=a, d=d)


def crps_from_counts(st: dict, c: np.ndarray) -> np.ndarray:
    """(draw, start, step) area-mean fair CRPS of the ensemble weighted by member counts ``c``:
    ``sum_m c_m a_m / M - sum_mn c_m c_n d_mn / (2 M (M - 1))``."""
    m = c.sum(-1)[:, None, None]
    return (np.einsum("slm,bm->bsl", st["a"], c) / m
            - np.einsum("bm,slmn,bn->bsl", c, st["d"], c) / (2 * m * (m - 1)))


def climatology_crps(y: np.ndarray, o: np.ndarray, w: np.ndarray) -> np.ndarray:
    """(start, step) area-mean fair CRPS of the climatological ensemble ``y`` (start, year, step, point)."""
    st = crps_stats(y, o, w)
    return crps_from_counts(st, np.ones((1, y.shape[1])))[0]


def temporal_acc_from_counts(x: np.ndarray, o: np.ndarray, w: np.ndarray, k: np.ndarray, c: np.ndarray,
                             batch: int = 50) -> np.ndarray:
    """(draw, step) area-weighted local temporal ACC: at every point the correlation across the starts
    (start counts ``k``) of the ensemble mean (member counts ``c``) with the observations."""
    xm = np.ascontiguousarray(np.moveaxis(x.astype("float32"), 1, 0))     # (member, start, step, point)
    m, s, l, p = xm.shape
    flat = np.nan_to_num(xm).reshape(m, -1)
    o = o.astype("float64")
    ok = np.isfinite(o).all(0) & np.isfinite(x).all((0, 1))               # (step, point)
    on = np.nan_to_num(o)
    out = np.empty((len(k), l))
    for b0 in range(0, len(k), batch):
        cb, kb = c[b0:b0 + batch], k[b0:b0 + batch]
        f = ((cb / cb.sum(-1, keepdims=True)).astype("float32") @ flat).reshape(len(cb), s, l, p).astype("float64")
        kw = (kb / kb.sum(-1, keepdims=True))[:, :, None, None]
        fa = f - (kw * f).sum(1, keepdims=True)
        oa = on[None] - (kw * on[None]).sum(1, keepdims=True)
        with np.errstate(invalid="ignore", divide="ignore"):
            r = (kw * fa * oa).sum(1) / np.sqrt((kw * fa ** 2).sum(1) * (kw * oa ** 2).sum(1))   # (draw, step, point)
        r = np.where(ok[None] & np.isfinite(r), r, np.nan)
        wr = np.where(np.isfinite(r), w, 0.0)
        out[b0:b0 + batch] = np.nansum(r * wr, -1) / wr.sum(-1)
    return out


# ---------------------------------------------------------------------------
# All scores of one verification set
# ---------------------------------------------------------------------------
def strategy_scores(anoms: dict[str, np.ndarray], obs: np.ndarray, clim: np.ndarray | None, w: np.ndarray,
                    contrasts: dict[str, tuple[str, str]], start_seasons, *, n_boot: int = 1000, seed: int = 42,
                    acc_types=("pattern", "temporal"), agreement: bool = True) -> dict:
    """Skill of every experiment and of the contrasts, with paired-bootstrap intervals.

    ``anoms[experiment]``: (start, member, step, point) forecast anomalies (members exchangeable, equal
    numbers); ``obs``: (start, step, point); ``clim``: (start, year, step, point) climatological ensemble
    for the CRPSS (None: no CRPSS); ``w``: (point,) area weights; ``contrasts``: {name: (a, b)} for a - b.
    Returns {metric: {"skill": {exp: (value, lo, hi)}, "contrast": {name: (value, lo, hi)}}}; values (step,).
    """
    exps = list(anoms)
    n_members = {a.shape[1] for a in anoms.values()}
    if len(n_members) != 1:
        raise ValueError(f"experiments have different ensemble sizes {n_members}; set max_members")
    m = n_members.pop()
    k, c = draws(start_seasons, m, n_boot, seed)
    full_k, full_c = np.ones((1, k.shape[1])), np.ones((1, m))
    pairs = {n: ab for n, ab in contrasts.items() if all(e in anoms for e in ab)}
    out = {}

    def summarize(metric, full, boot):
        res = {"skill": {}, "contrast": {}}
        for e in exps:
            res["skill"][e] = (full[e], *interval(boot[e]))
        for n, (a, b) in pairs.items():
            res["contrast"][n] = (full[a] - full[b], *interval(boot[a] - boot[b]))
        out[metric] = res

    def over_starts(v, kk):                                                # (draw, start, step) -> (draw, step)
        return np.einsum("bsl,bs->bl", v, kk) / kk.sum(-1)[:, None]

    pst = {e: pattern_stats(a, obs, w) for e, a in anoms.items()}
    if "pattern" in acc_types:
        summarize("pattern_acc", {e: over_starts(pattern_acc_from_counts(st, full_c), full_k)[0] for e, st in pst.items()},
                  {e: over_starts(pattern_acc_from_counts(st, c), k) for e, st in pst.items()})
    if "temporal" in acc_types:
        summarize("temporal_acc", {e: temporal_acc_from_counts(a, obs, w, full_k, full_c)[0] for e, a in anoms.items()},
                  {e: temporal_acc_from_counts(a, obs, w, k, c) for e, a in anoms.items()})
    if clim is not None:
        ref = climatology_crps(clim, obs, w)                               # (start, step)
        full, boot = {}, {}
        for e, a in anoms.items():
            st = crps_stats(a, obs, w)
            full[e] = 1 - over_starts(crps_from_counts(st, full_c), full_k)[0] / over_starts(ref[None], full_k)[0]
            boot[e] = 1 - over_starts(crps_from_counts(st, c), k) / over_starts(np.broadcast_to(ref, (len(k),) + ref.shape), k)
        summarize("crpss", full, boot)
    if agreement:                                                          # member index and starts resampled
        rng = np.random.default_rng(seed + 1)
        kk = rng.integers(0, m, (n_boot, m))
        per = {e: held_out_acc(st) for e, st in pst.items()}               # (member, start, step)
        full = {e: v.mean((0, 1)) for e, v in per.items()}
        boot = {e: np.stack([np.einsum("sl,s->l", v[kk[b]].mean(0), k[b]) / k[b].sum() for b in range(n_boot)])
                for e, v in per.items()}
        summarize("held_out_acc", full, boot)
    return out


def scores_table(scores: dict, steps, x_weeks, **labels) -> pd.DataFrame:
    """Long table: metric, kind (skill/contrast), name, step, x_weeks, value, ci_low, ci_high, and ``labels``."""
    rows = []
    for metric, res in scores.items():
        for kind, curves in res.items():
            for name, (v, lo, hi) in curves.items():
                for j, (s, x) in enumerate(zip(steps, x_weeks)):
                    rows.append(dict(**labels, metric=metric, kind=kind, name=name, step=int(s), x_weeks=float(x),
                                     value=float(v[j]), ci_low=float(lo[j]), ci_high=float(hi[j])))
    return pd.DataFrame(rows)


def from_table(table: pd.DataFrame, metric: str, kind: str) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    """{name: (x_weeks, value, ci_low, ci_high)} of one metric and kind from a ``scores_table``."""
    t = table[(table.metric == metric) & (table.kind == kind)]
    out = {}
    for name in dict.fromkeys(t.name):
        g = t[t.name == name].sort_values("step")
        out[name] = (g.x_weeks.to_numpy(), g.value.to_numpy(), g.ci_low.to_numpy(), g.ci_high.to_numpy())
    return out


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
METRIC_LABELS = {"pattern_acc": "Pattern ACC", "temporal_acc": "Area-averaged local temporal ACC",
                 "crpss": "CRPSS", "held_out_acc": "Held-out-member ACC"}


def _curves(ax, curves: dict, styles: dict, labels: dict | None = None, band_alpha: float = 0.15):
    for name, (x, v, lo, hi) in curves.items():
        st = dict(styles.get(name, {}))
        color = st.pop("color", None)
        ls = st.pop("linestyle", "-")
        line, = ax.plot(x, v, color=color, linestyle=ls, lw=st.pop("linewidth", 2.0),
                        label=(labels or {}).get(name, name))
        ax.fill_between(x, lo, hi, color=line.get_color(), alpha=band_alpha, lw=0)


def plot_strategy_figure(table: pd.DataFrame, acc_metric: str, exp_styles: dict, contrast_styles: dict,
                         contrast_labels: dict, *, title: str = "", xlabel: str = "Forecast lead (weeks)"):
    """(a) skill of every experiment, (b) ACC contrasts, (c) CRPSS contrasts; bands: bootstrap intervals."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6), constrained_layout=True)
    _curves(axes[0], from_table(table, acc_metric, "skill"), exp_styles, band_alpha=0.12)
    axes[0].set_title(f"(a) Forecast skill vs ERA5", loc="left", fontweight="bold")
    axes[0].set_ylabel(METRIC_LABELS[acc_metric])
    _curves(axes[1], from_table(table, acc_metric, "contrast"), contrast_styles, contrast_labels)
    axes[1].set_title("(b) Incremental skill of initialization strategies", loc="left", fontweight="bold")
    axes[1].set_ylabel(f"Δ {METRIC_LABELS[acc_metric]}")
    crpss = from_table(table, "crpss", "contrast")
    if crpss:
        _curves(axes[2], crpss, contrast_styles, contrast_labels)
    else:
        axes[2].text(0.5, 0.5, "not computed", ha="center", va="center", transform=axes[2].transAxes)
    axes[2].set_title("(c) Probabilistic skill differences", loc="left", fontweight="bold")
    axes[2].set_ylabel("Δ CRPSS (climatological reference)")
    for ax in axes:
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xlabel(xlabel)
        ax.grid(alpha=0.25)
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].legend(frameon=False, fontsize=8, loc="best")
    if title:
        fig.suptitle(title)
    return fig


def plot_agreement_figure(table: pd.DataFrame, exp_styles: dict, *, title: str = "",
                          xlabel: str = "Forecast lead (weeks)"):
    """Within-ensemble agreement: held-out-member ACC of every experiment, with bootstrap intervals."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.5, 4.4), constrained_layout=True)
    _curves(ax, from_table(table, "held_out_acc", "skill"), exp_styles, band_alpha=0.12)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Held-out-member ACC")
    ax.set_title("Within-ensemble agreement: held-out-member ACC", loc="left", fontweight="bold")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    if title:
        fig.suptitle(title)
    return fig


__all__ = ["draws", "season_strata", "pattern_stats", "pattern_acc_from_counts", "held_out_acc", "crps_stats",
           "crps_from_counts", "climatology_crps", "temporal_acc_from_counts", "strategy_scores", "scores_table",
           "from_table", "plot_strategy_figure", "plot_agreement_figure", "METRIC_LABELS"]
