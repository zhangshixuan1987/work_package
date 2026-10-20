"""Similarity evaluation and mosaic/cartopy spatial plotting for RTENSO2026 ICs.

Compares climatology-matched (``clim``) and observationally matched (``obs``)
restarts against the actual 2015 BruteForce restarts across key MPAS-Seaice,
ELM, and MOSART fields:

* Distance and rank checks against actual 2015 (``evaluate_distance_ranks``,
  ``plot_distance_vs_truth``).
* Field-level similarity metrics on lat-lon boxes and native model grids
  (``field_similarity_table``, ``plot_field_similarity_bars``,
  ``global_summary_table``).
* Spatial comparisons on lat-lon boxes (``plot_box_spatial``), on the native
  unstructured MPAS-Seaice mesh via ``mosaic`` (``build_mpassi_descriptors``,
  ``plot_mpassi_spatial``, ``plot_mpassi_four_seasons``), and on the native
  ELM ``grid1d`` land mesh via ``mosaic`` (``build_elm_descriptor``,
  ``plot_elm_spatial``, ``plot_elm_four_seasons``).
"""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

from clim_restart import (
    ALL_COMPONENTS,
    BRANCHES,
    ERA5_ROOT,
    KEY_FIELDS,
    OBS_FIELDS,
    SelectionPlan,
    _apply_box,
    _columns_to_gridcells,
    _to_boxes,
    comparison_points,
    compute_distances,
    distance_to_model_year,
    era5_available,
    load_key_fields,
    load_obs_fields,
    make_plan,
    select_years,
    source_file,
    summarize_restart,
)

FIELD_UNITS = {
    "NH ice concentration": "fraction",
    "SH ice concentration": "fraction",
    "NH ice thickness": "m",
    "SH ice thickness": "m",
    "NH snow on ice": "m",
    "SH snow on ice": "m",
    "snow water equivalent": "mm",
    "surface soil water": "kg m$^{-2}$",
    "subsurface soil water": "kg m$^{-2}$",
    "surface soil temperature": "K",
    "subsurface soil temperature": "K",
}

# Sequential (single-hue, light to dark) colormaps for each key field's state,
# and the diverging colormap for differences and standardized anomalies.
FIELD_CMAPS = {
    "NH ice concentration": "Blues",
    "SH ice concentration": "Blues",
    "NH ice thickness": "cmo.dense",
    "SH ice thickness": "cmo.dense",
    "NH snow on ice": "Blues",
    "SH snow on ice": "Blues",
    "snow water equivalent": "Blues",
    "surface soil water": "cmo.rain",
    "subsurface soil water": "cmo.rain",
    "surface soil temperature": "OrRd",
    "subsurface soil temperature": "OrRd",
}
DIFF_CMAP = "cmo.balance"
# Uneven bin edges for skewed fields (ice-sheet snow is ~1000 mm, seasonal snow ~10-200 mm).
FIELD_LEVELS = {
    "snow water equivalent": [0, 2, 5, 10, 25, 50, 100, 150, 200, 300, 500, 1000],
    "surface soil water": [0, 10, 20, 30, 40, 50, 60, 70, 80, 90],
    # Ice-sheet cells hold ~38,500 kg m-2 (ice-filled columns); other soils are ~400-3000.
    "subsurface soil water": [0, 250, 500, 750, 1000, 1250, 1500, 2000, 2500, 3000, 4000, 5000],
    "surface soil temperature": [240, 246, 252, 258, 264, 270, 276, 282, 288, 294, 300, 306, 312],
    "subsurface soil temperature": [240, 246, 252, 258, 264, 270, 276, 282, 288, 294, 300, 306, 312],
}
DIFF_PCT_LEVELS = {
    "snow water equivalent": [-100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100],
    "surface soil water": [-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50],
    "subsurface soil water": [-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50],
    "surface soil temperature": [-4, -3, -2, -1, 0, 1, 2, 3, 4],
    "subsurface soil temperature": [-4, -3, -2, -1, 0, 1, 2, 3, 4],
}
Z_LEVELS = [-3, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3]
NEUTRAL = "#f2f2f2"         # the two bins either side of zero on diverging scales
MISSING = {"mpassi": "0.78", "elm": "0.85"}   # background where there is no data


def discrete_cmap(cmap, lo: float, hi: float, data: np.ndarray | None = None, n: int = 10,
                  levels: Sequence[float] | None = None) -> dict:
    """``{"cmap", "norm"}`` for a discrete color scale with round bin edges over [lo, hi].

    ``levels`` gives explicit (possibly uneven) bin edges instead.

    A symmetric range (``lo == -hi``) is treated as diverging: zero is a bin edge
    and the two bins next to it are neutral gray. Out-of-range ``data`` (always,
    for diverging) adds end arrows.
    """
    import cmocean  # noqa: F401  (registers the cmo.* colormaps)
    from matplotlib import colormaps, ticker
    from matplotlib.colors import BoundaryNorm, ListedColormap, to_rgba

    diverging = bool(np.isclose(lo, -hi))
    edges = (np.asarray(levels, dtype=float) if levels is not None
             else ticker.MaxNLocator(nbins=n, symmetric=diverging).tick_values(lo, hi))
    if diverging:
        ext_lo = ext_hi = True
    else:
        finite = np.asarray(data)[np.isfinite(data)] if data is not None else np.array([])
        ext_lo = bool(finite.size and finite.min() < edges[0])
        ext_hi = bool(finite.size == 0 or finite.max() > edges[-1])
    extend = {(True, True): "both", (True, False): "min", (False, True): "max"}.get((ext_lo, ext_hi), "neither")
    nbins = len(edges) - 1
    total = nbins + ext_lo + ext_hi
    base = colormaps[cmap] if isinstance(cmap, str) else cmap
    colors = base(np.linspace(0.02, 0.98 if diverging else 0.9, total))
    if diverging:
        mid = ext_lo + nbins // 2
        colors[mid - 1] = colors[mid] = to_rgba(NEUTRAL)
    return {"cmap": ListedColormap(colors), "norm": BoundaryNorm(edges, total, extend=extend)}


BASE_FONTSIZE = 10.0


def font_sizes(fontsize: float = BASE_FONTSIZE) -> dict[str, float]:
    """Font sizes scaled from the panel-title size ``fontsize``."""
    return {
        "suptitle": 1.7 * fontsize,
        "title": 1.15 * fontsize,
        "label": 1.0 * fontsize,
        "tick": 0.88 * fontsize,
        "small": 0.8 * fontsize,
    }


def _panel(panel_size: tuple[float, float] | None, base: tuple[float, float], fontsize: float):
    """Panel (width, height) in inches: ``panel_size``, or ``base`` scaled with the font size."""
    if panel_size is not None:
        return panel_size
    k = fontsize / BASE_FONTSIZE
    return base[0] * k, base[1] * k


def _cbar(fig, mappable, fontsize: float = BASE_FONTSIZE, skip_ticks: bool = False, **kw):
    """Colorbar with customizable tick density and clean formatting."""
    from matplotlib.ticker import FuncFormatter
    f = font_sizes(fontsize)
    label = kw.pop("label", None)
    cb = fig.colorbar(mappable, **kw)
    edges = getattr(mappable.norm, "boundaries", None)
    if edges is not None:
        if skip_ticks and len(edges) > 8:
            cb.set_ticks(edges[::2])
        else:
            cb.set_ticks(edges)
            cb.ax.xaxis.set_major_formatter(
                FuncFormatter(lambda x, _: f"{int(round(x))}" if np.isclose(x, round(x)) else f"{x:g}")
            )
    cb.ax.tick_params(labelsize=f["tick"])
    if label:
        cb.set_label(label, fontsize=f["label"])
    return cb


def ensure_proj_data() -> None:
    """Ensure PROJ and GDAL data directories are set for Cartopy/mosaic."""
    proj_dir = Path(sys.prefix) / "share" / "proj"
    if proj_dir.is_dir():
        os.environ.setdefault("PROJ_DATA", str(proj_dir))
        os.environ.setdefault("PROJ_LIB", str(proj_dir))
    gdal_dir = Path(sys.prefix) / "share" / "gdal"
    if gdal_dir.is_dir():
        os.environ.setdefault("GDAL_DATA", str(gdal_dir))


def slugify(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def save_table(df: pd.DataFrame, name: str, diag_root: Path | None, enabled: bool = True) -> Path | None:
    if not enabled or diag_root is None:
        return None
    table_dir = Path(diag_root) / "data"
    table_dir.mkdir(parents=True, exist_ok=True)
    path = table_dir / f"{name}.csv"
    df.to_csv(path, index=False)
    print(f"[TABLE] {path}")
    return path


def save_fig(fig, name: str, diag_root: Path | None, enabled: bool = True, dpi: int = 180) -> Path | None:
    if not enabled or diag_root is None:
        return None
    fig_dir = Path(diag_root) / "figure"
    fig_dir.mkdir(parents=True, exist_ok=True)
    path = fig_dir / f"{name}.png"
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    print(f"[FIGURE] {path}")
    return path


# ---------------------------------------------------------------------------
# Selection & distance loading for test dates
# ---------------------------------------------------------------------------

def branches_in(table: Mapping[tuple[str, str], object]) -> list[str]:
    """Branch names in a ``{(branch, date): ...}`` mapping, in insertion order."""
    return list(dict.fromkeys(b for b, _ in table))


def _date_plan(plans: Mapping[tuple[str, str], SelectionPlan], date: str) -> SelectionPlan:
    """A plan for ``date`` (``clim`` if loaded); all branches share the candidate years and model fields."""
    return plans.get(("clim", date)) or next(p for (b, d), p in plans.items() if d == date)


def _has_obs(plan: SelectionPlan, component: str, year: int) -> bool:
    return bool(plan.obs_fields_for(component)) and era5_available(
        year, plan.mmdd, [component], obs_root=plan.obs_roots, obs_fields=plan.obs_fields
    )


def load_test_selections(
    test_dates: Sequence[str],
    src_root: Path,
    out_root: Path,
    *,
    era5_root: str | Path | Mapping[str, str | Path | None] | None = None,
    obs_root: str | Path | Mapping[str, str | Path | None] | None = None,
    obs_fields: Mapping[str, Mapping[str, str] | None] | None = None,
    obs_sets: Mapping[str, Mapping[str, object]] | None = None,
    branches: Sequence[str] | None = None,
    components: Sequence[str] = ALL_COMPONENTS,
    src_prefix: str = "20251024_s2d_spinup",
    out_prefix: str = "20251024_s2d_spinup",
    first_year: int | None = None,
    last_year: int | None = None,
    box: tuple[float, float] = (15.0, 30.0),
    box_shifts: int = 1,
    mode: str = "joint",
    component_weights: Mapping[str, float] | None = None,
    cache_dir: Path | None = None,
    log=lambda *_: None,
) -> tuple[dict, dict, dict]:
    """Load plans, distance tables, and selected years for ``test_dates``.

    ``obs_sets`` maps an obs branch name to ``{"root": ..., "fields": ...}``
    (as in ``02_build_climatology_ic.ipynb``); every branch other than ``clim``
    is an ``obs`` selection with that set. Without it, the ``obs`` branch uses
    ``obs_root``/``obs_fields``.
    """
    if obs_sets is None:
        obs_sets = {"obs": {"root": obs_root if obs_root is not None else era5_root, "fields": obs_fields}}
    if branches is None:
        branches = ["clim", *obs_sets]
    plans, distances, selections = {}, {}, {}
    for date in test_dates:
        for branch in branches:
            kind = "clim" if branch == "clim" else "obs"
            obs_set = obs_sets.get(branch, {})
            b_root, b_fields = obs_set.get("root"), obs_set.get("fields")
            if kind == "obs" and not era5_available(
                int(date[:4]), date[5:], components, obs_root=b_root, obs_fields=b_fields,
            ):
                continue
            plan = make_plan(
                date, src_root, out_root,
                obs_root=b_root, obs_fields=b_fields,
                branch=kind, components=components, label=branch,
                src_prefix=src_prefix, out_prefix=out_prefix, exclude_target_year=True,
                first_year=first_year, last_year=last_year,
            )
            table = compute_distances(
                plan, box=box, box_shifts=box_shifts, component_weights=component_weights,
                cache_dir=cache_dir, log=log,
            )
            key = (branch, date)
            plans[key], distances[key] = plan, table
            selections[key] = select_years(table, components, mode)
    return plans, distances, selections


def evaluate_distance_ranks(
    test_dates: Sequence[str],
    plans: Mapping[tuple[str, str], SelectionPlan],
    distances: Mapping[tuple[str, str], pd.DataFrame],
    *,
    target_year: int = 2015,
    box: tuple[float, float] = (15.0, 30.0),
    component_weights: Mapping[str, float] | None = None,
    cache_dir: Path | None = None,
    log=lambda *_: None,
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame], dict[str, pd.DataFrame]]:
    """Compare each branch's distances against the true distance to ``target_year``."""
    truth_obs, truth_all = {}, {}
    rows = []
    for date in test_dates:
        ref_plan = _date_plan(plans, date)
        truth_obs[date] = distance_to_model_year(
            ref_plan, target_year, box=box, observable_only=True,
            component_weights=component_weights, cache_dir=cache_dir, log=log,
        )
        truth_all[date] = distance_to_model_year(
            ref_plan, target_year, box=box, observable_only=False,
            component_weights=component_weights, cache_dir=cache_dir, log=log,
        )
        rank_obs = truth_obs[date]["joint"].rank()
        rank_all = truth_all[date]["joint"].rank()
        row = {
            "target_date": date,
            f"closest to {target_year} (obs fields)": int(truth_obs[date]["joint"].idxmin()),
            "min D_true (obs fields)": float(truth_obs[date]["joint"].min()),
            f"closest to {target_year} (all fields)": int(truth_all[date]["joint"].idxmin()),
            "min D_true (all fields)": float(truth_all[date]["joint"].min()),
        }
        for branch in branches_in(distances):
            if (branch, date) not in distances:
                continue
            year = int(distances[(branch, date)]["joint"].idxmin())
            row[f"{branch} year"] = year
            row[f"{branch}: D_true (obs fields)"] = float(truth_obs[date].loc[year, "joint"])
            row[f"{branch}: rank (obs fields)"] = int(rank_obs[year])
            row[f"{branch}: rank corr (obs)"] = float(distances[(branch, date)]["joint"].rank().corr(rank_obs))
            row[f"{branch}: D_true (all fields)"] = float(truth_all[date].loc[year, "joint"])
            row[f"{branch}: rank (all fields)"] = int(rank_all[year])
        rows.append(row)
    return pd.DataFrame(rows), truth_obs, truth_all


def plot_distance_vs_truth(
    test_dates: Sequence[str],
    distances: Mapping[tuple[str, str], pd.DataFrame],
    truth_obs: Mapping[str, pd.DataFrame],
    *,
    branches: Sequence[str] | None = None,
    fontsize: float = BASE_FONTSIZE,
    panel_size: tuple[float, float] | None = None,
):
    """Scatter plot of candidate-year branch distance vs true distance to 2015."""
    import matplotlib.pyplot as plt

    if branches is None:
        branches = branches_in(distances)
    f = font_sizes(fontsize)
    pw, ph = _panel(panel_size, (4.2, 3.6), fontsize)
    fig, axes = plt.subplots(
        len(branches), len(test_dates),
        figsize=(pw * len(test_dates), ph * len(branches)),
        sharey="row", squeeze=False,
    )
    for col, date in enumerate(test_dates):
        t_dist = truth_obs[date]["joint"]
        best_true_yr = int(t_dist.idxmin())
        for row, branch in enumerate(branches):
            ax = axes[row, col]
            if (branch, date) not in distances:
                ax.axis("off")
                continue
            b_dist = distances[(branch, date)]["joint"]
            sel_yr = int(b_dist.idxmin())
            other_yrs = [y for y in b_dist.index if y not in (sel_yr, best_true_yr)]

            ax.scatter(b_dist.loc[other_yrs], t_dist.loc[other_yrs], c="0.65", s=28, alpha=0.8, label="Candidates")
            ax.scatter(b_dist.loc[best_true_yr], t_dist.loc[best_true_yr], c="tab:green", marker="^", s=80, zorder=4,
                       label=f"Closest to 2015 ({best_true_yr})")
            ax.scatter(b_dist.loc[sel_yr], t_dist.loc[sel_yr], c="tab:red", s=85, zorder=5,
                       label=f"Selected {branch} ({sel_yr})")
            for y in (sel_yr, best_true_yr):
                ax.annotate(str(y), (b_dist.loc[y], t_dist.loc[y]),
                            textcoords="offset points", xytext=(5, 4), fontsize=f["small"])

            r_spearman = b_dist.rank().corr(t_dist.rank())
            true_rk = int(t_dist.rank()[sel_yr])
            ax.set_title(f"{branch} {date}\nsel={sel_yr} (rank {true_rk}/{len(b_dist)}), rank corr={r_spearman:.2f}",
                         fontsize=f["title"])
            if row == len(branches) - 1:
                ax.set_xlabel(f"{branch} distance $D$", fontsize=f["label"])
            ax.tick_params(labelsize=f["tick"])
            ax.grid(True, alpha=0.3)
        for row in range(len(branches)):
            axes[row, 0].set_ylabel(r"True distance to 2015 ($D_{\mathrm{true}}$, obs fields)", fontsize=f["label"])

    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="0.65", markersize=7, label="Candidates"),
        plt.Line2D([0], [0], marker="^", color="w", markerfacecolor="tab:green", markersize=9, label="Closest to 2015 truth"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="tab:red", markersize=9, label="Selected candidate year"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.985), ncol=3,
               fontsize=f["label"], frameon=True, facecolor="white", edgecolor="0.8")
    fig.suptitle("Candidate-year selection distance vs. true distance to actual 2015 IC",
                 fontsize=f["suptitle"], fontweight="bold", y=1.02)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    return fig


# ---------------------------------------------------------------------------
# Field-level similarity metrics and global summaries
# ---------------------------------------------------------------------------

def _weighted_corr(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> float:
    """Area-weighted anomaly pattern correlation between ``x`` and ``y``."""
    valid = np.isfinite(x) & np.isfinite(y) & (w > 0)
    if not np.any(valid):
        return float("nan")
    ww = w[valid] / w[valid].sum()
    dx = x[valid] - np.sum(ww * x[valid])
    dy = y[valid] - np.sum(ww * y[valid])
    denom = float(np.sqrt(np.sum(ww * dx ** 2) * np.sum(ww * dy ** 2)))
    return float(np.sum(ww * dx * dy) / denom) if denom > 0 else float("nan")


def _weighted_rms_diff(a: np.ndarray, b: np.ndarray, w: np.ndarray) -> float:
    """Area-weighted root-mean-square difference between ``a`` and ``b``."""
    valid = np.isfinite(a) & np.isfinite(b) & (w > 0)
    if not np.any(valid):
        return float("nan")
    ww = w[valid]
    return float(np.sqrt(np.sum(ww * (a[valid] - b[valid]) ** 2) / ww.sum()))


def _standardize_series(
    values: np.ndarray,
    weight: np.ndarray,
    target: np.ndarray | None = None,
    valid_mask: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray | None, np.ndarray, float, np.ndarray]:
    """Standardize candidate-year ``values`` (n_years, n_pts) and optional ``target`` (n_pts,).

    Returns ``(z_candidates, z_target, clim_mean, spread, valid_weight)``.
    """
    valid = np.all(np.isfinite(values), axis=0) & (weight > 0)
    if target is not None:
        valid &= np.isfinite(target)
    if valid_mask is not None:
        valid &= valid_mask
    w = np.where(valid, weight, 0.0)
    wsum = float(w.sum())
    v = np.where(valid, values, 0.0)
    mean = v.mean(axis=0)
    anom = v - mean
    spread = float(np.sqrt((w * (anom ** 2).mean(axis=0)).sum() / wsum)) if wsum > 0 else 0.0
    if spread <= 0:
        z = np.full_like(values, np.nan, dtype=np.float64)
        zt = None if target is None else np.full_like(target, np.nan, dtype=np.float64)
    else:
        z = np.where(valid, anom / spread, np.nan)
        zt = None if target is None else np.where(valid, (target - mean) / spread, np.nan)
    return z, zt, np.where(valid, mean, np.nan), spread, w


def field_similarity_table(
    plan: SelectionPlan,
    target_year: int,
    branch_years: Mapping[str, int],
    *,
    obs_plans: Mapping[str, SelectionPlan] | None = None,
    box: tuple[float, float] = (15.0, 30.0),
    cache_dir: Path | None = None,
    log=lambda *_: None,
) -> pd.DataFrame:
    """Detailed similarity metrics comparing selected years (and observations) against ``target_year``.

    ``obs_plans`` maps a ``source`` label to a plan whose observation set is
    compared with the actual anomaly (default ``{"ERA5": plan}``).
    """
    if obs_plans is None:
        obs_plans = {"ERA5": plan}
    cand_years = list(plan.candidate_years)
    year_to_idx = {y: i for i, y in enumerate(cand_years)}
    rows = []

    for comp in plan.components:
        if comp not in KEY_FIELDS:
            continue
        fields = load_key_fields(plan, comp, cache_dir, years=cand_years + [target_year], log=log)
        pts = comparison_points(plan, comp, cache_dir)
        obs_by_source = {
            label: load_obs_fields(op, comp, cand_years + [target_year], cache_dir, log=log)
            for label, op in obs_plans.items() if _has_obs(op, comp, target_year)
        }

        for fname, (vals_all, weight_raw) in fields.items():
            if comp == "elm":
                vals_grid, w_grid = _columns_to_gridcells(
                    vals_all, weight_raw, pts["col_grid"], pts["lat"].size
                )
            else:
                vals_grid, w_grid = vals_all, weight_raw
            lat, lon = pts["lat"], pts["lon"]

            cand_grid, targ_grid = vals_grid[:-1], vals_grid[-1]
            z_grid, zt_grid, _, _, wg_valid = _standardize_series(cand_grid, w_grid, targ_grid)

            box_vals, w_box, _ = _to_boxes(vals_grid, w_grid, lat, lon, box)
            z_box, zt_box, _, _, wb_valid = _standardize_series(box_vals[:-1], w_box, box_vals[-1])

            for branch, y in branch_years.items():
                if y not in year_to_idx:
                    continue
                idx = year_to_idx[y]
                diff = np.where(wg_valid > 0, cand_grid[idx] - targ_grid, 0.0)
                wsum_g = float(wg_valid.sum())
                rows.append({
                    "target_date": plan.target_date,
                    "component": comp,
                    "field": fname,
                    "units": FIELD_UNITS.get(fname, ""),
                    "source": branch,
                    "year": int(y),
                    "d_box": _weighted_rms_diff(z_box[idx], zt_box, wb_valid),
                    "r_box": _weighted_corr(z_box[idx], zt_box, wb_valid),
                    "d_grid": _weighted_rms_diff(z_grid[idx], zt_grid, wg_valid),
                    "r_grid": _weighted_corr(z_grid[idx], zt_grid, wg_valid),
                    "RMSE": _weighted_rms_diff(cand_grid[idx], targ_grid, wg_valid),
                    "PCC": _weighted_corr(cand_grid[idx], targ_grid, wg_valid),
                    "bias": float((wg_valid * diff).sum() / wsum_g) if wsum_g > 0 else np.nan,
                })

            for label, obs in obs_by_source.items():
                if fname not in obs:
                    continue
                o_all = obs[fname]
                valid_both = np.all(np.isfinite(vals_grid), axis=0) & np.all(np.isfinite(o_all), axis=0)
                _, zt_g_obs, _, _, wg_o = _standardize_series(
                    o_all[:-1], w_grid, o_all[-1], valid_mask=valid_both
                )
                _, zt_g_mod, _, _, _ = _standardize_series(
                    cand_grid, w_grid, targ_grid, valid_mask=valid_both
                )
                m_box, wb_o, o_box = _apply_box(vals_grid, w_grid, lat, lon, box, obs=o_all)
                _, zt_b_mod, _, _, wb_ov = _standardize_series(m_box[:-1], wb_o, m_box[-1])
                _, zt_b_obs, _, _, _ = _standardize_series(o_box[:-1], wb_o, o_box[-1])

                rows.append({
                    "target_date": plan.target_date,
                    "component": comp,
                    "field": fname,
                    "units": FIELD_UNITS.get(fname, ""),
                    "source": label,
                    "year": int(target_year),
                    "d_box": _weighted_rms_diff(zt_b_obs, zt_b_mod, wb_ov),
                    "r_box": _weighted_corr(zt_b_obs, zt_b_mod, wb_ov),
                    "d_grid": _weighted_rms_diff(zt_g_obs, zt_g_mod, wg_o),
                    "r_grid": _weighted_corr(zt_g_obs, zt_g_mod, wg_o),
                    "RMSE": np.nan,
                    "PCC": np.nan,
                    "bias": np.nan,
                })

    return pd.DataFrame(rows)


def plot_field_similarity_bars(
    sim_df: pd.DataFrame,
    *,
    fontsize: float = BASE_FONTSIZE,
    panel_size: tuple[float, float] | None = None,
):
    """Horizontal bar chart of four-season mean d_box and r_box by key field."""
    import matplotlib.pyplot as plt

    mean_sim = (
        sim_df.groupby(["component", "field", "source"], as_index=False)[["d_box", "r_box", "d_grid", "r_grid", "RMSE"]]
        .mean()
    )
    # Branches (clim first, then obs*) followed by the observation references.
    found = list(dict.fromkeys(sim_df["source"]))
    branch_src = [b for b in found if b == "clim" or b.startswith("obs")]
    sources = sorted(branch_src, key=lambda b: b != "clim") + [b for b in found if b not in branch_src]
    pivot_d = mean_sim.pivot(index=["component", "field"], columns="source", values="d_box")[sources]
    pivot_r = mean_sim.pivot(index=["component", "field"], columns="source", values="r_box")[sources]

    f = font_sizes(fontsize)
    pw, ph = _panel(panel_size, (6.5, 0.6 * len(pivot_d) * max(1, len(sources) / 3) + 1.0), fontsize)
    fig, axes = plt.subplots(1, 2, figsize=(2 * pw, ph), sharey=True, constrained_layout=True)
    labels = [f"{c}: {f}" for c, f in pivot_d.index]
    y_pos = np.arange(len(labels))
    width = 0.8 / len(sources)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for ax, piv, xlabel, vline in (
        (axes[0], pivot_d, r"Mean box distance $d_{\mathrm{box}}$ to actual 2015 ($\leftarrow$ smaller is closer)", 1.0),
        (axes[1], pivot_r, r"Mean box anomaly correlation $r_{\mathrm{box}}$ with actual 2015 (larger is closer $\rightarrow$)", 0.0),
    ):
        for k, src in enumerate(sources):
            label = f"{src} IC" if src in branch_src else f"{src} target anomaly"
            ax.barh(y_pos + (k - (len(sources) - 1) / 2) * width, piv[src], height=width,
                    label=label, color=colors[k % len(colors)], alpha=0.85 if src in branch_src else 0.7)
        ax.axvline(vline, color="k", lw=0.8, ls="--")
        ax.set_yticks(y_pos, labels)
        ax.invert_yaxis()
        ax.set_xlabel(xlabel, fontsize=f["label"])
        ax.tick_params(labelsize=f["tick"])
        ax.grid(True, axis="x", alpha=0.3)

    handles, leg_labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, leg_labels, loc="upper center", bbox_to_anchor=(0.5, 1.03),
               ncol=min(len(sources), 4), fontsize=f["label"], frameon=True, facecolor="white", edgecolor="0.8")
    fig.suptitle("Four-season mean similarity to actual 2015 IC by field (2015 sensitivity dates)",
                 fontsize=f["suptitle"], fontweight="bold", y=1.07)
    return fig


def global_summary_table(
    test_dates: Sequence[str],
    selections: Mapping[tuple[str, str], Mapping[str, int]],
    src_root: Path,
    src_prefix: str = "20251024_s2d_spinup",
    components: Sequence[str] = ALL_COMPONENTS,
    branches: Sequence[str] | None = None,
    target_year: int = 2015,
) -> pd.DataFrame:
    """Global integrated quantities comparing actual ``target_year`` against each branch."""
    if branches is None:
        branches = branches_in(selections)
    rows = []
    for date in test_dates:
        for comp in components:
            actual = summarize_restart(source_file(src_root, src_prefix, comp, target_year, date[5:]), comp)
            chosen = {}
            for branch in branches:
                if (branch, date) in selections:
                    year = selections[(branch, date)][comp]
                    chosen[branch] = (
                        year,
                        summarize_restart(source_file(src_root, src_prefix, comp, year, date[5:]), comp),
                    )
            for name in actual:
                row = {"target_date": date, "component": comp, "quantity": name, f"actual_{target_year}": actual[name]}
                for branch, (year, values) in chosen.items():
                    row[f"{branch} year"] = year
                    row[branch] = values[name]
                    row[f"{branch} - actual"] = values[name] - actual[name]
                rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Spatial field preparation and plotting (lat-lon boxes, MPAS-Seaice, ELM)
# ---------------------------------------------------------------------------

def prepare_spatial_field(
    plan: SelectionPlan,
    component: str,
    field_name: str,
    target_year: int,
    branch_years: Mapping[str, int],
    *,
    obs_plans: Mapping[str, SelectionPlan] | None = None,
    box: tuple[float, float] = (15.0, 30.0),
    cache_dir: Path | None = None,
    log=lambda *_: None,
) -> dict:
    """Prepare native-grid and box-averaged fields, anomalies, and differences for plotting.

    ``obs_plans`` maps a key to a plan whose observed anomaly is added under that
    key (default ``{"era5": plan}``); the keys are listed in ``"obs_refs"``.
    """
    if obs_plans is None:
        obs_plans = {"era5": plan}
    cand_years = list(plan.candidate_years)
    year_to_idx = {y: i for i, y in enumerate(cand_years)}
    fields = load_key_fields(plan, component, cache_dir, years=cand_years + [target_year], log=log)
    vals_all, weight_raw = fields[field_name]
    pts = comparison_points(plan, component, cache_dir)

    if component == "elm":
        vals_grid, w_grid = _columns_to_gridcells(
            vals_all, weight_raw, pts["col_grid"], pts["lat"].size
        )
    else:
        vals_grid, w_grid = vals_all, weight_raw
    lat, lon = pts["lat"], pts["lon"]

    cand_grid, targ_grid = vals_grid[:-1], vals_grid[-1]
    z_grid, zt_grid, clim_mean, spread, wg_valid = _standardize_series(cand_grid, w_grid, targ_grid)

    dlat, dlon = box
    nlat, nlon = int(round(180.0 / dlat)), int(round(360.0 / dlon))
    lat_edges = np.linspace(-90.0, 90.0, nlat + 1)
    lon_edges = np.linspace(0.0, 360.0, nlon + 1)

    box_vals, w_box, _ = _to_boxes(vals_grid, w_grid, lat, lon, box)
    z_box, zt_box, _, _, wb_valid = _standardize_series(box_vals[:-1], w_box, box_vals[-1])

    raw = {
        "actual": np.where(wg_valid > 0, targ_grid, np.nan),
        "climatology": clim_mean,
    }
    diff = {}
    z_map = {"actual": zt_grid}
    box_z = {"actual": np.where(wb_valid > 0, zt_box, np.nan).reshape(nlat, nlon)}
    box_diff_z = {}
    metrics = {}

    wsum_g = float(wg_valid.sum())

    for branch, y in branch_years.items():
        idx = year_to_idx[y]
        raw_b = np.where(wg_valid > 0, cand_grid[idx], np.nan)
        raw[branch] = raw_b
        diff[branch] = raw_b - raw["actual"]
        z_map[branch] = z_grid[idx]
        zb = np.where(wb_valid > 0, z_box[idx], np.nan)
        box_z[branch] = zb.reshape(nlat, nlon)
        box_diff_z[branch] = (zb - np.where(wb_valid > 0, zt_box, np.nan)).reshape(nlat, nlon)

        df_diff = np.where(wg_valid > 0, cand_grid[idx] - targ_grid, 0.0)
        metrics[branch] = {
            "year": int(y),
            "d_box": _weighted_rms_diff(z_box[idx], zt_box, wb_valid),
            "r_box": _weighted_corr(z_box[idx], zt_box, wb_valid),
            "d_grid": _weighted_rms_diff(z_grid[idx], zt_grid, wg_valid),
            "r_grid": _weighted_corr(z_grid[idx], zt_grid, wg_valid),
            "rmse": _weighted_rms_diff(cand_grid[idx], targ_grid, wg_valid),
            "pcc": _weighted_corr(cand_grid[idx], targ_grid, wg_valid),
            "bias": float((wg_valid * df_diff).sum() / wsum_g) if wsum_g > 0 else np.nan,
        }

    # Baseline: the candidate-year climatology against actual.
    clim_diff = np.where(wg_valid > 0, clim_mean - targ_grid, 0.0)
    metrics["climatology"] = {
        "rmse": _weighted_rms_diff(clim_mean, targ_grid, wg_valid),
        "pcc": _weighted_corr(clim_mean, targ_grid, wg_valid),
        "bias": float((wg_valid * clim_diff).sum() / wsum_g) if wsum_g > 0 else np.nan,
    }

    obs_refs = []
    vs_obs = {}     # {ref: {"actual" | branch: metrics against the observed anomaly}}
    for ref, op in obs_plans.items():
        if field_name not in op.obs_fields_for(component) or not _has_obs(op, component, target_year):
            continue
        obs_refs.append(ref)
        obs_all = load_obs_fields(op, component, cand_years + [target_year], cache_dir, log=log)[field_name]
        valid_both = np.all(np.isfinite(vals_grid), axis=0) & np.all(np.isfinite(obs_all), axis=0)
        _, zt_g_obs, _, _, wg_o = _standardize_series(obs_all[:-1], w_grid, obs_all[-1], valid_mask=valid_both)
        z_g_mod, zt_g_mod, _, _, _ = _standardize_series(cand_grid, w_grid, targ_grid, valid_mask=valid_both)
        z_map[ref] = zt_g_obs

        m_box, wb_o, o_box = _apply_box(vals_grid, w_grid, lat, lon, box, obs=obs_all)
        z_b_mod, zt_b_mod, _, _, wb_ov = _standardize_series(m_box[:-1], wb_o, m_box[-1])
        _, zt_b_obs, _, _, _ = _standardize_series(o_box[:-1], wb_o, o_box[-1])
        box_z[ref] = np.where(wb_ov > 0, zt_b_obs, np.nan).reshape(nlat, nlon)
        box_diff_z[ref] = np.where(wb_ov > 0, zt_b_obs - zt_b_mod, np.nan).reshape(nlat, nlon)

        metrics[ref] = {
            "year": int(target_year),
            "d_box": _weighted_rms_diff(zt_b_obs, zt_b_mod, wb_ov),
            "r_box": _weighted_corr(zt_b_obs, zt_b_mod, wb_ov),
            "d_grid": _weighted_rms_diff(zt_g_obs, zt_g_mod, wg_o),
            "r_grid": _weighted_corr(zt_g_obs, zt_g_mod, wg_o),
            "rmse": np.nan,
            "bias": np.nan,
        }

        def _vs(zg, zb):
            return {"d_box": _weighted_rms_diff(zb, zt_b_obs, wb_ov), "r_box": _weighted_corr(zb, zt_b_obs, wb_ov),
                    "d_grid": _weighted_rms_diff(zg, zt_g_obs, wg_o), "r_grid": _weighted_corr(zg, zt_g_obs, wg_o)}

        vs_obs[ref] = {"actual": _vs(zt_g_mod, zt_b_mod)}
        for branch, y in branch_years.items():
            vs_obs[ref][branch] = _vs(z_g_mod[year_to_idx[y]], z_b_mod[year_to_idx[y]])

    return {
        "component": component,
        "field_name": field_name,
        "units": FIELD_UNITS.get(field_name, ""),
        "target_date": plan.target_date,
        "target_year": int(target_year),
        "clim_years": (int(cand_years[0]), int(cand_years[-1])),
        "branch_years": {k: int(v) for k, v in branch_years.items()},
        "obs_refs": obs_refs,
        "lat": lat,
        "lon": lon,
        "weight": wg_valid,
        "clim_mean": clim_mean,
        "spread": spread,
        "raw": raw,
        "diff": diff,
        "z": z_map,
        "metrics": metrics,
        "vs_obs": vs_obs,
        "box": {
            "dlat": dlat,
            "dlon": dlon,
            "lat_edges": lat_edges,
            "lon_edges": lon_edges,
            "z": box_z,
            "diff_z": box_diff_z,
        },
    }


def plot_box_spatial(
    spatial: dict,
    *,
    z_lim: float = 2.0,
    fontsize: float = BASE_FONTSIZE,
    panel_size: tuple[float, float] | None = None,
):
    """Plot box-averaged standardized anomalies used in the distance metric."""
    ensure_proj_data()
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import matplotlib.path as mpath
    import matplotlib.pyplot as plt

    fname = spatial["field_name"]
    comp = spatial.get("component", "elm")
    tyear = spatial["target_year"]
    byears = spatial["branch_years"]
    metrics = spatial["metrics"]
    binfo = spatial["box"]
    lon_e, lat_e = binfo["lon_edges"], binfo["lat_edges"]

    date = spatial["target_date"]
    refs = [r for r in spatial.get("obs_refs", ["era5"]) if r in binfo["z"]]
    if refs:     # the observation is the truth; statistics of the model panels are against it
        metrics = spatial["vs_obs"][refs[0]]
        m = metrics["actual"]
        panels = [(_obs_title(refs[0], date, spatial.get("component")), binfo["z"][refs[0]]),
                  (f"{_actual_title(date)}\nvs truth: box PCC={m['r_box']:.2f}, RMSE={m['d_box']:.2f}",
                   binfo["z"]["actual"])]
    else:
        panels = [(f"Truth: E3SM actual {date} state", binfo["z"]["actual"])]
    for b in byears:
        if b in binfo["z"]:
            m = metrics[b]
            panels.append((f"{_branch_title(b, byears[b], date)}\n"
                           f"vs truth: box PCC={m['r_box']:.2f}, RMSE={m['d_box']:.2f}",
                           binfo["z"][b]))

    f = font_sizes(fontsize)
    pw, ph = _panel(panel_size, (4.4, 3.6), fontsize)   # room for three-line panel titles
    ncol = min(len(panels), 4)
    nrow = -(-len(panels) // ncol)

    circle_path = None
    if comp == "mpassi":
        if fname.startswith("NH"):
            proj = ccrs.NorthPolarStereo()
            extent = [-180, 180, 50, 90]
        else:
            proj = ccrs.SouthPolarStereo()
            extent = [-180, 180, -90, -50]
        theta = np.linspace(0, 2 * np.pi, 100)
        circle_path = mpath.Path(np.vstack([np.sin(theta), np.cos(theta)]).T * 0.5 + 0.5)
    else:
        proj = ccrs.Robinson()
        extent = [-180, 180, -60, 85]

    fig, axes = plt.subplots(
        nrow, ncol, figsize=(pw * ncol, ph * nrow),
        subplot_kw={"projection": proj},
        constrained_layout=True, squeeze=False,
    )
    for ax in axes.flat[len(panels):]:
        ax.axis("off")
    pm = None
    for ax, (title, grid_z) in zip(axes.flat, panels):
        ax.set_facecolor("#f5f7fa")
        ax.set_extent(extent, crs=ccrs.PlateCarree())
        if circle_path is not None:
            ax.set_boundary(circle_path, transform=ax.transAxes)
        ax.add_feature(cfeature.LAND, facecolor="#eceff1", zorder=0)
        pm = ax.pcolormesh(lon_e, lat_e, grid_z, transform=ccrs.PlateCarree(),
                           **discrete_cmap(DIFF_CMAP, -z_lim, z_lim), edgecolors="0.65", linewidth=0.25, zorder=2)
        ax.coastlines(linewidth=0.45, color="0.3", zorder=3)
        ax.set_title(title, fontsize=f["title"])
    _cbar(fig, pm, fontsize, ax=axes.ravel().tolist(), orientation="horizontal",
          fraction=0.045, pad=0.04, label=f"{int(binfo['dlat'])}°x{int(binfo['dlon'])}° box standardized anomaly z")
    fig.suptitle(f"Box-scale ({int(binfo['dlat'])}°x{int(binfo['dlon'])}°) anomalies: {comp} {fname} — {spatial['target_date']}",
                 fontsize=f["suptitle"], fontweight="bold")
    return fig


def build_mpassi_descriptors(mesh_path: Path) -> dict[str, object]:
    """Create ``mosaic.Descriptor`` objects for NH and SH polar stereographic projections."""
    ensure_proj_data()
    import cartopy.crs as ccrs
    import mosaic
    import xarray as xr

    ds = xr.open_dataset(mesh_path)
    transform = ccrs.PlateCarree()
    descriptors = {
        "NH": mosaic.Descriptor(ds, ccrs.NorthPolarStereo(), transform, use_latlon=True),
        "SH": mosaic.Descriptor(ds, ccrs.SouthPolarStereo(), transform, use_latlon=True),
    }
    for desc in descriptors.values():
        _ = desc.cell_patches
    return descriptors


def _ref_name(ref: str) -> str:
    return "ERA5" if ref == "era5" else ref


# Observation dataset behind each obs-set label, per component.
OBS_DATASETS = {"ERA5-Land/NSIDC": {"elm": "ERA5-Land", "mpassi": "NOAA/NSIDC CDR"}}


def _obs_title(ref: str, date: str, component: str | None = None) -> str:
    name = _ref_name(ref)
    dataset = OBS_DATASETS.get(name, {}).get(component, name)
    return f"Truth: observed {date} anomaly\n({dataset})"


def _actual_title(date: str) -> str:
    return f"Reference: E3SM actual {date} state"


def _branch_title(branch: str, year: int, date: str | None = None) -> str:
    when = f"{year}{date[4:]}" if date else str(year)
    return f"Selected initial condition: E3SM {when} state\n(year chosen by matching the observation)"


def _state_stats(m: dict) -> str:
    """PCC and RMSE of a field against the actual (reference) field."""
    return f"vs reference: PCC={m['pcc']:.2f}, RMSE={m['rmse']:.3g}"


def _state_titles(sp: dict, branch: str, year: int) -> dict:
    """Titles of the state-and-difference panels."""
    date = sp["target_date"]
    y0, y1 = sp.get("clim_years", ("", ""))
    return {"actual": f"Reference: E3SM actual {date} state",
            "branch": f"Selected initial condition: E3SM {year}{date[4:]} state",
            "clim": f"Baseline: E3SM {y0}-{y1} mean for {date[5:]}\n(no year selection)",
            "diff": "Error: selected initial condition - reference"}


def _z_stats(m: dict) -> str:
    """PCC and RMSE of a standardized anomaly against the actual one (grid), plus the box PCC."""
    return f"vs truth: PCC={m['r_grid']:.2f}, RMSE={m['d_grid']:.2f}, box PCC={m['r_box']:.2f}"


def _z_panels(spatial: dict, **_) -> list[tuple[str, np.ndarray]]:
    """``(title, z)`` panels: the observed anomaly (the truth), the actual model state (the reference), then
    each selected initial condition; the model panels' statistics are against the observed anomaly."""
    date, byears, metrics = spatial["target_date"], spatial["branch_years"], spatial["metrics"]
    comp = spatial.get("component")
    refs = [r for r in spatial.get("obs_refs", ["era5"]) if r in spatial["z"]]
    if refs:     # the observation is the truth; statistics of the model panels are against it
        truth = spatial["vs_obs"][refs[0]]
        panels = [(_obs_title(refs[0], date, comp), spatial["z"][refs[0]])]
        panels.append((f"{_actual_title(date)}\n{_z_stats(truth['actual'])}", spatial["z"]["actual"]))
        panels += [(f"{_branch_title(b, byears[b], date)}\n{_z_stats(truth[b])}", spatial["z"][b])
                   for b in byears if b in spatial["z"]]
    else:        # no observation for this field: compare with the actual state
        panels = [(f"Truth: E3SM actual {date} state", spatial["z"]["actual"])]
        panels += [(f"{_branch_title(b, byears[b], date)}\n{_z_stats(metrics[b])}", spatial["z"][b])
                   for b in byears if b in spatial["z"]]
    return panels


def _diff_limit(spatial: dict, pct: float, default: float) -> float:
    diffs = [d[np.isfinite(d)] for d in spatial["diff"].values()]
    all_d = np.concatenate(diffs) if diffs else np.array([default])
    return max(float(np.nanpercentile(np.abs(all_d), pct)) if all_d.size else default, 1e-3)


def _compute_diff_percent(fn: str, raw_b: np.ndarray, act: np.ndarray) -> np.ndarray:
    """Calculate percent difference 100 * (selected - reference) / reference safely."""
    diff = raw_b - act
    if "temperature" in fn:
        # In Kelvin (~240-315 K): reference is safely > 0
        with np.errstate(divide="ignore", invalid="ignore"):
            pct = 100.0 * diff / np.where(act > 0, act, np.nan)
            pct = np.clip(pct, -10.0, 10.0)
    elif "snow" in fn:
        # Snow water equivalent (mm)
        with np.errstate(divide="ignore", invalid="ignore"):
            pct = np.where(act >= 1.0, 100.0 * diff / act, np.nan)
            # Both act and candidate have no snow (<1 mm): 0% difference
            pct = np.where((act < 1.0) & (raw_b < 1.0), 0.0, pct)
            # Candidate has snow (>=1 mm) but reference had none (<1 mm): +100%
            pct = np.where((act < 1.0) & (raw_b >= 1.0), 100.0, pct)
            pct = np.clip(pct, -100.0, 100.0)
    else:
        # Soil moisture (kg m-2)
        with np.errstate(divide="ignore", invalid="ignore"):
            pct = np.where(act > 0.1, 100.0 * diff / act, np.nan)
            pct = np.where((act <= 0.1) & (raw_b <= 0.1), 0.0, pct)
            pct = np.where((act <= 0.1) & (raw_b > 0.1), 100.0, pct)
            pct = np.clip(pct, -100.0, 100.0)
    # Mask non-land or missing cells
    pct = np.where(np.isfinite(raw_b) & np.isfinite(act), pct, np.nan)
    return pct


def _obs_plans_for(
    plans: Mapping[tuple[str, str], SelectionPlan],
    date: str,
    obs_refs: Mapping[str, str] | None,
) -> dict[str, SelectionPlan] | None:
    """``{label: plan}`` of the obs branches named in ``obs_refs`` (``{branch: label}``) on ``date``."""
    if obs_refs is None:
        return None
    return {label: plans[(b, date)] for b, label in obs_refs.items() if (b, date) in plans}


def plot_mpassi_spatial(
    spatial: dict | Sequence[dict],
    descriptors: Mapping[str, object],
    *,
    kind: str = "state_and_diff",
    cmap: str | None = None,
    diff_cmap: str = DIFF_CMAP,
    vmin: float | None = None,
    vmax: float | None = None,
    diff_lim: float | None = None,
    z_lim: float = 2.5,
    fontsize: float = BASE_FONTSIZE,
    panel_size: tuple[float, float] | None = None,
):
    """Plot MPAS-Seaice fields on the native unstructured MPAS mesh using ``mosaic``.

    ``spatial`` is one ``prepare_spatial_field`` result or several (e.g. the NH
    and SH versions of a field), drawn one row each in their own polar
    projection. ``state_and_diff``: actual, each branch, and each branch minus
    actual. ``z``: actual, each observation reference, and each branch.
    """
    ensure_proj_data()
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import matplotlib.path as mpath
    import matplotlib.pyplot as plt
    import mosaic

    theta = np.linspace(0, 2 * np.pi, 100)
    circle_path = mpath.Path(np.vstack([np.sin(theta), np.cos(theta)]).T * 0.5 + 0.5)

    spatials = [spatial] if isinstance(spatial, Mapping) else list(spatial)
    rows = []   # (descriptor, extent, [(title, arr, cmap, lo, hi, colorbar label)])
    for sp in spatials:
        fname = sp["field_name"]
        hemi = "NH" if fname.startswith("NH") else "SH"
        prefix = f"{hemi}: " if len(spatials) > 1 else ""
        extent = [-180, 180, 45, 90] if hemi == "NH" else [-180, 180, -90, -45]
        tyear, units = sp["target_year"], sp["units"]
        byears, metrics = sp["branch_years"], sp["metrics"]
        branches = [b for b in byears if b in sp["raw"]]
        state_cmap = cmap or FIELD_CMAPS.get(fname, "cmo.ice")

        if kind == "state_and_diff":
            act = sp["raw"]["actual"]
            lo = 0.0 if vmin is None else vmin
            hi = vmax
            if hi is None:
                hi = 1.0 if "concentration" in fname else (
                    float(np.nanpercentile(act[np.isfinite(act)], 99)) if np.any(np.isfinite(act)) else 1.0)
            dl = _diff_limit(sp, 98, 0.5) if diff_lim is None else diff_lim
            titles = {b: _state_titles(sp, b, byears[b]) for b in branches}
            panels = [(f"{prefix}Reference: E3SM actual {sp['target_date']} state", act, state_cmap, lo, hi, units)]
            panels += [(f"{prefix}{titles[b]['branch']}\n{_state_stats(metrics[b])}", sp["raw"][b], state_cmap, lo,
                        hi, units) for b in branches]
            panels += [(
                f"{prefix}{titles[b]['diff']}\nbias={metrics[b]['bias']:+.3g}, "
                f"box PCC={metrics[b]['r_box']:.2f}",
                sp["diff"][b], diff_cmap, -dl, dl, f"diff ({units})",
            ) for b in branches]
        elif kind == "z":
            panels = [(f"{prefix}{t}", z, diff_cmap, -z_lim, z_lim, "z (std. anom.)") for t, z in _z_panels(sp)]
        else:
            raise ValueError(f"Unknown kind={kind!r}")

        desc = descriptors[hemi]
        if len(spatials) == 1 and kind == "z":
            rows += [(desc, extent, panels[i:i + 4]) for i in range(0, len(panels), 4)]
        else:
            rows.append((desc, extent, panels))

    ncol = max(len(p) for _, _, p in rows)
    f = font_sizes(fontsize)
    pw, ph = _panel(panel_size, (3.8, 3.8), fontsize)
    fig = plt.figure(figsize=(pw * ncol, ph * len(rows)), constrained_layout=True)
    last_coll = None
    for r, (desc, extent, panels) in enumerate(rows):
        state_axes, diff_axes = [], []
        state_coll, diff_coll = None, None
        for c, (title, arr, cm, lo, hi, label) in enumerate(panels):
            ax = fig.add_subplot(len(rows), ncol, r * ncol + c + 1, projection=desc.projection)
            ax.set_facecolor("#f5f7fa")
            coll = mosaic.polypcolor(
                ax, desc, np.ma.masked_invalid(arr), **discrete_cmap(cm, lo, hi, arr), antialiaseds=False
            )
            mosaic.coastlines(ax, desc, color="0.25", linewidth=0.5)
            ax.set_extent(extent, crs=ccrs.PlateCarree())
            ax.set_boundary(circle_path, transform=ax.transAxes)
            ax.set_title(title, fontsize=f["title"])
            last_coll = coll
            if "diff" in label:
                diff_axes.append(ax)
                diff_coll = coll
            else:
                state_axes.append(ax)
                state_coll = coll

        if kind == "state_and_diff":
            first_sp = spatials[r] if r < len(spatials) else spatials[0]
            r_name = first_sp["field_name"]
            if r_name[:3] in ("NH ", "SH "):
                r_name = r_name.split(" ", 1)[1]
            r_units = first_sp["units"]
            if state_axes and state_coll is not None:
                _cbar(fig, state_coll, fontsize, ax=state_axes, orientation="horizontal",
                      fraction=0.045, pad=0.05, aspect=36, shrink=0.92, skip_ticks=False,
                      label=f"{r_name} ({r_units})")
            if diff_axes and diff_coll is not None:
                _cbar(fig, diff_coll, fontsize, ax=diff_axes, orientation="horizontal",
                      fraction=0.045, pad=0.05, aspect=18, shrink=0.88, skip_ticks=False,
                      label=f"diff ({r_units})")
        elif kind == "z" and len(rows) == 1 and state_axes:
            _cbar(fig, last_coll, fontsize, ax=state_axes, orientation="horizontal",
                  fraction=0.045, pad=0.05, aspect=36, shrink=0.92, skip_ticks=False,
                  label="standardized anomaly z")

    if kind == "z" and len(rows) > 1 and last_coll is not None:
        _cbar(fig, last_coll, fontsize, ax=fig.axes, orientation="horizontal",
              fraction=0.045, pad=0.05, aspect=42, shrink=0.85, skip_ticks=False,
              label="standardized anomaly z")

    first = spatials[0]
    name = first["field_name"]
    if len(spatials) > 1:
        name = name.split(" ", 1)[1] if name[:3] in ("NH ", "SH ") else name
    fig.suptitle(f"MPAS-Seaice {name} ({first['units']}) — {first['target_date']}",
                 fontsize=f["suptitle"], fontweight="bold")
    return fig


def plot_mpassi_four_seasons(
    field_name: str,
    test_dates: Sequence[str],
    plans: Mapping[tuple[str, str], SelectionPlan],
    selections: Mapping[tuple[str, str], Mapping[str, int]],
    descriptors: Mapping[str, object],
    *,
    branches: Sequence[str] | None = None,
    box: tuple[float, float] = (15.0, 30.0),
    cache_dir: Path | None = None,
    cmap: str | None = None,
    vmin: float = 0.0,
    vmax: float = 1.0,
    diff_lim: float = 0.4,
    fontsize: float = BASE_FONTSIZE,
    panel_size: tuple[float, float] | None = None,
):
    """Four-season comparison of MPAS-Seaice state and each branch's difference using ``mosaic``."""
    ensure_proj_data()
    import cartopy.crs as ccrs
    import matplotlib.path as mpath
    import matplotlib.pyplot as plt
    import mosaic

    theta = np.linspace(0, 2 * np.pi, 100)
    circle_path = mpath.Path(np.vstack([np.sin(theta), np.cos(theta)]).T * 0.5 + 0.5)

    if branches is None:
        branches = branches_in(selections)
    cmap = cmap or FIELD_CMAPS.get(field_name, "cmo.ice")
    hemi = "NH" if field_name.startswith("NH") else "SH"
    extent = [-180, 180, 45, 90] if hemi == "NH" else [-180, 180, -90, -45]
    desc = descriptors[hemi]
    ncol = 1 + len(branches)
    f = font_sizes(fontsize)
    pw, ph = _panel(panel_size, (3.85, 3.3), fontsize)
    fig, axes = plt.subplots(
        len(test_dates), ncol, figsize=(pw * ncol, ph * len(test_dates)),
        subplot_kw={"projection": desc.projection}, constrained_layout=True, squeeze=False,
    )
    for row, date in enumerate(test_dates):
        byears = {b: selections[(b, date)]["mpassi"] for b in branches if (b, date) in selections}
        sp = prepare_spatial_field(
            _date_plan(plans, date), "mpassi", field_name, 2015, byears, obs_plans={},
            box=box, cache_dir=cache_dir, log=lambda *_: None,
        )
        panels = [(f"{date}: Actual 2015", sp["raw"]["actual"], cmap, vmin, vmax)]
        for b in branches:
            if b not in byears:
                panels.append(None)
                continue
            m = sp["metrics"][b]
            panels.append((f"{b} ({byears[b]}) - 2015\n{_state_stats(m)}",
                           sp["diff"][b], DIFF_CMAP, -diff_lim, diff_lim))
        for col, panel in enumerate(panels):
            ax = axes[row, col]
            if panel is None:
                ax.axis("off")
                continue
            title, arr, cm, lo, hi = panel
            ax.set_facecolor("#f5f7fa")
            coll = mosaic.polypcolor(ax, desc, np.ma.masked_invalid(arr), **discrete_cmap(cm, lo, hi, arr),
                                     antialiaseds=False)
            mosaic.coastlines(ax, desc, color="0.25", linewidth=0.5)
            ax.set_extent(extent, crs=ccrs.PlateCarree())
            ax.set_boundary(circle_path, transform=ax.transAxes)
            ax.set_title(title, fontsize=f["title"])
            if row == len(test_dates) - 1:
                lbl = "concentration" if col == 0 else "difference"
                _cbar(fig, coll, fontsize, ax=ax, orientation="horizontal", fraction=0.05, pad=0.04, label=lbl)
    fig.suptitle(f"MPAS-Seaice {field_name}: 4-season comparison across 2015 dates",
                 fontsize=f["suptitle"], fontweight="bold")
    return fig


def build_elm_descriptor(
    lat: np.ndarray,
    lon: np.ndarray,
    projection: object | None = None,
    res: float = 0.5,
) -> object:
    """Create a ``mosaic.Descriptor`` for the native ELM ``grid1d`` land mesh."""
    ensure_proj_data()
    import cartopy.crs as ccrs
    import mosaic
    import xarray as xr

    if projection is None:
        projection = ccrs.Robinson()

    nlat, nlon = int(round(180.0 / res)), int(round(360.0 / res))
    lon_deg = ((lon + 180.0) % 360.0) - 180.0
    ilat = np.clip(np.rint((lat + 90.0 - 0.5 * res) / res).astype(np.int64), 0, nlat - 1)
    ilon = np.clip(np.rint((lon_deg + 180.0 - 0.5 * res) / res).astype(np.int64), 0, nlon - 1)

    lat_v1d = np.linspace(-90.0, 90.0, nlat + 1)
    lon_v1d = np.linspace(-180.0, 180.0, nlon + 1)
    lon_v1d[0], lon_v1d[-1] = -180.0 + 1e-5, 180.0 - 1e-5
    lat_v1d[0], lat_v1d[-1] = -90.0 + 1e-5, 90.0 - 1e-5

    nv_col = nlon + 1
    v_sw = ilat * nv_col + ilon
    v_se = ilat * nv_col + (ilon + 1)
    v_ne = (ilat + 1) * nv_col + (ilon + 1)
    v_nw = (ilat + 1) * nv_col + ilon
    voc_raw = np.stack([v_sw, v_se, v_ne, v_nw], axis=1)
    u_verts, voc_0 = np.unique(voc_raw, return_inverse=True)
    voc_0 = voc_0.reshape(-1, 4)
    v_ilat, v_ilon = np.divmod(u_verts, nv_col)
    lat_vertex = np.radians(lat_v1d[v_ilat])
    lon_vertex = np.radians(lon_v1d[v_ilon])

    ncells = lat.size
    nverts = u_verts.size
    nedges = 2 * ncells
    cell_ids = np.arange(1, ncells + 1, dtype=np.int64)
    edge_ids = np.arange(1, nedges + 1, dtype=np.int64)

    cells_on_vert = np.zeros((nverts, 4), dtype=np.int64)
    edges_on_vert = np.zeros((nverts, 4), dtype=np.int64)
    for k in range(4):
        cells_on_vert[voc_0[:, k], k] = cell_ids
        edges_on_vert[voc_0[:, k], k] = edge_ids[k // 2 :: 2]

    cells_on_edge = np.column_stack([
        np.repeat(cell_ids, 2),
        np.zeros(nedges, dtype=np.int64),
    ])
    vertices_on_edge = np.vstack([
        np.column_stack([voc_0[:, 0] + 1, voc_0[:, 1] + 1]),
        np.column_stack([voc_0[:, 2] + 1, voc_0[:, 3] + 1]),
    ])
    lat_edge = np.repeat(np.radians(lat), 2)
    lon_edge = np.repeat(np.radians(lon_deg), 2)

    ds_elm = xr.Dataset(
        data_vars={
            "latCell": ("nCells", np.radians(lat)),
            "lonCell": ("nCells", np.radians(lon_deg)),
            "latVertex": ("nVertices", lat_vertex),
            "lonVertex": ("nVertices", lon_vertex),
            "latEdge": ("nEdges", lat_edge),
            "lonEdge": ("nEdges", lon_edge),
            "verticesOnCell": (("nCells", "maxEdges"), voc_0 + 1),
            "cellsOnEdge": (("nEdges", "TWO"), cells_on_edge),
            "cellsOnVertex": (("nVertices", "vertexDegree"), cells_on_vert),
            "verticesOnEdge": (("nEdges", "TWO"), vertices_on_edge),
            "edgesOnVertex": (("nVertices", "vertexDegree"), edges_on_vert),
        },
        attrs={"on_a_sphere": "YES"},
    )
    desc = mosaic.Descriptor(ds_elm, projection, ccrs.PlateCarree(), use_latlon=True)
    _ = desc.cell_patches
    return desc


def plot_elm_spatial(
    spatial: dict | Sequence[dict],
    elm_descriptor: object | None = None,
    *,
    kind: str = "state_and_diff",
    cmap: str | None = None,
    diff_cmap: str = DIFF_CMAP,
    vmin: float | None = None,
    vmax: float | None = None,
    diff_lim: float | None = None,
    diff_percent: bool = True,
    z_lim: float = 3.0,
    z_levels: Sequence[float] | None = None,
    fontsize: float = BASE_FONTSIZE,
    panel_size: tuple[float, float] | None = None,
):
    """Plot ELM fields on the native ``grid1d`` land mesh using ``mosaic.polypcolor``.

    ``state_and_diff``: one row with actual, each branch, the candidate-year
    climatology, and each branch minus actual. ``z``: actual, each observation
    reference, and each branch. ``spatial`` may be a list of fields, drawn one
    row each in a single figure.
    """
    ensure_proj_data()
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import matplotlib.pyplot as plt
    import mosaic

    spatials = [spatial] if isinstance(spatial, Mapping) else list(spatial)
    spatial = spatials[0]
    fname = spatial["field_name"]
    tyear = spatial["target_year"]
    byears = spatial["branch_years"]
    metrics = spatial["metrics"]
    lat, lon = spatial["lat"], spatial["lon"]

    if elm_descriptor is None:
        elm_descriptor = build_elm_descriptor(lat, lon)
    proj = elm_descriptor.projection

    f = font_sizes(fontsize)

    def draw(ax, arr, cm, lo, hi, title, levels=None):
        ax.set_facecolor("#f5f7fa")
        ax.add_feature(cfeature.LAND, facecolor="#eceff1", zorder=0)
        coll = mosaic.polypcolor(
            ax, elm_descriptor, np.ma.masked_invalid(arr), **discrete_cmap(cm, lo, hi, arr, levels=levels),
            antialiaseds=False, zorder=2,
        )
        ax.coastlines(linewidth=0.45, color="0.3", zorder=3)
        ax.set_extent([-180, 180, -60, 85], crs=ccrs.PlateCarree())
        ax.set_title(title, fontsize=f["title"])
        return coll

    if kind == "state_and_diff":
        # One row per field: actual, each branch, the candidate-year climatology,
        # and each branch minus actual.
        rows = []
        for sp in spatials:
            fn, units, by, met = sp["field_name"], sp["units"], sp["branch_years"], sp["metrics"]
            branches = [b for b in by if b in sp["raw"]]
            cm = cmap or FIELD_CMAPS.get(fn, "viridis")
            act = sp["raw"]["actual"]
            finite_act = act[np.isfinite(act)]
            levels = FIELD_LEVELS.get(fn)
            lo = vmin if vmin is not None else (
                float(np.nanpercentile(finite_act, 2)) if "temperature" in fn else 0.0)
            hi = vmax if vmax is not None else (
                float(np.nanpercentile(finite_act, 98)) if finite_act.size else 1.0)
            dl = diff_lim if diff_lim is not None else _diff_limit(sp, 97, 1.0)
            label = f"{fn} ({units})"
            titles = {b: _state_titles(sp, b, by[b]) for b in branches}
            clim_title = _state_titles(sp, "", 0)["clim"]
            clim_full = f"{clim_title}\n{_state_stats(met['climatology'])}"
            ref_title = f"Reference: E3SM actual {sp['target_date']} state"

            # Baseline climatology in column 1
            panels = [(f"[{fn}]\n{clim_full}" if len(spatials) > 1 else clim_full,
                       sp["raw"]["climatology"], cm, lo, hi, label, levels)]
            # Reference in column 2
            panels.append((ref_title, act, cm, lo, hi, label, levels))
            # Selected initial condition(s) in column 3
            panels += [(f"{titles[b]['branch']}\n{_state_stats(met[b])}", sp["raw"][b], cm, lo, hi, label, levels)
                       for b in branches]
            # Error (difference) in column 4
            for b in branches:
                m = met[b]
                raw_b = sp["raw"][b]
                if diff_percent:
                    diff_arr = _compute_diff_percent(fn, raw_b, act)
                    pct_levels = DIFF_PCT_LEVELS.get(fn, [-100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100])
                    dl_pct = float(pct_levels[-1])
                    finite_act = act[np.isfinite(act)]
                    finite_b = raw_b[np.isfinite(raw_b)]
                    mean_act = float(np.mean(finite_act)) if finite_act.size else 0.0
                    mean_b = float(np.mean(finite_b)) if finite_b.size else 0.0
                    if mean_act > 1e-4:
                        pct_bias = 100.0 * (mean_b - mean_act) / mean_act
                        bias_str = f"bias={pct_bias:+.1f}%"
                    else:
                        bias_str = f"bias={m['bias']:+.2g}"
                    panels.append((
                        f"{titles[b]['diff']}\n{bias_str}, box PCC={m['r_box']:.2f}",
                        diff_arr, diff_cmap, -dl_pct, dl_pct, "diff (%)", pct_levels,
                    ))
                else:
                    panels.append((
                        f"{titles[b]['diff']}\nbias={m['bias']:+.3g}, box PCC={m['r_box']:.2f}",
                        sp["diff"][b], diff_cmap, -dl, dl, f"diff ({units})", None,
                    ))
            rows.append(panels)

        ncol = max(len(r) for r in rows)
        pw, ph = _panel(panel_size, (4.4, 3.3), fontsize)
        fig, axes = plt.subplots(
            len(rows), ncol, figsize=(pw * ncol, ph * len(rows)),
            subplot_kw={"projection": proj}, constrained_layout=True, squeeze=False,
        )
        for r, panels in enumerate(rows):
            state_axes, diff_axes = [], []
            state_coll, diff_coll = None, None
            for c in range(ncol):
                if c >= len(panels):
                    axes[r, c].axis("off")
                    continue
                title, arr, cm, lo, hi, label, lev = panels[c]
                coll = draw(axes[r, c], arr, cm, lo, hi, title, levels=lev)
                if label.startswith("diff"):
                    diff_axes.append(axes[r, c])
                    diff_coll = coll
                else:
                    state_axes.append(axes[r, c])
                    state_coll = coll
            fn = spatials[r]["field_name"] if r < len(spatials) else spatials[0]["field_name"]
            units = spatials[r]["units"] if r < len(spatials) else spatials[0]["units"]
            if state_axes and state_coll is not None:
                _cbar(
                    fig, state_coll, fontsize, ax=state_axes, orientation="horizontal",
                    fraction=0.045, pad=0.05, aspect=36, shrink=0.92, skip_ticks=False,
                    label=f"{fn} ({units})",
                )
            if diff_axes and diff_coll is not None:
                diff_label = "diff (%)" if diff_percent else f"diff ({units})"
                _cbar(
                    fig, diff_coll, fontsize, ax=diff_axes, orientation="horizontal",
                    fraction=0.045, pad=0.05, aspect=18, shrink=0.88, skip_ticks=False,
                    label=diff_label,
                )
        what = "key fields" if len(spatials) > 1 else f"{fname} ({spatial['units']})"
        fig.suptitle(f"ELM {what} — {spatial['target_date']}", fontsize=f["suptitle"], fontweight="bold")
        return fig

    if kind == "z":
        if len(spatials) > 1:
            rows = []
            for sp in spatials:
                panels = _z_panels(sp, with_rgrid_ref=True)
                title, arr = panels[0]
                rows.append([(f"[{sp['field_name']}]\n{title}", arr)] + panels[1:])
        else:
            panels = _z_panels(spatial, with_rgrid_ref=True)
            ncol = min(len(panels), 3)
            rows = [panels[i:i + ncol] for i in range(0, len(panels), ncol)]
        ncol = max(len(r) for r in rows)
        nrow = len(rows)
        pw, ph = _panel(panel_size, (4.6, 3.4), fontsize)
        fig, axes = plt.subplots(
            nrow, ncol, figsize=(pw * ncol, ph * nrow),
            subplot_kw={"projection": proj}, constrained_layout=True, squeeze=False,
        )
        coll = None
        levs = z_levels if z_levels is not None else Z_LEVELS
        lo_z, hi_z = float(levs[0]), float(levs[-1])
        for r, panels in enumerate(rows):
            for c in range(ncol):
                if c >= len(panels):
                    axes[r, c].axis("off")
                    continue
                title, arr = panels[c]
                coll = draw(axes[r, c], arr, diff_cmap, lo_z, hi_z, title, levels=levs)
        _cbar(fig, coll, fontsize, ax=axes.ravel().tolist(), orientation="horizontal",
              fraction=0.045, pad=0.05, aspect=42, shrink=0.85, skip_ticks=False,
              label="standardized anomaly z")
        what = "key fields" if len(spatials) > 1 else fname
        fig.suptitle(f"ELM {what} standardized anomalies — {spatial['target_date']}",
                     fontsize=f["suptitle"], fontweight="bold")
        return fig

    raise ValueError(f"Unknown kind={kind!r}")


def plot_elm_four_seasons(
    field_name: str,
    test_dates: Sequence[str],
    plans: Mapping[tuple[str, str], SelectionPlan],
    selections: Mapping[tuple[str, str], Mapping[str, int]],
    elm_descriptor: object,
    *,
    branches: Sequence[str] | None = None,
    obs_refs: Mapping[str, str] | None = None,
    box: tuple[float, float] = (15.0, 30.0),
    cache_dir: Path | None = None,
    z_lim: float = 3.0,
    z_levels: Sequence[float] | None = None,
    fontsize: float = BASE_FONTSIZE,
    panel_size: tuple[float, float] | None = None,
):
    """Four-season comparison of ELM standardized anomalies on the native mesh via ``mosaic``.

    Columns: actual, each observation reference (``obs_refs`` maps an obs
    branch to its label; default ERA5), and each branch.
    """
    ensure_proj_data()
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import matplotlib.pyplot as plt
    import mosaic

    if branches is None:
        branches = branches_in(selections)
    rows = []
    for date in test_dates:
        byears = {b: selections[(b, date)]["elm"] for b in branches if (b, date) in selections}
        sp = prepare_spatial_field(
            _date_plan(plans, date), "elm", field_name, 2015, byears,
            obs_plans=_obs_plans_for(plans, date, obs_refs),
            box=box, cache_dir=cache_dir, log=lambda *_: None,
        )
        m = sp["metrics"]
        panels = _z_panels(sp)
        rows.append(panels)

    ncol = max(len(r) for r in rows)
    f = font_sizes(fontsize)
    pw, ph = _panel(panel_size, (4.1, 2.9), fontsize)
    fig, axes = plt.subplots(
        len(test_dates), ncol, figsize=(pw * ncol, ph * len(test_dates)),
        subplot_kw={"projection": elm_descriptor.projection}, constrained_layout=True, squeeze=False,
    )
    coll = None
    levs = z_levels if z_levels is not None else Z_LEVELS
    lo_z, hi_z = float(levs[0]), float(levs[-1])
    for row, panels in enumerate(rows):
        for col in range(ncol):
            ax = axes[row, col]
            if col >= len(panels):
                ax.axis("off")
                continue
            title, arr = panels[col]
            ax.set_facecolor("#f5f7fa")
            ax.add_feature(cfeature.LAND, facecolor="#eceff1", zorder=0)
            coll = mosaic.polypcolor(
                ax, elm_descriptor, np.ma.masked_invalid(arr),
                **discrete_cmap(DIFF_CMAP, lo_z, hi_z, levels=levs), antialiaseds=False, zorder=2,
            )
            ax.coastlines(linewidth=0.45, color="0.3", zorder=3)
            ax.set_extent([-180, 180, -60, 85], crs=ccrs.PlateCarree())
            ax.set_title(title, fontsize=f["title"])
    _cbar(fig, coll, fontsize, ax=axes.ravel().tolist(), orientation="horizontal",
          fraction=0.045, pad=0.05, aspect=42, shrink=0.85, skip_ticks=False,
          label="standardized anomaly z")
    fig.suptitle(f"ELM {field_name}: 4-season standardized anomalies",
                 fontsize=f["suptitle"], fontweight="bold")
    return fig
