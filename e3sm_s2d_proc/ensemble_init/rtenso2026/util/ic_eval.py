"""Similarity evaluation and mosaic/cartopy spatial plotting for RTENSO2026 ICs.

Compares climatology-matched (``clim``) and observationally matched (``obs``)
restarts against the actual 2015 BruteForce restarts across key MPAS-Seaice,
ELM, and MOSART fields:

* Distance and rank checks against actual 2015 (``evaluate_distance_ranks``,
  ``plot_distance_vs_truth``).
* Field-level similarity metrics on 30x60 boxes and native model grids
  (``field_similarity_table``, ``plot_field_similarity_bars``,
  ``global_summary_table``).
* Spatial comparisons on 30x60 boxes (``plot_box_spatial``), on the native
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
}


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

def load_test_selections(
    test_dates: Sequence[str],
    src_root: Path,
    out_root: Path,
    *,
    era5_root: str | Path | Mapping[str, str | Path | None] | None = None,
    obs_root: str | Path | Mapping[str, str | Path | None] | None = None,
    obs_fields: Mapping[str, Mapping[str, str] | None] | None = None,
    branches: Sequence[str] = BRANCHES,
    components: Sequence[str] = ALL_COMPONENTS,
    src_prefix: str = "20251024_s2d_spinup",
    out_prefix: str = "20251024_s2d_spinup",
    box: tuple[float, float] = (30.0, 60.0),
    mode: str = "joint",
    component_weights: Mapping[str, float] | None = None,
    cache_dir: Path | None = None,
    log=lambda *_: None,
) -> tuple[dict, dict, dict]:
    """Load plans, distance tables, and selected years for ``test_dates``."""
    plans, distances, selections = {}, {}, {}
    for date in test_dates:
        for branch in branches:
            if branch == "obs" and not era5_available(
                int(date[:4]), date[5:], components,
                era5_root=era5_root, obs_root=obs_root, obs_fields=obs_fields,
            ):
                continue
            plan = make_plan(
                date, src_root, out_root,
                era5_root=era5_root, obs_root=obs_root, obs_fields=obs_fields,
                branch=branch, components=components,
                src_prefix=src_prefix, out_prefix=out_prefix, exclude_target_year=True,
            )
            table = compute_distances(
                plan, box=box, component_weights=component_weights, cache_dir=cache_dir, log=log
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
    box: tuple[float, float] = (30.0, 60.0),
    component_weights: Mapping[str, float] | None = None,
    cache_dir: Path | None = None,
    log=lambda *_: None,
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame], dict[str, pd.DataFrame]]:
    """Compare ``clim`` and ``obs`` distances against true distance to ``target_year``."""
    truth_obs, truth_all = {}, {}
    rows = []
    for date in test_dates:
        ref_plan = plans.get(("clim", date)) or next(p for (b, d), p in plans.items() if d == date)
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
        for branch in BRANCHES:
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
    branches: Sequence[str] = BRANCHES,
):
    """Scatter plot of candidate-year branch distance vs true distance to 2015."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(
        len(branches), len(test_dates),
        figsize=(4.2 * len(test_dates), 3.6 * len(branches)),
        sharey="row", squeeze=False,
    )
    for col, date in enumerate(test_dates):
        t_dist = truth_obs[date]["joint"]
        best_true_yr = int(t_dist.idxmin())
        for row, branch in enumerate(branches):
            ax = axes[row, col]
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
                            textcoords="offset points", xytext=(5, 4), fontsize=8.5)

            r_spearman = b_dist.rank().corr(t_dist.rank())
            true_rk = int(t_dist.rank()[sel_yr])
            ax.set_title(f"{branch} {date}\nsel={sel_yr} (rank {true_rk}/{len(b_dist)}), rank corr={r_spearman:.2f}",
                         fontsize=10)
            ax.set_xlabel(f"{branch} distance $D$")
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=7.5, loc="lower right")
        for row in range(len(branches)):
            axes[row, 0].set_ylabel(r"True distance to 2015 ($D_{\mathrm{true}}$, obs fields)")

    fig.suptitle("Candidate-year selection distance vs. true distance to actual 2015 IC",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
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
    box: tuple[float, float] = (30.0, 60.0),
    cache_dir: Path | None = None,
    log=lambda *_: None,
) -> pd.DataFrame:
    """Detailed similarity metrics comparing selected years (and ERA5) against ``target_year``."""
    cand_years = list(plan.candidate_years)
    year_to_idx = {y: i for i, y in enumerate(cand_years)}
    rows = []

    for comp in plan.components:
        if comp not in KEY_FIELDS:
            continue
        fields = load_key_fields(plan, comp, cache_dir, years=cand_years + [target_year], log=log)
        pts = comparison_points(plan, comp, cache_dir)
        has_obs = comp in OBS_FIELDS and era5_available(target_year, plan.mmdd, [comp])
        obs = (load_obs_fields(plan, comp, cand_years + [target_year], cache_dir, log=log)
               if has_obs else {})

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
                    "bias": float((wg_valid * diff).sum() / wsum_g) if wsum_g > 0 else np.nan,
                })

            if fname in obs:
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
                    "source": "ERA5",
                    "year": int(target_year),
                    "d_box": _weighted_rms_diff(zt_b_obs, zt_b_mod, wb_ov),
                    "r_box": _weighted_corr(zt_b_obs, zt_b_mod, wb_ov),
                    "d_grid": _weighted_rms_diff(zt_g_obs, zt_g_mod, wg_o),
                    "r_grid": _weighted_corr(zt_g_obs, zt_g_mod, wg_o),
                    "RMSE": np.nan,
                    "bias": np.nan,
                })

    return pd.DataFrame(rows)


def plot_field_similarity_bars(sim_df: pd.DataFrame):
    """Horizontal bar chart of four-season mean d_box and r_box by key field."""
    import matplotlib.pyplot as plt

    mean_sim = (
        sim_df.groupby(["component", "field", "source"], as_index=False)[["d_box", "r_box", "d_grid", "r_grid", "RMSE"]]
        .mean()
    )
    pivot_d = mean_sim.pivot(index=["component", "field"], columns="source", values="d_box")[["clim", "obs", "ERA5"]]
    pivot_r = mean_sim.pivot(index=["component", "field"], columns="source", values="r_box")[["clim", "obs", "ERA5"]]

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True, constrained_layout=True)
    labels = [f"{c}: {f}" for c, f in pivot_d.index]
    y_pos = np.arange(len(labels))
    width = 0.26

    for ax, piv, xlabel, vline in (
        (axes[0], pivot_d, r"Mean box distance $d_{\mathrm{box}}$ to actual 2015 (smaller is closer)", 1.0),
        (axes[1], pivot_r, r"Mean box anomaly correlation $r_{\mathrm{box}}$ with actual 2015 (larger is closer)", 0.0),
    ):
        ax.barh(y_pos - width, piv["clim"], height=width, label="clim IC", color="tab:blue", alpha=0.85)
        ax.barh(y_pos, piv["obs"], height=width, label="obs IC", color="tab:orange", alpha=0.85)
        ax.barh(y_pos + width, piv["ERA5"], height=width, label="ERA5 target anomaly", color="tab:green", alpha=0.75)
        ax.axvline(vline, color="k", lw=0.8, ls="--")
        ax.set_yticks(y_pos, labels)
        ax.invert_yaxis()
        ax.set_xlabel(xlabel)
        ax.grid(True, axis="x", alpha=0.3)
        ax.legend(loc="best")

    fig.suptitle("Four-season mean similarity to actual 2015 IC by field (2015 sensitivity dates)",
                 fontsize=12.5, fontweight="bold")
    return fig


def global_summary_table(
    test_dates: Sequence[str],
    selections: Mapping[tuple[str, str], Mapping[str, int]],
    src_root: Path,
    src_prefix: str = "20251024_s2d_spinup",
    components: Sequence[str] = ALL_COMPONENTS,
    branches: Sequence[str] = BRANCHES,
    target_year: int = 2015,
) -> pd.DataFrame:
    """Global integrated quantities comparing actual ``target_year`` against ``clim`` and ``obs``."""
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
# Spatial field preparation and plotting (30x60 boxes, MPAS-Seaice, ELM)
# ---------------------------------------------------------------------------

def prepare_spatial_field(
    plan: SelectionPlan,
    component: str,
    field_name: str,
    target_year: int,
    branch_years: Mapping[str, int],
    *,
    box: tuple[float, float] = (30.0, 60.0),
    cache_dir: Path | None = None,
    log=lambda *_: None,
) -> dict:
    """Prepare native-grid and box-averaged fields, anomalies, and differences for plotting."""
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
            "bias": float((wg_valid * df_diff).sum() / wsum_g) if wsum_g > 0 else np.nan,
        }

    if component in OBS_FIELDS and field_name in OBS_FIELDS[component] and era5_available(target_year, plan.mmdd, [component]):
        obs_all = load_obs_fields(plan, component, cand_years + [target_year], cache_dir, log=log)[field_name]
        valid_both = np.all(np.isfinite(vals_grid), axis=0) & np.all(np.isfinite(obs_all), axis=0)
        _, zt_g_obs, _, _, wg_o = _standardize_series(obs_all[:-1], w_grid, obs_all[-1], valid_mask=valid_both)
        _, zt_g_mod, _, _, _ = _standardize_series(cand_grid, w_grid, targ_grid, valid_mask=valid_both)
        z_map["era5"] = zt_g_obs

        m_box, wb_o, o_box = _apply_box(vals_grid, w_grid, lat, lon, box, obs=obs_all)
        _, zt_b_mod, _, _, wb_ov = _standardize_series(m_box[:-1], wb_o, m_box[-1])
        _, zt_b_obs, _, _, _ = _standardize_series(o_box[:-1], wb_o, o_box[-1])
        box_z["era5"] = np.where(wb_ov > 0, zt_b_obs, np.nan).reshape(nlat, nlon)
        box_diff_z["era5"] = np.where(wb_ov > 0, zt_b_obs - zt_b_mod, np.nan).reshape(nlat, nlon)

        metrics["era5"] = {
            "year": int(target_year),
            "d_box": _weighted_rms_diff(zt_b_obs, zt_b_mod, wb_ov),
            "r_box": _weighted_corr(zt_b_obs, zt_b_mod, wb_ov),
            "d_grid": _weighted_rms_diff(zt_g_obs, zt_g_mod, wg_o),
            "r_grid": _weighted_corr(zt_g_obs, zt_g_mod, wg_o),
            "rmse": np.nan,
            "bias": np.nan,
        }

    return {
        "component": component,
        "field_name": field_name,
        "units": FIELD_UNITS.get(field_name, ""),
        "target_date": plan.target_date,
        "target_year": int(target_year),
        "branch_years": {k: int(v) for k, v in branch_years.items()},
        "lat": lat,
        "lon": lon,
        "weight": wg_valid,
        "clim_mean": clim_mean,
        "spread": spread,
        "raw": raw,
        "diff": diff,
        "z": z_map,
        "metrics": metrics,
        "box": {
            "dlat": dlat,
            "dlon": dlon,
            "lat_edges": lat_edges,
            "lon_edges": lon_edges,
            "z": box_z,
            "diff_z": box_diff_z,
        },
    }


def plot_box_spatial(spatial: dict, *, z_lim: float = 2.0):
    """Plot 30x60-degree box-averaged standardized anomalies used in the distance metric."""
    ensure_proj_data()
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt

    fname = spatial["field_name"]
    tyear = spatial["target_year"]
    byears = spatial["branch_years"]
    metrics = spatial["metrics"]
    binfo = spatial["box"]
    lon_e, lat_e = binfo["lon_edges"], binfo["lat_edges"]

    panels = [(f"Actual {tyear} box z", binfo["z"]["actual"])]
    if "era5" in binfo["z"]:
        m = metrics["era5"]
        panels.append((f"ERA5 {tyear} box z\nd_box={m['d_box']:.2f}, r_box={m['r_box']:.2f}",
                       binfo["z"]["era5"]))
    for b in ("obs", "clim"):
        if b in binfo["z"]:
            m = metrics[b]
            panels.append((f"{b} ({byears[b]}) box z\nd_box={m['d_box']:.2f}, r_box={m['r_box']:.2f}",
                           binfo["z"][b]))

    fig, axes = plt.subplots(
        1, len(panels), figsize=(4.2 * len(panels), 3.1),
        subplot_kw={"projection": ccrs.PlateCarree(central_longitude=180)},
        constrained_layout=True, squeeze=False,
    )
    pm = None
    for ax, (title, grid_z) in zip(axes[0], panels):
        ax.set_facecolor("0.92")
        pm = ax.pcolormesh(lon_e, lat_e, grid_z, transform=ccrs.PlateCarree(),
                           cmap="RdBu_r", vmin=-z_lim, vmax=z_lim, edgecolors="0.6", linewidth=0.3)
        ax.coastlines(linewidth=0.5, color="0.2")
        ax.set_global()
        ax.set_title(title, fontsize=9.5)
    fig.colorbar(pm, ax=axes[0, :].tolist(), orientation="horizontal",
                 fraction=0.06, pad=0.05, label=f"{int(binfo['dlat'])}°x{int(binfo['dlon'])}° box standardized anomaly z")
    fig.suptitle(f"Box-scale ({int(binfo['dlat'])}°x{int(binfo['dlon'])}°) anomalies: {spatial['component']} {fname} — {spatial['target_date']}",
                 fontsize=11.5, fontweight="bold")
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


def plot_mpassi_spatial(
    spatial: dict,
    descriptors: Mapping[str, object],
    *,
    kind: str = "state_and_diff",
    cmap: str = "YlGnBu_r",
    diff_cmap: str = "RdBu_r",
    vmin: float | None = None,
    vmax: float | None = None,
    diff_lim: float | None = None,
    z_lim: float = 2.5,
):
    """Plot MPAS-Seaice fields on the native unstructured MPAS mesh using ``mosaic``."""
    ensure_proj_data()
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt
    import mosaic

    fname = spatial["field_name"]
    hemi = "NH" if fname.startswith("NH") else "SH"
    desc = descriptors[hemi]
    proj = desc.projection
    extent = [-180, 180, 45, 90] if hemi == "NH" else [-180, 180, -90, -45]
    tyear = spatial["target_year"]
    units = spatial["units"]
    byears = spatial["branch_years"]
    metrics = spatial["metrics"]

    if kind == "state_and_diff":
        act = spatial["raw"]["actual"]
        if vmin is None:
            vmin = 0.0
        if vmax is None:
            vmax = float(np.nanpercentile(act[np.isfinite(act)], 99)) if np.any(np.isfinite(act)) else 1.0
            if "concentration" in fname:
                vmax = 1.0
        if diff_lim is None:
            diffs = [spatial["diff"][b] for b in ("clim", "obs") if b in spatial["diff"]]
            all_d = np.concatenate([d[np.isfinite(d)] for d in diffs]) if diffs else np.array([0.5])
            diff_lim = float(np.nanpercentile(np.abs(all_d), 98)) if all_d.size else 0.5
            diff_lim = max(diff_lim, 1e-3)

        panels = [(f"Actual {tyear}", act, cmap, vmin, vmax)]
        for b in ("clim", "obs"):
            if b in spatial["raw"]:
                panels.append((f"{b} ({byears[b]})", spatial["raw"][b], cmap, vmin, vmax))
        for b in ("clim", "obs"):
            if b in spatial["diff"]:
                m = metrics[b]
                panels.append((
                    f"{b} ({byears[b]}) - {tyear}\nRMSE={m['rmse']:.3g}, d_box={m['d_box']:.2f}, r_box={m['r_box']:.2f}",
                    spatial["diff"][b], diff_cmap, -diff_lim, diff_lim,
                ))
    elif kind == "z":
        panels = [(f"Actual {tyear} (z)", spatial["z"]["actual"], diff_cmap, -z_lim, z_lim)]
        if "era5" in spatial["z"]:
            m = metrics["era5"]
            panels.append((
                f"ERA5 {tyear} (z)\nd_box={m['d_box']:.2f}, r_box={m['r_box']:.2f}",
                spatial["z"]["era5"], diff_cmap, -z_lim, z_lim,
            ))
        for b in ("obs", "clim"):
            if b in spatial["z"]:
                m = metrics[b]
                panels.append((
                    f"{b} ({byears[b]}) (z)\nd_box={m['d_box']:.2f}, r_box={m['r_box']:.2f}, r_grid={m['r_grid']:.2f}",
                    spatial["z"][b], diff_cmap, -z_lim, z_lim,
                ))
    else:
        raise ValueError(f"Unknown kind={kind!r}")

    fig, axes = plt.subplots(
        1, len(panels), figsize=(3.8 * len(panels), 4.2),
        subplot_kw={"projection": proj}, constrained_layout=True, squeeze=False,
    )
    for ax, (title, arr, cm, lo, hi) in zip(axes[0], panels):
        ax.set_facecolor("0.85")
        coll = mosaic.polypcolor(
            ax, desc, np.ma.masked_invalid(arr), cmap=cm, vmin=lo, vmax=hi, antialiaseds=False
        )
        mosaic.coastlines(ax, desc, color="0.2", linewidth=0.5)
        ax.set_extent(extent, crs=ccrs.PlateCarree())
        ax.set_title(title, fontsize=9.5)
        fig.colorbar(coll, ax=ax, orientation="horizontal", fraction=0.055, pad=0.04,
                     label=units if cm == cmap else ("z (std. anom.)" if kind == "z" else f"diff ({units})"))
    fig.suptitle(f"MPAS-Seaice {fname} ({units}) — {spatial['target_date']}", fontsize=12, fontweight="bold")
    return fig


def plot_mpassi_four_seasons(
    field_name: str,
    test_dates: Sequence[str],
    plans: Mapping[tuple[str, str], SelectionPlan],
    selections: Mapping[tuple[str, str], Mapping[str, int]],
    descriptors: Mapping[str, object],
    *,
    box: tuple[float, float] = (30.0, 60.0),
    cache_dir: Path | None = None,
    cmap: str = "YlGnBu_r",
    vmin: float = 0.0,
    vmax: float = 1.0,
    diff_lim: float = 0.4,
):
    """Four-season comparison of MPAS-Seaice state and differences using ``mosaic``."""
    ensure_proj_data()
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt
    import mosaic

    hemi = "NH" if field_name.startswith("NH") else "SH"
    extent = [-180, 180, 45, 90] if hemi == "NH" else [-180, 180, -90, -45]
    desc = descriptors[hemi]
    fig, axes = plt.subplots(
        len(test_dates), 3, figsize=(11.5, 3.3 * len(test_dates)),
        subplot_kw={"projection": desc.projection}, constrained_layout=True,
    )
    for row, date in enumerate(test_dates):
        byears = {b: selections[(b, date)]["mpassi"] for b in BRANCHES if (b, date) in selections}
        sp = prepare_spatial_field(
            plans[("clim", date)], "mpassi", field_name, 2015, byears,
            box=box, cache_dir=cache_dir, log=lambda *_: None,
        )
        panels = [
            (f"{date}: Actual 2015", sp["raw"]["actual"], cmap, vmin, vmax),
            (f"clim ({byears['clim']}) - 2015 (d_box={sp['metrics']['clim']['d_box']:.2f}, r={sp['metrics']['clim']['r_box']:.2f})",
             sp["diff"]["clim"], "RdBu_r", -diff_lim, diff_lim),
            (f"obs ({byears['obs']}) - 2015 (d_box={sp['metrics']['obs']['d_box']:.2f}, r={sp['metrics']['obs']['r_box']:.2f})",
             sp["diff"]["obs"], "RdBu_r", -diff_lim, diff_lim),
        ]
        for col, (title, arr, cm, lo, hi) in enumerate(panels):
            ax = axes[row, col]
            ax.set_facecolor("0.85")
            coll = mosaic.polypcolor(ax, desc, np.ma.masked_invalid(arr), cmap=cm, vmin=lo, vmax=hi, antialiaseds=False)
            mosaic.coastlines(ax, desc, color="0.2", linewidth=0.5)
            ax.set_extent(extent, crs=ccrs.PlateCarree())
            ax.set_title(title, fontsize=9.5)
            if row == len(test_dates) - 1:
                fig.colorbar(coll, ax=ax, orientation="horizontal", fraction=0.05, pad=0.04)
    fig.suptitle(f"MPAS-Seaice {field_name}: 4-season comparison across 2015 dates", fontsize=13, fontweight="bold")
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
    spatial: dict,
    elm_descriptor: object | None = None,
    *,
    kind: str = "state_and_diff",
    cmap: str = "viridis",
    diff_cmap: str = "RdBu_r",
    vmin: float | None = None,
    vmax: float | None = None,
    diff_lim: float | None = None,
    z_lim: float = 2.5,
):
    """Plot ELM fields on the native ``grid1d`` land mesh using ``mosaic.polypcolor``."""
    ensure_proj_data()
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt
    import mosaic

    fname = spatial["field_name"]
    tyear = spatial["target_year"]
    units = spatial["units"]
    byears = spatial["branch_years"]
    metrics = spatial["metrics"]
    lat, lon = spatial["lat"], spatial["lon"]

    if elm_descriptor is None:
        elm_descriptor = build_elm_descriptor(lat, lon)
    proj = elm_descriptor.projection

    if kind == "state_and_diff":
        act = spatial["raw"]["actual"]
        finite_act = act[np.isfinite(act)]
        if vmin is None:
            vmin = float(np.nanpercentile(finite_act, 2)) if "temperature" in fname else 0.0
        if vmax is None:
            vmax = float(np.nanpercentile(finite_act, 98)) if finite_act.size else 1.0
        if diff_lim is None:
            diffs = [spatial["diff"][b] for b in ("clim", "obs") if b in spatial["diff"]]
            all_d = np.concatenate([d[np.isfinite(d)] for d in diffs]) if diffs else np.array([1.0])
            diff_lim = float(np.nanpercentile(np.abs(all_d), 97)) if all_d.size else 1.0
            diff_lim = max(diff_lim, 1e-3)

        fig, axes = plt.subplots(
            2, 3, figsize=(15, 7.2),
            subplot_kw={"projection": proj}, constrained_layout=True,
        )
        state_panels = [
            (axes[0, 0], f"Actual {tyear}", act),
            (axes[0, 1], f"clim ({byears.get('clim', '--')})", spatial["raw"].get("clim")),
            (axes[0, 2], f"obs ({byears.get('obs', '--')})", spatial["raw"].get("obs")),
        ]
        coll_state = None
        for ax, title, arr in state_panels:
            if arr is None:
                ax.axis("off")
                continue
            ax.set_facecolor("0.92")
            coll_state = mosaic.polypcolor(
                ax, elm_descriptor, np.ma.masked_invalid(arr),
                cmap=cmap, vmin=vmin, vmax=vmax, antialiaseds=False,
            )
            ax.coastlines(linewidth=0.5, color="0.25")
            ax.set_extent([-180, 180, -60, 85], crs=ccrs.PlateCarree())
            ax.set_title(title, fontsize=10)
        fig.colorbar(coll_state, ax=axes[0, :].tolist(), orientation="vertical",
                     fraction=0.018, pad=0.02, label=units)

        ax0 = axes[1, 0]
        ax0.set_facecolor("0.92")
        mosaic.polypcolor(
            ax0, elm_descriptor, np.ma.masked_invalid(spatial["raw"]["climatology"]),
            cmap=cmap, vmin=vmin, vmax=vmax, antialiaseds=False,
        )
        ax0.coastlines(linewidth=0.5, color="0.25")
        ax0.set_extent([-180, 180, -60, 85], crs=ccrs.PlateCarree())
        ax0.set_title("Candidate-year climatology", fontsize=10)

        coll_diff = None
        for col, b in ((1, "clim"), (2, "obs")):
            ax = axes[1, col]
            if b not in spatial["diff"]:
                ax.axis("off")
                continue
            m = metrics[b]
            ax.set_facecolor("0.92")
            coll_diff = mosaic.polypcolor(
                ax, elm_descriptor, np.ma.masked_invalid(spatial["diff"][b]),
                cmap=diff_cmap, vmin=-diff_lim, vmax=diff_lim, antialiaseds=False,
            )
            ax.coastlines(linewidth=0.5, color="0.25")
            ax.set_extent([-180, 180, -60, 85], crs=ccrs.PlateCarree())
            ax.set_title(
                f"{b} ({byears[b]}) - {tyear}\nRMSE={m['rmse']:.2f}, bias={m['bias']:+.2f}, d_box={m['d_box']:.2f}, r_box={m['r_box']:.2f}",
                fontsize=9.5,
            )
        if coll_diff is not None:
            fig.colorbar(coll_diff, ax=axes[1, 1:].tolist(), orientation="vertical",
                         fraction=0.018, pad=0.02, label=f"diff ({units})")
        fig.suptitle(f"ELM {fname} ({units}) on native grid (mosaic) — {spatial['target_date']}",
                     fontsize=12, fontweight="bold")
        return fig

    if kind == "z":
        panels = [(f"Actual {tyear} (z)", spatial["z"]["actual"])]
        if "era5" in spatial["z"]:
            m = metrics["era5"]
            panels.append((f"ERA5 {tyear} (z)\nd_box={m['d_box']:.2f}, r_box={m['r_box']:.2f}, r_grid={m['r_grid']:.2f}",
                           spatial["z"]["era5"]))
        for b in ("obs", "clim"):
            if b in spatial["z"]:
                m = metrics[b]
                panels.append((
                    f"{b} ({byears[b]}) (z)\nd_box={m['d_box']:.2f}, r_box={m['r_box']:.2f}, r_grid={m['r_grid']:.2f}",
                    spatial["z"][b],
                ))
        fig, axes = plt.subplots(
            1, len(panels), figsize=(4.3 * len(panels), 3.6),
            subplot_kw={"projection": proj}, constrained_layout=True, squeeze=False,
        )
        coll = None
        for ax, (title, arr) in zip(axes[0], panels):
            ax.set_facecolor("0.92")
            coll = mosaic.polypcolor(
                ax, elm_descriptor, np.ma.masked_invalid(arr),
                cmap=diff_cmap, vmin=-z_lim, vmax=z_lim, antialiaseds=False,
            )
            ax.coastlines(linewidth=0.5, color="0.25")
            ax.set_extent([-180, 180, -60, 85], crs=ccrs.PlateCarree())
            ax.set_title(title, fontsize=9.5)
        fig.colorbar(coll, ax=axes[0, :].tolist(), orientation="horizontal",
                     fraction=0.055, pad=0.04, label="standardized anomaly z")
        fig.suptitle(f"ELM {fname} standardized anomalies (native grid via mosaic) — {spatial['target_date']}",
                     fontsize=12, fontweight="bold")
        return fig

    raise ValueError(f"Unknown kind={kind!r}")


def plot_elm_four_seasons(
    field_name: str,
    test_dates: Sequence[str],
    plans: Mapping[tuple[str, str], SelectionPlan],
    selections: Mapping[tuple[str, str], Mapping[str, int]],
    elm_descriptor: object,
    *,
    box: tuple[float, float] = (30.0, 60.0),
    cache_dir: Path | None = None,
    z_lim: float = 2.5,
):
    """Four-season comparison of ELM standardized anomalies on the native mesh via ``mosaic``."""
    ensure_proj_data()
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt
    import mosaic

    fig, axes = plt.subplots(
        len(test_dates), 4, figsize=(16.5, 2.9 * len(test_dates)),
        subplot_kw={"projection": elm_descriptor.projection}, constrained_layout=True,
    )
    coll = None
    for row, date in enumerate(test_dates):
        byears = {b: selections[(b, date)]["elm"] for b in BRANCHES if (b, date) in selections}
        sp = prepare_spatial_field(
            plans[("clim", date)], "elm", field_name, 2015, byears,
            box=box, cache_dir=cache_dir, log=lambda *_: None,
        )
        m = sp["metrics"]
        row_panels = [
            (f"{date}: Actual 2015 (z)", sp["z"]["actual"]),
            (f"ERA5 2015 (r_box={m['era5']['r_box']:.2f}, r_grid={m['era5']['r_grid']:.2f})", sp["z"]["era5"]),
            (f"obs {byears['obs']} (d_box={m['obs']['d_box']:.2f}, r_box={m['obs']['r_box']:.2f})", sp["z"]["obs"]),
            (f"clim {byears['clim']} (d_box={m['clim']['d_box']:.2f}, r_box={m['clim']['r_box']:.2f})", sp["z"]["clim"]),
        ]
        for col, (title, arr) in enumerate(row_panels):
            ax = axes[row, col]
            ax.set_facecolor("0.92")
            coll = mosaic.polypcolor(
                ax, elm_descriptor, np.ma.masked_invalid(arr),
                cmap="RdBu_r", vmin=-z_lim, vmax=z_lim, antialiaseds=False,
            )
            ax.coastlines(linewidth=0.45, color="0.25")
            ax.set_extent([-180, 180, -60, 85], crs=ccrs.PlateCarree())
            ax.set_title(title, fontsize=9)
    fig.colorbar(coll, ax=axes.ravel().tolist(), orientation="horizontal",
                 fraction=0.025, pad=0.02, label="standardized anomaly z")
    fig.suptitle(f"ELM {field_name}: 4-season standardized anomalies (native grid via mosaic)",
                 fontsize=13, fontweight="bold")
    return fig

