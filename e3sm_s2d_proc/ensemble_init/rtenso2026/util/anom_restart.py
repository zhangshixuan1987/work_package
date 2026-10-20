"""Climatology plus observed anomaly ("clim_anom") restarts for RTENSO2026.

The ``clim`` branch of ``clim_restart`` selects the BruteForce year closest to
the candidate-year climatology. This module starts from that restart (the
*base*) and adds the observed (ERA5) anomaly on the target date to its key
fields, so the state is the model climatology plus the observed anomaly:

    target_i = xbar_i + a_i,      a_i = z_i * s_i,      z_i = (o_T,i - obar_i) / so_i

with ``xbar`` / ``obar`` the model / observed means over the candidate years at
point ``i`` (ELM grid cell or MPAS-Seaice cell), ``o_T`` the observed value on
the target date, and ``s_i`` / ``so_i`` the model / observed interannual
standard deviations at the point (``scale="point"``, the default): the observed
standardized anomaly is given the model's own variability, so model-ERA5 unit,
depth, and variance differences cancel. ``scale="field"`` uses the
field-wide ratio ``a_i = (o_T,i - obar_i) * S_model / S_obs`` (``S`` the
area-weighted RMS spread); ``scale="none"`` the raw observed anomaly (only
meaningful for matching units, such as ``ci`` and ``stl1``). In every mode
``|a_i| <= z_max * s_i``: this keeps out ERA5 artifacts such as the snow that
accumulates to its 10 m cap over high mountains and the Antarctic coast. With
``anom_box = (dlat, dlon)`` the anomaly is then averaged over lat-lon boxes.

The target is applied to the base restart as an increment that keeps the
restart internally consistent:

=================================  ======  ============================================
Key field                          Rule    Restart variables changed
=================================  ======  ============================================
NH/SH ice concentration            ratio   iceAreaCategory, iceVolumeCategory,
                                           snowVolumeCategory (all categories, same
                                           factor: category thickness, snow depth and
                                           the per-volume enthalpy/salinity tracers
                                           are unchanged)
snow water equivalent              ratio   H2OSNO, SNOW_DEPTH, INT_SNOW, snow-layer
                                           H2OSOI_LIQ/H2OSOI_ICE, DZSNO, ZSNO, ZISNO
                                           (snow density unchanged)
surface soil water (layers 1-3)    ratio   H2OSOI_LIQ, H2OSOI_ICE (liquid/ice split
                                           unchanged)
subsurface soil water (layers 4+)  ratio   H2OSOI_LIQ, H2OSOI_ICE
surface soil temperature           shift   T_SOISNO, soil layers 1..``soil_temp_layers``
subsurface soil temperature        shift   T_SOISNO, soil layers 4-10
=================================  ======  ============================================

Soil increments are continuous in depth (``continuous_profile=True``): the soil-water factor and the
soil-temperature shift vary linearly from their surface value to their deeper value (centred at about
1.49 m), solved per column so that each key field still gets exactly its target, and the temperature
shift decays to zero below 2.89 m (the bottom of ERA5 stl4/swvl4) with e-folding depth
``temp_decay_depth``. The base profile is continuous, so the adjusted profile is too.

``ratio`` fields are multiplied by ``target / base`` (per grid cell, bounded by
``ratio_bounds``); ``shift`` fields get ``target - base`` added. A ratio cannot
create snow or sea ice where the base has none, so those points keep the base
state (reported as ``frac unreachable``). With ``envelope=True`` the new value
of each column / cell is kept within the range of the candidate years at that
point (ELM has no porosity in its restart, so this is the bound that keeps soil
water physical). Only soil and crop columns (``cols1d_ityplun`` 1, 2) of ELM
are changed; glacier, lake, wetland, and urban columns keep the base state.
The water table, aquifer, river (MOSART, copied from the base land year), and
phase partition are not adjusted; the model adjusts them in its first steps.

Source files are only read; the base restart is copied to the output
directory and modified there.
"""

from __future__ import annotations

import json
import os
import shutil
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Mapping, Sequence

import netCDF4 as nc
import numpy as np
import pandas as pd

from clim_restart import (
    COMPONENTS,
    ELM_SOIL_INTERFACES,
    ELM_TSUB_WEIGHTS,
    KEY_FIELDS,
    SPVAL_THRESHOLD,
    SelectionPlan,
    _apply_box,
    _columns_to_gridcells,
    _model_points,
    _obs_index,
    _obs_slice,
    _read,
    _to_boxes,
    comparison_points,
    load_key_fields,
    load_obs_fields,
    read_key_fields,
)

# How each observable key field is applied to the base restart.
RULES = {
    "mpassi": {"NH ice concentration": "ratio", "SH ice concentration": "ratio"},
    "elm": {
        "snow water equivalent": "ratio",
        "surface soil water": "ratio",
        "subsurface soil water": "ratio",
        "surface soil temperature": "shift",
        "subsurface soil temperature": "shift",
    },
}
# Physical bounds of the key fields (applied to the target).
BOUNDS = {
    "NH ice concentration": (0.0, 1.0),
    "SH ice concentration": (0.0, 1.0),
    "snow water equivalent": (0.0, np.inf),
    "surface soil water": (0.0, np.inf),
    "subsurface soil water": (0.0, np.inf),
}
ELM_SOIL_LANDUNITS = (1, 2)      # istsoil, istcrop
EPS = 1.0e-6

# Vertical profile of the soil increments (continuous_profile=True). The increment goes linearly in depth
# from its surface value (at Z_TOP) to its deeper value (at Z_DEEP, the thickness-weighted centre of the
# subsurface-temperature layers, about 1.49 m) and stays there down to Z_OBS_BOTTOM (2.89 m, the bottom of
# ERA5 stl4/swvl4); the temperature increment then decays with depth to zero (e-folding TEMP_DECAY_DEPTH).
ELM_SOIL_NODES = 0.025 * (np.exp(0.5 * (np.arange(1, 16) - 0.5)) - 1)
Z_DEEP = float(ELM_TSUB_WEIGHTS @ ELM_SOIL_NODES[3:10] / ELM_TSUB_WEIGHTS.sum())
Z_OBS_BOTTOM = 2.89


@dataclass
class AnomalySettings:
    scale: str = "point"                         # "point", "field", or "none"
    z_max: float = 3.0                           # |anomaly| <= z_max model std devs at each point
    anom_box: tuple[float, float] | None = None  # smooth the anomaly over lat-lon boxes
    envelope: bool = True                        # keep each point within the candidate-year range
    ratio_bounds: tuple[float, float] = (0.25, 4.0)
    soil_temp_layers: int = 3                    # soil layers shifted by the surface increment (stepwise only)
    continuous_profile: bool = True              # smooth soil increments in depth (False: one value per layer group)
    temp_decay_depth: float = 3.0                # e-folding depth (m) of the temperature increment below 2.89 m
    fields: dict[str, list[str]] = field(default_factory=lambda: {c: list(f) for c, f in RULES.items()})

    def __post_init__(self):
        if self.scale not in ("point", "field", "none"):
            raise ValueError(f"scale must be 'point', 'field', or 'none', got {self.scale!r}")
        for comp, names in self.fields.items():
            unknown = set(names) - set(RULES.get(comp, {}))
            if unknown:
                raise ValueError(f"No anomaly rule for {comp} fields {sorted(unknown)}; use {list(RULES[comp])}")


def _wrms(x: np.ndarray, w: np.ndarray) -> float:
    """Area-weighted RMS of ``x`` over points."""
    return float(np.sqrt((w * x ** 2).sum() / w.sum())) if w.sum() > 0 else 0.0


def _spread(values: np.ndarray, w: np.ndarray) -> float:
    """Area-weighted RMS interannual spread S of ``values[n_years, n_points]``."""
    return _wrms(np.sqrt(((values - values.mean(axis=0)) ** 2).mean(axis=0)), w)


def _elm_soil_mask(plan: SelectionPlan, year: int) -> np.ndarray:
    with nc.Dataset(plan.source("elm", year), "r") as ds:
        return np.isin(_read(ds, "cols1d_ityplun"), ELM_SOIL_LANDUNITS) & (_read(ds, "cols1d_active") > 0)


def compute_increments(
    plan: SelectionPlan,
    base_years: Mapping[str, int],
    settings: AnomalySettings | None = None,
    *,
    cache_dir: Path | None = None,
    log=print,
) -> tuple[dict[str, dict], pd.DataFrame]:
    """Per-point factors / shifts that turn the base restart into clim + anomaly.

    ``plan`` must be an ``obs`` plan (observations and candidate years);
    ``base_years`` gives the base year per component (the ``clim`` selection).

    Returns ``(increments, stats)``: ``increments[comp][field]`` holds the
    rule, the per-point ``factor`` (ratio) or ``shift`` and the new key-field
    values ``new`` (model points: ELM columns, MPAS-Seaice cells); ``stats``
    summarizes each field (anomaly scale, RMS target anomaly, achieved
    correlation and RMS ratio of the realized anomaly, clipping).
    """
    settings = settings or AnomalySettings()
    increments: dict[str, dict] = {}
    rows = []
    for comp in plan.components:
        names = [f for f in settings.fields.get(comp, []) if f in plan.obs_fields_for(comp)]
        if not names:
            continue
        years = plan.candidate_years
        if base_years[comp] not in years:
            raise ValueError(f"Base year {base_years[comp]} for {comp} is not a candidate year")
        ib = years.index(base_years[comp])
        model = load_key_fields(plan, comp, cache_dir, log=log)
        obs = load_obs_fields(plan, comp, years + [plan.target_year], cache_dir, log)
        points = comparison_points(plan, comp, cache_dir)
        soil = _elm_soil_mask(plan, base_years[comp]) if comp == "elm" else None
        increments[comp] = {}
        for name in names:
            rule = RULES[comp][name]
            X, w = model[name]                           # [ny, npoint], model points
            if comp == "elm":
                col_grid = points["col_grid"]
                Xg, wg = _columns_to_gridcells(X, np.where(soil, w, 0.0), col_grid, points["lat"].size)
            else:
                col_grid = np.arange(X.shape[1])
                Xg, wg = X, w
            O = obs[name]
            valid = (wg > 0) & np.all(np.isfinite(Xg), axis=0) & np.all(np.isfinite(O), axis=0)
            wv = np.where(valid, wg, 0.0)
            Xg = np.where(valid, Xg, 0.0)
            O = np.where(valid, O, 0.0)

            clim = Xg.mean(axis=0)
            a_obs = O[-1] - O[:-1].mean(axis=0)
            s_mod, s_obs = Xg.std(axis=0), O[:-1].std(axis=0)
            S_mod, S_obs = _spread(Xg, wv), _spread(O[:-1], wv)
            if settings.scale == "point":
                ok = s_obs > 1e-3 * S_obs
                z = np.where(ok, a_obs / np.where(ok, s_obs, 1.0), 0.0)
                anom = np.clip(z, -settings.z_max, settings.z_max) * s_mod
                scale = S_mod / S_obs if S_obs > 0 else 0.0   # reported only
            else:
                scale = (S_mod / S_obs if S_obs > 0 else 0.0) if settings.scale == "field" else 1.0
                anom = np.clip(a_obs * scale, -settings.z_max * s_mod, settings.z_max * s_mod)
            if settings.anom_box is not None:
                boxed, _, index = _to_boxes(np.where(valid, anom, np.nan)[None, :], wv,
                                            points["lat"], points["lon"], settings.anom_box)
                anom = np.where(valid, np.nan_to_num(boxed[0][index]), 0.0)

            lo, hi = BOUNDS.get(name, (-np.inf, np.inf))
            target = np.clip(clim + anom, lo, hi)
            base_g = Xg[ib]
            if rule == "ratio":
                inc_g = np.where(valid & (base_g > EPS), target / np.where(base_g > EPS, base_g, 1.0), 1.0)
                inc_g = np.clip(inc_g, *settings.ratio_bounds)
            else:
                inc_g = np.where(valid, target - base_g, 0.0)

            # Model points: apply the grid-cell increment, then the envelope.
            x_b = X[ib]
            change = valid[col_grid] & np.isfinite(x_b)
            if soil is not None:
                change &= soil
            new = np.where(change, x_b * inc_g[col_grid] if rule == "ratio" else x_b + inc_g[col_grid], x_b)
            new = np.clip(new, lo, hi)
            clipped = np.zeros(new.shape, bool)
            if settings.envelope:
                env_lo, env_hi = np.nanmin(X, axis=0), np.nanmax(X, axis=0)
                bounded = np.clip(new, env_lo, env_hi)
                clipped = change & (np.abs(bounded - new) > 1e-9 * (1 + np.abs(new)))
                new = np.where(change, bounded, x_b)
            if rule == "ratio":
                inc = np.where(change & (x_b > EPS), new / np.where(x_b > EPS, x_b, 1.0), 1.0)
                new = np.where(change & (x_b > EPS), new, x_b)
            else:
                inc = np.where(change, new - x_b, 0.0)
            increments[comp][name] = {"rule": rule, "factor" if rule == "ratio" else "shift": inc,
                                      "new": new, "scale": scale}

            # Realized vs target anomaly on the grid cells.
            if comp == "elm":
                new_g, _ = _columns_to_gridcells(new[None, :], np.where(soil, w, 0.0), col_grid,
                                                 points["lat"].size)
                new_g = np.where(valid, new_g[0], 0.0)
            else:
                new_g = np.where(valid, new, 0.0)
            realized, wanted = new_g - clim, target - clim
            unreachable = valid & (rule == "ratio") & (base_g <= EPS) & (target > EPS)
            increments[comp][name]["maps"] = {      # comparison points (ELM grid cells, sea-ice cells)
                "lat": points["lat"], "lon": points["lon"], "valid": valid, "weight": wv,
                "target anomaly": np.where(valid, wanted, np.nan),
                "base anomaly": np.where(valid, base_g - clim, np.nan),
                "realized anomaly": np.where(valid, realized, np.nan),
            }
            rows.append({
                "component": comp, "field": name, "rule": rule, "base year": base_years[comp],
                "anomaly scale S_model/S_obs": scale,
                "RMS target anomaly": _wrms(wanted, wv),
                "RMS base anomaly": _wrms(base_g - clim, wv),
                "RMS realized anomaly": _wrms(realized, wv),
                "corr(realized, target)": float(np.corrcoef(realized[valid], wanted[valid])[0, 1])
                if valid.sum() > 2 and wanted[valid].std() > 0 and realized[valid].std() > 0 else np.nan,
                "frac changed": float(w[change & (np.abs(new - x_b) > 0)].sum() / max(w[change].sum(), 1e-30)),
                "frac clipped (envelope)": float(w[clipped].sum() / max(w[change].sum(), 1e-30)),
                "frac unreachable (no base ice/snow)": float(wv[unreachable].sum() / max(wv.sum(), 1e-30)),
            })
    return increments, pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------

def _scale_var(ds: nc.Dataset, name: str, factor: np.ndarray, layers: slice = slice(None)) -> None:
    """Multiply ``name`` by a per-point ``factor`` (points are the first axis after
    an optional ``Time`` axis; ``layers`` selects along the next axis)."""
    var = ds.variables[name]
    var.set_auto_maskandscale(False)
    data = np.asarray(var[...])
    view = data[0] if var.dimensions[0] == "Time" else data
    sub = view[:, layers] if view.ndim > 1 else view
    f = factor.reshape(factor.shape + (1,) * (sub.ndim - 1))
    sub[...] = np.where(np.abs(sub) < SPVAL_THRESHOLD, sub * f, sub)   # sub is a view of data
    var[...] = data


def _ramp(z_top: float) -> np.ndarray:
    """Weight of the deeper value at each soil node: 0 above ``z_top``, 1 below ``Z_DEEP``, linear between."""
    return np.clip((ELM_SOIL_NODES - z_top) / (Z_DEEP - z_top), 0.0, 1.0)


def _soil_water_factors(base: np.ndarray, r_top: np.ndarray, r_sub: np.ndarray,
                        bounds: tuple[float, float]) -> np.ndarray:
    """Per-layer factors (ncol, 15) for soil water that vary smoothly in depth.

    f_j = A (1 - phi_j) + B phi_j with A, B solved per column so that the totals over layers 1-3 and 4+
    change by ``r_top`` and ``r_sub`` (the key-field ratios). Falls back to the stepwise factors where the
    system is degenerate, and clips to ``bounds``.
    """
    phi = _ramp(ELM_SOIL_NODES[1])                       # surface value at the centre of layers 1-3
    top, sub = slice(0, 3), slice(3, None)
    b1, b2 = base[:, top].sum(1), base[:, sub].sum(1)
    p1, q1 = (base[:, top] * (1 - phi[top])).sum(1), (base[:, top] * phi[top]).sum(1)
    p2, q2 = (base[:, sub] * (1 - phi[sub])).sum(1), (base[:, sub] * phi[sub]).sum(1)
    det = p1 * q2 - q1 * p2
    ok = (b1 > EPS) & (b2 > EPS) & (np.abs(det) > EPS * np.maximum(b1 * b2, EPS))
    det_s = np.where(ok, det, 1.0)
    a = np.where(ok, (r_top * b1 * q2 - q1 * r_sub * b2) / det_s, 1.0)
    b = np.where(ok, (p1 * r_sub * b2 - p2 * r_top * b1) / det_s, 1.0)
    f = a[:, None] * (1 - phi) + b[:, None] * phi
    step = np.concatenate([np.repeat(r_top[:, None], 3, axis=1), np.repeat(r_sub[:, None], 12, axis=1)], axis=1)
    f = np.where(ok[:, None], f, step)
    return np.clip(f, *bounds)


def _soil_temp_shifts(d_top: np.ndarray, d_sub: np.ndarray, decay: float) -> np.ndarray:
    """Per-layer temperature shifts (ncol, 15) that vary smoothly in depth.

    g_j = A (1 - phi_j) + B phi_j above 2.89 m, decaying as B exp(-(z - 2.89) / decay) below, with A = d_top
    (soil layer 1 gets the surface increment) and B solved so that the thickness-weighted mean over layers
    4-10 (the subsurface key field) gets ``d_sub``.
    """
    phi = _ramp(ELM_SOIL_NODES[0])
    taper = np.where(ELM_SOIL_NODES > Z_OBS_BOTTOM, np.exp(-(ELM_SOIL_NODES - Z_OBS_BOTTOM) / decay), 1.0)
    w = np.zeros(15)
    w[3:10] = ELM_TSUB_WEIGHTS
    shape_b = phi * taper                       # contribution of B to each layer
    shape_a = (1 - phi) * taper
    a = d_top
    b = (d_sub * w.sum() - a * (w @ shape_a)) / (w @ shape_b)
    bound = 2.0 * np.maximum(np.abs(d_top), np.abs(d_sub))   # guard against extreme profiles
    b = np.clip(b, -bound, bound)
    return a[:, None] * shape_a + b[:, None] * shape_b


def _apply_elm(ds: nc.Dataset, incs: Mapping[str, dict], settings: "AnomalySettings") -> None:
    k = len(ds.dimensions["levsno"])
    ncol = len(ds.dimensions["column"])
    if "snow water equivalent" in incs:
        f = incs["snow water equivalent"]["factor"]
        for name in ("H2OSNO", "SNOW_DEPTH", "INT_SNOW"):
            if name in ds.variables:
                _scale_var(ds, name, f)
        for name in ("DZSNO", "ZSNO", "ZISNO"):
            _scale_var(ds, name, f)
        for name in ("H2OSOI_LIQ", "H2OSOI_ICE"):
            _scale_var(ds, name, f, layers=slice(0, k))

    one, zero = np.ones(ncol), np.zeros(ncol)
    r_top = incs.get("surface soil water", {}).get("factor", one)
    r_sub = incs.get("subsurface soil water", {}).get("factor", one)
    d_top = incs.get("surface soil temperature", {}).get("shift", zero)
    d_sub = incs.get("subsurface soil temperature", {}).get("shift", zero)

    def read(name):
        var = ds.variables[name]
        var.set_auto_maskandscale(False)
        return var, np.asarray(var[...])

    if "surface soil water" in incs or "subsurface soil water" in incs:
        (vl, liq), (vi, ice) = read("H2OSOI_LIQ"), read("H2OSOI_ICE")
        sl, si = liq[:, k:k + 15], ice[:, k:k + 15]
        gl, gi = np.abs(sl) < SPVAL_THRESHOLD, np.abs(si) < SPVAL_THRESHOLD
        if settings.continuous_profile:
            base = np.where(gl, sl, 0.0) + np.where(gi, si, 0.0)
            f = _soil_water_factors(base, r_top, r_sub, settings.ratio_bounds)
        else:
            f = np.concatenate([np.repeat(r_top[:, None], 3, 1), np.repeat(r_sub[:, None], 12, 1)], axis=1)
        liq[:, k:k + 15] = np.where(gl, sl * f, sl)
        ice[:, k:k + 15] = np.where(gi, si * f, si)
        vl[...] = liq
        vi[...] = ice

    if "surface soil temperature" in incs or "subsurface soil temperature" in incs:
        vt, t = read("T_SOISNO")
        st = t[:, k:k + 15]
        if settings.continuous_profile:
            g = _soil_temp_shifts(d_top, d_sub, settings.temp_decay_depth)
        else:
            g = np.zeros((ncol, 15))
            g[:, :settings.soil_temp_layers] = d_top[:, None]
            g[:, 3:10] = d_sub[:, None]
        t[:, k:k + 15] = np.where(np.abs(st) < SPVAL_THRESHOLD, st + g, st)
        vt[...] = t


def _apply_mpassi(ds: nc.Dataset, incs: Mapping[str, dict]) -> None:
    factor = np.ones(len(ds.dimensions["nCells"]))
    lat = _read(ds, "latCell")
    for name, inc in incs.items():
        hemi = lat > 0 if name.startswith("NH") else lat < 0
        factor = np.where(hemi, inc["factor"], factor)
    for name in ("iceAreaCategory", "iceVolumeCategory", "snowVolumeCategory"):
        _scale_var(ds, name, factor)


def write_anomaly_restarts(
    plan: SelectionPlan,
    base_years: Mapping[str, int],
    increments: Mapping[str, dict],
    stats: pd.DataFrame,
    settings: AnomalySettings,
    *,
    overwrite: bool = False,
    log=print,
) -> dict:
    """Copy each base restart to ``plan.out_dir`` and apply the increments.

    The river follows the land base year and is copied unchanged. Writes the
    rpointers and ``anomaly.json`` (plan, base years, settings, stats).
    """
    out_dir = Path(plan.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written = {}
    for comp in plan.components:
        year = base_years[comp]
        src, dst = plan.source(comp, year), plan.output(comp)
        if dst.exists() and not overwrite:
            raise FileExistsError(f"{dst} exists; pass overwrite=True to replace it")
        t0 = time.time()
        tmp = dst.with_name(f".{dst.name}.tmp")
        shutil.copyfile(src, tmp)
        if comp in increments:
            with nc.Dataset(tmp, "r+") as ds:
                if comp == "elm":
                    _apply_elm(ds, increments[comp], settings)
                elif comp == "mpassi":
                    _apply_mpassi(ds, increments[comp])
                ds.setncattr("rtenso2026_clim_anom",
                             f"base year {year} + observed anomaly on {plan.target_date} "
                             f"({', '.join(increments[comp])})")
        os.replace(tmp, dst)
        log(f"  {comp}: {src.name} {'+ anomaly ' if comp in increments else ''}-> {dst} "
            f"({time.time() - t0:.0f} s)")
        info = COMPONENTS[comp]
        pointer = f"./{dst.name}" if info["pointer"] == "file" else f"{plan.target_date}_00:00:00"
        (out_dir / info["rpointer"]).write_text(pointer + "\n")
        written[comp] = {"base year": int(year), "source": str(src), "output": str(dst),
                         "fields": list(increments.get(comp, {}))}

    summary = {
        "plan": asdict(plan),
        "method": "base restart (clim year) + observed anomaly on the target date",
        "settings": asdict(settings),
        "restarts": written,
        "stats": json.loads(stats.round(4).to_json(orient="records")),
    }
    (out_dir / "anomaly.json").write_text(json.dumps(summary, indent=2, default=str))
    return summary


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

def state_distance(
    plan: SelectionPlan,
    states: Mapping[str, Mapping[str, Mapping[str, np.ndarray]]],
    *,
    target: str | int = "obs",
    box: tuple[float, float] | None = None,
    cache_dir: Path | None = None,
    log=print,
) -> pd.DataFrame:
    """Distance d_f of arbitrary model states from the observed or a model target.

    ``states[label][comp][field]`` holds key-field values at model points (as
    from ``read_key_fields``). Standardization uses the candidate years, as in
    ``clim_restart.compute_distances`` (``target="obs"``) or
    ``distance_to_model_year`` (``target=<year>``), so the numbers compare
    directly with those tables. Only the observable fields are scored.
    """
    rows = {label: {} for label in states}
    years = plan.candidate_years
    for comp in plan.components:
        names = list(plan.obs_fields_for(comp))
        if comp not in KEY_FIELDS or not names:
            continue
        extra = [] if target == "obs" else [int(target)]
        model = load_key_fields(plan, comp, cache_dir, years=years + extra, log=log)
        points = comparison_points(plan, comp, cache_dir)
        if target == "obs":
            obs = load_obs_fields(plan, comp, years + [plan.target_year], cache_dir, log)
        for name in names:
            values, weight = model[name]
            stack = np.vstack([values[:len(years)]] + [states[l][comp][name][None, :] for l in states])
            if target == "obs":
                if comp == "elm":
                    stack, weight = _columns_to_gridcells(stack, weight, points["col_grid"], points["lat"].size)
                agg, wb, o = _apply_box(stack, weight, points["lat"], points["lon"], box, obs=obs[name])
            else:
                lat, lon = _model_points(plan, comp, cache_dir)
                stack = np.vstack([stack, values[-1][None, :]])
                agg, wb, _ = _apply_box(stack, weight, lat, lon, box)
            valid = np.all(np.isfinite(agg), axis=0)
            if target == "obs":
                valid &= np.all(np.isfinite(o), axis=0)
            w = np.where(valid, wb, 0.0)
            agg = np.where(valid, agg, 0.0)
            ny = len(years)
            mean = agg[:ny].mean(axis=0)
            S = _spread(agg[:ny], w)
            if target == "obs":
                o = np.where(valid, o, 0.0)
                obar = o[:-1].mean(axis=0)
                S_o = _spread(o[:-1], w)
                zt = (o[-1] - obar) / S_o if S_o > 0 else np.zeros_like(obar)
            else:
                zt = (agg[-1] - mean) / S if S > 0 else np.zeros_like(mean)
            for j, label in enumerate(states):
                z = (agg[ny + j] - mean) / S if S > 0 else np.zeros_like(mean)
                rows[label][f"{comp}: {name}"] = _wrms(z - zt, w)
    table = pd.DataFrame(rows).T
    for comp in plan.components:
        cols = [c for c in table.columns if c.startswith(f"{comp}: ")]
        if cols:
            table[comp] = np.sqrt((table[cols] ** 2).mean(axis=1))
    comps = [c for c in plan.components if c in table]
    table["joint"] = np.sqrt((table[comps] ** 2).mean(axis=1))
    return table


def read_states(plan: SelectionPlan, files: Mapping[str, Mapping[str, Path]]) -> dict:
    """Key fields ``{label: {comp: {field: values}}}`` of restart files ``{label: {comp: path}}``."""
    out = {}
    for label, paths in files.items():
        out[label] = {}
        for comp, path in paths.items():
            if comp in KEY_FIELDS:
                out[label][comp] = {k: v[0] for k, v in read_key_fields(comp, Path(path)).items()}
    return out


# ---------------------------------------------------------------------------
# Soil profiles: observations and the model before / after the adjustment
# ---------------------------------------------------------------------------

ERA5_SOIL_LEVELS = {     # ERA5 / ERA5-Land soil levels: (top, bottom) in m
    1: (0.0, 0.07), 2: (0.07, 0.28), 3: (0.28, 1.0), 4: (1.0, 2.89),
}


def _in_region(lat: np.ndarray, lon: np.ndarray, box: tuple[float, float, float, float]) -> np.ndarray:
    """Points inside ``box = (lat0, lat1, lon0, lon1)`` (lon in degrees east, -180..360)."""
    lat0, lat1, lon0, lon1 = box
    lon = np.mod(lon, 360.0)
    lo0, lo1 = np.mod(lon0, 360.0), np.mod(lon1, 360.0)
    in_lon = (lon >= lo0) & (lon <= lo1) if lo0 <= lo1 else (lon >= lo0) | (lon <= lo1)
    return (lat >= lat0) & (lat <= lat1) & in_lon


def soil_profiles(
    plan: SelectionPlan,
    base_year: int,
    new_file: Path,
    regions: Mapping[str, tuple[float, float, float, float]],
    *,
    nlayers: int = 10,
    log=print,
) -> dict:
    """Region-mean soil profiles for one date: model climatology, base (``clim``) year, adjusted restart
    (``new_file``), and observations (ERA5-Land or ERA5 ``swvl1-4`` / ``stl1-4``, from ``plan``).

    Model soil water is converted to volumetric water content, (liquid + ice) / (1000 kg m-3 * dz), on the
    standard ELM soil layers, so it compares directly with ERA5 ``swvl``. Region means use soil and crop
    columns weighted by ``cols1d_wtxy * cos(lat)``; observations use the ELM grid cells of the region,
    weighted by the same columns. Climatologies use the plan's candidate years.

    Returns ``{"model": {state: {"theta"|"tsoi": (nreg, nlayers)}}, "obs": {state: {...: (nreg, 4)}},
    "depth": layer nodes, "interfaces", "obs_levels", "regions"}`` with states ``clim``, ``base``, ``new``
    (model) and ``clim``, ``target`` (observations).
    """
    names = list(regions)
    dz = np.diff(ELM_SOIL_INTERFACES)[:nlayers]
    with nc.Dataset(plan.source("elm", base_year), "r") as ds:
        lat, lon = _read(ds, "cols1d_lat"), _read(ds, "cols1d_lon")
        w = (_read(ds, "cols1d_wtxy") * np.cos(np.radians(lat)) * (_read(ds, "cols1d_active") > 0)
             * np.isin(_read(ds, "cols1d_ityplun"), ELM_SOIL_LANDUNITS))
        col_grid = _read(ds, "cols1d_gridcell_index").astype(np.int64) - 1
    masks = np.stack([_in_region(lat, lon, regions[r]) for r in names]) * w        # (nreg, ncol)

    def model_profile(path: Path) -> dict:
        with nc.Dataset(path, "r") as ds:
            k = len(ds.dimensions["levsno"])
            t = _read(ds, "T_SOISNO")[:, k:k + nlayers]
            liq = _read(ds, "H2OSOI_LIQ")[:, k:k + nlayers]
            ice = _read(ds, "H2OSOI_ICE")[:, k:k + nlayers]
        good = (np.abs(t) < SPVAL_THRESHOLD) & (np.abs(liq) < SPVAL_THRESHOLD) & (np.abs(ice) < SPVAL_THRESHOLD)
        theta = np.where(good, (liq + ice) / (1000.0 * dz), 0.0)
        t = np.where(good, t, 0.0)
        out = {}
        for key, arr in (("theta", theta), ("tsoi", t)):
            num = masks @ arr
            den = masks @ good.astype(float)
            out[key] = num / np.where(den > 0, den, np.nan)
        return out

    years = list(plan.candidate_years)
    clim = []
    t0 = time.time()
    for n, y in enumerate(years):
        clim.append(model_profile(plan.source("elm", y)))
        log(f"  model profiles {y} ({n + 1}/{len(years)}, {time.time() - t0:.0f} s)")
    model = {"clim": {key: np.mean([c[key] for c in clim], axis=0) for key in ("theta", "tsoi")},
             "base": model_profile(plan.source("elm", base_year)),
             "new": model_profile(Path(new_file))}

    # Observations on the ELM grid cells of each region, weighted by the region's soil columns.
    root = plan.obs_root_for("elm")
    pts = comparison_points(plan, "elm", None)
    gw = np.stack([np.bincount(col_grid, weights=m, minlength=pts["lat"].size) for m in masks])  # (nreg, ngrid)
    i, j = _obs_index(pts["lat"], pts["lon"], root)

    def obs_profile(year: int) -> dict:
        out = {}
        for key, prefix in (("theta", "swvl"), ("tsoi", "stl")):
            levels = np.stack([_obs_slice(f"{prefix}{lev}", year, plan.mmdd, root)[i, j] for lev in ERA5_SOIL_LEVELS],
                              axis=1)                                                # (ngrid, 4)
            ok = np.isfinite(levels)
            num = gw @ np.where(ok, levels, 0.0)
            den = gw @ ok.astype(float)
            out[key] = num / np.where(den > 0, den, np.nan)
        return out

    obs_years = [obs_profile(y) for y in years]
    obs = {"clim": {key: np.mean([o[key] for o in obs_years], axis=0) for key in ("theta", "tsoi")},
           "target": obs_profile(plan.target_year)}
    log(f"  observation profiles read ({time.time() - t0:.0f} s)")
    return {"model": model, "obs": obs, "regions": names, "depth": ELM_SOIL_NODES[:nlayers],
            "interfaces": ELM_SOIL_INTERFACES[:nlayers + 1], "obs_levels": ERA5_SOIL_LEVELS}
