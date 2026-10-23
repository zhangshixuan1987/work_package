"""Saved outputs of the RTENSO2026 initial-condition workflows, for comparing the methods.

Each workflow notebook (01 climatology-match, 02 observation-match, 03
climatology plus observed anomaly) saves the same set of outputs under
``<diag_root>/<method>/data/``:

==========================  ===========================================================
File                        Content
==========================  ===========================================================
``config.json``             the notebook settings (years, box, mode, observations, ...)
``selection.csv``           per date and component: the year used (selected / base year)
``states_YYYYMMDD.npz``     the initial-condition key fields at model points (ELM
                            columns, MPAS-Seaice cells), one array per component/field
``<table>.csv``             tables written with ``ic_eval.save_table``, the same in every
                            notebook: ``input_summary``, ``distances``,
                            ``state_distance`` (D vs the observed anomaly and vs the
                            JRA55-FOSIRL 2015 IC, from ``evaluate_states``), ``field_similarity``,
                            ``global_summary``
==========================  ===========================================================

``load_table`` and ``load_states`` read them back for any set of methods, so
diagnostics comparing the methods can be built without rerunning the
workflows. The states are standardized only when compared: use
``clim_restart.make_plan`` with the same candidate and climatology years to
put them on a common scale.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

import ic_eval
from anom_restart import _spread, _wrms, state_distance
from clim_restart import (
    KEY_FIELDS,
    SelectionPlan,
    _apply_box,
    _columns_to_gridcells,
    comparison_points,
    era5_available,
    load_key_fields,
    load_obs_fields,
    source_file,
)

METHODS = ("clim", "obs", "hybrid_ic")   # default method directories, in workflow order


def data_dir(diag_root: str | Path, method: str) -> Path:
    path = Path(diag_root) / method / "data"
    path.mkdir(parents=True, exist_ok=True)
    return path


def year_states(plan: SelectionPlan, years: Mapping[str, int], cache_dir: Path | None = None,
                components: Sequence[str] = ("elm", "mpassi")) -> dict:
    """Key fields ``{component: {field: values}}`` of the archive year ``years[component]``."""
    out = {}
    for comp in components:
        if comp in plan.components and comp in years:
            fields = load_key_fields(plan, comp, cache_dir, years=[int(years[comp])], log=lambda *_: None)
            out[comp] = {name: values[0] for name, (values, _) in fields.items()}
    return out


def save_config(diag_root: str | Path, method: str, config: Mapping) -> Path:
    path = data_dir(diag_root, method) / "config.json"
    path.write_text(json.dumps(dict(config), indent=2, default=str))
    return path


def save_states(diag_root: str | Path, method: str, date: str, states: Mapping[str, Mapping[str, np.ndarray]],
                years: Mapping[str, int]) -> Path:
    """Write ``states_YYYYMMDD.npz`` (keys ``<component>|<field>``, plus the years used)."""
    arrays = {f"{comp}|{name}": np.asarray(v, dtype=np.float32)
              for comp, fields in states.items() for name, v in fields.items()}
    arrays.update({f"year|{comp}": np.int64(y) for comp, y in years.items()})
    path = data_dir(diag_root, method) / f"states_{date.replace('-', '')}.npz"
    np.savez_compressed(path, **arrays)
    return path


def load_states(diag_root: str | Path, method: str, date: str) -> tuple[dict, dict]:
    """``({component: {field: values}}, {component: year})`` saved by ``save_states``."""
    states, years = {}, {}
    with np.load(Path(diag_root) / method / "data" / f"states_{date.replace('-', '')}.npz") as z:
        for key in z.files:
            kind, name = key.split("|", 1)
            if kind == "year":
                years[name] = int(z[key])
            else:
                states.setdefault(kind, {})[name] = z[key].astype(np.float64)
    return states, years


def load_table(diag_root: str | Path, name: str, methods: Sequence[str] = METHODS) -> pd.DataFrame:
    """One saved table for several methods, with a ``method`` column (missing methods are skipped)."""
    parts = []
    for method in methods:
        path = Path(diag_root) / method / "data" / f"{name}.csv"
        if path.is_file():
            parts.append(pd.read_csv(path).assign(method=method))
    if not parts:
        raise FileNotFoundError(f"No {name}.csv under {diag_root}/{{{','.join(methods)}}}/data")
    return pd.concat(parts, ignore_index=True)


def evaluate_states(
    plan: SelectionPlan,
    states: Mapping[str, Mapping[str, Mapping[str, np.ndarray]]],
    *,
    box: tuple[float, float] | None,
    cache_dir: Path | None = None,
    source: str | None = None,
) -> pd.DataFrame:
    """Joint and field distances of initial-condition states, the common score of all methods.

    ``states[label][component][field]`` are key fields at model points. Each is
    scored (``anom_restart.state_distance``, observable fields) against the
    observed anomaly on the target date (when the plan's observations cover
    it) and against the actual model state of the target year (when that
    restart exists and is not a candidate or climatology year, i.e. 2015), at
    ``box`` and at model points. Lower is closer. The target-year reference
    is labelled ``"<source> <year>"`` (default ``ic_eval.IC_SOURCE``).
    """
    targets = []
    comps = [c for c in plan.components if plan.obs_fields_for(c)]
    if comps and era5_available(plan.target_year, plan.mmdd, comps, obs_root=plan.obs_roots,
                                obs_fields=plan.obs_fields):
        targets.append(("observed", "obs"))
    year = plan.target_year
    if year not in plan.years and all(source_file(plan.src_root, plan.src_prefix, c, year, plan.mmdd).is_file()
                                      for c in plan.components):
        targets.append((f"{source or ic_eval.IC_SOURCE} {year}", year))
    parts = []
    for versus, target in targets:
        for b in (box, None):
            table = state_distance(plan, states, target=target, box=b, cache_dir=cache_dir, log=lambda *_: None)
            parts.append(table.assign(target_date=plan.target_date, versus=versus,
                                      scale="points" if b is None else f"{b[0]:g}x{b[1]:g}"))
    if not parts:
        return pd.DataFrame()
    return pd.concat(parts).rename_axis("state").reset_index()


def reference_state(plan: SelectionPlan, cache_dir: Path | None = None) -> dict:
    """``{label: state}`` of the archive state of the target year (JRA55-FOSIRL 2015 for the 2015 dates);
    empty when that restart does not exist or is a candidate or climatology year."""
    year = plan.target_year
    if year in plan.years or not all(source_file(plan.src_root, plan.src_prefix, c, year, plan.mmdd).is_file()
                                     for c in plan.components):
        return {}
    return {f"{ic_eval.IC_SOURCE} {year}": year_states(plan, {"elm": year, "mpassi": year}, cache_dir)}


def observed_similarity(
    plan: SelectionPlan,
    states: Mapping[str, Mapping[str, Mapping[str, np.ndarray]]],
    *,
    box: tuple[float, float] | None,
    cache_dir: Path | None = None,
) -> pd.DataFrame:
    """Field-by-field similarity of initial-condition states to the observed anomaly on the target date.

    Standardized as ``evaluate_states`` (model and observed anomalies each from
    their own climatology-year mean and spread), at ``box`` and at model
    points: ``r`` is the area-weighted pattern correlation and ``d`` the RMS
    difference of the standardized anomalies (``d`` equals the field distance
    of ``evaluate_states``); ``rmse`` is ``d`` in model units (``d`` times the
    field's interannual spread ``S`` at that scale, the observed anomaly
    expressed in model units). One row per state, field, and scale; observable
    fields only.
    """
    rows = []
    years, ref = plan.years, plan.clim_rows
    ny = len(years)
    labels = list(states)
    for comp in plan.components:
        names = list(plan.obs_fields_for(comp))
        if comp not in KEY_FIELDS or not names:
            continue
        model = load_key_fields(plan, comp, cache_dir, years=years, log=lambda *_: None)
        points = comparison_points(plan, comp, cache_dir)
        obs = load_obs_fields(plan, comp, years + [plan.target_year], cache_dir, lambda *_: None)
        for name in names:
            values, weight = model[name]
            stack = np.vstack([values] + [np.asarray(states[l][comp][name])[None, :] for l in labels])
            if comp == "elm":
                stack, weight = _columns_to_gridcells(stack, weight, points["col_grid"], points["lat"].size)
            for b in (box, None):
                agg, wb, o = _apply_box(stack, weight, points["lat"], points["lon"], b, obs=obs[name])
                valid = np.all(np.isfinite(agg), axis=0) & np.all(np.isfinite(o), axis=0)
                w = np.where(valid, wb, 0.0)
                agg, o = np.where(valid, agg, 0.0), np.where(valid, o, 0.0)
                mean, S = agg[ref].mean(axis=0), _spread(agg[ref], w)
                obar, S_o = o[ref].mean(axis=0), _spread(o[ref], w)
                zt = (o[-1] - obar) / S_o if S_o > 0 else np.zeros_like(obar)
                for j, label in enumerate(labels):
                    z = (agg[ny + j] - mean) / S if S > 0 else np.zeros_like(mean)
                    rows.append({"target_date": plan.target_date, "state": label, "component": comp,
                                 "field": name, "scale": "points" if b is None else f"{b[0]:g}x{b[1]:g}",
                                 "r": ic_eval._weighted_corr(z, zt, w), "d": _wrms(z - zt, w),
                                 "S": S, "rmse": _wrms(z - zt, w) * S,
                                 "units": ic_eval.FIELD_UNITS.get(name, "")})
    return pd.DataFrame(rows)
