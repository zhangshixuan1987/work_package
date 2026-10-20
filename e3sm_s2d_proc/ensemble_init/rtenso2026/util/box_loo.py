"""Leave-one-out choice of the box size used by the ``obs`` year selection.

Step 1 of the RTENSO2026 workflow (``01_choose_box_size.ipynb``; before the
search in 02 and the 2015 evaluation in 04): the box size is chosen from the BruteForce archive alone,
so the 2015 test dates stay independent.

For a calendar date, each archive year ``y`` in turn is treated as the unknown
target: the ``obs`` selection (``clim_restart``) picks the year closest to the
observed (ERA5) anomaly of ``y`` from the other years, and the pick is scored
against the actual model state of ``y`` with the ``distance_to_model_year``
metric. Every box size is scored at the same fixed scales (``scores``), so the
comparison does not favour the box that matches the scoring scale.

Skill of a pick: ``(D_rms - D_pick) / (D_rms - D_best)``, with ``D_pick`` the
truth distance of the picked year, ``D_best`` the smallest over the candidates
(perfect pick, skill 1), and ``D_rms`` the RMS over the candidates (a random
pick, skill about 0).
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

import clim_restart as cr

BOXES = [None, (10, 20), (15, 30), (20, 40), (30, 60), (45, 90)]
SCORES = {"points": None, "30x60": (30, 60)}
FIELD_SETS = ("observable", "all")     # observable: fields with ERA5; all: every key field


def box_label(box) -> str:
    return "points" if box is None else f"{box[0]:g}x{box[1]:g}"


def _joint(tables: Mapping[str, list[np.ndarray]]) -> np.ndarray:
    comp = [np.sqrt(np.mean(np.square(v), axis=0)) for v in tables.values()]
    return np.sqrt(np.mean(np.square(comp), axis=0))


def loo_month(
    mmdd: str,
    src_root: str | Path,
    out_root: str | Path,
    cache_dir: str | Path,
    *,
    first_year: int | None = None,
    last_year: int | None = None,
    boxes: Sequence = BOXES,
    scores: Mapping[str, tuple[float, float] | None] = SCORES,
    box_shifts: int = 1,
    obs_root=None,
    obs_fields=None,
    components: Sequence[str] = ("elm", "mpassi"),
    log=print,
) -> pd.DataFrame:
    """One row per (held-out year, box) with the pick and its truth distance per score.

    ``box_shifts`` (1, 2, 4) selects with distances averaged over that many
    placements of the box edges, as ``clim_restart.compute_distances``.
    """
    plan = cr.make_plan(f"2000-{mmdd}", "obs", src_root=src_root, out_root=out_root,
                        components=components, first_year=first_year, last_year=last_year,
                        exclude_target_year=False, obs_root=obs_root, obs_fields=obs_fields)
    years = plan.candidate_years
    quiet = lambda *a: None
    data = {}
    for comp in components:
        data[comp] = (cr.load_key_fields(plan, comp, cache_dir, log=quiet),
                      cr.load_obs_fields(plan, comp, years, cache_dir, log=quiet),
                      cr.comparison_points(plan, comp, cache_dir),
                      *cr._model_points(plan, comp, cache_dir))
    rows = []
    for k, year in enumerate(years):
        cand = [i for i in range(len(years)) if i != k]
        truth = {}
        for sname, sbox in scores.items():
            for fset in FIELD_SETS:
                t = {}
                for comp, (model, _, _, lat, lon) in data.items():
                    names = list(plan.obs_fields_for(comp)) if fset == "observable" else list(model)
                    t[comp] = []
                    for name in names:
                        values, weight = model[name]
                        vb, wb, _ = cr._apply_box(values[cand + [k]], weight, lat, lon, sbox)
                        t[comp].append(cr.field_distance(vb[:-1], wb, target=vb[-1]))
                truth[sname, fset] = _joint(t)
        for box in boxes:
            t = {}
            for comp, (model, obs, pts, _, _) in data.items():
                t[comp] = []
                for name in plan.obs_fields_for(comp):
                    values, weight = model[name]
                    values = values[cand]
                    if comp == "elm":
                        values, weight = cr._columns_to_gridcells(values, weight, pts["col_grid"],
                                                                  pts["lat"].size)
                    d2 = []
                    for off in cr.box_offsets(box, box_shifts):
                        vb, wb, ob = cr._apply_box(values, weight, pts["lat"], pts["lon"], box,
                                                   obs=obs[name][cand + [k]], offset=off)
                        d2.append(cr.field_distance(vb, wb, obs_years=ob[:-1], obs_target=ob[-1]) ** 2)
                    t[comp].append(np.sqrt(np.mean(d2, axis=0)))
            pick = int(np.argmin(_joint(t)))
            row = {"month": mmdd, "year": year, "box": box_label(box), "picked": years[cand[pick]]}
            for (sname, fset), d in truth.items():
                key = f"{sname},{fset}"
                row[f"D[{key}]"] = d[pick]
                row[f"rank[{key}]"] = int((d < d[pick]).sum() + 1)
                row[f"best[{key}]"] = d.min()
                row[f"rms[{key}]"] = float(np.sqrt(np.mean(d ** 2)))
            rows.append(row)
        log(f"  {mmdd} {year}: done ({k + 1}/{len(years)})")
    return pd.DataFrame(rows)


def _loo_args(args):
    mmdd, kwargs = args
    return loo_month(mmdd, log=lambda *a: None, **kwargs)


def run_loo(months: Sequence[str], *, parallel: bool = True, **kwargs) -> pd.DataFrame:
    """``loo_month`` for several months (one process per month when ``parallel``)."""
    if parallel:
        with ProcessPoolExecutor(max_workers=len(months)) as pool:
            parts = list(pool.map(_loo_args, [(m, kwargs) for m in months]))
    else:
        parts = [loo_month(m, **kwargs) for m in months]
    return add_skill(pd.concat(parts, ignore_index=True))


def score_keys(df: pd.DataFrame) -> list[str]:
    return [c[2:-1] for c in df.columns if c.startswith("D[")]


def add_skill(df: pd.DataFrame) -> pd.DataFrame:
    for key in score_keys(df):
        df[f"skill[{key}]"] = (df[f"rms[{key}]"] - df[f"D[{key}]"]) / (df[f"rms[{key}]"] - df[f"best[{key}]"])
    return df


def skill_table(df: pd.DataFrame, key: str, boxes: Sequence[str]) -> pd.DataFrame:
    """Mean skill by box (rows) and month (columns, plus ``all``)."""
    t = df.pivot_table(index="box", columns="month", values=f"skill[{key}]", aggfunc="mean").loc[boxes]
    t["all"] = df.pivot_table(index="box", values=f"skill[{key}]", aggfunc="mean").loc[boxes].iloc[:, 0]
    return t


def paired_difference(df: pd.DataFrame, key: str, boxes: Sequence[str], reference: str) -> pd.DataFrame:
    """Mean and two standard errors of the skill difference from ``reference`` (paired by year)."""
    p = df.pivot_table(index=["month", "year"], columns="box", values=f"skill[{key}]")
    out = {}
    for month in [*sorted(df["month"].unique()), "all"]:
        sub = p if month == "all" else p.loc[month]
        d = sub[list(boxes)].sub(sub[reference], axis=0)
        out[month] = [f"{d[b].mean():+.2f} ± {2 * d[b].std() / np.sqrt(len(d)):.2f}" for b in boxes]
    return pd.DataFrame(out, index=list(boxes))


def choose_box(
    df: pd.DataFrame,
    boxes: Sequence[str],
    primary: Sequence[str] = ("points,observable", "points,all"),
    secondary: Sequence[str] = ("30x60,observable", "30x60,all"),
) -> tuple[str, pd.DataFrame]:
    """Box choice: tied with the best on the native grid, then best at large scale.

    Primary criterion (``primary`` scores, closeness at model points): a box is
    tied if, for every primary score, its skill is within two standard errors
    (paired by held-out year, all months) of the best box. Tie-breaker: the
    highest mean skill over the ``secondary`` scores (large-scale closeness).

    Returns ``(chosen box, table)``; the table has, per box, the skill and the
    paired difference from the best box for each primary score, ``tied``, and
    the mean secondary skill.
    """
    table = pd.DataFrame(index=list(boxes))
    tied = pd.Series(True, index=list(boxes))
    for key in primary:
        p = df.pivot_table(index=["month", "year"], columns="box", values=f"skill[{key}]")[list(boxes)]
        best = p.mean().idxmax()
        diff = p.sub(p[best], axis=0)
        se2 = 2 * diff.std() / np.sqrt(len(diff))
        table[f"skill [{key}]"] = p.mean()
        table[f"vs best [{key}]"] = [f"{diff[b].mean():+.3f} ± {se2[b]:.3f}" for b in boxes]
        tied &= (diff.mean() + se2 >= 0) | (diff.columns == best)
    table["tied with best on grid"] = tied
    table["large-scale skill (tie-breaker)"] = pd.concat(
        [df.pivot_table(index="box", values=f"skill[{k}]").iloc[:, 0] for k in secondary], axis=1
    ).mean(axis=1).loc[list(boxes)]
    chosen = table.loc[tied, "large-scale skill (tie-breaker)"].idxmax()
    return chosen, table
