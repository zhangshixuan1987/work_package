"""Atmospheric and land initialization benefits of the F20TR hindcasts (``hind_lead_eval/05``).

The benefit of an initialization change is a difference of skill between two experiments that
differ only in that change, e.g. DA_ATM - HYB_BASE (DART atmosphere instead of ERA5 + DART
perturbations, same DATM-forced land files member by member) or DA_COUPLED - DA_ATM (DART land
instead of the DATM-forced land, same DART atmosphere files). Experiments are referred to by
roles (``E00``, ``E10``, ``E01``, ``E11``: atmosphere, land initialized by DART when 1), and a
contrast is ``(roles added, roles subtracted)``.

Fields are the daily values of every member on a region's bounding box, (init, member, lead, lat,
lon), cached once; the skill is computed on weekly means (weeks 1-8) and, for the figure in the
style of Richter et al. (2024), on running means over ``smooth_days`` lead days at every lead day.
Anomalies are taken from the observed or the AMIP climatology (``util.fcst_init_data``).
The skill is the pattern ACC over the region (centred, cos-lat weighted) of the ensemble mean for
each start and week, averaged over the starts; the gridpoint ACC across the starts is also given
but rests on very few starts.

Uncertainty: a paired bootstrap that resamples the start dates (within each season) and the
member indices jointly, with the same draw for every experiment, so a contrast's interval reflects
both the few starts and the finite ensembles.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from configs.fcst_init_experiment_config import MEMBERS, init_dates
from util.fcst_init_data import model_climatology, open_campaign, open_reference, reference_climatology
from util.fcst_init_skill import WEEKS, weekly_means

ROLE_LABELS = {"E00": "reference atm + reference land", "E10": "DART atm + reference land",
               "E01": "reference atm + DART land", "E11": "DART atm + DART land"}


# ---------------------------------------------------------------------------
# Experiments and members
# ---------------------------------------------------------------------------
def role_configs(roles: dict) -> dict[str, str]:
    """``roles`` without the unset ones (value None)."""
    return {r: c for r, c in roles.items() if c}


def member_set(roles: dict, which: str = "matched") -> tuple[str, ...]:
    """Member indices used for every experiment.

    ``matched``: the indices every role's configuration has (for HYB_BASE, DA_ATM and DA_COUPLED the
    11 indices that HYB_BASE EN01-EN20 share with the DART member list), so member k of every
    experiment starts from the same land (or atmosphere) files and all ensembles have one size.
    ``all``: each configuration's own 20 members (only allowed when every role has the same list).
    """
    lists = [MEMBERS[c] for c in role_configs(roles).values()]
    if which == "matched":
        common = set(lists[0]).intersection(*lists[1:])
        return tuple(m for m in lists[0] if m in common)
    if which == "all":
        if any(tuple(x) != tuple(lists[0]) for x in lists):
            raise ValueError("MEMBER_SET='all' needs identical member lists; use 'matched'")
        return tuple(lists[0])
    raise ValueError("member set must be 'matched' or 'all'")


# ---------------------------------------------------------------------------
# Daily fields (cached) and their weekly and running means
# ---------------------------------------------------------------------------
def _box(da: xr.DataArray, region: dict) -> xr.DataArray:
    lo0, lo1 = region["lon"]
    da = da.sel(lat=slice(*region["lat"]))
    return da.sel(lon=slice(lo0, lo1)) if lo0 <= lo1 else da.where((da.lon >= lo0) | (da.lon <= lo1), drop=True)


def daily_members(config: str, season: str, variable: str, members, region: dict, *, max_lead: int,
                  data_root, land_support_dir, cache: Path | None = None, force: bool = False) -> xr.DataArray:
    """(init, member, lead, lat, lon) daily values of ``members`` on the region's box; cached in ``cache``."""
    key = ",".join(members)
    if cache is not None and cache.exists() and not force:
        da = xr.load_dataarray(cache)
        if da.attrs.get("members") == key and da.attrs.get("region") == str(region) and da.sizes["lead"] == max_lead:
            return da
    field = open_campaign(config, season, variable, members=list(members), max_lead=max_lead,
                          data_root=data_root, land_support_dir=land_support_dir)
    out = _box(field, region).drop_vars("time").load().transpose("init", "member", "lead", "lat", "lon")
    out = out.astype("float32").rename(variable).assign_attrs(config=config, season=season, members=key,
                                                              region=str(region), units=field.attrs.get("units", ""))
    if cache is not None:
        out.to_netcdf(cache)
    return out


def daily_references(variable: str, season: str, years, region: dict, *, max_lead: int, observation_root,
                     amip_root, amip_cases, land_support_dir, cache: Path | None = None,
                     force: bool = False) -> xr.Dataset:
    """Daily observations and the observed and AMIP climatologies, (init, lead, lat, lon), on the box."""
    tag = f"{min(years)}-{max(years)}|{amip_cases}|{region}|{max_lead}"
    if cache is not None and cache.exists() and not force:
        ds = xr.load_dataset(cache)
        if ds.attrs.get("tag") == tag:
            return ds
    inits = init_dates([season])[season]
    obs = open_reference(variable, inits, max_lead=max_lead, observation_root=observation_root)
    clim = reference_climatology(variable, inits, years, max_lead=max_lead, observation_root=observation_root)
    mclim = model_climatology(variable, inits, years, max_lead=max_lead, amip_root=amip_root, cases=amip_cases,
                              land_support_dir=land_support_dir)
    out = xr.Dataset({n: _box(d, region).drop_vars("time").load().transpose("init", "lead", "lat", "lon").astype("float32")
                      for n, d in (("obs", obs), ("clim_obs", clim), ("clim_amip", mclim))})
    out.attrs.update(tag=tag, clim_obs=clim.attrs.get("climatology_source", ""),
                     clim_amip=mclim.attrs.get("climatology_source", ""))
    if cache is not None:
        out.to_netcdf(cache)
    return out


def to_weeks(da: xr.DataArray) -> xr.DataArray:
    """Weekly means (weeks 1-8 = lead days 1-7, ..., 50-56), dimension ``week``."""
    return weekly_means(da)


def running_leads(da: xr.DataArray, days: int) -> xr.DataArray:
    """Centred running mean over ``days`` lead days (shorter windows at the ends), dimension ``lead``."""
    if days <= 1:
        return da
    return da.rolling(lead=days, center=True, min_periods=days // 2 + 1).mean()


def lead_in_weeks(lead) -> np.ndarray:
    """Lead day (1 = the start date, a 00-24Z mean) as elapsed weeks at the middle of the day."""
    return (np.asarray(lead, dtype=float) - 0.5) / 7.0


# ---------------------------------------------------------------------------
# Skill
# ---------------------------------------------------------------------------
def region_weights(mask: xr.DataArray) -> np.ndarray:
    """(lat, lon) cos-lat weights inside the region mask, zero outside."""
    w = np.cos(np.deg2rad(mask.lat.values))[:, None] * mask.values.astype(float)
    return w / w.sum()


def pattern_acc_np(f: np.ndarray, o: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Centred weighted pattern correlation over the last two axes (points missing in either are dropped)."""
    valid = np.isfinite(f) & np.isfinite(o) & (w > 0)
    ww = np.where(valid, w, 0.0)
    s = ww.sum((-2, -1))
    with np.errstate(invalid="ignore", divide="ignore"):
        mf = (np.where(valid, f, 0) * ww).sum((-2, -1)) / s
        mo = (np.where(valid, o, 0) * ww).sum((-2, -1)) / s
        df = np.where(valid, f - mf[..., None, None], 0)
        do = np.where(valid, o - mo[..., None, None], 0)
        cov = (ww * df * do).sum((-2, -1))
        return cov / np.sqrt((ww * df ** 2).sum((-2, -1)) * (ww * do ** 2).sum((-2, -1)))


def gridpoint_acc_np(f: np.ndarray, o: np.ndarray) -> np.ndarray:
    """Pearson correlation across the first axis (starts) at every point: (week, lat, lon)."""
    fa, oa = f - np.nanmean(f, 0), o - np.nanmean(o, 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.nanmean(fa * oa, 0) / np.sqrt(np.nanmean(fa ** 2, 0) * np.nanmean(oa ** 2, 0))


def contrast_values(skill: dict[str, np.ndarray], contrasts: dict) -> dict[str, np.ndarray]:
    """{name: sum(skill of added roles) - sum(skill of subtracted roles)} for the contrasts whose roles exist."""
    out = {}
    for name, (plus, minus) in contrasts.items():
        if all(r in skill for r in (*plus, *minus)):
            out[name] = sum(skill[r] for r in plus) - sum(skill[r] for r in minus)
    return out


def _init_strata(inits_season: list[str]) -> list[np.ndarray]:
    """Indices of the starts of each season (resampled within the season)."""
    seasons = np.asarray(inits_season)
    return [np.flatnonzero(seasons == s) for s in dict.fromkeys(seasons)]


def era5_benefits(anoms: dict[str, np.ndarray], obs: np.ndarray, w: np.ndarray, contrasts: dict,
                  inits_season: list[str], *, n_boot: int = 1000, seed: int = 42) -> dict:
    """Panel (b): skill of each role and the contrasts against the observations, with a paired bootstrap.

    ``anoms[role]``: (init, member, week, lat, lon) forecast anomalies, members in matched order;
    ``obs``: (init, week, lat, lon) observed anomalies. Returns per-start pattern ACC, the mean
    curves, bootstrap percentiles (2.5, 97.5) and the gridpoint-ACC maps.
    """
    rng = np.random.default_rng(seed)
    nmem = next(iter(anoms.values())).shape[1]
    strata = _init_strata(inits_season)
    per_start = {r: pattern_acc_np(a.mean(1), obs, w) for r, a in anoms.items()}          # (init, week)
    curves = {r: v.mean(0) for r, v in per_start.items()}
    boot = {name: [] for name in contrast_values(curves, contrasts)}
    boot_skill = {r: [] for r in anoms}
    for _ in range(n_boot):
        ii = np.concatenate([rng.choice(s, s.size) for s in strata])
        mm = rng.integers(0, nmem, nmem)
        sk = {r: pattern_acc_np(a[ii][:, mm].mean(1), obs[ii], w).mean(0) for r, a in anoms.items()}
        for r, v in sk.items():
            boot_skill[r].append(v)
        for name, v in contrast_values(sk, contrasts).items():
            boot[name].append(v)
    maps = {r: gridpoint_acc_np(a.mean(1), obs) for r, a in anoms.items()}
    return dict(per_start=per_start, skill=curves, contrasts=contrast_values(curves, contrasts),
                skill_ci={r: np.nanpercentile(v, (2.5, 97.5), 0) for r, v in boot_skill.items()},
                contrast_ci={n: np.nanpercentile(v, (2.5, 97.5), 0) for n, v in boot.items()},
                maps=maps, map_contrasts=contrast_values(maps, contrasts))


def model_member_benefits(anoms: dict[str, np.ndarray], truth_role: str, w: np.ndarray, contrasts: dict,
                          inits_season: list[str], *, n_boot: int = 1000, seed: int = 42) -> dict:
    """Panel (c): verification against a held-out member of ``truth_role``, rotated over all members.

    For held-out member k, every experiment's forecast is the mean of its other members (index k is
    dropped everywhere: in DA_ATM and DA_COUPLED member k shares the atmosphere initial state, so
    keeping it would favour one experiment), and the truth is member k of ``truth_role``. Anomalies
    of forecasts and truth come from the same (AMIP) climatology. Skill is averaged over k.
    """
    rng = np.random.default_rng(seed)
    nmem = anoms[truth_role].shape[1]
    strata = _init_strata(inits_season)
    # per_start[r][k]: (init, week) pattern ACC with member k held out
    per_start = {r: np.stack([pattern_acc_np(np.delete(a, k, 1).mean(1), anoms[truth_role][:, k], w)
                              for k in range(nmem)]) for r, a in anoms.items()}
    curves = {r: v.mean((0, 1)) for r, v in per_start.items()}
    boot = {name: [] for name in contrast_values(curves, contrasts)}
    boot_skill = {r: [] for r in anoms}
    for _ in range(n_boot):
        ii = np.concatenate([rng.choice(s, s.size) for s in strata])
        kk = rng.integers(0, nmem, nmem)
        sk = {r: v[kk][:, ii].mean((0, 1)) for r, v in per_start.items()}
        for r, v in sk.items():
            boot_skill[r].append(v)
        for name, v in contrast_values(sk, contrasts).items():
            boot[name].append(v)
    return dict(per_start=per_start, skill=curves, contrasts=contrast_values(curves, contrasts),
                skill_ci={r: np.nanpercentile(v, (2.5, 97.5), 0) for r, v in boot_skill.items()},
                contrast_ci={n: np.nanpercentile(v, (2.5, 97.5), 0) for n, v in boot.items()})


def results_table(result: dict, coords, panel: str, label: str, axis: str = "week") -> pd.DataFrame:
    """Long table of the skill and contrast curves with their bootstrap intervals."""
    rows = []
    for kind, curves, cis in (("skill", result["skill"], result["skill_ci"]),
                              ("contrast", result["contrasts"], result["contrast_ci"])):
        for name, v in curves.items():
            for j, x in enumerate(coords):
                rows.append(dict(panel=panel, set=label, kind=kind, name=name, axis=axis, step=int(x), value=float(v[j]),
                                 ci_low=float(cis[name][0][j]), ci_high=float(cis[name][1][j])))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
COLORS = {"atmosphere": "#c9b93d", "land": "#348a43", "interaction": "0.4", "total": "k"}


def plot_benefit_figure(weeks, panel_b: dict, panel_c: dict | None, *, title: str, ylabel: str = "ΔACC",
                        names=("atmosphere", "land")):
    """Three panels: (a) schematic hypothesis, (b) ERA5 verification, (c) model-member verification."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2), constrained_layout=True)
    x = np.linspace(1, max(weeks), 200)
    ax = axes[0]
    ax.plot(x, np.exp(-(x - 1) / 1.6), color=COLORS["atmosphere"], lw=2.4, label="Atmosphere")
    ax.plot(x, 0.45 * (1 - np.exp(-(x - 0.5) / 2.0)) * np.exp(-(x - 1) / 9.0), color=COLORS["land"], lw=2.4,
            label="Land")
    ax.set_title("(a) Hypothesis (schematic, not data)", loc="left")
    ax.set_ylabel("Relative initialization benefit")
    ax.set_yticks([])
    lims = []
    for ax, res, ttl in ((axes[1], panel_b, "(b) Verification: ERA5"),
                         (axes[2], panel_c, "(c) Verification: held-out model member")):
        if res is None:
            ax.text(0.5, 0.5, "not computed", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(ttl, loc="left")
            continue
        for name in names:
            if name not in res["contrasts"]:
                continue
            v, (lo, hi) = res["contrasts"][name], res["contrast_ci"][name]
            ax.plot(weeks, v, color=COLORS.get(name, "k"), lw=2.2, marker="o", label=name.capitalize())
            ax.fill_between(weeks, lo, hi, color=COLORS.get(name, "k"), alpha=0.18, lw=0)
            lims += [np.nanmin(lo), np.nanmax(hi)]
        ax.axhline(0, color="k", lw=0.8)
        ax.set_title(ttl, loc="left")
        ax.set_ylabel(ylabel)
    if lims:
        pad = 0.05 * (max(lims) - min(lims) or 1)
        for ax in axes[1:]:
            ax.set_ylim(min(lims) - pad, max(lims) + pad)
    for ax in axes:
        ax.set_xlabel("Forecast lead week")
        ax.set_xlim(min(weeks) - 0.2, max(weeks) + 0.2)
        ax.set_xticks(list(weeks))
        ax.grid(alpha=0.2)
        ax.legend(frameon=False, fontsize=9)
    fig.suptitle(title)
    return fig


def plot_skill_curves(weeks, panel_b: dict, panel_c: dict | None, roles: dict, styles: dict, *, title: str):
    """Absolute pattern ACC of every role, ERA5 and model-member verification."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4), constrained_layout=True, sharey=True)
    for ax, res, ttl in ((axes[0], panel_b, "(a) Verification: ERA5"), (axes[1], panel_c, "(b) Held-out model member")):
        if res is not None:
            for r, v in res["skill"].items():
                cfg = roles[r]
                st = styles.get(cfg, {})
                lo, hi = res["skill_ci"][r]
                ax.plot(weeks, v, marker="o", label=f"{r} {cfg}", color=st.get("color"), lw=2)
                ax.fill_between(weeks, lo, hi, color=st.get("color"), alpha=0.12, lw=0)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_title(ttl, loc="left")
        ax.set_xlabel("Forecast lead week")
        ax.set_xticks(list(weeks))
        ax.grid(alpha=0.2)
        ax.legend(frameon=False, fontsize=8)
    axes[0].set_ylabel("Pattern ACC (ensemble mean)")
    fig.suptitle(title)
    return fig


def plot_paper_figure(x_weeks, panel_b: dict, panel_c: dict | None, *, names=("atmosphere", "land"),
                      title: str = "", ylabel: str = "ΔACC", window=(2, 4), model_label: str = "E3SM"):
    """The three-panel figure in the style of Richter et al. (2024, Fig. 1): curves at every lead day,
    the area down to zero shaded, the weeks-2-4 window in grey, the legend in panel (b) only."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.0), constrained_layout=True)
    xmax = float(np.nanmax(x_weeks))
    x = np.linspace(0, xmax, 300)
    # schematic shapes of Richter et al. (2024, Fig. 1a): atmospheric benefit high at week 0 and gone after
    # about 2 weeks; land benefit rising to a peak near 1.5 weeks and decaying slowly
    sketch = {"atmosphere": 0.56 * np.exp(-(x / 1.25) ** 1.3),
              "land": 0.226 * (1 - np.exp(-x / 0.5)) * np.exp(-np.maximum(x - 1.5, 0) / 5.0)}
    ax = axes[0]
    for name in names:
        if name in sketch:
            c = COLORS.get(name, "k")
            ax.fill_between(x, 0, sketch[name], color=c, alpha=0.15, lw=0)
            ax.plot(x, sketch[name], color=c, lw=2)
    ax.set_title("(a) Hypothesis", fontweight="bold")
    ax.set_ylabel("Predictability", fontweight="bold")
    ax.set_yticks([])
    lims = [0.0]
    for ax, res, ttl in ((axes[1], panel_b, f"(b) {model_label} (Verif = ERA5)"),
                         (axes[2], panel_c, f"(c) {model_label} (Verif = held-out member)")):
        ax.set_title(ttl, fontweight="bold")
        ax.set_ylabel(ylabel, fontweight="bold")
        if res is None:
            ax.text(0.5, 0.5, "not computed", ha="center", va="center", transform=ax.transAxes)
            continue
        for name in names:
            if name not in res["contrasts"]:
                continue
            v, c = np.asarray(res["contrasts"][name]), COLORS.get(name, "k")
            ax.fill_between(x_weeks, 0, v, color=c, alpha=0.15, lw=0)
            ax.plot(x_weeks, v, color=c, lw=2, label=name)
            lims += [np.nanmin(v), np.nanmax(v)]
        ax.axhline(0, color="k", lw=0.6)
    hi, lo = max(lims), min(lims)
    pad = 0.06 * ((hi - lo) or 1)
    for ax in axes:
        ax.axvspan(*window, color="0.5", alpha=0.12, lw=0)
        ax.set_xlim(0, xmax)
        ax.set_xlabel("Week", fontweight="bold")
        ax.grid(alpha=0.25)
    for ax in axes[1:]:
        ax.set_ylim(lo - pad, hi + pad)
    axes[0].set_ylim(0, None)
    axes[1].legend(frameon=True, fontsize=9, loc="upper right")
    if title:
        fig.suptitle(title)
    return fig


def week_index(weeks, wanted) -> list[int]:
    return [list(weeks).index(w) for w in wanted if w in list(weeks)]


__all__ = ["ROLE_LABELS", "COLORS", "WEEKS", "role_configs", "member_set", "daily_members", "daily_references",
           "to_weeks", "running_leads", "lead_in_weeks", "plot_paper_figure",
           "region_weights", "pattern_acc_np", "gridpoint_acc_np", "contrast_values", "era5_benefits",
           "model_member_benefits", "results_table", "plot_benefit_figure", "plot_skill_curves", "week_index"]
