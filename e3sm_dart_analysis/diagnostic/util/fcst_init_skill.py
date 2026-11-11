"""Lead-day and weekly S2S skill of the F20TR initialization hindcasts.

Inputs are the stacked fields of ``util.fcst_init_data``: model (init, member, lead,
lat, lon), reference and climatology (init, lead, lat, lon), on the 1-degree grid.

Daily ensemble performance (``lead_scores``), pooled over the start dates of a season:
ensemble-mean RMSE, bias, ensemble spread (unbiased member variance, corrected by
sqrt((N+1)/N)), spread/RMSE ratio, and fair CRPS, cos-lat weighted over each region
on the points where every source is valid.

Weekly S2S verification follows the protocol of the ESP-Lab ``s2s_skill`` notebooks
(``esp_lab/diagnostics/s2s_core.py``; the window definition, the gridpoint ACC, and its
one-sided t-test are ported from there):

* weeks 1-8 are non-overlapping 7-day blocks of lead days 1-7, ..., 50-56;
* observed anomalies are taken from the reference climatology on the verification days,
  forecast anomalies from the climatology of the free-running E3SMv3 AMIP runs
  (``bias="amip"``; ``util.fcst_init_data.model_climatology``). A season has only 4
  start dates, too few for the hindcast's own lead-dependent climatology; the AMIP
  climatology removes the model's mean bias without using observations of the verified
  period. ``bias="loo"`` (the ensemble-mean error of each week averaged over the other
  starts of the season) is kept for comparison only: the other starts verify on
  overlapping days, so it leaks their observed anomalies into the forecast anomaly;
* ``pattern_acc``: centred, cos-lat weighted spatial anomaly correlation over a region
  for each start and week (the robust score with few starts), and its mean;
* ``gridpoint_acc``: temporal ACC at each gridpoint over the start dates (pooled over
  seasons when several are given), with the critical correlation of the t-test for
  that sample size (8 starts -> r_crit = 0.62 at 95%);
* ``msss``: 1 - MSE(forecast anomaly) / MSE(climatology forecast), regional;
* ``regional_index`` / ``index_skill``: ACC and nRMSE over the starts of the cos-lat weighted
  regional-mean anomaly (the regional index of ESP-Lab ``3a_sst_skill_ts``), with the one-sided
  t-test; the area means of the gridpoint skill maps (``window_skill`` + ``window_region_summary``)
  are the regional curves of ESP-Lab ``2a_regional_acc_skill_ts``;
* ``window_skill``: gridpoint ACC and normalized RMSE of multi-week means (weeks 1-2, 3-4,
  5-6, 7-8) over the start dates, the S2S analogue of the ESP-Lab lead-time skill maps
  (``1a_*_leadtime_acc_skill_map``, ``1b_*_leadtime_rmse_skill_map``).
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.stats as stats
import xarray as xr

# ---------------------------------------------------------------------------
# Weekly windows (ESP-Lab s2s_core.build_weekly_window)
# ---------------------------------------------------------------------------
WEEKS = {w: (7 * (w - 1) + 1, 7 * w) for w in range(1, 9)}
WEEK_LABELS = {w: f"Week {w} (days {a}-{b})" for w, (a, b) in WEEKS.items()}


def weekly_means(field: xr.DataArray, weeks=tuple(WEEKS), lead_dim: str = "lead") -> xr.DataArray:
    """Mean over the lead days of each week (dim ``week``); a week needs all its days present."""
    leads = set(field[lead_dim].values.tolist())
    out, kept = [], []
    for w in weeks:
        a, b = WEEKS[w]
        if not set(range(a, b + 1)) <= leads:
            continue
        mean = field.sel({lead_dim: list(range(a, b + 1))}).mean(lead_dim, skipna=False)
        out.append(mean.drop_vars([c for c in ("time",) if c in mean.coords]))
        kept.append(w)
    if not out:
        raise ValueError("No complete week in the lead range")
    return xr.concat(out, pd.Index(kept, name="week"))


# ---------------------------------------------------------------------------
# Daily ensemble performance
# ---------------------------------------------------------------------------
def _weights(lat, valid, masks):
    cos = np.cos(np.deg2rad(lat)).clip(min=0)[:, None]
    return [np.where(valid & masks.sel(region=r).values, cos, 0.0) for r in masks.region.values]


def _area(field, w):
    """cos-lat weighted mean over the last two axes with weights ``w`` (broadcast)."""
    num = np.nansum(field * w, axis=(-2, -1))
    den = np.sum(np.broadcast_to(w, field.shape), axis=(-2, -1))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num / den, np.nan)


def lead_scores(model: xr.DataArray, reference: xr.DataArray, masks: xr.DataArray) -> xr.Dataset:
    """(region, lead) RMSE, bias, spread, spread ratio and fair CRPS pooled over ``init``,
    plus the per-start RMSE (region, init, lead)."""
    from util.forecast_daily_skill import _fair_crps

    model = model.transpose("init", "member", "lead", "lat", "lon")
    reference = reference.transpose("init", "lead", "lat", "lon")
    f = np.asarray(model.values, dtype=np.float64)
    o = np.asarray(reference.values, dtype=np.float64)
    n = f.shape[1]
    valid = np.isfinite(o) & np.isfinite(f).all(axis=1)              # (init, lead, lat, lon)
    w = [np.where(valid, wr, 0.0) for wr in _weights(model.lat.values, True, masks)]
    mean = f.mean(axis=1)
    err = np.where(valid, mean - o, np.nan)
    var = np.where(valid, f.var(axis=1, ddof=1), np.nan)
    crps = np.stack([np.where(valid[i], _fair_crps(f[i], o[i]), np.nan) for i in range(f.shape[0])])
    out = {}
    for name, fld in (("mse", err ** 2), ("bias", err), ("var", var), ("crps_fair", crps)):
        out[name] = np.stack([_area(fld, wr) for wr in w])        # (region, init, lead)
    regions = masks.region.values
    coords = {"region": regions, "init": model.init.values, "lead": model.lead.values}
    ds = xr.Dataset({
        "rmse": (("region", "lead"), np.sqrt(np.nanmean(out["mse"], axis=1))),
        "bias": (("region", "lead"), np.nanmean(out["bias"], axis=1)),
        "spread": (("region", "lead"), np.sqrt(np.nanmean(out["var"], axis=1) * (n + 1) / n)),
        "crps_fair": (("region", "lead"), np.nanmean(out["crps_fair"], axis=1)),
        "rmse_init": (("region", "init", "lead"), np.sqrt(out["mse"])),
    }, coords=coords)
    ds["spread_ratio"] = ds.spread / ds.rmse
    ds.attrs.update(n_members=n, n_inits=f.shape[0])
    return ds


# ---------------------------------------------------------------------------
# Weekly anomalies and S2S skill
# ---------------------------------------------------------------------------
def weekly_anomalies(model: xr.DataArray, reference: xr.DataArray, climatology: xr.DataArray, *,
                     bias: str = "amip", model_climatology: xr.DataArray | None = None
                     ) -> tuple[xr.DataArray, xr.DataArray]:
    """Weekly ensemble-mean and reference anomalies (init, week, lat, lon).

    The reference anomaly is always ``reference - climatology``. The forecast anomaly by ``bias``:

    * ``"amip"``: ``model - model_climatology`` (the AMIP climatology of
      ``util.fcst_init_data.model_climatology``). Removes the model's mean bias without
      using any observation of the verified period; the initialization drift stays in
      the forecast error.
    * ``"none"``: ``model - climatology`` (the reference climatology; the model bias counts as error).
    * ``"loo"``: as ``"none"`` minus the mean ensemble-mean error of the same week over the
      other start dates (needs >= 2 starts). With starts 5 days apart in one season the
      other starts verify on overlapping days, so this leaks their observed anomalies into
      the forecast anomaly; kept for comparison.
    """
    fw = weekly_means(model.mean("member"))
    ow = weekly_means(reference)
    cw = weekly_means(climatology)
    fa, oa = fw - cw, ow - cw
    if bias == "amip":
        if model_climatology is None:
            raise ValueError("bias='amip' needs model_climatology")
        fa = fw - weekly_means(model_climatology)
    elif bias == "loo":
        n = fa.sizes["init"]
        if n < 2:
            raise ValueError("bias='loo' needs at least 2 start dates")
        err = fw - ow
        ok = err.notnull()
        others = ok.sum("init") - ok                       # valid other starts at each point
        loo = (err.sum("init") - err.fillna(0)) / others   # mean error of the other starts
        fa = fa - loo.where(others > 0)
    elif bias != "none":
        raise ValueError("bias must be 'amip', 'none' or 'loo'")
    fa.attrs, oa.attrs = {"bias": bias}, {}
    return fa, oa


def pattern_acc(fa: xr.DataArray, oa: xr.DataArray, masks: xr.DataArray) -> xr.DataArray:
    """(region, init, week) centred cos-lat weighted pattern ACC."""
    f = fa.transpose("init", "week", "lat", "lon").values.astype(np.float64)
    o = oa.transpose("init", "week", "lat", "lon").values.astype(np.float64)
    valid = np.isfinite(f) & np.isfinite(o)
    out = []
    for wr in _weights(fa.lat.values, True, masks):
        w = np.where(valid, wr, 0.0)
        mf, mo = _area(np.where(valid, f, 0), w), _area(np.where(valid, o, 0), w)
        df, do = np.where(valid, f - mf[..., None, None], 0), np.where(valid, o - mo[..., None, None], 0)
        cov, vf, vo = _area(df * do, w), _area(df ** 2, w), _area(do ** 2, w)
        with np.errstate(invalid="ignore", divide="ignore"):
            out.append(cov / np.sqrt(vf * vo))
    return xr.DataArray(np.stack(out), dims=("region", "init", "week"),
                        coords={"region": masks.region.values, "init": fa.init.values, "week": fa.week.values},
                        name="pattern_acc")


def msss(fa: xr.DataArray, oa: xr.DataArray, masks: xr.DataArray) -> xr.Dataset:
    """(region, week) RMSE of the anomaly forecast, of the climatology forecast, and MSSS."""
    f = fa.transpose("init", "week", "lat", "lon").values.astype(np.float64)
    o = oa.transpose("init", "week", "lat", "lon").values.astype(np.float64)
    valid = np.isfinite(f) & np.isfinite(o)
    mse_f, mse_c = [], []
    for wr in _weights(fa.lat.values, True, masks):
        w = np.where(valid, wr, 0.0)
        with warnings.catch_warnings():  # weeks without any valid start stay NaN
            warnings.simplefilter("ignore", RuntimeWarning)
            mse_f.append(np.nanmean(_area(np.where(valid, (f - o) ** 2, 0), w), axis=0))
            mse_c.append(np.nanmean(_area(np.where(valid, o ** 2, 0), w), axis=0))
    coords = {"region": masks.region.values, "week": fa.week.values}
    mf, mc = np.array(mse_f), np.array(mse_c)
    return xr.Dataset({"rmse_anom": (("region", "week"), np.sqrt(mf)),
                       "rmse_clim": (("region", "week"), np.sqrt(mc)),
                       "msss": (("region", "week"), 1 - mf / mc)}, coords=coords)


def gridpoint_acc(fa: xr.DataArray, oa: xr.DataArray, sample_dim: str = "init") -> xr.DataArray:
    """(week, lat, lon) Pearson correlation over ``sample_dim`` (ESP-Lab compute_weekly_acc)."""
    mp, op = fa - fa.mean(sample_dim), oa - oa.mean(sample_dim)
    cov = (mp * op).mean(sample_dim)
    den = np.sqrt((mp ** 2).mean(sample_dim) * (op ** 2).mean(sample_dim))
    acc = xr.where(den > 1e-12, cov / den, np.nan).rename("acc")
    acc.attrs.update(n_samples=int(fa.sizes[sample_dim]))
    return acc


WINDOWS = ((1, 2), (3, 4), (5, 6), (7, 8))   # two-week lead windows of the lead-time skill maps


def window_label(weeks) -> str:
    """``(1, 2)`` -> "Weeks 1-2" (days 1-14)."""
    weeks = tuple(weeks)
    return f"Week {weeks[0]}" if len(weeks) == 1 else f"Weeks {weeks[0]}-{weeks[-1]}"


def window_means(field: xr.DataArray, windows=WINDOWS) -> xr.DataArray:
    """(window, ...) mean of a weekly field (dim ``week``) over the weeks of each window; a window
    needs all its weeks (a missing week gives NaN)."""
    parts = [field.reindex(week=list(w)).mean("week", skipna=False) for w in windows]
    return xr.concat(parts, pd.Index([window_label(w) for w in windows], name="window"))


def window_skill(fa: xr.DataArray, oa: xr.DataArray, windows=WINDOWS, *, sample_dim: str = "init",
                 min_obs_fraction: float = 0.01) -> xr.Dataset:
    """(window, lat, lon) gridpoint skill of the window means over ``sample_dim`` (the start dates).

    * ``acc``: Pearson correlation of the ensemble-mean and observed anomalies (``gridpoint_acc``);
    * ``rmse``: RMSE of the ensemble-mean anomaly; ``rmse_clim``: RMSE of the climatology forecast
      (the RMS of the observed anomaly);
    * ``nrmse`` = ``rmse / rmse_clim``: below 1 the forecast beats climatology (nRMSE^2 = 1 - MSSS at
      the gridpoint). ESP-Lab divides by the standard deviation of its detrended observed anomalies; with
      anomalies from a multi-year climatology and only a few starts the RMS is the matching scale (it keeps
      the observed mean anomaly of the starts). Cells where ``rmse_clim`` is below ``min_obs_fraction`` x
      its median over the map are left out (ESP-Lab ``min_obs_std_fraction``), e.g. dry desert precipitation.
    """
    fw, ow = window_means(fa, windows), window_means(oa, windows)
    acc = gridpoint_acc(fw, ow, sample_dim)
    rmse = np.sqrt(((fw - ow) ** 2).mean(sample_dim))
    clim = np.sqrt((ow ** 2).mean(sample_dim))
    floor = min_obs_fraction * clim.median(("lat", "lon"))
    nrmse = (rmse / clim).where(clim > floor)
    n = int(fa.sizes[sample_dim])
    return xr.Dataset(dict(acc=acc.drop_attrs(), rmse=rmse, rmse_clim=clim, nrmse=nrmse)).assign_attrs(
        n_samples=n, acc_critical_95=critical_acc(n) if n > 2 else np.nan,
        windows=";".join(window_label(w) for w in windows), nrmse_definition="RMSE / RMS of the observed anomaly")


def area_fraction(mask: xr.DataArray, valid: xr.DataArray | None = None) -> float:
    """cos-lat weighted fraction (%) of the valid cells of a (lat, lon) map where ``mask`` is true."""
    valid = mask.notnull() if valid is None else valid
    w = np.cos(np.deg2rad(mask.lat)).broadcast_like(mask).where(valid)
    tot = float(w.sum())
    return 100.0 * float(w.where(mask.fillna(False).astype(bool)).sum()) / tot if tot > 0 else np.nan


def window_region_summary(skill: xr.Dataset, masks: xr.DataArray) -> pd.DataFrame:
    """Per region and window, from the ``window_skill`` maps: cos-lat weighted mean ACC, the area (%) with ACC
    above the 95% critical value, the regional nRMSE sqrt(sum w RMSE^2 / sum w RMSE_clim^2) and the area (%)
    with nRMSE < 1 (beating climatology)."""
    rc = skill.attrs.get("acc_critical_95", np.nan)
    w = (masks.astype(float) * np.cos(np.deg2rad(skill.lat))).transpose("region", "lat", "lon")   # (region, lat, lon)
    dims = ("lat", "lon")

    def wsum(x, valid):
        return xr.dot(w, (x.where(valid, 0.0) if x.dtype.kind == "f" else x.astype(float).where(valid, 0.0)), dim=dims)

    ok_a = skill.acc.notnull()
    ok_r = skill.nrmse.notnull()
    na, nr = xr.dot(w, ok_a.astype(float), dim=dims), xr.dot(w, ok_r.astype(float), dim=dims)
    out = xr.Dataset(dict(
        acc_mean=wsum(skill.acc, ok_a) / na,
        acc_sig_area_pct=100 * wsum(skill.acc > rc, ok_a) / na if np.isfinite(rc) else np.nan * na,
        nrmse_region=np.sqrt(wsum(skill.rmse ** 2, ok_r) / wsum(skill.rmse_clim ** 2, ok_r)),
        nrmse_lt1_area_pct=100 * wsum(skill.nrmse < 1, ok_r) / nr,
    ))
    df = out.to_dataframe().reset_index()[["region", "window", "acc_mean", "acc_sig_area_pct", "nrmse_region",
                                           "nrmse_lt1_area_pct"]]
    return df.astype({"region": str, "window": str})

def regional_index(field: xr.DataArray, masks: xr.DataArray) -> xr.DataArray:
    """(region, ...) cos-lat weighted mean of a (..., lat, lon) field over each region's valid cells."""
    w = (masks.astype(float) * np.cos(np.deg2rad(field.lat))).transpose("region", "lat", "lon")
    valid = field.notnull()
    num = xr.dot(w, field.where(valid, 0.0), dim=("lat", "lon"))
    return num / xr.dot(w, valid.astype(float), dim=("lat", "lon"))


def index_skill(fi: xr.DataArray, oi: xr.DataArray, sample_dim: str = "init") -> xr.Dataset:
    """Skill over ``sample_dim`` of a regional index (``regional_index`` of the ensemble-mean and observed
    anomalies): ``acc`` (Pearson), ``acc_sig`` (above the one-sided 95% critical value), ``nrmse`` = RMSE /
    RMS of the observed index (below 1 beats climatology), ``rmse``, and ``n``."""
    n = int(fi.sizes[sample_dim])
    acc = gridpoint_acc(fi, oi, sample_dim).drop_attrs()
    rmse = np.sqrt(((fi - oi) ** 2).mean(sample_dim))
    clim = np.sqrt((oi ** 2).mean(sample_dim))
    nrmse = xr.where(clim > 0, rmse / clim, np.nan)
    rc = critical_acc(n) if n > 2 else np.nan
    return xr.Dataset(dict(acc=acc, acc_sig=acc > rc, rmse=rmse, rmse_clim=clim, nrmse=nrmse)).assign_attrs(
        n_samples=n, acc_critical_95=rc)


def critical_acc(n_samples: int, alpha: float = 0.05) -> float:
    """One-sided t-test critical correlation for ``n_samples`` (ESP-Lab compute_weekly_acc_significance)."""
    df = max(1, n_samples - 2)
    t = stats.t.ppf(1 - alpha, df)
    return float(t / np.sqrt(df + t ** 2))


def pool_seasons(per_season: dict[str, xr.DataArray]) -> xr.DataArray:
    """Concatenate per-season anomalies along ``init`` (each season's anomalies are already
    relative to its own climatology and drift)."""
    return xr.concat(list(per_season.values()), "init")


# ---------------------------------------------------------------------------
# Land memory and coupling (ESP-Lab 4a_soil_moisture_memory, ensemble-based)
# ---------------------------------------------------------------------------
def member_anomalies(field: xr.DataArray) -> xr.DataArray:
    """(init, member, week, ...) weekly anomalies from the mean over starts and members of each week."""
    w = weekly_means(field)
    return w - w.mean(("init", "member"))


def _corr(a: xr.DataArray, b: xr.DataArray, dims) -> xr.DataArray:
    a, b = a - a.mean(dims), b - b.mean(dims)
    den = np.sqrt((a ** 2).mean(dims) * (b ** 2).mean(dims))
    return xr.where(den > 1e-12, (a * b).mean(dims) / den, np.nan)


def memory_and_coupling(sm: xr.DataArray, other: dict[str, xr.DataArray]) -> xr.Dataset:
    """Per gridpoint and week, over members and starts:

    * ``memory``: r(SM'(week 1), SM'(week w)), persistence of the initial soil-moisture anomaly;
    * ``coupling_<name>``: r(SM'(w), X'(w)) for each field in ``other`` (e.g. TSA, LHFLX).
    """
    smw = member_anomalies(sm)
    dims = ("init", "member")
    out = {"memory": _corr(smw.sel(week=1, drop=True), smw, dims)}
    for name, field in other.items():
        out[f"coupling_{name}"] = _corr(smw, member_anomalies(field), dims)
    ds = xr.Dataset(out)
    ds.attrs.update(n_samples=int(smw.sizes["init"] * smw.sizes["member"]))
    return ds


# ---------------------------------------------------------------------------
# RMM skill pooled over start dates
# ---------------------------------------------------------------------------
def rmm_on_leads(rmm: xr.Dataset, init, max_lead: int) -> xr.Dataset:
    """RMM1/RMM2 ([member,] lead) of one start: lead day 1 = the start date."""
    days = pd.date_range(pd.Timestamp(np.datetime64(init, "D")), periods=max_lead, freq="D")
    out = rmm[["RMM1", "RMM2"]].reindex(time=days.values)
    return out.assign_coords(lead=("time", np.arange(1, max_lead + 1))).swap_dims(time="lead").drop_vars("time")


def rmm_lead_skill(fcst: dict, obs: xr.Dataset, *, max_lead: int = 56) -> xr.Dataset:
    """(lead) skill pooled over start dates (Rashid et al. 2011): ``fcst`` = {init: RMM (member, time)}.

    ``cor_ensmean``/``rmse_ensmean``: bivariate COR/RMSE of the ensemble mean over the starts;
    ``spread``: sqrt of the mean unbiased member variance of the RMM vector;
    ``amp_error_ensmean``: mean ensemble-mean amplitude minus observed amplitude.
    """
    from util.mjo_rmm_projection import bivariate_cor, bivariate_rmse

    f = xr.concat([rmm_on_leads(r, i, max_lead) for i, r in fcst.items()], pd.Index(list(fcst), name="init"))
    o = xr.concat([rmm_on_leads(obs, i, max_lead) for i in fcst], pd.Index(list(fcst), name="init"))
    m1, m2 = f.RMM1.mean("member"), f.RMM2.mean("member")
    n = f.sizes["member"]
    var = (((f.RMM1 - m1) ** 2 + (f.RMM2 - m2) ** 2).sum("member") / max(n - 1, 1)).mean("init")
    cor, rmse = [], []
    for L in f.lead.values:
        e1, e2 = m1.sel(lead=L).values, m2.sel(lead=L).values
        o1, o2 = o.RMM1.sel(lead=L).values, o.RMM2.sel(lead=L).values
        ok = np.isfinite(e1) & np.isfinite(o1)
        cor.append(bivariate_cor(e1[ok], e2[ok], o1[ok], o2[ok]) if ok.sum() >= 3 else np.nan)
        rmse.append(bivariate_rmse(e1[ok], e2[ok], o1[ok], o2[ok]) if ok.any() else np.nan)
    amp = (np.hypot(m1, m2) - np.hypot(o.RMM1, o.RMM2)).mean("init")
    return xr.Dataset({"cor_ensmean": ("lead", np.array(cor)), "rmse_ensmean": ("lead", np.array(rmse)),
                       "spread": np.sqrt(var), "amp_error_ensmean": amp},
                      coords={"lead": f.lead.values}).assign_attrs(n_inits=len(fcst), n_members=n)


# ---------------------------------------------------------------------------
# Scorecard tables
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ScoreSpec:
    metric: str
    higher_is_better: bool


SCORECARD_METRICS = {
    "pattern_acc": ScoreSpec("pattern ACC", True),
    "msss": ScoreSpec("MSSS", True),
    "rmse_anom": ScoreSpec("weekly RMSE", False),
}


def relative_change(table: pd.DataFrame, baseline: str, metric: str) -> pd.DataFrame:
    """Change of ``metric`` relative to ``baseline`` (index config): % for RMSE, difference otherwise."""
    base = table.loc[baseline]
    if metric.startswith("rmse"):
        return 100 * (table - base) / base
    return table - base
