"""Shared inputs of the 9a-9h stratosphere notebooks (polar vortex and SSWs in the 2026 starts).

Every notebook builds a ``StratConfig`` from its run controls (data roots are notebook inputs) and
loads only what it needs:

* ``load_era5``: ERA5 daily zonal-mean indices, their smoothed daily climatologies and the
  calendar-month climatology on all levels (lazy reads: only reductions are loaded);
* ``load_references``: Reanalysis monthly zonal means of every start month (May/Nov 1981-2010;
  Feb/Aug the starts of every fifth year, 1985-2010) for the lead-time climatologies;
* ``load_model``: the pooled 2026 starts, monthly and daily zonal means, and the daily indices with
  their climatology (``cfg.clim_method``; the native daily Z gets the daily-minus-monthly offset);
* ``ssw_catalogue``: SSW central dates and per-event properties, ERA5 and every model member.

The caches themselves are built by ``9a_strat_build_caches`` (``mode="auto"``); the other notebooks
read them with ``mode="require"``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import xarray as xr

from util import rtenso_case as rc
from util import strat_case as sc

MONTH_NAMES = {1: "Jan", 2: "Feb", 3: "Mar", 4: "Apr", 5: "May", 6: "Jun", 7: "Jul", 8: "Aug",
               9: "Sep", 10: "Oct", 11: "Nov", 12: "Dec"}
INIT_COLORS = {2: "#1f77b4", 5: "#2ca02c", 8: "#9467bd", 11: "#d62728"}      # as notebook 8
DJF = ((0, 12), (1, 1), (1, 2))                                   # (years after the winter year, month)
DJF_DAYS = np.array([31.0, 31.0, 28.0])                          # day weights (no-leap model calendar)


def default_pool(prefix=rc.MODEL_PREFIX):
    return {"Reanalysis-CLIM": f"{prefix}_RTENSO2026_clim",
            "Reanalysis-OBS": f"{prefix}_RTENSO2026_obs",
            "Reanalysis-NSIDC": f"{prefix}_RTENSO2026_obs_ice_nsidc"}


@dataclass
class StratConfig:
    """Run controls shared by the 9a-9h notebooks."""
    era5_zonal_dir: Path                       # ERA5/1x1/daily_zonal
    era5_daily_dir: Path                       # ERA5/1x1/daily (PSL)
    simulation_root: Path                      # raw eam.h1 history (cache keys; read only by 9a)
    pool: dict = field(default_factory=default_pool)
    year: int = 2026
    init_months: tuple = (2, 5, 8, 11)
    clim_method: str = "direct"                # "direct", "lead_drift" or "target_month"
    clim_years: tuple = rc.CLIM_YEARS
    era5_first_year: int = 1979
    ref_level: float = 10.0
    mode: str = "require"                      # "auto" builds missing caches (9a); "require" reads only
    workers: int = 4
    ssw: dict = field(default_factory=lambda: {"sep_days": 20, "recover_days": 10})


# ---------------------------------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------------------------------
def doy_of(da):
    return xr.DataArray(sc.noleap_doy(da["time"].values), dims="time", coords={"time": da["time"]})


def on_days(clim, da):
    """A (doy, ...) climatology on the days of ``da``."""
    return clim.sel(doy=doy_of(da)).drop_vars("doy")


def mid_month(times):
    """The 15th of each time's month (datetime64[ns])."""
    return sc.to_datetime64(times).astype("datetime64[M]").astype("datetime64[ns]") + np.timedelta64(14, "D")


def case_of(member):
    return str(member).split(":")[0]


def source_of(prefix):
    return prefix.split(f"{rc.MODEL_PREFIX}_")[-1]


def target_time(year, init, nlead=rc.NLEAD):
    """Mid-month target dates of each lead of a start."""
    out = []
    for lead in range(nlead):
        m = init - 1 + lead
        out.append(np.datetime64(f"{year + m // 12:04d}-{m % 12 + 1:02d}-15"))
    return np.asarray(out, dtype="datetime64[ns]")


def djf_mean(series, w):
    """Day-weighted DJF mean of a monthly series stamped on the 15th (dims other than time kept);
    None unless Dec, Jan and Feb are all present."""
    times = [np.datetime64(f"{w + dy}-{m:02d}-15") for dy, m in DJF]
    if not all(t in series["time"].values for t in times):
        return None
    v = series.sel(time=times)
    return (v * xr.DataArray(DJF_DAYS, dims="time", coords={"time": v["time"]})).sum("time", skipna=False) / DJF_DAYS.sum()


def winter_label(w):
    return f"{w}/{str(w + 1)[2:]}"


# ---------------------------------------------------------------------------------------------------
# ERA5
# ---------------------------------------------------------------------------------------------------
def _era5_stamp(era5_zonal_dir):
    files = list(Path(era5_zonal_dir).glob("ERA5.daily_zonal.en00.*.nc"))
    return int(max(f.stat().st_mtime for f in files)) if files else 0


def load_era5(cfg: StratConfig) -> SimpleNamespace:
    """ERA5 daily zonal-mean indices and climatologies (``cfg.clim_years``).

    u60 (time): U(60N, ref level); pcz (time, plev): polar-cap height; pct: polar-cap T at the ref
    level; vt: 100-hPa eddy heat flux 45-75N; qbo (time, plev): monthly U 5S-5N; pch: polar-cap
    height anomaly / daily sigma; month_clim (month, plev, lat) of U, T, Z; clim_* and std_* (doy, ...);
    fields: the lazy daily zonal means. The reductions are cached under ``S2D_DIAG_ROOT/observations/strat``
    and rebuilt when the archive files change (``mode="rebuild"`` forces it).
    """
    era5 = sc.era5_daily_zonal(cfg.era5_zonal_dir, years=(cfg.era5_first_year, cfg.year), load=False)
    o = SimpleNamespace(fields=era5, obs_end=pd.Timestamp(era5["time"].values[-1]))
    path = rc.diag_dir("observations", "strat", f"ERA5_daily_zonal_indices_y{cfg.era5_first_year}-{cfg.year}"
                       f"_clim{cfg.clim_years[0]}-{cfg.clim_years[1]}_p{cfg.ref_level:g}_v{sc.CACHE_VERSION}.nc")
    expected = {"source_mtime": _era5_stamp(cfg.era5_zonal_dir), "ref_level": cfg.ref_level,
                "clim_years": f"{cfg.clim_years[0]}-{cfg.clim_years[1]}"}
    if cfg.mode != "rebuild" and sc._cache_ok(path, expected):
        ds = xr.load_dataset(path)
        o.u60, o.pct, o.vt, o.pcz = ds["u60"], ds["pct"], ds["vt"], ds["pcz"]
        o.qbo = ds["qbo"].rename(time_month="time")
        o.month_clim = xr.Dataset({k: ds[f"month_clim_{k}"] for k in ("U", "T", "Z")})
    else:
        o.u60 = sc.lat_mean(era5["U"].sel(plev=cfg.ref_level), sc.U60).load()
        o.pcz = sc.lat_mean(era5["Z"], sc.PCAP).load()
        o.pct = sc.lat_mean(era5["T"].sel(plev=cfg.ref_level), sc.PCAP).load()
        o.vt = sc.lat_mean(era5["VT"].sel(plev=100.0), sc.VT_BAND).load()
        o.qbo = sc.monthly_means(sc.lat_mean(era5["U"], sc.EQ)).load()
        o.month_clim = sc.obs_month_climatology(era5[["U", "T", "Z"]], cfg.clim_years).load()
        ds = xr.Dataset({"u60": o.u60, "pct": o.pct, "vt": o.vt, "pcz": o.pcz, "qbo": o.qbo.rename(time="time_month"),
                         **{f"month_clim_{k}": o.month_clim[k] for k in ("U", "T", "Z")}})
        ds.attrs.update(expected)
        sc.atomic_to_netcdf(ds, path)
    o.clim_u60, o.std_u60 = sc.daily_climatology(o.u60, cfg.clim_years)
    o.clim_pcz, o.std_pcz = sc.daily_climatology(o.pcz, cfg.clim_years)
    o.clim_vt, _ = sc.daily_climatology(o.vt, cfg.clim_years)
    o.clim_pct, _ = sc.daily_climatology(o.pct, cfg.clim_years)
    o.pch = (o.pcz - on_days(o.clim_pcz, o.pcz)) / on_days(o.std_pcz, o.pcz)
    o.u60_month = sc.monthly_means(o.u60)
    print(f"ERA5 daily zonal means {pd.Timestamp(era5['time'].values[0]):%d %b %Y} - {o.obs_end:%d %b %Y}")
    return o


def observed_nino(cfg: StratConfig, product="HadISST1.1"):
    """Observed monthly Nino3.4 anomaly (mid-month stamps) and its day-weighted DJF mean per winter (Series)."""
    nino = rc.box_mean(rc.obs_anomaly("TS", years=(cfg.era5_first_year, cfg.year), mode="auto", product=product),
                       rc.NINO_REGIONS["Nino3.4"]).load()
    nino = nino.assign_coords(time=mid_month(nino["time"].values))
    djf = {}
    for w in range(cfg.era5_first_year, cfg.year):
        v = djf_mean(nino, w)
        djf[w] = np.nan if v is None else float(v)
    return nino, pd.Series(djf)


# ---------------------------------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------------------------------
def load_references(cfg: StratConfig) -> SimpleNamespace:
    """Reanalysis monthly zonal means (year, member, lead, plev, lat) of every start month."""
    years = range(cfg.clim_years[0], cfg.clim_years[1] + 1)
    ref = {m: sc.monthly_zonal(rc.REFERENCE["case_prefix"], "Reanalysis", years, m, mode=cfg.mode, workers=cfg.workers,
                               verbose=cfg.mode != "require")
           for m in rc.REFERENCE["init_months"]}
    direct = {m: sc.monthly_zonal(rc.REFERENCE["case_prefix"], "Reanalysis", years, m, mode=cfg.mode, workers=cfg.workers,
                                  verbose=cfg.mode != "require")
              for m in rc.DIRECT_INIT_MONTHS}
    out = SimpleNamespace(ref=ref, direct=direct, all={**ref, **direct})
    out.plev = ref[rc.REFERENCE["init_months"][0]]["plev"].values
    out.for_method = lambda method: out.all if method == "direct" else out.ref
    return out


def pooled_monthly(cfg: StratConfig, init):
    """Pooled monthly zonal means of one start with target-month ``time`` on ``lead``."""
    out = sc.pooled_monthly_zonal(cfg.pool, cfg.year, init, mode=cfg.mode, workers=cfg.workers, verbose=False)
    return out.assign_coords(time=("lead", target_time(cfg.year, init, out.sizes["lead"])))


def load_model(cfg: StratConfig, era5: SimpleNamespace, refs: SimpleNamespace) -> SimpleNamespace:
    """Pooled 2026 starts: monthly and daily zonal means, lead-time climatologies and daily/monthly indices.

    Daily (per start, (member, time)): u60 raw, u60_clim (the start's climatology on the same days),
    u60_adj (raw minus the model-minus-ERA5 climatological difference), pct (raw polar-cap T at the ref
    level), vt (100-hPa heat flux), pch (polar-cap height anomaly / ERA5 daily sigma, (member, time, plev)).
    Monthly (per start): mm_u60 raw and mm_u60_anom (member, lead), mm_qbo (member, lead, plev).
    The native daily Z gets the ensemble-mean daily-minus-monthly offset (``sc.daily_monthly_offset``).
    """
    month_clim = era5.month_clim.sel(plev=refs.plev)
    m = SimpleNamespace(monthly={}, daily={}, lead_clim={}, z_offset={}, u60={}, u60_clim={}, u60_adj={}, pct={},
                        vt={}, pch={}, mm_u60={}, mm_u60_anom={}, mm_qbo={})
    std_pcz = sc.interp_logp(era5.std_pcz, sc.DAILY_LEVELS)
    for init in cfg.init_months:
        mo = pooled_monthly(cfg, init)
        d = sc.pooled_daily_zonal(cfg.simulation_root, cfg.pool, cfg.year, init, mode=cfg.mode, workers=cfg.workers,
                                  verbose=False)
        lc = sc.lead_climatology(init, refs.for_method(cfg.clim_method), month_clim, method=cfg.clim_method)
        m.monthly[init], m.daily[init], m.lead_clim[init] = mo, d, lc
        lcd = sc.interp_logp(lc[["U", "Z", "T"]], sc.DAILY_LEVELS)
        m.z_offset[init] = sc.daily_monthly_offset(d, mo, cfg.year, init)
        lcd["Z"] = lcd["Z"] + m.z_offset[init]
        cd = sc.daily_from_lead(lcd, cfg.year, init, d["time"].values).assign_coords(time=d["time"])
        m.u60[init] = sc.lat_mean(d["U"].sel(plev=cfg.ref_level), sc.U60)
        m.u60_clim[init] = sc.lat_mean(cd["U"].sel(plev=cfg.ref_level), sc.U60)
        m.u60_adj[init] = m.u60[init] - (m.u60_clim[init] - on_days(era5.clim_u60, m.u60_clim[init]))
        m.pct[init] = sc.lat_mean(d["T"].sel(plev=cfg.ref_level), sc.PCAP)
        m.vt[init] = sc.lat_mean(d["VT"].sel(plev=100.0), sc.VT_BAND)
        pcz = sc.lat_mean(d["Z"], sc.PCAP) - sc.lat_mean(cd["Z"], sc.PCAP)
        m.pch[init] = (pcz / on_days(std_pcz, pcz)).dropna("plev", how="all")
        m.mm_u60[init] = sc.lat_mean(mo["U"].sel(plev=cfg.ref_level), sc.U60)
        m.mm_u60_anom[init] = m.mm_u60[init] - sc.lat_mean(lc["U"].sel(plev=cfg.ref_level), sc.U60)
        m.mm_qbo[init] = sc.lat_mean(mo["U"], sc.EQ)
        print(f"[init {MONTH_NAMES[init]} {cfg.year}] {d.sizes['member']} members, daily "
              f"{pd.Timestamp(d['time'].values[0]):%d %b %Y} - {pd.Timestamp(d['time'].values[-1]):%d %b %Y}; "
              f"climatology: {lc.attrs.get('climatology_method', cfg.clim_method)}")
    return m


def offset_table(model: SimpleNamespace, levels=(10.0, 50.0, 100.0, 500.0, 1000.0)):
    """Polar-cap daily-minus-monthly Z offset by start: mean and standard deviation over leads (m)."""
    rows = []
    for init, off in model.z_offset.items():
        for p in levels:
            v = sc.lat_mean(off.sel(plev=p), sc.PCAP)
            rows.append({"start": MONTH_NAMES[init], "plev (hPa)": p, "mean (m)": float(v.mean()),
                         "std over leads": float(v.std()), "min": float(v.min()), "max": float(v.max())})
    return pd.DataFrame(rows).pivot_table(index="plev (hPa)", columns="start", values=["mean (m)", "std over leads"])


# ---------------------------------------------------------------------------------------------------
# SSW catalogue
# ---------------------------------------------------------------------------------------------------
def ssw_catalogue(cfg: StratConfig, era5: SimpleNamespace, model: SimpleNamespace) -> SimpleNamespace:
    """SSW central dates of ERA5 and of every member (raw and bias-adjusted wind), with per-event details.

    obs_events / model_ssw[init] / model_ssw_adj[init]: DataFrames (series, winter, central_date);
    events (all, with min_u, easterly_days, recovery_date, days_to_recovery, source, wind, case).
    """
    kw = cfg.ssw
    c = SimpleNamespace(model_ssw={}, model_ssw_adj={})
    c.obs_events = sc.ssw_table(era5.u60.values, era5.u60["time"].values, "ERA5", **kw)
    c.obs_winters = sc.covered_winters(era5.u60["time"].values)
    c.obs_counts = pd.Series(0, index=c.obs_winters).add(c.obs_events.groupby("winter").size(), fill_value=0).astype(int)
    cat = [sc.ssw_event_details(era5.u60.values, era5.u60["time"].values, c.obs_events["central_date"], kw["recover_days"])
           .assign(source="ERA5", series="ERA5", wind="ERA5")]
    for init in cfg.init_months:
        for kind, src, store in (("raw", model.u60[init], c.model_ssw), ("bias-adjusted", model.u60_adj[init], c.model_ssw_adj)):
            tabs = [sc.ssw_table(src.sel(member=mem).values, src["time"].values, str(mem), **kw) for mem in src["member"].values]
            store[init] = pd.concat(tabs, ignore_index=True)
            for mem, ev in store[init].groupby("series"):
                cat.append(sc.ssw_event_details(src.sel(member=mem).values, src["time"].values, ev["central_date"],
                                                kw["recover_days"]).assign(source=f"E3SM init{init:02d}", series=mem, wind=kind))
    ev = pd.concat(cat, ignore_index=True)
    ev["winter"] = [sc.winter_of(d) for d in ev["central_date"]]
    ev["case"] = [case_of(s) for s in ev["series"]]
    c.events = ev
    return c


def member_counts(cfg: StratConfig, model: SimpleNamespace, catalogue: SimpleNamespace, kind="raw"):
    """{(init, winter): Series of SSW counts per member} for every fully covered winter of each start."""
    tabs = catalogue.model_ssw if kind == "raw" else catalogue.model_ssw_adj
    out = {}
    for init in cfg.init_months:
        members = [str(mm) for mm in model.u60[init]["member"].values]
        for w in sc.covered_winters(model.u60[init]["time"].values):
            t = tabs[init]
            out[(init, w)] = t[t["winter"] == w].groupby("series").size().reindex(members, fill_value=0)
    return out


# ---------------------------------------------------------------------------------------------------
# DJF vortex strength
# ---------------------------------------------------------------------------------------------------
def djf_vortex(cfg: StratConfig, model: SimpleNamespace, refs: SimpleNamespace, era5: SimpleNamespace, method=None):
    """{(init, winter): DJF U(60N, ref level) anomaly per member (monthly product)} with the ``method``
    climatology (default ``cfg.clim_method``)."""
    method = method or cfg.clim_method
    month_clim = era5.month_clim.sel(plev=refs.plev)
    out = {}
    for init in cfg.init_months:
        raw = model.mm_u60[init].assign_coords(time=model.monthly[init]["time"]).swap_dims(lead="time")
        lc = sc.lead_climatology(init, refs.for_method(method), month_clim, method=method)
        clim = sc.lat_mean(lc["U"].sel(plev=cfg.ref_level), sc.U60).assign_coords(
            time=("lead", model.monthly[init]["time"].values)).swap_dims(lead="time")
        for w in (cfg.year, cfg.year + 1):
            r, c = djf_mean(raw, w), djf_mean(clim, w)
            if r is not None and c is not None:
                out[(init, w)] = r - c
    return out
