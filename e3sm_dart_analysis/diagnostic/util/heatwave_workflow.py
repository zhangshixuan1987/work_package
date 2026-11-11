"""Workflow of the summer 2012 heat-wave notebooks (``hind_lead_eval/06a-06h``).

``HeatwaveCase`` holds one set of ``HeatwaveSettings`` and gives every notebook the same products,
each computed once and cached on disk with the settings it depends on:

* ``obs``: ERA5/GPCP fields on the map box, both climatologies and the calendar-day thresholds
  (cached, ``observed_fields.nc``), and from them the observed events, event window, metrics, coverage;
* ``field(config, variable)``: daily member fields on the map box (cached, ``box_<CONFIG>_<SEASON>_<VAR>.nc``);
* ``series``: regional series of every configuration, start and member (from the two above, in memory);
* ``metrics``: heat metrics of every member on the event window (cached, ``heat_metrics.csv``).

A product is read when its file exists with the same settings, unless ``force_compute`` (this case's
results) or ``force_inputs`` (also the model and observation reads) is set. Methods: ``util/heatwave_case.py``.
"""
from __future__ import annotations

from functools import cached_property
from types import SimpleNamespace

import numpy as np
import pandas as pd
import xarray as xr

from configs.fcst_init_experiment_config import balanced_runs, case_stamp, check_configs, init_dates, readiness
from configs.fcst_init_experiment_config import members as config_members
from configs.hind_lead_eval_cases import HeatwaveSettings
from configs.output_paths import workflow_output_dirs
from util.fcst_init_data import region_masks
from util.heatwave_case import (bootstrap_contrast, calendar, climatologies, coverage_table, crps_ensemble,
                                heat_metrics, lead_window_scores, model_box, observed_box, observed_events,
                                percentile_threshold, region_mean, region_weights, runs)
from util.hind_case_io import cached_dataset, cached_table

SUBDIR = "heatwave_2012"
VARIABLES = ("TREFHT", "PRECT")
REFERENCE_LABELS = {"none": "observed-reference anomalies, no mean-bias adjustment (differences = absolute error)",
                    "amip": "AMIP-reference anomalies (AMIP mean bias removed; drift kept)"}
METRIC_UNITS = dict(event_mean_T_error="K", event_mean_anomaly="K", peak_anomaly="K", peak_date_error_days="days",
                    exceedance_days="days", longest_run_days="days", accumulated_exceedance="K day",
                    area_fraction_mean="0-1", event_precip_mm="mm", peak_anomaly_tol="K",
                    peak_date_error_tol_days="days", longest_run_tol_days="days", crps_event_mean="K")
PRIMARY_METRICS = ["event_mean_T_error", "event_mean_anomaly", "peak_anomaly", "peak_date_error_days", "exceedance_days",
                   "longest_run_days", "accumulated_exceedance", "area_fraction_mean", "event_precip_mm"]
TOLERANT_METRICS = ["peak_anomaly_tol", "peak_date_error_tol_days", "longest_run_tol_days"]
START_METRICS = ["event_mean_T_error", "event_mean_anomaly", "exceedance_days", "accumulated_exceedance",
                 "area_fraction_mean", "event_precip_mm"]
CONTRAST_METRICS = ["event_mean_anomaly", "exceedance_days", "event_precip_mm", "crps_event_mean"]
SERIES_YLABELS = {"TREFHT": "TREFHT anomaly (K)", "PRECT": "precipitation since start (mm)"}

# settings each product depends on
_OBS_KEYS = ("map_box", "region", "obs_period", "clim_years", "exclude_years", "amip_cases", "percentile",
             "half_window")
_EVENT_KEYS = _OBS_KEYS + ("min_event_days", "event_window", "peak_tolerance_days")
_MODEL_KEYS = ("configs", "season", "n_days", "max_members", "map_box")


class HeatwaveCase:
    """Products of the heat-wave case for one set of settings (see the module docstring)."""

    def __init__(self, settings: HeatwaveSettings, *, force_compute: bool = False, force_inputs: bool = False):
        self.cfg = settings
        self.force_inputs = force_inputs
        self.force = force_compute or force_inputs
        self._fields = {}
        check_configs(list(settings.configs))
        p = settings.paths
        self.data_dir, self.figure_dir = workflow_output_dirs(p.workflow_name, SUBDIR, output_root=p.diagnostic_output_root)
        self.scored = [c for c in settings.configs if c not in settings.ic_defect]   # scores, rankings, contrasts
        self.primary = [c for c in settings.primary if c in self.scored]
        self.inits = list(init_dates([settings.season])[settings.season])
        self.members = {c: tuple(config_members(c, settings.season))[:settings.max_members] for c in settings.configs}
        self.years = [y for y in range(settings.clim_years[0], settings.clim_years[1] + 1)
                      if y not in settings.exclude_years]
        r = settings.region
        self.label = (f"{r['name']} land ({r['lat'][0]:g}-{r['lat'][1]:g}N, "
                      f"{360 - r['lon'][1]:g}-{360 - r['lon'][0]:g}W)")
        print(f"[DATA] {self.data_dir}\n[FIGURE] {self.figure_dir}\n[SCORED] {self.scored}\n"
              f"[QC ONLY] {list(settings.ic_defect)}")

    def stamp(self, *keys: str, extra: str = "") -> str:
        return self.cfg.stamp(*keys) + extra + case_stamp(self.cfg.configs)

    @property
    def metrics_stamp(self) -> str:
        return self.stamp(*sorted(set(_EVENT_KEYS + _MODEL_KEYS + ("biases",))),
                          extra=str({c: m for c, m in self.members.items()}))

    # ------------------------------------------------------------------
    # Inputs
    # ------------------------------------------------------------------
    @cached_property
    def runs(self) -> list:
        """(config, season) pairs whose members are post-processed (raises when one is missing)."""
        status = readiness(list(self.cfg.configs), [self.cfg.season], streams=(("atm", "daily"), ("atm", "6hourly")),
                           data_root=self.cfg.paths.input_data_root)
        out = balanced_runs(status, self.cfg.require_all_members)
        print(f"[RUNS] {out}\n[STARTS] {[str(i) for i in self.inits]}\n"
              f"[MEMBERS] {({c: len(m) for c, m in self.members.items()})}")
        return out

    def _observed_fields(self) -> xr.Dataset:
        cfg, p = self.cfg, self.cfg.paths
        obs_T = observed_box("TREFHT", *cfg.obs_period, cfg.map_box, observation_root=p.observation_root)
        obs_P = observed_box("PRECT", *cfg.obs_period, cfg.map_box, observation_root=p.observation_root)
        mask = region_masks([cfg.region], obs_T.lat.values, obs_T.lon.values,
                            land_support_dir=p.land_support_dir).isel(region=0, drop=True)
        clim = climatologies("TREFHT", *cfg.obs_period, self.years, cfg.map_box, observation_root=p.observation_root,
                             amip_root=p.amip_root, amip_cases=cfg.amip_cases, land_support_dir=p.land_support_dir)
        thr = percentile_threshold("TREFHT", *cfg.obs_period, self.years, cfg.map_box, region_weights(mask),
                                   q=cfg.percentile, half_window=cfg.half_window, observation_root=p.observation_root)
        ds = xr.merge([obs_T.rename("obs_T"), obs_P.rename("obs_P"), mask.rename("mask").astype("int8"),
                       clim.clim_obs, clim.clim_amip, thr.thr_grid, thr.thr_region], join="exact",
                      combine_attrs="drop")
        ds.attrs.update(clim_obs=clim.attrs["clim_obs"], clim_amip=clim.attrs["clim_amip"])
        return ds

    @cached_property
    def obs(self) -> SimpleNamespace:
        """Observed fields (cached) and everything derived from them: regional series, events, window,
        observed metrics and outcome, simulation coverage."""
        cfg = self.cfg
        ds = cached_dataset(self.data_dir / "observed_fields.nc", self._observed_fields, self.stamp(*_OBS_KEYS),
                            force=self.force_inputs)
        mask = ds.mask.astype(bool)
        w = region_weights(mask)
        o = SimpleNamespace(T=ds.obs_T, P=ds.obs_P, clim=ds[["clim_obs", "clim_amip"]], mask=mask, w=w,
                            times=ds.time.values)
        o.s = region_mean(ds.obs_T, w)                                  # absolute regional-mean T (K)
        o.c_reg = {"none": region_mean(ds.clim_obs, w), "amip": region_mean(ds.clim_amip, w)}
        o.a = o.s - o.c_reg["none"]
        o.thr_region = ds.thr_region
        o.thr_anom = ds.thr_region - o.c_reg["none"]                    # regional threshold as an anomaly
        o.thr_anom_grid = ds.thr_grid - ds.clim_obs                     # gridpoint threshold as an anomaly
        o.area = region_mean((ds.obs_T > ds.thr_grid).astype("float64"), w)
        o.P_reg = region_mean(ds.obs_P, w)
        o.events = observed_events(o.s, ds.thr_region, o.a, cfg.min_event_days)
        if cfg.event_window is None:
            cand = o.events[pd.to_datetime(o.events.onset) >= pd.Timestamp("2012-06-01")]
            best = cand.loc[cand.accumulated_exceedance.idxmax()]
            o.window, o.peak = (str(best.onset), str(best.end)), str(best.peak_date)
        else:
            o.window = tuple(cfg.event_window)
            seg = o.a.sel(time=slice(*o.window))
            o.peak = str(pd.Timestamp(seg.time.values[int(np.argmax(seg.values))]).date())
        m = heat_metrics(o.a.values, o.times, o.thr_anom.values, o.area.values, o.window, o.peak,
                         cfg.peak_tolerance_days)
        m["event_mean_T"] = float(o.s.sel(time=slice(*o.window)).mean())
        m["event_mean_T_error"] = 0.0
        m["event_precip_mm"] = float(o.P_reg.sel(time=slice(*o.window)).sum())
        o.metrics = m
        o.y = int(m["longest_run_days"] >= cfg.min_event_days)          # observed outcome of the event definition
        o.y_tol = int(m["longest_run_tol_days"] >= cfg.min_event_days)
        o.coverage = coverage_table(self.inits, cfg.n_days, o.window)
        o.event_leads = {r.init: (r.lead_at_onset, r.lead_at_end) for r in o.coverage.itertuples()}
        o.clim_sources = (ds.attrs.get("clim_obs", ""), ds.attrs.get("clim_amip", ""))
        print(f"[EVENT] window {o.window[0]} to {o.window[1]}, observed peak {o.peak}")
        return o

    def field(self, config: str, variable: str) -> xr.DataArray:
        """(init, member, lead, lat, lon) daily field of ``config`` on the map box (cached)."""
        key = (config, variable)
        if key not in self._fields:
            cfg = self.cfg
            cache = self.data_dir / f"box_{config}_{cfg.season}_{variable}.nc"
            if self.force_inputs or not cache.exists():
                self.runs   # every configuration must be post-processed before a read
            da = model_box(config, cfg.season, variable, self.members[config], cfg.map_box, n_days=cfg.n_days,
                           cache=cache, force=self.force_inputs, workers=cfg.workers,
                           data_root=cfg.paths.input_data_root, land_support_dir=cfg.paths.land_support_dir)
            if not (np.allclose(da.lat, self.obs.T.lat) and np.allclose(da.lon, self.obs.T.lon)):
                raise ValueError(f"model ({config} {variable}) and observation grids differ on the map box")
            print(f"[MODEL] {config} {variable}: {dict(da.sizes)}")
            self._fields[key] = da
        return self._fields[key]

    def on_calendar(self, da: xr.DataArray, init) -> xr.DataArray:
        """(member, time, lat, lon) of one start on its calendar days, within the observed period."""
        x = da.sel(init=np.datetime64(init, "ns")).drop_vars("init", errors="ignore")
        x = x.assign_coords(lead=calendar(init, x.lead.values).astype("datetime64[ns]")).rename(lead="time")
        return x.sel(time=slice(*self.cfg.obs_period))

    # ------------------------------------------------------------------
    # Derived products
    # ------------------------------------------------------------------
    @cached_property
    def series(self) -> SimpleNamespace:
        """Regional means of every member and start, all configurations: absolute temperature ``absT``,
        anomalies against both references ``anom[bias, c, init]``, area fraction above the gridpoint threshold
        ``area[bias, c, init]``, precipitation ``prect`` (daily) and ``precip`` (accumulated from the start)."""
        o, s = self.obs, SimpleNamespace(absT={}, anom={}, area={}, prect={}, precip={})
        for c in self.cfg.configs:
            for init in self.inits:
                T = self.on_calendar(self.field(c, "TREFHT"), init)
                t = T.time.values
                s.absT[c, init] = region_mean(T, o.w)
                s.prect[c, init] = region_mean(self.on_calendar(self.field(c, "PRECT"), init), o.w)
                s.precip[c, init] = s.prect[c, init].cumsum("time")
                for b in self.cfg.biases:
                    anom_grid = T - (o.clim.clim_amip if b == "amip" else o.clim.clim_obs).sel(time=t)
                    s.anom[b, c, init] = region_mean(anom_grid, o.w)
                    s.area[b, c, init] = region_mean((anom_grid > o.thr_anom_grid.sel(time=t)).astype("float64"), o.w)
        s.obs_cum = {init: o.P_reg.sel(time=slice(np.datetime64(init, "ns"), None)).cumsum("time") for init in self.inits}
        print(f"[SERIES] {len(s.absT)} (config, start) pairs, {len(s.anom)} anomaly series")
        return s

    def series_dataset(self) -> xr.Dataset:
        """The regional series on the calendar days of the observed period (NaN where a start has no
        simulated day), dimensions (config, init, member_index, time); member labels in ``member``."""
        days = pd.date_range(*self.cfg.obs_period, freq="D").values
        s, configs = self.series, list(self.cfg.configs)

        def stack(d, prefix=()):
            per_c = []
            for c in configs:
                per_i = [d[(*prefix, c, i)].reindex(time=days).drop_vars("member")
                         .rename(member="member_index") for i in self.inits]
                per_c.append(xr.concat(per_i, pd.Index(np.array(self.inits, dtype="datetime64[ns]"), name="init")))
            return xr.concat(per_c, pd.Index(configs, name="config"), join="outer")

        ds = xr.Dataset({"absT": stack(s.absT), "prect": stack(s.prect), "precip": stack(s.precip),
                         **{f"anom_{b}": stack(s.anom, (b,)) for b in self.cfg.biases},
                         **{f"area_{b}": stack(s.area, (b,)) for b in self.cfg.biases}})
        n = max(len(m) for m in self.members.values())
        labels = np.array([[m[k] if k < len(m) else "" for k in range(n)] for m in
                           (self.members[c] for c in configs)])
        ds["member"] = (("config", "member_index"), labels)
        ds["obs_anom"], ds["obs_precip"] = self.obs.a, self.obs.P_reg
        ds.attrs.update(region=self.label, units_T="K", units_P="mm/day; precip accumulated from the start (mm)")
        return ds

    def _metrics(self) -> pd.DataFrame:
        o, s, rows = self.obs, self.series, []
        for c in self.cfg.configs:
            for init in self.inits:
                t = s.absT[c, init].time.values
                th = o.thr_anom.sel(time=t).values
                T_event = s.absT[c, init].sel(time=slice(*o.window)).mean("time")
                p_event = s.prect[c, init].sel(time=slice(*o.window)).sum("time")
                for b in self.cfg.biases:
                    for k, m in enumerate(s.absT[c, init].member.values):
                        r = heat_metrics(s.anom[b, c, init].isel(member=k).values, t, th,
                                         s.area[b, c, init].isel(member=k).values, o.window, o.peak,
                                         self.cfg.peak_tolerance_days)
                        rows.append(dict(bias=b, config=c, init=str(init), member=str(m), **r,
                                         event_mean_T=float(T_event.isel(member=k)),
                                         event_precip_mm=float(p_event.isel(member=k))))
        metrics = pd.DataFrame(rows)
        metrics["event_mean_T_error"] = metrics.event_mean_T - o.metrics["event_mean_T"]
        return metrics

    @cached_property
    def metrics(self) -> pd.DataFrame:
        """Heat metrics of every member, start, configuration and bias treatment on the event window (cached)."""
        return cached_table(self.data_dir / "heat_metrics.csv", self._metrics, self.metrics_stamp, force=self.force)

    @property
    def scored_metrics(self) -> pd.DataFrame:
        return self.metrics[self.metrics.config.isin(self.scored)]

    # ------------------------------------------------------------------
    # Scores of the individual notebooks (cached)
    # ------------------------------------------------------------------
    def init_qc(self) -> pd.DataFrame:
        """Ensemble-mean error of regional-mean TREFHT on lead day 1 and in week 1, every configuration."""
        def build():
            rows = []
            for c in self.cfg.configs:
                for init in self.inits:
                    x = self.series.absT[c, init]
                    e = x.mean("member") - self.obs.s.sel(time=x.time)
                    rows.append(dict(config=c, init=str(init), lead1_error=float(e.isel(time=0)),
                                     week1_error=float(e.isel(time=slice(0, 7)).mean()),
                                     status="IC defect: QC only" if c in self.cfg.ic_defect else "scored"))
            return pd.DataFrame(rows)
        stamp = self.stamp(*sorted(set(_OBS_KEYS + _MODEL_KEYS + ("ic_defect",))), extra=str(self.members))
        return cached_table(self.data_dir / "init_qc.csv", build, stamp, force=self.force)

    def lead_scores(self) -> pd.DataFrame:
        """Ensemble-mean error, spread and CRPS of the lead-window means of every configuration and start."""
        def build():
            rows = []
            for b in self.cfg.biases:
                for c in self.cfg.configs:
                    for init in self.inits:
                        for r in lead_window_scores(self.series.anom[b, c, init], self.obs.a, init, self.cfg.lead_windows):
                            rows.append(dict(bias=b, config=c, init=str(init), **r))
            return pd.DataFrame(rows)
        stamp = self.metrics_stamp + self.cfg.stamp("lead_windows")
        return cached_table(self.data_dir / "lead_window_scores.csv", build, stamp, force=self.force)

    def event_scores(self) -> pd.DataFrame:
        """Per start: CRPS of the event-window mean anomaly, event probability and Brier score (exact and
        timing-tolerant)."""
        o = self.obs

        def build():
            rows = []
            for (b, c, i), g in self.metrics.groupby(["bias", "config", "init"], sort=False):
                p = float((g.longest_run_days >= self.cfg.min_event_days).mean())
                p_tol = float((g.longest_run_tol_days >= self.cfg.min_event_days).mean())
                rows.append(dict(bias=b, config=c, init=i, n_members=len(g),
                                 crps_event_mean=crps_ensemble(g.event_mean_anomaly, o.metrics["event_mean_anomaly"]),
                                 ens_mean_error=float(g.event_mean_anomaly.mean() - o.metrics["event_mean_anomaly"]),
                                 p_event=p, brier=(p - o.y) ** 2, p_event_tol=p_tol, brier_tol=(p_tol - o.y_tol) ** 2))
            return pd.DataFrame(rows)
        return cached_table(self.data_dir / "event_scores.csv", build, self.metrics_stamp, force=self.force)

    def contrasts(self) -> pd.DataFrame:
        """Differences between two scored experiments (a minus b) per start, with a 90 % bootstrap interval
        from resampling the members of each ensemble independently."""
        pairs = [(a, b) for a, b in self.cfg.contrasts if a in self.scored and b in self.scored]
        o = self.obs.metrics["event_mean_anomaly"]

        def build():
            rows = []
            for bias in self.cfg.biases:
                mb = self.metrics[self.metrics.bias == bias]
                for a, b in pairs:
                    for i in map(str, self.inits):
                        A, B = (mb[(mb.config == x) & (mb.init == i)] for x in (a, b))
                        for m in CONTRAST_METRICS:
                            if m == "crps_event_mean":
                                r = bootstrap_contrast(A.event_mean_anomaly, B.event_mean_anomaly,
                                                       stat=lambda x: crps_ensemble(x, o), n_boot=self.cfg.n_boot)
                            else:
                                r = bootstrap_contrast(A[m], B[m], n_boot=self.cfg.n_boot)
                            rows.append(dict(bias=bias, contrast=f"{a} - {b}", init=i, metric=m, **r))
            return pd.DataFrame(rows)
        stamp = self.metrics_stamp + self.cfg.stamp("contrasts", "n_boot", "ic_defect")
        print(f"[PAIRS] {pairs}")
        return cached_table(self.data_dir / "contrasts.csv", build, stamp, force=self.force)

    def threshold_sensitivity(self) -> pd.DataFrame:
        """Observed hot days and longest run, and the event probability of each scored experiment and start,
        with the threshold recomputed from other reference periods and percentiles (fixed event window)."""
        cfg, o = self.cfg, self.obs

        def build():
            rows = []
            for y0, y1 in cfg.sens_threshold_years:
                ys = [y for y in range(y0, y1 + 1) if y not in cfg.exclude_years]
                for q in cfg.sens_percentiles:
                    th = percentile_threshold("TREFHT", *cfg.obs_period, ys, cfg.map_box, o.w, q=q,
                                              half_window=cfg.half_window,
                                              observation_root=cfg.paths.observation_root).thr_region
                    hot = (o.s > th).sel(time=slice(*o.window)).values
                    base = dict(threshold_years=f"{y0}-{y1}", percentile=q, obs_hot_days=int(hot.sum()),
                                obs_longest_run=max((e - s + 1 for s, e in runs(hot)), default=0))
                    for c in self.scored:
                        for init in self.inits:
                            x = self.series.absT[c, init]
                            hot_m = (x > th.sel(time=x.time)).sel(time=slice(*o.window)).values   # (member, time)
                            longest = [max((e - s + 1 for s, e in runs(h)), default=0) for h in hot_m]
                            rows.append(dict(base, config=c, init=str(init),
                                             p_event=float(np.mean(np.array(longest) >= cfg.min_event_days))))
            return pd.DataFrame(rows)
        stamp = self.metrics_stamp + cfg.stamp("sens_threshold_years", "sens_percentiles", "ic_defect")
        return cached_table(self.data_dir / "threshold_sensitivity.csv", build, stamp, force=self.force)


__all__ = ["HeatwaveCase", "REFERENCE_LABELS", "METRIC_UNITS", "PRIMARY_METRICS", "TOLERANT_METRICS",
           "START_METRICS", "CONTRAST_METRICS", "SERIES_YLABELS"]
