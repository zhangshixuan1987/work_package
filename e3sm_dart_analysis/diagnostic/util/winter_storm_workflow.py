"""Workflow of the January 2012 Pacific Northwest storm notebooks (``hind_lead_eval/07a-07d``).

``WinterStormCase`` holds one set of ``WinterStormSettings`` and gives every notebook the same products,
each computed once and cached on disk with the settings it depends on:

* ``obs``: ERA5 TREFHT, GPCP precipitation, both TREFHT climatologies and the region masks on the map box
  (cached, ``observed_fields.nc``), and from them the regional observations, the observed peak, the
  observed metrics and the coverage;
* ``field(config, variable)``: daily member fields on the map or environment box
  (cached, ``box_{map,env}_<CONFIG>_DJF_<VAR>.nc``), ``pooled``: the members of all starts on calendar days;
* ``series``: regional series of every configuration, region and bias treatment (in memory);
* ``metrics``: storm metrics of every member, region and bias treatment (cached, ``storm_metrics.csv``).

A product is read when its file exists with the same settings, unless ``force_compute`` (this case's
results) or ``force_inputs`` (also the model and observation reads) is set. Methods:
``util/winter_storm_case.py`` and ``util/heatwave_case.py``.
"""
from __future__ import annotations

from functools import cached_property
from types import SimpleNamespace

import numpy as np
import pandas as pd
import xarray as xr

from configs.fcst_init_experiment_config import balanced_runs, case_stamp, check_configs, init_dates, readiness
from configs.fcst_init_experiment_config import members as config_members
from configs.hind_lead_eval_cases import WinterStormSettings
from configs.output_paths import workflow_output_dirs
from util.fcst_init_data import region_masks
from util.heatwave_case import (calendar, climatologies, coverage_table, model_box, observed_box, region_mean,
                                region_weights)
from util.hind_case_io import cached_dataset, cached_table
from util.winter_storm_case import KELVIN, pattern_scores, storm_metrics

SUBDIR = "winter_storm_2012"
MAP_VARS, ENV_VARS = ("TREFHT", "PRECT", "H2OSNO"), ("PSL", "Z500", "TUQ", "TVQ")
SERIES_VARS = ("TREFHT", "PRECT", "H2OSNO")
SERIES_YLABELS = {"TREFHT": "daily-mean T2m (C)", "PRECT": "precipitation (mm/day)", "H2OSNO": "snow water (mm, model)"}
METRICS = ["event_precip_mm", "max_daily_precip_mm", "peak_precip_date_error_days", "event_mean_t_c", "event_min_t_c",
           "precip_freezing_days_mm", "precip_near_freezing_days_mm", "precip_pattern_r", "t_pattern_r"]
OBS_PATTERN = dict(precip_pattern_r=1.0, precip_pattern_rmse_mm=0.0, t_pattern_r=1.0, t_pattern_rmse_c=0.0)

_OBS_KEYS = ("map_box", "regions", "obs_period", "clim_years", "amip_cases")
_MODEL_KEYS = ("configs", "season", "n_days", "max_members", "map_box", "env_box")


class WinterStormCase:
    """Products of the winter-storm case for one set of settings (see the module docstring)."""

    def __init__(self, settings: WinterStormSettings, *, force_compute: bool = False, force_inputs: bool = False):
        self.cfg = settings
        self.force_inputs = force_inputs
        self.force = force_compute or force_inputs
        self._fields, self._pooled = {}, {}
        check_configs(list(settings.configs))
        p = settings.paths
        self.data_dir, self.figure_dir = workflow_output_dirs(p.workflow_name, SUBDIR, output_root=p.diagnostic_output_root)
        self.inits = list(init_dates([settings.season])[settings.season])
        self.members = {c: tuple(config_members(c, settings.season))[:settings.max_members] for c in settings.configs}
        self.region_names = [r["name"] for r in settings.regions]
        print(f"[DATA] {self.data_dir}\n[FIGURE] {self.figure_dir}\n[IC DEFECT] {list(settings.ic_defect)}")

    def stamp(self, *keys: str, extra: str = "") -> str:
        return self.cfg.stamp(*keys) + extra + case_stamp(self.cfg.configs)

    @property
    def metrics_stamp(self) -> str:
        return self.stamp(*sorted(set(_OBS_KEYS + _MODEL_KEYS + ("event_window", "peak_region", "peak_tolerance_days",
                                                                 "biases"))), extra=str(self.members))

    # ------------------------------------------------------------------
    # Inputs
    # ------------------------------------------------------------------
    @cached_property
    def runs(self) -> list:
        """(config, season) pairs whose members are post-processed (raises when one is missing)."""
        status = readiness(list(self.cfg.configs), [self.cfg.season],
                           streams=(("atm", "daily"), ("atm", "6hourly"), ("lnd", "daily")),
                           data_root=self.cfg.paths.input_data_root)
        out = balanced_runs(status, self.cfg.require_all_members)
        print(f"[RUNS] {out}\n[STARTS] {[str(i) for i in self.inits]}\n"
              f"[MEMBERS] {({c: len(m) for c, m in self.members.items()})}")
        return out

    def _observed_fields(self) -> xr.Dataset:
        cfg, p = self.cfg, self.cfg.paths
        years = range(cfg.clim_years[0], cfg.clim_years[1] + 1)
        obs_T = observed_box("TREFHT", *cfg.obs_period, cfg.map_box, observation_root=p.observation_root)
        obs_P = observed_box("PRECT", *cfg.obs_period, cfg.map_box, observation_root=p.observation_root)
        clim = climatologies("TREFHT", *cfg.obs_period, years, cfg.map_box, observation_root=p.observation_root,
                             amip_root=p.amip_root, amip_cases=cfg.amip_cases, land_support_dir=p.land_support_dir)
        masks = region_masks(list(cfg.regions), obs_T.lat.values, obs_T.lon.values, land_support_dir=p.land_support_dir)
        ds = xr.merge([obs_T.rename("obs_T"), obs_P.rename("obs_P"), masks.rename("masks").astype("int8"),
                       clim.clim_obs, clim.clim_amip], join="exact", combine_attrs="drop")
        ds.attrs.update(clim_obs=clim.attrs["clim_obs"], clim_amip=clim.attrs["clim_amip"])
        return ds

    @cached_property
    def obs(self) -> SimpleNamespace:
        """Observed fields (cached) and everything derived from them: regional weights and series, observed
        precipitation peak, event-mean maps, observed metrics, simulation coverage."""
        cfg = self.cfg
        ds = cached_dataset(self.data_dir / "observed_fields.nc", self._observed_fields, self.stamp(*_OBS_KEYS),
                            force=self.force_inputs)
        masks = ds.masks.astype(bool)
        o = SimpleNamespace(T=ds.obs_T - KELVIN, P=ds.obs_P, clim=ds[["clim_obs", "clim_amip"]], times=ds.time.values)
        o.bias_map = ds.clim_amip - ds.clim_obs                        # model mean bias of each day (K)
        o.W = {r: region_weights(masks.sel(region=r, drop=True)) for r in self.region_names}
        o.W_ALL = region_weights(masks.any("region"))
        print({k: int((v > 0).sum()) for k, v in o.W.items()}, "land points")
        o.reg = {r: dict(TREFHT=region_mean(o.T, w), PRECT=region_mean(o.P, w), H2OSNO=None) for r, w in o.W.items()}
        pk = o.reg[cfg.peak_region]["PRECT"].sel(time=slice(*cfg.event_window))
        o.peak = str(pd.Timestamp(pk.time.values[int(np.argmax(pk.values))]).date())
        o.event_T = o.T.sel(time=slice(*cfg.event_window)).mean("time")
        o.event_P = o.P.sel(time=slice(*cfg.event_window)).sum("time")
        o.metrics = pd.DataFrame([dict(region=r, **storm_metrics(o.reg[r]["TREFHT"].values, o.reg[r]["PRECT"].values,
                                                                 o.times, cfg.event_window, o.peak,
                                                                 tolerance_days=cfg.peak_tolerance_days))
                                  for r in o.W]).set_index("region")
        o.coverage = coverage_table(self.inits, cfg.n_days, cfg.event_window)
        print(f"[EVENT] {cfg.event_window[0]} to {cfg.event_window[1]}; observed {cfg.peak_region} "
              f"precipitation peak {o.peak}")
        return o

    def field(self, config: str, variable: str) -> xr.DataArray:
        """(init, member, lead, lat, lon) daily field of ``config``: map variables on the map box,
        environment variables on the environment box (cached)."""
        key = (config, variable)
        if key not in self._fields:
            cfg = self.cfg
            tag, box = ("map", cfg.map_box) if variable in MAP_VARS else ("env", cfg.env_box)
            cache = self.data_dir / f"box_{tag}_{config}_{cfg.season}_{variable}.nc"
            if self.force_inputs or not cache.exists():
                self.runs
            da = model_box(config, cfg.season, variable, self.members[config], box, n_days=cfg.n_days, cache=cache,
                           force=self.force_inputs, workers=cfg.workers, data_root=cfg.paths.input_data_root,
                           land_support_dir=cfg.paths.land_support_dir)
            if tag == "map" and not (np.allclose(da.lat, self.obs.T.lat) and np.allclose(da.lon, self.obs.T.lon)):
                raise ValueError(f"model ({config} {variable}) and observation grids differ on the map box")
            print(f"[MODEL] {config} {variable}: {dict(da.sizes)}")
            self._fields[key] = da
        return self._fields[key]

    def pooled(self, config: str, variable: str) -> xr.DataArray:
        """(member, time, lat, lon) on the calendar days of the observed period, the members of all starts
        stacked (labels made unique per start: ``MM-DD_ENxx``); days a start does not cover are NaN."""
        key = (config, variable)
        if key not in self._pooled:
            da = self.field(config, variable)
            days = pd.date_range(*self.cfg.obs_period, freq="D").values
            out = []
            for init in self.inits:
                x = da.sel(init=np.datetime64(init, "ns")).drop_vars("init", errors="ignore")
                x = x.assign_coords(lead=calendar(init, x.lead.values).astype("datetime64[ns]")).rename(lead="time")
                x = x.reindex(time=days).assign_coords(member=[f"{str(init)[5:10]}_{m}" for m in x.member.values])
                out.append(x)
            self._pooled[key] = xr.concat(out, "member")
        return self._pooled[key]

    def temperature(self, config: str, bias: str) -> xr.DataArray:
        """Pooled daily-mean TREFHT (C), minus the day's AMIP-minus-observed climatology with ``amip``."""
        T = self.pooled(config, "TREFHT") - KELVIN
        return T - self.obs.bias_map if bias == "amip" else T

    # ------------------------------------------------------------------
    # Derived products
    # ------------------------------------------------------------------
    @cached_property
    def series(self) -> dict:
        """``series[bias, region, config]``: regional TREFHT (C), PRECT, H2OSNO of every pooled member."""
        out = {}
        for b in self.cfg.biases:
            for c in self.cfg.configs:
                T, P, S = self.temperature(c, b), self.pooled(c, "PRECT"), self.pooled(c, "H2OSNO")
                for r, w in self.obs.W.items():
                    out[b, r, c] = dict(TREFHT=region_mean(T, w), PRECT=region_mean(P, w), H2OSNO=region_mean(S, w))
        return out

    def series_dataset(self) -> xr.Dataset:
        """The regional series, dimensions (bias, region, config, member_index, time); labels in ``member``."""
        cfg, configs = self.cfg, list(self.cfg.configs)

        def stack(v):
            return xr.concat([xr.concat([xr.concat([self.series[b, r, c][v].drop_vars("member")
                                                    .rename(member="member_index") for c in configs],
                                                   pd.Index(configs, name="config"), join="outer")
                                         for r in self.region_names], pd.Index(self.region_names, name="region"))
                              for b in cfg.biases], pd.Index(list(cfg.biases), name="bias"))

        ds = xr.Dataset({v: stack(v) for v in SERIES_VARS})
        n = ds.sizes["member_index"]
        ds["member"] = (("config", "member_index"),
                        np.array([[str(m) for m in self.pooled(c, "TREFHT").member.values][:n] for c in configs]))
        for v in ("TREFHT", "PRECT"):
            ds[f"obs_{v}"] = xr.concat([self.obs.reg[r][v] for r in self.region_names],
                                       pd.Index(self.region_names, name="region"))
        ds.attrs.update(units_TREFHT="C", units_PRECT="mm/day", units_H2OSNO="mm")
        return ds

    def _metrics(self) -> pd.DataFrame:
        cfg, o, rows = self.cfg, self.obs, []
        for b in cfg.biases:
            for c in cfg.configs:
                T = self.temperature(c, b)
                start_of = np.array([m.split("_")[0] for m in T.member.values])
                ev_T = T.sel(time=slice(*cfg.event_window)).mean("time")
                ev_P = self.pooled(c, "PRECT").sel(time=slice(*cfg.event_window)).sum("time")
                for k, m in enumerate(T.member.values):
                    spat = dict(zip(("precip_pattern_r", "precip_pattern_rmse_mm"),
                                    pattern_scores(ev_P.isel(member=k), o.event_P, o.W_ALL)))
                    spat.update(zip(("t_pattern_r", "t_pattern_rmse_c"),
                                    pattern_scores(ev_T.isel(member=k), o.event_T, o.W_ALL)))
                    for r in o.W:
                        s = self.series[b, r, c]
                        met = storm_metrics(s["TREFHT"].isel(member=k).values, s["PRECT"].isel(member=k).values,
                                            o.times, cfg.event_window, o.peak, tolerance_days=cfg.peak_tolerance_days)
                        rows.append(dict(bias=b, config=c, init=start_of[k], member=str(m), region=r, **met, **spat))
        return pd.DataFrame(rows)

    @cached_property
    def metrics(self) -> pd.DataFrame:
        """Storm metrics of every pooled member, region, configuration and bias treatment (cached)."""
        return cached_table(self.data_dir / "storm_metrics.csv", self._metrics, self.metrics_stamp, force=self.force)

    def event_maps(self) -> xr.Dataset:
        """Event-total precipitation and event-mean T2m anomaly (vs each climatology): observations and the
        ensemble-and-start mean of every configuration (cached, ``storm_event_maps.nc``)."""
        cfg, o = self.cfg, self.obs

        def build():
            win = slice(*cfg.event_window)
            configs = pd.Index(list(cfg.configs), name="config")
            ds = xr.Dataset({"obs_precip": o.event_P,
                             "obs_t2m_anom": (o.T + KELVIN - o.clim.clim_obs).sel(time=win).mean("time"),
                             "precip": xr.concat([self.pooled(c, "PRECT").sel(time=win).sum("time").mean("member")
                                                  for c in cfg.configs], configs)})
            for b in cfg.biases:
                ref = o.clim.clim_amip if b == "amip" else o.clim.clim_obs
                ds[f"t2m_anom_{b}"] = xr.concat([(self.pooled(c, "TREFHT") - ref).sel(time=win).mean(("time", "member"))
                                                 for c in cfg.configs], configs)
            ds.attrs.update(units_precip="event total (mm)", units_t2m="anomaly (K)")
            return ds
        return cached_dataset(self.data_dir / "storm_event_maps.nc", build, self.metrics_stamp, force=self.force)

    def environment(self) -> xr.Dataset:
        """PSL, Z500, TUQ, TVQ on the environment days: ERA5 (``source = OBS``) and the ensemble-and-start
        mean of every configuration (cached, ``storm_environment.nc``). Z500 is the 00Z state ending each day."""
        cfg = self.cfg
        days = pd.to_datetime(list(cfg.env_days))

        def build():
            root = cfg.paths.observation_root
            obs = {}
            for v in ENV_VARS:
                if v == "Z500":   # 00Z state ending each day: read the next day's 00Z and label it with the day
                    x = observed_box(v, str((days.min() + pd.Timedelta(days=1)).date()),
                                     str((days.max() + pd.Timedelta(days=1)).date()), cfg.env_box, observation_root=root)
                    x = x.assign_coords(time=x.time.values - np.timedelta64(1, "D"))
                else:
                    x = observed_box(v, str(days.min().date()), str(days.max().date()), cfg.env_box, observation_root=root)
                obs[v] = x.sel(time=days.values)
            sources = [xr.Dataset(obs)]
            for c in cfg.configs:
                sources.append(xr.Dataset({v: self.pooled(c, v).sel(time=days.values).mean("member") for v in ENV_VARS}))
            return xr.concat(sources, pd.Index(["OBS", *cfg.configs], name="source"), join="override",
                             combine_attrs="drop")
        stamp = self.stamp(*sorted(set(_MODEL_KEYS + ("env_days",))), extra=str(self.members))
        return cached_dataset(self.data_dir / "storm_environment.nc", build, stamp, force=self.force)


__all__ = ["WinterStormCase", "MAP_VARS", "ENV_VARS", "SERIES_VARS", "SERIES_YLABELS", "METRICS", "OBS_PATTERN"]
