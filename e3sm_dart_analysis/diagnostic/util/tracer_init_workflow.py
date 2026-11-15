"""Workflow of the tracer-initialization notebook (``hind_lead_eval/08``).

``TracerInitCase`` holds one set of ``TracerInitSettings`` and gives the notebook its products, each cached on
disk with the settings it depends on:

* ``summary(config, season, variable)``: every member reduced as it is read (``util.tracer_init.reduce_start``):
  regional and zonal means of each member, ensemble mean and variance of the lead-window means, and the regional
  ensemble-mean error against the reference where one exists (``summary_<CONFIG>_<SEASON>_<VAR>.nc``);
* ``reference(season, variable)``: regional and zonal means of the reference (``reference_<SEASON>_<VAR>.nc``);
* ``contrasts()``: regional difference, standard error, member spread and z of every contrast by lead day
  (``tracer_init_contrasts.csv``); ``recovery()``: day-1 and peak difference and the recovery days
  (``tracer_init_recovery.csv``); ``skill()``: regional bias, RMSE and spread against the reference
  (``tracer_init_skill.csv``);
* ``zonal(season)``: zonal-mean contrasts, latitude x lead (``tracer_init_zonal_<SEASON>.nc``); ``maps(season)``:
  lead-window contrast maps (``tracer_init_maps_<SEASON>.nc``).

A product is read when its file exists with the same settings, unless ``force_compute`` (this case's results) or
``force_inputs`` (also the model and observation reads) is set. Methods: ``util/tracer_init.py``.
"""
from __future__ import annotations

from functools import cached_property

import numpy as np
import pandas as pd
import xarray as xr

from configs.fcst_init_experiment_config import (balanced_runs, case_names, check_configs, experiments_table,
                                                  init_dates, readiness)
from configs.fcst_init_experiment_config import members as config_members
from configs.hind_lead_eval_cases import TracerInitSettings
from configs.output_paths import workflow_output_dirs
from util.fcst_init_data import open_reference, region_label, region_masks, select_regions, spec
from util.heatwave_case import _load_member_box, _parallel_map
from util.hind_case_io import cached_dataset, cached_table
from util.tracer_init import (GLOBE, contrast_stats, contrast_stats_windows, pooled_skill, recovery, reduce_reference,
                              reduce_start)

SUBDIR = "tracer_init"

_READ_KEYS = ("regions", "custom_regions", "max_lead", "lead_windows")
_CONTRAST_KEYS = ("contrasts",)


class TracerInitCase:
    """Products of the tracer-initialization analysis for one set of settings (see the module docstring)."""

    def __init__(self, settings: TracerInitSettings, *, force_compute: bool = False, force_inputs: bool = False):
        self.cfg = settings
        self.force_inputs = force_inputs
        self.force = force_compute or force_inputs
        configs = list(settings.configs)
        check_configs(configs)
        for name, (a, b) in settings.contrasts.items():
            if a not in configs or b not in configs:
                raise ValueError(f"contrast {name!r} uses {a!r}/{b!r}; both must be in configs {configs}")
        if settings.test not in configs or settings.control not in configs:
            raise ValueError(f"test {settings.test!r} and control {settings.control!r} must be in configs {configs}")
        for v in settings.variables:
            spec(v)
        p = settings.paths
        self.data_dir, self.figure_dir = workflow_output_dirs(p.workflow_name, SUBDIR,
                                                              output_root=p.diagnostic_output_root)
        self.regions = select_regions(settings.regions, settings.custom_regions)
        self.region_labels = {r["name"]: region_label(r) for r in self.regions}
        self.members = {(c, s): tuple(config_members(c, s))[:settings.max_members]
                        for c in configs for s in settings.seasons}
        self._masks = None
        print(f"[DATA] {self.data_dir}\n[FIGURE] {self.figure_dir}")
        print(f"[CONTRASTS] " + "; ".join(f"{n}: {a} - {b}" for n, (a, b) in settings.contrasts.items()))

    def experiments(self) -> pd.DataFrame:
        return experiments_table(list(self.cfg.configs))

    @property
    def units(self) -> dict[str, str]:
        return {v: spec(v).units for v in self.cfg.variables}

    def verifiable(self, variable: str) -> bool:
        return bool(spec(variable).reference)

    def _stamp(self, *keys: str, configs=None) -> str:
        configs = self.cfg.configs if configs is None else configs
        members = {f"{c} {s}": m for (c, s), m in self.members.items() if c in configs}
        return self.cfg.stamp(*keys) + str(members) + str({c: case_names()[c] for c in configs})

    # ------------------------------------------------------------------
    # Inputs
    # ------------------------------------------------------------------
    @cached_property
    def seasons(self) -> list[str]:
        """Seasons of the settings, every configuration post-processed in each (raises otherwise). Without a
        forced re-read, existing summaries are trusted and the model output is not scanned."""
        cfg = self.cfg
        cached = all(self._summary_path(c, s, v).exists() for c in cfg.configs for s in cfg.seasons
                     for v in cfg.variables)
        if not cached or self.force_inputs:
            status = readiness(list(cfg.configs), cfg.seasons, streams=(("atm", "daily"), ("atm", "6hourly")),
                               data_root=cfg.paths.input_data_root)
            balanced_runs(status, cfg.require_all_members)
        print(f"[SEASONS] {list(cfg.seasons)}; members: "
              + ", ".join(f"{c} {s} {len(m)}" for (c, s), m in self.members.items()))
        return list(cfg.seasons)

    def masks(self, lat, lon) -> xr.DataArray:
        if self._masks is None:
            self._masks = region_masks(self.regions, np.asarray(lat), np.asarray(lon),
                                       land_support_dir=self.cfg.paths.land_support_dir)
        return self._masks

    def _summary_path(self, config, season, variable):
        return self.data_dir / f"summary_{config}_{season}_{variable}.nc"

    def _reference_field(self, season: str, variable: str) -> xr.DataArray:
        inits = init_dates([season])[season]
        return open_reference(variable, inits, max_lead=self.cfg.max_lead,
                              observation_root=self.cfg.paths.observation_root).drop_vars("time").load()

    def summary(self, config: str, season: str, variable: str) -> xr.Dataset:
        """(init, member, ...) reductions of every member of ``config`` in ``season`` (cached)."""
        cfg = self.cfg

        def build():
            inits = list(init_dates([season])[season])
            members = list(self.members[config, season])
            names = {config: case_names()[config]}
            ref = self._reference_field(season, variable) if self.verifiable(variable) else None
            parts = []
            for i in inits:
                jobs = [(config, i, m, variable, cfg.max_lead, GLOBE, cfg.paths.input_data_root,
                         cfg.paths.land_support_dir, names) for m in members]
                fields = _parallel_map(_load_member_box, jobs, cfg.workers)
                ens = xr.concat(fields, pd.Index(members, name="member"), coords="minimal", compat="override")
                masks = self.masks(ens.lat.values, ens.lon.values)
                r = None
                if ref is not None:
                    r = ref.sel(init=np.datetime64(i, "ns"))
                    if not (np.allclose(r.lat, ens.lat) and np.allclose(r.lon, ens.lon)):
                        raise ValueError(f"{variable}: model and reference grids differ")
                parts.append(reduce_start(ens, masks, cfg.lead_windows, r))
                print(f"[READ] {config} {season} {variable} {i}: {len(members)} members", flush=True)
            out = xr.concat(parts, pd.Index(np.array(inits, dtype="datetime64[ns]"), name="init"))
            out["win_mean"] = out.win_mean.astype("float32")
            out["win_var"] = out.win_var.astype("float32")
            return out.assign_attrs(config=config, season=season, variable=variable, units=spec(variable).units,
                                    n_members=len(members))
        return cached_dataset(self._summary_path(config, season, variable), build,
                              self._stamp(*_READ_KEYS, configs=[config]), force=self.force_inputs)

    def reference(self, season: str, variable: str) -> xr.Dataset | None:
        """Regional and zonal means of the reference of ``variable`` (cached); None without a reference."""
        if not self.verifiable(variable):
            return None

        def build():
            ref = self._reference_field(season, variable)
            out = reduce_reference(ref, self.masks(ref.lat.values, ref.lon.values))
            return out.assign_attrs(reference=spec(variable).reference, variable=variable,
                                    units=spec(variable).units)
        return cached_dataset(self.data_dir / f"reference_{season}_{variable}.nc", build,
                              self.cfg.stamp(*_READ_KEYS), force=self.force_inputs)

    def series(self, variables, region: str):
        """Regional series of every configuration ({(season, config, variable): (init, member, lead)}) and of
        the references ({(season, variable): (init, lead)})."""
        model, obs = {}, {}
        for s in self.seasons:
            for v in variables:
                for c in self.cfg.configs:
                    model[s, c, v] = self.summary(c, s, v).reg.sel(region=region)
                ref = self.reference(s, v)
                if ref is not None:
                    obs[s, v] = ref.reg.sel(region=region)
        return model, obs

    # ------------------------------------------------------------------
    # Results (cached)
    # ------------------------------------------------------------------
    def contrasts(self) -> pd.DataFrame:
        """Long table (contrast, season, variable, region, lead): delta, se, spread, z."""
        cfg = self.cfg

        def build():
            rows = []
            for s in self.seasons:
                for v in cfg.variables:
                    for name, (a, b) in cfg.contrasts.items():
                        st = contrast_stats(self.summary(a, s, v).reg, self.summary(b, s, v).reg)
                        rows.append(st.to_dataframe().reset_index().assign(contrast=name, season=s, variable=v))
            cols = ["contrast", "season", "variable", "region", "lead", "delta", "se", "spread", "z"]
            return pd.concat(rows, ignore_index=True)[cols]
        return cached_table(self.data_dir / "tracer_init_contrasts.csv", build,
                            self._stamp(*_READ_KEYS, *_CONTRAST_KEYS, "variables", "seasons"), force=self.force)

    def recovery(self) -> pd.DataFrame:
        """(contrast, season, variable, region): day-1 and peak difference, the member spread at the peak, the
        last lead days outside k x spread and significant (0: never; max_lead: not recovered)."""
        cfg = self.cfg

        def build():
            t = self.contrasts()
            rows = []
            for key, g in t.groupby(["contrast", "season", "variable", "region"], sort=False):
                st = xr.Dataset.from_dataframe(g.set_index("lead")[["delta", "se", "spread", "z"]])
                rows.append(dict(zip(("contrast", "season", "variable", "region"), key),
                                 **recovery(st, cfg.recovery_k, cfg.z_crit)))
            out = pd.DataFrame(rows)
            out["units"] = out.variable.map(self.units)
            return out
        return cached_table(self.data_dir / "tracer_init_recovery.csv", build,
                            self._stamp(*_READ_KEYS, *_CONTRAST_KEYS, "variables", "seasons", "recovery_k", "z_crit"),
                            force=self.force)

    def skill(self) -> pd.DataFrame:
        """Long table (config, season, variable, region, lead): bias, rmse, spread, spread_ratio of the
        variables with a reference."""
        cfg = self.cfg

        def build():
            rows = []
            for s in self.seasons:
                for v in cfg.variables:
                    if not self.verifiable(v):
                        continue
                    for c in cfg.configs:
                        sk = pooled_skill(self.summary(c, s, v))
                        rows.append(sk.to_dataframe().reset_index().assign(config=c, season=s, variable=v))
            cols = ["config", "season", "variable", "region", "lead", "bias", "rmse", "spread", "spread_ratio"]
            return pd.concat(rows, ignore_index=True)[cols]
        return cached_table(self.data_dir / "tracer_init_skill.csv", build,
                            self._stamp(*_READ_KEYS, "variables", "seasons"), force=self.force)

    def skill_datasets(self, table: pd.DataFrame) -> dict:
        """{(season, config, variable): (region, lead) Dataset} of a ``skill`` table, as ``plot_lead_metric`` takes."""
        out = {}
        for (c, s, v), g in table.groupby(["config", "season", "variable"], sort=False):
            ds = xr.Dataset.from_dataframe(g.set_index(["region", "lead"])[["bias", "rmse", "spread", "spread_ratio"]])
            out[s, c, v] = ds.assign_attrs(units=self.units[v])
        return out

    def zonal(self, season: str) -> xr.Dataset:
        """(contrast, variable, lead, lat) zonal-mean contrast statistics (cached)."""
        cfg = self.cfg

        def build():
            per_c = []
            for name, (a, b) in cfg.contrasts.items():
                per_v = [contrast_stats(self.summary(a, season, v).zonal, self.summary(b, season, v).zonal)
                         for v in cfg.variables]
                per_c.append(xr.concat(per_v, pd.Index(list(cfg.variables), name="variable")))
            out = xr.concat(per_c, pd.Index(list(cfg.contrasts), name="contrast"))
            return out.astype("float32").transpose("contrast", "variable", "lead", "lat")
        return cached_dataset(self.data_dir / f"tracer_init_zonal_{season}.nc", build,
                              self._stamp(*_READ_KEYS, *_CONTRAST_KEYS, "variables"), force=self.force)

    def maps(self, season: str, variables) -> xr.Dataset:
        """(contrast, variable, window, lat, lon) lead-window contrast maps of ``variables`` (cached)."""
        cfg = self.cfg
        variables = list(variables)

        def build():
            per_c = []
            for name, (a, b) in cfg.contrasts.items():
                per_v = [contrast_stats_windows(self.summary(a, season, v), self.summary(b, season, v))
                         for v in variables]
                per_c.append(xr.concat(per_v, pd.Index(variables, name="variable")))
            out = xr.concat(per_c, pd.Index(list(cfg.contrasts), name="contrast"))
            return out.astype("float32").transpose("contrast", "variable", "window", "lat", "lon")
        return cached_dataset(self.data_dir / f"tracer_init_maps_{season}.nc", build,
                              self._stamp(*_READ_KEYS, *_CONTRAST_KEYS) + str(variables), force=self.force)


__all__ = ["TracerInitCase", "SUBDIR"]
