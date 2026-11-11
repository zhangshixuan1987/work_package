"""Workflow of the initialization-benefit notebooks (``hind_lead_eval/05a-05c``).

``InitBenefitCase`` holds one set of ``InitBenefitSettings`` and gives every notebook the same products:

* ``daily`` / ``refs``: daily member values of every role and the observations and both climatologies on
  the region's box, per usable season (cached, ``daily_<CONFIG>_<SEASON>_<VAR>.nc``,
  ``daily_ref_<SEASON>_<VAR>.nc``);
* ``curves(scale)``: skill and contrast curves with paired-bootstrap intervals on weekly means or daily
  running means, panel (b) per bias treatment and season set and panel (c) (cached,
  ``init_benefit_curves_<scale>.csv``); ``panel(...)`` turns them back into the plotting dictionaries;
* ``maps(bias, set)``: gridpoint-ACC contrasts on weekly means (cached, ``init_benefit_maps_<BIAS>_<SET>.nc``).

A product is read when its file exists with the same settings, unless ``force_compute`` (this case's
results) or ``force_inputs`` (also the model and observation reads) is set. Methods: ``util/init_benefit.py``.
"""
from __future__ import annotations

from functools import cached_property

import numpy as np
import pandas as pd
import xarray as xr

from configs.fcst_init_experiment_config import CONFIGS, case_stamp, check_configs, readiness
from configs.hind_lead_eval_cases import InitBenefitSettings
from configs.output_paths import workflow_output_dirs
from util.fcst_init_data import region_masks
from util.hind_case_io import cached_dataset, cached_table
from util.init_benefit import (contrast_values, daily_members, daily_references, era5_benefits, gridpoint_acc_np,
                               lead_in_weeks, member_set, model_member_benefits, region_weights, results_table,
                               role_configs, running_leads, to_weeks)

SUBDIR = "init_benefit"
SCALES = ("weekly", "daily")

_DATA_KEYS = ("roles", "max_members", "variable", "region", "max_lead", "clim_years", "amip_cases")
_SCORE_KEYS = _DATA_KEYS + ("contrasts", "biases", "model_truth_role", "n_boot", "seed", "smooth_days")


class InitBenefitCase:
    """Products of the initialization-benefit analysis for one set of settings (see the module docstring)."""

    def __init__(self, settings: InitBenefitSettings, *, force_compute: bool = False, force_inputs: bool = False):
        self.cfg = settings
        self.force_inputs = force_inputs
        self.force = force_compute or force_inputs
        self.roles = role_configs(settings.roles)
        check_configs(list(self.roles.values()))
        p = settings.paths
        self.data_dir, self.figure_dir = workflow_output_dirs(p.workflow_name, SUBDIR, output_root=p.diagnostic_output_root)
        # {season: {role: member directories}}; the DA labels differ by season
        self.members = {s: {r: tuple(m)[:settings.max_members] for r, m in member_set(settings.roles, s).items()}
                        for s in settings.seasons}
        self.n_members = len(next(iter(self.members[settings.seasons[0]].values())))
        print(f"[DATA] {self.data_dir}\n[FIGURE] {self.figure_dir}")
        for s, by_role in self.members.items():
            for r, m in by_role.items():
                print(f"[MEMBERS] {s} {r}: {len(m)}: {', '.join(m)}")

    def role_table(self) -> pd.DataFrame:
        return pd.DataFrame([dict(role=r, config=c, atmosphere=CONFIGS[c][0], land=CONFIGS[c][1])
                             for r, c in self.roles.items()])

    # ------------------------------------------------------------------
    # Inputs
    # ------------------------------------------------------------------
    def _cache(self, config: str, season: str) -> object:
        return self.data_dir / f"daily_{config}_{season}_{self.cfg.variable}.nc"

    @cached_property
    def seasons(self) -> list[str]:
        """Seasons in which every role's members are post-processed. When every daily cache already exists
        (and no re-read is forced), the cached seasons are used without scanning the model output."""
        cfg = self.cfg
        if not self.force_inputs and all(self._cache(c, s).exists() for s in cfg.seasons for c in self.roles.values()):
            print(f"[SEASONS] from the daily caches: {list(cfg.seasons)}")
            return list(cfg.seasons)
        status = readiness(list(self.roles.values()), cfg.seasons, streams=(("atm", "6hourly"),),
                           data_root=cfg.paths.input_data_root)
        summary = status.groupby(["config", "season"]).agg(ready=("ready", "sum"), members=("members", "sum"))
        print(summary.unstack("season"))
        ready = {(c, s) for (c, s), row in summary.iterrows()
                 if row.ready >= (row.members if cfg.require_all_members else 1)}
        seasons = [s for s in cfg.seasons if all((c, s) in ready for c in self.roles.values())]
        print(f"[SEASONS] usable (every role ready): {seasons or 'none'}")
        if not seasons:
            raise RuntimeError("No season with every role post-processed")
        return seasons

    @cached_property
    def refs(self) -> dict[str, xr.Dataset]:
        """Daily observations and both climatologies per season, (init, lead, lat, lon) on the box (cached)."""
        cfg, p = self.cfg, self.cfg.paths
        years = list(range(cfg.clim_years[0], cfg.clim_years[1] + 1))
        out = {s: daily_references(cfg.variable, s, years, cfg.region, max_lead=cfg.max_lead,
                                   observation_root=p.observation_root, amip_root=p.amip_root,
                                   amip_cases=cfg.amip_cases, land_support_dir=p.land_support_dir,
                                   cache=self.data_dir / f"daily_ref_{s}_{cfg.variable}.nc", force=self.force_inputs)
               for s in self.seasons}
        d0 = out[self.seasons[0]]
        print(f"[CLIM] observed: {d0.attrs['clim_obs']}\n[CLIM] model: {d0.attrs['clim_amip']}")
        return out

    @cached_property
    def daily(self) -> dict[tuple[str, str], xr.DataArray]:
        """(init, member, lead, lat, lon) daily values of every role's members per season (cached)."""
        cfg, out = self.cfg, {}
        for s in self.seasons:
            for r, c in self.roles.items():
                out[r, s] = daily_members(c, s, cfg.variable, self.members[s][r], cfg.region, max_lead=cfg.max_lead,
                                          cache=self._cache(c, s), force=self.force_inputs,
                                          data_root=cfg.paths.input_data_root,
                                          land_support_dir=cfg.paths.land_support_dir)
                print(f"[DAILY] {r} {c} {s}: {dict(out[r, s].sizes)}")
        return out

    # ------------------------------------------------------------------
    # Grid, axes and anomalies
    # ------------------------------------------------------------------
    @cached_property
    def grid(self) -> xr.Dataset:
        """Lat/lon of the region's box and the region mask (from the first season's references)."""
        d0 = self.refs[self.seasons[0]]
        mask = region_masks([self.cfg.region], d0.lat.values, d0.lon.values,
                            land_support_dir=self.cfg.paths.land_support_dir).isel(region=0, drop=True)
        return xr.Dataset({"mask": mask})

    @property
    def weights(self) -> np.ndarray:
        return region_weights(self.grid.mask)

    @property
    def leads(self) -> np.ndarray:
        return self.refs[self.seasons[0]].lead.values

    @cached_property
    def weeks(self) -> np.ndarray:
        return to_weeks(self.refs[self.seasons[0]].obs).week.values

    @property
    def x_weeks(self) -> np.ndarray:
        return lead_in_weeks(self.leads)

    @property
    def sets(self) -> dict[str, list[str]]:
        """Season sets: each season alone (when there are several) and all together (``ALL``)."""
        return {**({s: [s] for s in self.seasons} if len(self.seasons) > 1 else {}), "ALL": list(self.seasons)}

    def axis(self, scale: str) -> tuple[str, np.ndarray]:
        return ("week", self.weeks) if scale == "weekly" else ("lead", self.leads)

    def _prepare(self, da: xr.DataArray, scale: str) -> np.ndarray:
        """Weekly means or ``smooth_days`` running means of an anomaly field, dimensions in a fixed order."""
        out = to_weeks(da) if scale == "weekly" else running_leads(da, self.cfg.smooth_days)
        dims = [d for d in ("init", "member", self.axis(scale)[0], "lat", "lon") if d in out.dims]
        return out.transpose(*dims).values

    def stack(self, set_seasons, scale: str, forecast_clim: str):
        """Forecast anomalies per role, observed anomalies, and the season of each start."""
        fa = {r: np.concatenate([self._prepare(self.daily[r, s] - self.refs[s][forecast_clim], scale)
                                 for s in set_seasons]) for r in self.roles}
        oa = np.concatenate([self._prepare(self.refs[s].obs - self.refs[s].clim_obs, scale) for s in set_seasons])
        labels = [s for s in set_seasons for _ in range(self.refs[s].sizes["init"])]
        return fa, oa, labels

    def _stamp(self, *keys: str) -> str:
        return self.cfg.stamp(*keys) + str(self.members) + str(self.seasons) + case_stamp(self.roles.values())

    # ------------------------------------------------------------------
    # Scores (cached)
    # ------------------------------------------------------------------
    def curves(self, scale: str) -> pd.DataFrame:
        """Skill and contrast curves with 2.5-97.5 % paired-bootstrap intervals on the ``scale`` time scale
        (``weekly``: weeks 1-8; ``daily``: ``smooth_days`` running means at every lead day) for panel (b)
        (every bias treatment and season set) and panel (c) (AMIP climatology)."""
        cfg = self.cfg
        axis, coords = self.axis(scale)

        def build():
            tables = []
            for name, set_seasons in self.sets.items():
                for b in cfg.biases:
                    fa, oa, labels = self.stack(set_seasons, scale, "clim_amip" if b == "amip" else "clim_obs")
                    res = era5_benefits(fa, oa, self.weights, cfg.contrasts, labels, n_boot=cfg.n_boot, seed=cfg.seed)
                    tables.append(results_table(res, coords, "ERA5", name, axis).assign(bias=b, scale=scale))
                if cfg.model_truth_role in self.roles:
                    fa, _, labels = self.stack(set_seasons, scale, "clim_amip")
                    res = model_member_benefits(fa, cfg.model_truth_role, self.weights, cfg.contrasts, labels,
                                                n_boot=cfg.n_boot, seed=cfg.seed)
                    tables.append(results_table(res, coords, "model member", name, axis).assign(bias="amip", scale=scale))
                print(f"[SCORES] {scale} {name}")
            return pd.concat(tables, ignore_index=True)
        return cached_table(self.data_dir / f"init_benefit_curves_{scale}.csv", build, self._stamp(*_SCORE_KEYS),
                            force=self.force)

    @staticmethod
    def panel(table: pd.DataFrame, panel: str, set_name: str, bias: str = "amip") -> dict | None:
        """The ``skill``, ``skill_ci``, ``contrasts``, ``contrast_ci`` dictionaries of one panel
        (``ERA5`` or ``model member``) from a ``curves`` table, as the plotting functions take them."""
        t = table[(table.panel == panel) & (table.set == set_name) & (table.bias == bias)]
        if t.empty:
            return None
        out = {}
        for kind, key in (("skill", "skill"), ("contrast", "contrasts")):
            sub = t[t.kind == kind]
            names = list(dict.fromkeys(sub.name))
            g = {n: sub[sub.name == n].sort_values("step") for n in names}
            out[key] = {n: g[n].value.to_numpy() for n in names}
            out[f"{kind}_ci"] = {n: np.vstack([g[n].ci_low.to_numpy(), g[n].ci_high.to_numpy()]) for n in names}
        return out

    def maps(self, bias: str, set_name: str) -> xr.Dataset:
        """Gridpoint ACC across the starts of each week, per role, and its contrasts (cached)."""
        cfg = self.cfg

        def build():
            set_seasons = self.sets[set_name]
            fa, oa, labels = self.stack(set_seasons, "weekly", "clim_amip" if bias == "amip" else "clim_obs")
            acc = {r: gridpoint_acc_np(a.mean(1), oa) for r, a in fa.items()}
            g = self.grid
            ds = xr.Dataset({k: (("week", "lat", "lon"), v) for k, v in contrast_values(acc, cfg.contrasts).items()},
                            coords={"week": self.weeks, "lat": g.lat, "lon": g.lon}).where(g.mask)
            ds.attrs.update(n_starts=len(labels), note="gridpoint ACC across starts; very few starts")
            return ds
        return cached_dataset(self.data_dir / f"init_benefit_maps_{bias}_{set_name}.nc", build,
                              self._stamp(*_DATA_KEYS, "contrasts"), force=self.force)


__all__ = ["InitBenefitCase", "SCALES"]
