"""Collect observation and per-member model file lists for land diagnostics.

``S2SFileCollector`` resolves, for each experiment of a forecast group, the
observation files (via the observation registry) and the post-processed
land time-series files of every ensemble member.

Output of ``collect`` (per model variable)::

    {experiment: {"obs": [obs files], "model": {member: [files]}, "ts_dir": ..., "nens": ...}}
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple


@dataclass(frozen=True)
class FileCollectionConfig:
    group: str
    freq: str
    run: str
    obs: str
    period: str
    model_template: str
    ens_prefix: str = "en"
    ens_width: int = 2
    ens_start: int = 1

    def parse_period(self) -> Tuple[int, int, int, int]:
        m = re.match(r"^(\d{4})(\d{2})-(\d{4})(\d{2})$", self.period)
        if not m:
            raise ValueError(f"Bad period '{self.period}', expected 'YYYYMM-YYYYMM'")
        y0, m0, y1, m1 = map(int, m.groups())
        if (y1, m1) < (y0, m0):
            raise ValueError(f"Bad period '{self.period}' (end < start)")
        return y0, m0, y1, m1

    def years(self) -> List[int]:
        y0, _, y1, _ = self.parse_period()
        return list(range(y0, y1 + 1))

    def ens_labels(self, nens: int) -> List[str]:
        return [
            f"{self.ens_prefix}{i:0{self.ens_width}d}"
            for i in range(self.ens_start, nens + self.ens_start)
        ]


class S2SFileCollector:
    """
    Collect obs + model files for S2S / diagnostics.

    s2s_var_dict: dict[model_var] -> obs_var
    """

    def __init__(
        self,
        *,
        exp_list: dict,
        exp_dict: dict,
        obs_registry,
        s2s_var_dict: Dict[str, str],
        get_obs_file_func,
    ):
        self.exp_list = exp_list
        self.exp_dict = exp_dict
        self.obs_registry = obs_registry
        self.s2s_var_dict = s2s_var_dict
        self.get_obs_file = get_obs_file_func

    def resolve_obs_file(self, obs: str, freq: str, year: int, var: Optional[str] = None) -> str:
        return self.get_obs_file(self.obs_registry, obs, freq=freq, year=year, var=var)

    def model_ts_dir(self, run_meta, freq: str, ens: Optional[str] = None) -> str:
        """Land time-series directory; per-member ``<run_id>/<ens>/<lnd>/`` when that exists.

        Ensemble cases keep each member's archive under ``<run_id>/EN##/``; the
        registry's ``lnd_path`` points at ``<run_id>/<lnd>/`` (single-instance
        layout), which is used as the fallback.
        """
        shared = os.path.join(run_meta.lnd_path, "ts", freq)
        if ens is None:
            return shared
        case_dir = run_meta.lnd_path[: -len(run_meta.lnd)].rstrip(os.sep)
        member = os.path.join(case_dir, ens, run_meta.lnd, "ts", freq)
        return member if os.path.isdir(member) else shared

    def model_files_for_ens(
        self,
        cfg: FileCollectionConfig,
        *,
        run_meta,
        var: str,
        ens: str,
    ) -> List[str]:
        ts_dir = self.model_ts_dir(run_meta, cfg.freq, ens)
        years = cfg.years()
        return [
            os.path.join(
                ts_dir,
                cfg.model_template % {"var": var, "ens": ens, "year": y, "period": cfg.period},
            )
            for y in years
        ]

    def collect_one_var(self, cfg: FileCollectionConfig, *, var: str, obs_var: str, verbose: bool = True):
        years = cfg.years()
        models = self.exp_list[cfg.group]["models"]

        # obs files
        obs_files = [self.resolve_obs_file(cfg.obs, cfg.freq, y, obs_var) for y in years]
        missing_obs = [f for f in obs_files if not os.path.exists(f)]
        if missing_obs:
            if verbose:
                print(f"[MISSING OBS] {var} (obs_var={obs_var})")
                for f in missing_obs:
                    print("  ", f)
            return None

        out_var = {}

        for exp in models:
            run_meta = self.exp_dict[exp]["runs"][cfg.run]
            if run_meta is None:
                if verbose:
                    print(f"SKIP {exp}: no run='{cfg.run}'")
                continue

            nens = int(self.exp_dict[exp]["nens"])
            ts_dir = self.model_ts_dir(run_meta, cfg.freq)

            model_by_ens = {}
            missing_any = False

            for ens in cfg.ens_labels(nens):
                files = self.model_files_for_ens(cfg, run_meta=run_meta, var=var, ens=ens)

                missing = [f for f in files if not os.path.exists(f)]
                if missing:
                    missing_any = True
                    if verbose:
                        print(f"[MISSING MOD] {exp} {var} {ens} (ts_dir={ts_dir})")
                        for f in missing[:5]:
                            print("   ", f)
                        if len(missing) > 5:
                            print(f"   ... {len(missing)-5} more")

                model_by_ens[ens] = files

            out_var[exp] = {"obs": obs_files, "model": model_by_ens, "ts_dir": ts_dir, "nens": nens}

            if verbose and missing_any:
                print(f"[WARN] {exp} {var}: some ensemble members missing files")

        return out_var

    def collect(
        self,
        cfg: FileCollectionConfig,
        *,
        vars_to_process: Optional[List[str]] = None,
        verbose: bool = True,
    ) -> Dict[str, dict]:
        out = {}
        for var, obs_var in self.s2s_var_dict.items():
            if vars_to_process is not None and var not in vars_to_process:
                continue
            res = self.collect_one_var(cfg, var=var, obs_var=obs_var, verbose=verbose)
            if res is None:
                continue
            out[var] = res
        return out
