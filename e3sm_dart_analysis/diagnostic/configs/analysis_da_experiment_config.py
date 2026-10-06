"""Generic experiment setup helpers for analysis_da observation notebooks."""

from __future__ import annotations

import copy
import os
import re
from pathlib import Path
from typing import Dict, Mapping, MutableMapping, Optional, Sequence


DEFAULT_RESOLUTION = "ne30pg2_r05_IcoswISC30E3r5"
DEFAULT_MODEL_RESOLUTION = f"F20TR_{DEFAULT_RESOLUTION}"
DEFAULT_MACHINE = "compy"
DEFAULT_ATM_SUBDIR = "archive/post/atm/180x360_aave"
DEFAULT_LND_SUBDIR = "archive/post/lnd/180x360_aave"
DEFAULT_OBS_DIAG_SETS = ("obs_seq", "obs_diag", "obs_common", "closest_member", "cam6_common")

_SEASON_RE = re.compile(r"(?:-|_)(S\d+)\b")
_RUN_PERIOD_RE = re.compile(r"^\d{6}-\d{6}$")
_OBS_PERIOD_RE = re.compile(r"^\d{10}-\d{10}$")


def _deep_update(base: MutableMapping, updates: Mapping) -> MutableMapping:
    for key, value in updates.items():
        if isinstance(value, Mapping) and isinstance(base.get(key), MutableMapping):
            _deep_update(base[key], value)
        else:
            base[key] = copy.deepcopy(value)
    return base


def merge_experiments(
    experiments: Optional[Mapping[str, Mapping]],
    overrides: Optional[Mapping[str, Optional[Mapping]]] = None,
) -> Dict[str, dict]:
    """Return a deep-copied experiment mapping with optional add/update/remove overrides."""

    specs = copy.deepcopy(dict(experiments or {}))
    for name, override in (overrides or {}).items():
        if override is None:
            specs.pop(name, None)
        elif name in specs:
            _deep_update(specs[name], override)
        else:
            specs[name] = copy.deepcopy(override)
    return specs


def _require_experiments(experiments: Optional[Mapping[str, Mapping]]) -> Mapping[str, Mapping]:
    if not experiments:
        raise ValueError(
            "No experiments were provided. Define USER_EXPERIMENTS in the notebook "
            "and pass it to the experiment setup helper."
        )
    return experiments


def _season_from_name(name: str) -> Optional[str]:
    match = _SEASON_RE.search(name)
    return match.group(1) if match else None


def _valid_period(period: str, pattern: re.Pattern, expected: str, exp_name: str) -> str:
    if not pattern.match(period):
        raise ValueError(f"Invalid period '{period}' for {exp_name}; expected '{expected}'.")
    return period


def _build_run(
    exp_name: str,
    spec: Optional[Mapping],
    *,
    data_path: str,
    resolution: str,
    machine: str,
    atm_subdir: str,
    lnd_subdir: str,
    group_key: Optional[str],
) -> Optional[dict]:
    if spec is None:
        return None

    period = _valid_period(str(spec["period"]), _RUN_PERIOD_RE, "YYYYMM-YYYYMM", exp_name)
    name = str(spec["name"])
    compset = str(spec["compset"])
    run_id = str(spec.get("run_id") or f"{name}_{compset}_{resolution}_{machine}")
    alias = spec.get("alias", spec.get("alia", name))

    out = {
        "run_id": run_id,
        "name": name,
        "alias": alias,
        "key": spec.get("key", alias),
        "group_key": group_key,
        "compset": compset,
        "period": period,
        "atm": atm_subdir,
        "lnd": lnd_subdir,
        "atm_path": os.path.join(data_path, run_id, atm_subdir),
        "lnd_path": os.path.join(data_path, run_id, lnd_subdir),
    }
    if "obs_diag" in spec:
        out["obs_diag"] = list(spec["obs_diag"])
    return out


def build_run_experiment_info(
    data_path: str | os.PathLike,
    experiments: Mapping[str, Mapping],
    *,
    resolution: str = DEFAULT_RESOLUTION,
    machine: str = DEFAULT_MACHINE,
    atm_subdir: str = DEFAULT_ATM_SUBDIR,
    lnd_subdir: str = DEFAULT_LND_SUBDIR,
    default_run_order: Sequence[str] = ("fc", "wc", "da"),
    experiment_overrides: Optional[Mapping[str, Optional[Mapping]]] = None,
) -> Dict[str, dict]:
    """Build run-oriented metadata from a notebook-provided experiment mapping."""

    specs = merge_experiments(_require_experiments(experiments), experiment_overrides)
    data_path = str(data_path)

    exp_dict: Dict[str, dict] = {}
    for exp_name, meta in sorted(specs.items()):
        group_key = meta.get("key")
        runs = {
            "da": _build_run(
                exp_name,
                meta.get("da_run"),
                data_path=data_path,
                resolution=resolution,
                machine=machine,
                atm_subdir=atm_subdir,
                lnd_subdir=lnd_subdir,
                group_key=group_key,
            ),
            "fc": _build_run(
                exp_name,
                meta.get("fc_run"),
                data_path=data_path,
                resolution=resolution,
                machine=machine,
                atm_subdir=atm_subdir,
                lnd_subdir=lnd_subdir,
                group_key=group_key,
            ),
            "wc": _build_run(
                exp_name,
                meta.get("wc_run"),
                data_path=data_path,
                resolution=resolution,
                machine=machine,
                atm_subdir=atm_subdir,
                lnd_subdir=lnd_subdir,
                group_key=group_key,
            ),
        }
        default_run = next((runs[k] for k in default_run_order if runs.get(k)), None)
        exp_dict[exp_name] = {
            "nens": meta["nens"],
            "season": _season_from_name(exp_name),
            "group_key": group_key,
            "runs": runs,
            "default_run": default_run,
            "key": (default_run or {}).get("key"),
            "period": (default_run or {}).get("period"),
        }

    return exp_dict


# Backward-friendly name for notebooks that already call extract_exp_info(...).
extract_exp_info = build_run_experiment_info


def build_obs_diag_config(
    experiments: Mapping[str, Mapping],
    exp_name: Optional[str] = None,
    *,
    data_path: str | os.PathLike | None = None,
    case_name: str = "JAN2011",
    resolution: str = DEFAULT_MODEL_RESOLUTION,
    machine: str = DEFAULT_MACHINE,
    diag_key: str = "obs_diag_output",
    frequency: str = "6hourly",
    region: str = "Northern Hemisphere",
    path_template: Optional[str] = None,
    file_template: str = "%(RUNNAME).dart.e.eam_%(KEY).%(TIME).nc",
    diag_sets: Sequence[str] = DEFAULT_OBS_DIAG_SETS,
    experiment_overrides: Optional[Mapping[str, Optional[Mapping]]] = None,
) -> dict:
    """Build the obs-diagnostic config shape from notebook-provided experiments."""

    data_root = Path(data_path or os.environ.get("E3SM_DART_DATA_PATH", "/compyfs/zhan391/v3_dart_cda_scratch"))
    if path_template is None:
        path_template = str(
            data_root / "%(RUNNAME)" / "%(CASENAME)" / "%(COMP)" / "dart_diagnostics" / "%(DIAG)"
        )

    specs = merge_experiments(_require_experiments(experiments), experiment_overrides)
    for name, meta in specs.items():
        _valid_period(str(meta["period"]), _OBS_PERIOD_RE, "YYYYMMDDHH-YYYYMMDDHH", name)
        sets = tuple(meta.get("diag_sets", diag_sets))
        for index, diag_set in enumerate(sets, start=1):
            meta.setdefault(f"diag{index}", diag_set)
        meta["diag_sets"] = [meta[f"diag{index}"] for index in range(1, len(sets) + 1)]

    if exp_name:
        specs = {exp_name: specs[exp_name]} if exp_name in specs else {}

    return {
        "global": {
            "case_name": case_name,
            "resolution": resolution,
            "machine": machine,
            "diag_key": diag_key,
            "frequency": frequency,
            "region": region,
            "path_template": path_template,
            "file_template": file_template,
        },
        "experiments": specs,
    }


def extract_obs_diag_config(exp_name: Optional[str] = None, **kwargs) -> dict:
    """Backward-friendly wrapper with ``exp_name`` as the first argument."""

    return build_obs_diag_config(exp_name=exp_name, **kwargs)
