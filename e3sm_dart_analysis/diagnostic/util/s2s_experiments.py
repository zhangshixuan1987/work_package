from __future__ import annotations
import os
import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple


# -------------------------
# Types (portable specs + resolved metadata)
# -------------------------

@dataclass(frozen=True)
class RunSpec:
    """Input spec (portable, no paths resolved)."""
    name: str
    compset: str
    period: str
    run_id: Optional[str] = None
    obs_diag: Optional[List[str]] = None
    extra: Dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class ExperimentSpec:
    """Input spec for a single experiment entry."""
    nens: int
    key: str
    da_run: Optional[RunSpec] = None
    fc_run: Optional[RunSpec] = None
    wc_run: Optional[RunSpec] = None


@dataclass(frozen=True)
class RunMeta:
    """Resolved run metadata (paths computed)."""
    run_id: str
    name: str
    compset: str
    period: str
    atm: str
    lnd: str
    atm_path: str
    lnd_path: str
    obs_diag: Optional[List[str]] = None
    extra: Dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class ExperimentMeta:
    """Resolved experiment metadata."""
    nens: int
    season: Optional[str]
    group_key: str
    runs: Dict[str, Optional[RunMeta]]  # keys: da|fc|wc
    default_run: Optional[RunMeta]
    key: str
    period: Optional[str]


# -------------------------
# Registry (edit/extend for other projects)
# -------------------------

DEFAULT_EXPERIMENTS: Dict[str, ExperimentSpec] = {
    "CTRL": ExperimentSpec(
        nens=1,
        key="ctrl",
        da_run=None,
        fc_run=RunSpec(compset="F20TR", name="CTRL", period="201201-201212"),
        wc_run=RunSpec(compset="WCYCL20TR", name="CTRL", period="201201-201212"),
    ),
    "CTRL10-S0": ExperimentSpec(
        nens=10,
        key="dart_en10",
        da_run=RunSpec(
            compset="F20TR",
            name="CTRLEN10",
            period="201112-201112",
            obs_diag=["obs_seq", "obs_diag", "obs_common", "closest_member"],
        ),
        fc_run=RunSpec(compset="F20TR", name="CTRLEN10_15day", period="201201-201202"),
        wc_run=RunSpec(compset="WCYCL20TR", name="CTRLEN10_15day", period="201201-201202"),
    ),
    "CAPT10-S0": ExperimentSpec(
        nens=10,
        key="dart_en10",
        da_run=None,
        fc_run=RunSpec(compset="F20TR", name="CAPTEN10_15day", period="201201-201202"),
        wc_run=RunSpec(compset="WCYCL20TR", name="CAPTEN10_15day", period="201201-201202"),
    ),
    "DART10-S0": ExperimentSpec(
        nens=10,
        key="dart_en10",
        da_run=RunSpec(
            compset="F20TR",
            name="DARTEN10",
            period="201112-201112",
            obs_diag=["obs_seq", "obs_diag", "obs_common", "closest_member"],
        ),
        fc_run=None,
        wc_run=None,
    ),
    "DART20-S0": ExperimentSpec(
        nens=20,
        key="dart_en20",
        da_run=RunSpec(
            compset="F20TR",
            name="DARTEN20",
            period="201112-201112",
            obs_diag=["obs_seq", "obs_diag", "obs_common", "closest_member"],
        ),
        fc_run=RunSpec(compset="F20TR", name="DARTEN20_15day", period="201201-201202"),
        wc_run=RunSpec(compset="WCYCL20TR", name="DARTEN20_15day", period="201201-201202"),
    ),
    "DART40-S0": ExperimentSpec(
        nens=40,
        key="dart_en40",
        da_run=RunSpec(
            compset="F20TR",
            name="DARTEN40",
            period="201112-201112",
            obs_diag=["obs_seq", "obs_diag", "obs_common", "closest_member"],
        ),
        fc_run=RunSpec(compset="F20TR", name="DARTEN40_15day", period="201201-201202"),
        wc_run=RunSpec(compset="WCYCL20TR", name="DARTEN40_15day", period="201201-201202"),
    ),
    "CAM80-S0": ExperimentSpec(
        nens=80,
        key="dart_en80",
        da_run=RunSpec(
            compset="F20TR",
            name="DARTEN40",
            period="201112-201112",
            run_id="f.e21.FHIST_BGC.f09_025.CAM6assim.011",
            obs_diag=["obs_seq", "obs_diag", "obs_common", "closest_member"],
        ),
        fc_run=None,
        wc_run=None,
    ),
    "DART40INF0p6-S0": ExperimentSpec(
        nens=40,
        key="dart_en40",
        da_run=RunSpec(
            compset="F20TR",
            name="DARTEN40_INF0p6",
            period="201112-201112",
            obs_diag=["obs_seq", "obs_diag", "obs_common", "closest_member"],
            extra={"alia": "DARTEN40"},
        ),
        fc_run=None,
        wc_run=None,
    ),
    "CTRL10-S1": ExperimentSpec(
        nens=10,
        key="ctrl_en10",
        da_run=None,
        fc_run=RunSpec(compset="F20TR", name="CTRLEN10s1_15day", period="201206-201207"),
        wc_run=RunSpec(compset="WCYCL20TR", name="CTRLEN10s1_15day", period="201206-201207"),
    ),
    "CAPT10-S1": ExperimentSpec(
        nens=10,
        key="capt_en10",
        da_run=None,
        fc_run=RunSpec(compset="F20TR", name="CAPTEN10S1_15day", period="201206-201207"),
        wc_run=RunSpec(compset="WCYCL20TR", name="CAPTEN10S1_15day", period="201206-201207"),
    ),
    "DART40-S1": ExperimentSpec(
        nens=40,
        key="dart_en40",
        da_run=RunSpec(
            compset="F20TR",
            name="DARTEN40S1",
            period="201205-201205",
            obs_diag=["obs_seq", "obs_diag", "obs_common", "closest_member"],
        ),
        fc_run=RunSpec(compset="F20TR", name="DARTEN40S1_15day", period="201206-201207"),
        wc_run=RunSpec(compset="WCYCL20TR", name="DARTEN40S1_15day", period="201206-201207"),
    ),
    "CAM80-S1": ExperimentSpec(
        nens=80,
        key="dart_en80",
        da_run=RunSpec(
            compset="F20TR",
            name="DARTEN40",
            period="201205-201205",
            run_id="f.e22.FHIST_BGC.f09_025.CAM6assim.011",
            obs_diag=["obs_seq", "obs_diag", "obs_common", "closest_member"],
        ),
        fc_run=None,
        wc_run=None,
    ),
}


# -------------------------
# Builder / helpers
# -------------------------

_SEASON_RE = re.compile(r"(?:-|_)(S\d+)\b")
_PERIOD_RE = re.compile(r"^\d{6}-\d{6}$")


def season_from_name(exp_name: str) -> Optional[str]:
    m = _SEASON_RE.search(exp_name)
    return m.group(1) if m else None


def require_valid_period(p: str) -> str:
    if not _PERIOD_RE.match(p):
        raise ValueError(f"Invalid period '{p}' for experiment; expected 'YYYYMM-YYYYMM'.")
    return p


def build_run(
    data_path: str,
    spec: Optional[RunSpec],
    *,
    resolution: str,
    machine: str,
    atm_subdir: str,
    lnd_subdir: str,
) -> Optional[RunMeta]:
    if spec is None:
        return None

    period = require_valid_period(spec.period)
    run_id = spec.run_id or f"{spec.name}_{spec.compset}_{resolution}_{machine}"

    atm_path = os.path.join(data_path, run_id, atm_subdir)
    lnd_path = os.path.join(data_path, run_id, lnd_subdir)

    return RunMeta(
        run_id=run_id,
        name=spec.name,
        compset=spec.compset,
        period=period,
        atm=atm_subdir,
        lnd=lnd_subdir,
        atm_path=atm_path,
        lnd_path=lnd_path,
        obs_diag=list(spec.obs_diag) if spec.obs_diag else None,
        extra=dict(spec.extra) if spec.extra else {},
    )


def build_experiments(
    data_path: str,
    experiments: Dict[str, ExperimentSpec] = DEFAULT_EXPERIMENTS,
    *,
    resolution: str = "ne30pg2_r05_IcoswISC30E3r5",
    machine: str = "compy",
    atm_subdir: str = "archive/post/atm/180x360_aave",
    lnd_subdir: str = "archive/post/lnd/180x360_aave",
    default_preference: Tuple[str, ...] = ("fc", "wc", "da"),
    return_mode: str = "legacy",
) -> Dict[str, dict]:
    """
    Build resolved experiment metadata.

    Parameters
    ----------
    data_path
        Root path that contains run_id directories.
    experiments
        Registry mapping exp_name -> ExperimentSpec.
    return_mode
        - "legacy": return dict matching your original shape (runs/default_run as objects)
        - "typed":  return dict of ExperimentMeta objects

    Returns
    -------
    Dict[str, dict] (legacy) or Dict[str, ExperimentMeta] (typed)
    """
    typed_out: Dict[str, ExperimentMeta] = {}

    for exp_name in sorted(experiments.keys()):
        spec = experiments[exp_name]

        runs = {
            "da": build_run(
                data_path, spec.da_run,
                resolution=resolution, machine=machine,
                atm_subdir=atm_subdir, lnd_subdir=lnd_subdir
            ),
            "fc": build_run(
                data_path, spec.fc_run,
                resolution=resolution, machine=machine,
                atm_subdir=atm_subdir, lnd_subdir=lnd_subdir
            ),
            "wc": build_run(
                data_path, spec.wc_run,
                resolution=resolution, machine=machine,
                atm_subdir=atm_subdir, lnd_subdir=lnd_subdir
            ),
        }

        default_run = None
        for k in default_preference:
            if runs.get(k) is not None:
                default_run = runs[k]
                break

        typed_out[exp_name] = ExperimentMeta(
            nens=spec.nens,
            season=season_from_name(exp_name),
            group_key=spec.key,
            runs=runs,
            default_run=default_run,
            key=spec.key,
            period=default_run.period if default_run else None,
        )

    if return_mode == "typed":
        # returns: Dict[str, ExperimentMeta]
        return typed_out  # type: ignore[return-value]

    if return_mode != "legacy":
        raise ValueError("return_mode must be one of: 'legacy', 'typed'.")

    # Legacy compatibility: Dict[str, dict] with same keys you used before
    legacy_out: Dict[str, dict] = {}
    for exp_name, em in typed_out.items():
        legacy_out[exp_name] = {
            "nens": em.nens,
            "season": em.season,
            "group_key": em.group_key,
            "runs": em.runs,
            "default_run": em.default_run,
            "key": em.key,
            "period": em.period,
        }

    return legacy_out


# -------------------------
# Optional convenience helpers (small, useful, non-invasive)
# -------------------------

def select_experiments(
    exp_dict: Dict[str, dict],
    *,
    season: Optional[str] = None,
    key: Optional[str] = None,
) -> Dict[str, dict]:
    """
    Filter a legacy exp_dict returned by build_experiments(return_mode='legacy').
    """
    out: Dict[str, dict] = {}
    for name, meta in exp_dict.items():
        if season is not None and meta.get("season") != season:
            continue
        if key is not None and meta.get("key") != key and meta.get("group_key") != key:
            continue
        out[name] = meta
    return out


def get_run(
    exp_dict: Dict[str, dict],
    exp_name: str,
    which: str = "default",
) -> Optional[RunMeta]:
    """
    Convenience accessor for legacy exp_dict.
      which: 'default' | 'da' | 'fc' | 'wc'
    """
    meta = exp_dict[exp_name]
    if which == "default":
        return meta.get("default_run")
    return meta["runs"].get(which)
