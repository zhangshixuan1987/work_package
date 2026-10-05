"""Reusable model, observation, and mode registries for MOV workflows."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Mapping, Sequence


DEFAULT_PERIOD = (1985, 2014)
DEFAULT_MODEL_DATA_ROOT = Path("/lcrc/group/e3sm2/ac.qtang/E3SMv3")
DEFAULT_MODEL_SUBDIR = Path("post/atm/180x360_aave/cmip_ts/monthly")


@dataclass(frozen=True)
class MOVExperiment:
    """One explicitly configured experiment usable by either MOV backend."""

    key: str
    label: str
    model_name: str
    pcmdi_tag: str
    data_root: Path = DEFAULT_MODEL_DATA_ROOT
    data_subdir: Path = DEFAULT_MODEL_SUBDIR

    def raw_dataset(self) -> Dict[str, Any]:
        """Return the dataset entry expected by :class:`ModeAnalyzer`."""
        return {
            "name": self.model_name,
            "dir": Path(self.data_root),
            "subdir": Path(self.data_subdir),
        }


def experiments_to_raw_datasets(
    experiments: Sequence[MOVExperiment],
) -> Dict[str, Dict[str, Any]]:
    """Build a raw-analysis dataset registry from explicit experiments."""
    return {experiment.key: experiment.raw_dataset() for experiment in experiments}

# Named preset; other projects can supply the same mapping interface.
V3_RRM_MODELS: Dict[str, Dict[str, Any]] = {
    "LR": {"name": "v3.LR.amip_0101"},
    "NARRM": {"name": "v3.NARRM.amip_0101"},
    "NARRM_r0125": {"name": "v3.NARRM_r0125.amip_0101"},
    "EARRM": {"name": "v3.EARRM.amip_0101"},
    "AMZRRM": {"name": "v3.AMZRRM.amip_0101"},
}
V3_RRM_MODEL_ORDER = tuple(V3_RRM_MODELS)
V3_RRM_PCMDI_MODEL_TAGS = tuple(cfg["name"] for cfg in V3_RRM_MODELS.values())

OBS_SOURCES = {
    "psl": {
        "name": "NOAA-20C",
        "data": Path("/lcrc/group/e3sm/diagnostics/observations/obs4MIPs_PCMDI_monthly/NOAA-20C/psl_183601_201512.nc"),
    },
    "ts": {
        "name": "HadISST2",
        "data": Path("/lcrc/group/e3sm/diagnostics/observations/obs4MIPs_PCMDI_monthly/HadISST2/ts_186901_202212.nc"),
    },
}

MODE_CONFIG = {
    "AMO": {"eof": "eof1", "eof_num": 1, "var": "ts", "obs": "HadISST2", "lat_bnds": (0, 70), "lon_bnds": (-80, 0), "season": "annual"},
    "PDO": {"eof": "eof1", "eof_num": 1, "var": "ts", "obs": "HadISST2", "lat_bnds": (20, 70), "lon_bnds": (110, -100), "season": "annual"},
    "NPGO": {"eof": "eof2", "eof_num": 2, "var": "ts", "obs": "HadISST2", "lat_bnds": (20, 70), "lon_bnds": (110, -100), "season": "annual"},
    "NAM": {"eof": "eof1", "eof_num": 1, "var": "psl", "obs": "NOAA-20C", "lat_bnds": (20, 90), "lon_bnds": (-180, 180), "season": "DJF"},
    "NAO": {"eof": "eof1", "eof_num": 1, "var": "psl", "obs": "NOAA-20C", "lat_bnds": (20, 80), "lon_bnds": (-80, 40), "season": "DJF"},
    "NPO": {"eof": "eof2", "eof_num": 2, "var": "psl", "obs": "NOAA-20C", "lat_bnds": (20, 80), "lon_bnds": (120, -120), "season": "DJF"},
    "PNA": {"eof": "eof1", "eof_num": 1, "var": "psl", "obs": "NOAA-20C", "lat_bnds": (20, 80), "lon_bnds": (120, -60), "season": "DJF"},
    "SAM": {"eof": "eof1", "eof_num": 1, "var": "psl", "obs": "NOAA-20C", "lat_bnds": (-90, -20), "lon_bnds": (-180, 180), "season": "annual"},
    "PSA1": {"eof": "eof2", "eof_num": 2, "var": "psl", "obs": "NOAA-20C", "lat_bnds": (-90, -20), "lon_bnds": (120, -60), "season": "annual"},
    "PSA2": {"eof": "eof3", "eof_num": 3, "var": "psl", "obs": "NOAA-20C", "lat_bnds": (-90, -20), "lon_bnds": (120, -60), "season": "annual"},
}


def model_order(models: Mapping[str, Mapping[str, Any]]) -> tuple[str, ...]:
    """Return stable display order for any model registry."""
    return tuple(models)


def pcmdi_model_tags(models: Mapping[str, Mapping[str, Any]]) -> tuple[str, ...]:
    """Return PCMDI directory/model tags for any model registry."""
    return tuple(str(cfg["name"]) for cfg in models.values())


def build_model_datasets(
    models: Mapping[str, Mapping[str, Any]],
    *,
    root: Path = DEFAULT_MODEL_DATA_ROOT,
    subdir: Path = DEFAULT_MODEL_SUBDIR,
) -> Dict[str, Dict[str, Any]]:
    """Attach a data root and subdirectory to an arbitrary model registry."""
    return {
        key: {"name": cfg["name"], "dir": Path(root), "subdir": Path(subdir)}
        for key, cfg in models.items()
    }


def v3_rrm_model_datasets(root: Path = DEFAULT_MODEL_DATA_ROOT) -> Dict[str, Dict[str, Any]]:
    """Convenience preset for the five v3 RRM configurations."""
    return build_model_datasets(V3_RRM_MODELS, root=root)
