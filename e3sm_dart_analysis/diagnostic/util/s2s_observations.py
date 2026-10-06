# ---------------------------------------------------------------------
# Observations registry + helpers (drop into DataUtil.py)
# ---------------------------------------------------------------------

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import Dict, Optional, Tuple, Iterable, List, Any


@dataclass(frozen=True)
class ObsDatasetSpec:
    path: str
    template: str
    period: str  # e.g. "1979-2018"


@dataclass(frozen=True)
class ObsSourceSpec:
    monthly: Optional[ObsDatasetSpec] = None
    daily: Optional[ObsDatasetSpec] = None


def parse_year_range(period: str) -> Tuple[int, int]:
    """
    Parse 'YYYY-YYYY' -> (start_year, end_year) inclusive.
    """
    m = re.match(r"^\s*(\d{4})\s*-\s*(\d{4})\s*$", period)
    if not m:
        raise ValueError(f"Invalid obs period '{period}' (expected 'YYYY-YYYY').")
    y0, y1 = int(m.group(1)), int(m.group(2))
    if y1 < y0:
        raise ValueError(f"Invalid obs period '{period}' (end < start).")
    return y0, y1


def year_in_period(year: int, period: str) -> bool:
    y0, y1 = parse_year_range(period)
    return y0 <= year <= y1


def _format_template(template: str, *, year: int, var: Optional[str] = None, **kwargs) -> str:
    """
    Supports templates written like:
      'file_%(year).nc'
      '%(var)_%(year)01_%(year)12.nc'
    """
    mapping: Dict[str, Any] = {"year": year}
    if var is not None:
        mapping["var"] = var
    mapping.update(kwargs)
    try:
        return template % mapping
    except KeyError as e:
        raise KeyError(
            f"Missing template key {e} for template='{template}'. "
            f"Provided keys={sorted(mapping.keys())}"
        ) from e


def build_obs_registry(raw: Dict[str, dict]) -> Dict[str, ObsSourceSpec]:
    """
    Convert your nested dict into typed specs (optional but helpful).
    """
    out: Dict[str, ObsSourceSpec] = {}
    for name, spec in raw.items():
        monthly = None
        daily = None
        if "monthly" in spec:
            monthly = ObsDatasetSpec(**spec["monthly"])
        if "daily" in spec:
            daily = ObsDatasetSpec(**spec["daily"])
        out[name] = ObsSourceSpec(monthly=monthly, daily=daily)
    return out


def get_obs_file(
    obs_registry: Dict[str, ObsSourceSpec],
    source: str,
    *,
    freq: str,
    year: int,
    var: Optional[str] = None,
    check_period: bool = True,
    **template_kwargs,
) -> str:
    """
    Return full file path for a given obs source/frequency/year.

    Examples:
      get_obs_file(reg, "ERA5", freq="daily", year=2012, var="PRECT")
      get_obs_file(reg, "ERA5", freq="monthly", year=2012)
      get_obs_file(reg, "GPCP", freq="daily", year=2009)
    """
    if source not in obs_registry:
        raise KeyError(f"Unknown obs source '{source}'. Available={list(obs_registry)}")

    src = obs_registry[source]
    if freq not in ("daily", "monthly"):
        raise ValueError("freq must be 'daily' or 'monthly'")

    ds: Optional[ObsDatasetSpec] = getattr(src, freq)
    if ds is None:
        raise ValueError(f"Obs source '{source}' does not provide freq='{freq}'")

    if check_period and not year_in_period(year, ds.period):
        raise ValueError(f"{source} {freq} does not cover year={year} (period={ds.period})")

    fname = _format_template(ds.template, year=year, var=var, **template_kwargs)
    return os.path.join(ds.path, fname)


def list_obs_sources(obs_registry: Dict[str, ObsSourceSpec], *, freq: Optional[str] = None) -> List[str]:
    """
    List available obs sources, optionally filtered by frequency.
    """
    if freq is None:
        return sorted(obs_registry.keys())
    if freq not in ("daily", "monthly"):
        raise ValueError("freq must be 'daily' or 'monthly'")
    out = []
    for k, v in obs_registry.items():
        if getattr(v, freq) is not None:
            out.append(k)
    return sorted(out)


def obs_coverage(obs_registry: Dict[str, ObsSourceSpec], source: str, *, freq: str) -> Tuple[int, int]:
    """
    Return (start_year, end_year) inclusive for a given source/freq.
    """
    src = obs_registry[source]
    ds: Optional[ObsDatasetSpec] = getattr(src, freq)
    if ds is None:
        raise ValueError(f"Obs source '{source}' does not provide freq='{freq}'")
    return parse_year_range(ds.period)


# ---------------------------------------------------------------------
# Your raw obs group as a registry constant
# ---------------------------------------------------------------------

RAW_OBS_GROUP = {
    "ERA5": {
        "monthly": {
            "path": "/compyfs/zhan391/acme_init/Observations/ERA5/monthly",
            "template": "ERA5_analysis_monthly_%(year)d.nc",
            "period": "1979-2018",
        },
        "daily": {
            "path": "/compyfs/zhan391/acme_init/Observations/ERA5/daily",
            "template": "%(var)s_%(year)d01_%(year)d12.nc",
            "period": "1979-2018",
        },
    },
    "CERES-OAFlux": {
        "monthly": {
            "path": "/compyfs/zhan391/acme_init/Observations/CERES-OAFlux/monthly",
            "template": "CERES-OAFlux_%(year)d.nc",
            "period": "2001-2018",
        }
    },
    "NOAA-OLR": {
        "daily": {
            "path": "/compyfs/zhan391/acme_init/Observations/NOAA-OLR/daily",
            "template": "FLUT_%(year)d01_%(year)d12.nc",
            "period": "2009-2015",
        }
    },
    "GPCP": {
        "monthly": {
            "path": "/compyfs/zhan391/acme_init/Observations/GPCP/monthly",
            "template": "PRECT.monthly.%(year)d.nc",
            "period": "1979-2017",
        },
        "daily": {
            "path": "/compyfs/zhan391/acme_init/Observations/GPCP/daily",
            "template": "PRECT.daily.%(year)d.nc",
            "period": "2010-2015",
        },
    },
    "GPM": {
        "daily": {
            "path": "/compyfs/zhan391/acme_init/Observations/IMERG/daily",
            "template": "PRECT.daily.%(year)d.nc",
            "period": "2001-2020",
        }
    },
}

OBS_REGISTRY = build_obs_registry(RAW_OBS_GROUP)
