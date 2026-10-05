"""Per-region settings for the climo-based bias and budget notebooks.

Settings come from ``config/regions.json`` (``defaults`` merged with the
region's entry); region bounds come from ``exp_info.REGION_CATALOG``::

    RS = region_settings("Greenland")      # sfc_mask_kind, refexp, ...
    RUN_CATALOG, VARIABLE_CATALOG = region_catalogs("Antarctic")
"""

from __future__ import annotations

import copy
import json
from dataclasses import replace
from pathlib import Path
from typing import Any, Dict, Tuple

import exp_info

CONFIG_FILE = Path(__file__).resolve().parents[1] / "config" / "regions.json"

with CONFIG_FILE.open() as _f:
    _CONFIG = json.load(_f)

REGIONS = tuple(_CONFIG["regions"])


def region_settings(region: str) -> Dict[str, Any]:
    """Return the merged settings for ``region`` plus its ``lat_range``/``lon_range``."""
    if region not in _CONFIG["regions"]:
        raise KeyError(f"Unknown region {region!r}; configured regions: {list(REGIONS)} "
                       f"({CONFIG_FILE})")
    if region not in exp_info.REGION_CATALOG:
        raise KeyError(f"Region {region!r} has no bounds in exp_info.REGION_CATALOG")
    settings = copy.deepcopy(_CONFIG["defaults"])
    settings.update(copy.deepcopy(_CONFIG["regions"][region]))
    bounds = exp_info.REGION_CATALOG[region]
    settings.update(region=region, lat_range=bounds.lat_range, lon_range=bounds.lon_range)
    for key in ("refexp", "sfc_mask_ref"):
        if settings[key] not in exp_info.RUN_CATALOG:
            raise KeyError(f"{region}: {key}={settings[key]!r} is not in exp_info.RUN_CATALOG")
    return settings


def region_catalogs(region: str) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Return (RUN_CATALOG, VARIABLE_CATALOG) with the region's exclusions and overrides."""
    settings = region_settings(region)
    runs = {k: v for k, v in exp_info.RUN_CATALOG.items() if k not in settings["exclude_runs"]}
    variables = dict(exp_info.VARIABLE_CATALOG)
    for name, fields in settings["variable_overrides"].items():
        variables[name] = replace(variables[name], **fields)
    return runs, variables


__all__ = ["REGIONS", "region_settings", "region_catalogs"]
