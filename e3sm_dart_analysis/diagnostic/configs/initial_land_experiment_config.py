"""Configuration defaults for initial-land diagnostic notebooks and scripts."""

from __future__ import annotations

import os
from pathlib import Path

from configs.output_paths import OUTPUT_DATA_ROOT, OUTPUT_FIGURE_ROOT
from configs.input_paths import OBSERVATION_ROOT, LAND_SUPPORT_DIR, REGRID_MAP_DIR, PRODUCT_DIRS
from configs.initial_land_experiment_registry import build_experiments


DEFAULT_DATA_ROOT = Path(os.environ.get("E3SM_DART_INITIAL_LAND_DATA_ROOT", "/compyfs/zhan391/v3_dart_cda_scratch"))
DEFAULT_OBS_ROOT = OBSERVATION_ROOT
DEFAULT_FIGURE_DIR = OUTPUT_FIGURE_ROOT / "analysis_lnd_init"
DEFAULT_DIAG_DIR = OUTPUT_DATA_ROOT / "analysis_lnd_init"
DEFAULT_REGRID_MAP_DIR = REGRID_MAP_DIR
DEFAULT_LAND_REFERENCE_DIR = LAND_SUPPORT_DIR
DEFAULT_LANDMASK_FILE = DEFAULT_LAND_REFERENCE_DIR / "landmask_1x1.nc"
DEFAULT_SOIL_LAYER_FILE = DEFAULT_LAND_REFERENCE_DIR / "dzsoi_elm.nc"

DA_SOIL_MOISTURE_EXPERIMENTS = {
    "S0": ("CTRL10-S0", "DART10-S0", "DART20-S0", "DART40-S0"),
    "S1": ("DART40-S1",),
}

DA_TARGET_DATES = {"S0": "2012-01-01", "S1": "2012-06-01"}

DEFAULT_OBSERVATIONS = {
    "CPC_SOM": {
        "path": DEFAULT_OBS_ROOT / "CPC_SOM/monthly",
        "template": "SOILWATER_10CM.daily.%(year)s.nc",
        "frequency": "monthly",
        "nens": 1,
    },
    "ESA_CCI": {
        "path": DEFAULT_OBS_ROOT / PRODUCT_DIRS["ESA_CCI"] / "daily",
        "template": "H2OSOI.daily.%(year)s.nc",
        "frequency": "daily",
        "nens": 1,
    },
    "GPCP": {
        "path": DEFAULT_OBS_ROOT / "GPCP/daily",
        "template": "PRECT.daily.%(year)s.nc",
        "frequency": "daily",
        "nens": 1,
    },
}

DEFAULT_SOIL_MOISTURE_PLOT = {
    "model_keys": list(DA_SOIL_MOISTURE_EXPERIMENTS["S0"]),
    "obs_key": "ESA_CCI",
    "target_date": "2012-01-01",
    "depth_cm": 5,
    "mask_land": True,
    "confidence_level": 0.05,
    "use_mass_units": True,
    "bias_levels": [-8, -6, -4, -2, -1, 0, 1, 2, 4, 6, 8],
    "spread_levels": [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5],
    "figsize": (14.0, 6.6),
    "fontz": 9,
}

RESTART_OUTPUT_DIR = DEFAULT_DIAG_DIR / "restart"

DEFAULT_RESTART_REGRID = {
    "out_agg": RESTART_OUTPUT_DIR / "agg_gridcell.nc",
    "dst_scrip": DEFAULT_REGRID_MAP_DIR / "cmip6_180x360_scrip.20181001.nc",
    "src_scrip": RESTART_OUTPUT_DIR / "src_unstruct_from_agg.scrip.nc",
    "map_file": RESTART_OUTPUT_DIR / "map_unstruct_r05_to_cmip6_180x360_aave.c251029.nc",
    "regrid_out": RESTART_OUTPUT_DIR / "elm_180x360_aave.nc",
}


def stringify_paths(mapping):
    """Return a shallow copy with Path values converted to strings for notebook code."""
    return {
        key: {inner_key: str(value) if isinstance(value, Path) else value for inner_key, value in spec.items()}
        for key, spec in mapping.items()
    }


def get_soil_moisture_data_dict(
    *,
    data_root=DEFAULT_DATA_ROOT,
    obs_root=DEFAULT_OBS_ROOT,
    season="S0",
    run_segment="da",
    include_observations=True,
    include_models=True,
):
    """Return DA land-state inputs rebased beneath notebook-configured roots.

    Model files come from the cycled DA runs, never the 15-day hindcasts. The
    remapped ELM ``h0`` file represents the final six-hour DA cycle: 1 January
    for S0 and 1 June for S1.
    """
    data_root = Path(data_root)
    obs_root = Path(obs_root)
    if run_segment != "da":
        raise ValueError("Initial-land soil-moisture diagnostics require run_segment='da'.")
    if season not in DA_SOIL_MOISTURE_EXPERIMENTS:
        raise ValueError(
            f"season must be one of {sorted(DA_SOIL_MOISTURE_EXPERIMENTS)}"
        )
    data = {}
    if include_observations:
        for key, spec in DEFAULT_OBSERVATIONS.items():
            item = dict(spec)
            item["path"] = obs_root / Path(spec["path"]).relative_to(DEFAULT_OBS_ROOT)
            data[key] = item
    if include_models:
        registry = build_experiments(str(data_root))
        for experiment in DA_SOIL_MOISTURE_EXPERIMENTS[season]:
            metadata = registry[experiment]
            run = metadata["runs"][run_segment]
            if run is None:
                raise ValueError(f"{experiment} has no {run_segment!r} run.")
            year_month = f"{run.period[:4]}-{run.period[4:6]}"
            data[experiment] = {
                "path": data_root / run.run_id,
                "template": (
                    f"%(ensemble)s/archive/post/lnd/180x360_aave/clim/"
                    f"{run.run_id}.%(ensemble)s.elm.h0.{year_month}.nc"
                ),
                "frequency": "final_da_cycle",
                "nens": metadata["nens"],
            }
    return stringify_paths(data)
