"""Configuration defaults for initial-land diagnostic notebooks and scripts."""

from __future__ import annotations

import os
from pathlib import Path

from configs.output_paths import OUTPUT_DATA_ROOT, OUTPUT_FIGURE_ROOT


DEFAULT_DATA_ROOT = Path(os.environ.get("E3SM_DART_INITIAL_LAND_DATA_ROOT", "/compyfs/zhan391/v3_dart_cda_scratch"))
DEFAULT_OBS_ROOT = Path(os.environ.get("E3SM_DART_OBS_ROOT", "/compyfs/zhan391/v3_dart_cda_scratch/reference"))
DEFAULT_FIGURE_DIR = OUTPUT_FIGURE_ROOT / "analysis_lnd_init"
DEFAULT_DIAG_DIR = OUTPUT_DATA_ROOT / "analysis_lnd_init"
DEFAULT_REGRID_MAP_DIR = Path(os.environ.get("E3SM_DART_REGRID_MAP_DIR", "/compyfs/zhan391/v3_dart_cda_scratch/reference/regrid_maps"))
DEFAULT_LAND_REFERENCE_DIR = Path(os.environ.get("E3SM_DART_LAND_REFERENCE_DIR", "/compyfs/zhan391/v3_dart_cda_scratch/reference/lnd_sea_mask"))
DEFAULT_LANDMASK_FILE = DEFAULT_LAND_REFERENCE_DIR / "landmask_1x1.nc"
DEFAULT_SOIL_LAYER_FILE = DEFAULT_LAND_REFERENCE_DIR / "dzsoi_elm.nc"

DEFAULT_MODEL_EXPERIMENTS = {
    "CTRLEN10": {
        "path": DEFAULT_DATA_ROOT / "CTRLEN10_15day_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy/archive/post/lnd/180x360_aave/ts/daily",
        "template": "H2OSOI.%(ensemble)s.%(year)s.nc",
        "frequency": "daily",
        "nens": 10,
    },
    "CAPTEN10": {
        "path": DEFAULT_DATA_ROOT / "CAPTEN10_15day_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy/archive/post/lnd/180x360_aave/ts/daily",
        "template": "H2OSOI.%(ensemble)s.%(year)s.nc",
        "frequency": "daily",
        "nens": 10,
    },
    "DARTEN20": {
        "path": DEFAULT_DATA_ROOT / "DARTEN20_15day_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy/archive/post/lnd/180x360_aave/ts/daily",
        "template": "H2OSOI.%(ensemble)s.%(year)s.nc",
        "frequency": "daily",
        "nens": 20,
    },
    "DARTEN40": {
        "path": DEFAULT_DATA_ROOT / "DARTEN40_15day_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy/archive/post/lnd/180x360_aave/ts/daily",
        "template": "H2OSOI.%(ensemble)s.%(year)s.nc",
        "frequency": "daily",
        "nens": 40,
    },
}

DEFAULT_OBSERVATIONS = {
    "CPC_SOM": {
        "path": DEFAULT_OBS_ROOT / "CPC_SOM/monthly",
        "template": "SOILWATER_10CM.daily.%(year)s.nc",
        "frequency": "monthly",
        "nens": 1,
    },
    "ESA_CCI": {
        "path": DEFAULT_OBS_ROOT / "ESA_CCI/daily",
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
    "model_keys": ["CTRLEN10", "CAPTEN10", "DARTEN20", "DARTEN40"],
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
    include_observations=True,
    include_models=True,
):
    """Return data definitions rebased beneath notebook-configured input roots."""
    data_root = Path(data_root)
    obs_root = Path(obs_root)
    data = {}
    if include_observations:
        for key, spec in DEFAULT_OBSERVATIONS.items():
            item = dict(spec)
            item["path"] = obs_root / Path(spec["path"]).relative_to(DEFAULT_OBS_ROOT)
            data[key] = item
    if include_models:
        for key, spec in DEFAULT_MODEL_EXPERIMENTS.items():
            item = dict(spec)
            item["path"] = data_root / Path(spec["path"]).relative_to(DEFAULT_DATA_ROOT)
            data[key] = item
    return stringify_paths(data)
