import os
import glob
import xarray as xr
from typing import Dict, List, Tuple, Optional

import numpy as np
import cftime

# === Base Experiment Path ===
base_path = "/lcrc/group/e3sm/ac.szhang/acme_scratch/e3sm_project/E3SMv3_testings"


# === Helper: Create Entry for One AMIP Experiment ===
def make_run_entry(
    name: str,
    period: str,
    *,
    path: str = base_path,
    glb: Optional[str] = None,
    atm: Optional[str] = None,
    lnd: Optional[str] = None,
    rof: Optional[str] = None,
    cmip6: Optional[str] = None,
    **kwargs,
) -> Dict:
    entry = {"name": name, "period": period, "path": path}
    if glb is not None:
        entry["glb"] = glb
    if atm is not None:
        entry["atm"] = atm
    if lnd is not None:
        entry["lnd"] = lnd
    if rof is not None:
        entry["rof"] = rof
    if cmip6 is not None:
        entry["cmip6"] = cmip6
    entry.update(kwargs)
    return entry


# === AMIP Experiments ===
run_dict: Dict[str, Dict] = {
    k: make_run_entry(
        name=v,
        period="2012-2013",
        glb="post/atm/glb/ts/monthly/1yr",
        atm="post/atm/180x360_aave",
        lnd="post/lnd/glb/ts/monthly/1yr",
        rof="post/rof/native/ts/monthly/1yr",
        # CMIP6-style variables live here on disk
        cmip6="post/atm/180x360_aave/cmip_ts",
    )
    for k, v in {
        "CTRL": "20260114.v3-LR-AMIP-CTRL1.ne30pg2_r05_IcoswISC30E3r5",
        "DARTOP": "20260114.v3-LR-AMIP-DART1.ne30pg2_r05_IcoswISC30E3r5",
        "ERA5OP": "20260114.v3-LR-AMIP-ERA51.ne30pg2_r05_IcoswISC30E3r5",
        "ERA5DL": "20260114.v3-LR-AMIP-ERA5DL.ne30pg2_r05_IcoswISC30E3r5",
        "ERA5BL": "20260114.v3-LR-AMIP-ERA5BL.ne30pg2_r05_IcoswISC30E3r5",
        "DARTBL": "20260114.v3-LR-AMIP-DARTBL.ne30pg2_r05_IcoswISC30E3r5",
        "DARTBA": "20260114.v3-LR-AMIP-DARTBA.ne30pg2_r05_IcoswISC30E3r5",
        "DARTAE": "20260114.v3-LR-AMIP-DARTAE.ne30pg2_r05_IcoswISC30E3r5",
    }.items()
}


# === Utility: Normalize Aliases ===
def resolve_aliases(var_dict: Dict[str, Dict[str, Dict]]) -> Dict[str, Dict[str, Dict]]:
    for comp_vars in var_dict.values():
        for meta in comp_vars.values():
            if "alia" in meta and not isinstance(meta["alia"], list):
                meta["alia"] = [meta["alia"]]
    return var_dict


# === Utility: Flatten (component, variable) ===
def flatten_variable_dict(var_dict: Dict[str, Dict[str, Dict]]) -> Dict[Tuple[str, str], Dict]:
    return {(comp, var): meta for comp, vars_ in var_dict.items() for var, meta in vars_.items()}


# === Component → run_dict key mapping ===
COMP_TO_RUNKEY = {
    "EAM": "atm",
    "EAM_CMIP6": "cmip6",
    "ELM": "lnd",
    "MOSART": "rof",
}


# === Path helpers ===
def get_comp_dir(run: Dict, run_key: str) -> str:
    if run_key not in run:
        raise KeyError(f"Run '{run['name']}' missing '{run_key}'. Keys={list(run.keys())}")
    return os.path.join(run["path"], run["name"], run[run_key])


def find_nc_files(
    run: Dict,
    run_key: str,
    var: str,
    *,
    is_cmip6: bool = False,
) -> List[str]:
    """
    - Native EAM files: TS_201201_201212.nc (or similar)
    - CMIP6-style files: tas_Amon_*_201201-201212.nc (or similar)
    """
    base = get_comp_dir(run, run_key)

    # Search a few common nesting patterns (harmless if they don't exist).
    search_dirs = [base]

    if is_cmip6:
        # Typical CMIP6 postproc structures vary:
        #   base/*.nc
        #   base/Amon/*.nc, base/AERmon/*.nc
        #   base/monthly/*.nc
        #   base/monthly/Amon/*.nc, base/monthly/AERmon/*.nc
        search_dirs += [
            os.path.join(base, "Amon"),
            os.path.join(base, "AERmon"),
            os.path.join(base, "monthly"),
            os.path.join(base, "monthly", "Amon"),
            os.path.join(base, "monthly", "AERmon"),
        ]

        # CMIP6 pattern: <var>_<table>_..._<YYYYMM-YYYYMM>.nc
        patterns = [f"{var}_*.nc"]
    else:
        # Native postproc pattern: VAR_201201_201212.nc etc.
        patterns = [f"{var}_*.nc", f"*{var}*.nc", f"*.{var}.*.nc"]

    files: List[str] = []
    for d in search_dirs:
        for p in patterns:
            files.extend(glob.glob(os.path.join(d, p)))

    return sorted(set(files))


def open_var_from_files(files: List[str], var: str) -> xr.DataArray:
    if not files:
        raise FileNotFoundError(f"No files found for variable '{var}'")
    ds = xr.open_mfdataset(files, combine="by_coords")
    if var not in ds:
        raise KeyError(f"Variable '{var}' not in dataset. Example vars: {list(ds.data_vars)[:30]}")
    return ds[var]


# === Convenience: find/open by (component, var) ===
def find_files_for_compvar(run: Dict, component: str, var: str) -> List[str]:
    run_key = COMP_TO_RUNKEY.get(component)
    if run_key is None:
        raise KeyError(f"Unknown component '{component}'. Known={list(COMP_TO_RUNKEY.keys())}")
    return find_nc_files(run, run_key, var, is_cmip6=(component == "EAM_CMIP6"))


def open_compvar(run: Dict, component: str, var: str) -> xr.DataArray:
    files = find_files_for_compvar(run, component, var)
    return open_var_from_files(files, var)


def set_subdir(run_key: str, subdir: str, runs: Optional[List[str]] = None) -> None:
    """
    Convenience: override a path key (atm/lnd/rof/cmip6) for selected runs or all runs.

    Example:
        set_subdir("atm", "post/atm/180x360_aave/ts/monthly/1yr")
        set_subdir("cmip6", "post/atm/180x360_aave/cmip_ts/monthly", runs=["CTRL"])
    """
    targets = runs if runs is not None else list(run_dict.keys())
    for k in targets:
        run_dict[k][run_key] = subdir


# === Variable Metadata by Component (your version + added ELM/MOSART variables) ===
def collect_variable_list(components: List[str]) -> Dict[str, Dict[str, Dict]]:
    all_dict = {
        "EAM": {
            "TS": {"unit": "K", "min": -6, "max": 6, "nlev": 13},
            "TREFHT": {"unit": "K", "min": -6, "max": 6, "nlev": 13},
            "PSL": {"unit": "hPa", "min": -6, "max": 6, "nlev": 13},
            "PRECT": {"unit": "mm day$^{-1}$", "min": -6, "max": 6, "nlev": 13},
            "PRECST": {"unit": "mm day$^{-1}$", "min": -6, "max": 6, "nlev": 13},
            "EVAP": {"unit": "mm day$^{-1}$", "min": -6, "max": 6, "nlev": 13},
            "FLUT": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "ICEFRAC": {"unit": "%", "min": -6, "max": 6, "nlev": 13},
            "RESTOA": {"unit": "W m$^{-2}$", "min": -8, "max": 8, "nlev": 13},
            "RESSURF": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "RESNET": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "PS": {"unit": "Pa", "min": -6, "max": 6, "nlev": 13},
            "U": {"unit": "m s$^{-1}$", "min": -6, "max": 6, "nlev": 13},
            "CLDLOW": {"unit": "%", "min": -6, "max": 6, "nlev": 13},
            "CLDMED": {"unit": "%", "min": -6, "max": 6, "nlev": 13},
            "CLDHGH": {"unit": "%", "min": -6, "max": 6, "nlev": 13},
            "CLDTOT": {"unit": "%", "min": -6, "max": 6, "nlev": 13},
            "QFLX": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "SHFLX": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "TAUX": {"unit": "N m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "TAUY": {"unit": "N m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "FSNT": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "FSNTOA": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "FLNT": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "FSNS": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "FLNS": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "PRECC": {"unit": "mm day$^{-1}$", "min": -6, "max": 6, "nlev": 13},
            "PRECL": {"unit": "mm day$^{-1}$", "min": -6, "max": 6, "nlev": 13},
            "PRECSC": {"unit": "mm day$^{-1}$", "min": -6, "max": 6, "nlev": 13},
            "PRECSL": {"unit": "mm day$^{-1}$", "min": -6, "max": 6, "nlev": 13},
            "LANDFRAC": {"unit": "1", "min": -6, "max": 6, "nlev": 13},
            "OCNFRAC": {"unit": "1", "min": -6, "max": 6, "nlev": 13},
        },

        "EAM_CMIP6": {
            "abs550aer": {"unit": "1", "min": -6, "max": 6, "nlev": 13},
            "od550aer": {"unit": "1", "min": -6, "max": 6, "nlev": 13},
            "clt": {"unit": "%", "min": -6, "max": 6, "nlev": 13},
            "clivi": {"unit": "kg m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "clwvi": {"unit": "kg m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "evspsbl": {"unit": "kg m$^{-2}$ s$^{-1}$", "min": -6, "max": 6, "nlev": 13},
            "hfls": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "hfss": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "huss": {"unit": "1", "min": -6, "max": 6, "nlev": 13},
            "pr": {"unit": "kg m$^{-2}$ s$^{-1}$", "min": -6, "max": 6, "nlev": 13},
            "prc": {"unit": "kg m$^{-2}$ s$^{-1}$", "min": -6, "max": 6, "nlev": 13},
            "prsn": {"unit": "kg m$^{-2}$ s$^{-1}$", "min": -6, "max": 6, "nlev": 13},
            "prw": {"unit": "kg m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "ps": {"unit": "Pa", "min": -6, "max": 6, "nlev": 13},
            "psl": {"unit": "Pa", "min": -6, "max": 6, "nlev": 13},
            "rlds": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "rldscs": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "rlus": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "rlut": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "rlutcs": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "rsds": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "rsdscs": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "rsdt": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "rsus": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "rsuscs": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "rsut": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "rsutcs": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "rtmt": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "sfcWind": {"unit": "m s$^{-1}$", "min": -6, "max": 6, "nlev": 13},
            "ta": {"unit": "K", "min": -6, "max": 6, "nlev": 13},
            "tas": {"unit": "K", "min": -6, "max": 6, "nlev": 13},
            "tasmax": {"unit": "K", "min": -6, "max": 6, "nlev": 13},
            "tasmin": {"unit": "K", "min": -6, "max": 6, "nlev": 13},
            "tauu": {"unit": "N m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "tauv": {"unit": "N m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "ts": {"unit": "K", "min": -6, "max": 6, "nlev": 13},
            "ua": {"unit": "m s$^{-1}$", "min": -6, "max": 6, "nlev": 13},
            "va": {"unit": "m s$^{-1}$", "min": -6, "max": 6, "nlev": 13},
            "zg": {"unit": "m", "min": -6, "max": 6, "nlev": 13},
        },

        "ELM": {
            # Hydrology
            "AR": {"unit": "mm day$^{-1}$", "min": -6, "max": 6, "nlev": 13},
            "QOVER": {"unit": "mm day$^{-1}$", "min": -6, "max": 6, "nlev": 13},
            "QRUNOFF": {"unit": "mm day$^{-1}$", "min": -6, "max": 6, "nlev": 13},
            "QSOIL": {"unit": "mm day$^{-1}$", "min": -6, "max": 6, "nlev": 13},
            "QINTR": {"unit": "mm day$^{-1}$", "min": -6, "max": 6, "nlev": 13},
            "QVEGE": {"unit": "mm day$^{-1}$", "min": -6, "max": 6, "nlev": 13},
            "QVEGT": {"unit": "mm day$^{-1}$", "min": -6, "max": 6, "nlev": 13},

            "SOILWATER_10CM": {"unit": "mm", "min": -6, "max": 6, "nlev": 13},
            "SOILICE": {"unit": "kg m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "SOILLIQ": {"unit": "kg m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "SOIL4C": {"unit": "gC m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "SOIL3C": {"unit": "gC m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "SOIL2C": {"unit": "gC m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "SOIL1C": {"unit": "gC m$^{-2}$", "min": -6, "max": 6, "nlev": 13},

            "H2OSNO": {"unit": "kg m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "SNOW_DEPTH": {"unit": "m", "min": -6, "max": 6, "nlev": 13},

            # Fluxes / state
            "FSH": {"unit": "W m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "RH2M": {"unit": "%", "min": -6, "max": 6, "nlev": 13},
            "TSA": {"unit": "K", "min": -6, "max": 6, "nlev": 13},
            "TSOI": {"unit": "K", "min": -6, "max": 6, "nlev": 13},

            # LAI
            "LAISHA": {"unit": "m$^{2}$ m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "LAISUN": {"unit": "m$^{2}$ m$^{-2}$", "min": -6, "max": 6, "nlev": 13},

            # Carbon fluxes (leave units generic-ish; file attrs are truth)
            "GPP": {"unit": "unknown", "min": -6, "max": 6, "nlev": 13},
            "NBP": {"unit": "unknown", "min": -6, "max": 6, "nlev": 13},
            "HR": {"unit": "unknown", "min": -6, "max": 6, "nlev": 13},

            # Carbon pools
            "TOTVEGC": {"unit": "gC m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "TOTLITC": {"unit": "gC m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "CWDC": {"unit": "gC m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
            "WOOD_HARVESTC": {"unit": "gC m$^{-2}$", "min": -6, "max": 6, "nlev": 13},
        },

        "MOSART": {
            "RIVER_DISCHARGE_OVER_LAND_LIQ": {"unit": "m$^{3}$ s$^{-1}$", "min": -6, "max": 6, "nlev": 13},
            "RUNOFF": {"unit": "mm day$^{-1}$", "min": -6, "max": 6, "nlev": 13},
        },
    }

    var_dict: Dict[str, Dict[str, Dict]] = {}
    for comp in components:
        if comp not in all_dict:
            raise ValueError(f"Component '{comp}' is not recognized. Must be one of {list(all_dict.keys())}")
        var_dict[comp] = all_dict[comp]
    return var_dict

def set_monthly_time_to_midmonth_cftime(da, time_name="time", day=16, hour=0, minute=0, second=0, calendar_cls=None):
    """
    Replace da[time_name] with cftime dates at YYYY-MM-<day> <hour>:<minute>:<second>.
    calendar_cls: e.g., cftime.DatetimeNoLeap. If None, infer from existing time values;
                  if existing is numpy datetime64, default to DatetimeNoLeap.
    """
    tvals = da[time_name].values

    # infer calendar class
    if calendar_cls is None:
        t0 = tvals[0]
        if isinstance(t0, cftime.datetime):
            calendar_cls = t0.__class__
        else:
            # ERA5 (numpy datetime64) -> choose model-like calendar explicitly
            calendar_cls = cftime.DatetimeNoLeap

    # get year/month robustly for either cftime or numpy datetime64
    years = []
    months = []
    for tt in tvals:
        if isinstance(tt, cftime.datetime):
            years.append(tt.year); months.append(tt.month)
        else:
            # numpy datetime64 -> convert to YYYY-MM string without pandas
            s = np.datetime_as_string(tt, unit="D")  # "YYYY-MM-DD"
            years.append(int(s[0:4])); months.append(int(s[5:7]))

    new_t = np.array(
        [calendar_cls(y, m, day, hour, minute, second) for y, m in zip(years, months)],
        dtype=object,
    )

    da = da.assign_coords({time_name: (time_name, new_t)})
    # optional but sometimes helpful metadata
    da[time_name].attrs["calendar"] = getattr(calendar_cls, "__name__", "cftime")
    return da

__all__ = [
    "base_path",
    "run_dict",
    "collect_variable_list",
    "resolve_aliases",
    "flatten_variable_dict",
    "COMP_TO_RUNKEY",
    "get_comp_dir",
    "find_nc_files",
    "open_var_from_files",
    "find_files_for_compvar",
    "open_compvar",
    "set_subdir",
    "set_monthly_time_to_midmonth_cftime",
    "standardize_spacetime"
]
