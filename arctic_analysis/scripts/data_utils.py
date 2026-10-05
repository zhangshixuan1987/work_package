from __future__ import annotations

import os
from typing import Dict, Tuple, Iterable, Optional

import numpy as np
import pandas as pd
import xarray as xr
from xcdat import open_dataset as xcdat_open

# ----------------------------
# Experiment / run metadata
# ----------------------------

RUN_DICT = {
    "v3.LR.historical.en00": {"name": "v3.LR.historical.0051", "period": "185001-202412"},
    "v3.LR.historical.en01": {"name": "v3.LR.historical.0091", "period": "185001-202412"},
    "v3.LR.historical.en02": {"name": "v3.LR.historical.0101", "period": "185001-202412"},
    "v3.LR.historical.en03": {"name": "v3.LR.historical.0111", "period": "185001-202412"},
    "v3.LR.historical.en04": {"name": "v3.LR.historical.0121", "period": "185001-202412"},
    "v3.LR.historical.en05": {"name": "v3.LR.historical.0131", "period": "185001-202412"},
    "v3.LR.historical.en06": {"name": "v3.LR.historical.0141", "period": "185001-202412"},
    "v3.LR.historical.en07": {"name": "v3.LR.historical.0151", "period": "185001-202412"},
    "v3.LR.historical.en08": {"name": "v3.LR.historical.0161", "period": "185001-202412"},
    "v3.LR.historical.en09": {"name": "v3.LR.historical.0171", "period": "185001-202412"},
    "v3.LR.historical.en10": {"name": "v3.LR.historical.0181", "period": "185001-202412"},
    "v3.LR.historical.en11": {"name": "v3.LR.historical.0191", "period": "185001-202412"},
    "v3.LR.historical.en12": {"name": "v3.LR.historical.0201", "period": "185001-202412"},
    "v3.LR.historical.en13": {"name": "v3.LR.historical.0211", "period": "185001-202412"},
    "v3.LR.historical.en14": {"name": "v3.LR.historical.0221", "period": "185001-202412"},
    "v3.LR.historical.en15": {"name": "v3.LR.historical.0231", "period": "185001-202412"},
    "v3.LR.historical.en16": {"name": "v3.LR.historical.0241", "period": "185001-202412"},
    "v3.LR.historical.en17": {"name": "v3.LR.historical.0251", "period": "185001-202412"},
    "v3.LR.historical.en18": {"name": "v3.LR.historical.0261", "period": "185001-202412"},
    "v3.LR.historical.en19": {"name": "v3.LR.historical.0271", "period": "185001-202412"},
    "v3.LR.historical.en20": {"name": "v3.LR.historical.0281", "period": "185001-202412"},
    "v3.LR.historical.en21": {"name": "v3.LR.historical.0291", "period": "185001-202412"},
    "v3.LR.historical.en22": {"name": "v3.LR.historical.0301", "period": "185001-202412"},
    "v3.LR.historical.en23": {"name": "v3.LR.historical.0311", "period": "185001-202412"},
    "v3.LR.historical.en24": {"name": "v3.LR.historical.0321", "period": "185001-202412"},
}

def get_run_meta(run: str) -> dict:
    """
    Return metadata for a given run (ensemble member).
    Example keys: 'name', 'period'.
    Raises KeyError if run not found.
    """
    if run not in RUN_DICT:
        raise KeyError(f"Unknown run '{run}'. Options: {list(RUN_DICT.keys())[:5]} ...")
    return RUN_DICT[run]

def list_runs(prefix: str = "") -> list[str]:
    """
    List all runs, optionally filtering by a prefix string.
    """
    return [k for k in RUN_DICT if k.startswith(prefix)]

# ----------------------------
# Reference variable metadata
# ----------------------------

REFERENCE_VAR_DICT = {
    "SST"          : {"name": "Sea Surface Temperature", "unit": r"$^\circ$C",   "min": -6,   "max": 6,   "nlev": 11},
    "SICONC"       : {"name": "Sea Ice Concentration",   "unit": "%",            "min": -6,   "max": 6,   "nlev": 11},
    "SIMASS"       : {"name": "Sea Ice Mass",            "unit": r"kg m$^{-2}$", "min": -6,   "max": 6,   "nlev": 11},
    "SITHICK"      : {"name": "Sea Ice Thickness",       "unit": "mm",           "min": -6,   "max": 6,   "nlev": 11},
    "iceArea"      : {"name": "Ice Area",                "unit": "m$^2$",        "min": -6,   "max": 6,   "nlev": 11},
    "iceVolume"    : {"name": "Ice Volume",              "unit": "m$^3$",        "min": -0.5, "max": 0.5, "nlev": 11},
    "iceThickness" : {"name": "Ice Thickness",           "unit": "m",            "min": -0.5, "max": 0.5, "nlev": 11},
    "TS"           : {"name": "Skin Temperature",        "unit": "K",            "min": -6,   "max": 6,   "nlev": 11},
    "TURFLX"       : {"name": "Turbulent Flux(SH + LH)", "unit": r"W m$^{-2}$",  "min": -20,  "max": 20,  "nlev": 11},
    "LHFLX"        : {"name": "Latent Heat Flux",        "unit": r"W m$^{-2}$",  "min": -20,  "max": 20,  "nlev": 11},
    "SHFLX"        : {"name": "Sensible Heat Flux",      "unit": r"W m$^{-2}$",  "min": -20,  "max": 20,  "nlev": 11},
    "RADNET"       : {"name": "Net Radiative Flux",      "unit": r"W m$^{-2}$",  "min": -20,  "max": 20,  "nlev": 11},
    "FSNS"         : {"name": "Net SW Flux (surf)",      "unit": r"W m$^{-2}$",  "min": -20,  "max": 20,  "nlev": 11},
    "FSDS"         : {"name": "Downward SW Flux (surf)", "unit": r"W m$^{-2}$",  "min": -20,  "max": 20,  "nlev": 11},
    "FSUS"         : {"name": "Upward SW Flux (surf)",   "unit": r"W m$^{-2}$",  "min": -20,  "max": 20,  "nlev": 11},
    "FSNSC"        : {"name": "Net SW Flux(surf,clr)",   "unit": r"W m$^{-2}$",  "min": -20,  "max": 20,  "nlev": 11},
    "FSDSC"        : {"name": "Downward SW Flux (clr)",  "unit": r"W m$^{-2}$",  "min": -20,  "max": 20,  "nlev": 11},
    "FSUSC"        : {"name": "Upward SW Flux (clr)",    "unit": r"W m$^{-2}$",  "min": -20,  "max": 20,  "nlev": 11},
    "FLDS"         : {"name": "Downward LW Flux (surf)", "unit": r"W m$^{-2}$",  "min": -20,  "max": 20,  "nlev": 11},
    "FLNS"         : {"name": "Net LW Flux (surf)",      "unit": r"W m$^{-2}$",  "min": -20,  "max": 20,  "nlev": 11},
    "FLUS"         : {"name": "Upward LW Flux (surf)",   "unit": r"W m$^{-2}$",  "min": -20,  "max": 20,  "nlev": 11},
    "FLUSC"        : {"name": "Upward LW Flux (clr)",    "unit": r"W m$^{-2}$",  "min": -20,  "max": 20,  "nlev": 11},
    "FLDSC"        : {"name": "Downward LW Flux (clr)",  "unit": r"W m$^{-2}$",  "min": -20,  "max": 20,  "nlev": 11},
    "FLNSC"        : {"name": "Net LW Flux (clr)",       "unit": r"W m$^{-2}$",  "min": -20,  "max": 20,  "nlev": 11},
    "SWCF_SRF"     : {"name": "SW Cloud Radiative Effect","unit": r"W m$^{-2}$", "min": -20,  "max": 20,  "nlev": 11},
    "LWCF_SRF"     : {"name": "LW Cloud Radiative Effect","unit": r"W m$^{-2}$", "min": -20,  "max": 20,  "nlev": 11},
    "ALBEDOC_SRF"  : {"name": "Surface Albedo (clr)",    "unit": "%",            "min": -6,   "max": 6,   "nlev": 13},
    "ALBEDO_SRF"   : {"name": "Surface Albedo",          "unit": "%",            "min": -6,   "max": 6,   "nlev": 13},
}


def get_var_meta(var: str, *,
                 default_unit: str = "",
                 default_min: float = -1.0,
                 default_max: float = 1.0,
                 default_nlev: int = 11) -> dict:
    """
    Return plotting metadata for `var`.
    Falls back to provided defaults if `var` not in REFERENCE_VAR_DICT.
    """
    meta = REFERENCE_VAR_DICT.get(var, None)
    if meta is None:
        return {"unit": default_unit, "min": default_min, "max": default_max, "nlev": default_nlev}
    # basic validation & copy
    unit = meta.get("unit", default_unit)
    vmin = float(meta.get("min", default_min))
    vmax = float(meta.get("max", default_max))
    nlev = int(meta.get("nlev", default_nlev))
    if nlev < 2:
        nlev = max(2, nlev)
    if vmax <= vmin:
        vmax = vmin + 1.0
    return {"unit": unit, "min": vmin, "max": vmax, "nlev": nlev}

def get_levels(var: str, *, symmetric: bool = False) -> np.ndarray:
    """
    Build contour/colormap levels for `var` using its metadata.
    If symmetric=True, forces levels to be symmetric about zero.
    """
    meta = get_var_meta(var)
    vmin, vmax, nlev = meta["min"], meta["max"], meta["nlev"]
    if symmetric:
        bound = max(abs(vmin), abs(vmax))
        vmin, vmax = -bound, bound
    return np.linspace(vmin, vmax, nlev)

def format_unit(var: str) -> str:
    """Return the LaTeX-ready unit string for `var`."""
    return get_var_meta(var)["unit"]


# ----------------------------
# Regions & Longitude Helpers
# ----------------------------

def define_region(regnam: str = "global") -> Tuple[Tuple[float, float], Tuple[float, float]]:
    """
    Return (lat_min, lat_max), (lon_min, lon_max) for a named region.
    Longitudes are in the [-180, 180] convention.
    """
    reg_dict: Dict[str, Tuple[Tuple[float, float], Tuple[float, float]]] = {
        "global":    ((-90, 90),  (-180, 180)),
        "Atlantic":  ((0, 60),    (265 - 360, 345 - 360)),  # convert to [-180,180]: [ -95, -15 ]
        "Arctic":    ((50, 90),   (-180, 180)),
        "CONUS":     ((25, 50),   (235 - 360, 295 - 360)),  # [-125, -65]
        "Antarctic": ((-90, -50), (-180, 180)),
        "PolarN":    ((50, 90),   (-180, 180)),
        "Greenland": ((55, 90),   (-75, -10)),
    }
    if regnam not in reg_dict:
        raise KeyError(f"Unknown region '{regnam}'. Options: {list(reg_dict)}")
    return reg_dict[regnam]

def list_regions() -> list[str]:
    """
    Return the list of supported region names for `define_region`.

    Note:
        This mirrors the keys inside `define_region`. If you add a new region
        there, please update this list accordingly.
    """
    return [
        "global",
        "Atlantic",
        "Arctic",
        "CONUS",
        "Antarctic",
        "PolarN",
        "Greenland",
    ]


def normalize_longitude(ds, lon_name="lon"):
    """
    Convert 0..360-style longitudes (e.g., 0.5, 1.5, ..., 359.5) to -180..180.
    Keeps -180..180 grids unchanged. Sorts by longitude afterward.
    """
    lons = ds[lon_name]
    # Detect a 0..360 grid (allowing half-degree centers like 359.5)
    if np.all(np.isfinite(lons)) and float(lons.min()) >= 0.0 and float(lons.max()) <= 360.0 - 1e-6:
        new = (lons + 180.0) % 360.0 - 180.0
        # If exactly -180 sneaks in on other grids, map to +180 to avoid a duplicate with existing +180
        new = xr.where(np.isclose(new, -180.0), 180.0, new)
        ds = ds.assign_coords({lon_name: new}).sortby(lon_name)
    return ds

# ----------------------------
# Time Utilities
# ----------------------------

def _ensure_valid_time(
    ds: xr.Dataset,
    period: str,
    frequency: str = "monthly",
    time_name: str = "time",
) -> xr.Dataset:
    """
    Ensure `ds[time_name]` is a proper datetime64 index. If missing/invalid,
    reconstruct from a period string like '185001-202412'.

    frequency:
      - "monthly" -> month start dates "MS"
    """
    if time_name in ds.coords:
        try:
            t = ds[time_name]
            if hasattr(t, "dt") and np.issubdtype(t.dtype, np.datetime64):
                return ds  # already valid
        except Exception:
            pass  # fall through to reconstruct

    if frequency != "monthly":
        raise NotImplementedError("Only 'monthly' reconstruction is implemented.")

    try:
        start_str, end_str = period.split("-")
        start = pd.to_datetime(f"{start_str[:4]}-{start_str[4:]}-01")
        end   = pd.to_datetime(f"{end_str[:4]}-{end_str[4:]}-01")
    except Exception as e:
        raise ValueError(f"Failed to parse period '{period}'. Expected 'YYYYMM-YYYYMM'.") from e

    new_time = pd.date_range(start=start, end=end, freq="MS")

    # If dataset already has a time dim length, ensure consistency
    nt = ds.dims.get(time_name, None)
    if nt is not None and nt != len(new_time):
        raise ValueError(
            f"Expected {len(new_time)} time steps from period {period}, but dataset has {nt} along '{time_name}'."
        )

    print(f"[INFO] Reconstructing '{time_name}' using monthly steps {start.date()} → {end.date()}")
    return ds.assign_coords({time_name: new_time})

# ----------------------------
# Mask Generators
# ----------------------------

def _add_bounds_if_missing(ds: xr.Dataset) -> xr.Dataset:
    """
    Use xcdat's bounds accessor if available to add missing bounds. No-op if not available.
    """
    try:
        if hasattr(ds, "bounds"):
            return ds.bounds.add_missing_bounds()
    except Exception:
        pass
    return ds


def _subset_region(
    ds: xr.Dataset | xr.DataArray,
    lat_slice: Tuple[float, float],
    lon_slice: Tuple[float, float],
    lat_name: str = "lat",
    lon_name: str = "lon",
) -> xr.Dataset | xr.DataArray:
    """Subset by simple lat/lon slices (expects lon in [-180, 180])."""
    if lat_name in ds.coords:
        ds = ds.sel({lat_name: slice(lat_slice[0], lat_slice[1])})
    if lon_name in ds.coords:
        ds = ds.sel({lon_name: slice(lon_slice[0], lon_slice[1])})
    return ds


def gen_land_mask(
    data_dir: str,
    exp: str,
    group: str,
    period: str,
    var: str,
    lat_slice: Tuple[float, float],
    lon_slice: Tuple[float, float],
    *,
    lon_name: str = "lon",
    lat_name: str = "lat",
    positive_is_land: bool = True,
    land_threshold: float = 0.0,
    inclusive: bool = False,
) -> xr.DataArray:
    """
    Generate a binary land mask aligned with a reference field:
    1 where `refds[var] > land_threshold` (if positive_is_land), else 0.

    Expected file layout: {data_dir}/{group}/{exp}.{var}.{period}.nc
    (Adjust as needed if your layout differs.)
    """
    refpath = os.path.join(data_dir, group, f"{exp}.{var}.{period}.nc")
    if not os.path.exists(refpath):
        raise FileNotFoundError(f"Reference file not found: {refpath}")

    refds = xcdat_open(refpath, decode_times=True)
    refds = _ensure_valid_time(refds, period, frequency="monthly", time_name="time")
    refds = normalize_longitude(refds, lon_name=lon_name)
    refds = _add_bounds_if_missing(refds)
    refds = _subset_region(refds, lat_slice, lon_slice, lat_name=lat_name, lon_name=lon_name)

    if var not in refds:
        raise KeyError(f"Variable '{var}' not found in {refpath}.")

    da = refds[var]
    # robust comparison with NaNs
    if positive_is_land:
        cmp = da >= land_threshold if inclusive else da > land_threshold
        mask = xr.where(cmp, 1, 0)
    else:
        cmp = da <= land_threshold if inclusive else da < land_threshold
        mask = xr.where(cmp, 1, 0)

    mask = mask.where(np.isfinite(da), 0)  # NaN -> 0
    mask = mask.astype(np.int8).rename("land_mask")
    mask.attrs.update({"long_name": "Binary land mask", "source_file": os.path.basename(refpath)})
    return mask


def gen_ice_mask(
    data_dir: str,
    exp: str,
    group: str,
    period: str,
    var: str,
    lat_slice: Tuple[float, float],
    lon_slice: Tuple[float, float],
    *,
    lon_name: str = "lon",
    lat_name: str = "lat",
    ice_threshold: float = 15.0,
    inclusive: bool = False,
) -> xr.DataArray:
    """
    Generate a binary sea-ice mask (1 where ice concentration exceeds threshold, else 0).

    Expected file layout: {data_dir}/{group}/{exp}.{var}.{period}.nc
    var = ice concentration (%), default threshold = 15%
    """
    refpath = os.path.join(data_dir, group, f"{exp}.{var}.{period}.nc")
    if not os.path.exists(refpath):
        raise FileNotFoundError(f"Reference file not found: {refpath}")

    refds = xcdat_open(refpath, decode_times=True)
    refds = _ensure_valid_time(refds, period, frequency="monthly", time_name="time")
    refds = normalize_longitude(refds, lon_name=lon_name)
    refds = _add_bounds_if_missing(refds)
    refds = _subset_region(refds, lat_slice, lon_slice, lat_name=lat_name, lon_name=lon_name)

    if var not in refds:
        raise KeyError(f"Variable '{var}' not found in {refpath}.")

    da = refds[var]
    mask = xr.where(da >= ice_threshold, 1, 0) if inclusive else xr.where(da > ice_threshold, 1, 0)
    mask = mask.where(np.isfinite(da), 0)
    mask = mask.astype(np.int8).rename("ice_mask")
    mask.attrs.update(
        {
            "long_name": f"Binary sea-ice mask (> {ice_threshold}%)",
            "source_file": os.path.basename(refpath),
        }
    )
    return mask

__all__ = [
  "RUN_DICT", "get_run_meta", "list_runs",
  "REFERENCE_VAR_DICT", "get_var_meta", "get_levels", "format_unit",
  "define_region", "list_regions", "normalize_longitude",
  "gen_land_mask", "gen_ice_mask",
]