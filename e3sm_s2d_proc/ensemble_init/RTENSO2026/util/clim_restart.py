"""Select a BruteForce year for the RTENSO2026 ELM, MOSART, and MPAS-Seaice restarts.

For one initialization date (for example Feb 01), the restarts for the same
calendar date are read for every candidate year of an existing
initial-condition archive, and one year is selected. Its restart files are
copied unchanged (renamed to the target date), so each component state is a
real, self-consistent model state. Two branches choose the year differently:

``clim``  closest to the climatology of the candidate years;
``obs``   closest to the observed (ERA5) anomaly on the target date.

Distance metric
---------------
For key field ``f``, the candidate-year values ``x_y`` at each point are turned
into standardized anomalies::

    z_y,i = (x_y,i - xbar_i) / S,   S**2 = sum_i w_i s_i**2 / sum_i w_i

with ``xbar_i`` and ``s_i**2`` the mean and interannual variance at point ``i``
and ``w_i`` an area weight (``S`` is the area-weighted RMS interannual spread).
The distance of year ``y`` from a target ``t`` is::

    d_f(y) = sqrt( sum_i w_i (z_y,i - t_i)**2 / sum_i w_i )

* ``clim``: ``t = 0``. ``d_f`` is the area-weighted RMS anomaly relative to the
  RMS spread; its mean square over the candidate years is exactly 1.
* ``obs``: ``t`` is the ERA5 standardized anomaly on the target date,
  ``(o_T - obar) / S_obs``, with the ERA5 mean and spread taken over the same
  candidate years and calendar date. Model and ERA5 anomalies are compared, so
  model-ERA5 biases and unit or layer-depth differences cancel.

Points missing in any year (model or ERA5) get zero weight. The component
distance ``D_c`` is the RMS of its field distances, and the joint distance is
``D = sqrt(sum_c a_c D_c**2 / sum_c a_c)`` with component weights ``a_c``
(default 1). ``mode="joint"`` takes the year with the smallest ``D`` for all
components (one coupled state); ``mode="per_component"`` takes each
component's smallest ``D_c``.

The river (MOSART) is not part of the distance: its storage is driven by land
runoff, so it would count land twice. Its restart is still written, always
from the land (ELM) year, so river and land stay consistent.

Key fields (``KEY_FIELDS``) and their ERA5 counterparts (``OBS_FIELDS``):

=========  ==========================  =================================  =====================
Component  Field                       Model definition                   ERA5 (``obs`` branch)
=========  ==========================  =================================  =====================
mpassi     NH/SH ice concentration     sum of iceAreaCategory             ci
           NH/SH ice thickness         sum of iceVolumeCategory (m)       --
           NH/SH snow on ice           sum of snowVolumeCategory (m)      --
elm        snow water equivalent       H2OSNO (mm)                        sd
           surface soil water          liquid + ice, soil layers 1-3      swvl1 (0-7 cm)
           subsurface soil water       liquid + ice, soil layers 4+       swvl2-4, depth-weighted (7-289 cm)
           surface soil temperature    T_SOISNO, soil layer 1 (K)         stl1
           subsurface soil temperature T_SOISNO, soil layers 4-10,        stl2-4, depth-weighted (7-289 cm)
                                       depth-weighted to 2.89 m (K)
=========  ==========================  =================================  =====================

Weights: sea ice ``areaCell`` (split by hemisphere); land
``cols1d_wtxy * cos(lat)`` on active columns (summed to 0.5-degree grid cells
for the ``obs`` comparison).
ERA5 values are taken at the nearest 0.25-degree grid point at 00Z.

Other observation sets are recognized from the layout of the observation root
(``obs_source``): ERA5-Land (``<root>/1x1/3hourly/ERA5_Land_<date>_00UTC.nc``,
same variable names as ERA5) for land and the NOAA/NSIDC sea-ice CDR
(``<root>/1x1/daily/ICEFRAC.daily.<year>.nc``, read as ``ci``) for sea ice,
both on a 1-degree grid and sampled at the nearest cell.

Source files are only ever opened read-only, and output may not be written
inside the source archive.
"""

from __future__ import annotations

import glob
import json
import os
import shutil
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import netCDF4 as nc
import numpy as np
import pandas as pd

SPVAL_THRESHOLD = 1.0e35
ALL_COMPONENTS = ("elm", "mosart", "mpassi")
BRANCHES = ("clim", "obs")

COMPONENTS = {
    "elm": {"stream": "elm.r", "sep": "-", "rpointer": "rpointer.lnd", "pointer": "file"},
    "mosart": {"stream": "mosart.r", "sep": "-", "rpointer": "rpointer.rof", "pointer": "file"},
    "mpassi": {"stream": "mpassi.rst", "sep": "_", "rpointer": "rpointer.ice", "pointer": "date"},
}

# Field name -> (model definition, domain and weight).
KEY_FIELDS = {
    "mpassi": {
        "NH ice concentration": ("sum of iceAreaCategory", "lat > 0, areaCell"),
        "SH ice concentration": ("sum of iceAreaCategory", "lat < 0, areaCell"),
        "NH ice thickness": ("sum of iceVolumeCategory (m)", "lat > 0, areaCell"),
        "SH ice thickness": ("sum of iceVolumeCategory (m)", "lat < 0, areaCell"),
        "NH snow on ice": ("sum of snowVolumeCategory (m)", "lat > 0, areaCell"),
        "SH snow on ice": ("sum of snowVolumeCategory (m)", "lat < 0, areaCell"),
    },
    "elm": {
        "snow water equivalent": ("H2OSNO (mm)", "active columns, cols1d_wtxy*cos(lat)"),
        "surface soil water": ("H2OSOI_LIQ + H2OSOI_ICE, soil layers 1-3 (kg m-2)",
                               "active columns, cols1d_wtxy*cos(lat)"),
        "subsurface soil water": ("H2OSOI_LIQ + H2OSOI_ICE, soil layers 4 and below (kg m-2)",
                                  "active columns, cols1d_wtxy*cos(lat)"),
        "surface soil temperature": ("T_SOISNO, soil layer 1 (K)", "active columns, cols1d_wtxy*cos(lat)"),
        "subsurface soil temperature": ("T_SOISNO, soil layers 4-10, depth-weighted over 0.09-2.89 m (K)",
                                        "active columns, cols1d_wtxy*cos(lat)"),
    },
}

# Components whose restart follows another component's selected year instead
# of having a distance of their own: the river follows the land.
FOLLOWS = {"mosart": "elm"}

# ERA5 counterpart of each observable key field.
OBS_FIELDS = {
    "mpassi": {"NH ice concentration": "ci", "SH ice concentration": "ci"},
    "elm": {
        "snow water equivalent": "sd",
        "surface soil water": "swvl1",
        "subsurface soil water": "swvl2-4",
        "surface soil temperature": "stl1",
        "subsurface soil temperature": "stl2-4",
    },
}

ERA5_ROOT = "/global/cfs/projectdirs/m3522/cmip6/ERA5/e5.oper.an.sfc"
ERA5_CODES = {
    "ci": "128_031_ci", "sd": "128_141_sd", "stl1": "128_139_stl1",
    "stl2": "128_170_stl2", "stl3": "128_183_stl3", "stl4": "128_236_stl4",
    "swvl1": "128_039_swvl1", "swvl2": "128_040_swvl2",
    "swvl3": "128_041_swvl3", "swvl4": "128_042_swvl4",
}
ERA5_SWVL_THICKNESS = {"swvl1": 0.07, "swvl2": 0.21, "swvl3": 0.72, "swvl4": 1.89}   # m
ERA5_STL_THICKNESS = {"stl2": 0.21, "stl3": 0.72, "stl4": 1.89}                       # m
# Depth-weighted ERA5 layer groups: name -> {variable: layer thickness (m)}.
ERA5_LAYER_GROUPS = {
    "swvl2-4": {v: ERA5_SWVL_THICKNESS[v] for v in ("swvl2", "swvl3", "swvl4")},
    "stl2-4": ERA5_STL_THICKNESS,
}

# ELM standard soil layers (nlevgrnd = 15): node depths z_j = 0.025 (exp(0.5 (j - 0.5)) - 1) m and
# interfaces halfway between nodes. The restart does not store them.
_ELM_Z = 0.025 * (np.exp(0.5 * (np.arange(1, 16) - 0.5)) - 1)
ELM_SOIL_INTERFACES = np.concatenate([[0.0], 0.5 * (_ELM_Z[:-1] + _ELM_Z[1:]), [1.5 * _ELM_Z[-1] - 0.5 * _ELM_Z[-2]]])
# Subsurface soil temperature: soil layers 4-10, each weighted by its thickness above 2.89 m (the
# bottom of ERA5 stl4), so the model layers cover the depth range of ERA5 stl2-4.
ELM_TSUB_WEIGHTS = np.clip(np.minimum(ELM_SOIL_INTERFACES[4:11], 2.89) - ELM_SOIL_INTERFACES[3:10], 0.0, None)

# Observation sets other than ERA5 (0.25 degree), recognized from the root
# layout. Both are on a 1-degree grid, lat -89.5..89.5, lon 0.5..359.5.
ERA5_LAND_FILE = "1x1/3hourly/ERA5_Land_{date}_00UTC.nc"     # one 00Z snapshot per file
NSIDC_FILE = "1x1/daily/ICEFRAC.daily.{year:04d}.nc"         # daily means, one file per year
NSIDC_VARS = {"ci": "ICEFRAC"}


def metric_table(obs_fields: Any | None = None) -> pd.DataFrame:
    """Key fields, their weights, and ERA5 counterparts as a table."""
    fmap = _resolve_obs_fields(obs_fields)
    return pd.DataFrame(
        [{"component": c, "field": f, "model definition": d, "domain, weight": w,
          "ERA5 (obs branch)": fmap.get(c, {}).get(f, "--")}
         for c, fields in KEY_FIELDS.items() for f, (d, w) in fields.items()]
    )


# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------

def restart_name(prefix: str, component: str, year: int, mmdd: str) -> str:
    info = COMPONENTS[component]
    return f"{prefix}.{info['stream']}.{year:04d}-{mmdd}{info['sep']}00000.nc"


def source_file(src_root, prefix: str, component: str, year: int, mmdd: str) -> Path:
    return Path(src_root) / f"{year:04d}-{mmdd}-00000" / restart_name(prefix, component, year, mmdd)


def available_years(src_root, prefix: str, component: str, mmdd: str) -> list[int]:
    """Years whose ``<year>-<mmdd>-00000`` directory holds the component restart."""
    years = []
    for path in sorted(Path(src_root).glob(f"????-{mmdd}-00000")):
        year = int(path.name[:4])
        if source_file(src_root, prefix, component, year, mmdd).is_file():
            years.append(year)
    return years


def _resolve_obs_roots(
    obs_root: str | Path | Mapping[str, str | Path | None] | None = None,
    era5_root: str | Path | Mapping[str, str | Path | None] | None = None,
) -> dict[str, str]:
    """Resolve observation root directories for ``elm`` and ``mpassi``.

    If ``obs_root`` / ``era5_root`` is ``None`` (or a component entry in a
    mapping is ``None``), falls back automatically to the default ``ERA5_ROOT``.
    """
    spec = obs_root if obs_root is not None else era5_root
    if spec is None:
        default_str = str(Path(ERA5_ROOT).resolve())
        return {c: default_str for c in OBS_FIELDS}
    if isinstance(spec, Mapping):
        out = {}
        for c in OBS_FIELDS:
            val = spec.get(c)
            out[c] = str(Path(ERA5_ROOT if val is None else val).resolve())
        return out
    resolved = str(Path(spec).resolve())
    return {c: resolved for c in OBS_FIELDS}


def _resolve_obs_fields(
    obs_fields: Any | None = None,
) -> dict[str, dict[str, str]]:
    """Resolve observable field mappings per component, defaulting to ``OBS_FIELDS``."""
    if obs_fields is None:
        return {c: dict(v) for c, v in OBS_FIELDS.items()}
    out = {}
    for c, default_map in OBS_FIELDS.items():
        custom = obs_fields.get(c) if isinstance(obs_fields, Mapping) else None
        if custom is None:
            out[c] = dict(default_map)
        elif isinstance(custom, Mapping):
            out[c] = dict(custom)
        elif isinstance(custom, (list, tuple, set)):
            # If user provided a list of ERA5 variable names (e.g. ["sd", "swvl1", ...]),
            # match them back to the corresponding model key field names
            out[c] = {f: v for f, v in default_map.items() if v in custom}
        else:
            out[c] = dict(default_map)
    return out


def obs_source(root: str | Path | None = None) -> str:
    """Observation set under ``root``: ``"era5"``, ``"era5_land"``, or ``"nsidc"``."""
    root = Path(ERA5_ROOT if root is None else root)
    if (root / "1x1" / "3hourly").is_dir() and any((root / "1x1" / "3hourly").glob("ERA5_Land_*.nc")):
        return "era5_land"
    if (root / "1x1" / "daily").is_dir() and any((root / "1x1" / "daily").glob("ICEFRAC.daily.*.nc")):
        return "nsidc"
    return "era5"


def era5_file(var: str, year: int, month: int, era5_root: str | Path | None = None) -> Path | None:
    root = ERA5_ROOT if era5_root is None else era5_root
    hits = glob.glob(f"{Path(root)}/{year:04d}{month:02d}/"
                     f"e5.oper.an.sfc.{ERA5_CODES[var]}.ll025sc.{year:04d}{month:02d}0100_*.nc")
    return Path(hits[0]) if hits else None


def _nsidc_time_index(path: Path, year: int, mmdd: str) -> int | None:
    """Index of <year>-<mmdd> in a NSIDC daily file (time in days since 1970-01-01)."""
    day = (np.datetime64(f"{year:04d}-{mmdd}") - np.datetime64("1970-01-01")).astype(int)
    with nc.Dataset(path, "r") as ds:
        hits = np.flatnonzero(np.floor(ds.variables["time"][:]) == day)
    return int(hits[0]) if hits.size else None


def obs_file(var: str, year: int, mmdd: str, root: str | Path | None = None) -> Path | None:
    """File holding observed ``var`` on <year>-<mmdd> under ``root``, or None."""
    source = obs_source(root)
    if source == "era5":
        return era5_file(var, year, int(mmdd[:2]), era5_root=root)
    if source == "era5_land":
        path = Path(root) / ERA5_LAND_FILE.format(date=f"{year:04d}-{mmdd}")
        if not path.is_file():
            return None
        with nc.Dataset(path, "r") as ds:
            return path if var in ds.variables else None
    if var not in NSIDC_VARS:
        return None
    path = Path(root) / NSIDC_FILE.format(year=year)
    return path if path.is_file() and _nsidc_time_index(path, year, mmdd) is not None else None


def era5_available(
    target_date_or_year: int | str,
    mmdd_or_components: str | Sequence[str] | None = None,
    components: Sequence[str] | None = None,
    *,
    era5_root: str | Path | Mapping[str, str | Path | None] | None = None,
    obs_root: str | Path | Mapping[str, str | Path | None] | None = None,
    obs_fields: Any | None = None,
) -> bool:
    if isinstance(target_date_or_year, str):
        # target_date like "2015-02-01"
        year = int(target_date_or_year[:4])
        mmdd = target_date_or_year[5:]
        if isinstance(mmdd_or_components, (list, tuple)):
            actual_components = mmdd_or_components
        else:
            actual_components = components or ALL_COMPONENTS
    else:
        year = int(target_date_or_year)
        mmdd = str(mmdd_or_components)
        actual_components = components or ALL_COMPONENTS

    roots = _resolve_obs_roots(obs_root=obs_root, era5_root=era5_root)
    fmap = _resolve_obs_fields(obs_fields)
    for c in actual_components:
        if c not in fmap:
            continue
        c_root = roots.get(c, ERA5_ROOT)
        needed = {
            v
            for f in fmap[c].values()
            for v in (list(ERA5_LAYER_GROUPS[f]) if f in ERA5_LAYER_GROUPS else [f])
        }
        if not all(obs_file(v, year, mmdd, root=c_root) is not None for v in needed):
            return False
    return True


@dataclass
class SelectionPlan:
    target_date: str              # YYYY-MM-DD
    branch: str                   # "clim" or "obs"
    components: list[str]
    candidate_years: list[int]
    src_root: str
    src_prefix: str
    out_dir: str
    out_prefix: str
    excluded_years: list[int] = field(default_factory=list)
    era5_root: str = ERA5_ROOT
    obs_roots: dict[str, str] = field(default_factory=lambda: _resolve_obs_roots())
    obs_fields: dict[str, dict[str, str]] = field(default_factory=lambda: _resolve_obs_fields())
    label: str = ""               # output subdirectory; defaults to the branch
    clim_years: list[int] = field(default_factory=list)   # climatology (mean, spread); [] = candidate years

    def __post_init__(self):
        if not self.clim_years:
            self.clim_years = list(self.candidate_years)

    @property
    def years(self) -> list[int]:
        """Every archive year read: the candidate and climatology years."""
        return sorted(set(self.candidate_years) | set(self.clim_years))

    @property
    def clim_rows(self) -> list[int]:
        """Rows of the climatology years in ``years``."""
        return [self.years.index(y) for y in self.clim_years]

    @property
    def candidate_rows(self) -> list[int]:
        """Rows of the candidate years in ``years``."""
        return [self.years.index(y) for y in self.candidate_years]

    @property
    def mmdd(self) -> str:
        return self.target_date[5:]

    @property
    def target_year(self) -> int:
        return int(self.target_date[:4])

    def obs_root_for(self, component: str) -> str:
        return self.obs_roots.get(component, self.era5_root or ERA5_ROOT)

    def obs_fields_for(self, component: str) -> dict[str, str]:
        return self.obs_fields.get(component, OBS_FIELDS.get(component, {}))

    def source(self, component: str, year: int) -> Path:
        return source_file(self.src_root, self.src_prefix, component, year, self.mmdd)

    def output(self, component: str) -> Path:
        return Path(self.out_dir) / restart_name(self.out_prefix, component, self.target_year, self.mmdd)


def make_plan(
    target_date: str,
    *args,
    src_root: str | Path | None = None,
    out_root: str | Path | None = None,
    era5_root: str | Path | Mapping[str, str | Path | None] | None = None,
    obs_root: str | Path | Mapping[str, str | Path | None] | None = None,
    obs_fields: Any | None = None,
    branch: str = "clim",
    components: Sequence[str] = ALL_COMPONENTS,
    src_prefix: str = "20251024_s2d_spinup",
    out_prefix: str = "20251024_s2d_spinup",
    first_year: int | None = None,
    last_year: int | None = None,
    exclude_years: Iterable[int] = (),
    exclude_target_year: bool = True,
    label: str | None = None,
    clim_first_year: int | None = None,
    clim_last_year: int | None = None,
) -> SelectionPlan:
    """Candidate years (all components available) and output directory
    ``<out_root>/YYYYMMDD/<label>`` (``label`` defaults to the branch).

    The candidate years searched are ``first_year``-``last_year`` without
    ``exclude_years``. The climatology (the mean and interannual spread that
    standardize every field, model and observed) uses the available years
    ``clim_first_year``-``clim_last_year``; with both ``None`` it is the
    candidate years. The target year is left out of both when
    ``exclude_target_year``.

    If ``obs_root`` / ``era5_root`` or ``obs_fields`` are ``None`` (or omitted
    for a component), the default ERA5 setup is used automatically.
    """
    if len(args) == 1:
        if args[0] in BRANCHES:
            branch = args[0]
        else:
            src_root = args[0]
    elif len(args) == 2:
        if args[0] in BRANCHES:
            branch = args[0]
            src_root = args[1]
        else:
            src_root = args[0]
            out_root = args[1]
    elif len(args) >= 3:
        if args[0] in BRANCHES:
            branch = args[0]
            src_root = args[1]
            out_root = args[2]
        else:
            src_root = args[0]
            out_root = args[1]
            branch = args[2]

    if src_root is None or out_root is None:
        raise ValueError(
            f"make_plan requires both src_root and out_root, got src_root={src_root!r}, out_root={out_root!r}"
        )

    if branch not in BRANCHES:
        raise ValueError(f"branch must be one of {BRANCHES}, got {branch!r}")
    unknown = set(components) - set(COMPONENTS)
    if unknown:
        raise ValueError(f"Unknown components {sorted(unknown)}; use {list(COMPONENTS)}")
    src_root = Path(src_root).resolve()
    out_root = Path(out_root).resolve()
    resolved_roots = _resolve_obs_roots(obs_root=obs_root, era5_root=era5_root)
    resolved_fields = _resolve_obs_fields(obs_fields)
    if out_root == src_root or src_root in out_root.parents:
        raise ValueError(f"Output root {out_root} is inside the source archive {src_root}")

    target_year, mmdd = int(target_date[:4]), target_date[5:]
    if branch == "obs" and not era5_available(
        target_year, mmdd, components, obs_root=resolved_roots, obs_fields=resolved_fields
    ):
        raise FileNotFoundError(f"ERA5 for {target_date} is not available under {resolved_roots}")

    years = None
    for comp in components:
        found = set(available_years(src_root, src_prefix, comp, mmdd))
        years = found if years is None else years & found
    years = sorted(years)
    clim_years = []
    if clim_first_year is not None or clim_last_year is not None:
        clim_years = [y for y in years if (clim_first_year is None or y >= clim_first_year)
                      and (clim_last_year is None or y <= clim_last_year)
                      and not (exclude_target_year and y == target_year)]
        if len(clim_years) < 3:
            raise ValueError(f"Need at least three climatology years for {mmdd}; found {clim_years}")
    if first_year is not None:
        years = [y for y in years if y >= first_year]
    if last_year is not None:
        years = [y for y in years if y <= last_year]
    excluded = set(exclude_years)
    if exclude_target_year:
        excluded.add(target_year)
    years = [y for y in years if y not in excluded]
    if len(years) < 3:
        raise ValueError(f"Need at least three candidate years for {mmdd}; found {years}")

    return SelectionPlan(
        target_date=target_date,
        branch=branch,
        components=list(components),
        candidate_years=years,
        src_root=str(src_root),
        src_prefix=src_prefix,
        out_dir=str(out_root / target_date.replace("-", "") / (label or branch)),
        out_prefix=out_prefix,
        excluded_years=sorted(excluded),
        era5_root=resolved_roots.get("elm", ERA5_ROOT),
        obs_roots=resolved_roots,
        obs_fields=resolved_fields,
        label=label or branch,
        clim_years=clim_years,
    )


# ---------------------------------------------------------------------------
# Model key fields
# ---------------------------------------------------------------------------

def _read(ds: nc.Dataset, name: str) -> np.ndarray:
    var = ds.variables[name]
    var.set_auto_maskandscale(False)
    return np.asarray(var[...])


def _clean(data: np.ndarray) -> np.ndarray:
    data = np.asarray(data, dtype=np.float64)
    return np.where(np.abs(data) < SPVAL_THRESHOLD, data, np.nan)


def read_key_fields(component: str, path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Key fields of one restart as ``{name: (values, weights)}``, one value per point."""
    out = {}
    with nc.Dataset(path, "r") as ds:
        if component == "mpassi":
            lat = _read(ds, "latCell")
            cell = _read(ds, "areaCell")
            totals = {
                "ice concentration": _read(ds, "iceAreaCategory")[0, :, :, 0].sum(axis=1),
                "ice thickness": _read(ds, "iceVolumeCategory")[0, :, :, 0].sum(axis=1),
                "snow on ice": _read(ds, "snowVolumeCategory")[0, :, :, 0].sum(axis=1),
            }
            for hemi, mask in (("NH", lat > 0), ("SH", lat < 0)):
                for name, data in totals.items():
                    out[f"{hemi} {name}"] = (_clean(data), np.where(mask, cell, 0.0))
        elif component == "elm":
            k = len(ds.dimensions["levsno"])
            water = _clean(_read(ds, "H2OSOI_LIQ")[:, k:]) + _clean(_read(ds, "H2OSOI_ICE")[:, k:])
            weight = (_read(ds, "cols1d_wtxy") * np.cos(np.radians(_read(ds, "cols1d_lat")))
                      * (_read(ds, "cols1d_active") > 0))
            out["snow water equivalent"] = (_clean(_read(ds, "H2OSNO")), weight)
            out["surface soil water"] = (water[:, :3].sum(axis=1), weight)
            out["subsurface soil water"] = (water[:, 3:].sum(axis=1), weight)
            tsoi = _clean(_read(ds, "T_SOISNO")[:, k:])
            out["surface soil temperature"] = (tsoi[:, 0], weight)
            out["subsurface soil temperature"] = (
                tsoi[:, 3:10] @ ELM_TSUB_WEIGHTS / ELM_TSUB_WEIGHTS.sum(), weight)
        else:
            raise ValueError(f"No key fields for {component}")
    return out


def load_key_fields(
    plan: SelectionPlan,
    component: str,
    cache_dir: Path | None = None,
    years: Sequence[int] | None = None,
    log=print,
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Key fields as ``{name: (values[n_years, n_points], weights)}``.

    ``years`` defaults to the plan's ``years`` (candidate and climatology years). With ``cache_dir``,
    fields for every year read are kept (as float32) in
    ``<cache_dir>/<component>_<MM-DD>.npz`` and reused.
    """
    years = list(plan.years if years is None else years)
    names = list(KEY_FIELDS[component])
    cache = Path(cache_dir) / f"{component}_{plan.mmdd}.npz" if cache_dir else None
    stored: dict[int, list[np.ndarray]] = {}
    weights = None
    if cache is not None and cache.is_file():
        with np.load(cache) as z:
            if "names" in z.files and list(z["names"]) == names:   # else rebuild
                weights = [z[f"w{j}"] for j in range(len(names))]
                arrays = [z[f"f{j}"] for j in range(len(names))]   # each access loads the array
                for i, y in enumerate(z["years"]):
                    stored[int(y)] = [a[i] for a in arrays]

    missing = [y for y in years if y not in stored]
    t0 = time.time()
    for n, year in enumerate(missing):
        fields = read_key_fields(component, plan.source(component, year))
        stored[year] = [fields[name][0] for name in names]
        weights = [fields[name][1] for name in names]
        log(f"  {component} {plan.mmdd} {year}: key fields read "
            f"({n + 1}/{len(missing)}, {time.time() - t0:.0f} s)")

    if cache is not None and missing:
        cache.parent.mkdir(parents=True, exist_ok=True)
        all_years = sorted(stored)
        arrays = {f"f{j}": np.stack([stored[y][j] for y in all_years]).astype(np.float32)
                  for j in range(len(names))}
        arrays.update({f"w{j}": weights[j] for j in range(len(names))})
        tmp = cache.with_name(f".{cache.stem}.{os.getpid()}.npz")
        np.savez(tmp, years=np.array(all_years), names=np.array(names), **arrays)
        os.replace(tmp, cache)

    return {name: (np.stack([stored[y][j] for y in years]).astype(np.float64), weights[j])
            for j, name in enumerate(names)}


# ---------------------------------------------------------------------------
# Observation (ERA5) fields at model comparison points
# ---------------------------------------------------------------------------

def comparison_points(plan: SelectionPlan, component: str, cache_dir: Path | None = None) -> dict:
    """Coordinates of the model points.

    Sea ice: MPAS cells. Land: ELM grid cells, with the column-to-grid-cell map
    (``col_grid``).
    """
    cache = Path(cache_dir) / f"points_{component}.npz" if cache_dir else None
    if cache is not None and cache.is_file():
        with np.load(cache) as z:
            return {k: z[k] for k in z.files}
    with nc.Dataset(plan.source(component, plan.candidate_years[0]), "r") as ds:
        if component == "mpassi":
            pts = {"lat": np.degrees(_read(ds, "latCell")), "lon": np.degrees(_read(ds, "lonCell"))}
        elif component == "elm":
            pts = {"lat": _read(ds, "grid1d_lat"), "lon": _read(ds, "grid1d_lon"),
                   "col_grid": _read(ds, "cols1d_gridcell_index").astype(np.int64) - 1}
        else:
            raise ValueError(component)
    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        tmp = cache.with_name(f".{cache.stem}.{os.getpid()}.npz")
        np.savez(tmp, **pts)
        os.replace(tmp, cache)
    return pts


def _era5_slice(var: str, year: int, mmdd: str, era5_root: str | Path = ERA5_ROOT) -> np.ndarray:
    """ERA5 field at 00Z on <year>-<mmdd>, shape (721, 1440), NaN where missing."""
    path = era5_file(var, year, int(mmdd[:2]), era5_root=era5_root)
    if path is None:
        raise FileNotFoundError(f"ERA5 {var} for {year}-{mmdd} not found under {era5_root}")
    with nc.Dataset(path, "r") as ds:
        name = next(n for n, v in ds.variables.items() if v.ndim == 3)
        data = ds.variables[name][(int(mmdd[3:]) - 1) * 24]
    return np.ma.filled(data.astype(np.float64), np.nan)


def _era5_index(lat: np.ndarray, lon: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    i = np.clip(np.rint((90.0 - lat) / 0.25).astype(np.int64), 0, 720)
    j = np.rint(np.mod(lon, 360.0) / 0.25).astype(np.int64) % 1440
    return i, j


def _obs_slice(var: str, year: int, mmdd: str, root: str | Path) -> np.ndarray:
    """Observed ``var`` on <year>-<mmdd> on the native grid of ``root``, NaN where missing."""
    source = obs_source(root)
    if source == "era5":
        return _era5_slice(var, year, mmdd, era5_root=root)
    path = obs_file(var, year, mmdd, root=root)
    if path is None:
        raise FileNotFoundError(f"{source} {var} for {year}-{mmdd} not found under {root}")
    with nc.Dataset(path, "r") as ds:
        if source == "era5_land":
            data = ds.variables[var][0]
        else:
            data = ds.variables[NSIDC_VARS[var]][_nsidc_time_index(path, year, mmdd)]
    data = np.ma.filled(np.ma.asarray(data).astype(np.float64), np.nan)
    data[np.abs(data) > SPVAL_THRESHOLD] = np.nan
    return data


def _obs_index(lat: np.ndarray, lon: np.ndarray, root: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """Nearest-cell indices of (lat, lon) on the grid of ``root``."""
    if obs_source(root) == "era5":
        return _era5_index(lat, lon)
    i = np.clip(np.floor(lat + 90.0).astype(np.int64), 0, 179)     # lat -89.5 .. 89.5
    j = np.floor(np.mod(lon, 360.0)).astype(np.int64) % 360        # lon 0.5 .. 359.5
    return i, j


def read_obs_fields(plan: SelectionPlan, component: str, year: int, points: dict) -> list[np.ndarray]:
    """Observed values of the component's observable fields at the comparison points."""
    c_root = plan.obs_root_for(component)
    i, j = _obs_index(points["lat"], points["lon"], c_root)
    out = []
    cache: dict[str, np.ndarray] = {}

    def era5(var):
        if var not in cache:
            cache[var] = _obs_slice(var, year, plan.mmdd, c_root)[i, j]
        return cache[var]

    for var in plan.obs_fields_for(component).values():
        if var == "swvl2-4":      # total water over the layers (m water)
            out.append(sum(era5(v) * dz for v, dz in ERA5_LAYER_GROUPS[var].items()))
        elif var in ERA5_LAYER_GROUPS:     # depth-weighted mean (temperature)
            group = ERA5_LAYER_GROUPS[var]
            out.append(sum(era5(v) * dz for v, dz in group.items()) / sum(group.values()))
        else:
            out.append(era5(var))
    return out


def load_obs_fields(
    plan: SelectionPlan,
    component: str,
    years: Sequence[int],
    cache_dir: Path | None = None,
    log=print,
) -> dict[str, np.ndarray]:
    """ERA5 fields ``{name: values[n_years, n_points]}`` at the comparison points.

    With ``cache_dir``, every year read is kept in
    ``<cache_dir>/obs_<component>_<MM-DD>.npz`` (ERA5) or
    ``obs_<component>_<source>_<MM-DD>.npz`` (other observation sets) and reused.
    """
    names = list(plan.obs_fields_for(component))
    points = comparison_points(plan, component, cache_dir)
    source = obs_source(plan.obs_root_for(component))
    tag = "" if source == "era5" else f"{source}_"
    cache = Path(cache_dir) / f"obs_{component}_{tag}{plan.mmdd}.npz" if cache_dir else None
    stored: dict[int, list[np.ndarray]] = {}
    if cache is not None and cache.is_file():
        with np.load(cache) as z:
            if "names" in z.files and list(z["names"]) == names:   # else rebuild
                arrays = [z[f"f{j}"] for j in range(len(names))]
                for i, y in enumerate(z["years"]):
                    stored[int(y)] = [a[i] for a in arrays]

    missing = [y for y in years if y not in stored]
    t0 = time.time()
    for n, year in enumerate(missing):
        stored[year] = read_obs_fields(plan, component, year, points)
        log(f"  {source} {component} {year}-{plan.mmdd} read ({n + 1}/{len(missing)}, {time.time() - t0:.0f} s)")

    if cache is not None and missing:
        cache.parent.mkdir(parents=True, exist_ok=True)
        all_years = sorted(stored)
        arrays = {f"f{j}": np.stack([stored[y][j] for y in all_years]).astype(np.float32)
                  for j in range(len(names))}
        tmp = cache.with_name(f".{cache.stem}.{os.getpid()}.npz")
        np.savez(tmp, years=np.array(all_years), names=np.array(names), **arrays)
        os.replace(tmp, cache)

    return {name: np.stack([stored[y][j] for y in years]).astype(np.float64)
            for j, name in enumerate(names)}


def _columns_to_gridcells(values: np.ndarray, weight: np.ndarray, col_grid: np.ndarray,
                          ngrid: int) -> tuple[np.ndarray, np.ndarray]:
    """Weighted mean of column values per grid cell; returns (values, grid-cell weights)."""
    valid = np.all(np.isfinite(values), axis=0)
    w = np.where(valid, weight, 0.0)
    wg = np.bincount(col_grid, weights=w, minlength=ngrid)
    out = np.stack([np.bincount(col_grid, weights=w * np.where(valid, v, 0.0), minlength=ngrid)
                    for v in values])
    out = out / np.where(wg > 0, wg, 1.0)
    out[:, wg <= 0] = np.nan
    return out, wg


def box_offsets(box: tuple[float, float] | None, shifts: int = 1) -> list[tuple[float, float]]:
    """Grid offsets (dlat, dlon degrees) for ``shifts`` placements of the box edges.

    1: the standard grid; 2: also shifted by half a box in both directions;
    4: also shifted by half a box in latitude only and in longitude only.
    """
    if box is None or shifts == 1:
        return [(0.0, 0.0)]
    dlat, dlon = box
    if shifts == 2:
        return [(0.0, 0.0), (dlat / 2, dlon / 2)]
    if shifts == 4:
        return [(0.0, 0.0), (dlat / 2, 0.0), (0.0, dlon / 2), (dlat / 2, dlon / 2)]
    raise ValueError(f"box shifts must be 1, 2, or 4, got {shifts}")


def _to_boxes(values: np.ndarray, weight: np.ndarray, lat: np.ndarray, lon: np.ndarray,
              box: tuple[float, float], offset: tuple[float, float] = (0.0, 0.0),
              ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Weighted means over lat-lon boxes of ``box = (dlat, dlon)`` degrees.

    ``offset = (olat, olon)`` shifts the box edges by that many degrees (the
    polar bands are then partial boxes). Returns (box values, box weights, box
    index of each point). Points missing in any row are left out.
    """
    dlat, dlon = box
    olat, olon = offset
    nlat = int(np.ceil((180 + olat) / dlat - 1e-9))
    nlon = int(round(360 / dlon))
    index = (np.minimum(np.floor((lat + 90 + olat) / dlat), nlat - 1).astype(np.int64) * nlon
             + np.floor(np.mod(lon + olon, 360) / dlon).astype(np.int64) % nlon)
    nbox = nlat * nlon
    valid = np.all(np.isfinite(values), axis=0)
    w = np.where(valid, weight, 0.0)
    wb = np.bincount(index, weights=w, minlength=nbox)
    out = np.stack([np.bincount(index, weights=w * np.where(valid, v, 0.0), minlength=nbox)
                    for v in values])
    out = out / np.where(wb > 0, wb, 1.0)
    out[:, wb <= 0] = np.nan
    return out, wb, index


def _apply_box(values, weight, lat, lon, box, obs=None, offset=(0.0, 0.0)):
    """Aggregate model (and observed) values to boxes with the model weights.

    With ``obs``, both are averaged over the same points: those valid in every
    model and observed row. ``offset`` shifts the box edges (see ``_to_boxes``).
    """
    if box is None:
        return values, weight, obs
    if obs is not None:
        valid = np.all(np.isfinite(values), axis=0) & np.all(np.isfinite(obs), axis=0)
        values = np.where(valid, values, np.nan)
        obs = np.where(valid, obs, np.nan)
    agg, wb, _ = _to_boxes(values, weight, lat, lon, box, offset)
    if obs is not None:
        obs, _, _ = _to_boxes(obs, weight, lat, lon, box, offset)
    return agg, wb, obs


# ---------------------------------------------------------------------------
# Distances and selection
# ---------------------------------------------------------------------------

def field_distance(
    values: np.ndarray,
    weight: np.ndarray,
    *,
    target: np.ndarray | None = None,
    obs_years: np.ndarray | None = None,
    obs_target: np.ndarray | None = None,
    ref: Sequence[int] | None = None,
) -> np.ndarray:
    """Distance d_f of each year (row of ``values``) from a target.

    No target: the climatology (``t = 0``). ``target``: a model state at the
    same points (standardized with the climatological mean and spread).
    ``obs_years``/``obs_target``: observed values for the same years and
    the target date (standardized with their own climatological mean and spread).
    ``ref`` are the rows that form the climatology (default: all rows).
    """
    valid = np.all(np.isfinite(values), axis=0)
    if target is not None:
        valid &= np.isfinite(target)
    if obs_years is not None:
        valid &= np.all(np.isfinite(obs_years), axis=0) & np.isfinite(obs_target)
    w = np.where(valid, weight, 0.0)
    wsum = w.sum()
    if wsum <= 0:
        return np.zeros(values.shape[0])

    rows = slice(None) if ref is None else list(ref)

    def standardize(v, t=None):
        v = np.where(valid, v, 0.0)
        mean = v[rows].mean(axis=0)
        anom = v - mean
        spread = np.sqrt((w * (anom[rows] ** 2).mean(axis=0)).sum() / wsum)
        if spread <= 0:
            return np.zeros_like(anom), None if t is None else np.zeros(v.shape[1])
        zt = None if t is None else (np.where(valid, t, 0.0) - mean) / spread
        return anom / spread, zt

    z, zt = standardize(values, target)
    if obs_years is not None:
        _, zt = standardize(obs_years, obs_target)
    if zt is None:
        zt = 0.0
    return np.sqrt((w * (z - zt) ** 2).sum(axis=1) / wsum)


def _combine(table: pd.DataFrame, components: Sequence[str],
             component_weights: Mapping[str, float] | None) -> pd.DataFrame:
    a = {c: 1.0 for c in components}
    a.update(component_weights or {})
    total = sum(a[c] for c in components)
    table["joint"] = np.sqrt(sum(a[c] * table[c] ** 2 for c in components) / total)
    return table


def _model_points(plan: SelectionPlan, component: str, cache_dir) -> tuple[np.ndarray, np.ndarray]:
    """Latitude and longitude of the model key-field points (ELM: columns)."""
    pts = comparison_points(plan, component, cache_dir)
    if component == "elm":
        return pts["lat"][pts["col_grid"]], pts["lon"][pts["col_grid"]]
    return pts["lat"], pts["lon"]


def compute_distances(
    plan: SelectionPlan,
    *,
    box: tuple[float, float] | None = None,
    box_shifts: int = 1,
    component_weights: Mapping[str, float] | None = None,
    cache_dir: Path | None = None,
    log=print,
) -> pd.DataFrame:
    """Distance table for the plan's branch, indexed by candidate year.

    ``box = (dlat, dlon)`` compares area-weighted means over lat-lon boxes of
    that size instead of individual model points (ERA5 is averaged over the
    same points and boxes). With ``box_shifts`` 2 or 4, each field distance is
    the RMS over that many placements of the box edges (``box_offsets``), so
    the result does not hinge on where the edges fall.

    Columns: ``<component>: <field>`` (d_f), ``<component>`` (D_c), ``joint`` (D).
    """
    table = pd.DataFrame(index=pd.Index(plan.candidate_years, name="year"))
    used = []
    ref, cand = plan.clim_rows, plan.candidate_rows
    for comp in plan.components:
        obs_comp_fields = plan.obs_fields_for(comp)
        if comp not in KEY_FIELDS or (plan.branch == "obs" and not obs_comp_fields):
            continue
        model = load_key_fields(plan, comp, cache_dir, years=plan.years, log=log)
        cols = []
        if plan.branch == "clim":
            lat, lon = _model_points(plan, comp, cache_dir)
            for name, (values, weight) in model.items():
                d2 = [field_distance(*_apply_box(values, weight, lat, lon, box, offset=off)[:2], ref=ref)[cand] ** 2
                      for off in box_offsets(box, box_shifts)]
                table[f"{comp}: {name}"] = np.sqrt(np.mean(d2, axis=0))
                cols.append(f"{comp}: {name}")
        else:
            obs = load_obs_fields(plan, comp, plan.years + [plan.target_year], cache_dir, log)
            points = comparison_points(plan, comp, cache_dir)
            for name in obs_comp_fields:
                values, weight = model[name]
                if comp == "elm":
                    values, weight = _columns_to_gridcells(values, weight, points["col_grid"],
                                                           points["lat"].size)
                d2 = []
                for off in box_offsets(box, box_shifts):
                    vb, wb, o = _apply_box(values, weight, points["lat"], points["lon"], box,
                                           obs=obs[name], offset=off)
                    d2.append(field_distance(vb, wb, obs_years=o[:-1], obs_target=o[-1], ref=ref)[cand] ** 2)
                table[f"{comp}: {name}"] = np.sqrt(np.mean(d2, axis=0))
                cols.append(f"{comp}: {name}")
        table[comp] = np.sqrt((table[cols] ** 2).mean(axis=1))
        used.append(comp)
    return _combine(table, used, component_weights)


def distance_to_model_year(
    plan: SelectionPlan,
    target_year: int,
    *,
    box: tuple[float, float] | None = None,
    observable_only: bool = False,
    component_weights: Mapping[str, float] | None = None,
    cache_dir: Path | None = None,
    log=print,
) -> pd.DataFrame:
    """Distance of each candidate year from the model state of ``target_year``.

    Used to check a selection when the actual restart of the target year
    exists (the 2015 dates): same standardization as ``clim``, with the
    target year's anomaly as the target. ``observable_only`` restricts it to
    the fields used by the ``obs`` branch.
    """
    table = pd.DataFrame(index=pd.Index(plan.candidate_years, name="year"))
    used = []
    for comp in plan.components:
        obs_comp_fields = plan.obs_fields_for(comp)
        if comp not in KEY_FIELDS or (observable_only and not obs_comp_fields):
            continue
        fields = load_key_fields(plan, comp, cache_dir, years=plan.years + [target_year], log=log)
        lat, lon = _model_points(plan, comp, cache_dir)
        cols = []
        for name, (values, weight) in fields.items():
            if observable_only and name not in obs_comp_fields:
                continue
            values, weight, _ = _apply_box(values, weight, lat, lon, box)
            table[f"{comp}: {name}"] = field_distance(values[:-1], weight, target=values[-1],
                                                      ref=plan.clim_rows)[plan.candidate_rows]
            cols.append(f"{comp}: {name}")
        table[comp] = np.sqrt((table[cols] ** 2).mean(axis=1))
        used.append(comp)
    return _combine(table, used, component_weights)


def input_summary(plan: SelectionPlan, cache_dir: Path | None = None, log=print) -> pd.DataFrame:
    """Read (and cache) the key fields of ``plan`` and summarize them, one row per field and source.

    Model: every archive year read (``plan.years``). Observations (fields with an
    observed counterpart in the plan's observation set, if available for the
    target date): the same years and the target date. Columns give the years
    and points read, the area-weighted climatological mean, and the
    interannual spread ``S`` over the climatology years.
    """
    ref = plan.clim_rows
    obs_ok = era5_available(plan.target_year, plan.mmdd,
                            [c for c in plan.components if c in KEY_FIELDS and plan.obs_fields_for(c)],
                            obs_root=plan.obs_roots, obs_fields=plan.obs_fields)
    rows = []
    for comp in plan.components:
        if comp not in KEY_FIELDS:
            continue
        model = load_key_fields(plan, comp, cache_dir, years=plan.years, log=log)
        points = comparison_points(plan, comp, cache_dir)
        obs_names = plan.obs_fields_for(comp) if obs_ok else {}
        obs = load_obs_fields(plan, comp, plan.years + [plan.target_year], cache_dir, log) if obs_names else {}
        for name, (values, weight) in model.items():
            sources = [("model", "", values, weight)]
            if name in obs:
                wg = weight
                if comp == "elm":
                    _, wg = _columns_to_gridcells(values, weight, points["col_grid"], points["lat"].size)
                sources.append(("observed", obs_names[name], obs[name][:-1], wg))
            for source, var, v, w in sources:
                valid = np.all(np.isfinite(v), axis=0) & (w > 0)
                wv = np.where(valid, w, 0.0)
                v = np.where(valid, v, 0.0)[ref]
                mean = v.mean(axis=0)
                rows.append({
                    "target_date": plan.target_date, "component": comp, "field": name, "source": source,
                    "variable": var or "restart", "years read": f"{plan.years[0]}-{plan.years[-1]} ({len(plan.years)})",
                    "points": int(valid.sum()),
                    "clim mean": float((wv * mean).sum() / wv.sum()) if wv.sum() > 0 else np.nan,
                    "spread S": float(np.sqrt((wv * ((v - mean) ** 2).mean(axis=0)).sum() / wv.sum()))
                    if wv.sum() > 0 else np.nan,
                })
    return pd.DataFrame(rows)


def select_years(distances: pd.DataFrame, components: Sequence[str], mode: str = "joint") -> dict[str, int]:
    """Selected year per component.

    ``joint``: the year with the smallest joint distance for every component.
    ``per_component``: each component's smallest D_c. A component in
    ``FOLLOWS`` (the river) takes the year of the component it follows (the
    land); any other component without a distance takes the joint year.
    """
    if mode not in ("joint", "per_component"):
        raise ValueError(f"mode must be 'joint' or 'per_component', got {mode!r}")
    joint = int(distances["joint"].idxmin())
    if mode == "joint":
        return {comp: joint for comp in components}
    own = {c: int(distances[c].idxmin()) for c in components if c in distances}
    return {c: own.get(c, own.get(FOLLOWS.get(c), joint)) for c in components}


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------

def write_selected_restarts(
    plan: SelectionPlan,
    selection: Mapping[str, int],
    distances: pd.DataFrame,
    *,
    mode: str = "joint",
    box: tuple[float, float] | None = None,
    box_shifts: int = 1,
    component_weights: Mapping[str, float] | None = None,
    overwrite: bool = False,
    log=print,
) -> dict:
    """Copy each selected restart to the target-date file name and write rpointers.

    ``selection.json`` next to the outputs records the branch, candidate
    years, selected year per component, and the full distance table.
    """
    out_dir = Path(plan.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    copied = {}
    previous = {}
    if (out_dir / "selection.json").is_file():
        previous = json.loads((out_dir / "selection.json").read_text()).get("selection", {})
    for comp, year in selection.items():
        src, dst = plan.source(comp, year), plan.output(comp)
        info = COMPONENTS[comp]
        unchanged = (dst.is_file() and previous.get(comp, {}).get("source") == str(src)
                     and dst.stat().st_size == src.stat().st_size)
        if unchanged:
            log(f"  {comp}: {dst} already holds {year}; not copied")
        else:
            if dst.exists() and not overwrite:
                raise FileExistsError(f"{dst} exists; pass overwrite=True to replace it")
            t0 = time.time()
            tmp = dst.with_name(f".{dst.name}.tmp")
            shutil.copyfile(src, tmp)
            os.replace(tmp, dst)
            log(f"  {comp}: {src.name} -> {dst} ({time.time() - t0:.0f} s)")
        pointer = f"./{dst.name}" if info["pointer"] == "file" else f"{plan.target_date}_00:00:00"
        (out_dir / info["rpointer"]).write_text(pointer + "\n")
        copied[comp] = {"year": int(year), "source": str(src), "output": str(dst)}

    summary = {
        "plan": asdict(plan),
        "mode": mode,
        "box_degrees": list(box) if box else None,
        "box_shifts": box_shifts,
        "component_weights": dict(component_weights or {}),
        "selection": copied,
        "fields": [c for c in distances.columns if ": " in c],
        "distances": json.loads(distances.round(4).to_json(orient="index")),
    }
    (out_dir / "selection.json").write_text(json.dumps(summary, indent=2))
    return summary


def run_selection(
    plan: SelectionPlan,
    *,
    box: tuple[float, float] | None = None,
    box_shifts: int = 1,
    mode: str = "joint",
    component_weights: Mapping[str, float] | None = None,
    cache_dir: Path | None = None,
    write: bool = True,
    overwrite: bool = False,
    log=print,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Compute distances, select the year(s), and optionally copy the files."""
    log(f"{plan.target_date} [{plan.branch}]: {len(plan.candidate_years)} candidate years "
        f"({plan.candidate_years[0]}-{plan.candidate_years[-1]}, excluded {plan.excluded_years}), "
        f"climatology {plan.clim_years[0]}-{plan.clim_years[-1]} ({len(plan.clim_years)} years)")
    distances = compute_distances(plan, box=box, box_shifts=box_shifts,
                                  component_weights=component_weights, cache_dir=cache_dir, log=log)
    selection = select_years(distances, plan.components, mode)
    log(f"  selected ({mode}): {selection}")
    if write:
        write_selected_restarts(plan, selection, distances, mode=mode, box=box, box_shifts=box_shifts,
                                component_weights=component_weights, overwrite=overwrite, log=log)
    return distances, selection


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

def summarize_restart(path: Path, component: str) -> dict[str, float]:
    """Global summary numbers for comparing restarts."""
    out = {}
    with nc.Dataset(path, "r") as ds:
        if component == "mpassi":
            cell = _read(ds, "areaCell")
            lat = _read(ds, "latCell")
            area = _read(ds, "iceAreaCategory")[0, :, :, 0].sum(axis=1)
            vol = _read(ds, "iceVolumeCategory")[0, :, :, 0].sum(axis=1)
            snow = _read(ds, "snowVolumeCategory")[0, :, :, 0].sum(axis=1)
            for hemi, mask in (("NH", lat > 0), ("SH", lat < 0)):
                out[f"{hemi} ice area (1e6 km2)"] = float((area * cell)[mask].sum() / 1e12)
                out[f"{hemi} ice extent >15% (1e6 km2)"] = float(cell[mask & (area > 0.15)].sum() / 1e12)
                out[f"{hemi} ice volume (1e3 km3)"] = float((vol * cell)[mask].sum() / 1e12)
                out[f"{hemi} snow volume (1e3 km3)"] = float((snow * cell)[mask].sum() / 1e12)
        elif component == "elm":
            w = np.where(_read(ds, "cols1d_active") > 0, _read(ds, "cols1d_wtxy"), 0.0)

            def wmean(x):
                good = np.abs(x) < SPVAL_THRESHOLD
                ww = np.where(good, w, 0.0)
                return float((np.where(good, x, 0.0) * ww).sum() / ww.sum())

            k = len(ds.dimensions["levsno"])
            out["mean H2OSNO (mm)"] = wmean(_read(ds, "H2OSNO"))
            out["snow-covered column fraction"] = float(w[_read(ds, "SNLSNO") < 0].sum() / w.sum())
            out["top-3-layer soil liquid (kg m-2)"] = wmean(_read(ds, "H2OSOI_LIQ")[:, k:k + 3].sum(axis=1))
            out["top-3-layer soil ice (kg m-2)"] = wmean(_read(ds, "H2OSOI_ICE")[:, k:k + 3].sum(axis=1))
            out["top soil layer T (K)"] = wmean(_read(ds, "T_SOISNO")[:, k])
            out["mean water table depth ZWT (m)"] = wmean(_read(ds, "ZWT"))
        elif component == "mosart":
            for name in ("RTM_VOLR_LIQ", "RTM_WR_LIQ", "RTM_WT_LIQ"):
                out[f"total {name} (km3)"] = float(_read(ds, name).sum() / 1e9)
    return out
