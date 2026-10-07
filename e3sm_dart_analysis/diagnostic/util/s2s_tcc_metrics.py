from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Sequence, Mapping, Tuple, Dict, Optional, List

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from matplotlib.colors import Normalize
from matplotlib.colors import BoundaryNorm, ListedColormap, LinearSegmentedColormap
from matplotlib.colors import TwoSlopeNorm
from matplotlib.ticker import FixedLocator

from matplotlib.ticker import FixedLocator

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
    HAS_CARTOPY = True
except ImportError:
    HAS_CARTOPY = False

from util.s2s_experiments import (
    build_experiments,
    DEFAULT_EXPERIMENTS
)

from util.s2s_observations import (
    OBS_REGISTRY,
    get_obs_file,
    list_obs_sources,
    obs_coverage,
)

from util.s2s_time_windows import (
    build_windows_from_start
)
from util.array_stats import nanquantile

@dataclass(frozen=True)
class FileCollectionConfig:
    group: str
    freq: str
    run: str
    obs: str
    period: str
    model_template: str            # <-- REQUIRED, explicit
    ens_prefix: str = "EN"
    ens_width: int = 2
    ens_start: int = 1

    def parse_period(self) -> Tuple[int, int, int, int]:
        m = re.match(r"^(\d{4})(\d{2})-(\d{4})(\d{2})$", self.period)
        if not m:
            raise ValueError(f"Bad period '{self.period}', expected 'YYYYMM-YYYYMM'")
        y0, m0, y1, m1 = map(int, m.groups())
        if (y1, m1) < (y0, m0):
            raise ValueError(f"Bad period '{self.period}' (end < start)")
        return y0, m0, y1, m1

    def years(self) -> List[int]:
        y0, _, y1, _ = self.parse_period()
        return list(range(y0, y1 + 1))

    def ens_labels(self, nens: int) -> List[str]:
        # FIX: use self.ens_start
        return [
            f"{self.ens_prefix}{i:0{self.ens_width}d}"
            for i in range(self.ens_start, nens + self.ens_start)
        ]


class S2SFileCollector:
    """
    Collect obs + model files for S2S / ACC diagnostics.
    """

    def __init__(
        self,
        *,
        exp_list: dict,
        exp_dict: dict,
        obs_registry,
        s2s_var_dict: Dict[str, str],
        get_obs_file_func,
    ):
        self.exp_list = exp_list
        self.exp_dict = exp_dict
        self.obs_registry = obs_registry
        self.s2s_var_dict = s2s_var_dict
        self.get_obs_file = get_obs_file_func

    # ---------- path resolvers ----------

    def resolve_obs_file(
        self,
        obs: str,
        freq: str,
        year: int,
        var: Optional[str] = None,
    ) -> str:
        return self.get_obs_file(
            self.obs_registry,
            obs,
            freq=freq,
            year=year,
            var=var,
        )

    def model_ts_dir(self, run_meta, freq: str, ens: Optional[str] = None) -> str:
        atm_path = run_meta.atm_path
        if ens:
            marker = os.sep + "archive" + os.sep
            if marker not in atm_path:
                raise ValueError(f"Cannot insert member directory into {atm_path}")
            atm_path = atm_path.replace(marker, os.sep + ens + marker, 1)
        return os.path.join(atm_path, "ts", freq)

    def resolve_model_file(
        self,
        *,
        run_meta,
        freq: str,
        year: int,
        var: str,
        ens: str,
        template: str,
    ) -> str:
        ts_dir = self.model_ts_dir(run_meta, freq, ens)
        fname = template % {
            "var": var,
            "ens": ens,
            "year": year,
            "period": run_meta.period,
        }
        return os.path.join(ts_dir, fname)

    # ---------- model file discovery ----------

    def candidates_model_files(
        self,
        *,
        ts_dir: str,
        var: str,
        ens: str,
        years: List[int],
        period: str,
        templates: List[str],
    ) -> List[str]:
        out: List[str] = []
        for tpl in templates:
            for y in years:
                name = tpl % {
                    "var": var,
                    "ens": ens,
                    "year": y,
                    "period": period,
                }
                path = os.path.join(ts_dir, name)
                if os.path.exists(path):
                    out.append(path)

        # de-duplicate while preserving order
        seen = set()
        uniq: List[str] = []
        for p in out:
            if p not in seen:
                uniq.append(p)
                seen.add(p)
        return uniq

    def model_files_for_ens(
        self,
        cfg: FileCollectionConfig,
        *,
        run_meta,
        var: str,
        ens: str,
    ) -> List[str]:
        ts_dir = self.model_ts_dir(run_meta, cfg.freq, ens)
        years = cfg.years()
        return [
            os.path.join(
                ts_dir,
                cfg.model_template % {"var": var, "ens": ens, "year": y, "period": cfg.period},
            )
            for y in years
        ]

    # ---------- main collection ----------

    def collect_one_var(
        self,
        cfg: FileCollectionConfig,
        *,
        var: str,
        obs_var: str,
        verbose: bool = True,
    ):
        years = cfg.years()
        models = self.exp_list[cfg.group]["models"]

        # --- obs files ---
        obs_files = [self.resolve_obs_file(cfg.obs, cfg.freq, y, obs_var) for y in years]
        missing_obs = [f for f in obs_files if not os.path.exists(f)]
        if missing_obs:
            if verbose:
                print(f"[MISSING OBS] {var} (obs_var={obs_var})")
                for f in missing_obs:
                    print("  ", f)
            return None

        out_var = {}

        for exp in models:
            run_meta = self.exp_dict[exp]["runs"][cfg.run]
            if run_meta is None:
                if verbose:
                    print(f"SKIP {exp}: no run='{cfg.run}'")
                continue

            nens = int(self.exp_dict[exp]["nens"])
            ts_dir = self.model_ts_dir(run_meta, cfg.freq)

            model_by_ens = {}
            missing_any = False

            for ens in cfg.ens_labels(nens):
                files = self.model_files_for_ens(cfg, run_meta=run_meta, var=var, ens=ens)

                missing = [f for f in files if not os.path.exists(f)]
                if missing:
                    missing_any = True
                    if verbose:
                        print(f"[MISSING MOD] {exp} {var} {ens} (ts_dir={ts_dir})")
                        for f in missing[:5]:
                            print("   ", f)
                        if len(missing) > 5:
                            print(f"   ... {len(missing)-5} more")

                model_by_ens[ens] = files  # keep full list for later loading

            out_var[exp] = {
                "obs": obs_files,
                "model": model_by_ens,
                "ts_dir": ts_dir,
                "nens": nens,
            }

            if verbose and missing_any:
                print(f"[WARN] {exp} {var}: some ensemble members missing files")

        return out_var

    def collect(
        self,
        cfg: FileCollectionConfig,
        *,
        vars_to_process: Optional[List[str]] = None,
        verbose: bool = True,
    ) -> Dict[str, dict]:
        out: Dict[str, dict] = {}
        for var, obs_var in self.s2s_var_dict.items():
            if vars_to_process is not None and var not in vars_to_process:
                continue

            res = self.collect_one_var(cfg, var=var, obs_var=obs_var, verbose=verbose)
            if res is None:
                continue
            out[var] = res

        return out


# ============================================================
# Physics normalization registries (unchanged)
# ============================================================

FLUX_SIGN_REGISTRY = {
    ("ERA5",  "FLNS"):  +1,
    ("ERA5",  "FLNSC"): +1,
    ("ERA5",  "FSNS"):  +1,
    ("ERA5",  "FSNSC"): +1,
    ("ERA5",  "FLDS"):  +1,
    ("ERA5",  "FSDS"):  +1,
    ("ERA5",  "FLUT"):  +1,

    ("MODEL", "FLNS"):  +1,
    ("MODEL", "FLNSC"): +1,
    ("MODEL", "FSNS"):  +1,
    ("MODEL", "FSNSC"): +1,
    ("MODEL", "FLDS"):  +1,
    ("MODEL", "FSDS"):  +1,
    ("MODEL", "FLUT"):  +1,
}


# ============================================================
# Window configs
# ============================================================

@dataclass(frozen=True)
class WindowedTCCBaseConfig:
    windows: Dict[str, Tuple[str, str]]  # label -> (start_date, end_date), ISO yyyy-mm-dd
    regions: Optional[List[dict]] = None  # [{"name":..., "lat":(lo,hi), "lon":(lo,hi)}]
    lat_range: Optional[Tuple[float, float]] = None
    lon_range: Optional[Tuple[float, float]] = None

    # grid/time checks
    require_same_grid: bool = True
    grid_tol: float = 0.0
    require_nonempty_time: bool = True

    # output
    out_dir: str = "./base_tcc_rmse"
    overwrite: bool = False
    verbose: bool = True


@dataclass(frozen=True)
class WindowedTCCAggregateConfig:
    member_quantiles: Tuple[float, ...] = (0.1, 0.5, 0.9)
    overwrite: bool = True
    verbose: bool = True


# ============================================================
# Stage 1: base builder (UPDATED to compute/save BIAS)
# ============================================================

class S2SWindowedTCCBaseBuilder:
    """
    Stage 1: build base files per (exp, var, ens).
      output file contains:
        - TCC(window, lat, lon)
        - RMSE(window, lat, lon)
        - BIAS(window, lat, lon)   <-- NEW
        - coords: window_start/window_end
        - attrs: exp/ens/var/region/period
    """

    def __init__(self, cfg: WindowedTCCBaseConfig):
        self.cfg = cfg
        os.makedirs(self.cfg.out_dir, exist_ok=True)

    # ---------- helpers ----------

    @staticmethod
    def _parse_yyyymm_period(period: str) -> Tuple[int, int, int, int]:
        m = re.match(r"^(\d{4})(\d{2})-(\d{4})(\d{2})$", period)
        if not m:
            raise ValueError(f"Bad period '{period}', expected 'YYYYMM-YYYYMM'")
        y0, m0, y1, m1 = map(int, m.groups())
        return y0, m0, y1, m1

    @classmethod
    def _select_period_time(cls, da: xr.DataArray, period: str) -> xr.DataArray:
        y0, m0, y1, m1 = cls._parse_yyyymm_period(period)
        t0 = np.datetime64(f"{y0:04d}-{m0:02d}-01")
        if m1 == 12:
            t1 = np.datetime64(f"{y1 + 1:04d}-01-01")
        else:
            t1 = np.datetime64(f"{y1:04d}-{m1 + 1:02d}-01")
        return da.sel(time=slice(t0, t1))

    @staticmethod
    def _to_daily_date_time(da: xr.DataArray) -> xr.DataArray:
        if "time" not in da.coords:
            raise KeyError("Missing time coordinate.")
        t = da["time"]

        if np.issubdtype(t.dtype, np.datetime64):
            t_daily = t.astype("datetime64[D]")
        else:
            t_daily = xr.DataArray(
                np.array([np.datetime64(str(v)[:10]) for v in t.values], dtype="datetime64[D]"),
                dims=t.dims, coords=t.coords, name=t.name
            )

        out = da.assign_coords(time=t_daily)

        # collapse duplicates after daily casting
        if out.indexes["time"].has_duplicates:
            out = out.groupby("time").mean("time")

        return out

    @staticmethod
    def _open_var(files: List[str], var: str) -> xr.DataArray:
        ds = xr.open_mfdataset(files, combine="by_coords")
        if var in ds:
            return ds[var]
        for k in ds.data_vars:
            if k.lower() == var.lower():
                return ds[k]
        raise KeyError(f"Variable '{var}' not found in files. Available={list(ds.data_vars)}")

    @staticmethod
    def _guess_lat_lon_names(da: xr.DataArray) -> Tuple[str, str]:
        lat_candidates = ["lat", "latitude", "LAT", "nlat"]
        lon_candidates = ["lon", "longitude", "LON", "nlon"]
        lat = next((n for n in lat_candidates if (n in da.dims or n in da.coords)), None)
        lon = next((n for n in lon_candidates if (n in da.dims or n in da.coords)), None)
        if lat is None or lon is None:
            raise ValueError(f"Cannot infer lat/lon. dims={da.dims} coords={list(da.coords)}")
        return lat, lon

    @staticmethod
    def _wrap_lon_0_360(da: xr.DataArray, lon_name: str) -> xr.DataArray:
        if lon_name not in da.coords:
            return da
        lon = da[lon_name]
        if not np.issubdtype(lon.dtype, np.number):
            return da
        lon2 = (lon % 360.0)
        return da.assign_coords({lon_name: lon2}).sortby(lon_name)

    @staticmethod
    def _subset_region(
        da: xr.DataArray,
        lat: str,
        lon: str,
        lat_range: Optional[Tuple[float, float]],
        lon_range: Optional[Tuple[float, float]],
    ) -> xr.DataArray:
        out = da
        if lat_range is not None:
            lo, hi = lat_range
            latv = out[lat]
            if latv.size > 1 and (latv[0] > latv[-1]):
                out = out.sel({lat: slice(hi, lo)})
            else:
                out = out.sel({lat: slice(lo, hi)})
        if lon_range is not None:
            lo, hi = lon_range
            out = out.sel({lon: slice(lo, hi)})
        return out

    @staticmethod
    def _same_1d_values(a: xr.DataArray, b: xr.DataArray, tol: float) -> bool:
        if a.size != b.size:
            return False
        av = np.asarray(a.values)
        bv = np.asarray(b.values)
        if tol == 0.0:
            return np.array_equal(av, bv)
        return np.allclose(av, bv, rtol=0.0, atol=tol, equal_nan=True)

    @staticmethod
    def _rename_dim_if_present(da: xr.DataArray, old: str, new: str) -> xr.DataArray:
        if old == new:
            return da
        if (new in da.dims) or (new in da.coords):
            return da
        rename_map = {}
        if old in da.dims:
            rename_map[old] = new
        if old in da.coords:
            rename_map[old] = new
        return da.rename(rename_map) if rename_map else da

    def _ensure_names_match(
        self,
        mod: xr.DataArray,
        obs: xr.DataArray,
        mod_lat: str,
        mod_lon: str,
        obs_lat: str,
        obs_lon: str,
    ) -> xr.DataArray:
        if (mod_lat == obs_lat) and (mod_lon == obs_lon):
            return mod
        lat_ok = self._same_1d_values(mod[mod_lat], obs[obs_lat], self.cfg.grid_tol)
        lon_ok = self._same_1d_values(mod[mod_lon], obs[obs_lon], self.cfg.grid_tol)
        if self.cfg.require_same_grid and not (lat_ok and lon_ok):
            raise ValueError("Model/obs grids differ; regrid explicitly or set require_same_grid=False.")
        mod = self._rename_dim_if_present(mod, mod_lat, obs_lat)
        mod = self._rename_dim_if_present(mod, mod_lon, obs_lon)
        return mod

    def _align_and_check(self, mod: xr.DataArray, obs: xr.DataArray) -> tuple[xr.DataArray, xr.DataArray]:
        mod, obs = xr.align(mod, obs, join="inner")
        mod = mod.sortby("time")
        obs = obs.sortby("time")
        if self.cfg.require_nonempty_time and mod.sizes.get("time", 0) == 0:
            raise ValueError("No overlapping times between model and obs after alignment.")
        return mod, obs

    # ---------- kernels ----------
    @staticmethod
    def _maybe_convert_precip_to_mmday(da: xr.DataArray, *, var_key: str, verbose: bool = False) -> xr.DataArray:
        k = (var_key or "").lower()
        name = (da.name or "").lower()

        is_pr = (
            k in ("pr", "precip", "prect", "precc", "precl", "precsl", "tp")
            or "precip" in k or "precip" in name
            or k.startswith("pr") or name.startswith("pr")
            or name in ("prect", "precc", "precl", "precsl", "tp")
        )
        if not is_pr:
            return da

        u = (da.attrs.get("units") or "").strip().lower()
        u0 = u.replace(" ", "")

        def _set(out, note):
            out.attrs = dict(da.attrs)
            out.attrs["units"] = "mm/day"
            out.attrs["note_unit_fix"] = note
            return out

        if u0 in ("mm/day", "mm/d", "mmd-1", "mmdy-1") or ("mm" in u0 and ("day" in u0 or "d-1" in u0)):
            return da

        if ("kg" in u0 and "m-2" in u0 and "s-1" in u0) or u0 in ("kgm-2s-1", "kg/m2/s", "kgm**-2s**-1"):
            if verbose:
                print(f"[UNIT FIX] {da.name or var_key}: kg m-2 s-1 -> mm/day (×86400)")
            return _set(da * 86400.0, "converted from kg m-2 s-1 via ×86400")

        if u0 in ("m/s", "ms-1") or ("m" in u0 and "s-1" in u0 and "kg" not in u0 and "mm" not in u0):
            if verbose:
                print(f"[UNIT FIX] {da.name or var_key}: m/s -> mm/day (×86400×1000)")
            return _set(da * 86400.0 * 1000.0, "converted from m s-1 via ×86400×1000")

        if u0 in ("mm/s", "mms-1"):
            if verbose:
                print(f"[UNIT FIX] {da.name or var_key}: mm/s -> mm/day (×86400)")
            return _set(da * 86400.0, "converted from mm s-1 via ×86400")

        if u0 in ("m/day", "m/d", "md-1"):
            if verbose:
                print(f"[UNIT FIX] {da.name or var_key}: m/day -> mm/day (×1000)")
            return _set(da * 1000.0, "converted from m day-1 via ×1000")

        return da

    @staticmethod
    def _maybe_convert_geopotential_to_height(
        da: xr.DataArray,
        *,
        var_key: str,
        g: float = 9.80665,
        verbose: bool = False,
    ) -> xr.DataArray:
        k = (var_key or "").lower()
        name = (da.name or "").lower()

        looks_like_z = (
            k.startswith("z")
            or "geopot" in k
            or "geopot" in name
            or name.startswith("z")
            or name in ("zg", "phi")
        )
        if not looks_like_z:
            return da

        units = (da.attrs.get("units") or "").strip().lower()
        u0 = units.replace(" ", "").replace("^", "").replace("**", "")

        def is_geopot(u: str) -> bool:
            return ("m2s-2" in u) or ("m2/s2" in u)

        def is_height(u: str) -> bool:
            return u in ("m", "meter", "metre", "meters", "metres")

        if u0:
            if is_geopot(u0):
                if verbose:
                    print(f"[UNIT FIX] {da.name or var_key}: geopotential -> height (÷g)")
                out = da / g
                out.attrs = dict(da.attrs)
                out.attrs["units"] = "m"
                out.attrs["note_unit_fix"] = "converted from geopotential (m2 s-2) via /g"
                return out
            if is_height(u0):
                return da

        try:
            vmax = float(da.max(skipna=True))
        except Exception:
            return da

        if 2.0e4 < vmax < 5.0e6:
            if verbose:
                print(f"[UNIT FIX] {da.name or var_key}: heuristic geopotential -> height (÷g), vmax={vmax:.3g}")
            out = da / g
            out.attrs = dict(da.attrs)
            out.attrs["units"] = "m"
            out.attrs["note_unit_fix"] = "heuristic geopotential->height via /g (units missing/unknown)"
            return out

        return da

    def _maybe_fix_flux_sign_by_registry(
        self,
        da: xr.DataArray,
        *,
        var_key: str,
        source: str,   # "ERA5" / "MODEL" / "NOAA-OLR" ...
        registry: dict,
    ) -> xr.DataArray:
        v = (var_key or da.name or "").upper()

        s = registry.get((source, v), registry.get(v, +1))
        s = int(s)

        if s == -1:
            out = -da
            out.attrs = dict(da.attrs)
            out.attrs["note_sign_fix"] = f"flipped sign by registry: ({source},{v})"
            return out
        return da

    def _normalize_physics(self, da: xr.DataArray, *, var_key: str, source: str) -> xr.DataArray:
        """
        Normalize per-field physics once (no window loop cost):
          - lon to [0,360)
          - daily time coordinate
          - precip to mm/day (if needed)
          - geopotential -> height (if needed)
          - flux sign by registry (if needed)
        """
        lat, lon = self._guess_lat_lon_names(da)
        da = self._wrap_lon_0_360(da, lon)

        da = self._to_daily_date_time(da)

        da = self._maybe_convert_precip_to_mmday(da, var_key=var_key, verbose=self.cfg.verbose)
        da = self._maybe_convert_geopotential_to_height(da, var_key=var_key, verbose=self.cfg.verbose)

        da = self._maybe_fix_flux_sign_by_registry(
            da, var_key=var_key, source=source, registry=FLUX_SIGN_REGISTRY
        )
        return da

    @staticmethod
    def _tcc_over_time(
        f: xr.DataArray,
        o: xr.DataArray,
        *,
        min_count: int = 2,
    ) -> xr.DataArray:
        valid = f.notnull() & o.notnull()
        f2 = f.where(valid)
        o2 = o.where(valid)

        fa = f2 - f2.mean("time", skipna=True)
        oa = o2 - o2.mean("time", skipna=True)

        num = (fa * oa).sum("time", skipna=True)
        den = xr.ufuncs.sqrt(
            (fa * fa).sum("time", skipna=True)
            * (oa * oa).sum("time", skipna=True)
        )

        n = valid.sum("time")
        out = (num / den).where((den > 0) & (n >= min_count))
        return out

    @staticmethod
    def _rmse_over_time(f: xr.DataArray, o: xr.DataArray, *, min_count: int = 1) -> xr.DataArray:
        d2 = (f - o) ** 2
        mse = d2.mean("time", skipna=True)
        if min_count > 1:
            n = d2.notnull().sum("time")
            mse = mse.where(n >= min_count)
        return xr.ufuncs.sqrt(mse)

    # NEW: Bias = mean(f - o)
    @staticmethod
    def _bias_over_time(f: xr.DataArray, o: xr.DataArray, *, min_count: int = 1) -> xr.DataArray:
        valid = f.notnull() & o.notnull()
        d = (f - o).where(valid)

        bias = d.mean("time", skipna=True)

        if min_count > 1:
            n = valid.sum("time")
            bias = bias.where(n >= min_count)

        return bias

    # ---------- filenames ----------

    def base_path(
        self,
        *,
        group: str,
        freq: str,
        run: str,
        obs: str,
        period: str,
        region: str,
        exp: str,
        var: str,
        ens: str,
    ) -> str:
        fname = f"s2s_base_tcc_rmse_{group}_{freq}_{run}_{obs}_{period}_{region}_{exp}_{ens}_{var}.nc"
        return os.path.join(self.cfg.out_dir, fname)

    def ensmean_path(
        self,
        *,
        group: str,
        freq: str,
        run: str,
        obs: str,
        period: str,
        region: str,
        exp: str,
        var: str,
    ) -> str:
        fname = f"s2s_base_tcc_rmse_{group}_{freq}_{run}_{obs}_{period}_{region}_{exp}_ENSMEAN_{var}.nc"
        return os.path.join(self.cfg.out_dir, fname)

    @staticmethod
    def _print_da_extent(tag: str, da: xr.DataArray):
        def _rng(x):
            try:
                return float(x.min()), float(x.max()), x.size
            except Exception:
                return None

        msg = [f"[{tag}]"]

        if "time" in da.coords:
            t = da["time"]
            msg.append(f"time={str(t.min().values)[:10]} → {str(t.max().values)[:10]} (n={t.size})")

        for dim in ("lat", "latitude"):
            if dim in da.coords:
                lo, hi, n = _rng(da[dim])
                msg.append(f"{dim}={lo:.2f} → {hi:.2f} (n={n})")
                break

        for dim in ("lon", "longitude"):
            if dim in da.coords:
                lo, hi, n = _rng(da[dim])
                msg.append(f"{dim}={lo:.2f} → {hi:.2f} (n={n})")
                break

        print(" | ".join(msg))

    # ---------- main: build per member ----------
    def build_member_file(
        self,
        all_files: Dict[str, dict],
        *,
        group: str,
        freq: str,
        run: str,
        obs_name: str,
        period: str,
        var: str,
        exp: str,
        ens: str,
        obs_var_override: Optional[Dict[str, str]] = None,
    ) -> str:
        """
        Build ONE base file for a single (exp, var, ens) across all configured regions.
        Saves: TCC, RMSE, BIAS.
        """
        if var not in all_files or exp not in all_files[var]:
            raise KeyError(f"Missing var/exp in all_files: var={var} exp={exp}")

        regions = self.cfg.regions or [dict(name="DEFAULT", lat=self.cfg.lat_range, lon=self.cfg.lon_range)]

        # open obs once
        obs_files = all_files[var][exp]["obs"]
        obs_var = (obs_var_override or {}).get(var, var)
        obs0 = self._open_var(obs_files, obs_var)
        obs0 = self._select_period_time(obs0, period)
        obs0 = self._normalize_physics(obs0, var_key=var, source=obs_name)
        if self.cfg.verbose:
            self._print_da_extent(f"OBS normalized ({var})", obs0)

        lat_obs, lon_obs = self._guess_lat_lon_names(obs0)

        # open member once
        model_files = all_files[var][exp]["model"].get(ens, [])
        if not model_files:
            raise FileNotFoundError(f"No model files for exp={exp} var={var} ens={ens}")

        mem = self._open_var(model_files, var)
        mem = self._select_period_time(mem, period)
        mem = self._normalize_physics(mem, var_key=var, source="MODEL")
        if self.cfg.verbose:
            self._print_da_extent(f"MODEL normalized ({exp} {ens} {var})", mem)

        lat_m, lon_m = self._guess_lat_lon_names(mem)

        # per-region compute and pack
        reg_ds_list = []
        for reg in regions:
            rname = str(reg.get("name", "REGION"))
            rlat = reg.get("lat", None)
            rlon = reg.get("lon", None)

            obs_reg = self._subset_region(obs0, lat_obs, lon_obs, rlat, rlon)
            mem_reg = self._subset_region(mem, lat_m, lon_m, rlat, rlon)
            mem_reg = self._ensure_names_match(mem_reg, obs_reg, lat_m, lon_m, lat_obs, lon_obs)
            mem_a, obs_a = self._align_and_check(mem_reg, obs_reg)

            if self.cfg.verbose:
                self._print_da_extent(f"ALIGNED OBS region={rname}", obs_a)
                self._print_da_extent(f"ALIGNED MODEL region={rname}", mem_a)

            wlabels, wstart, wend = [], [], []
            tcc_list, rmse_list, bias_list = [], [], []  # NEW bias_list

            for wlab in sorted(self.cfg.windows.keys()):
                s, e = self.cfg.windows[wlab]
                ts = slice(np.datetime64(s), np.datetime64(e))
                f = mem_a.sel(time=ts)
                o = obs_a.sel(time=ts)
                if self.cfg.verbose:
                    self._print_da_extent(f"WINDOW {wlab} OBS", o)
                    self._print_da_extent(f"WINDOW {wlab} MODEL", f)
                if f.sizes.get("time", 0) == 0:
                    continue

                tcc_list.append(self._tcc_over_time(f, o))
                rmse_list.append(self._rmse_over_time(f, o))
                bias_list.append(self._bias_over_time(f, o))  # NEW

                wlabels.append(wlab)
                wstart.append(np.datetime64(s))
                wend.append(np.datetime64(e))

            if not tcc_list:
                if self.cfg.verbose:
                    print(f"[SKIP] no windows produced for {exp} {ens} {var} region={rname}")
                continue

            TCC = xr.concat(tcc_list, dim=xr.IndexVariable("window", wlabels))
            RMSE = xr.concat(rmse_list, dim=xr.IndexVariable("window", wlabels))
            BIAS = xr.concat(bias_list, dim=xr.IndexVariable("window", wlabels))  # NEW

            TCC = TCC.assign_coords(
                window_start=("window", np.array(wstart, dtype="datetime64[D]")),
                window_end=("window", np.array(wend, dtype="datetime64[D]")),
            )
            RMSE = RMSE.assign_coords(window_start=TCC["window_start"], window_end=TCC["window_end"])
            BIAS = BIAS.assign_coords(window_start=TCC["window_start"], window_end=TCC["window_end"])  # NEW

            ds_reg = xr.Dataset(
                data_vars=dict(
                    TCC=TCC,
                    RMSE=RMSE,
                    BIAS=BIAS,  # NEW
                ),
                coords=dict(window=TCC["window"]),
                attrs=dict(region=rname),
            ).assign_coords(region=rname).expand_dims("region")

            reg_ds_list.append(ds_reg)

        if not reg_ds_list:
            raise ValueError(f"No region results produced for exp={exp} ens={ens} var={var}")

        ds_out = xr.concat(reg_ds_list, dim="region")
        ds_out.attrs.update(
            dict(
                group=group,
                freq=freq,
                run=run,
                obs=obs_name,
                period=period,
                exp=exp,
                ens=ens,
                var=var,
                note="Base file: windowed TCC/RMSE/BIAS maps for a single ensemble member.",
            )
        )

        path = self.base_path(
            group=group, freq=freq, run=run, obs=obs_name, period=period,
            region="ALL", exp=exp, var=var, ens=ens
        )
        if (not self.cfg.overwrite) and os.path.exists(path):
            if self.cfg.verbose:
                print(f"[SKIP EXISTS] {path}")
            return path

        ds_out.to_netcdf(path)
        if self.cfg.verbose:
            print(f"[SAVED] {path}")
        return path

    def build_ensmean_file(
        self,
        all_files: Dict[str, dict],
        *,
        group: str,
        freq: str,
        run: str,
        obs_name: str,
        period: str,
        var: str,
        exp: str,
        obs_var_override: Optional[Dict[str, str]] = None,
    ) -> str:
        """
        Build ONE base file for ensemble-mean-based windowed TCC/RMSE/BIAS maps: (exp, var).
        """
        if var not in all_files or exp not in all_files[var]:
            raise KeyError(f"Missing var/exp in all_files: var={var} exp={exp}")

        regions = self.cfg.regions or [dict(name="DEFAULT", lat=self.cfg.lat_range, lon=self.cfg.lon_range)]

        # ---- open obs once ----
        obs_files = all_files[var][exp]["obs"]
        obs_var = (obs_var_override or {}).get(var, var)

        obs0 = self._open_var(obs_files, obs_var)
        obs0 = self._select_period_time(obs0, period)
        obs0 = self._normalize_physics(obs0, var_key=var, source=obs_name)
        if self.cfg.verbose:
            self._print_da_extent(f"OBS normalized ({var})", obs0)

        lat_obs, lon_obs = self._guess_lat_lon_names(obs0)

        # ---- open + build ensemble mean once ----
        model_by_ens = all_files[var][exp]["model"]
        members = []
        ens_names = []
        for ens, files in model_by_ens.items():
            if not files:
                continue
            da = self._open_var(files, var)
            da = self._select_period_time(da, period)
            da = self._normalize_physics(da, var_key=var, source="MODEL")
            members.append(da)
            ens_names.append(ens)

        if not members:
            raise ValueError(f"No member data for exp={exp} var={var}")

        mod_ens = xr.concat(members, dim=xr.IndexVariable("ens", ens_names))
        mod_mean = mod_ens.mean("ens")
        lat_m, lon_m = self._guess_lat_lon_names(mod_mean)

        # ---- per-region compute ----
        reg_ds_list = []
        for reg in regions:
            rname = str(reg.get("name", "REGION"))
            rlat = reg.get("lat", None)
            rlon = reg.get("lon", None)

            obs_reg = self._subset_region(obs0, lat_obs, lon_obs, rlat, rlon)
            mod_reg = self._subset_region(mod_mean, lat_m, lon_m, rlat, rlon)

            mod_reg = self._ensure_names_match(mod_reg, obs_reg, lat_m, lon_m, lat_obs, lon_obs)
            mod_a, obs_a = self._align_and_check(mod_reg, obs_reg)

            wlabels, wstart, wend = [], [], []
            tcc_list, rmse_list, bias_list = [], [], []  # NEW bias_list

            for wlab in sorted(self.cfg.windows.keys()):
                s, e = self.cfg.windows[wlab]
                ts = slice(np.datetime64(s), np.datetime64(e))
                f = mod_a.sel(time=ts)
                o = obs_a.sel(time=ts)
                if f.sizes.get("time", 0) == 0:
                    continue

                tcc_list.append(self._tcc_over_time(f, o))
                rmse_list.append(self._rmse_over_time(f, o))
                bias_list.append(self._bias_over_time(f, o))  # NEW

                wlabels.append(wlab)
                wstart.append(np.datetime64(s))
                wend.append(np.datetime64(e))

            if not tcc_list:
                if self.cfg.verbose:
                    print(f"[SKIP] ensmean: no windows produced for {exp} {var} region={rname}")
                continue

            TCC = xr.concat(tcc_list, dim=xr.IndexVariable("window", wlabels))
            RMSE = xr.concat(rmse_list, dim=xr.IndexVariable("window", wlabels))
            BIAS = xr.concat(bias_list, dim=xr.IndexVariable("window", wlabels))  # NEW

            TCC = TCC.assign_coords(
                window_start=("window", np.array(wstart, dtype="datetime64[D]")),
                window_end=("window", np.array(wend, dtype="datetime64[D]")),
            )
            RMSE = RMSE.assign_coords(window_start=TCC["window_start"], window_end=TCC["window_end"])
            BIAS = BIAS.assign_coords(window_start=TCC["window_start"], window_end=TCC["window_end"])  # NEW

            ds_reg = xr.Dataset(
                data_vars=dict(
                    TCC_ensmean=TCC,
                    RMSE_ensmean=RMSE,
                    BIAS_ensmean=BIAS,  # NEW
                ),
                coords=dict(window=TCC["window"]),
                attrs=dict(region=rname),
            ).assign_coords(region=rname).expand_dims("region")

            reg_ds_list.append(ds_reg)

        if not reg_ds_list:
            raise ValueError(f"No ensmean region results produced for exp={exp} var={var}")

        ds_out = xr.concat(reg_ds_list, dim="region")
        ds_out.attrs.update(
            dict(
                group=group, freq=freq, run=run, obs=obs_name, period=period,
                exp=exp, var=var,
                note="Ensemble-mean-based base file: windowed TCC/RMSE/BIAS maps computed from ensemble-mean time series.",
            )
        )

        path = self.ensmean_path(
            group=group, freq=freq, run=run, obs=obs_name, period=period,
            region="ALL", exp=exp, var=var
        )
        if (not self.cfg.overwrite) and os.path.exists(path):
            if self.cfg.verbose:
                print(f"[SKIP EXISTS] {path}")
            return path

        ds_out.to_netcdf(path)
        if self.cfg.verbose:
            print(f"[SAVED] {path}")
        return path


# ============================================================
# Stage 2: aggregator (UPDATED to aggregate BIAS)
# ============================================================

class S2SWindowedTCCAggregator:
    """
    Stage 2: aggregate member base files into ensemble summary for each (exp,var).
    Produces:
      - TCC_mean/std/quantile over ens
      - RMSE_mean/std/quantile over ens
      - BIAS_mean/std/quantile over ens   <-- NEW
      - optional: *_ensmean from ensmean base file
    """

    def __init__(self, cfg: WindowedTCCAggregateConfig):
        self.cfg = cfg

    @staticmethod
    def _open_many(paths: List[str]) -> xr.Dataset:
        return xr.open_mfdataset(paths, combine="by_coords")

    def aggregate_exp_var(
        self,
        member_paths: List[str],
        *,
        out_nc: str,
        ensmean_path: str | None = None,
    ) -> str:

        ds_list = []
        ens_names = []
        for p in member_paths:
            ds = xr.open_dataset(p)
            ens = str(ds.attrs.get("ens", "ens"))
            ds_list.append(ds.load())
            ds.close()
            ens_names.append(ens)

        ds_ens = xr.concat(ds_list, dim=xr.IndexVariable("ens", ens_names))

        out = xr.Dataset()
        out["TCC_member"] = ds_ens["TCC"]
        out["RMSE_member"] = ds_ens["RMSE"]
        out["BIAS_member"] = ds_ens["BIAS"]  # NEW

        out["TCC_mean"] = out["TCC_member"].mean("ens")
        out["RMSE_mean"] = out["RMSE_member"].mean("ens")
        out["BIAS_mean"] = out["BIAS_member"].mean("ens")  # NEW

        out["TCC_std"] = out["TCC_member"].std("ens")
        out["RMSE_std"] = out["RMSE_member"].std("ens")
        out["BIAS_std"] = out["BIAS_member"].std("ens")  # NEW

        qs = list(self.cfg.member_quantiles or ())
        if qs:
            out["TCC_quantile"] = nanquantile(out["TCC_member"], qs, "ens").rename({"quantile": "q"})
            out["RMSE_quantile"] = nanquantile(out["RMSE_member"], qs, "ens").rename({"quantile": "q"})
            out["BIAS_quantile"] = nanquantile(out["BIAS_member"], qs, "ens").rename({"quantile": "q"})  # NEW

        # add ensmean diagnostics if provided
        if ensmean_path is not None and os.path.exists(ensmean_path):
            ds_em = xr.open_dataset(ensmean_path)
            if "TCC_ensmean" in ds_em:
                out["TCC_ensmean"] = ds_em["TCC_ensmean"]
            if "RMSE_ensmean" in ds_em:
                out["RMSE_ensmean"] = ds_em["RMSE_ensmean"]
            if "BIAS_ensmean" in ds_em:  # NEW
                out["BIAS_ensmean"] = ds_em["BIAS_ensmean"]
            ds_em.close()

        # attrs
        for k, v in ds_list[0].attrs.items():
            out.attrs[k] = v
        out.attrs["note"] = (
            "Ensemble aggregation from per-member base files "
            "(mean/std/quantiles) + optional ensmean-based diagnostics. "
            "Includes TCC/RMSE/BIAS."
        )

        if (not self.cfg.overwrite) and os.path.exists(out_nc):
            raise FileExistsError(f"Output exists: {out_nc}")

        out.to_netcdf(out_nc)
        if self.cfg.verbose:
            print(f"[SAVED] {out_nc}")
        return out_nc

class RegionSpec:
    name: str
    lat: tuple[float, float]
    lon: tuple[float, float]
    mask: str = "none"   # "none" | "land" | "ocean"


class RegionalReduceConfig:
    out_dir: Path
    overwrite: bool = False
    verbose: bool = False

    # ens naming (members)
    ens_start: int = 1
    ens_prefix: str = "EN"
    ens_width: int = 2

    # ensemble-mean (skill of ensemble-mean forecast field)
    include_ensmean: bool = True
    ensmean_tag: str = "ENSMEAN"
    compute_ensmean_if_missing: bool = False  # keep False for your case

    # coords
    lat_names: tuple[str, ...] = ("lat", "latitude", "nlat")
    lon_names: tuple[str, ...] = ("lon", "longitude", "nlon")

    # mask
    landmask_file: Path | None = None
    landmask_var: str = "landfrac"
    land_threshold: float = 0.5

    # quantiles for member distribution (e.g., for shading)
    member_quantiles: tuple[float, ...] = (0.10, 0.25, 0.50, 0.75, 0.90)

    # Optional: require at least this fraction of region weights to be valid; otherwise set mean to NaN
    min_valid_frac: float = 0.0  # set to e.g. 0.05 if desired

    # NEW: save both quantile definitions
    save_member_quantiles_post: bool = True  # quantile(regional_mean over ens)
    save_member_quantiles_pre: bool = True   # regional_mean(quantile map over ens)


class S2SRegionalSkillReducer:
    """
    Reads per-member base files (and optional ENSMEAN base file) containing spatial skill maps
    (e.g., TCC/RMSE) on the model grid and computes regional area-weighted means.

    ENSMEAN is treated as a separate product: "skill of ensemble-mean forecast field",
    and is NOT included in member ensemble statistics.

    Quantiles:
      - qpostXX: quantile across ens of the *regional-mean* member skill (current behavior)
      - qpreXX : regional-mean of the *gridpoint* ensemble quantile map (new behavior)

    Key robustness:
      - Canonicalizes longitude to 0..360 for BOTH data and landmask
      - Region selection uses wrap-aware logic and treats (0,360) as full globe
      - Handles NaNs by computing weighted mean only over finite values and
        outputs valid-weight sums + valid fractions (optionally masks low coverage)
      - Removes input 'region' dim/coord before creating output region dim
    """

    def __init__(self, cfg: RegionalReduceConfig, regions: tuple[RegionSpec, ...]):
        self.cfg = cfg
        self.regions = regions
        self.cfg.out_dir.mkdir(parents=True, exist_ok=True)

    # ---------- utilities ----------
    def _pick_coord(self, ds: xr.Dataset, names: tuple[str, ...], kind: str) -> str:
        for n in names:
            if n in ds.coords or n in ds.dims:
                return n
        raise KeyError(
            f"Cannot find {kind} coord/dim. Tried {names}. "
            f"Found coords={list(ds.coords)}, dims={list(ds.dims)}"
        )

    def _lon_to_0_360_1d(self, lon: xr.DataArray) -> xr.DataArray:
        """Convert 1D lon coordinate to [0, 360) and ensure increasing order."""
        lon2 = lon.copy()
        if float(lon2.min()) < 0.0:
            lon2 = lon2 % 360.0

        if lon2.ndim == 1:
            vals = lon2.values
            if not np.all(np.diff(vals) >= 0):
                order = np.argsort(vals)
                lon2 = lon2.isel({lon2.dims[0]: order})
        return lon2

    def _canonicalize_lon_ds(self, ds: xr.Dataset, lon_name: str) -> xr.Dataset:
        """Canonicalize ds lon coordinate to 0..360 and sort ds by lon if necessary. Assumes lon is 1D."""
        if lon_name not in ds.coords and lon_name not in ds.dims:
            return ds

        lon_old = ds[lon_name]
        lon_new = self._lon_to_0_360_1d(lon_old)
        ds2 = ds.assign_coords({lon_name: lon_new})
        ds2 = ds2.sortby(lon_name)
        return ds2

    def _normalize_region_lon_0_360(self, lon0: float, lon1: float) -> tuple[float, float, bool]:
        """Normalize region lon bounds to 0..360. Returns (lonmin, lonmax, is_full_globe)."""
        span = lon1 - lon0
        if span >= 360.0 - 1e-6:
            return 0.0, 360.0, True
        lonmin = lon0 % 360.0
        lonmax = lon1 % 360.0
        return lonmin, lonmax, False

    def _clean_input_region(self, ds: xr.Dataset) -> xr.Dataset:
        """
        Clean the input 'region' axis from base files:
          - if region dim exists and is singleton: squeeze it (also removes index coord)
          - if region dim exists and is non-singleton: rename to region_in to avoid collision
          - drop any leftover non-index coord/var named 'region'
        """
        if "region" in ds.dims:
            n = ds.sizes.get("region", 0)
            if n == 1:
                ds = ds.squeeze("region", drop=True)
            else:
                ds = ds.rename({"region": "region_in"})

        if "region" in ds.coords:
            try:
                ds = ds.reset_coords("region", drop=True)
            except ValueError:
                ds = ds.rename({"region": "region_in_coord"})

        if "region" in ds.variables and "region" not in ds.dims:
            ds = ds.drop_vars("region", errors="ignore")

        return ds

    def _sanitize_region_for_output(self, ds: xr.Dataset) -> xr.Dataset:
        """Ensure there is NO 'region' dim/coord/var before creating OUTPUT region dim."""
        if "region" in ds.dims:
            n = ds.sizes.get("region", 0)
            if n == 1:
                ds = ds.squeeze("region", drop=True)
            else:
                ds = ds.rename({"region": "region_in"})

        if "region" in ds.coords:
            ds = ds.reset_coords("region", drop=True)

        if "region" in ds.variables and "region" not in ds.dims:
            ds = ds.drop_vars("region", errors="ignore")

        return ds

    def _squeeze_region_da(self, da: xr.DataArray) -> xr.DataArray:
        if "region" in da.dims:
            n = da.sizes.get("region", 0)
            if n == 1:
                da = da.squeeze("region", drop=True)
            else:
                raise ValueError(f"Unexpected non-singleton 'region' dim: {da.sizes['region']}")
        return da

    def _area_w2d(self, lat: xr.DataArray, lon: xr.DataArray) -> xr.DataArray:
        """Area weight ~ cos(lat); broadcast to (lat, lon)."""
        w1 = np.cos(np.deg2rad(lat))
        w1 = xr.where(w1 < 0, 0.0, w1)
        return w1.broadcast_like(
            xr.DataArray(
                np.zeros((lat.size, lon.size)),
                coords={lat.name: lat, lon.name: lon},
                dims=(lat.name, lon.name),
            )
        )

    def _ens_tag(self, k: int) -> str:
        return f"{self.cfg.ens_prefix}{k:0{self.cfg.ens_width}d}"

    def _canonical_metric_name(self, vn: str) -> str:
        base = vn
        for suf in ("_ensmean", "_member"):
            if base.endswith(suf):
                base = base[: -len(suf)]
        return base

    def _validate_quantiles(self) -> tuple[float, ...]:
        qs = tuple(sorted(set(self.cfg.member_quantiles)))
        for q in qs:
            if not (0.0 <= q <= 1.0):
                raise ValueError(f"member_quantiles must be in [0,1]. Got {q}")
        return qs

    # ---------- file collection ----------
    def collect_member_and_ensmean_files(
        self,
        input_dir: Path,
        *,
        group: str,
        freq: str,
        run: str,
        obs: str,
        period: str,
        exp: str,
        var: str,
        nens: int,
    ) -> tuple[list[Path], list[int], Path | None]:
        member_files: list[Path] = []
        member_ens: list[int] = []

        for k in range(self.cfg.ens_start, self.cfg.ens_start + nens):
            tag = self._ens_tag(k)
            f = input_dir / f"s2s_base_tcc_rmse_{group}_{freq}_{run}_{obs}_{period}_ALL_{exp}_{tag}_{var}.nc"
            if f.exists():
                member_files.append(f)
                member_ens.append(k)

        ensmean_file = None
        if self.cfg.include_ensmean:
            tag = self.cfg.ensmean_tag
            f = input_dir / f"s2s_base_tcc_rmse_{group}_{freq}_{run}_{obs}_{period}_ALL_{exp}_{tag}_{var}.nc"
            if f.exists():
                ensmean_file = f

        return member_files, member_ens, ensmean_file

    # ---------- region + landmask ----------
    def _region_bool(self, lat: xr.DataArray, lon: xr.DataArray, reg: RegionSpec) -> xr.DataArray:
        """Region mask on canonical lon (expected 0..360). Wrap-aware; (0,360) treated as full globe."""
        lat2d, lon2d = xr.broadcast(lat, lon)

        latmin, latmax = reg.lat
        in_lat = (lat2d >= latmin) & (lat2d <= latmax)

        lon0, lon1 = reg.lon
        lonmin, lonmax, full = self._normalize_region_lon_0_360(lon0, lon1)

        if full:
            in_lon = xr.ones_like(lon2d, dtype=bool)
        else:
            if lonmin <= lonmax:
                in_lon = (lon2d >= lonmin) & (lon2d <= lonmax)
            else:
                in_lon = (lon2d >= lonmin) | (lon2d <= lonmax)

        return (in_lat & in_lon).transpose(lat.name, lon.name)

    def _load_landmask_bool(self, target_lat: xr.DataArray, target_lon: xr.DataArray, want: str) -> xr.DataArray:
        """want: 'none' | 'land' | 'ocean'"""
        if want == "none":
            return xr.ones_like(self._area_w2d(target_lat, target_lon), dtype=bool)

        if self.cfg.landmask_file is None:
            raise ValueError("Region requests land/ocean mask but cfg.landmask_file is None")

        ds_m = xr.open_dataset(self.cfg.landmask_file, decode_times=False)
        if self.cfg.landmask_var not in ds_m:
            raise KeyError(
                f"landmask_var={self.cfg.landmask_var} not in {self.cfg.landmask_file}. "
                f"Found {list(ds_m.data_vars)}"
            )

        latm = self._pick_coord(ds_m, self.cfg.lat_names, "lat")
        lonm = self._pick_coord(ds_m, self.cfg.lon_names, "lon")

        ds_m = self._canonicalize_lon_ds(ds_m, lonm)
        m = ds_m[self.cfg.landmask_var]

        target_lon = self._lon_to_0_360_1d(target_lon)

        if (m.sizes.get(latm) != target_lat.size) or (m.sizes.get(lonm) != target_lon.size):
            m = m.interp({latm: target_lat, lonm: target_lon}, method="nearest")

        is_land = (m >= self.cfg.land_threshold)
        out = is_land if want == "land" else (~is_land)

        out = out.rename({latm: target_lat.name, lonm: target_lon.name})
        return out.astype(bool).transpose(target_lat.name, target_lon.name)

    # ---------- math helpers ----------
    def _weighted_mean_latlon(self, da: xr.DataArray, w2d: xr.DataArray, lat_name: str, lon_name: str) -> xr.DataArray:
        """Weighted mean over lat/lon using weights w2d, ignoring NaNs in da."""
        w = w2d
        for d in da.dims:
            if d not in w.dims:
                w = w.expand_dims({d: da.sizes[d]})
        w = w.transpose(*da.dims)

        valid = xr.ufuncs.isfinite(da)
        w_eff = w.where(valid, 0.0)

        num = (da.where(valid) * w_eff).sum(dim=(lat_name, lon_name), skipna=True)
        den = w_eff.sum(dim=(lat_name, lon_name), skipna=True)

        out = num / den
        out = out.where(den > 0)
        return out

    def _valid_weight_sum(self, da: xr.DataArray, w2d: xr.DataArray, lat_name: str, lon_name: str) -> xr.DataArray:
        """Sum of weights in region where da is finite."""
        w = w2d
        for d in da.dims:
            if d not in w.dims:
                w = w.expand_dims({d: da.sizes[d]})
        w = w.transpose(*da.dims)

        valid = xr.ufuncs.isfinite(da)
        return w.where(valid, 0.0).sum(dim=(lat_name, lon_name), skipna=True)

    # ---------- open helpers ----------
    def _open_stack_members(self, member_files: list[Path], member_ens: list[int]) -> xr.Dataset:
        dsets = []
        for f, ensv in zip(member_files, member_ens):
            ds = xr.open_dataset(f)
            ds = self._clean_input_region(ds)
            ds = ds.expand_dims({"ens": [ensv]})
            dsets.append(ds)
        return xr.concat(dsets, dim="ens")

    def _open_ensmean(self, ensmean_file: Path) -> xr.Dataset:
        ds = xr.open_dataset(ensmean_file)
        ds = self._clean_input_region(ds)
        return ds

    # ---------- main API ----------
    def reduce_one_exp_var(
        self,
        input_dir: Path,
        *,
        group: str,
        freq: str,
        run: str,
        obs: str,
        period: str,
        exp: str,
        var: str,
        nens: int,
    ) -> Path | None:
        member_files, member_ens, ensmean_file = self.collect_member_and_ensmean_files(
            input_dir,
            group=group, freq=freq, run=run, obs=obs, period=period, exp=exp, var=var, nens=nens
        )

        if not member_files and ensmean_file is None:
            if self.cfg.verbose:
                print(f"[miss] {group} {obs} {exp} {var}: no member or ensmean files")
            return None

        return self.reduce_members_to_one_output(
            member_files=member_files,
            member_ens=member_ens,
            ensmean_file=ensmean_file,
            group=group, obs=obs, period=period, exp=exp, var=var,
        )

    def reduce_members_to_one_output(
        self,
        member_files: list[Path],
        member_ens: list[int],
        ensmean_file: Path | None,
        *,
        group: str,
        obs: str,
        period: str,
        exp: str,
        var: str,
    ) -> Path:
        out_name = f"s2s_regional_tcc_rmse_{group}_{obs}_{period}_{exp}_{var}.nc"
        out_path = self.cfg.out_dir / out_name

        if out_path.exists() and not self.cfg.overwrite:
            if self.cfg.verbose:
                print(f"[skip] exists: {out_path}")
            return out_path

        qs = self._validate_quantiles()

        ds_mem = self._open_stack_members(member_files, member_ens) if member_files else None
        ds_em = self._open_ensmean(ensmean_file) if ensmean_file is not None else None

        ds_ref = ds_mem if ds_mem is not None else ds_em
        if ds_ref is None:
            raise RuntimeError("Internal: no reference dataset available.")

        lat_name = self._pick_coord(ds_ref, self.cfg.lat_names, "lat")
        lon_name = self._pick_coord(ds_ref, self.cfg.lon_names, "lon")

        if ds_mem is not None:
            ds_mem = self._canonicalize_lon_ds(ds_mem, lon_name)
        if ds_em is not None:
            ds_em = self._canonicalize_lon_ds(ds_em, lon_name)

        ds_ref2 = ds_mem if ds_mem is not None else ds_em
        if ds_ref2 is None:
            raise RuntimeError("Internal: no reference dataset available after lon canonicalization.")

        lat = ds_ref2[lat_name]
        lon = ds_ref2[lon_name]

        # map vars that are maps (have lat/lon)
        map_vars_mem = []
        if ds_mem is not None:
            map_vars_mem = [vn for vn, da in ds_mem.data_vars.items() if (lat_name in da.dims and lon_name in da.dims)]
        map_vars_em = []
        if ds_em is not None:
            map_vars_em = [vn for vn, da in ds_em.data_vars.items() if (lat_name in da.dims and lon_name in da.dims)]

        mem_by_metric: dict[str, list[str]] = {}
        for vn in map_vars_mem:
            mem_by_metric.setdefault(self._canonical_metric_name(vn), []).append(vn)

        em_by_metric: dict[str, list[str]] = {}
        for vn in map_vars_em:
            em_by_metric.setdefault(self._canonical_metric_name(vn), []).append(vn)

        metrics = sorted(set(mem_by_metric.keys()) | set(em_by_metric.keys()))
        if self.cfg.verbose:
            print(f"[reduce] {group} {obs} {exp} {var}: metrics={metrics}")
            if ds_em is not None:
                print(f"         ensmean map vars: {map_vars_em}")

        w_area = self._area_w2d(lat, lon)

        lo_cache: dict[str, xr.DataArray] = {"none": xr.ones_like(w_area, dtype=bool)}
        if any(r.mask == "land" for r in self.regions):
            lo_cache["land"] = self._load_landmask_bool(lat, lon, "land")
        if any(r.mask == "ocean" for r in self.regions):
            lo_cache["ocean"] = self._load_landmask_bool(lat, lon, "ocean")

        out_list: list[xr.Dataset] = []
        for reg in self.regions:
            m_reg = self._region_bool(lat, lon, reg)
            m_lo = lo_cache.get(reg.mask)
            if m_lo is None:
                m_lo = self._load_landmask_bool(lat, lon, reg.mask)
            m = (m_reg & m_lo)

            w = w_area.where(m, 0.0)

            reg_ds = xr.Dataset()

            region_wsum = w.sum(dim=(lat_name, lon_name), skipna=True)
            region_wsum = self._squeeze_region_da(region_wsum)
            reg_ds["region_area_weight_sum"] = region_wsum

            if ds_mem is not None:
                for metric in metrics:
                    vns = mem_by_metric.get(metric, [])
                    if not vns:
                        continue
                    vn = vns[0]  # assume one map var per metric in member dataset

                    reg_mean_mem = self._weighted_mean_latlon(ds_mem[vn], w, lat_name, lon_name)
                    reg_mean_mem = self._squeeze_region_da(reg_mean_mem)

                    vw_mem = self._valid_weight_sum(ds_mem[vn], w, lat_name, lon_name)
                    vw_mem = self._squeeze_region_da(vw_mem)

                    valid_frac_mem = (vw_mem / region_wsum).where(region_wsum > 0)

                    if self.cfg.min_valid_frac > 0:
                        reg_mean_mem = reg_mean_mem.where(valid_frac_mem >= self.cfg.min_valid_frac)

                    reg_ds[f"{metric}_regional_member"] = reg_mean_mem
                    reg_ds[f"{metric}_valid_wsum_member"] = vw_mem
                    reg_ds[f"{metric}_valid_frac_member"] = valid_frac_mem

                    reg_ds[f"{metric}_member_ens_mean"] = reg_mean_mem.mean("ens", skipna=True)
                    reg_ds[f"{metric}_member_ens_std"] = reg_mean_mem.std("ens", skipna=True)

                    # ---- Quantiles (POST): quantile of regional mean across members ----
                    if self.cfg.save_member_quantiles_post:
                        with np.errstate(all="ignore"):
                            for q in qs:
                                qlab = int(round(q * 100))
                                reg_ds[f"{metric}_member_qpost{qlab:02d}"] = reg_mean_mem.quantile(
                                    q, dim="ens", skipna=True
                                )

                    # ---- Quantiles (PRE): regional mean of gridpoint quantile maps ----
                    if self.cfg.save_member_quantiles_pre:
                        # qmap dims: (q, ... , lat, lon)  (q coordinate name: "q")
                        qmap = nanquantile(ds_mem[vn], qs, "ens").rename({"quantile": "q"})
                        qreg = self._weighted_mean_latlon(qmap, w, lat_name, lon_name)
                        qreg = self._squeeze_region_da(qreg)

                        with np.errstate(all="ignore"):
                            for q in qs:
                                qlab = int(round(q * 100))
                                reg_ds[f"{metric}_member_qpre{qlab:02d}"] = qreg.sel(q=q).drop_vars("q", errors="ignore")

            if ds_em is not None:
                for metric in metrics:
                    vns = em_by_metric.get(metric, [])
                    if not vns:
                        continue
                    vn = vns[0]  # assume one map var per metric in ensmean dataset

                    reg_mean_em = self._weighted_mean_latlon(ds_em[vn], w, lat_name, lon_name)
                    reg_mean_em = self._squeeze_region_da(reg_mean_em)

                    vw_em = self._valid_weight_sum(ds_em[vn], w, lat_name, lon_name)
                    vw_em = self._squeeze_region_da(vw_em)

                    valid_frac_em = (vw_em / region_wsum).where(region_wsum > 0)

                    if self.cfg.min_valid_frac > 0:
                        reg_mean_em = reg_mean_em.where(valid_frac_em >= self.cfg.min_valid_frac)

                    reg_ds[f"{metric}_regional_ensmean"] = reg_mean_em
                    reg_ds[f"{metric}_valid_wsum_ensmean"] = vw_em
                    reg_ds[f"{metric}_valid_frac_ensmean"] = valid_frac_em

            reg_ds = self._sanitize_region_for_output(reg_ds)
            reg_ds = reg_ds.expand_dims({"region": [reg.name]})
            reg_ds = reg_ds.assign_attrs(region_mask=reg.mask)
            out_list.append(reg_ds)

        ds_out = xr.concat(out_list, dim="region")
        ds_out = ds_out.assign_attrs(
            group=group,
            obs=obs,
            exp=exp,
            var=var,
            period=period,
            member_quantiles=",".join(f"{q:.3f}" for q in qs),
            min_valid_frac=str(self.cfg.min_valid_frac),
            source_member_files=";".join(str(p) for p in member_files) if member_files else "",
            source_ensmean_file=str(ensmean_file) if ensmean_file is not None else "",
            member_quantile_post_definition="qpostXX = quantile over ens of regional-mean member skill",
            member_quantile_pre_definition="qpreXX = regional-mean of gridpoint quantile map over ens",
        )

        ds_out.to_netcdf(out_path)
        return out_path
