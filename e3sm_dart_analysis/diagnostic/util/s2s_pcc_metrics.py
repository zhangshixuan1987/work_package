# === Core Python ===
import os
import glob
import re
import collections
import cftime
from datetime import datetime
from typing import Tuple, Dict, Optional, List

# === Numerical & Data Handling ===
import numpy as np
import numpy.ma as ma
import pandas as pd
import xarray as xr
import xcdat as xcd
import xskillscore as xs
from xskillscore import rmse, pearson_r

from dataclasses import dataclass

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
        # FIX: use self.ens_start (was ens_start undefined)
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
# Skill assessment class
# ============================================================
@dataclass(frozen=True)
class SkillMetricConfig:
    # default single-region behavior (backward compatible if regions=None)
    lat_range: Optional[Tuple[float, float]] = (-20.0, 20.0)
    lon_range: Optional[Tuple[float, float]] = None

    require_same_grid: bool = True
    grid_tol: float = 0.0
    require_nonempty_time: bool = True

    compute_members: bool = True
    member_quantiles: Tuple[float, ...] = (0.1, 0.5, 0.9)  # q10/median/q90

    # NEW: multi-region evaluation
    # Each region is a dict: {"name": str, "lat": (lo,hi), "lon": (lo,hi), "mask": Optional["land"]}
    regions: Optional[List[dict]] = None

    # NEW: landmask support (for regions with mask="land")
    landmask_file: Optional[str] = None
    landmask_var: Optional[str] = None   # if None, auto-detect
    land_threshold: float = 0.5

    out_nc: str = "s2s_skill_acc_rmse.nc"
    overwrite: bool = True
    verbose: bool = True


class S2SSkillAssessor:
    """
    Compute ACC and RMSE lead-time curves from S2SFileCollector output.

    Output dims:
      - ACC(var, region, exp, time)
      - RMSE(var, region, exp, time)
      - optional: ACC_member(var, region, exp, ens, time) and RMSE_member(...)

    Regions:
      If cfg.regions is None -> use a single region named "DEFAULT" using cfg.lat_range/lon_range.
      If cfg.regions is provided -> loop over those region definitions.
    """

    def __init__(self, cfg: SkillMetricConfig):
        self.cfg = cfg
        self._landmask_cache: dict = {}  # cache landmask interpolated onto obs grid

    # ----------------------------
    # Helpers: time/IO
    # ----------------------------

    @staticmethod
    def _parse_yyyymm_period(period: str) -> Tuple[int, int, int, int]:
        m = re.match(r"^(\d{4})(\d{2})-(\d{4})(\d{2})$", period)
        if not m:
            raise ValueError(f"Bad period '{period}', expected 'YYYYMM-YYYYMM'")
        y0, m0, y1, m1 = map(int, m.groups())
        if (y1, m1) < (y0, m0):
            raise ValueError(f"Bad period '{period}' (end < start)")
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
    def _open_var(files: List[str], var: str) -> xr.DataArray:
        ds = xr.open_mfdataset(files, combine="by_coords")
        if var in ds:
            return ds[var]
        for k in ds.data_vars:
            if k.lower() == var.lower():
                return ds[k]
        raise KeyError(f"Variable '{var}' not found. Available={list(ds.data_vars)}")

    def _open_member(self, files: List[str], *, var: str, period: str) -> xr.DataArray:
        da = self._open_var(files, var)
        da = self._select_period_time(da, period)
        da = self._to_daily_date_time(da)
        return da

    @staticmethod
    def _to_daily_date_time(obj: xr.DataArray | xr.Dataset) -> xr.DataArray | xr.Dataset:
        """Normalize time coordinate to daily datetime64[D] (floor to date)."""
        if "time" not in obj.coords:
            raise KeyError("Object has no 'time' coordinate to normalize.")

        t = obj["time"]
        if np.issubdtype(t.dtype, np.datetime64):
            t_daily = t.astype("datetime64[D]")
        else:
            t_daily = xr.DataArray(
                np.array([np.datetime64(str(v)[:10]) for v in t.values], dtype="datetime64[D]"),
                dims=t.dims,
                coords=t.coords,
                name=t.name,
            )
        return obj.assign_coords(time=t_daily)

    # ----------------------------
    # Helpers: coords / region / weights
    # ----------------------------

    @staticmethod
    def _guess_lat_lon_names(da: xr.DataArray) -> Tuple[str, str]:
        lat_candidates = ["lat", "latitude", "LAT", "nlat"]
        lon_candidates = ["lon", "longitude", "LON", "nlon"]
        lat = next((n for n in lat_candidates if (n in da.dims or n in da.coords)), None)
        lon = next((n for n in lon_candidates if (n in da.dims or n in da.coords)), None)
        if lat is None or lon is None:
            raise ValueError(f"Cannot infer lat/lon from dims={da.dims}, coords={list(da.coords)}")
        return lat, lon

    @staticmethod
    def _wrap_lon_0_360(da: xr.DataArray, lon_name: str) -> xr.DataArray:
        """Normalize lon coord to [0, 360) and sort by lon."""
        if lon_name not in da.coords:
            return da
        lon = da[lon_name]
        if not np.issubdtype(lon.dtype, np.number):
            return da
        lon2 = (lon % 360.0)
        out = da.assign_coords({lon_name: lon2}).sortby(lon_name)
        return out

    @staticmethod
    def _subset_region(
        da: xr.DataArray,
        lat_name: str,
        lon_name: str,
        lat_range: Optional[Tuple[float, float]],
        lon_range: Optional[Tuple[float, float]],
    ) -> xr.DataArray:
        out = da
        if lat_range is not None:
            lo, hi = lat_range
            latv = out[lat_name]
            if latv.size > 1 and (latv[0] > latv[-1]):
                out = out.sel({lat_name: slice(hi, lo)})
            else:
                out = out.sel({lat_name: slice(lo, hi)})

        if lon_range is not None:
            lo, hi = lon_range
            out = out.sel({lon_name: slice(lo, hi)})

        return out

    @staticmethod
    def _coslat_weights(lat: xr.DataArray) -> xr.DataArray:
        w = np.cos(np.deg2rad(lat))
        return xr.DataArray(w, coords={lat.dims[0]: lat}, dims=lat.dims).clip(min=0.0)

    @staticmethod
    def _weighted_mean_space(a: xr.DataArray, w2d: xr.DataArray) -> xr.DataArray:
        dims = list(w2d.dims)
        return (a * w2d).sum(dims) / w2d.sum(dims)

    # ----------------------------
    # Helpers: grid matching
    # ----------------------------

    @staticmethod
    def _same_1d_values(a: xr.DataArray, b: xr.DataArray, *, tol: float = 0.0) -> bool:
        if a.size != b.size:
            return False
        av = np.asarray(a.values)
        bv = np.asarray(b.values)
        if tol == 0.0:
            return np.array_equal(av, bv)
        return np.allclose(av, bv, rtol=0.0, atol=tol, equal_nan=True)

    @classmethod
    def _ensure_latlon_names_match(
        cls,
        mod: xr.DataArray,
        obs: xr.DataArray,
        *,
        mod_lat: str,
        mod_lon: str,
        obs_lat: str,
        obs_lon: str,
        tol: float = 0.0,
    ) -> xr.DataArray:
        """Rename model lat/lon dim/coord names to obs names if coord values match."""
        if (mod_lat == obs_lat) and (mod_lon == obs_lon):
            return mod

        if mod_lat not in mod.coords or mod_lon not in mod.coords:
            raise ValueError(f"Model missing coords '{mod_lat}'/'{mod_lon}'. coords={list(mod.coords)}")
        if obs_lat not in obs.coords or obs_lon not in obs.coords:
            raise ValueError(f"Obs missing coords '{obs_lat}'/'{obs_lon}'. coords={list(obs.coords)}")

        lat_ok = cls._same_1d_values(mod[mod_lat], obs[obs_lat], tol=tol)
        lon_ok = cls._same_1d_values(mod[mod_lon], obs[obs_lon], tol=tol)

        if not (lat_ok and lon_ok):
            raise ValueError(
                "Model/obs lat/lon coord values differ, cannot safely rename.\n"
                f"  model names: {mod_lat}, {mod_lon}\n"
                f"  obs names:   {obs_lat}, {obs_lon}\n"
                f"  model sizes: lat={mod[mod_lat].size} lon={mod[mod_lon].size}\n"
                f"  obs sizes:   lat={obs[obs_lat].size} lon={obs[obs_lon].size}\n"
                "You likely need explicit regridding."
            )

        mod = cls._rename_dim_if_present(mod, mod_lat, obs_lat)
        mod = cls._rename_dim_if_present(mod, mod_lon, obs_lon)
        return mod

    @staticmethod
    def _require_1d_coord(da: xr.DataArray, name: str) -> xr.DataArray:
        if name not in da.coords:
            raise ValueError(f"Missing coord '{name}' in {da.name or 'DataArray'} coords={list(da.coords)}")
        c = da[name]
        if c.ndim != 1:
            raise ValueError(f"Coord '{name}' is not 1D (ndim={c.ndim}). dims={c.dims}")
        return c

    def _align_and_check(
        self,
        mod: xr.DataArray,
        obs: xr.DataArray,
        *,
        lat_name: str,
        lon_name: str,
    ) -> tuple[xr.DataArray, xr.DataArray]:
        if "time" not in mod.coords or "time" not in obs.coords:
            raise ValueError("Both model and obs must have 'time' coordinate.")

        mod, obs = xr.align(mod, obs, join="inner")
        mod = mod.sortby("time")
        obs = obs.sortby("time")

        if self.cfg.require_nonempty_time and mod.sizes.get("time", 0) == 0:
            raise ValueError("No overlapping times between model and obs after alignment.")

        mlat = self._require_1d_coord(mod, lat_name)
        mlon = self._require_1d_coord(mod, lon_name)
        olat = self._require_1d_coord(obs, lat_name)
        olon = self._require_1d_coord(obs, lon_name)

        if self.cfg.require_same_grid:
            lat_ok = self._same_1d_values(mlat, olat, tol=self.cfg.grid_tol)
            lon_ok = self._same_1d_values(mlon, olon, tol=self.cfg.grid_tol)
            if not (lat_ok and lon_ok):
                raise ValueError(
                    "Model/obs grids differ (lat/lon not identical). "
                    "Regrid explicitly (recommended) or set require_same_grid=False.\n"
                    f"  lat sizes: mod={mlat.size} obs={olat.size}\n"
                    f"  lon sizes: mod={mlon.size} obs={olon.size}\n"
                    f"  lat[0:3]: mod={np.asarray(mlat.values[:3])} obs={np.asarray(olat.values[:3])}\n"
                    f"  lon[0:3]: mod={np.asarray(mlon.values[:3])} obs={np.asarray(olon.values[:3])}"
                )

        return mod, obs

    # ----------------------------
    # Helpers: metrics
    # ----------------------------

    @classmethod
    def _pattern_corr(cls, x: xr.DataArray, y: xr.DataArray, w2d: xr.DataArray) -> xr.DataArray:
        x, y = xr.align(x, y, join="inner")
        m = xr.ufuncs.isfinite(x) & xr.ufuncs.isfinite(y)
        x = x.where(m)
        y = y.where(m)

        mx = cls._weighted_mean_space(x, w2d)
        my = cls._weighted_mean_space(y, w2d)
        xa = x - mx
        ya = y - my

        cov = cls._weighted_mean_space(xa * ya, w2d)
        vx = cls._weighted_mean_space(xa * xa, w2d)
        vy = cls._weighted_mean_space(ya * ya, w2d)

        return cov / xr.ufuncs.sqrt(vx * vy)

    @classmethod
    def _weighted_rmse(cls, x: xr.DataArray, y: xr.DataArray, w2d: xr.DataArray) -> xr.DataArray:
        x, y = xr.align(x, y, join="inner")
        m = xr.ufuncs.isfinite(x) & xr.ufuncs.isfinite(y)
        d2 = (x - y) ** 2
        d2 = d2.where(m)
        mse = cls._weighted_mean_space(d2, w2d)
        return xr.ufuncs.sqrt(mse)

    @classmethod
    def _concat_ens_mean(
        cls,
        model_files_by_ens: Dict[str, List[str]],
        *,
        var: str,
        period: str,
    ) -> xr.DataArray:
        members: List[xr.DataArray] = []
        ens_names: List[str] = []
        for ens, files in model_files_by_ens.items():
            if not files:
                continue
            da = cls._open_var(files, var)
            da = cls._select_period_time(da, period)
            members.append(da)
            ens_names.append(ens)

        if not members:
            raise ValueError(f"No model files found for var={var} across ensembles.")

        ens_da = xr.concat(members, dim=xr.IndexVariable("ens", ens_names))
        return ens_da.mean("ens")

    # ----------------------------
    # Helpers: unit normalization
    # ----------------------------

    @staticmethod
    def _maybe_convert_geopotential_to_height(
        da: xr.DataArray,
        *,
        var_key: str,
        g: float = 9.80665,
        verbose: bool = False,
    ) -> xr.DataArray:
        """Convert geopotential (m2/s2) to geopotential height (m) when appropriate."""
        k = (var_key or "").lower()
        name = (da.name or "").lower()
        if not (
            k.startswith("z")
            or "geopot" in k
            or "geopot" in name
            or name.startswith("z")
            or name in ("zg", "phi")
        ):
            return da

        units = (da.attrs.get("units") or "").strip().lower().replace("^", "").replace("**", "")
        units = units.replace(" ", "")

        def _is_geopotential_units(u: str) -> bool:
            return ("m2s-2" in u) or ("m2/s2" in u) or ("m2s^-2" in u) or ("m2s**-2" in u)

        def _is_height_units(u: str) -> bool:
            return u in ("m", "meter", "metre", "meters", "metres")

        if units:
            if _is_geopotential_units(units):
                if verbose:
                    print(f"[UNIT FIX] {da.name or var_key}: geopotential -> height (÷g)")
                out = da / g
                out.attrs = dict(da.attrs)
                out.attrs["units"] = "m"
                out.attrs["note_unit_fix"] = "converted from geopotential (m2 s-2) via /g"
                return out
            if _is_height_units(units):
                return da

        # heuristic fallback
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
            out.attrs["note_unit_fix"] = "heuristic: treated as geopotential and converted via /g"
            return out

        return da

    @staticmethod
    def _maybe_convert_precip_to_mmday(
        da: xr.DataArray,
        *,
        var_key: str,
        verbose: bool = False,
    ) -> xr.DataArray:
        """Normalize precipitation-like variables to mm/day."""
        k = (var_key or "").lower()
        name = (da.name or "").lower()

        is_pr = (
            k in ("pr", "precip", "precsl", "prect", "precc", "precl", "tp")
            or "precip" in k
            or "precip" in name
            or k.startswith("pr")
            or name.startswith("pr")
            or name in ("prect", "precc", "precsl", "precl", "tp")
        )
        if not is_pr:
            return da

        u = (da.attrs.get("units") or "").strip().lower()
        u0 = u.replace(" ", "")

        def _set_units(out: xr.DataArray, note: str) -> xr.DataArray:
            out.attrs = dict(da.attrs)
            out.attrs["units"] = "mm/day"
            out.attrs["note_unit_fix"] = note
            return out

        if u0 in ("mm/day", "mmd-1", "mm/d", "mmdy-1") or ("mm" in u0 and ("day" in u0 or "d-1" in u0)):
            return da

        if ("kg" in u0 and "m-2" in u0 and "s-1" in u0) or u0 in ("kgm-2s-1", "kg/m2/s", "kgm**-2s**-1"):
            if verbose:
                print(f"[UNIT FIX] {da.name or var_key}: kg m-2 s-1 -> mm/day (×86400)")
            return _set_units(da * 86400.0, "converted from kg m-2 s-1 via ×86400")

        if u0 in ("m/s", "ms-1", "msec-1") or ("m" in u0 and "s-1" in u0 and "kg" not in u0 and "mm" not in u0):
            if verbose:
                print(f"[UNIT FIX] {da.name or var_key}: m/s -> mm/day (×86400×1000)")
            return _set_units(da * 86400.0 * 1000.0, "converted from m s-1 via ×86400×1000")

        if u0 in ("mm/s", "mms-1"):
            if verbose:
                print(f"[UNIT FIX] {da.name or var_key}: mm/s -> mm/day (×86400)")
            return _set_units(da * 86400.0, "converted from mm s-1 via ×86400")

        if u0 in ("m/day", "md-1", "m/d"):
            if verbose:
                print(f"[UNIT FIX] {da.name or var_key}: m/day -> mm/day (×1000)")
            return _set_units(da * 1000.0, "converted from m day-1 via ×1000")

        # heuristic fallback if units missing
        if not u0:
            try:
                vmax = float(da.max(skipna=True))
            except Exception:
                return da
            if 0 < vmax < 1e-2:
                if vmax < 5e-5:
                    if verbose:
                        print(f"[UNIT FIX] {da.name or var_key}: heuristic m/s -> mm/day (×86400×1000), vmax={vmax:.3g}")
                    return _set_units(da * 86400.0 * 1000.0, "heuristic: treated as m s-1 and converted to mm/day")
                else:
                    if verbose:
                        print(f"[UNIT FIX] {da.name or var_key}: heuristic kg m-2 s-1 -> mm/day (×86400), vmax={vmax:.3g}")
                    return _set_units(da * 86400.0, "heuristic: treated as kg m-2 s-1 and converted to mm/day")

        return da

    # ----------------------------
    # Helpers: landmask
    # ----------------------------

    @staticmethod
    def _pick_landmask_var(ds: xr.Dataset, prefer: Optional[str] = None) -> str:
        if prefer and prefer in ds.data_vars:
            return prefer
        candidates = [
            "landmask", "LANDMASK", "land_mask",
            "sftlf", "SFTLF",
            "landfrac", "LANDFRAC", "LANDFRAC_PFT",
            "mask", "MASK",
        ]
        for v in candidates:
            if v in ds.data_vars:
                return v
        # fallback: first data var
        if ds.data_vars:
            return list(ds.data_vars)[0]
        raise ValueError("Landmask dataset has no data variables.")

    def _get_land_mask_on_obs_grid(
        self,
        *,
        obs_ref: xr.DataArray,
        lat_obs: str,
        lon_obs: str,
    ) -> xr.DataArray:
        """
        Return boolean land mask on obs grid (lat_obs/lon_obs), cached.
        """
        path = self.cfg.landmask_file
        if not path:
            raise ValueError("Region requested mask='land' but cfg.landmask_file is not set.")

        # cache key: file + obs grid signature + threshold-var choice
        key = (
            os.path.abspath(path),
            self.cfg.landmask_var,
            float(self.cfg.land_threshold),
            lat_obs,
            lon_obs,
            int(obs_ref[lat_obs].size),
            int(obs_ref[lon_obs].size),
            float(np.asarray(obs_ref[lat_obs].values).min()),
            float(np.asarray(obs_ref[lat_obs].values).max()),
            float(np.asarray(obs_ref[lon_obs].values).min()),
            float(np.asarray(obs_ref[lon_obs].values).max()),
        )
        if key in self._landmask_cache:
            return self._landmask_cache[key]

        ds = xr.open_dataset(path)
        vname = self._pick_landmask_var(ds, prefer=self.cfg.landmask_var)
        lm = ds[vname]

        # squeeze extra dims if any
        for d in list(lm.dims):
            if d not in lm.coords:
                continue
        for d in list(lm.dims):
            if d not in ("lat", "latitude", "LAT", "nlat", "lon", "longitude", "LON", "nlon"):
                lm = lm.isel({d: 0})

        lat_lm, lon_lm = self._guess_lat_lon_names(lm)
        lm = self._rename_dim_if_present(lm, lat_lm, "lat")
        lm = self._rename_dim_if_present(lm, lon_lm, "lon")
        lm = self._wrap_lon_0_360(lm, "lon")

        # convert to fraction [0,1] if looks like percent
        lm_max = float(lm.max(skipna=True))
        lm_frac = lm / 100.0 if lm_max > 1.5 else lm

        # wrap obs lon to 0..360 for consistent interp
        obs_ref2 = self._wrap_lon_0_360(obs_ref, lon_obs)
        obs_lat = obs_ref2[lat_obs]
        obs_lon = obs_ref2[lon_obs]

        lm_on_obs = lm_frac.interp(lat=obs_lat, lon=obs_lon, method="nearest")

        # rename to obs coord names
        lm_on_obs = self._rename_dim_if_present(lm_on_obs, "lat", lat_obs)
        lm_on_obs = self._rename_dim_if_present(lm_on_obs, "lon", lon_obs)

        mask_bool = (lm_on_obs >= float(self.cfg.land_threshold))

        self._landmask_cache[key] = mask_bool
        return mask_bool

    @staticmethod
    def _rename_dim_if_present(da: xr.DataArray, old: str, new: str) -> xr.DataArray:
        if old == new:
            return da
        # if target name already exists, don't rename to avoid conflicts
        if (new in da.dims) or (new in da.coords):
            return da
        rename_map = {}
        if old in da.dims:
            rename_map[old] = new
        if old in da.coords:
            rename_map[old] = new
        if rename_map:
            return da.rename(rename_map)
        return da

    # ----------------------------
    # Main API
    # ----------------------------

    def compute(
        self,
        all_files: Dict[str, dict],
        *,
        period: str,
        obs_var_override: Optional[Dict[str, str]] = None,
    ) -> xr.Dataset:
        var_list = sorted(all_files.keys())
        if not var_list:
            raise ValueError("all_files is empty; nothing to compute.")

        first_var = var_list[0]
        exp_list = sorted(all_files[first_var].keys())
        if not exp_list:
            raise ValueError(f"No experiments found under var='{first_var}' in all_files.")

        # Regions: either user-specified or a single default region
        regions = self.cfg.regions
        if not regions:
            regions = [dict(name="DEFAULT", lat=self.cfg.lat_range, lon=self.cfg.lon_range, mask=None)]

        acc_out: List[xr.DataArray] = []
        rmse_out: List[xr.DataArray] = []
        acc_members_out: List[xr.DataArray] = []
        rmse_members_out: List[xr.DataArray] = []

        for var in var_list:
            if not all_files[var]:
                continue

            any_exp = next(iter(all_files[var].keys()))
            obs_files = all_files[var][any_exp]["obs"]
            obs_var = (obs_var_override or {}).get(var, var)

            # ---- open obs once per var ----
            obs_da0 = self._open_var(obs_files, obs_var)
            obs_da0 = self._maybe_convert_geopotential_to_height(obs_da0, var_key=var, verbose=self.cfg.verbose)
            obs_da0 = self._maybe_convert_precip_to_mmday(obs_da0, var_key=var, verbose=self.cfg.verbose)

            obs_da0 = self._select_period_time(obs_da0, period)
            obs_da0 = self._to_daily_date_time(obs_da0)

            lat_obs0, lon_obs0 = self._guess_lat_lon_names(obs_da0)
            obs_da0 = self._wrap_lon_0_360(obs_da0, lon_obs0)

            # after wrapping lon, re-detect names (still same, but safe)
            lat_obs0, lon_obs0 = self._guess_lat_lon_names(obs_da0)

            # ---- per-region loop ----
            acc_regions_var: List[xr.DataArray] = []
            rmse_regions_var: List[xr.DataArray] = []
            acc_members_regions_var: List[xr.DataArray] = []
            rmse_members_regions_var: List[xr.DataArray] = []

            for reg in regions:
                reg_name = str(reg.get("name", "REGION"))
                reg_lat = reg.get("lat", None)
                reg_lon = reg.get("lon", None)
                reg_mask = reg.get("mask", None)

                obs_reg = self._subset_region(obs_da0, lat_obs0, lon_obs0, reg_lat, reg_lon)
                if obs_reg.sizes.get("time", 0) == 0:
                    if self.cfg.verbose:
                        print(f"[SKIP] {var} region={reg_name}: obs has no data after subsetting.")
                    continue

                # weights on region obs grid
                wlat = self._coslat_weights(obs_reg[lat_obs0])
                w2d = wlat.broadcast_like(obs_reg.isel(time=0))

                # optional land mask (also apply to weights so denom is land-only)
                if reg_mask and str(reg_mask).lower() == "land":
                    land2d = self._get_land_mask_on_obs_grid(obs_ref=obs_reg.isel(time=0), lat_obs=lat_obs0, lon_obs=lon_obs0)
                    obs_reg = obs_reg.where(land2d)
                    w2d = w2d.where(land2d)

                # collect exp results for this region
                acc_exp_list: List[xr.DataArray] = []
                rmse_exp_list: List[xr.DataArray] = []
                acc_mem_exp_list: List[xr.DataArray] = []
                rmse_mem_exp_list: List[xr.DataArray] = []

                for exp in exp_list:
                    if exp not in all_files[var]:
                        continue

                    model_by_ens = all_files[var][exp]["model"]

                    # ---- ensemble mean ----
                    mod_mean = self._concat_ens_mean(model_by_ens, var=var, period=period)
                    mod_mean = self._maybe_convert_geopotential_to_height(mod_mean, var_key=var, verbose=self.cfg.verbose)
                    mod_mean = self._maybe_convert_precip_to_mmday(mod_mean, var_key=var, verbose=self.cfg.verbose)

                    mod_mean = self._select_period_time(mod_mean, period)
                    mod_mean = self._to_daily_date_time(mod_mean)

                    lat_mod, lon_mod = self._guess_lat_lon_names(mod_mean)
                    mod_mean = self._wrap_lon_0_360(mod_mean, lon_mod)
                    lat_mod, lon_mod = self._guess_lat_lon_names(mod_mean)

                    mod_reg = self._subset_region(mod_mean, lat_mod, lon_mod, reg_lat, reg_lon)

                    # rename model lat/lon to obs names if safe
                    mod_reg = self._ensure_latlon_names_match(
                        mod_reg, obs_reg,
                        mod_lat=lat_mod, mod_lon=lon_mod,
                        obs_lat=lat_obs0, obs_lon=lon_obs0,
                        tol=self.cfg.grid_tol,
                    )

                    # apply same land mask (on obs grid) if requested
                    if reg_mask and str(reg_mask).lower() == "land":
                        land2d = self._get_land_mask_on_obs_grid(obs_ref=obs_reg.isel(time=0), lat_obs=lat_obs0, lon_obs=lon_obs0)
                        mod_reg = mod_reg.where(land2d)

                    mod_aligned, obs_aligned = self._align_and_check(mod_reg, obs_reg, lat_name=lat_obs0, lon_name=lon_obs0)

                    if self.cfg.verbose:
                        print(
                            f"[DEBUG] {var} region={reg_name} exp={exp} aligned ranges:\n"
                            f"  model: min={float(mod_aligned.min()):.4g}, max={float(mod_aligned.max()):.4g}\n"
                            f"  obs:   min={float(obs_aligned.min()):.4g}, max={float(obs_aligned.max()):.4g}"
                        )

                    acc_em = self._pattern_corr(mod_aligned, obs_aligned, w2d)
                    rmse_em = self._weighted_rmse(mod_aligned, obs_aligned, w2d)

                    acc_exp_list.append(acc_em.assign_coords(exp=exp).expand_dims("exp"))
                    rmse_exp_list.append(rmse_em.assign_coords(exp=exp).expand_dims("exp"))

                    # ---- per-member ----
                    if self.cfg.compute_members:
                        acc_list: List[xr.DataArray] = []
                        rmse_list: List[xr.DataArray] = []
                        ens_names: List[str] = []

                        for ens, files in model_by_ens.items():
                            if not files:
                                continue

                            mem = self._open_member(files, var=var, period=period)
                            mem = self._maybe_convert_geopotential_to_height(mem, var_key=var, verbose=False)
                            mem = self._maybe_convert_precip_to_mmday(mem, var_key=var, verbose=False)

                            lat_m, lon_m = self._guess_lat_lon_names(mem)
                            mem = self._wrap_lon_0_360(mem, lon_m)
                            lat_m, lon_m = self._guess_lat_lon_names(mem)

                            mem_reg = self._subset_region(mem, lat_m, lon_m, reg_lat, reg_lon)
                            mem_reg = self._ensure_latlon_names_match(
                                mem_reg, obs_reg,
                                mod_lat=lat_m, mod_lon=lon_m,
                                obs_lat=lat_obs0, obs_lon=lon_obs0,
                                tol=self.cfg.grid_tol,
                            )

                            if reg_mask and str(reg_mask).lower() == "land":
                                land2d = self._get_land_mask_on_obs_grid(obs_ref=obs_reg.isel(time=0), lat_obs=lat_obs0, lon_obs=lon_obs0)
                                mem_reg = mem_reg.where(land2d)

                            mem_aligned, obs_mem = self._align_and_check(mem_reg, obs_reg, lat_name=lat_obs0, lon_name=lon_obs0)

                            acc_m = self._pattern_corr(mem_aligned, obs_mem, w2d)
                            rmse_m = self._weighted_rmse(mem_aligned, obs_mem, w2d)

                            acc_list.append(acc_m)
                            rmse_list.append(rmse_m)
                            ens_names.append(ens)

                        if acc_list:
                            acc_members = xr.concat(acc_list, dim=xr.IndexVariable("ens", ens_names))
                            rmse_members = xr.concat(rmse_list, dim=xr.IndexVariable("ens", ens_names))

                            acc_members = acc_members.assign_coords(exp=exp).expand_dims("exp")
                            rmse_members = rmse_members.assign_coords(exp=exp).expand_dims("exp")

                            acc_mem_exp_list.append(acc_members)
                            rmse_mem_exp_list.append(rmse_members)

                if not acc_exp_list:
                    continue

                acc_reg = xr.concat(acc_exp_list, dim="exp").assign_coords(region=reg_name).expand_dims("region")
                rmse_reg = xr.concat(rmse_exp_list, dim="exp").assign_coords(region=reg_name).expand_dims("region")

                acc_regions_var.append(acc_reg)
                rmse_regions_var.append(rmse_reg)

                if self.cfg.compute_members and acc_mem_exp_list:
                    accm_reg = xr.concat(acc_mem_exp_list, dim="exp").assign_coords(region=reg_name).expand_dims("region")
                    rmsem_reg = xr.concat(rmse_mem_exp_list, dim="exp").assign_coords(region=reg_name).expand_dims("region")
                    acc_members_regions_var.append(accm_reg)
                    rmse_members_regions_var.append(rmsem_reg)

            if not acc_regions_var:
                continue

            # (region, exp, time) -> add var
            ACC_var = xr.concat(acc_regions_var, dim="region").assign_coords(var=var).expand_dims("var")
            RMSE_var = xr.concat(rmse_regions_var, dim="region").assign_coords(var=var).expand_dims("var")
            acc_out.append(ACC_var)
            rmse_out.append(RMSE_var)

            if self.cfg.compute_members and acc_members_regions_var:
                ACCm_var = xr.concat(acc_members_regions_var, dim="region").assign_coords(var=var).expand_dims("var")
                RMSEm_var = xr.concat(rmse_members_regions_var, dim="region").assign_coords(var=var).expand_dims("var")
                acc_members_out.append(ACCm_var)
                rmse_members_out.append(RMSEm_var)

        if not acc_out:
            raise ValueError("No ACC/RMSE results produced (check inputs, regions, and period).")

        ACC = xr.concat(acc_out, dim="var").astype("float32")
        RMSE = xr.concat(rmse_out, dim="var").astype("float32")

        member_vars: Dict[str, xr.DataArray] = {}
        if self.cfg.compute_members and acc_members_out:
            ACCm = xr.concat(acc_members_out, dim="var").astype("float32")  # (var, region, exp, ens, time)
            RMSEm = xr.concat(rmse_members_out, dim="var").astype("float32")

            member_vars["ACC_member"] = ACCm
            member_vars["RMSE_member"] = RMSEm
            member_vars["ACC_member_mean"] = ACCm.mean("ens")
            member_vars["RMSE_member_mean"] = RMSEm.mean("ens")
            member_vars["ACC_member_std"] = ACCm.std("ens")
            member_vars["RMSE_member_std"] = RMSEm.std("ens")

            qs = self.cfg.member_quantiles
            if qs:
                member_vars["ACC_member_quantile"] = ACCm.quantile(list(qs), dim="ens").rename({"quantile": "q"})
                member_vars["RMSE_member_quantile"] = RMSEm.quantile(list(qs), dim="ens").rename({"quantile": "q"})

        # lead_days based on ACC time coordinate
        t = ACC["time"]
        lead = ((t - t.isel(time=0)) / np.timedelta64(1, "D")).astype("float32").rename("lead_days")

        ds = xr.Dataset(
            data_vars=dict(
                ACC=ACC,
                RMSE=RMSE,
                **member_vars,
            ),
            coords=dict(
                time=t,
                lead_days=lead,
                exp=ACC["exp"],
                var=ACC["var"],
                region=ACC["region"],
            ),
            attrs=dict(
                period=period,
                require_same_grid=str(self.cfg.require_same_grid),
                grid_tol=str(self.cfg.grid_tol),
                compute_members=str(getattr(self.cfg, "compute_members", False)),
                regions=",".join([str(r.get("name", "REGION")) for r in (self.cfg.regions or [])]) or "DEFAULT",
                landmask_file=str(self.cfg.landmask_file),
                note="ACC is weighted spatial pattern correlation; RMSE is weighted area RMSE; supports multi-region + optional land masking; unit fixes: geopotential->height and precip->mm/day.",
            ),
        )
        return ds

    def save(self, ds: xr.Dataset, out_nc: Optional[str] = None) -> str:
        path = out_nc or self.cfg.out_nc
        if (not self.cfg.overwrite) and os.path.exists(path):
            raise FileExistsError(f"Output exists: {path}")
        ds.to_netcdf(path)
        return path

class BiweeklyAccMapConfig:
    lat_range: Optional[Tuple[float, float]] = (-20.0, 20.0)
    lon_range: Optional[Tuple[float, float]] = None

    init_time: str = "2012-01-01"
    lead_bins: Tuple[Tuple[int, int], ...] = ((1, 14), (15, 28), (29, 42), (43, 56))
    lead_bin_labels: Tuple[str, ...] = ("wk1-2", "wk3-4", "wk5-6", "wk7-8")
    min_samples: int = 15

    require_same_grid: bool = True
    grid_tol: float = 0.0
    require_nonempty_time: bool = True

    # new
    compute_member_stats: bool = True
    save_member_maps: bool = False   # WARNING: can be large

    out_nc: str = "s2s_acc_map_biweekly.nc"
    overwrite: bool = True
    verbose : bool = True


class BiweeklyAccMapAssessor:
    """
    Biweekly ACC-map assessor.

    Output variables (always):
      - ACC_EMEAN(var, exp, lead_bin, lat, lon)
      - N_SAMPLES_EMEAN(var, exp, lead_bin, lat, lon)

    If cfg.compute_member_stats:
      - ACC_MEMBER_MEAN(var, exp, lead_bin, lat, lon)
      - ACC_MEMBER_STD(var, exp, lead_bin, lat, lon)
      - P_ACC_POS(var, exp, lead_bin, lat, lon)      # fraction of members with ACC>0
      - N_ENS_USED(var, exp)                         # count of members included (non-missing)
      - N_SAMPLES_MEMBER_MEAN(var, exp, lead_bin, lat, lon)  # mean paired samples across members

    If cfg.save_member_maps:
      - ACC_MEMBER(var, exp, ens, lead_bin, lat, lon)
      - N_SAMPLES_MEMBER(var, exp, ens, lead_bin, lat, lon)
    """

    def __init__(self, cfg: BiweeklyAccMapConfig):
        self.cfg = cfg
        self._validate_cfg()

    # ----------------------------
    # Helpers (mostly unchanged)
    # ----------------------------
    def _validate_cfg(self) -> None:
        if len(self.cfg.lead_bins) != len(self.cfg.lead_bin_labels):
            raise ValueError("lead_bins and lead_bin_labels must have same length.")
        if self.cfg.min_samples < 2:
            raise ValueError("min_samples must be >= 2.")
        for a, b in self.cfg.lead_bins:
            if a > b:
                raise ValueError(f"Invalid lead bin ({a},{b}).")
        _ = np.datetime64(self.cfg.init_time)

    @staticmethod
    def _parse_yyyymm_period(period: str) -> Tuple[int, int, int, int]:
        m = re.match(r"^(\d{4})(\d{2})-(\d{4})(\d{2})$", period)
        if not m:
            raise ValueError(f"Bad period '{period}', expected 'YYYYMM-YYYYMM'")
        y0, m0, y1, m1 = map(int, m.groups())
        if (y1, m1) < (y0, m0):
            raise ValueError(f"Bad period '{period}' (end < start)")
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
    def _open_var(files: List[str], var: str) -> xr.DataArray:
        ds = xr.open_mfdataset(files, combine="by_coords")
        if var in ds:
            return ds[var]
        for k in ds.data_vars:
            if k.lower() == var.lower():
                return ds[k]
        raise KeyError(f"Variable '{var}' not found. Available={list(ds.data_vars)}")

    @staticmethod
    def _guess_lat_lon_names(da: xr.DataArray) -> Tuple[str, str]:
        lat_candidates = ["lat", "latitude", "LAT", "nlat"]
        lon_candidates = ["lon", "longitude", "LON", "nlon"]
        lat = next((n for n in lat_candidates if (n in da.dims or n in da.coords)), None)
        lon = next((n for n in lon_candidates if (n in da.dims or n in da.coords)), None)
        if lat is None or lon is None:
            raise ValueError(f"Cannot infer lat/lon from dims={da.dims}, coords={list(da.coords)}")
        return lat, lon

    @staticmethod
    def _subset_region(
        da: xr.DataArray,
        lat_name: str,
        lon_name: str,
        lat_range: Optional[Tuple[float, float]],
        lon_range: Optional[Tuple[float, float]],
    ) -> xr.DataArray:
        out = da
        if lat_range is not None:
            lo, hi = lat_range
            latv = out[lat_name]
            if latv.size > 1 and (latv[0] > latv[-1]):
                out = out.sel({lat_name: slice(hi, lo)})
            else:
                out = out.sel({lat_name: slice(lo, hi)})
        if lon_range is not None:
            lo, hi = lon_range
            out = out.sel({lon_name: slice(lo, hi)})
        return out

    @staticmethod
    def _same_1d_values(a: xr.DataArray, b: xr.DataArray, *, tol: float = 0.0) -> bool:
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
        rename_map = {}
        if old in da.dims:
            rename_map[old] = new
        if old in da.coords:
            rename_map[old] = new
        if rename_map:
            return da.rename(rename_map)
        return da

    @classmethod
    def _ensure_latlon_names_match(
        cls,
        mod: xr.DataArray,
        obs: xr.DataArray,
        *,
        mod_lat: str,
        mod_lon: str,
        obs_lat: str,
        obs_lon: str,
        tol: float = 0.0,
    ) -> xr.DataArray:
        if (mod_lat == obs_lat) and (mod_lon == obs_lon):
            return mod

        if mod_lat not in mod.coords or mod_lon not in mod.coords:
            raise ValueError(f"Model missing coords '{mod_lat}'/'{mod_lon}'. coords={list(mod.coords)}")
        if obs_lat not in obs.coords or obs_lon not in obs.coords:
            raise ValueError(f"Obs missing coords '{obs_lat}'/'{obs_lon}'. coords={list(obs.coords)}")

        lat_ok = cls._same_1d_values(mod[mod_lat], obs[obs_lat], tol=tol)
        lon_ok = cls._same_1d_values(mod[mod_lon], obs[obs_lon], tol=tol)

        if not (lat_ok and lon_ok):
            raise ValueError("Model/obs lat/lon values differ; regrid explicitly.")

        mod = cls._rename_dim_if_present(mod, mod_lat, obs_lat)
        mod = cls._rename_dim_if_present(mod, mod_lon, obs_lon)
        return mod

    @staticmethod
    def _to_daily_date_time(obj: xr.DataArray | xr.Dataset) -> xr.DataArray | xr.Dataset:
        if "time" not in obj.coords:
            raise KeyError("Object has no 'time' coordinate.")
        t = obj["time"]
        if np.issubdtype(t.dtype, np.datetime64):
            t_daily = t.astype("datetime64[D]")
        else:
            t_daily = xr.DataArray(
                np.array([np.datetime64(str(v)[:10]) for v in t.values], dtype="datetime64[D]"),
                dims=t.dims,
                coords=t.coords,
                name=t.name,
            )
        return obj.assign_coords(time=t_daily)

    def _align_and_check(self, mod: xr.DataArray, obs: xr.DataArray, *, lat_name: str, lon_name: str):
        mod, obs = xr.align(mod, obs, join="inner")
        mod = mod.sortby("time")
        obs = obs.sortby("time")
        if self.cfg.require_nonempty_time and mod.sizes.get("time", 0) == 0:
            raise ValueError("No overlapping times between model and obs.")
        if self.cfg.require_same_grid:
            if not self._same_1d_values(mod[lat_name], obs[lat_name], tol=self.cfg.grid_tol):
                raise ValueError("lat grids differ; regrid or disable require_same_grid.")
            if not self._same_1d_values(mod[lon_name], obs[lon_name], tol=self.cfg.grid_tol):
                raise ValueError("lon grids differ; regrid or disable require_same_grid.")
        return mod, obs

    @classmethod
    def _concat_ens_mean(cls, model_files_by_ens: Dict[str, List[str]], *, var: str, period: str) -> xr.DataArray:
        members = []
        ens_names = []
        for ens, files in model_files_by_ens.items():
            if not files:
                continue
            da = cls._open_var(files, var)
            da = cls._select_period_time(da, period)
            members.append(da)
            ens_names.append(ens)
        if not members:
            raise ValueError(f"No model files found for var={var} across ensembles.")
        ens_da = xr.concat(members, dim=xr.IndexVariable("ens", ens_names))
        return ens_da.mean("ens")

    # ----------------------------
    # NEW: ACC-map core
    # ----------------------------
    def _add_lead_day_coord(self, da: xr.DataArray) -> xr.DataArray:
        init = np.datetime64(self.cfg.init_time)
        t = da["time"].values
        lead = ((t - init) / np.timedelta64(1, "D")).astype("float32")
        return da.assign_coords(lead_day=("time", lead))

    @staticmethod
    def _maybe_convert_geopotential_to_height(
        da: xr.DataArray,
        *,
        var_key: str,
        prefer: str = "m",   # normalize to meters
        g: float = 9.80665,
        verbose: bool = False,
    ) -> xr.DataArray:
        """
        Ensure geopotential-type variables are in geopotential height (m).

        Handles:
          - ERA-style geopotential: units like "m2 s-2", "m^2/s^2", "m**2 s**-2"
          - Height already: units like "m", "meter"
          - Heuristic fallback based on magnitude if units missing/garbled

        var_key: your loop variable name (e.g., "Z500", "Z850", "Z", etc.)
        """

        # Only attempt for Z-like variables.
        # Covers: Z, Zxxx, Z500, Z850, zg, geopotential, etc.
        k = (var_key or "").lower()
        name = (da.name or "").lower()
        if not (
            k.startswith("z") or
            "geopot" in k or "geopot" in name or
            name.startswith("z") or name in ("zg", "phi")
        ):
            return da

        units = (da.attrs.get("units") or "").strip().lower().replace("^", "").replace("**", "")
        units = units.replace(" ", "")

        def _is_geopotential_units(u: str) -> bool:
            # common encodings of m^2 s^-2
            return (
                "m2s-2" in u or
                "m2/s2" in u or
                "m2s**-2" in u or
                "m2s^-2" in u
            )

        def _is_height_units(u: str) -> bool:
            return u in ("m", "meter", "metre", "meters", "metres")

        # If units tell us clearly
        if units:
            if _is_geopotential_units(units):
                if verbose:
                    print(f"[UNIT FIX] {da.name or var_key}: geopotential -> height (divide by g)")
                out = da / g
                out.attrs = dict(da.attrs)
                out.attrs["units"] = "m"
                out.attrs["note_unit_fix"] = "converted from geopotential (m2 s-2) via /g"
                return out

            if _is_height_units(units):
                return da

        # Fallback heuristic (units missing or untrusted):
        # geopotential ~ O(1e4-1e5); height ~ O(1e3-1e4)
        try:
            vmax = float(da.max(skipna=True))
        except Exception:
            return da

        # These thresholds are conservative and work across pressure levels.
        # Heights at low levels can be small, but won't exceed ~2e4 m.
        if vmax > 2.0e4 and vmax < 5.0e6:
            if verbose:
                print(f"[UNIT FIX] {da.name or var_key}: heuristic geopotential -> height (vmax={vmax:.3g})")
            out = da / g
            out.attrs = dict(da.attrs)
            out.attrs["units"] = "m"
            out.attrs["note_unit_fix"] = "heuristic: treated as geopotential and converted via /g"
            return out

        return da

    @staticmethod
    def _maybe_convert_precip_to_mmday(
        da: xr.DataArray,
        *,
        var_key: str,
        verbose: bool = False,
    ) -> xr.DataArray:
        """
        Normalize precipitation-like variables to mm/day.

        Handles common units:
          - m/s
          - kg m-2 s-1  (or kg/m^2/s)
          - mm/day, mm d-1
          - mm/s, m/day, etc.
        Falls back to a magnitude heuristic if units missing.
        """
        k = (var_key or "").lower()
        name = (da.name or "").lower()

        # decide if this looks like precipitation
        is_pr = (
            k in ("pr", "precip", "prect", "precc", "precl", "tp") or
            "precip" in k or "precip" in name or
            k.startswith("pr") or name.startswith("pr") or
            name in ("prect", "precc", "precl", "tp")
        )
        if not is_pr:
            return da

        u = (da.attrs.get("units") or "").strip().lower()
        u0 = u.replace(" ", "")

        def _set_units(out: xr.DataArray, note: str) -> xr.DataArray:
            out.attrs = dict(da.attrs)
            out.attrs["units"] = "mm/day"
            out.attrs["note_unit_fix"] = note
            return out

        # already mm/day-ish
        if u0 in ("mm/day", "mmd-1", "mm/d", "mmdy-1") or ("mm" in u0 and ("day" in u0 or "d-1" in u0)):
            return da

        # kg m-2 s-1  -> mm/day
        if ("kg" in u0 and "m-2" in u0 and "s-1" in u0) or u0 in ("kgm-2s-1", "kg/m2/s", "kgm**-2s**-1"):
            if verbose:
                print(f"[UNIT FIX] {da.name or var_key}: kg m-2 s-1 -> mm/day (×86400)")
            out = da * 86400.0
            return _set_units(out, "converted from kg m-2 s-1 via ×86400")

        # m/s -> mm/day
        if u0 in ("m/s", "ms-1", "msec-1") or ("m" in u0 and "s-1" in u0 and "kg" not in u0 and "mm" not in u0):
            if verbose:
                print(f"[UNIT FIX] {da.name or var_key}: m/s -> mm/day (×86400×1000)")
            out = da * 86400.0 * 1000.0
            return _set_units(out, "converted from m s-1 via ×86400×1000")

        # mm/s -> mm/day
        if u0 in ("mm/s", "mms-1"):
            if verbose:
                print(f"[UNIT FIX] {da.name or var_key}: mm/s -> mm/day (×86400)")
            out = da * 86400.0
            return _set_units(out, "converted from mm s-1 via ×86400")

        # m/day -> mm/day
        if u0 in ("m/day", "md-1", "m/d"):
            if verbose:
                print(f"[UNIT FIX] {da.name or var_key}: m/day -> mm/day (×1000)")
            out = da * 1000.0
            return _set_units(out, "converted from m day-1 via ×1000")

        # Fallback heuristic (if units missing): typical precip rates are small in SI.
        # If values are ~1e-7 to 1e-4 it might be m/s; if ~1e-6 to 1e-3 might be kg/m2/s.
        # We DO NOT auto-convert without some confidence; use a conservative rule.
        if not u0:
            try:
                vmax = float(da.max(skipna=True))
            except Exception:
                return da

            # If vmax is tiny (<1e-2) it's likely a rate in per-second units (m/s or kg/m2/s).
            # Distinguish by checking if it's closer to m/s scale (very tiny) vs kg/m2/s (slightly larger),
            # but both can overlap. We'll choose kg/m2/s unless extremely tiny.
            if 0 < vmax < 1e-2:
                # extremely tiny -> treat as m/s
                if vmax < 5e-5:
                    if verbose:
                        print(f"[UNIT FIX] {da.name or var_key}: heuristic m/s -> mm/day (vmax={vmax:.3g})")
                    out = da * 86400.0 * 1000.0
                    return _set_units(out, "heuristic: treated as m s-1 and converted to mm/day")
                else:
                    if verbose:
                        print(f"[UNIT FIX] {da.name or var_key}: heuristic kg m-2 s-1 -> mm/day (vmax={vmax:.3g})")
                    out = da * 86400.0
                    return _set_units(out, "heuristic: treated as kg m-2 s-1 and converted to mm/day")

        return da

    @staticmethod
    def _nan_corr_1d(x: np.ndarray, y: np.ndarray) -> float:
        m = np.isfinite(x) & np.isfinite(y)
        if m.sum() < 2:
            return np.nan
        x0 = x[m] - np.mean(x[m])
        y0 = y[m] - np.mean(y[m])
        den = np.sqrt(np.sum(x0 * x0) * np.sum(y0 * y0))
        if den == 0:
            return np.nan
        return float(np.sum(x0 * y0) / den)

    def _select_lead_bin(self, da: xr.DataArray, a: int, b: int) -> xr.DataArray:
        m = (da["lead_day"] >= float(a)) & (da["lead_day"] <= float(b))
        return da.sel(time=da.time.where(m, drop=True))

    def _corr_map_time(self, fc: xr.DataArray, obs: xr.DataArray) -> tuple[xr.DataArray, xr.DataArray]:
        """
        Correlation across time at each (lat,lon). Returns (acc_map, n_samples_map).
        Handles dask chunking on time core dim.
        """
        tname = "time"

        # Fix for your error: core dim time must be single chunk for dask='parallelized'
        if hasattr(fc.data, "chunks"):
            fc = fc.chunk({tname: -1})
        if hasattr(obs.data, "chunks"):
            obs = obs.chunk({tname: -1})

        paired_valid = xr.apply_ufunc(
            lambda a, b: np.isfinite(a) & np.isfinite(b),
            fc, obs,
            input_core_dims=[[tname], [tname]],
            output_core_dims=[[tname]],
            vectorize=True,
            dask="parallelized",
            output_dtypes=[bool],
        )
        n_samples = paired_valid.sum(tname).astype("int32")

        acc = xr.apply_ufunc(
            self._nan_corr_1d,
            fc, obs,
            input_core_dims=[[tname], [tname]],
            output_core_dims=[[]],
            vectorize=True,
            dask="parallelized",
            output_dtypes=[float],
        ).astype("float32")

        acc = acc.where(n_samples >= self.cfg.min_samples)
        return acc, n_samples

    def _compute_biweekly_maps(self, fc: xr.DataArray, obs: xr.DataArray) -> tuple[xr.DataArray, xr.DataArray]:
        acc_list = []
        ns_list = []
        for (a, b), lab in zip(self.cfg.lead_bins, self.cfg.lead_bin_labels):
            fc_bin = self._select_lead_bin(fc, a, b)
            obs_bin = self._select_lead_bin(obs, a, b)

            if fc_bin.sizes.get("time", 0) == 0 or obs_bin.sizes.get("time", 0) == 0:
                template = fc.isel(time=0, drop=True)
                acc = xr.full_like(template, np.nan).astype("float32")
                ns = xr.zeros_like(template, dtype="int32")
            else:
                acc, ns = self._corr_map_time(fc_bin, obs_bin)

            acc = acc.expand_dims(lead_bin=[lab])
            ns = ns.expand_dims(lead_bin=[lab])
            acc_list.append(acc)
            ns_list.append(ns)

        acc_all = xr.concat(acc_list, dim="lead_bin")
        ns_all = xr.concat(ns_list, dim="lead_bin")
        return acc_all, ns_all

    # ----------------------------
    # NEW: member stats
    # ----------------------------
    def _compute_member_products(
        self,
        model_by_ens: Dict[str, List[str]],
        *,
        var: str,
        period: str,
        obs_da: xr.DataArray,
        lat_obs: str,
        lon_obs: str,
    ):
        """
        Returns a dict of member summary fields (+ optional full member maps).
        All outputs are lead_bin x lat x lon (except N_ENS_USED scalar).
        """
        acc_members = []
        ns_members = []
        ens_names = []

        for ens, files in model_by_ens.items():
            if not files:
                continue
            try:
                mem = self._open_var(files, var)
                mem = self._maybe_convert_geopotential_to_height(
                    mem, var_key=var, verbose=self.cfg.verbose
                )
                mem = self._maybe_convert_precip_to_mmday(
                    mem, var_key=var, verbose=self.cfg.verbose
                )
            except Exception:
                continue
            mem = self._select_period_time(mem, period)
            mem = self._to_daily_date_time(mem)

            lat_mod, lon_mod = self._guess_lat_lon_names(mem)
            mem = self._subset_region(mem, lat_mod, lon_mod, self.cfg.lat_range, self.cfg.lon_range)
            mem = self._ensure_latlon_names_match(
                mem, obs_da,
                mod_lat=lat_mod, mod_lon=lon_mod,
                obs_lat=lat_obs, obs_lon=lon_obs,
                tol=self.cfg.grid_tol,
            )

            mem_al, obs_al = self._align_and_check(mem, obs_da, lat_name=lat_obs, lon_name=lon_obs)

            mem_al = self._add_lead_day_coord(mem_al)
            obs_al = self._add_lead_day_coord(obs_al)

            acc_map, ns_map = self._compute_biweekly_maps(mem_al, obs_al)

            acc_members.append(acc_map)
            ns_members.append(ns_map)
            ens_names.append(ens)

        n_ens_used = len(acc_members)
        if n_ens_used == 0:
            # Create NaN placeholders based on obs grid
            template = obs_da.isel(time=0, drop=True)
            acc_nan = xr.full_like(template, np.nan).astype("float32")
            ns_zero = xr.zeros_like(template, dtype="int32")
            acc_nan = xr.concat([acc_nan.expand_dims(lead_bin=[lab]) for lab in self.cfg.lead_bin_labels], dim="lead_bin")
            ns_zero = xr.concat([ns_zero.expand_dims(lead_bin=[lab]) for lab in self.cfg.lead_bin_labels], dim="lead_bin")
            return {
                "ACC_MEMBER_MEAN": acc_nan,
                "ACC_MEMBER_STD": acc_nan,
                "P_ACC_POS": acc_nan,
                "N_SAMPLES_MEMBER_MEAN": ns_zero.astype("float32"),
                "N_ENS_USED": xr.DataArray(0, dims=()),
                "ACC_MEMBER": None,
                "N_SAMPLES_MEMBER": None,
            }

        acc_ens = xr.concat(acc_members, dim=xr.IndexVariable("ens", ens_names))
        ns_ens = xr.concat(ns_members, dim=xr.IndexVariable("ens", ens_names))

        acc_member_mean = acc_ens.mean("ens", skipna=True).astype("float32")
        acc_member_std = acc_ens.std("ens", skipna=True).astype("float32")
        p_acc_pos = (acc_ens > 0).mean("ens", skipna=True).astype("float32")
        ns_member_mean = ns_ens.mean("ens", skipna=True).astype("float32")

        out = {
            "ACC_MEMBER_MEAN": acc_member_mean,
            "ACC_MEMBER_STD": acc_member_std,
            "P_ACC_POS": p_acc_pos,
            "N_SAMPLES_MEMBER_MEAN": ns_member_mean,
            "N_ENS_USED": xr.DataArray(n_ens_used, dims=()),
            "ACC_MEMBER": acc_ens if self.cfg.save_member_maps else None,
            "N_SAMPLES_MEMBER": ns_ens if self.cfg.save_member_maps else None,
        }
        return out

    # ----------------------------
    # Main API (same signature style)
    # ----------------------------
    def compute(
        self,
        all_files: Dict[str, dict],
        *,
        period: str,
        obs_var_override: Optional[Dict[str, str]] = None,
    ) -> xr.Dataset:

        var_list = sorted(all_files.keys())
        if not var_list:
            raise ValueError("all_files is empty; nothing to compute.")

        first_var = var_list[0]
        exp_list = sorted(all_files[first_var].keys())
        if not exp_list:
            raise ValueError(f"No experiments found under var='{first_var}'.")

        acc_emean_out = []
        ns_emean_out = []

        # member outputs (optional)
        acc_mem_mean_out = []
        acc_mem_std_out = []
        p_pos_out = []
        n_ens_used_out = []
        ns_mem_mean_out = []
        acc_member_full_out = []
        ns_member_full_out = []

        for var in var_list:
            if not all_files[var]:
                continue

            any_exp = next(iter(all_files[var].keys()))
            obs_files = all_files[var][any_exp]["obs"]
            obs_var = (obs_var_override or {}).get(var, var)

            obs_da = self._open_var(obs_files, obs_var)
            # ---- unit fix for Z500: geopotential (m2/s2) -> geopotential height (m)
            obs_da = self._maybe_convert_geopotential_to_height(
                obs_da, var_key=var, verbose=self.cfg.verbose
            )
            obs_da = self._maybe_convert_precip_to_mmday(
                obs_da, var_key=var, verbose=self.cfg.verbose
            )

            obs_da = self._select_period_time(obs_da, period)
            obs_da = self._to_daily_date_time(obs_da)

            lat_obs, lon_obs = self._guess_lat_lon_names(obs_da)
            obs_da = self._subset_region(obs_da, lat_obs, lon_obs, self.cfg.lat_range, self.cfg.lon_range)

            if obs_da.sizes.get("time", 0) == 0:
                raise ValueError(f"OBS has no data after period selection for var={var} (obs_var={obs_var}).")

            acc_var = []
            ns_var = []

            # member stat containers per var
            acc_mem_mean_var = []
            acc_mem_std_var = []
            p_pos_var = []
            n_ens_used_var = []
            ns_mem_mean_var = []
            acc_member_full_var = []
            ns_member_full_var = []

            for exp in exp_list:
                if exp not in all_files[var]:
                    continue

                model_by_ens = all_files[var][exp]["model"]

                # ---- ensemble mean product (main)
                mod_mean = self._concat_ens_mean(model_by_ens, var=var, period=period)
                mod_mean = self._to_daily_date_time(mod_mean)
                mod_mean = self._maybe_convert_geopotential_to_height(
                        mod_mean, var_key=var, verbose=self.cfg.verbose
                )
                mod_mean = self._maybe_convert_precip_to_mmday(
                        mod_mean, var_key=var, verbose=self.cfg.verbose
                )

                lat_mod, lon_mod = self._guess_lat_lon_names(mod_mean)
                mod_mean = self._subset_region(mod_mean, lat_mod, lon_mod, self.cfg.lat_range, self.cfg.lon_range)
                mod_mean = self._ensure_latlon_names_match(
                    mod_mean, obs_da,
                    mod_lat=lat_mod, mod_lon=lon_mod,
                    obs_lat=lat_obs, obs_lon=lon_obs,
                    tol=self.cfg.grid_tol,
                )

                mod_al, obs_al = self._align_and_check(mod_mean, obs_da, lat_name=lat_obs, lon_name=lon_obs)
                mod_al = self._add_lead_day_coord(mod_al)
                obs_al = self._add_lead_day_coord(obs_al)

                acc_maps, ns_maps = self._compute_biweekly_maps(mod_al, obs_al)

                acc_maps = acc_maps.assign_coords(exp=exp).expand_dims("exp")
                ns_maps = ns_maps.assign_coords(exp=exp).expand_dims("exp")
                acc_var.append(acc_maps)
                ns_var.append(ns_maps)

                # ---- member stats (optional)
                if self.cfg.compute_member_stats:
                    mem_out = self._compute_member_products(
                        model_by_ens,
                        var=var,
                        period=period,
                        obs_da=obs_da,
                        lat_obs=lat_obs,
                        lon_obs=lon_obs,
                    )

                    acc_mem_mean_var.append(mem_out["ACC_MEMBER_MEAN"].assign_coords(exp=exp).expand_dims("exp"))
                    acc_mem_std_var.append(mem_out["ACC_MEMBER_STD"].assign_coords(exp=exp).expand_dims("exp"))
                    p_pos_var.append(mem_out["P_ACC_POS"].assign_coords(exp=exp).expand_dims("exp"))
                    ns_mem_mean_var.append(mem_out["N_SAMPLES_MEMBER_MEAN"].assign_coords(exp=exp).expand_dims("exp"))
                    n_ens_used_var.append(mem_out["N_ENS_USED"].assign_coords(exp=exp).expand_dims("exp"))

                    if self.cfg.save_member_maps and (mem_out["ACC_MEMBER"] is not None):
                        acc_member_full_var.append(mem_out["ACC_MEMBER"].assign_coords(exp=exp).expand_dims("exp"))
                        ns_member_full_var.append(mem_out["N_SAMPLES_MEMBER"].assign_coords(exp=exp).expand_dims("exp"))

            if not acc_var:
                continue

            acc_emean = xr.concat(acc_var, dim="exp").assign_coords(var=var).expand_dims("var")
            ns_emean = xr.concat(ns_var, dim="exp").assign_coords(var=var).expand_dims("var")
            acc_emean_out.append(acc_emean)
            ns_emean_out.append(ns_emean)

            if self.cfg.compute_member_stats:
                acc_mem_mean_out.append(xr.concat(acc_mem_mean_var, dim="exp").assign_coords(var=var).expand_dims("var"))
                acc_mem_std_out.append(xr.concat(acc_mem_std_var, dim="exp").assign_coords(var=var).expand_dims("var"))
                p_pos_out.append(xr.concat(p_pos_var, dim="exp").assign_coords(var=var).expand_dims("var"))
                ns_mem_mean_out.append(xr.concat(ns_mem_mean_var, dim="exp").assign_coords(var=var).expand_dims("var"))
                n_ens_used_out.append(xr.concat(n_ens_used_var, dim="exp").assign_coords(var=var).expand_dims("var"))

                if self.cfg.save_member_maps and acc_member_full_var:
                    acc_member_full_out.append(xr.concat(acc_member_full_var, dim="exp").assign_coords(var=var).expand_dims("var"))
                    ns_member_full_out.append(xr.concat(ns_member_full_var, dim="exp").assign_coords(var=var).expand_dims("var"))

        if not acc_emean_out:
            raise ValueError("No ACC map results produced (check inputs and period).")

        ACC_EMEAN = xr.concat(acc_emean_out, dim="var").astype("float32")
        N_SAMPLES_EMEAN = xr.concat(ns_emean_out, dim="var").astype("int32")

        data_vars = {
            "ACC_EMEAN": ACC_EMEAN,
            "N_SAMPLES_EMEAN": N_SAMPLES_EMEAN,
        }

        if self.cfg.compute_member_stats:
            data_vars.update({
                "ACC_MEMBER_MEAN": xr.concat(acc_mem_mean_out, dim="var").astype("float32"),
                "ACC_MEMBER_STD": xr.concat(acc_mem_std_out, dim="var").astype("float32"),
                "P_ACC_POS": xr.concat(p_pos_out, dim="var").astype("float32"),
                "N_SAMPLES_MEMBER_MEAN": xr.concat(ns_mem_mean_out, dim="var").astype("float32"),
                "N_ENS_USED": xr.concat(n_ens_used_out, dim="var").astype("int32"),
            })

        if self.cfg.save_member_maps and acc_member_full_out:
            data_vars["ACC_MEMBER"] = xr.concat(acc_member_full_out, dim="var").astype("float32")
            data_vars["N_SAMPLES_MEMBER"] = xr.concat(ns_member_full_out, dim="var").astype("int32")

        ds = xr.Dataset(
            data_vars=data_vars,
            coords=dict(
                var=ACC_EMEAN["var"],
                exp=ACC_EMEAN["exp"],
                lead_bin=ACC_EMEAN["lead_bin"],
            ),
            attrs=dict(
                period=period,
                init_time=self.cfg.init_time,
                lead_bins=",".join([f"{a}-{b}" for a, b in self.cfg.lead_bins]),
                lead_bin_labels=",".join(self.cfg.lead_bin_labels),
                lat_range=str(self.cfg.lat_range),
                lon_range=str(self.cfg.lon_range),
                require_same_grid=str(self.cfg.require_same_grid),
                grid_tol=str(self.cfg.grid_tol),
                min_samples=str(self.cfg.min_samples),
                compute_member_stats=str(self.cfg.compute_member_stats),
                save_member_maps=str(self.cfg.save_member_maps),
                note="ACC_* are correlation across TIME at each grid point within each lead-bin. "
                     "EMEAN uses ensemble-mean first. MEMBER_* summarizes per-member ACC maps.",
            ),
        )
        return ds

    def save(self, ds: xr.Dataset, out_nc: Optional[str] = None) -> str:
        path = out_nc or self.cfg.out_nc
        if (not self.cfg.overwrite) and os.path.exists(path):
            raise FileExistsError(f"Output exists: {path}")
        ds.to_netcdf(path)
        return path

    def save_per_var_exp(
        self,
        ds: xr.Dataset,
        *,
        out_dir: str,
        group: str,
        freq: str,
        run: str,
        obs: str,
        overwrite: Optional[bool] = None,
    ) -> List[str]:
        """
        Split ds into per (var, exp) netcdf files:
          s2s_accmap_{group}_{freq}_{run}_{obs}_{exp}_{var}.nc
        """
        out = []
        ow = self.cfg.overwrite if overwrite is None else overwrite
        os.makedirs(out_dir, exist_ok=True)

        for var in ds["var"].values.tolist():
            for exp in ds["exp"].values.tolist():
                ds_one = ds.sel(var=var, exp=exp)
                fn = f"s2s_accmap_{group}_{freq}_{run}_{obs}_{exp}_{var}.nc"
                path = os.path.join(out_dir, fn)
                if (not ow) and os.path.exists(path):
                    out.append(path)
                    continue
                ds_one.to_netcdf(path)
                out.append(path)
        return out
