"""Surface-flux time-series extraction from MPAS-O (accumulated fluxes -> monthly/annual series).

v1 classes: surface ocean and TKE/AMOC analyses. *V2 names at the end: the ocean-heat-uptake version.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple, Sequence, Any, Union, List

import os
import re
import glob
import numpy as np
import xarray as xr
import cftime


YearRange = Tuple[int, int]


# ==============================
# Requests / Options
# ==============================
@dataclass
class FluxExtractOptions:
    engine: Optional[str] = "netcdf4"
    chunks: Optional[dict] = None

    # time handling
    calendar: str = "noleap"
    # if your files have time decoding issues, we "repair" time attrs + reassign
    reassign_monthly_time: bool = True
    length_policy: str = "auto"  # "auto" | "error" | "trim_ds" | "truncate_new_time"
    verbose: bool = False


@dataclass
class FluxSeriesRequest:
    """
    A single extracted series.

    key:
      - short key in var_dict (e.g., "Qnet", "Mass_net")
      - OR a raw stem (e.g., "netEnergyFlux")
      - OR a derived key (defined in derived_recipes)

    scale:
      multiply final series (e.g., convert per-second to per-year)
    per_area:
      if True, divide by global area A (must be provided via global_area_m2 or mesh_path or exemplar)
    """
    key: str
    scale: float = 1.0
    per_area: bool = False


# ==============================
# Batch extraction config (matches your __main__ usage)
# ==============================
@dataclass
class FluxExtractConfig:
    """
    Batch-extraction config for segmented flux files.

    ts_window: (year_min, year_max) in simulation years (integer years)
    annual_mean: if True return annual series; else monthly
    """
    ts_window: YearRange
    annual_mean: bool = False
    annual_weighted: bool = True
    annual_drop_incomplete: bool = True
    annual_min_days: int = 360

    calendar: str = "noleap"
    reassign_monthly_time: bool = True
    length_policy: str = "auto"


# ==============================
# Extractor
# ==============================
class FluxTimeSeriesExtractor:
    """
    Flux-style MPAS-O extractor for segmented files:
      <stem>_YYYYMM-YYYYMM.nc
    """

    _date_re = re.compile(r"_(\d{6})-(\d{6})\.nc$")
    _ym_re   = re.compile(r"_(\d{6})-(\d{6})")

    def __init__(
        self,
        *,
        options: Optional[FluxExtractOptions] = None,
        var_dict: Optional[Dict[str, Dict[str, Any]]] = None,
        derived_recipes: Optional[Dict[str, Dict[str, Any]]] = None,
        # area normalization (optional)
        global_area_m2: Optional[float] = None,
        mesh_path: Optional[str] = None,
        exemplar: Optional[Tuple[str, str, str]] = None,  # (base_dir, exp, sub_dir) containing mpasTimeSeriesOcean.nc
    ):
        self.opt = options or FluxExtractOptions()

        # var_dict schema: {shortkey: {"name": stem, "index": i, "glob": f"{stem}_*.nc"}}
        self._var_dict = var_dict or self.build_default_var_dict()
        self._stem_to_short = {v["name"]: k for k, v in self._var_dict.items()}

        # derived recipes schema:
        #   derived[key] = {"terms": [(subkey, coef), ...]} OR {"expr": "Qnet - ..."}
        self._derived = derived_recipes or {}

        self._A = None
        if (global_area_m2 is not None) or mesh_path or exemplar:
            self._A = self._init_global_area_once(global_area_m2, mesh_path, exemplar)

        self._series_cache: Dict[Tuple, xr.DataArray] = {}

    # --------------------------
    # var_dict
    # --------------------------
    @staticmethod
    def build_default_var_dict() -> Dict[str, Dict[str, Any]]:
        shortname_map = {
            # Energy
            "netEnergyFlux": "Qnet_mpas",   # 0 in the first month after a restart; Qnet is derived (alt_global_timeseries)
            "accumulatedShortWaveHeatFlux": "Qsw",
            "accumulatedLongWaveHeatFluxDown": "Qlw_dn",
            "accumulatedLongWaveHeatFluxUp": "Qlw_up",
            "accumulatedLatentHeatFlux": "Qlh",
            "accumulatedSensibleHeatFlux": "Qsh",
            "accumulatedSeaIceHeatFlux": "Qsi",
            "accumulatedFrazilHeatFlux": "Qfrazil",
            "accumulatedMeltingSnowHeatFlux": "Qsnow",
            "accumulatedMeltingIceRunoffHeatFlux": "Qirunoff",
            "accumulatedIcebergHeatFlux": "Qberg",
            "accumulatedLandIceHeatFlux": "Qli",
            "accumulatedRainTemperatureFlux": "Qrain",
            "accumulatedRiverRunoffTemperatureFlux": "Qriver",
            "accumulatedEvapTemperatureFlux": "Qevap",
            "accumulatedSeaIceTemperatureFlux": "QsiT",
            "accumulatedIcebergTemperatureFlux": "QbergT",
            # Freshwater
            "netFreshwaterInput": "FWnet",
            "accumulatedEvaporationFlux": "FW_evap",
            "accumulatedRainFlux": "FW_rain",
            "accumulatedRiverRunoffFlux": "FW_river",
            "accumulatedSnowFlux": "FW_snow",
            "accumulatedSeaIceFlux": "FW_si",
            "accumulatedIcebergFlux": "FW_berg",
            "accumulatedLandIceFlux": "FW_li",
            "accumulatedIceRunoffFlux": "FW_irun",
            "accumulatedFrazilFlux": "FW_frazil",
            # Salt
            "netSaltFlux": "Salt_net",
            "salinityRestoringFluxAvg": "Salt_rest",
            # Mass/Volume
            "netMassFlux": "Mass_net",
            "massChange": "Mass_chg",
            "totalVolumeChange": "Vol_chg",
            "volumeCellGlobal": "Vol_glb",
            # Diagnostics
            "absoluteEnergyError": "Err_energy",
            "absoluteSaltError": "Err_salt",
            "absoluteFreshWaterConservation": "Err_fw",
            "kineticEnergyCellAvg": "KE",
            "enstrophyAvg": "Enstrophy",
            "pressureAvg": "Pressure",
            "layerThicknessAvg": "Layer_thk",
            # State
            "temperatureAvg": "Temp",
            "temperatureFluxAvg": "Temp_flux",
            "salinityAvg": "Salt",
            # Flux averages
            "evaporationFluxAvg": "FW_evap_avg",
            "rainFluxAvg": "FW_rain_avg",
            "riverRunoffFluxAvg": "FW_river_avg",
            "snowFluxAvg": "FW_snow_avg",
            "icebergFreshWaterFluxAvg": "FW_berg_avg",
            "iceRunoffFluxAvg": "FW_irun_avg",
            "seaIceFreshWaterFluxAvg": "FW_si_avg",
        }

        var_dict: Dict[str, Dict[str, Any]] = {}
        for i, (stem, shortkey) in enumerate(shortname_map.items()):
            var_dict[shortkey] = {"name": stem, "index": i, "glob": f"{stem}_*.nc"}
        return var_dict

    # --------------------------
    # global area (optional)
    # --------------------------
    def _init_global_area_once(
        self,
        global_area_m2: Optional[float],
        mesh_path: Optional[str],
        exemplar: Optional[Tuple[str, str, str]],
    ) -> float:
        if global_area_m2 is not None:
            A = float(global_area_m2)
            if not np.isfinite(A) or A <= 0:
                raise ValueError(f"Invalid global_area_m2={A}")
            return A

        eng = self.opt.engine or "netcdf4"
        ch  = self.opt.chunks

        if mesh_path:
            with xr.open_dataset(mesh_path, engine=eng, chunks=ch) as ds:
                if "areaCell" not in ds:
                    raise KeyError(f"'areaCell' not found in mesh: {mesh_path}")
                A = float(ds["areaCell"].sum().compute().item())
            if not np.isfinite(A) or A <= 0:
                raise ValueError(f"Invalid global area from mesh: {A}")
            return A

        if exemplar:
            base_dir, exp, sub_dir = exemplar
            candidates = [
                os.path.join(base_dir, exp, sub_dir, "mpasTimeSeriesOcean.nc"),
                os.path.join(base_dir, exp, "post", "time_series", "mpasTimeSeriesOcean.nc"),
                os.path.join(base_dir, exp, sub_dir, "post", "time_series", "mpasTimeSeriesOcean.nc"),
            ]
            fn = next((p for p in candidates if os.path.exists(p)), None)
            if fn is None:
                raise FileNotFoundError(f"mpasTimeSeriesOcean.nc not found via exemplar={exemplar}")

            with xr.open_dataset(fn, engine=eng, chunks=ch) as ts:
                var = "timeMonthly_avg_areaCellGlobal"
                if var not in ts:
                    raise KeyError(f"{var} not found in {fn}")
                A = float(ts[var].median(skipna=True).compute().item())

            if not np.isfinite(A) or A <= 0:
                raise ValueError(f"Invalid global area from exemplar: {A}")
            return A

        raise ValueError("Provide one of: global_area_m2, mesh_path, or exemplar=(base_dir, exp, sub_dir).")

    # --------------------------
    # helpers: file discovery
    # --------------------------
    def _resolve_base(self, base_dir: str, exp: str, sub_dir: str) -> str:
        if os.path.isabs(sub_dir):
            return sub_dir
        if sub_dir.startswith(exp):
            return os.path.join(base_dir, sub_dir)
        return os.path.join(base_dir, exp, sub_dir)

    def _files_for_key(
        self,
        *,
        base_dir: str,
        exp: str,
        sub_dir: str,
        short_or_stem: str,
        year_min: int,
        year_max: int,
    ) -> List[str]:
        base = self._resolve_base(base_dir, exp, sub_dir)

        if short_or_stem in self._var_dict:
            stem = self._var_dict[short_or_stem]["name"]
            pat = self._var_dict[short_or_stem]["glob"]
        else:
            stem = short_or_stem
            pat = f"{stem}_*.nc"

        fns = glob.glob(os.path.join(base, pat))
        if not fns:
            return []

        win_lo = year_min * 100 + 1
        win_hi = year_max * 100 + 12

        kept: List[str] = []
        for p in fns:
            bn = os.path.basename(p)
            m = self._ym_re.search(bn)
            if not m:
                kept.append(p)
                continue
            start_ym = int(m.group(1))
            end_ym   = int(m.group(2))
            if (start_ym <= win_hi) and (end_ym >= win_lo):
                kept.append(p)

        def _start_key(pth: str):
            mm = self._ym_re.search(os.path.basename(pth))
            return (int(mm.group(1)) if mm else 10**12, os.path.basename(pth))

        kept.sort(key=_start_key)
        return kept

    # --------------------------
    # helpers: time normalization
    # --------------------------
    def _normalize_calendar(self, cal: Optional[str]) -> Optional[str]:
        if not isinstance(cal, str):
            return cal
        cal_norm = cal.strip().lower().replace("-", "_")
        if cal_norm in {"no_leap", "noleap"}:
            return "noleap"
        if cal_norm in {"gregorian", "standard"}:
            return "gregorian"
        if cal_norm in {"proleptic_gregorian", "proleptic-gregorian"}:
            return "proleptic_gregorian"
        return cal

    def _normalize_units(self, units: Optional[str]) -> Optional[str]:
        _TIME_RE = re.compile(r"(\d{2}):(\d{2}):(\d{1,2})(\.\d+)?$")
        if not isinstance(units, str):
            return units
        if " since " not in units:
            return units
        head, rest = units.split(" since ", 1)
        if " " in rest:
            date_part, time_part = rest.split(" ", 1)
            m = _TIME_RE.search(time_part.strip())
            if m:
                hh, mm, ss, frac = m.groups()
                ss = ss.zfill(2)
                frac = frac or ""
                time_norm = f"{hh}:{mm}:{ss}{frac}"
                rest = f"{date_part} {time_norm}"
            return f"{head} since {rest}"
        return units

    def _fix_time_attrs(self, var: xr.DataArray) -> dict:
        attrs = dict(var.attrs)
        if "calendar" in attrs:
            attrs["calendar"] = self._normalize_calendar(attrs["calendar"])
        if "units" in attrs:
            attrs["units"] = self._normalize_units(attrs["units"])
        return attrs

    def reassign_monthly_time_with_bounds(
        self,
        ds: xr.Dataset,
        *,
        time: str = "time",
        year_min: int,
        year_max: int,
        calendar: str,
        length_policy: str,
        verbose: bool,
    ) -> xr.Dataset:
        if time not in ds.coords and time not in ds.dims:
            raise KeyError(f"'{time}' not found in dataset dims/coords.")

        y0, y1 = int(year_min), int(year_max)
        if y1 < y0:
            y0, y1 = y1, y0

        cal = (calendar or "noleap").lower()
        cal_map = {
            "noleap": cftime.DatetimeNoLeap,
            "365_day": cftime.DatetimeNoLeap,
            "360_day": cftime.Datetime360Day,
            "gregorian": cftime.DatetimeGregorian,
            "standard": cftime.DatetimeGregorian,
            "proleptic_gregorian": cftime.DatetimeProlepticGregorian,
            "julian": cftime.DatetimeJulian,
        }
        cf_cls = cal_map.get(cal, cftime.DatetimeNoLeap)

        n_expected = (y1 - y0 + 1) * 12
        new_time_full = xr.date_range(
            start=cf_cls(y0, 1, 1),
            periods=n_expected,
            freq="MS",
            use_cftime=True,
        ).tolist()

        n_ds = ds.sizes[time]

        if n_ds == n_expected:
            new_time = new_time_full
            ds_out = ds
        elif length_policy == "error":
            raise ValueError(f"Length mismatch: dataset has {n_ds} but expected {n_expected} months.")
        elif length_policy in ("auto", "trim_ds"):
            if n_ds > n_expected:
                ds_out = ds.isel({time: slice(0, n_expected)})
                new_time = new_time_full
                if verbose:
                    print(f"[diag] Trimmed ds {n_ds} -> {n_expected}")
            else:
                if length_policy == "trim_ds":
                    raise ValueError(f"Dataset has {n_ds} months but needs {n_expected}")
                ds_out = ds
                new_time = new_time_full[:n_ds]
                if verbose:
                    print(f"[diag] Truncated new_time to {n_ds} months")
        elif length_policy == "truncate_new_time":
            if n_expected >= n_ds:
                ds_out = ds
                new_time = new_time_full[:n_ds]
            else:
                raise ValueError(f"Dataset longer than requested window: ds={n_ds}, expected={n_expected}")
        else:
            raise ValueError(f"Unknown length_policy='{length_policy}'")

        ds_out = ds_out.assign_coords({time: (time, np.asarray(new_time))}).sortby(time)
        ds_out[time].encoding = dict(ds_out[time].encoding)
        ds_out[time].encoding.setdefault("calendar", cal)

        if verbose:
            print(f"[diag] Reassigned '{time}' monthly: {y0}-01 -> {y1}-12 (calendar={cal})")
        return ds_out

    def _open_mfdataset_robust(self, fns: Sequence[str], *, year_min: int, year_max: int) -> xr.Dataset:
        ds = xr.open_mfdataset(
            list(fns),
            combine="by_coords",
            decode_times=False,
            engine=self.opt.engine or "netcdf4",
            chunks=self.opt.chunks,
        )

        time_name = "time" if "time" in ds.variables else next((n for n in ("Time", "t") if n in ds.variables), None)
        bnd_name  = next((n for n in ("time_bnds", "time_bounds") if n in ds.variables), None)

        if time_name is not None:
            ds = ds.copy(deep=False)
            ds[time_name].attrs = self._fix_time_attrs(ds[time_name])
        if bnd_name is not None:
            ds[bnd_name].attrs = self._fix_time_attrs(ds[bnd_name])
            if time_name is not None:
                ds[time_name].attrs.setdefault("bounds", bnd_name)

        coder = xr.coders.CFDatetimeCoder(use_cftime=True)
        ds = xr.decode_cf(ds, decode_times=coder, mask_and_scale=True, decode_coords="all")

        if self.opt.reassign_monthly_time:
            if time_name and time_name != "time" and time_name in ds.dims:
                ds = ds.rename({time_name: "time"})
            if "time" in ds.dims:
                ds = self.reassign_monthly_time_with_bounds(
                    ds,
                    time="time",
                    year_min=year_min,
                    year_max=year_max,
                    calendar=self.opt.calendar,
                    length_policy=self.opt.length_policy,
                    verbose=self.opt.verbose,
                )
        return ds

    # --------------------------
    # helpers: make 1D time series
    # --------------------------
    @staticmethod
    def _normalize_time_1d(da: xr.DataArray) -> xr.DataArray:
        tdim = "time" if "time" in da.dims else next((d for d in da.dims if d.lower().startswith("time")), da.dims[0])
        if tdim != "time":
            da = da.rename({tdim: "time"})
        for d in list(da.dims):
            if d != "time":
                da = da.mean(d, skipna=True)
        return da.sortby("time")

    # --------------------------
    # annual mean helper (weighted)
    # --------------------------
    @staticmethod
    def to_annual_mean(
        da: xr.DataArray,
        *,
        weighted: bool = True,
        drop_incomplete: bool = True,
        min_days: int = 360,
        ymin: Optional[int] = None,
        ymax: Optional[int] = None,
    ) -> xr.DataArray:
        if "time" not in da.dims:
            raise ValueError("to_annual_mean expects time dimension named 'time'")
        if not hasattr(da["time"], "dt"):
            raise ValueError("to_annual_mean requires decoded datetime/CF-time on 'time'")

        if weighted:
            days = xr.DataArray(da["time"].dt.days_in_month, dims="time")
            num = (da * days).groupby("time.year").sum("time", skipna=True)
            den = days.groupby("time.year").sum("time")

            if drop_incomplete:
                # den has dim 'year' and usually coords, but be robust
                keep = den.where(den >= min_days, drop=True).coords.get("year", None)
                if keep is None:
                    try:
                        keep = den.get_index("year").values
                    except Exception:
                        keep = np.arange(den.sizes["year"])
                num = num.sel(year=keep)
                den = den.sel(year=keep)

            ann = num / den
        else:
            ann = da.groupby("time.year").mean("time", skipna=True)

            if drop_incomplete:
                nmon = xr.ones_like(da).groupby("time.year").sum("time")
                keep = nmon.where(nmon >= 12, drop=True).coords.get("year", None)
                if keep is None:
                    try:
                        keep = nmon.get_index("year").values
                    except Exception:
                        keep = np.arange(nmon.sizes["year"])
                ann = ann.sel(year=keep)

        # --- robustly get the year labels BEFORE renaming ---
        year_vals = ann.coords.get("year", None)
        if year_vals is None:
            try:
                year_vals = ann.get_index("year").values
            except Exception:
                year_vals = np.arange(ann.sizes["year"])
        year_vals = np.asarray(year_vals, dtype=int)

        # rename dim -> time and explicitly assign coordinate values
        ann = ann.rename({"year": "time"}).assign_coords(time=("time", year_vals))

        if ymin is not None or ymax is not None:
            lo = int(ymin) if ymin is not None else int(ann["time"].values.min())
            hi = int(ymax) if ymax is not None else int(ann["time"].values.max())
            ann = ann.sel(time=slice(lo, hi))

        return ann


    # --------------------------
    # core: load one series (monthly)
    # --------------------------
    def _load_monthly_series(
        self,
        *,
        base_dir: str,
        exp: str,
        sub_dir: str,
        stem: str,
        year_min: int,
        year_max: int,
        per_area: bool,
    ) -> xr.DataArray:
        cache_key = (base_dir, exp, sub_dir, stem, year_min, year_max, per_area, self._A)
        if cache_key in self._series_cache:
            return self._series_cache[cache_key]
    
        fns = self._files_for_key(
            base_dir=base_dir, exp=exp, sub_dir=sub_dir,
            short_or_stem=stem, year_min=year_min, year_max=year_max
        )
        if not fns:
            raise FileNotFoundError(
                f"No files found for stem='{stem}' under {self._resolve_base(base_dir, exp, sub_dir)}"
            )
    
        ds = self._open_mfdataset_robust(fns, year_min=year_min, year_max=year_max)
        try:
            if stem not in ds.data_vars:
                raise KeyError(f"'{stem}' not found in opened dataset (stem='{stem}')")
    
            da = self._normalize_time_1d(ds[stem])
    
            if per_area:
                if self._A is None or (not np.isfinite(self._A)) or self._A <= 0:
                    raise ValueError("per_area=True requested but global area A is not initialized.")
                da = da / float(self._A)
    
            try:
                years = da["time"].dt.year.astype(int)
                da = da.where((years >= int(year_min)) & (years <= int(year_max)), drop=True)
            except Exception:
                pass
    
            # ✅ critical: materialize before closing underlying files (robust for dask/netCDF backends)
            da = da.load()
    
            self._series_cache[cache_key] = da
            return da
    
        finally:
            # ✅ always close even if an exception occurs
            ds.close()

    # --------------------------
    # derived variable evaluation
    # --------------------------
    def _eval_derived(
        self,
        *,
        base_dir: str,
        exp: str,
        sub_dir: str,
        key: str,
        year_min: int,
        year_max: int,
    ) -> xr.DataArray:
        recipe = self._derived.get(key)
        if recipe is None:
            raise KeyError(f"No derived recipe for '{key}'")

        if "terms" in recipe:
            out = None
            for subkey, coef in recipe["terms"]:
                da = self.extract_monthly_da(
                    base_dir=base_dir, exp=exp, sub_dir=sub_dir,
                    key=subkey, year_min=year_min, year_max=year_max,
                    per_area=False
                )
                term = float(coef) * da
                out = term if out is None else (out + term)
            return self._normalize_time_1d(out)

        if "expr" in recipe:
            allowed: Dict[str, xr.DataArray] = {}
            for sk in re.findall(r"[A-Za-z_]\w*", recipe["expr"]):
                if sk not in allowed:
                    allowed[sk] = self.extract_monthly_da(
                        base_dir=base_dir, exp=exp, sub_dir=sub_dir,
                        key=sk, year_min=year_min, year_max=year_max,
                        per_area=False
                    )
            out = eval(recipe["expr"], {"__builtins__": {}}, allowed)
            return self._normalize_time_1d(out)

        raise ValueError(f"Bad derived recipe for '{key}'. Use 'terms' or 'expr'.")

    # --------------------------
    # public API: return xarray DataArray (monthly)
    # --------------------------
    def extract_monthly_da(
        self,
        *,
        base_dir: str,
        exp: str,
        sub_dir: str,
        key: str,
        year_min: int,
        year_max: int,
        per_area: bool = False,
    ) -> xr.DataArray:
        if key in self._derived:
            da = self._eval_derived(base_dir=base_dir, exp=exp, sub_dir=sub_dir, key=key, year_min=year_min, year_max=year_max)
            if per_area:
                if self._A is None or (not np.isfinite(self._A)) or self._A <= 0:
                    raise ValueError("per_area=True requested but global area A is not initialized.")
                da = da / float(self._A)
            return da

        stem = self._var_dict[key]["name"] if key in self._var_dict else key
        return self._load_monthly_series(
            base_dir=base_dir, exp=exp, sub_dir=sub_dir,
            stem=stem, year_min=year_min, year_max=year_max,
            per_area=per_area
        )

    # --------------------------
    # public API: return numpy (monthly or annual)
    # --------------------------
    def extract_series(
        self,
        *,
        base_dir: str,
        exp: str,
        sub_dir: str,
        req: FluxSeriesRequest,
        year_min: int,
        year_max: int,
        annual_mean: bool = False,
        annual_weighted: bool = True,
        annual_drop_incomplete: bool = True,
        annual_min_days: int = 360,
    ) -> np.ndarray:
        da = self.extract_monthly_da(
            base_dir=base_dir, exp=exp, sub_dir=sub_dir,
            key=req.key, year_min=year_min, year_max=year_max,
            per_area=req.per_area
        )

        if annual_mean:
            ann = self.to_annual_mean(
                da,
                weighted=annual_weighted,
                drop_incomplete=annual_drop_incomplete,
                min_days=annual_min_days,
                ymin=year_min,
                ymax=year_max,
            )
            out = np.asarray(ann.values, float).reshape(-1)
        else:
            out = np.asarray(da.values, float).reshape(-1)

        out = out * float(req.scale)
        if not np.isfinite(out).any():
            raise ValueError(f"All-NaN series for exp='{exp}', key='{req.key}'")
        return out

    def extract_group(
        self,
        *,
        base_dir: str,
        exp: str,
        sub_dir: str,
        requests: Dict[str, FluxSeriesRequest],
        year_min: int,
        year_max: int,
        annual_mean: bool = False,
        annual_weighted: bool = True,
        annual_drop_incomplete: bool = True,
        annual_min_days: int = 360,
    ) -> Dict[str, np.ndarray]:
        out: Dict[str, np.ndarray] = {}
        for alias, req in requests.items():
            out[alias] = self.extract_series(
                base_dir=base_dir, exp=exp, sub_dir=sub_dir, req=req,
                year_min=year_min, year_max=year_max,
                annual_mean=annual_mean,
                annual_weighted=annual_weighted,
                annual_drop_incomplete=annual_drop_incomplete,
                annual_min_days=annual_min_days,
            )
        return out

    def extract_many(
        self,
        *,
        base_dir: str,
        exps: Sequence[str],
        subdirs: Dict[str, str],
        requests: Dict[str, FluxSeriesRequest],
        year_min: int,
        year_max: int,
        annual_mean: bool = False,
        annual_weighted: bool = True,
        annual_drop_incomplete: bool = True,
        annual_min_days: int = 360,
    ) -> Dict[str, Dict[str, np.ndarray]]:
        out: Dict[str, Dict[str, np.ndarray]] = {}
        for exp in exps:
            out[exp] = self.extract_group(
                base_dir=base_dir,
                exp=exp,
                sub_dir=subdirs[exp],
                requests=requests,
                year_min=year_min,
                year_max=year_max,
                annual_mean=annual_mean,
                annual_weighted=annual_weighted,
                annual_drop_incomplete=annual_drop_incomplete,
                annual_min_days=annual_min_days,
            )
        return out


# ==============================
# Registry + convenience wrappers (to match your __main__)
# ==============================
class FluxRegistry:
    """
    Minimal registry:
      - var_dict in extractor schema (shortkey -> stem + glob)
      - derived recipes store
      - config builder used by your __main__
    """

    def __init__(self):
        self._var_dict = FluxTimeSeriesExtractor.build_default_var_dict()
        self._derived: Dict[str, Dict[str, Any]] = {}

    def get_var_dict(self) -> Dict[str, Dict[str, Any]]:
        return self._var_dict

    def get_derived(self) -> Dict[str, Dict[str, Any]]:
        return self._derived

    def add_derived_terms(self, key: str, *, terms: Sequence[Tuple[str, float]]) -> None:
        self._derived[str(key)] = {"terms": [(str(k), float(c)) for k, c in terms]}

    def add_derived_expr(self, key: str, *, expr: str) -> None:
        self._derived[str(key)] = {"expr": str(expr)}

    # ---- config builder (FIXED) ----
    def build_extract_config(
        self,
        *,
        ts_window: YearRange,
        annual_mean: bool = False,
        annual_weighted: bool = True,
        annual_drop_incomplete: bool = True,
        annual_min_days: int = 360,
        calendar: str = "noleap",
        reassign_monthly_time: bool = True,
        length_policy: str = "auto",
        # keep API compatibility (ignored)
        year_token: Optional[str] = None,
        region_key: Optional[str] = None,
        region_dim: Optional[str] = None,
        region_names_var: Optional[str] = None,
    ) -> FluxExtractConfig:
        return FluxExtractConfig(
            ts_window=ts_window,
            annual_mean=annual_mean,
            annual_weighted=annual_weighted,
            annual_drop_incomplete=annual_drop_incomplete,
            annual_min_days=annual_min_days,
            calendar=calendar,
            reassign_monthly_time=reassign_monthly_time,
            length_policy=length_policy,
        )

def _build_requests_from_var_dict_and_derived(
    *,
    var_dict: Dict[str, Dict[str, Any]],
    derived: Dict[str, Dict[str, Any]],
    per_area_vars: Optional[Sequence[str]] = None,
) -> Dict[str, FluxSeriesRequest]:
    """
    Build extraction requests:
      - include ALL var_dict keys passed in
      - include ALL derived keys passed in
      - per_area=True for keys in per_area_vars
    """
    per_area_set = set(per_area_vars or ())
    reqs: Dict[str, FluxSeriesRequest] = {}

    for k in var_dict.keys():
        reqs[k] = FluxSeriesRequest(key=k, scale=1.0, per_area=(k in per_area_set))

    for k in derived.keys():
        reqs[k] = FluxSeriesRequest(key=k, scale=1.0, per_area=(k in per_area_set))

    return reqs

def extract_all_flux(
    *,
    base_dir: str,
    exps: Sequence[str],
    exp_subdirs: Dict[str, str],
    var_dict: Dict[str, Dict[str, Any]],
    derived: Optional[Dict[str, Dict[str, Any]]],
    cfg: FluxExtractConfig,
    engine: Optional[str] = "netcdf4",
    chunks: Optional[dict] = None,
    verbose: bool = False,
    per_area_vars: Optional[Sequence[str]] = None,   # <--- ADD
    global_area_m2: Optional[float] = None,
    mesh_path: Optional[str] = None,
    exemplar: Optional[Tuple[str, str, str]] = None,
) -> Dict[str, Dict[str, np.ndarray]]:
    """
    Batch extract for your __main__ pattern.

    Returns:
      data[exp][key] -> 1D numpy
        - annual if cfg.annual_mean True
        - monthly otherwise
    """
    year_min, year_max = cfg.ts_window
    derived = derived or {}

    opt = FluxExtractOptions(
        engine=engine,
        chunks=chunks,
        calendar=cfg.calendar,
        reassign_monthly_time=cfg.reassign_monthly_time,
        length_policy=cfg.length_policy,
        verbose=bool(verbose),
    )

    X = FluxTimeSeriesExtractor(
        options=opt,
        var_dict=var_dict,
        derived_recipes=derived,
        global_area_m2=global_area_m2,
        mesh_path=mesh_path,
        exemplar=exemplar,
    )

    requests = _build_requests_from_var_dict_and_derived(
        var_dict=var_dict,
        derived=derived,
        per_area_vars=per_area_vars,
    )
    out: Dict[str, Dict[str, np.ndarray]] = {}
    for exp in exps:
        out[exp] = X.extract_group(
            base_dir=base_dir,
            exp=exp,
            sub_dir=exp_subdirs[exp],
            requests=requests,
            year_min=year_min,
            year_max=year_max,
            annual_mean=cfg.annual_mean,
            annual_weighted=cfg.annual_weighted,
            annual_drop_incomplete=cfg.annual_drop_incomplete,
            annual_min_days=cfg.annual_min_days,
        )
    return out


# =============================================================================
# v2, used by
# the ocean-heat-uptake analysis. Differences from v1: renamed temperature-flux keys
# (QrainT/QriverT/QevapT), month-length fix (dt.daysinmonth), per-term per_area
# handling in derived quantities, verbose diagnostics.
# =============================================================================
class FluxTimeSeriesExtractorV2(FluxTimeSeriesExtractor):
    @staticmethod
    def build_default_var_dict() -> Dict[str, Dict[str, Any]]:
        shortname_map = {
            # Energy
            "netEnergyFlux": "Qnet_mpas",   # 0 in the first month after a restart; Qnet is derived (alt_global_timeseries)
            "accumulatedShortWaveHeatFlux": "Qsw",
            "accumulatedLongWaveHeatFluxDown": "Qlw_dn",
            "accumulatedLongWaveHeatFluxUp": "Qlw_up",
            "accumulatedLatentHeatFlux": "Qlh",
            "accumulatedSensibleHeatFlux": "Qsh",
            "accumulatedSeaIceHeatFlux": "Qsi",
            "accumulatedFrazilHeatFlux": "Qfrazil",
            "accumulatedMeltingSnowHeatFlux": "Qsnow",
            "accumulatedMeltingIceRunoffHeatFlux": "Qirunoff",
            "accumulatedIcebergHeatFlux": "Qberg",
            "accumulatedLandIceHeatFlux": "Qli",
            "accumulatedRainTemperatureFlux": "QrainT",
            "accumulatedRiverRunoffTemperatureFlux": "QriverT",
            "accumulatedEvapTemperatureFlux": "QevapT",
            "accumulatedSeaIceTemperatureFlux": "QsiT",
            "accumulatedIcebergTemperatureFlux": "QbergT",
            # Freshwater
            "netFreshwaterInput": "FWnet",
            "accumulatedEvaporationFlux": "FW_evap",
            "accumulatedRainFlux": "FW_rain",
            "accumulatedRiverRunoffFlux": "FW_river",
            "accumulatedSnowFlux": "FW_snow",
            "accumulatedSeaIceFlux": "FW_si",
            "accumulatedIcebergFlux": "FW_berg",
            "accumulatedLandIceFlux": "FW_li",
            "accumulatedIceRunoffFlux": "FW_irun",
            "accumulatedFrazilFlux": "FW_frazil",
            # Salt
            "netSaltFlux": "Salt_net",
            "salinityRestoringFluxAvg": "Salt_rest",
            # Mass/Volume
            "netMassFlux": "Mass_net",
            "massChange": "Mass_chg",
            "totalVolumeChange": "Vol_chg",
            "volumeCellGlobal": "Vol_glb",
            # Diagnostics
            "absoluteEnergyError": "Err_energy",
            "absoluteSaltError": "Err_salt",
            "absoluteFreshWaterConservation": "Err_fw",
            "kineticEnergyCellAvg": "KE",
            "enstrophyAvg": "Enstrophy",
            "pressureAvg": "Pressure",
            "layerThicknessAvg": "Layer_thk",
            # State
            "temperatureAvg": "Temp",
            "temperatureFluxAvg": "Temp_flux",
            "salinityAvg": "Salt",
            # Flux averages
            "evaporationFluxAvg": "FW_evap_avg",
            "rainFluxAvg": "FW_rain_avg",
            "riverRunoffFluxAvg": "FW_river_avg",
            "snowFluxAvg": "FW_snow_avg",
            "icebergFreshWaterFluxAvg": "FW_berg_avg",
            "iceRunoffFluxAvg": "FW_irun_avg",
            "seaIceFreshWaterFluxAvg": "FW_si_avg",
        }

        var_dict: Dict[str, Dict[str, Any]] = {}
        for i, (stem, shortkey) in enumerate(shortname_map.items()):
            var_dict[shortkey] = {"name": stem, "index": i, "glob": f"{stem}_*.nc"}
        return var_dict

    def _load_monthly_series(
        self,
        *,
        base_dir: str,
        exp: str,
        sub_dir: str,
        stem: str,
        year_min: int,
        year_max: int,
        per_area: bool,
    ) -> xr.DataArray:
        cache_key = (base_dir, exp, sub_dir, stem, year_min, year_max, per_area, self._A)
        if cache_key in self._series_cache:
            return self._series_cache[cache_key]

        fns = self._files_for_key(
            base_dir=base_dir, exp=exp, sub_dir=sub_dir,
            short_or_stem=stem, year_min=year_min, year_max=year_max
        )
        if not fns:
            raise FileNotFoundError(
                f"No files found for stem='{stem}' under {self._resolve_base(base_dir, exp, sub_dir)}"
            )

        ds = self._open_mfdataset_robust(fns, year_min=year_min, year_max=year_max)
        try:
            if stem not in ds.data_vars:
                raise KeyError(f"'{stem}' not found in opened dataset (stem='{stem}')")

            da = self._normalize_time_1d(ds[stem])

            if self.opt.verbose:
                try:
                    y0 = int(da["time"].dt.year.min().compute().item())
                    y1 = int(da["time"].dt.year.max().compute().item())
                    print(f"[diag] '{stem}' time years after open: {y0}..{y1} (n={da.sizes['time']})")
                except Exception:
                    print(f"[diag] '{stem}' time not datetime-like; dims={da.dims}, n={da.sizes.get('time', None)}")

            if per_area:
                if self._A is None or (not np.isfinite(self._A)) or self._A <= 0:
                    raise ValueError("per_area=True requested but global area A is not initialized.")
                da = da / float(self._A)

            # Window trimming by year is only meaningful if the time coordinate is year-like.
            # Keep behavior but make it safer: only apply if dt.year works.
            try:
                years = da["time"].dt.year.astype(int)
                da = da.where((years >= int(year_min)) & (years <= int(year_max)), drop=True)
            except Exception:
                pass

            # critical: materialize before closing underlying files (robust for dask/netCDF backends)
            da = da.load()

            self._series_cache[cache_key] = da
            return da

        finally:
            # always close even if an exception occurs
            ds.close()

    def _eval_derived(
        self,
        *,
        base_dir: str,
        exp: str,
        sub_dir: str,
        key: str,
        year_min: int,
        year_max: int,
        per_area: bool,
    ) -> xr.DataArray:
        recipe = self._derived.get(key)
        if recipe is None:
            raise KeyError(f"No derived recipe for '{key}'")
    
        if self.opt.verbose:
            print(f"\n[INFO] _eval_derived: key='{key}' exp='{exp}' window={year_min}..{year_max}")
            print(f"  recipe keys: {sorted(recipe.keys())}")
            print(f"  per_area: {per_area}")
    
        if "terms" in recipe:
            if self.opt.verbose:
                print("  mode: terms")
                print("  terms:")
                for subkey, coef in recipe["terms"]:
                    print(f"    {coef:+g} * {subkey}")
    
            out = None
            for subkey, coef in recipe["terms"]:
                stem = self._var_dict[subkey]["name"] if subkey in self._var_dict else subkey
                da = self._load_monthly_series(
                    base_dir=base_dir, exp=exp, sub_dir=sub_dir,
                    stem=stem, year_min=year_min, year_max=year_max,
                    per_area=per_area
                )
            
                term = float(coef) * da
                vmin = float(np.nanmin(term.values))
                vmax = float(np.nanmax(term.values))
                
                print(f"  term '{subkey} {coef}' result: dims={term.dims}, sizes={dict(term.sizes)} range: vmin={vmin}, vmax={vmax}")
                
                # only these are temperature-flux terms converted with rho*cp
                #if subkey in {"QriverT", "QrainT", "QevapT", "QsiT", "QbergT"}:
                #    if self._A is None or (not np.isfinite(self._A)) or self._A <= 0:
                #        raise ValueError(f"Need global area A to area-normalize rho*cp term '{subkey}' in derived '{key}'")
                #    term = term / float(self._A)
                #if self.opt.verbose:
                #        print(f"    term '{subkey}': applied rho*cp then /A (A={self._A:.3e} m^2)")
            
                out = term if out is None else (out + term)
                
            out = self._normalize_time_1d(out)
    
        elif "expr" in recipe:
            expr = recipe["expr"]
            toks = re.findall(r"[A-Za-z_]\w*", expr)
    
            if self.opt.verbose:
                print("  mode: expr")
                print(f"  expr: {expr}")
                print(f"  tokens: {sorted(set(toks))}")
    
            allowed: Dict[str, xr.DataArray] = {}
            for sk in toks:
                if sk not in allowed:
                    stem = self._var_dict[sk]["name"] if sk in self._var_dict else sk
                    allowed[sk] = self._load_monthly_series(
                        base_dir=base_dir, exp=exp, sub_dir=sub_dir,
                        stem=stem, year_min=year_min, year_max=year_max,
                        per_area=per_area
                    )
                    if self.opt.verbose:
                        da = allowed[sk]
                        try:
                            y0 = int(da["time"].dt.year.min().compute().item())
                            y1 = int(da["time"].dt.year.max().compute().item())
                            print(f"    loaded '{sk}': dims={da.dims}, ntime={da.sizes.get('time')}, years={y0}..{y1}")
                        except Exception:
                            print(f"    loaded '{sk}': dims={da.dims}, sizes={dict(da.sizes)}")
    
            out = eval(expr, {"__builtins__": {}}, allowed)
            out = self._normalize_time_1d(out)
    
        else:
            raise ValueError(f"Bad derived recipe for '{key}'. Use 'terms' or 'expr'.")
    
        # --- final sanity summary (verbose only) ---
        if self.opt.verbose:
            try:
                n = int(out.sizes.get("time", -1))
                y0 = int(out["time"].dt.year.min().compute().item())
                y1 = int(out["time"].dt.year.max().compute().item())
                frac_nan = float(np.isnan(out.values).mean())
                vmin = float(np.nanmin(out.values))
                vmax = float(np.nanmax(out.values))
                print(f"  derived '{key}' result: dims={out.dims}, ntime={n}, years={y0}..{y1}, "
                      f"nan_frac={frac_nan:.3f}, min={vmin:.3e}, max={vmax:.3e}")
            except Exception:
                print(f"  derived '{key}' result: dims={out.dims}, sizes={dict(out.sizes)}")
    
        return out

    def extract_monthly_da(
        self,
        *,
        base_dir: str,
        exp: str,
        sub_dir: str,
        key: str,
        year_min: int,
        year_max: int,
        per_area: bool = False,
    ) -> xr.DataArray:
        if key in self._derived:
            da = self._eval_derived(
                base_dir=base_dir, exp=exp, sub_dir=sub_dir,
                key=key, year_min=year_min, year_max=year_max,
                per_area=per_area
            )
            return da

        stem = self._var_dict[key]["name"] if key in self._var_dict else key
        return self._load_monthly_series(
            base_dir=base_dir, exp=exp, sub_dir=sub_dir,
            stem=stem, year_min=year_min, year_max=year_max,
            per_area=per_area
        )

    @staticmethod
    def to_annual_mean(
        da: xr.DataArray,
        *,
        weighted: bool = True,
        drop_incomplete: bool = True,
        min_days: int = 360,
        ymin: Optional[int] = None,
        ymax: Optional[int] = None,
    ) -> xr.DataArray:
        if "time" not in da.dims:
            raise ValueError("to_annual_mean expects time dimension named 'time'")
        if not hasattr(da["time"], "dt"):
            raise ValueError("to_annual_mean requires decoded datetime/CF-time on 'time'")

        if weighted:
            # BUGFIX: xarray accessor is daysinmonth (not days_in_month)
            days = da["time"].dt.daysinmonth.astype("float64")

            num = (da * days).groupby("time.year").sum("time", skipna=True)
            den = days.groupby("time.year").sum("time")

            if drop_incomplete:
                keep = den.where(den >= min_days, drop=True).coords.get("year", None)
                if keep is None:
                    try:
                        keep = den.get_index("year").values
                    except Exception:
                        keep = np.arange(den.sizes["year"])
                num = num.sel(year=keep)
                den = den.sel(year=keep)

            ann = num / den
        else:
            ann = da.groupby("time.year").mean("time", skipna=True)

            if drop_incomplete:
                nmon = xr.ones_like(da).groupby("time.year").sum("time")
                keep = nmon.where(nmon >= 12, drop=True).coords.get("year", None)
                if keep is None:
                    try:
                        keep = nmon.get_index("year").values
                    except Exception:
                        keep = np.arange(nmon.sizes["year"])
                ann = ann.sel(year=keep)

        # --- robustly get the year labels BEFORE renaming ---
        year_vals = ann.coords.get("year", None)
        if year_vals is None:
            try:
                year_vals = ann.get_index("year").values
            except Exception:
                year_vals = np.arange(ann.sizes["year"])
        year_vals = np.asarray(year_vals, dtype=int)

        # rename dim -> time and explicitly assign coordinate values
        ann = ann.rename({"year": "time"}).assign_coords(time=("time", year_vals)).sortby("time")

        if ymin is not None or ymax is not None:
            lo = int(ymin) if ymin is not None else int(ann["time"].values.min())
            hi = int(ymax) if ymax is not None else int(ann["time"].values.max())
            ann = ann.sel(time=slice(lo, hi))

        return ann


class FluxRegistryV2(FluxRegistry):
    def __init__(self):
        self._var_dict = FluxTimeSeriesExtractorV2.build_default_var_dict()
        self._derived: Dict[str, Dict[str, Any]] = {}


def extract_all_flux_v2(
    *,
    base_dir: str,
    exps: Sequence[str],
    exp_subdirs: Dict[str, str],
    var_dict: Dict[str, Dict[str, Any]],
    derived: Optional[Dict[str, Dict[str, Any]]],
    cfg: FluxExtractConfig,
    engine: Optional[str] = "netcdf4",
    chunks: Optional[dict] = None,
    verbose: bool = False,
    per_area_vars: Optional[Sequence[str]] = None,
    global_area_m2: Optional[float] = None,
    mesh_path: Optional[str] = None,
    exemplar: Optional[Tuple[str, str, str]] = None,
) -> Dict[str, Dict[str, np.ndarray]]:
    """
    Batch extract for your __main__ pattern.

    Returns:
      data[exp][key] -> 1D numpy
        - annual if cfg.annual_mean True
        - monthly otherwise
    """
    year_min, year_max = cfg.ts_window
    derived = derived or {}

    opt = FluxExtractOptions(
        engine=engine,
        chunks=chunks,
        calendar=cfg.calendar,
        reassign_monthly_time=cfg.reassign_monthly_time,
        length_policy=cfg.length_policy,
        verbose=bool(verbose),
    )

    X = FluxTimeSeriesExtractorV2(
        options=opt,
        var_dict=var_dict,
        derived_recipes=derived,
        global_area_m2=global_area_m2,
        mesh_path=mesh_path,
        exemplar=exemplar,
    )

    requests = _build_requests_from_var_dict_and_derived(
        var_dict=var_dict,
        derived=derived,
        per_area_vars=per_area_vars,
    )

    out: Dict[str, Dict[str, np.ndarray]] = {}
    for exp in exps:
        out[exp] = X.extract_group(
            base_dir=base_dir,
            exp=exp,
            sub_dir=exp_subdirs[exp],
            requests=requests,
            year_min=year_min,
            year_max=year_max,
            annual_mean=cfg.annual_mean,
            annual_weighted=cfg.annual_weighted,
            annual_drop_incomplete=cfg.annual_drop_incomplete,
            annual_min_days=cfg.annual_min_days,
        )
    return out


FluxExtractOptionsV2 = FluxExtractOptions      # identical in v1 and v2
FluxSeriesRequestV2 = FluxSeriesRequest        # identical in v1 and v2
