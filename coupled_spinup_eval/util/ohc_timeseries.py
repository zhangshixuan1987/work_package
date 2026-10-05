"""Ocean heat content time series and depth-time (Hovmoller) extraction from MPAS-Analysis output."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple, Sequence, Any, List, Union

import os
import re
import glob
import numpy as np
import xarray as xr
import cftime

YearRange = Tuple[int, int]
DepthRange = Tuple[float, float]  # (zmin, zmax) meters, positive-down

# ==============================
# Options / Requests
# ==============================
@dataclass
class OHCExtractOptions:
    engine: Optional[str] = "netcdf4"
    chunks: Optional[dict] = None

    calendar: str = "noleap"
    reassign_monthly_time: bool = True
    length_policy: str = "auto"  # "auto" | "error" | "trim_ds" | "truncate_new_time"
    verbose: bool = False


@dataclass
class OHCSeriesRequest:
    """
    One extracted series.

    kind:
      - "OHC": ocean heat content anomaly/absolute integrated over depth window (J or ZJ)
      - "OHU": inferred heat uptake = d(OHC)/dt / area (W m-2) (derived from OHC tendency)

    layer:
      depth window in meters (positive-down), [zmin, zmax)
    region_index:
      integer MPAS ocean-region index along region_dim (usually nOceanRegions)

    anomaly:
      if True, subtract baseline mean (first baseline_months)
    to_ZJ:
      if True, return OHC in 1e22 J
    rho, cp:
      constants

    use_mask_in_volume:
      True: V ~ (Aavg * Nmask) * overlap
      False: V ~ Aavg * overlap
      (defaults True; also carried consistently into OHU area normalization)
    """
    kind: str = "OHC"
    layer: DepthRange = (0.0, 700.0)

    # region selection
    region_index: int = 0
    region_dim: Optional[str] = "nOceanRegions"
    region_name: Optional[str] = None
    region_names_var: Optional[str] = "regionNames"

    # anomaly / units
    anomaly: bool = False
    baseline_months: int = 12
    to_ZJ: bool = True

    # physics
    rho: float = 1026.0
    cp: float = 3996.0

    # volume/area interpretation
    use_mask_in_volume: bool = True

    # OHU options
    smooth_months: Optional[int] = None
    post_smooth_months: Optional[int] = None
    area_tol_rel: float = 1e-6


@dataclass
class OHCHovmollerRequest:
    """
    Extract a 2D OHC matrix for Hovmöller plots.

    mode:
      - "per_area": J m-2 per layer  (rho*cp*T*Hlayer)
      - "total":    J per layer      (per_area * effective_layer_area)

    to_units:
      - "J" (for total) or "J m^-2" (for per_area)  [no scaling]
      - "ZJ" / "1e22 J"      (for total only)
      - "GJ/m^2"             (for per_area only)

    anomaly:
      subtract baseline mean (calendar_year or first_N months)

    normalize_by_surface_area:
      If mode="total", divide by surface (lev=0) region area to return J m^-2.

    apply_window_overlap:
      If True and window is not None, use overlap thickness for the window [zmin,zmax)
      (partial boundary layers) before forming per-layer OHC.
      This is mainly for closure tests against the 1D window-integrated OHC.
    """
    mode: str = "total"                  # "total" or "per_area"
    anomaly: bool = True
    baseline_mode: str = "calendar_year" # "calendar_year" or "first_N"
    first_N: int = 12

    to_units: str = "ZJ"                 # "ZJ"/"1e22 J" or "GJ/m^2" or "J"

    # physics
    rho: float = 1026.0
    cp: float = 3996.0

    # region selection
    region_index: int = 0
    region_dim: Optional[str] = "nOceanRegions"
    region_name: Optional[str] = None
    region_names_var: Optional[str] = "regionNames"

    # depth coordinate choice
    depth_coord: str = "bottom"          # "bottom" | "mid"

    # cumulative option
    return_cumulative: bool = False
    from_top: bool = True

    # legacy-compatible option
    normalize_by_surface_area: bool = False

    # optional closure option vs 1D window-integrated OHC
    window: Optional[DepthRange] = None
    apply_window_overlap: bool = False

    # option to use mask in area 
    use_mask_in_area: bool = True


# ==============================
# Batch extraction config
# ==============================
@dataclass
class OHCExtractConfig:
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
class OHCTimeSeriesExtractor:
    """
    MPAS-O OHC extractor.
    """

    _ym_re = re.compile(r"_(\d{6})-(\d{6})")

    def __init__(
        self,
        *,
        options: Optional[OHCExtractOptions] = None,
        file_key: str = "mpasTimeSeriesOcean",
        fdepth: Optional[str] = None,   # must contain refBottomDepth
    ):
        self.opt = options or OHCExtractOptions()
        self.file_key = str(file_key)
        self.fdepth = fdepth

        self._series_cache: Dict[Tuple, Union[xr.DataArray, xr.Dataset]] = {}
        self._depth_cache: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}

    # --------------------------
    # file discovery
    # --------------------------
    def _resolve_base(self, base_dir: str, exp: str, sub_dir: str) -> str:
        if os.path.isabs(sub_dir):
            return sub_dir
        if sub_dir.startswith(exp):
            return os.path.join(base_dir, sub_dir)
        return os.path.join(base_dir, exp, sub_dir)

    def _find_timeseries_files(self, *, base_dir: str, exp: str, sub_dir: str) -> List[str]:
        base = self._resolve_base(base_dir, exp, sub_dir)

        single = glob.glob(os.path.join(base, f"{self.file_key}.nc"))
        if single:
            return single

        segs = glob.glob(os.path.join(base, f"{self.file_key}_*.nc"))
        if not segs:
            return []

        def start_stamp(path: str) -> int:
            m = re.search(
                rf"{re.escape(self.file_key)}_(\d{{6}})-(\d{{6}})\.nc$",
                os.path.basename(path),
            )
            return int(m.group(1)) if m else 0

        segs.sort(key=start_stamp)
        return segs

    def _files_for_window(
        self,
        *,
        base_dir: str,
        exp: str,
        sub_dir: str,
        year_min: int,
        year_max: int,
    ) -> List[str]:
        fns = self._find_timeseries_files(base_dir=base_dir, exp=exp, sub_dir=sub_dir)
        if not fns:
            return []

        if len(fns) == 1 and fns[0].endswith(f"{self.file_key}.nc"):
            return fns

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
            end_ym = int(m.group(2))
            if (start_ym <= win_hi) and (end_ym >= win_lo):
                kept.append(p)

        def _start_key(pth: str):
            mm = self._ym_re.search(os.path.basename(pth))
            return (int(mm.group(1)) if mm else 10**12, os.path.basename(pth))

        kept.sort(key=_start_key)
        return kept

    def _depth_coord_from_fdepth(self, *, lev_dim: str, kind: str = "bottom") -> xr.DataArray:
        if not self.fdepth or (not os.path.exists(self.fdepth)):
            raise FileNotFoundError(f"fdepth not available: {self.fdepth}")

        z_top_vals, z_bot_vals = self._get_ref_edges(self.fdepth)

        if str(kind).lower() == "mid":
            z = 0.5 * (z_top_vals + z_bot_vals)
        else:
            z = z_bot_vals

        depth = xr.DataArray(z, dims=[lev_dim], name="depth")
        depth.attrs.update(units="m", positive="down", long_name="depth below sea surface")
        return depth

    # --------------------------
    # time repair helpers
    # --------------------------
    @staticmethod
    def _normalize_calendar(cal: Optional[str]) -> Optional[str]:
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

    @staticmethod
    def _normalize_units(units: Optional[str]) -> Optional[str]:
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
            ds_out = ds
            new_time = new_time_full
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
        bnd_name = next((n for n in ("time_bnds", "time_bounds") if n in ds.variables), None)

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
    # annual mean helper (robust year-coord fix)
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

        year_vals = ann.coords.get("year", None)
        if year_vals is None:
            try:
                year_vals = ann.get_index("year").values
            except Exception:
                year_vals = np.arange(ann.sizes["year"])
        year_vals = np.asarray(year_vals, dtype=int)

        ann = ann.rename({"year": "time"}).assign_coords(time=("time", year_vals))

        if ymin is not None or ymax is not None:
            lo = int(ymin) if ymin is not None else int(ann["time"].values.min())
            hi = int(ymax) if ymax is not None else int(ann["time"].values.max())
            ann = ann.sel(time=slice(lo, hi))

        return ann

    # --------------------------
    # depth geometry helpers
    # --------------------------
    def _get_ref_edges(self, fdepth: str) -> Tuple[np.ndarray, np.ndarray]:
        """
        Returns (z_top, z_bot) arrays of length nlev (meters, positive-down).
        Cached by fdepth path.

        NOTE: We only sort if the file is not already monotonic increasing, to avoid
        accidentally breaking alignment with lev indexing.
        """
        if fdepth in self._depth_cache:
            return self._depth_cache[fdepth]

        if not os.path.exists(fdepth):
            raise FileNotFoundError(f"fdepth not found: {fdepth}")

        with xr.open_dataset(fdepth) as fds:
            if "refBottomDepth" not in fds:
                raise KeyError(f"refBottomDepth not found in {fdepth}")
            z_bot_ref = fds["refBottomDepth"].squeeze()

        z_bot_vals = np.asarray(z_bot_ref.values, dtype=float)

        # ensure positive-down
        if np.nanmax(z_bot_vals) <= 0:
            z_bot_vals = -z_bot_vals

        # only sort if non-monotonic
        if np.any(np.diff(z_bot_vals) < 0):
            z_bot_vals = np.sort(z_bot_vals)

        z_top_vals = np.concatenate(([0.0], z_bot_vals[:-1]))

        self._depth_cache[fdepth] = (z_top_vals, z_bot_vals)
        return z_top_vals, z_bot_vals

    @staticmethod
    def _overlap_from_edges(z_top: xr.DataArray, z_bot: xr.DataArray, zmin: float, zmax: float) -> xr.DataArray:
        dz = (z_bot - z_top)
        out = xr.apply_ufunc(
            lambda zt, zb, dz_: np.clip(np.minimum(zb, zmax) - np.maximum(zt, zmin), 0.0, dz_),
            z_top, z_bot, dz,
            dask="parallelized",
            output_dtypes=[float],
        )
        return out.where((z_bot > z_top))

    def _overlap_from_fdepth(
        self,
        Hlike: xr.DataArray,
        *,
        lev_dim: str,
        zmin: float,
        zmax: float,
    ) -> xr.DataArray:
        if not self.fdepth or (not os.path.exists(self.fdepth)):
            raise FileNotFoundError(f"fdepth not available: {self.fdepth}")

        z_top_vals, z_bot_vals = self._get_ref_edges(self.fdepth)

        z_top_ref = xr.DataArray(z_top_vals, dims=[lev_dim])
        z_bot_ref = xr.DataArray(z_bot_vals, dims=[lev_dim])

        if lev_dim in Hlike.coords:
            tgt = Hlike[lev_dim]
            z_top_ref = z_top_ref.assign_coords({lev_dim: tgt})
            z_bot_ref = z_bot_ref.assign_coords({lev_dim: tgt})

        z_top = z_top_ref.broadcast_like(Hlike).astype(float)
        z_bot = z_bot_ref.broadcast_like(Hlike).astype(float)

        return self._overlap_from_edges(z_top, z_bot, zmin, zmax)

    @staticmethod
    def _smooth(da: xr.DataArray, time_dim: str, win: Optional[int]) -> xr.DataArray:
        if win and int(win) > 1:
            w = int(win)
            return da.rolling({time_dim: w}, center=True, min_periods=max(2, w // 2)).mean()
        return da

    # --------------------------
    # core OHC computation
    # --------------------------
    def _compute_ohc_monthly(
        self,
        ds: xr.Dataset,
        *,
        region_index: int,
        region_dim: Optional[str],
        layer: DepthRange,
        rho: float,
        cp: float,
        anomaly: bool,
        baseline_months: int,
        to_ZJ: bool,
        temp_name: str = "timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature",
        thk_name: str  = "timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerThickness",
        area_name: str = "timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerArea",
        mask_name: str = "timeMonthly_avg_avgValueWithinOceanLayerRegion_sumLayerMaskValue",
        use_mask_in_volume: bool = True,
    ) -> xr.DataArray:
        for req in (temp_name, thk_name, area_name, mask_name):
            if req not in ds:
                raise KeyError(f"Required variable '{req}' not found in dataset.")

        T0 = ds[temp_name]
        time_dim = "time" if "time" in T0.dims else T0.dims[0]

        lev_dim = next((d for d in ("nVertLevels", "z", "nz", "nLevels", "lev") if d in ds.dims), None)
        if lev_dim is None:
            raise ValueError(f"Cannot find vertical dim among {tuple(ds.dims)}")

        if region_dim and region_dim in ds.dims:
            reg_dim = region_dim
        else:
            reg_dim = next((d for d in ("nOceanRegionsTmp", "nOceanRegions", "region", "nRegions") if d in ds.dims), None)
            if reg_dim is None:
                raise ValueError(f"Cannot find region dim among {tuple(ds.dims)}")

        zmin, zmax = float(layer[0]), float(layer[1])
        if not np.isfinite(zmax) or zmax <= zmin:
            raise ValueError(f"Invalid depth window layer={layer}")

        # region slice
        T    = ds[temp_name].isel({reg_dim: int(region_index)})
        Hraw = ds[thk_name].isel({reg_dim: int(region_index)})
        Aavg = ds[area_name].isel({reg_dim: int(region_index)})
        Nmsk = ds[mask_name].isel({reg_dim: int(region_index)})  
        fac = float(rho) * float(cp)

        # no wet-mask of T/H/A (trust the MPAS region aggregation)
        if self.fdepth and os.path.exists(self.fdepth):
            overlap = self._overlap_from_fdepth(Hraw, lev_dim=lev_dim, zmin=zmin, zmax=zmax)
            geom_src = "refBottomDepth"
        else:
            H = Hraw.clip(min=0.0)
            z_bot = H.cumsum(dim=lev_dim)
            z_top = (z_bot - H).clip(min=0.0)
            overlap = self._overlap_from_edges(z_top, z_bot, zmin, zmax)
            geom_src = "layerThickness"

        vol_fac = (Aavg * Nmsk) if use_mask_in_volume else Aavg
        Vlayer = (vol_fac * overlap)
        ts = (fac * T * Vlayer).sum(dim=lev_dim, skipna=True).rename("OHC")

        if anomaly:
            n0 = int(max(1, min(int(baseline_months), ts.sizes.get(time_dim, ts.shape[0]))))
            base = ts.isel({time_dim: slice(0, n0)}).mean(time_dim, skipna=True)
            ts = (ts - base).rename("OHC")

        if to_ZJ:
            ts = ts / 1.0e22
            units = "1e22 J"
        else:
            units = "J"

        ts = ts.assign_attrs(
            units=units,
            long_name=f"OHC{' anomaly' if anomaly else ''} ({zmin:g}–{zmax:g} m)",
            rho0=float(rho),
            cp0=float(cp),
            factor_used=fac,
            factor_meaning="rho*cp",
            depth_min_m=zmin,
            depth_max_m=zmax,
            region_index=int(region_index),
            use_mask_in_volume=bool(use_mask_in_volume),
            geometry_source=str(geom_src),
        )

        if time_dim != "time" and time_dim in ts.dims:
            ts = ts.rename({time_dim: "time"})
        return ts

    def _infer_ohu_from_ohc(
        self,
        ohc: xr.DataArray,
        *,
        ds_for_area: xr.Dataset,
        region_index: int,
        region_dim: Optional[str],
        lev_dim: str,
        smooth_months: Optional[int],
        post_smooth_months: Optional[int],
        area_tol_rel: float,
        use_mask_in_area: bool = True,
        area_name: str = "timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerArea",
        mask_name: str = "timeMonthly_avg_avgValueWithinOceanLayerRegion_sumLayerMaskValue",
    ) -> xr.Dataset:
    
        # --- region dim discovery (unchanged) ---
        if region_dim and region_dim in ds_for_area.dims:
            reg_dim = region_dim
        else:
            reg_dim = next((d for d in ("nOceanRegionsTmp", "nOceanRegions", "region", "nRegions") if d in ds_for_area.dims), None)
            if reg_dim is None:
                raise ValueError("Cannot find region dim for area computation.")
    
        if area_name not in ds_for_area or mask_name not in ds_for_area:
            raise KeyError(f"Need {area_name} and {mask_name} in dataset for OHU area normalization.")
    
        Aavg  = ds_for_area[area_name].isel({reg_dim: int(region_index)})
        Nmask = ds_for_area[mask_name].isel({reg_dim: int(region_index)})
        
        # --- INFO: region & dimension sanity check ---
        print("\n[INFO] _infer_ohu_from_ohc: region / dimension check")
        print(f"  Requested region_dim : {region_dim}")
        print(f"  Detected  region_dim : {reg_dim}")
        print(f"  Using region_index  : {int(region_index)}")
        
        print(f"  area_name dims      : {ds_for_area[area_name].dims}")
        print(f"  mask_name dims      : {ds_for_area[mask_name].dims}")
        
        print(f"  Aavg dims after sel : {Aavg.dims}")
        print(f"  Nmask dims after sel: {Nmask.dims}")
        
        # time sanity
        if "time" in Aavg.dims:
            print(f"  time size (area)    : {Aavg.sizes['time']}")
        else:
            print("  WARNING: 'time' not found in area variable dims")
        
        # vertical sanity
        if lev_dim in Aavg.dims:
            print(f"  lev_dim '{lev_dim}' size : {Aavg.sizes[lev_dim]}")
        else:
            print(f"  WARNING: lev_dim '{lev_dim}' not in area variable dims")
        
        # quick magnitude check (surface layer, first time)
        try:
            A0 = (Aavg * Nmask).isel({lev_dim: 0, "time": 0})
            print(f"  Surface area sample : {float(A0.mean().compute().item()):.3e} m^2")
        except Exception as e:
            print(f"  NOTE: could not compute surface area sample ({e})")

        tcoord = ohc["time"]
        
        # robust time-in-seconds for np.gradient (handles cftime safely)
        is_cftime = (tcoord.size > 0) and isinstance(tcoord.values[0], cftime.datetime)
        
        if (not is_cftime):
            try:
                t0 = tcoord.astype("datetime64[ns]")
                tsec = (t0 - t0.isel(time=0)) / np.timedelta64(1, "s")
                tsec = tsec.astype("float64")
            except Exception:
                is_cftime = True
        
        if is_cftime:
            try:
                days = tcoord.dt.daysinmonth.astype("float64")
            except Exception:
                days = xr.full_like(tcoord, 30, dtype="float64")
            # month-center or month-start; either is fine for scale, keep simple:
            tsec = (days.cumsum("time") - days.isel(time=0)) * 86400.0
        
        tsec = xr.DataArray(tsec, dims=["time"], coords={"time": tcoord})

        # --- ensure OHC is in Joules before smoothing ---
        ohc_J = ohc
        units = str(ohc.attrs.get("units", "")).strip()
        if units.startswith("1e"):
            # e.g. "1e22 J"
            scale = float(units.split()[0])
            ohc_J = ohc * scale
            ohc_J.attrs["units"] = "J"
        elif units.lower() in ("j", "joule", "joules"):
            pass  # already in J
        else:
            raise ValueError(f"Unrecognized OHC units: '{units}'")
    
        # --- smooth pre-derivative (same idea as original) ---
        E = self._smooth(ohc_J, "time", smooth_months)
    
        # Ensure core dim is single chunk for dask-parallelized gufunc
        if getattr(E.data, "chunks", None) is not None:
            E = E.chunk({"time": -1})
        if getattr(tsec.data, "chunks", None) is not None:
            tsec = tsec.chunk({"time": -1})
    
        dEdt_W = xr.apply_ufunc(
            np.gradient, E, tsec,
            input_core_dims=[["time"], ["time"]],
            output_core_dims=[["time"]],
            vectorize=True,
            dask="parallelized",
            output_dtypes=[E.dtype],
        )  # W (region total if OHC is J)
    
        # --- AREA NORMALIZATION: make identical to infer_flux_from_layer_ohc ---
        A_base = (Aavg * Nmask) if use_mask_in_area else Aavg
        A_t = A_base.isel({lev_dim: 0})  # (time, ...)
    
        # Collapse over ALL dims (exactly like original)
        A_mean_da = A_t.mean(dim=list(A_t.dims), skipna=True)
        A_min_da  = A_t.min(dim=list(A_t.dims),  skipna=True)
        A_max_da  = A_t.max(dim=list(A_t.dims),  skipna=True)
    
        # robust scalar extraction (handles dask and non-dask)
        def _to_float_scalar(x: xr.DataArray) -> float:
            x2 = x.compute() if hasattr(x.data, "compute") else x
            return float(x2.item())
    
        A_mean = _to_float_scalar(A_mean_da)
        A_min  = _to_float_scalar(A_min_da)
        A_max  = _to_float_scalar(A_max_da)
    
        A_den = A_mean if (A_mean > 0.0) else max(A_max, 1.0)
        rel_range = (A_max - A_min) / A_den if A_den > 0 else np.inf
    
        if np.isfinite(rel_range) and rel_range <= float(area_tol_rel):
            # match original: scalar from time=0
            A_region = _to_float_scalar(A_t.isel(time=0))
            OHU_raw = (dEdt_W / A_region).rename("OHU_raw")
            area_mode = "constant"
        else:
            OHU_raw = (dEdt_W / A_t).rename("OHU_raw")
            area_mode = "time_varying"
    
        OHU = self._smooth(OHU_raw, "time", post_smooth_months).rename("OHU")
    
        OHU_raw = OHU_raw.assign_attrs(
            units="W m-2",
            long_name="Inferred OHU (raw)",
            area_mode=area_mode,
            area_rel_range=float(rel_range),
            use_mask_in_area=bool(use_mask_in_area),
            area_tol_rel=float(area_tol_rel),
        )
        OHU = OHU.assign_attrs(
            units="W m-2",
            long_name="Inferred OHU (smoothed)" if (post_smooth_months and post_smooth_months > 1) else "Inferred OHU",
            area_mode=area_mode,
            area_rel_range=float(rel_range),
            use_mask_in_area=bool(use_mask_in_area),
            area_tol_rel=float(area_tol_rel),
        )
    
        return xr.Dataset(dict(OHU=OHU, OHU_raw=OHU_raw, OHC=ohc.rename("OHC")))

    def _compute_ohc_hovmoller_monthly(
        self,
        ds: xr.Dataset,
        *,
        region_index: int,
        region_dim: Optional[str],
        rho: float,
        cp: float,
        anomaly: bool,
        baseline_mode: str,
        first_N: int,
        mode: str,
        to_units: str,
        depth_coord: str,
        return_cumulative: bool,
        from_top: bool,
        temp_name: str = "timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature",
        thk_name: str  = "timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerThickness",
        area_name: str = "timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerArea",
        mask_name: str = "timeMonthly_avg_avgValueWithinOceanLayerRegion_sumLayerMaskValue",
        use_mask_in_area: bool = True,
        # --- NEW (optional) ---
        normalize_by_surface_area: bool = False,
        window: Optional[DepthRange] = None,
        apply_window_overlap: bool = False,
    ) -> Union[xr.DataArray, xr.Dataset]:
    
        for req in (temp_name, thk_name, area_name, mask_name):
            if req not in ds:
                raise KeyError(f"Required variable '{req}' not found in dataset.")
    
        # robust dims
        T0 = ds[temp_name]
        time_dim = "time" if "time" in T0.dims else T0.dims[0]
    
        lev_dim = next((d for d in ("nVertLevels", "z", "nz", "nLevels", "lev") if d in ds.dims), None)
        if lev_dim is None:
            raise ValueError(f"Cannot find vertical dim among {tuple(ds.dims)}")
    
        if region_dim and region_dim in ds.dims:
            reg_dim = region_dim
        else:
            reg_dim = next((d for d in ("nOceanRegionsTmp", "nOceanRegions", "region", "nRegions") if d in ds.dims), None)
            if reg_dim is None:
                raise ValueError(f"Cannot find region dim among {tuple(ds.dims)}")
                
        # region slice
        T    = ds[temp_name].isel({reg_dim: int(region_index)})
        Hraw = ds[thk_name].isel({reg_dim: int(region_index)})
        Aavg = ds[area_name].isel({reg_dim: int(region_index)})
        Nmsk = ds[mask_name].isel({reg_dim: int(region_index)})  
        fac = rho * cp 

        if self.opt.verbose:
            nreg = ds.sizes.get(reg_dim, None)
            if nreg is not None and not (0 <= int(region_index) < nreg):
                print(f"[diag][region] WARNING: region_index={region_index} out of bounds for '{reg_dim}' (size={nreg})")
            else:
                print(f"[diag][region] using region_index={region_index} (dim='{reg_dim}', size={nreg})")
        
            try:
                wet_frac = float((Nmsk > 0).mean().compute().item())
            except Exception:
                wet_frac = np.nan
        
            try:
                surf_area = float(Aavg.isel({lev_dim: 0}).mean().compute().item())
            except Exception:
                surf_area = np.nan
        
            print(f"[diag][region] wet_fraction={wet_frac:.3f} mean_surface_area={surf_area:.3e}")

        if self.fdepth and os.path.exists(self.fdepth):
            geom_src = "refBottomDepth"
        else:
            geom_src = "layerThickness"
          
        # ohc per layer
        ohc_pa = (fac * T * Hraw).rename("OHC_layer")  # J/m^2 per layer
        ohc_pa.attrs.update(units="J m^-2", geometry_source=str(geom_src))

        # total OHC per layer (J): multiply by effective wet area
        Aeff = (Aavg * Nmsk) if use_mask_in_area else Aavg 
        ohc_total  = (ohc_pa * Aeff).rename("OHC_layer") # J per layer
        ohc_total.attrs.update(units="J", geometry_source=str(geom_src))
    
        mode_l = str(mode).lower()
        if mode_l == "per_area":
            ohc = ohc_pa
            intended = "J m^-2"
        elif mode_l == "total":
            ohc = ohc_total
            intended = "J"
        else:
            raise ValueError("mode must be 'per_area' or 'total'")
            
        # optional: convert total J -> J m^-2 by dividing by surface (lev=0) area
        if (mode_l == "total") and bool(normalize_by_surface_area):
            A_region = Aeff.isel({lev_dim: 0})
            A_region = xr.where(A_region > 0, A_region, np.nan)
            ohc = (ohc / A_region).rename("OHC_layer")
            intended = "J m^-2"
            ohc.attrs["units"] = intended
            
        # anomaly
        if anomaly:
            bm = str(baseline_mode).lower()
            if (bm == "calendar_year") and hasattr(ohc[time_dim], "dt"):
                yr0 = int(ohc[time_dim].dt.year.min().item())
                base = ohc.where(ohc[time_dim].dt.year == yr0, drop=True).mean(time_dim, skipna=True)
            else:
                n0 = int(max(1, min(int(first_N), ohc.sizes.get(time_dim, ohc.shape[0]))))
                base = ohc.isel({time_dim: slice(0, n0)}).mean(time_dim, skipna=True)
            ohc = (ohc - base).rename("OHC_layer")
            ohc.attrs["units"] = intended
    
        # unit scaling
        ureq = str(to_units).strip().lower()
        if mode_l == "total":
            if ureq in ("zj", "1e22 j", "10^22 j"):
                ohc = ohc / 1e22
                ohc.attrs["units"] = r"10$^22$ J"
        else:
            if ureq in ("gj/m^2", "gj m^-2", "gj per m^2"):
                ohc = ohc / 1e9
                ohc.attrs["units"] = r"GJ m$^-2$"
    
        # normalize time dim name
        if time_dim != "time" and time_dim in ohc.dims:
            ohc = ohc.rename({time_dim: "time"})
            time_dim = "time"
    
        # --------------------------
        # Depth coordinate (axis only)
        # --------------------------
        if self.fdepth and os.path.exists(self.fdepth):
            depth = self._depth_coord_from_fdepth(lev_dim=lev_dim, kind=depth_coord)  # meters

            if self.opt.verbose:
                dvals = np.asarray(depth.values, dtype=float)
                finite = np.isfinite(dvals)

                if finite.any():
                    zmin_ = float(dvals[finite].min())
                    zmax_ = float(dvals[finite].max())
                    mono_ = bool(np.all(np.diff(dvals[finite]) >= 0))
                else:
                    zmin_, zmax_, mono_ = np.nan, np.nan, False

                print(
                    f"[diag][depth] source=fdepth "
                    f"path={os.path.basename(self.fdepth)} "
                    f"kind={depth_coord} "
                    f"lev_dim={lev_dim} "
                    f"nlev={dvals.size} "
                    f"range=({zmin_:.1f},{zmax_:.1f}) m "
                    f"monotonic={mono_}"
                )

        else:
            # fallback: derive depth axis from thickness (no fdepth available)
            H = Hraw.clip(min=0.0)
            Hmean = H.mean(dim=time_dim, skipna=True)

            z_bot = Hmean.cumsum(dim=lev_dim)
            z_top = z_bot.shift({lev_dim: 1}, fill_value=0.0)

            if str(depth_coord).lower() == "bottom":
                depth = z_bot.rename("depth")
            else:
                depth = (0.5 * (z_top + z_bot)).rename("depth")

            depth.attrs.update(units="m", positive="down", long_name="Depth below sea surface")

            if self.opt.verbose:
                dvals = np.asarray(depth.values, dtype=float)
                finite = np.isfinite(dvals)

                if finite.any():
                    zmin_ = float(dvals[finite].min())
                    zmax_ = float(dvals[finite].max())
                    mono_ = bool(np.all(np.diff(dvals[finite]) >= 0))
                else:
                    zmin_, zmax_, mono_ = np.nan, np.nan, False

                print(
                    f"[diag][depth] source=thickness_fallback "
                    f"kind={depth_coord} "
                    f"lev_dim={lev_dim} "
                    f"nlev={dvals.size} "
                    f"range=({zmin_:.1f},{zmax_:.1f}) m "
                    f"monotonic={mono_}"
                )
                
        # --- ROBUST: force depth to become the DIMENSION (not just a coord) ---
        # Ensure depth is 1D over lev_dim
        if depth.dims != (lev_dim,):
            if len(depth.dims) == 1:
                depth = depth.rename({depth.dims[0]: lev_dim})
            else:
                raise ValueError(f"depth must be 1D over '{lev_dim}', got dims={depth.dims}")
    
        # attach depth as coordinate tied to lev_dim
        ohc = ohc.assign_coords(depth=(lev_dim, np.asarray(depth.values, dtype=float)))
    
        # swap lev_dim -> depth
        ohc = ohc.swap_dims({lev_dim: "depth"})
        ohc["depth"].attrs.update(
            units="m",
            positive="down",
            long_name="Depth below sea surface",
        )
        
        if self.opt.verbose:
            print(
                f"[diag][hov] dims={dict(ohc.sizes)} "
                f"depth_dim={'depth' in ohc.dims} "
                f"sorted_depth={bool(np.all(np.diff(np.asarray(ohc['depth'].values)) >= 0))}"
            )
            
        # drop lev_dim coord to avoid confusion with 0..nLevels
        if lev_dim in ohc.coords:
            ohc = ohc.drop_vars(lev_dim)
    
        ohc = ohc.sortby("depth")
    
        ohc = ohc.assign_attrs(
            long_name=f"{'Δ' if anomaly else ''}OHC per layer ({mode_l})",
            region_index=int(region_index),
            rho0=float(rho),
            cp0=float(cp),
            factor_used=fac,
            geometry_source=str(geom_src),
        )

        if not return_cumulative:
            return xr.Dataset(dict(OHC_layer=ohc))
    
        if from_top:
            cum = ohc.cumsum("depth").rename("OHC_cumulative")
        else:
            cum = ohc.sortby(
                "depth",
                ascending=False
            ).cumsum("depth").sortby("depth").rename("OHC_cumulative")
            
        cum.attrs["units"] = ohc.attrs.get("units", "")
    
        return xr.Dataset(dict(OHC_layer=ohc, OHC_cumulative=cum))

    # --------------------------
    # public API
    # --------------------------
    def extract_series(
        self,
        *,
        base_dir: str,
        exp: str,
        sub_dir: str,
        year_min: int,
        year_max: int,
        req: OHCSeriesRequest,
        annual_mean: bool = False,
        annual_weighted: bool = True,
        annual_drop_incomplete: bool = True,
        annual_min_days: int = 360,
    ) -> np.ndarray:
        obj = self.extract_monthly_da(
            base_dir=base_dir, exp=exp, sub_dir=sub_dir,
            year_min=year_min, year_max=year_max, req=req
        )
    
        if isinstance(obj, xr.Dataset):
            da = obj["OHU"] if (req.kind.upper() == "OHU") else obj["OHC"]
        else:
            da = obj
    
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
    
        if not np.isfinite(out).any():
            raise ValueError(f"All-NaN series for exp='{exp}', req={req}")
        return out

    def extract_monthly_da(
        self,
        *,
        base_dir: str,
        exp: str,
        sub_dir: str,
        year_min: int,
        year_max: int,
        req: OHCSeriesRequest,
    ) -> Union[xr.DataArray, xr.Dataset]:
        cache_key = (
            base_dir, exp, sub_dir, year_min, year_max,
            req.kind, req.layer, req.region_index, req.region_dim,
            req.anomaly, req.baseline_months, req.to_ZJ,
            req.rho, req.cp, req.use_mask_in_volume,
            req.smooth_months, req.post_smooth_months, req.area_tol_rel,
            self.file_key, self.fdepth, self.opt.calendar, self.opt.reassign_monthly_time, self.opt.length_policy,
        )
        if cache_key in self._series_cache:
            return self._series_cache[cache_key]

        fns = self._files_for_window(base_dir=base_dir, exp=exp, sub_dir=sub_dir, year_min=year_min, year_max=year_max)
        if not fns:
            raise FileNotFoundError(f"No {self.file_key} files found under {self._resolve_base(base_dir, exp, sub_dir)}")

        ds = self._open_mfdataset_robust(fns, year_min=year_min, year_max=year_max)

        if "time" in ds.dims and hasattr(ds["time"], "dt"):
            yrs = ds["time"].dt.year.astype(int)
            ds = ds.where((yrs >= int(year_min)) & (yrs <= int(year_max)), drop=True)

        ohc = self._compute_ohc_monthly(
            ds,
            region_index=req.region_index,
            region_dim=req.region_dim,
            layer=req.layer,
            rho=req.rho,
            cp=req.cp,
            anomaly=req.anomaly,
            baseline_months=req.baseline_months,
            to_ZJ=req.to_ZJ,
            use_mask_in_volume=bool(getattr(req, "use_mask_in_volume", True)),  # <-- ADD THIS
        )
        
        if req.kind.upper() == "OHC":
            ds.close()
            self._series_cache[cache_key] = ohc
            return ohc

        lev_dim = next((d for d in ("nVertLevels", "z", "nz", "nLevels", "lev") if d in ds.dims), None)
        if lev_dim is None:
            ds.close()
            raise ValueError("Cannot infer lev_dim for OHU computation.")

        out = self._infer_ohu_from_ohc(
            ohc,
            ds_for_area=ds,
            region_index=req.region_index,
            region_dim=req.region_dim,
            lev_dim=lev_dim,
            smooth_months=req.smooth_months,
            post_smooth_months=req.post_smooth_months,
            area_tol_rel=req.area_tol_rel,
            use_mask_in_area=req.use_mask_in_volume,  # keep consistent with OHC volume proxy
        )
        ds.close()
        self._series_cache[cache_key] = out
        return out

    def extract_hovmoller_monthly(
        self,
        *,
        base_dir: str,
        exp: str,
        sub_dir: str,
        year_min: int,
        year_max: int,
        req: OHCHovmollerRequest,
    ) -> Union[xr.DataArray, xr.Dataset]:
        fns = self._files_for_window(base_dir=base_dir, exp=exp, sub_dir=sub_dir, year_min=year_min, year_max=year_max)
        if not fns:
            raise FileNotFoundError(f"No {self.file_key} files found under {self._resolve_base(base_dir, exp, sub_dir)}")

        ds = self._open_mfdataset_robust(fns, year_min=year_min, year_max=year_max)

        if "time" in ds.dims and hasattr(ds["time"], "dt"):
            yrs = ds["time"].dt.year.astype(int)
            ds = ds.where((yrs >= int(year_min)) & (yrs <= int(year_max)), drop=True)

        out = self._compute_ohc_hovmoller_monthly(
            ds,
            region_index=req.region_index,
            region_dim=req.region_dim,
            rho=req.rho,
            cp=req.cp,
            anomaly=req.anomaly,
            baseline_mode=req.baseline_mode,
            first_N=req.first_N,
            mode=req.mode,
            to_units=req.to_units,
            depth_coord=req.depth_coord,
            return_cumulative=req.return_cumulative,
            from_top=req.from_top,
            normalize_by_surface_area=req.normalize_by_surface_area,
            window=req.window,
            apply_window_overlap=req.apply_window_overlap,
            use_mask_in_area=req.use_mask_in_area,
        )

        ds.close()
        return out
        
    def extract_hovmoller_many(
        self,
        *,
        base_dir: str,
        exps: Sequence[str],
        subdirs: Dict[str, str],
        year_min: int,
        year_max: int,
        req: OHCHovmollerRequest,
        annual_mean: bool = False,
        annual_weighted: bool = True,
        annual_drop_incomplete: bool = True,
        annual_min_days: int = 360,
    ) -> Dict[str, Union[xr.DataArray, xr.Dataset]]:
        """
        Batch extract Hovmöller OHC for many experiments.
    
        Returns:
          out[exp] -> DataArray (time, depth) or Dataset with OHC_layer/OHC_cumulative
          If annual_mean=True, annualizes along time for each returned field.
        """
        out: Dict[str, Union[xr.DataArray, xr.Dataset]] = {}
    
        for exp in exps:
            obj = self.extract_hovmoller_monthly(
                base_dir=base_dir,
                exp=exp,
                sub_dir=subdirs[exp],
                year_min=year_min,
                year_max=year_max,
                req=req,
            )
    
            if not annual_mean:
                out[exp] = obj
                continue
    
            # annualize
            if isinstance(obj, xr.Dataset):
                obj2 = obj.copy()
                if "OHC_layer" in obj2:
                    obj2["OHC_layer"] = self.to_annual_mean(
                        obj2["OHC_layer"],
                        weighted=annual_weighted,
                        drop_incomplete=annual_drop_incomplete,
                        min_days=annual_min_days,
                        ymin=year_min,
                        ymax=year_max,
                    )
                if "OHC_cumulative" in obj2:
                    obj2["OHC_cumulative"] = self.to_annual_mean(
                        obj2["OHC_cumulative"],
                        weighted=annual_weighted,
                        drop_incomplete=annual_drop_incomplete,
                        min_days=annual_min_days,
                        ymin=year_min,
                        ymax=year_max,
                    )
                out[exp] = obj2
            else:
                out[exp] = self.to_annual_mean(
                    obj,
                    weighted=annual_weighted,
                    drop_incomplete=annual_drop_incomplete,
                    min_days=annual_min_days,
                    ymin=year_min,
                    ymax=year_max,
                )
    
        return out

    def extract_many(
        self,
        *,
        base_dir: str,
        exps: Sequence[str],
        subdirs: Dict[str, str],
        year_min: int,
        year_max: int,
        requests: Dict[str, OHCSeriesRequest],
        annual_mean: bool = False,
        annual_weighted: bool = True,
        annual_drop_incomplete: bool = True,
        annual_min_days: int = 360,
    ) -> Dict[str, Dict[str, np.ndarray]]:
        out: Dict[str, Dict[str, np.ndarray]] = {}
        for exp in exps:
            sub_dir = subdirs[exp]
            out[exp] = {}
            for alias, req in requests.items():
                out[exp][alias] = self.extract_series(
                    base_dir=base_dir,
                    exp=exp,
                    sub_dir=sub_dir,
                    year_min=year_min,
                    year_max=year_max,
                    req=req,
                    annual_mean=annual_mean,
                    annual_weighted=annual_weighted,
                    annual_drop_incomplete=annual_drop_incomplete,
                    annual_min_days=annual_min_days,
                )
        return out



# ==============================
# Convenience wrappers
# ==============================
def extract_ohc(
    *,
    base_dir: str,
    exps: Sequence[str],
    exp_subdirs: Dict[str, str],
    cfg: OHCExtractConfig,
    requests: Dict[str, OHCSeriesRequest],
    file_key: str = "mpasTimeSeriesOcean",
    fdepth: str,
    engine: Optional[str] = "netcdf4",
    chunks: Optional[dict] = None,
    verbose: bool = False,
) -> Dict[str, Dict[str, np.ndarray]]:
    year_min, year_max = cfg.ts_window

    opt = OHCExtractOptions(
        engine=engine,
        chunks=chunks,
        calendar=cfg.calendar,
        reassign_monthly_time=cfg.reassign_monthly_time,
        length_policy=cfg.length_policy,
        verbose=bool(verbose),
    )

    X = OHCTimeSeriesExtractor(
        options=opt,
        file_key=file_key,
        fdepth=fdepth,
    )

    return X.extract_many(
        base_dir=base_dir,
        exps=exps,
        subdirs=exp_subdirs,
        year_min=year_min,
        year_max=year_max,
        requests=requests,
        annual_mean=cfg.annual_mean,
        annual_weighted=cfg.annual_weighted,
        annual_drop_incomplete=cfg.annual_drop_incomplete,
        annual_min_days=cfg.annual_min_days,
    )

def extract_all_ohc(
    *,
    root_dir: str,
    exps: Sequence[str],
    exp_subdirs: Dict[str, str],
    var_dict: Dict[str, Dict[str, Any]],
    cfg: OHCExtractConfig,
    region_index: int,
    region_dim: str = "nOceanRegions",
    anomaly: bool = True,
    baseline_months: int = 12,
    to_ZJ: bool = True,
    rho: float = 1026.0,
    cp: float = 3996.0,
    use_mask_in_volume: bool = True,   # <-- add this
    file_key: str = "mpasTimeSeriesOcean",
    fdepth: str = "",
    engine: str = "netcdf4",
    chunks: dict | None = None,
    verbose: bool = False,
) -> Dict[str, Dict[str, np.ndarray]]:

    requests: Dict[str, OHCSeriesRequest] = {}
    for key, meta in var_dict.items():
        kind = meta.get("kind", "OHC")
        layer_meta = meta.get("layer", {})
        zmin = float(layer_meta.get("min", 0.0))
        zmax = float(layer_meta.get("max", 700.0))
        requests[key] = OHCSeriesRequest(
            kind=kind,
            layer=(zmin, zmax),
            region_index=int(region_index),
            region_dim=region_dim,
            anomaly=bool(anomaly),
            baseline_months=int(baseline_months),
            to_ZJ=bool(to_ZJ),
            rho=float(rho),
            cp=float(cp),
            use_mask_in_volume=bool(use_mask_in_volume),  # <-- add this
        )

    return extract_ohc(
        base_dir=root_dir,
        exps=exps,
        exp_subdirs=exp_subdirs,
        cfg=cfg,
        requests=requests,
        file_key=file_key,
        fdepth=fdepth,
        engine=engine,
        chunks=chunks,
        verbose=verbose,
    )

def extract_all_ohc_hovmoller(
    *,
    root_dir: str,
    exps: Sequence[str],
    exp_subdirs: Dict[str, str],
    cfg: OHCExtractConfig,
    req: OHCHovmollerRequest,
    file_key: str = "mpasTimeSeriesOcean",
    fdepth: str = "",
    engine: str = "netcdf4",
    chunks: dict | None = None,
    verbose: bool = False,
) -> Dict[str, Union[xr.DataArray, xr.Dataset]]:
    year_min, year_max = cfg.ts_window

    opt = OHCExtractOptions(
        engine=engine,
        chunks=chunks,
        calendar=cfg.calendar,
        reassign_monthly_time=cfg.reassign_monthly_time,
        length_policy=cfg.length_policy,
        verbose=bool(verbose),
    )
    X = OHCTimeSeriesExtractor(options=opt, file_key=file_key, fdepth=fdepth)

    return X.extract_hovmoller_many(
        base_dir=root_dir,
        exps=exps,
        subdirs=exp_subdirs,
        year_min=year_min,
        year_max=year_max,
        req=req,
        annual_mean=cfg.annual_mean,
        annual_weighted=cfg.annual_weighted,
        annual_drop_incomplete=cfg.annual_drop_incomplete,
        annual_min_days=cfg.annual_min_days,
    )
    
# ==============================
# Registry
# ==============================
class OHCRegistry:
    """
    FluxRegistry-style convenience wrapper:
      - var_dict: canonical OHC keys -> {kind, layer{min,max}}
      - region_dict: region_key -> region_index
      - build_extract_config(): produces OHCExtractConfig
    """

    def __init__(self):
        self._var_dict: Dict[str, Dict[str, Any]] = self.build_default_var_dict()
        # IMPORTANT: fill with YOUR MPAS region indices
        self._region_dict = {
            "Arctic": {
                "name": "arctic",
                "index": 0,
            },
            "Equator": {
                "name": "equatorial",
                "index": 1,
            },
            "Southern Ocean": {
                "name": "so",
                "index": 2,
            },
            "Nino3": {
                "name": "nino3",
                "index": 3,
            },
            "Nino4": {
                "name": "nino4",
                "index": 4,
            },
            "Nino3.4": {
                "name": "nino3.4",
                "index": 5,
            },
            "Global": {
                "name": "global",
                "index": 6,
            },
        }

    def get_var_dict(self) -> Dict[str, Dict[str, Any]]:
        return self._var_dict

    def get_region_index(self, region_key: str) -> int:
        if region_key not in self._region_dict:
            raise KeyError(
                f"Unknown region_key='{region_key}'. "
                f"Available: {list(self._region_dict.keys())}"
            )
        return int(self._region_dict[region_key]["index"])

    @staticmethod
    def build_default_var_dict() -> Dict[str, Dict[str, Any]]:
        return {
            # --- OHC ---
            "OHC_0_700":       {"kind": "OHC", "layer": {"min": 0.0,    "max": 700.0,  "unit": "m"}},
            "OHC_700_2000":    {"kind": "OHC", "layer": {"min": 700.0,  "max": 2000.0, "unit": "m"}},
            "OHC_2000_bottom": {"kind": "OHC", "layer": {"min": 2000.0, "max": 6000.0, "unit": "m"}},
            "OHC_0_bottom":    {"kind": "OHC", "layer": {"min": 0.0,    "max": 6000.0, "unit": "m"}},
            # --- OHU (derived from OHC tendency / area) ---
            # Units are W m-2, depth window applies to the underlying OHC you differentiate.
            "OHU_0_700":       {"kind": "OHU", "layer": {"min": 0.0,    "max": 700.0,  "unit": "m"}},
            "OHU_700_2000":    {"kind": "OHU", "layer": {"min": 700.0,  "max": 2000.0, "unit": "m"}},
            "OHU_2000_bottom": {"kind": "OHU", "layer": {"min": 2000.0, "max": 6000.0, "unit": "m"}},
            "OHU_0_bottom":    {"kind": "OHU", "layer": {"min": 0.0,    "max": 6000.0, "unit": "m"}},
        }
    
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
    ) -> OHCExtractConfig:
        return OHCExtractConfig(
            ts_window=ts_window,
            annual_mean=annual_mean,
            annual_weighted=annual_weighted,
            annual_drop_incomplete=annual_drop_incomplete,
            annual_min_days=annual_min_days,
            calendar=calendar,
            reassign_monthly_time=reassign_monthly_time,
            length_policy=length_policy,
        )