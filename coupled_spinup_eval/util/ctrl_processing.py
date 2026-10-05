"""Processing for the ctrl analyses (one experiment: spin-up followed by its piControl):
annual global time series written to data/ caches.

  MPASDiagnosticsBuilder (shared base) and one subclass per analysis:
    SSTBuilder, SSSBuilder, SSHBuilder, SeaIceBuilder  (14_ocn_ts_ctrl_analysis)
    OHCBuilder                                        (15_ohc_ts_ctrl_analysis)
    OceanBudgetBuilder                                (13_ocn_conserve_ctrl_analysis)
  AMOCDiagnosticsBuilder                              (14_ocn_ts_ctrl_analysis)
  AtmosphereBalanceAnalysis                           (12_atm_flux_ctrl_analysis)
"""
import os, re, glob, warnings, string
import numpy as np
import xarray as xr
import cftime
import contextlib
from scipy.stats import linregress
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter, FixedLocator
import pandas as pd
try:
    from cftime import num2date
except Exception:
    num2date = None
import os
import glob
import string
import warnings


class MPASDiagnosticsBuilder:
    """Shared reader/annualizer for MPAS-O / MPAS-SI global time series (spin-up + piControl). Use a subclass."""

    def __init__(
        self,
        base_dir,
        spinup_name,
        picontrol_name,
        subdir="post/time_series",
        file_stem="mpasTimeSeriesOcean",          # NEW: stem for discovery (prefix of files)
        stamp_regex=None,                         # NEW: how to parse yyyymm[-dd]-yyyymm[-dd]
        file_suffix=".nc",
        xr_engine="h5netcdf",
        xr_chunks=None,
        var_name="timeMonthly_avg_avgValueWithinOceanRegion_avgSurfaceTemperature",
        default_units=None,                        # e.g., {"timeMonthly_avg_*SurfaceTemperature": "degC"}
        force_reprocess=False,
        region=None,     # e.g., 6 for "Global Ocean"
        region_dim_candidates=("nOceanRegions", "nOceanRegionsTmp", "nRegions", "region", "oceanRegions"),
        level_index=None,      # e.g., 0 for surface/first layer
        level_dim_candidates=("nVertLevels", "nVertLevelsP1", "nLevels", "z", "depth"),
    ):
        self.base_dir = base_dir
        self.spinup_name = spinup_name
        self.picontrol_name = picontrol_name
        self.subdir = subdir

        self.file_stem = file_stem
        self.file_suffix = file_suffix
        self._STAMP_RE = re.compile(stamp_regex) if stamp_regex else None

        self.xr_engine = xr_engine
        self.xr_chunks = xr_chunks
        self.var_name = var_name

        self.default_units = default_units or {}
        self.force_reprocess = force_reprocess
        
        self.ocn_region_map = {
            0: {"name": "Arctic", "short": "arctic"}, 
            1: {"name": "Equatorial (15S-15N)", "short": "equatorial"}, 
            2: {"name": "Southern Ocean", "short": "so"}, 
            3: {"name": "Nino 3",   "short": "nino3"}, 
            4: {"name": "Nino 4",   "short": "nino4"}, 
            5: {"name": "Nino 3.4", "short": "nino3.4"}, 
            6: {"name": "Global Ocean", "short": "global"}, 
        }
        
        # Resolve region_short → region_index if provided
        if region is not None:
            region_index = self.region_index_from_short(region.lower())
            
        self.region_index = region_index
        self.region_dim_candidates = tuple(region_dim_candidates)
        self.level_index = level_index
        self.level_dim_candidates = tuple(level_dim_candidates)

    def region_index_from_short(self, short):
        short = str(short).lower()
        for idx, info in self.ocn_region_map.items():
            if info["short"].lower() == short:
                return idx
        raise ValueError(
            f"Unknown ocean region short name '{short}'. "
            f"Valid options: {[v['short'] for v in self.ocn_region_map.values()]}"
        )

    def _is_fresh(self, target, sources):
        if not os.path.exists(target):
            return False
        t_mtime = os.path.getmtime(target)
        for s in sources:
            try:
                if os.path.getmtime(s) > t_mtime:
                    return False
            except FileNotFoundError:
                return False
        return True

    def _is_fresh_multi(self, target, sources):
        if not os.path.exists(target):
            return False
        t_mtime = os.path.getmtime(target)
        for s in sources:
            if (not os.path.exists(s)) or (os.path.getmtime(s) > t_mtime):
                return False
        return True

    def _extract_stamps(self, fname):
        if not self._STAMP_RE:
            return None, None
        m = self._STAMP_RE.search(os.path.basename(fname))
        if not m:
            return None, None
        return int(m.group(1)), int(m.group(2))

    def _find_ts_files(self, case_root):
        # Try exact file first (common when there is a single merged file)
        exact = os.path.join(case_root, self.subdir, f"{self.file_stem}{self.file_suffix}")
        if os.path.isfile(exact):
            return [exact]

        # Otherwise, glob for any matching files
        pattern = os.path.join(case_root, self.subdir, f"{self.file_stem}*{self.file_suffix}")
        files = glob.glob(pattern)
        if not files:
            pattern = os.path.join(case_root, self.subdir, "**", f"{self.file_stem}*{self.file_suffix}")
            files = glob.glob(pattern, recursive=True)

        if not files:
            return []

        # If a stamp regex is provided, sort by the start stamp; otherwise keep lexicographic order
        if self._STAMP_RE:
            with_stamps, no_stamps = [], []
            for f in files:
                s0, s1 = self._extract_stamps(f)
                (with_stamps if s0 is not None else no_stamps).append((f, s0))
            with_stamps.sort(key=lambda t: t[1])
            return [f for f, _ in with_stamps] + [f for f in no_stamps]
        else:
            # No stamps — just sort path names so ordering is deterministic
            files.sort()
            return files

    def _get_time_name(self, ds):
        if "time" in ds.dims:
            return "time"
        if "Time" in ds.dims:
            return "Time"

        def looks_like_time(n):
            n = n.lower()
            return ("time" in n) or ("xtime" in n)

        for name in list(ds.coords) + list(ds.variables):
            da = ds[name]
            if da.ndim == 1 and looks_like_time(name):
                return da.dims[0]
        for name in list(ds.coords) + list(ds.variables):
            da = ds[name]
            if da.ndim == 1:
                return da.dims[0]
        raise KeyError("No recognizable time dimension in dataset.")

    def _compute_year_coord(self, ds, tname):
        v = ds[tname]

        # char arrays / string time (e.g., xtime; or "YYYY-MM-DD_*")
        if (v.dtype.kind in ("S", "U")) or ("str" in str(v.dtype).lower()):
            arr = v.values
            if arr.ndim == 2:
                strings = []
                for row in arr:
                    s = "".join([(c.decode("utf-8") if isinstance(c, (bytes, np.bytes_)) else str(c)) for c in row]).strip()
                    strings.append(s)
            else:
                strings = [(s.decode("utf-8") if isinstance(s, (bytes, np.bytes_)) else str(s)) for s in v.values]
            years = [int(s[:4]) if len(s) >= 4 and s[:4].isdigit() else np.nan for s in strings]
            return np.array(years, dtype=float)

        units = (v.attrs.get("units") or "").lower()
        if "months since" in units:
            m = re.search(r"months since\s+(\d{1,4})", units)
            base_year = int(m.group(1)) if m else 1
            months = np.asarray(v.values, dtype=float)
            return np.floor(base_year + months / 12.0).astype(int)
        if "days since" in units:
            m = re.search(r"days since\s+(\d{1,4})", units)
            base_year = int(m.group(1)) if m else 1
            days = np.asarray(v.values, dtype=float)
            return np.floor(base_year + days / 365.0).astype(int)
        if np.issubdtype(v.dtype, np.datetime64):
            return v.dt.year.values

        # fallback: simple index (rare)
        n = v.sizes[v.dims[0]]
        return np.arange(n)

    def _aggregate_annual_mean_stream(self, files):
        if not files:
            raise FileNotFoundError("No MPAS time series files found.")

        year_all, val_all = [], []
        vname = self.var_name
        chosen_units = None

        for f in files:
            ds = xr.open_dataset(f, engine=self.xr_engine, decode_times=False, chunks=self.xr_chunks)
            try:
                tname = self._get_time_name(ds)
                years = self._compute_year_coord(ds, tname)
                da = ds[vname]
                
                # --- Select region (e.g., Global Ocean = 6) if variable has a region dim ---                
                if self.region_index is not None:
                    print(self.region_dim_candidates)
                    print(da.dims)
                    
                    for rd in self.region_dim_candidates:
                        if rd in da.dims:
                            da = da.isel({rd: int(self.region_index)})
                            break
                
                # --- Select vertical level (e.g., first layer = 0) if variable has a level dim ---
                if self.level_index is not None:
                    for zd in self.level_dim_candidates:
                        if zd in da.dims:
                            da = da.isel({zd: int(self.level_index)})
                            break
            
                # try to remember/propagate units (fallback to provided default_units)
                if chosen_units is None:
                    chosen_units = da.attrs.get("units") or self.default_units.get(vname)

                vals = da.values
                if vals.ndim != 1:
                    # move time to front, then mean over the rest
                    time_ax = da.get_axis_num(tname) if tname in da.dims else 0
                    if time_ax != 0:
                        vals = np.moveaxis(vals, time_ax, 0)
                    while vals.ndim > 1:
                        vals = np.nanmean(vals, axis=-1)

                year_all.append(years.astype(np.int64, copy=False))
                val_all.append(np.asarray(vals, dtype=float))
            finally:
                ds.close()

        years = np.concatenate(year_all)
        vals  = np.concatenate(val_all)
        good = np.isfinite(vals) & np.isfinite(years)
        years = years[good]
        vals  = vals[good]
        if years.size == 0:
            raise ValueError("All annual values are NaN after screening.")

        y0, y1 = int(np.nanmin(years)), int(np.nanmax(years))
        bins = y1 - y0 + 1
        idx = years - y0
        ssum = np.bincount(idx, weights=vals, minlength=bins)
        scount = np.bincount(idx, minlength=bins)
        with np.errstate(invalid="ignore", divide="ignore"):
            annual = ssum / scount

        year_coord = np.arange(y0, y1 + 1, dtype=int)
        ds_out = xr.Dataset(
            data_vars={vname: (("year",), annual, {"units": chosen_units or self.default_units.get(vname, "")})},
            coords={"year": ("year", year_coord)},
            attrs={},
        )
        return ds_out

    def build_sources(self, out_dir, suffix="annual", var_name="diag"):
        """
        Returns dict with paths to annual NetCDFs for spinup & piControl.

        Output files:
          <spinup_name>_<var_name>_timeseries_<suffix>.nc
          <picontrol_name>_<var_name>_timeseries_<suffix>.nc
        """
        os.makedirs(out_dir, exist_ok=True)
        spinup_root = os.path.join(self.base_dir, self.spinup_name)
        pctl_root   = os.path.join(self.base_dir, self.picontrol_name)

        spinup_files = self._find_ts_files(spinup_root)
        pctl_files   = self._find_ts_files(pctl_root)

        if not spinup_files:
            raise FileNotFoundError(f"No MPAS files found under {spinup_root}/{self.subdir} with stem '{self.file_stem}'")
        if not pctl_files:
            raise FileNotFoundError(f"No MPAS files found under {pctl_root}/{self.subdir} with stem '{self.file_stem}'")

        spinup_nc = os.path.join(out_dir, f"{self.spinup_name}_{var_name}_timeseries_{suffix}.nc")
        pctl_nc   = os.path.join(out_dir, f"{self.picontrol_name}_{var_name}_timeseries_{suffix}.nc")

        # Respect force_reprocess
        if self.force_reprocess:
            for path in (spinup_nc, pctl_nc):
                try:
                    os.remove(path)
                except FileNotFoundError:
                    pass

        if not self._is_fresh(spinup_nc, spinup_files):
            ds_s = self._aggregate_annual_mean_stream(spinup_files)
            enc = {self.var_name: {"zlib": True, "complevel": 4}}
            ds_s.to_netcdf(spinup_nc, engine=self.xr_engine, encoding=enc)

        if not self._is_fresh(pctl_nc, pctl_files):
            ds_p = self._aggregate_annual_mean_stream(pctl_files)
            enc = {self.var_name: {"zlib": True, "complevel": 4}}
            ds_p.to_netcdf(pctl_nc, engine=self.xr_engine, encoding=enc)

        return {"spinup": spinup_nc, "piControl": pctl_nc}

    def _centered_filter(self, y, window=11):
        w = max(1, int(window))
        if w % 2 == 0:
            w += 1
        pad = w // 2
        y = np.asarray(y, dtype=float)
        isfin = np.isfinite(y).astype(float)
        y_filled = np.where(np.isfinite(y), y, 0.0)
        k = np.ones(w, dtype=float)
        num = np.convolve(y_filled, k, mode="same")
        den = np.convolve(isfin, k, mode="same")
        out = np.full_like(y, np.nan, dtype=float)
        good = den > 0
        out[good] = num[good] / den[good]
        out[:pad] = np.nan
        out[-pad:] = np.nan
        return out

    def _sliding_trend(self, y, x, window):
        n = len(x)
        w = max(2, int(window))
        if w % 2 == 0:
            w += 1
        half = w // 2
        out = np.full(n, np.nan, dtype=float)
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        for i in range(n):
            i0 = max(0, i - half)
            i1 = min(n, i + half + 1)
            xx = x[i0:i1]
            yy = y[i0:i1]
            m = np.isfinite(xx) & np.isfinite(yy)
            if m.sum() >= 3:
                coef = np.polyfit(xx[m], yy[m], 1)  # per-year
                out[i] = coef[0]
        return out

    def _load_year_val(self, ncpath):
        """Load (years, values) for self.var_name from a per-source annual file."""
        vname = self.var_name
        ds = xr.open_dataset(ncpath, decode_times=False, engine=self.xr_engine)
        try:
            if "year" in ds.coords:
                years = ds["year"].values
                vals  = ds[vname].values
            else:
                tname = self._get_time_name(ds)
                years = self._compute_year_coord(ds, tname)
                vals  = ds[vname].values
            vals  = np.asarray(vals, dtype=float)
            years = np.asarray(years, dtype=int)
            n = min(len(years), len(vals))
            return years[:n], vals[:n]
        finally:
            ds.close()

    def make_years(self, sequence, y_spin, y_pi):
        if sequence == "calendar":
            return np.concatenate([y_spin, y_pi])

        if sequence == "append":
            # Always continue sequentially after spin-up (1..N + N+1..)
            offset = int(np.nanmax(y_spin)) if len(y_spin) else 0
            y_pi_seq = (np.arange(len(y_pi)) + offset + 1).astype(int)
            return np.concatenate([y_spin.astype(int, copy=False), y_pi_seq])

        # "auto": try to detect overlap/reset and choose a sensible behavior
        y_spin = np.asarray(y_spin, dtype=int)
        y_pi   = np.asarray(y_pi, dtype=int)
        overlap = (y_pi.min() <= y_spin.max())
        reset_like = (np.median(y_pi) <= np.median(y_spin)) or (y_pi.min() <= 5)
        continues = (y_pi.min() > y_spin.max()) and (y_pi.min() - y_spin.max() <= 5)

        if overlap or reset_like:
            offset = int(np.nanmax(y_spin)) if len(y_spin) else 0
            y_pi_seq = (np.arange(len(y_pi)) + offset + 1).astype(int)
            return np.concatenate([y_spin, y_pi_seq])
        elif continues:
            return np.concatenate([y_spin, y_pi])
        else:
            offset = int(np.nanmax(y_spin)) if len(y_spin) else 0
            y_pi_seq = (np.arange(len(y_pi)) + offset + 1).astype(int)
            return np.concatenate([y_spin, y_pi_seq])

    def build_combined_tidy(
        self,
        spinup_nc,
        pctl_nc,
        combined_tidy_nc,
        filter_window=11,
        trend_windows=(10,),
        add_per_decade=True,
        sequence="auto",   # "auto" | "calendar" | "append"
        force=False,
        strict=True,
    ):
        deps = [spinup_nc, pctl_nc]
        if (not force) and self._is_fresh_multi(combined_tidy_nc, deps):
            return combined_tidy_nc

        vname = self.var_name

        # Load sources
        y_s, v_s = self._load_year_val(spinup_nc)
        y_p, v_p = self._load_year_val(pctl_nc)
        n_s, n_p = len(y_s), len(y_p)

        # Construct combined 'year'
        years_combined = self.make_years(sequence, y_s, y_p).astype(int, copy=False)
        values_combined = np.concatenate([v_s, v_p])

        # Drop immediate duplicates (defensive)
        if years_combined.size > 1:
            keep = np.ones_like(years_combined, dtype=bool)
            for i in range(len(years_combined) - 1):
                if years_combined[i] == years_combined[i + 1]:
                    keep[i] = False
            if not np.all(keep):
                years_combined  = years_combined[keep]
                values_combined = values_combined[keep]

        n_comb = len(years_combined)
        expected = n_s + n_p
        if strict and n_comb != expected:
            raise ValueError(
                f"[combine] Expected {expected} samples (spinup={n_s} + piControl={n_p}) "
                f"but got {n_comb}. sequence='{sequence}'. "
                f"Spinup years {int(y_s.min())}-{int(y_s.max()) if len(y_s) else -1}, "
                f"piControl years {int(y_p.min())}-{int(y_p.max()) if len(y_p) else -1}."
            )

        # Coords & data
        model_year = np.arange(n_comb, dtype=int)
        rgn = np.zeros(n_comb, dtype=int)

        # Units propagation
        units = self.default_units.get(vname, "")
        # Try reading units from one of the per-source files (cheap)
        try:
            ds_tmp = xr.open_dataset(spinup_nc, engine=self.xr_engine, decode_times=False)
            units = ds_tmp[vname].attrs.get("units", units)
            ds_tmp.close()
        except Exception:
            pass

        data_vars = {vname: (("time",), values_combined, {"units": units})}

        if filter_window and filter_window > 1:
            data_vars[f"filtered_{int(filter_window)}yr"] = (
                ("time",),
                self._centered_filter(values_combined, filter_window)
            )

        for W in (trend_windows or ()):
            W = int(W)
            tname = f"sliding_trend_{W}yr"
            trend = self._sliding_trend(values_combined, years_combined, W)
            data_vars[tname] = (("time",), trend)
            if add_per_decade:
                data_vars[f"{tname}_per_decade"] = (("time",), trend * 10.0)

        ds_out = xr.Dataset(
            data_vars=data_vars,
            coords={
                "time": ("time", np.arange(n_comb, dtype=int)),
                "year": ("time", years_combined),
                "model_year": ("time", model_year),
                "rgn": ("time", rgn),
            },
            attrs={
                "series": f"{self.spinup_name} + {self.picontrol_name} (sequence={sequence})",
                "filter_window": int(filter_window) if filter_window else 0,
                "trend_windows": ",".join(str(int(w)) for w in (trend_windows or [])),
                "units": units,
                "spinup_len": int(n_s),
                "picontrol_len": int(n_p),
                "combined_len": int(n_comb),
                "spinup_year_min": int(np.nanmin(y_s)) if n_s else -1,
                "spinup_year_max": int(np.nanmax(y_s)) if n_s else -1,
                "picontrol_year_min": int(np.nanmin(y_p)) if n_p else -1,
                "picontrol_year_max": int(np.nanmax(y_p)) if n_p else -1,
                "combined_year_min": int(np.nanmin(years_combined)) if n_comb else -1,
                "combined_year_max": int(np.nanmax(years_combined)) if n_comb else -1,
            },
        )
        enc = {k: {"zlib": True, "complevel": 4} for k in ds_out.data_vars}
        ds_out.to_netcdf(combined_tidy_nc, engine=self.xr_engine, encoding=enc)
        return combined_tidy_nc

    def _derive_ssh(self, ds):
        factor = 100.0  # m -> cm, if inputs are m^3 and m^2 this yields m, then *100 => cm
        v_area = "timeMonthly_avg_areaCellGlobal"
        v_vol  = "timeMonthly_avg_volumeCellGlobal"

        for vname in (v_area, v_vol):
            if vname not in ds:
                raise KeyError(f"Variable '{vname}' not found in dataset.")

        tdim = self._get_time_name(ds)

        area = ds[v_area].astype("float64")
        vol  = ds[v_vol ].astype("float64")
        ssh  = (vol / area) * factor  # cm

        if not np.isfinite(ssh.values).all():
            raise ValueError("Non-finite values in derived SSH proxy.")

        ssh.name = "SSH"
        ssh.attrs.update({"long_name": "Global mean sea surface height (proxy)", "units": "cm"})
        # keep dims with time first
        if ssh.dims[0] != tdim:
            ssh = ssh.transpose(tdim, ...)
        return ssh


class SSTBuilder(MPASDiagnosticsBuilder):
    """Global-mean SST annual series."""


class SSSBuilder(MPASDiagnosticsBuilder):
    """Global-mean SSS annual series."""


class SSHBuilder(MPASDiagnosticsBuilder):
    """Global-mean SSH annual series."""

    def __init__(
        self,
        base_dir,
        spinup_name,
        picontrol_name,
        subdir="post/time_series",
        file_stem="mpasTimeSeriesOcean",          # NEW: stem for discovery (prefix of files)
        stamp_regex=None,                         # NEW: how to parse yyyymm[-dd]-yyyymm[-dd]
        file_suffix=".nc",
        xr_engine="h5netcdf",
        xr_chunks=None,
        var_name="timeMonthly_avg_avgValueWithinOceanRegion_avgSurfaceTemperature",
        default_units=None,                        # e.g., {"timeMonthly_avg_*SurfaceTemperature": "degC"}
        force_reprocess=False,
    ):
        self.base_dir = base_dir
        self.spinup_name = spinup_name
        self.picontrol_name = picontrol_name
        self.subdir = subdir

        self.file_stem = file_stem
        self.file_suffix = file_suffix
        self._STAMP_RE = re.compile(stamp_regex) if stamp_regex else None

        self.xr_engine = xr_engine
        self.xr_chunks = xr_chunks
        self.var_name = var_name

        self.default_units = default_units or {}
        self.force_reprocess = force_reprocess

    def _compute_year_coord(self, ds, tname):
        v = ds[tname]

        # 1) CFTime / numpy datetime64 (most robust): use xarray's datetime accessor
        try:
            years = v.dt.year.values  # works for CFTime and numpy datetime64
            return years.astype(int, copy=False)
        except Exception:
            pass

        # 2) char arrays / string time (e.g., xtime* with YYYY-MM-DD..)
        if (v.dtype.kind in ("S", "U")) or ("str" in str(v.dtype).lower()):
            arr = v.values
            if arr.ndim == 2:
                strings = [
                    "".join([(c.decode("utf-8") if isinstance(c, (bytes, np.bytes_)) else str(c)) for c in row]).strip()
                    for row in arr
                ]
            else:
                strings = [(s.decode("utf-8") if isinstance(s, (bytes, np.bytes_)) else str(s)) for s in arr]
            years = [int(s[:4]) if len(s) >= 4 and s[:4].isdigit() else np.nan for s in strings]
            return np.asarray(years, dtype=float)

        # 3) numeric time with units
        units = (v.attrs.get("units") or "").lower()
        if "months since" in units:
            m = re.search(r"months since\s+(\d{1,4})", units)
            base_year = int(m.group(1)) if m else 1
            months = np.asarray(v.values, dtype=float)
            return np.floor(base_year + months / 12.0).astype(int)
        if "days since" in units:
            m = re.search(r"days since\s+(\d{1,4})", units)
            base_year = int(m.group(1)) if m else 1
            days = np.asarray(v.values, dtype=float)
            return np.floor(base_year + days / 365.0).astype(int)

        # 4) already numeric years
        try:
            return np.asarray(v.values, dtype=int)
        except Exception:
            # fallback: simple index (rare)
            n = v.sizes[v.dims[0]]
            return np.arange(n, dtype=int)

    def _aggregate_annual_mean_stream(self, files):
        if not files:
            raise FileNotFoundError("No MPAS time series files found.")

        year_all, val_all = [], []
        chosen_units = None
        vname = self.var_name
        
        for f in files:
            ds = xr.open_dataset(f, engine=self.xr_engine, decode_times=True, chunks=self.xr_chunks)
            try:
                if vname == "SSH":
                    da = self._derive_ssh(ds)  # expects to create a DataArray named vname
                elif vname not in ds:
                    raise KeyError(f"Variable '{vname}' not found in {f}")
                else:
                    da = ds[vname]
                    
                tname = self._get_time_name(ds)
                years = self._compute_year_coord(ds, tname)
                
                # Remember/propagate units once (fallback to default_units)
                if chosen_units is None:
                    chosen_units = da.attrs.get("units") or self.default_units.get(vname)

                # Reduce over all non-time dims safely (works with dask; keeps attrs)
                non_time_dims = [d for d in da.dims if d != tname]
                if non_time_dims:
                    da = da.mean(dim=non_time_dims, skipna=True, keep_attrs=True)

                # Align time axis name just in case
                if tname not in da.dims and "time" in ds.dims:
                    # if reduction removed tname, something's off
                    raise ValueError(f"After reduction, '{vname}' lost its time dimension in {f}")

                # Materialize arrays at the last moment
                vals = np.asarray(da.values, dtype=float)

                # Basic sanity
                if vals.ndim != 1:
                    raise ValueError(f"Expected 1D time series for '{vname}' from {f}, got shape {vals.shape}")

                year_all.append(np.asarray(years, dtype=np.int64))
                val_all.append(vals)

            finally:
                ds.close()
                
        # Concatenate across files, sort by year, and aggregate duplicates (mean)
        if not year_all:
            raise FileNotFoundError(f"No usable chunks for variable '{vname}' in: {files}")

        # --- Concatenate all monthly points ---
        years = np.concatenate(year_all)
        vals  = np.concatenate(val_all)
        good = np.isfinite(vals) & np.isfinite(years)
        years, vals = years[good], vals[good]

        if years.size == 0:
            raise ValueError("All values are NaN after screening.")

        # --- Round or floor to integer years to ensure grouping works ---
        years_int = np.floor(years).astype(int)

        # --- Group by integer year (true annual mean) ---
        df = pd.DataFrame({vname: vals, "year": years_int})
        df_annual = df.groupby("year", as_index=False)[vname].mean(numeric_only=True)
        
        # --- Construct final dataset ---
        ds_out = xr.Dataset(
            data_vars={
                vname: (
                    ("year",),
                    df_annual[vname].to_numpy(dtype=float),
                    {"units": chosen_units or self.default_units.get(vname, "")},
                )
            },
            coords={"year": ("year", df_annual["year"].to_numpy(dtype=int))},
            attrs={
                "aggregation": "annual mean (12 months averaged)",
                "n_source_files": len(files),
            },
        )
        return ds_out

    def build_combined_tidy(
        self,
        spinup_nc,
        pctl_nc,
        combined_tidy_nc,
        filter_window=11,
        trend_windows=(10,),
        add_per_decade=True,
        sequence="auto",   # "auto" | "calendar" | "append"
        force=False,
        strict=True,
    ):
        deps = [spinup_nc, pctl_nc]
        if (not force) and self._is_fresh_multi(combined_tidy_nc, deps):
            return combined_tidy_nc

        vname = self.var_name

        # Load sources
        y_s, v_s = self._load_year_val(spinup_nc)
        y_p, v_p = self._load_year_val(pctl_nc)
        n_s, n_p = len(y_s), len(y_p)

        # Construct combined 'year'
        years_combined = self.make_years(sequence, y_s, y_p).astype(int, copy=False)
        values_all = np.concatenate([v_s, v_p])
        
        # --- Derive anomaly relative to first valid year ---
        values_all = np.asarray(values_all, dtype=float)
        if np.any(np.isfinite(values_all)):
            first_idx = int(np.argmax(np.isfinite(values_all)))  # first finite
            baseline = values_all[first_idx]
            values_combined = values_all - baseline
        else:
            values_combined = np.full_like(values_all, np.nan)

        # Drop immediate duplicates (defensive)
        if years_combined.size > 1:
            keep = np.ones_like(years_combined, dtype=bool)
            for i in range(len(years_combined) - 1):
                if years_combined[i] == years_combined[i + 1]:
                    keep[i] = False
            if not np.all(keep):
                years_combined  = years_combined[keep]
                values_combined = values_combined[keep]

        n_comb = len(years_combined)
        expected = n_s + n_p
        if strict and n_comb != expected:
            raise ValueError(
                f"[combine] Expected {expected} samples (spinup={n_s} + piControl={n_p}) "
                f"but got {n_comb}. sequence='{sequence}'. "
                f"Spinup years {int(y_s.min())}-{int(y_s.max()) if len(y_s) else -1}, "
                f"piControl years {int(y_p.min())}-{int(y_p.max()) if len(y_p) else -1}."
            )

        # Coords & data
        model_year = np.arange(n_comb, dtype=int)
        rgn = np.zeros(n_comb, dtype=int)

        # Units propagation
        units = self.default_units.get(vname, "")
        # Try reading units from one of the per-source files (cheap)
        try:
            ds_tmp = xr.open_dataset(spinup_nc, engine=self.xr_engine, decode_times=False)
            units = ds_tmp[vname].attrs.get("units", units)
            ds_tmp.close()
        except Exception:
            pass

        data_vars = {vname: (("time",), values_combined, {"units": units})}

        if filter_window and filter_window > 1:
            data_vars[f"filtered_{int(filter_window)}yr"] = (
                ("time",),
                self._centered_filter(values_combined, filter_window)
            )

        for W in (trend_windows or ()):
            W = int(W)
            tname = f"sliding_trend_{W}yr"
            trend = self._sliding_trend(values_combined, years_combined, W)
            data_vars[tname] = (("time",), trend)
            if add_per_decade:
                data_vars[f"{tname}_per_decade"] = (("time",), trend * 10.0)

        ds_out = xr.Dataset(
            data_vars=data_vars,
            coords={
                "time": ("time", np.arange(n_comb, dtype=int)),
                "year": ("time", years_combined),
                "model_year": ("time", model_year),
                "rgn": ("time", rgn),
            },
            attrs={
                "series": f"{self.spinup_name} + {self.picontrol_name} (sequence={sequence})",
                "filter_window": int(filter_window) if filter_window else 0,
                "trend_windows": ",".join(str(int(w)) for w in (trend_windows or [])),
                "units": units,
                "spinup_len": int(n_s),
                "picontrol_len": int(n_p),
                "combined_len": int(n_comb),
                "spinup_year_min": int(np.nanmin(y_s)) if n_s else -1,
                "spinup_year_max": int(np.nanmax(y_s)) if n_s else -1,
                "picontrol_year_min": int(np.nanmin(y_p)) if n_p else -1,
                "picontrol_year_max": int(np.nanmax(y_p)) if n_p else -1,
                "combined_year_min": int(np.nanmin(years_combined)) if n_comb else -1,
                "combined_year_max": int(np.nanmax(years_combined)) if n_comb else -1,
            },
        )   
        enc = {k: {"zlib": True, "complevel": 4} for k in ds_out.data_vars}
        ds_out.to_netcdf(combined_tidy_nc, engine=self.xr_engine, encoding=enc)
        return combined_tidy_nc


class SeaIceBuilder(MPASDiagnosticsBuilder):
    """NH/SH sea-ice volume annual series."""

    def __init__(
        self,
        base_dir,
        spinup_name,
        picontrol_name,
        subdir="post/time_series",
        file_stem="mpasTimeSeriesOcean",
        stamp_regex=None,
        file_suffix=".nc",
        xr_engine="netcdf4",
        xr_chunks=None,
        var_name="iceVolume",
        var_unit=None,
        var_fact=None,
        default_units=None,
        force_reprocess=False,
        year_mode="calendar",  # "calendar" | "zero_based" (applies to per-source output)
    ):
        self.base_dir = base_dir
        self.spinup_name = spinup_name
        self.picontrol_name = picontrol_name
        self.subdir = subdir

        self.file_stem = file_stem
        self.file_suffix = file_suffix
        self._STAMP_RE = re.compile(stamp_regex) if stamp_regex else None

        self.xr_engine = xr_engine
        self.xr_chunks = xr_chunks
        self.var_name = var_name
        self.var_unit = var_unit
        self.var_fact = var_fact

        self.default_units = default_units or {}
        self.force_reprocess = force_reprocess
        self.year_mode = year_mode

    def _find_ts_files(self, case_root):
        exact = os.path.join(case_root, self.subdir, f"{self.file_stem}{self.file_suffix}")
        if os.path.isfile(exact):
            return [exact]
        pattern = os.path.join(case_root, self.subdir, f"{self.file_stem}*{self.file_suffix}")
        files = glob.glob(pattern)
        if not files:
            pattern = os.path.join(case_root, self.subdir, "**", f"{self.file_stem}*{self.file_suffix}")
            files = glob.glob(pattern, recursive=True)
        if not files:
            return []
        if self._STAMP_RE:
            with_stamps, no_stamps = [], []
            for f in files:
                s0, _ = self._extract_stamps(f)
                (with_stamps if s0 is not None else no_stamps).append((f, s0))
            with_stamps.sort(key=lambda t: t[1])
            return [f for f, _ in with_stamps] + [f for f in no_stamps]
        files.sort()
        return files

    def _normalize_units_text(self, u):
        if not u: return u
        if isinstance(u, bytes):
            try: u = u.decode("utf-8", errors="ignore")
            except Exception: u = str(u)
        return str(u).replace("$^", "^").replace("$", "").replace("−", "-").strip()

    def _seaice_calendar(self, ds):
        """
        Return calendar string, preferring attrs. If missing, infer from the day counts.
        - 'noleap' if total_days ≈ N*365 exactly
        - '360_day' if total_days ≈ N*360 exactly
        - otherwise 'proleptic_gregorian'
        """
        # If provided anywhere, trust it
        cal = (
            ds["startTime"].attrs.get("calendar")
            or ds["endTime"].attrs.get("calendar")
            or ds.attrs.get("calendar")
        )
        if cal:
            return str(cal)

        # Infer from endTime span
        try:
            total_days = float(ds["endTime"].values[-1])
        except Exception:
            return "proleptic_gregorian"

        # Check closeness to an integer number of years in each calendar
        def is_close_integer(x, tol=1e-8):
            return abs(x - round(x)) < tol

        if is_close_integer(total_days / 365.0):
            return "noleap"
        if is_close_integer(total_days / 360.0):
            return "360_day"
        return "proleptic_gregorian"

    def _seaice_origin(self, ds):
        for key in ("startTime", "endTime"):
            if key in ds.variables:
                u = ds[key].attrs.get("units")
                if u:
                    txt = str(u)
                    m = re.search(r"(days|hours)\s+since\s+(\d{3,4}-\d{2}-\d{2})", txt, flags=re.I)
                    if m: return m.group(1).lower(), m.group(2)
                    m2 = re.search(r"(days|hours)\s+since\s+(\d{3,4})[/-](\d{2})[/-](\d{2})", txt, flags=re.I)
                    if m2: return m2.group(1).lower(), f"{m2.group(2)}-{m2.group(3)}-{m2.group(4)}"
        return "days", "0001-01-01"

    def _construct_time_coord(self, ds):
        # Use midpoint of [startTime, endTime] if available
        if ("startTime" in ds.variables) and ("endTime" in ds.variables):
            st = np.asarray(ds["startTime"].values, dtype="float64")
            et = np.asarray(ds["endTime"].values, dtype="float64")
            mid = 0.5 * (st + et)
            unit_kind, origin = self._seaice_origin(ds)
            calendar = self._seaice_calendar(ds)

            if num2date is not None:
                try:
                    units_str = f"{unit_kind} since {origin} 00:00:00"
                    dates = num2date(mid, units=units_str, calendar=calendar, only_use_cftime_datetimes=True)
                    time = xr.DataArray(dates, dims=("Time",), name="time")
                    time.attrs.update(units=units_str, calendar=calendar)
                    return time
                except Exception:
                    pass

            # numpy fallback
            scale = 86400.0 if unit_kind == "days" else 3600.0
            base = np.datetime64(origin)
            td = (mid.astype("float64") * scale).astype("timedelta64[s]").astype("timedelta64[ns]")
            time = xr.DataArray(base + td, dims=("Time",), name="time")
            time.attrs.update(units=f"{unit_kind} since {origin} 00:00:00", calendar="standard")
            return time

        # Fallback: try any existing reasonable time coord
        for cand in ("time", "Time", "xtime", "xTime"):
            if cand in ds:
                return ds[cand]
        raise KeyError("No startTime/endTime and no recognizable time coordinate found.")

    def _ensure_timecoord(self, ds):
        time_da = self._construct_time_coord(ds)

        # Use ds.sizes instead of ds.dims (future-proof)
        if "Time" in ds.sizes and ds.sizes["Time"] == time_da.size:
            ds = ds.assign_coords(time=("Time", time_da.values)).swap_dims({"Time": "time"})
        else:
            found = None
            for dname, dlen in ds.sizes.items():   # <-- sizes, not dims
                if dlen == time_da.size:
                    found = dname
                    break
            if found is None:
                raise ValueError("Cannot align constructed 'time' to any existing dimension.")
            ds = ds.assign_coords(time=(found, time_da.values)).swap_dims({found: "time"})

        return ds.sortby("time")

    def _aggregate_annual_mean_stream(self, files):
        if not files:
            raise FileNotFoundError("No MPAS time series files found.")

        vname = self.var_name
        vunit = self.var_unit
        vfact = self.var_fact
        chosen_units = None
        ann_list = []

        for f in files:
            ds = xr.open_dataset(f, engine=self.xr_engine, decode_times=False, chunks=self.xr_chunks)
            try:
                # Optional: normalize common sea-ice units text
                for vv in ("iceArea", "iceVolume", "iceThickness"):
                    if vv in ds.variables and "units" in ds[vv].attrs:
                        ds[vv].attrs["units"] = self._normalize_units_text(ds[vv].attrs["units"])

                # Ensure CF-like time coordinate
                ds = self._ensure_timecoord(ds)

                if vname not in ds.variables:
                    raise KeyError(f"Variable '{vname}' not found in {f}")

                da = ds[vname]

                # Apply conversion factor / units override if provided
                if vfact is not None:
                    da = da.astype("float64") * float(vfact)
                if vunit is not None:
                    da = da.copy()  # avoid mutating on-disk attrs
                    da.attrs["units"] = vunit

                # Record units once
                if chosen_units is None:
                    chosen_units = (vunit if vunit is not None else
                                    da.attrs.get("units") or self.default_units.get(vname, ""))

                # Reduce extra dims by mean (keep time)
                extra_dims = [d for d in da.dims if d != "time"]
                if extra_dims:
                    da = da.mean(dim=extra_dims, skipna=True, keep_attrs=True)

                # --- Robust annualization: group by calendar year (no edge-year drop)
                try:
                    year_index = da.time.dt.year
                except Exception:
                    year_index = xr.DataArray(pd.DatetimeIndex(np.array(da.time.values)).year, dims=("time",))

                ann = da.groupby(year_index).mean(keep_attrs=True)

                # Get the group dimension name and coerce it to 'year'
                gdim = ann.dims[0]                    # usually 'group'
                yrs  = ann[gdim].values.astype(int)
                ann  = ann.rename({gdim: "year"}).assign_coords(year=("year", yrs))

                ann_list.append(ann)
            finally:
                ds.close()

        if not ann_list:
            raise ValueError("No annual data produced.")

        # Concatenate all annual blocks and collapse overlaps by year (mean)
        all_ann = xr.concat(ann_list, dim="year")
        grouped = all_ann.groupby("year").mean(skipna=True, keep_attrs=True)

        # Optionally renumber to zero-based per-source years
        years = grouped["year"].values.astype(int)
        if self.year_mode == "zero_based":
            years = np.arange(years.size, dtype=int)
            grouped = grouped.assign_coords(year=("year", years))

        ds_out = xr.Dataset(
            data_vars={vname: (("year",), grouped.values, {"units": chosen_units})},
            coords={"year": ("year", years)},
            attrs={}
        )
        return ds_out

    def build_sources(self, out_dir, suffix="annual", var_name="diag"):
        os.makedirs(out_dir, exist_ok=True)
        spinup_root = os.path.join(self.base_dir, self.spinup_name)
        pctl_root   = os.path.join(self.base_dir, self.picontrol_name)

        spinup_files = self._find_ts_files(spinup_root)
        pctl_files   = self._find_ts_files(pctl_root)
        if not spinup_files:
            raise FileNotFoundError(f"No files under {spinup_root}/{self.subdir} with stem '{self.file_stem}'")
        if not pctl_files:
            raise FileNotFoundError(f"No files under {pctl_root}/{self.subdir} with stem '{self.file_stem}'")

        spinup_nc = os.path.join(out_dir, f"{self.spinup_name}_{var_name}_timeseries_{suffix}.nc")
        pctl_nc   = os.path.join(out_dir, f"{self.picontrol_name}_{var_name}_timeseries_{suffix}.nc")

        if self.force_reprocess:
            for p in (spinup_nc, pctl_nc):
                try: os.remove(p)
                except FileNotFoundError: pass

        if not self._is_fresh(spinup_nc, spinup_files):
            ds_s = self._aggregate_annual_mean_stream(spinup_files)
            enc = {self.var_name: {"zlib": True, "complevel": 4}}
            ds_s.to_netcdf(spinup_nc, engine=self.xr_engine, encoding=enc)

        if not self._is_fresh(pctl_nc, pctl_files):
            ds_p = self._aggregate_annual_mean_stream(pctl_files)
            enc = {self.var_name: {"zlib": True, "complevel": 4}}
            ds_p.to_netcdf(pctl_nc, engine=self.xr_engine, encoding=enc)

        return {"spinup": spinup_nc, "piControl": pctl_nc}

    def _centered_filter(self, y, window=11):
        w = max(1, int(window))
        if w % 2 == 0: w += 1
        pad = w // 2
        y = np.asarray(y, dtype=float)
        isfin = np.isfinite(y).astype(float)
        y_fill = np.where(np.isfinite(y), y, 0.0)
        k = np.ones(w, dtype=float)
        num = np.convolve(y_fill, k, mode="same")
        den = np.convolve(isfin, k, mode="same")
        out = np.full_like(y, np.nan, dtype=float)
        m = den > 0
        out[m] = num[m] / den[m]
        out[:pad] = np.nan; out[-pad:] = np.nan
        return out

    def _sliding_trend(self, y, x, window):
        n = len(x)
        w = max(3, int(window) | 1)  # odd, >=3
        h = w // 2
        out = np.full(n, np.nan, dtype=float)
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        for i in range(n):
            i0 = max(0, i - h); i1 = min(n, i + h + 1)
            xx = x[i0:i1]; yy = y[i0:i1]
            m = np.isfinite(xx) & np.isfinite(yy)
            if m.sum() >= 3:
                out[i] = np.polyfit(xx[m], yy[m], 1)[0]  # per-year slope
        return out

    def _load_year_val(self, ncpath):
        vname = self.var_name
        ds = xr.open_dataset(ncpath, decode_times=False, engine=self.xr_engine)
        try:
            years = ds["year"].values
            vals  = ds[vname].values
            n = min(len(years), len(vals))
            return years[:n], np.asarray(vals[:n], dtype=float)
        finally:
            ds.close()

    def build_combined_tidy(
        self,
        spinup_nc,
        pctl_nc,
        combined_tidy_nc,
        filter_window=11,
        trend_windows=(10,),
        add_per_decade=True,
        sequence="auto",   # "auto" | "calendar" | "append"
        force=False,
        strict=True,
        normalize_year0=False,       # rebase combined 'year' to start at 0
    ):
        deps = [spinup_nc, pctl_nc]
        if (not force) and self._is_fresh_multi(combined_tidy_nc, deps):
            return combined_tidy_nc

        vname = self.var_name
        y_s, v_s = self._load_year_val(spinup_nc)
        y_p, v_p = self._load_year_val(pctl_nc)
        n_s, n_p = len(y_s), len(y_p)

        years_combined = self.make_years(sequence, y_s, y_p).astype(int, copy=False)
        if normalize_year0:
            years_combined = years_combined - int(np.nanmin(years_combined))

        values_combined = np.concatenate([v_s, v_p])

        # simple adjacent-duplicate drop
        if years_combined.size > 1:
            keep = np.ones_like(years_combined, dtype=bool)
            for i in range(len(years_combined) - 1):
                if years_combined[i] == years_combined[i + 1]:
                    keep[i] = False
            if not np.all(keep):
                years_combined  = years_combined[keep]
                values_combined = values_combined[keep]

        n_comb = len(years_combined)
        expected = n_s + n_p
        if strict and n_comb != expected:
            raise ValueError(
                f"[combine] Expected {expected} samples (spinup={n_s} + piControl={n_p}) "
                f"but got {n_comb}. sequence='{sequence}'. "
                f"Spinup years {int(y_s.min())}-{int(y_s.max()) if len(y_s) else -1}, "
                f"piControl years {int(y_p.min())}-{int(y_p.max()) if len(y_p) else -1}."
            )

        # units
        units = self.default_units.get(vname, "")
        try:
            ds_tmp = xr.open_dataset(spinup_nc, engine=self.xr_engine, decode_times=False)
            units = ds_tmp[vname].attrs.get("units", units); ds_tmp.close()
        except Exception:
            pass

        data_vars = {vname: (("time",), values_combined, {"units": units})}
        if filter_window and filter_window > 1:
            data_vars[f"filtered_{int(filter_window)}yr"] = (
                ("time",), self._centered_filter(values_combined, filter_window)
            )
        for W in (trend_windows or ()):
            W = int(W)
            tname = f"sliding_trend_{W}yr"
            trend = self._sliding_trend(values_combined, years_combined, W)
            data_vars[tname] = (("time",), trend)
            if add_per_decade:
                data_vars[f"{tname}_per_decade"] = (("time",), trend * 10.0)

        ds_out = xr.Dataset(
            data_vars=data_vars,
            coords={
                "time": ("time", np.arange(n_comb, dtype=int)),
                "year": ("time", years_combined),
                "model_year": ("time", np.arange(n_comb, dtype=int)),
                "rgn": ("time", np.zeros(n_comb, dtype=int)),
            },
            attrs={
                "series": f"{self.spinup_name} + {self.picontrol_name} (sequence={sequence})",
                "filter_window": int(filter_window) if filter_window else 0,
                "trend_windows": ",".join(str(int(w)) for w in (trend_windows or [])),
                "units": units,
                "spinup_len": int(n_s),
                "picontrol_len": int(n_p),
                "combined_len": int(n_comb),
                "spinup_year_min": int(np.nanmin(y_s)) if n_s else -1,
                "spinup_year_max": int(np.nanmax(y_s)) if n_s else -1,
                "picontrol_year_min": int(np.nanmin(y_p)) if n_p else -1,
                "picontrol_year_max": int(np.nanmax(y_p)) if n_p else -1,
                "combined_year_min": int(np.nanmin(years_combined)) if n_comb else -1,
                "combined_year_max": int(np.nanmax(years_combined)) if n_comb else -1,
            },
        )
        enc = {k: {"zlib": True, "complevel": 4} for k in ds_out.data_vars}
        ds_out.to_netcdf(combined_tidy_nc, engine=self.xr_engine, encoding=enc)
        return combined_tidy_nc


class OHCBuilder(MPASDiagnosticsBuilder):
    """OHC by depth layer annual series."""

    def __init__(
        self,
        base_dir,
        spinup_name,
        picontrol_name,
        subdir="post/time_series",
        file_stem="mpasTimeSeriesOcean",
        stamp_regex=None,
        file_suffix=".nc",
        xr_engine="h5netcdf",
        xr_chunks=None,
        var_name="timeMonthly_avg_avgValueWithinOceanRegion_avgSurfaceTemperature",
        var_unit=None,
        var_fact=None,
        default_units=None,
        force_reprocess=False,
        file_depth="depth.nc",
        region_id=None,
        # OHC-only params (used ONLY when var_name == "OHC")
        ohc_layer_lab="0-700m",
        ohc_layer_int="0-700",
        ohc_ref_temp=None,
        # --- SSH handling ---
        ssh_mode="raw",                 # "raw" (default) | "derived"
        ssh_candidates=None,            # possible SSH variable names in MPAS files
    ):
        self.base_dir = base_dir
        self.spinup_name = spinup_name
        self.picontrol_name = picontrol_name
        self.subdir = subdir

        self.file_stem = file_stem
        self.file_suffix = file_suffix
        self._STAMP_RE = re.compile(stamp_regex) if stamp_regex else None

        self.xr_engine = xr_engine
        self.xr_chunks = xr_chunks
        self.var_name  = var_name
        self.var_unit  = var_unit
        self.var_fact  = var_fact

        self.default_units   = default_units or {}
        self.force_reprocess = force_reprocess
        self.file_depth      = file_depth
        self.region_id       = region_id

        # OHC params (enforce "single band" only when deriving OHC)
        self.ohc_layer_lab = ohc_layer_lab
        self.ohc_layer_int = ohc_layer_int
        self.ohc_ref_temp  = ohc_ref_temp

        # SSH params
        self.ssh_mode = ssh_mode
        self.ssh_candidates = ssh_candidates or [
            # add/adjust to your datasets as needed
            "timeMonthly_avg_ssh",
            "timeMonthly_avg_surfaceHeight",
            "timeMonthly_avg_avgSSH",
            "timeMonthly_avg_avgValueWithinOceanRegion_avgSurfaceHeight",
        ]

    def _find_ts_files(self, case_root):
        exact = os.path.join(case_root, self.subdir, f"{self.file_stem}{self.file_suffix}")
        if os.path.isfile(exact):
            return [exact]
        pattern = os.path.join(case_root, self.subdir, f"{self.file_stem}*{self.file_suffix}")
        files = glob.glob(pattern)
        if not files:
            pattern = os.path.join(case_root, self.subdir, "**", f"{self.file_stem}*{self.file_suffix}")
            files = glob.glob(pattern, recursive=True)
        if not files:
            return []
        if self._STAMP_RE:
            with_stamps, no_stamps = [], []
            for f in files:
                s0, s1 = self._extract_stamps(f)
                (with_stamps if s0 is not None else no_stamps).append((f, s0))
            with_stamps.sort(key=lambda t: t[1])
            return [f for f, _ in with_stamps] + [f for f in no_stamps]
        files.sort()
        return files

    def _compute_year_coord(self, ds, tname):
        v = ds[tname]
        try:
            years = v.dt.year.values
            return years.astype(int, copy=False)
        except Exception:
            pass
        if (v.dtype.kind in ("S", "U")) or ("str" in str(v.dtype).lower()):
            arr = v.values
            if arr.ndim == 2:
                strings = ["".join([(c.decode("utf-8") if isinstance(c, (bytes, np.bytes_)) else str(c)) for c in row]).strip() for row in arr]
            else:
                strings = [(s.decode("utf-8") if isinstance(s, (bytes, np.bytes_)) else str(s)) for s in arr]
            years = [int(s[:4]) if len(s) >= 4 and s[:4].isdigit() else np.nan for s in strings]
            return np.asarray(years, dtype=float)
        units = (v.attrs.get("units") or "").lower()
        if "months since" in units:
            m = re.search(r"months since\s+(\d{1,4})", units)
            base_year = int(m.group(1)) if m else 1
            months = np.asarray(v.values, dtype=float)
            return np.floor(base_year + months / 12.0).astype(int)
        if "days since" in units:
            m = re.search(r"days since\s+(\d{1,4})", units)
            base_year = int(m.group(1)) if m else 1
            days = np.asarray(v.values, dtype=float)
            return np.floor(base_year + days / 365.0).astype(int)
        try:
            return np.asarray(v.values, dtype=int)
        except Exception:
            n = v.sizes[v.dims[0]]
            return np.arange(n, dtype=int)

    def _select_ssh_da(self, ds, tdim):
        """Pick SSH variable directly from dataset."""
        for name in self.ssh_candidates:
            if name in ds:
                da = ds[name]
                if da.dims[0] != tdim:
                    da = da.transpose(tdim, ...)
                return da
        raise KeyError(f"No SSH variable found. Tried: {self.ssh_candidates}. "
                       f"Provide ssh_candidates=[...] or set ssh_mode='derived'.")

    def _derive_ssh(self, ds):
        """Legacy proxy SSH from volume/area (use only if ssh_mode='derived')."""
        # optional scaling (e.g., convert m→cm for SSH)
        if self.var_fact is not None:
            factor = float(self.var_fact)
        else:
            factor = 100.0 # m→cm (proxy is in meters), adjust as needed
        if self.var_unit is not None:
            units = str(self.var_unit)
        else:
            units = "cm"
            
        v_area = "timeMonthly_avg_areaCellGlobal"
        v_vol  = "timeMonthly_avg_volumeCellGlobal"
        for vname in (v_area, v_vol):
            if vname not in ds:
                raise KeyError(f"Variable '{vname}' not found in dataset.")
        tdim = self._get_time_name(ds)
        area = ds[v_area].astype("float64")
        vol  = ds[v_vol ].astype("float64")
        ssh  = (vol / area) * factor
        if not np.isfinite(ssh.values).all():
            raise ValueError("Non-finite values in derived SSH proxy.")
        ssh.name = "SSH"
        ssh.attrs.update({"long_name": "Global mean sea surface height (proxy)", "units": units})
        if ssh.dims[0] != tdim:
            ssh = ssh.transpose(tdim, ...)
        return ssh

    def _parse_band(self, s, bottom_val):
        z0s, z1s = s.split("-")
        z0 = float(z0s.replace("m",""))
        z1 = bottom_val if (z1s.lower() in ["bottom","lower"]) else float(z1s.replace("m",""))
        z0 = max(0.0, z0)
        z1 = max(z0, z1)
        return z0, z1

    def _bottom_depth_value(self):
        if not os.path.exists(self.file_depth):
            raise FileNotFoundError(f"Depth file not found: {self.file_depth}")
        with xr.open_dataset(self.file_depth, engine=self.xr_engine, decode_times=False) as dx:
            if "refBottomDepth" not in dx:
                raise KeyError(f"Depth variable 'refBottomDepth' not found in {self.file_depth}")
            return float(dx["refBottomDepth"].max().values)

    def _derive_ohc_single_band(self, ds):
        """
        Used ONLY when self.var_name == 'OHC' (raw → single-band OHC).
        Returns 1D DataArray named 'OHC_<safe_label>'.
        """
        if isinstance(self.ohc_layer_lab, (list, tuple)):
            if len(self.ohc_layer_lab) != 1:
                raise ValueError("Provide exactly one layer label for OHC derivation.")
            lab = self.ohc_layer_lab[0]
        else:
            lab = str(self.ohc_layer_lab)

        if isinstance(self.ohc_layer_int, (list, tuple)):
            if len(self.ohc_layer_int) != 1:
                raise ValueError("Provide exactly one layer interval for OHC derivation.")
            intr = self.ohc_layer_int[0]
        else:
            intr = str(self.ohc_layer_int)

        rho = 1026.0
        cp  = 3996.0

        v_acel = "timeMonthly_avg_areaCellGlobal"
        v_area = "timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerArea"
        v_lthk = "timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerThickness"
        v_ltmp = "timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature"
        v_mask = "timeMonthly_avg_avgValueWithinOceanLayerRegion_sumLayerMaskValue"

        for v in (v_area, v_lthk, v_ltmp, v_mask, v_acel):
            if v not in ds:
                raise KeyError(f"Variable '{v}' not found in dataset.")
        tdim = self._get_time_name(ds)

        areaCellGlobal = ds[v_acel].astype("float64")
        area = ds[v_area].astype("float64")[:, self.region_id, :]
        thk  = ds[v_lthk].astype("float64")[:, self.region_id, :]
        tmp  = ds[v_ltmp].astype("float64")[:, self.region_id, :]
        msk  = ds[v_mask].astype("float64")[:, self.region_id, :] if v_mask in ds else 1.0
        T = tmp - float(self.ohc_ref_temp) if self.ohc_ref_temp is not None else tmp

        ohc_layer = rho * cp * area * thk * T * msk
        if self.var_fact is not None:
            ohc_layer = ohc_layer * float(self.var_fact)

        vdim_candidates = [d for d in ohc_layer.dims if d != tdim]
        if not vdim_candidates:
            da = ohc_layer
            if da.dims and da.dims[0] != tdim:
                da = da.transpose(tdim, ...)
            out = xr.DataArray(
                da.values, dims=da.dims, coords=da.coords, name=self.var_name,
                attrs={"units": (self.var_unit or "J")}
            )
            if not np.isfinite(out.values).all():
                raise ValueError("Non-finite values in OHC single-band computation.")
            return out
        vdim = vdim_candidates[0]

        bottom_depth = self._bottom_depth_value()
        with xr.open_dataset(self.file_depth, engine=self.xr_engine, decode_times=False) as dx:
            refBot = dx["refBottomDepth"].astype("float64")
        if refBot.dims != (vdim,):
            try:
                refBot = refBot.rename({refBot.dims[0]: vdim})
            except Exception:
                raise ValueError(f"Depth vector dims {refBot.dims} do not match layer dim '{vdim}'.")

        top = xr.concat([xr.zeros_like(refBot.isel({vdim: 0})),
                         refBot.isel({vdim: slice(None, -1)})], dim=vdim)
        bot = refBot

        with np.errstate(divide="ignore", invalid="ignore"):
            heat_per_m = ohc_layer / thk
            heat_per_m = heat_per_m.where(np.isfinite(heat_per_m), 0.0)

        z0, z1 = self._parse_band(intr, bottom_depth)
        overlap = xr.apply_ufunc(np.minimum, bot, xr.DataArray(z1)) - \
                  xr.apply_ufunc(np.maximum, top, xr.DataArray(z0))
        overlap = overlap.clip(min=0.0)

        band_series = (heat_per_m * overlap).sum(dim=vdim)
        if band_series.dims and band_series.dims[0] != tdim:
            band_series = band_series.transpose(tdim, ...)
        out = xr.DataArray(
            band_series.values, dims=band_series.dims, coords=band_series.coords,
            name=self.var_name, attrs={"units": (self.var_unit or "J")}
        )
        if not np.isfinite(out.values).all():
            raise ValueError("Non-finite values in OHC single-band computation.")
        return out

    def _aggregate_annual_mean_stream(self, files):
        if not files:
            raise FileNotFoundError("No MPAS time series files found.")

        year_all, val_all = [], []
        chosen_units = None
        vname = self.var_name
        last_main_name = vname

        for f in files:
            ds = xr.open_dataset(f, engine=self.xr_engine, decode_times=True, chunks=self.xr_chunks)
            try:
                tname = self._get_time_name(ds)

                if vname == "SSH":
                    if self.ssh_mode == "derived":
                        da = self._derive_ssh(ds)
                    else:
                        da = self._select_ssh_da(ds, tname)
                elif vname == "OHC":
                    da = self._derive_ohc_single_band(ds)
                elif vname not in ds:
                    raise KeyError(f"Variable '{vname}' not found in {f}")
                else:
                    da = ds[vname]

                # reduce extra dims by averaging
                non_time_dims = [d for d in da.dims if d != tname]
                if non_time_dims:
                    da = da.mean(dim=non_time_dims, skipna=True, keep_attrs=True)
                    
                vals = np.asarray(da.values, dtype=float)
                if vals.ndim != 1:
                    raise ValueError(f"Expected 1D time series for '{da.name or vname}' from {f}, got {vals.shape}")

                if chosen_units is None:
                    chosen_units = da.attrs.get("units") or \
                                   self.default_units.get(da.name or vname) or \
                                   self.var_unit

                years = self._compute_year_coord(ds, tname)
                year_all.append(np.asarray(years, dtype=np.int64))
                val_all.append(vals)
                last_main_name = da.name or vname
            finally:
                ds.close()

        years = np.concatenate(year_all)
        vals  = np.concatenate(val_all)
        good = np.isfinite(vals) & np.isfinite(years)
        years, vals = years[good], vals[good]
        if years.size == 0:
            raise ValueError("All values are NaN after screening.")

        years_int = np.floor(years).astype(int)
        main_name = last_main_name

        df = pd.DataFrame({main_name: vals, "year": years_int})
        df_annual = df.groupby("year", as_index=False)[main_name].mean(numeric_only=True)

        ds_out = xr.Dataset(
            data_vars={
                main_name: (
                    ("year",),
                    df_annual[main_name].to_numpy(dtype=float),
                    {"units": (chosen_units or self.default_units.get(main_name, self.var_unit or ""))}
                )
            },
            coords={"year": ("year", df_annual["year"].to_numpy(dtype=int))},
            attrs={"aggregation": "annual mean (12 months averaged)", "n_source_files": len(files)},
        )
        return ds_out

    def build_sources(self, out_dir, suffix="annual", var_name="diag"):
        os.makedirs(out_dir, exist_ok=True)
        spinup_files = self._find_ts_files(os.path.join(self.base_dir, self.spinup_name))
        pctl_files   = self._find_ts_files(os.path.join(self.base_dir, self.picontrol_name))

        if not spinup_files:
            raise FileNotFoundError(
                f"No MPAS files under {self.base_dir}/{self.spinup_name}/{self.subdir} with stem '{self.file_stem}'")
        if not pctl_files:
            raise FileNotFoundError(
                f"No MPAS files under {self.base_dir}/{self.picontrol_name}/{self.subdir} with stem '{self.file_stem}'")

        spinup_nc = os.path.join(out_dir, f"{self.spinup_name}_{var_name}_timeseries_{suffix}.nc")
        pctl_nc   = os.path.join(out_dir, f"{self.picontrol_name}_{var_name}_timeseries_{suffix}.nc")

        if self.force_reprocess:
            for path in (spinup_nc, pctl_nc):
                try: os.remove(path)
                except FileNotFoundError: pass

        if not self._is_fresh(spinup_nc, spinup_files):
            ds_s = self._aggregate_annual_mean_stream(spinup_files)
            enc = {k: {"zlib": True, "complevel": 4} for k in ds_s.data_vars}
            ds_s.to_netcdf(spinup_nc, engine=self.xr_engine, encoding=enc)

        if not self._is_fresh(pctl_nc, pctl_files):
            ds_p = self._aggregate_annual_mean_stream(pctl_files)
            enc = {k: {"zlib": True, "complevel": 4} for k in ds_p.data_vars}
            ds_p.to_netcdf(pctl_nc, engine=self.xr_engine, encoding=enc)

        return {"spinup": spinup_nc, "piControl": pctl_nc}

    def _resolve_main_var_name_from_nc(self, ncpath):
        with xr.open_dataset(ncpath, decode_times=False, engine=self.xr_engine) as ds:
            if self.var_name in ds:
                return self.var_name
            ohc_vars = [k for k in ds.data_vars if k.startswith("OHC_")]
            if len(ohc_vars) == 1:
                return ohc_vars[0]
            if self.var_name.startswith("OHC_") and self.var_name in ohc_vars:
                return self.var_name
            if ds.data_vars:
                return list(ds.data_vars)[0]
        raise KeyError(f"No suitable variable found in {ncpath}")

    def _load_year_val(self, ncpath, var_name=None):
        with xr.open_dataset(ncpath, decode_times=False, engine=self.xr_engine) as ds:
            vname = var_name or (self.var_name if self.var_name in ds else None)
            if vname is None:
                ohc_vars = [k for k in ds.data_vars if k.startswith("OHC_")]
                if self.var_name.startswith("OHC_") and self.var_name in ohc_vars:
                    vname = self.var_name
                elif len(ohc_vars) == 1:
                    vname = ohc_vars[0]
                else:
                    raise KeyError(f"Variable '{self.var_name}' not found in {ncpath}. Candidates: {ohc_vars}")
            years = ds["year"].values if "year" in ds.coords else \
                    self._compute_year_coord(ds, self._get_time_name(ds))
            vals  = np.asarray(ds[vname].values, dtype=float)
            years = np.asarray(years, dtype=int)
            n = min(len(years), len(vals))
            return years[:n], vals[:n]

    def build_combined_tidy(
        self,
        spinup_nc,
        pctl_nc,
        combined_tidy_nc,
        filter_window=11,
        trend_windows=(10,),
        add_per_decade=True,
        sequence="auto",
        force=False,
        strict=True,
    ):
        deps = [spinup_nc, pctl_nc]
        if (not force) and self._is_fresh_multi(combined_tidy_nc, deps):
            return combined_tidy_nc

        y_s, v_s = self._load_year_val(spinup_nc)
        y_p, v_p = self._load_year_val(pctl_nc)
        n_s, n_p = len(y_s), len(y_p)

        years_combined = self.make_years(sequence, y_s, y_p).astype(int, copy=False)
        values_all = np.concatenate([v_s, v_p]).astype(float)

        if np.any(np.isfinite(values_all)):
            first_idx = int(np.argmax(np.isfinite(values_all)))
            baseline = values_all[first_idx]
            values_combined = values_all - baseline
        else:
            values_combined = np.full_like(values_all, np.nan)

        if years_combined.size > 1:
            keep = np.ones_like(years_combined, dtype=bool)
            for i in range(len(years_combined) - 1):
                if years_combined[i] == years_combined[i + 1]:
                    keep[i] = False
            if not np.all(keep):
                years_combined  = years_combined[keep]
                values_combined = values_combined[keep]

        n_comb = len(years_combined)
        expected = n_s + n_p
        if strict and n_comb != expected:
            raise ValueError(
                f"[combine] Expected {expected} samples (spinup={n_s} + piControl={n_p}) "
                f"but got {n_comb}. sequence='{sequence}'. "
                f"Spinup years {int(y_s.min())}-{int(y_s.max()) if len(y_s) else -1}, "
                f"piControl years {int(y_p.min())}-{int(y_p.max()) if len(y_p) else -1}."
            )

        main_key_out = self._resolve_main_var_name_from_nc(spinup_nc)
        units = self.default_units.get(main_key_out, self.var_unit or "")
        try:
            with xr.open_dataset(spinup_nc, engine=self.xr_engine, decode_times=False) as ds_tmp:
                if main_key_out in ds_tmp:
                    units = ds_tmp[main_key_out].attrs.get("units", units)
        except Exception:
            pass

        data_vars = {main_key_out: (("time",), values_combined, {"units": units})}

        if filter_window and filter_window > 1:
            data_vars[f"filtered_{int(filter_window)}yr"] = (
                ("time",), self._centered_filter(values_combined, filter_window)
            )

        for W in (trend_windows or ()):
            W = int(W)
            tname = f"sliding_trend_{W}yr"
            trend = self._sliding_trend(values_combined, years_combined, W)
            data_vars[tname] = (("time",), trend)
            if add_per_decade:
                data_vars[f"{tname}_per_decade"] = (("time",), trend * 10.0)

        ds_out = xr.Dataset(
            data_vars=data_vars,
            coords={
                "time": ("time", np.arange(n_comb, dtype=int)),
                "year": ("time", years_combined),
                "model_year": ("time", np.arange(n_comb, dtype=int)),
                "rgn": ("time", np.zeros(n_comb, dtype=int)),
            },
            attrs={
                "series": f"{self.spinup_name} + {self.picontrol_name} (sequence={sequence})",
                "filter_window": int(filter_window) if filter_window else 0,
                "trend_windows": ",".join(str(int(w)) for w in (trend_windows or [])),
                "units": units,
                "spinup_len": int(n_s), "picontrol_len": int(n_p), "combined_len": int(n_comb),
                "spinup_year_min": int(np.nanmin(y_s)) if n_s else -1,
                "spinup_year_max": int(np.nanmax(y_s)) if n_s else -1,
                "picontrol_year_min": int(np.nanmin(y_p)) if n_p else -1,
                "picontrol_year_max": int(np.nanmax(y_p)) if n_p else -1,
                "combined_year_min": int(np.nanmin(years_combined)) if n_comb else -1,
                "combined_year_max": int(np.nanmax(years_combined)) if n_comb else -1,
            },
        )
        enc = {k: {"zlib": True, "complevel": 4} for k in ds_out.data_vars}
        ds_out.to_netcdf(combined_tidy_nc, engine=self.xr_engine, encoding=enc)
        return combined_tidy_nc


class OceanBudgetBuilder(MPASDiagnosticsBuilder):
    """Ocean volume/salt/freshwater conservation diagnostics."""

    def __init__(
        self,
        base_dir,
        spinup_name,
        picontrol_name,
        subdir="post/time_series",
        file_stem="mpasTimeSeriesOcean",
        stamp_regex=None,
        file_suffix=".nc",
        xr_engine="h5netcdf",
        xr_chunks=None,
        var_name="timeMonthly_avg_avgValueWithinOceanRegion_avgSurfaceTemperature",
        default_units=None,
        force_reprocess=False,

        # --- selection for mode="mean" (and generally useful) ---
        ocn_region_map=None,
        region_short=None,  # e.g., "global", "nino3.4"
        region_index=None,  # e.g., 6 for "Global Ocean" (overrides region_short)
        region_dim_candidates=("nOceanRegions", "nRegions", "region", "oceanRegions", "nOceanRegionsTmp"),
        level_index=None,   # e.g., 0 for surface/first layer
        level_dim_candidates=("nVertLevels", "nVertLevelsP1", "nLevels", "z", "depth"),
    ):
        self.base_dir = base_dir
        self.spinup_name = spinup_name
        self.picontrol_name = picontrol_name
        self.subdir = subdir

        self.file_stem = file_stem
        self.file_suffix = file_suffix
        self._STAMP_RE = re.compile(stamp_regex) if stamp_regex else None

        self.xr_engine = xr_engine
        self.xr_chunks = xr_chunks
        self.var_name = var_name

        self.default_units = default_units or {}
        self.force_reprocess = force_reprocess

        # region map (order must match MPAS region indices)
        self.ocn_region_map = ocn_region_map or {
            0: {"name": "Arctic", "short": "arctic"},
            1: {"name": "Equatorial (15S-15N)", "short": "equatorial"},
            2: {"name": "Southern Ocean", "short": "so"},
            3: {"name": "Nino 3",   "short": "nino3"},
            4: {"name": "Nino 4",   "short": "nino4"},
            5: {"name": "Nino 3.4", "short": "nino3.4"},
            6: {"name": "Global Ocean", "short": "global"},
        }

        self.region_dim_candidates = tuple(region_dim_candidates)
        self.level_dim_candidates = tuple(level_dim_candidates)
        self.level_index = level_index

        # Resolve region_index from short name (if needed)
        self.region_short = region_short
        self.region_index = region_index
        if self.region_index is None and self.region_short is not None:
            key = str(self.region_short).lower()
            ridx = None
            for i, info in self.ocn_region_map.items():
                if str(info.get("short", "")).lower() == key:
                    ridx = int(i)
                    break
            if ridx is None:
                raise ValueError(
                    f"Unknown ocean region short='{self.region_short}'. "
                    f"Valid: {[v['short'] for v in self.ocn_region_map.values()]}"
                )
            self.region_index = ridx

    def _find_ts_files(self, case_root):
        exact = os.path.join(case_root, self.subdir, f"{self.file_stem}{self.file_suffix}")
        if os.path.isfile(exact):
            return [exact]

        pattern = os.path.join(case_root, self.subdir, f"{self.file_stem}*{self.file_suffix}")
        files = glob.glob(pattern)
        if not files:
            pattern = os.path.join(case_root, self.subdir, "**", f"{self.file_stem}*{self.file_suffix}")
            files = glob.glob(pattern, recursive=True)

        if not files:
            return []

        if self._STAMP_RE:
            with_stamps, no_stamps = [], []
            for f in files:
                s0, _ = self._extract_stamps(f)
                (with_stamps if s0 is not None else no_stamps).append((f, s0))
            with_stamps.sort(key=lambda t: t[1])
            return [f for f, _ in with_stamps] + [f for f in no_stamps]
        else:
            files.sort()
            return files

    def _compute_year_coord(self, ds, tname):
        v = ds[tname]

        # CFTime / numpy datetime64
        try:
            return v.dt.year.values.astype(int, copy=False)
        except Exception:
            pass

        # string time
        if (v.dtype.kind in ("S", "U")) or ("str" in str(v.dtype).lower()):
            arr = v.values
            if arr.ndim == 2:
                strings = []
                for row in arr:
                    s = "".join(
                        [(c.decode("utf-8") if isinstance(c, (bytes, np.bytes_)) else str(c)) for c in row]
                    ).strip()
                    strings.append(s)
            else:
                strings = [(s.decode("utf-8") if isinstance(s, (bytes, np.bytes_)) else str(s)) for s in v.values]
            years = [int(s[:4]) if len(s) >= 4 and s[:4].isdigit() else np.nan for s in strings]
            return np.asarray(years, dtype=float)

        # numeric time with units
        units = (v.attrs.get("units") or "").lower()
        if "months since" in units:
            m = re.search(r"months since\s+(\d{1,4})", units)
            base_year = int(m.group(1)) if m else 1
            months = np.asarray(v.values, dtype=float)
            return np.floor(base_year + months / 12.0).astype(int)
        if "days since" in units:
            m = re.search(r"days since\s+(\d{1,4})", units)
            base_year = int(m.group(1)) if m else 1
            days = np.asarray(v.values, dtype=float)
            return np.floor(base_year + days / 365.0).astype(int)

        # fallback: simple index
        n = v.sizes[v.dims[0]]
        return np.arange(n, dtype=int)

    def _maybe_select_region_level(self, da, tname):
        # region select
        if self.region_index is not None:
            for rd in self.region_dim_candidates:
                if rd in da.dims:
                    da = da.isel({rd: int(self.region_index)})
                    break

        # level select
        if self.level_index is not None:
            for zd in self.level_dim_candidates:
                if zd in da.dims:
                    da = da.isel({zd: int(self.level_index)})
                    break

        return da

    def _get_ocean_area_m2(self, files, area_cell_name="areaCell"):
        if not files:
            raise FileNotFoundError("No MPAS time series files found (needed for ocean area).")
        f0 = files[0]
        ds = xr.open_dataset(f0, engine=self.xr_engine, decode_times=False, chunks=self.xr_chunks)
        try:
            if area_cell_name not in ds:
                raise KeyError(
                    f"'{area_cell_name}' not found in {os.path.basename(f0)}. "
                    f"Provide ocean_area_m2 explicitly or ensure a file with '{area_cell_name}' is present."
                )
            A = ds[area_cell_name]
            return float(A.sum().values)
        finally:
            ds.close()

    @staticmethod
    def _first_finite_index(arr):
        arr = np.asarray(arr, dtype=float)
        m = np.isfinite(arr)
        if not m.any():
            return None
        return int(np.argmax(m))

    def _aggregate_annual_mean_stream(self, files):
        if not files:
            raise FileNotFoundError("No MPAS time series files found.")

        year_all, val_all = [], []
        vname = self.var_name
        chosen_units = None

        for f in files:
            ds = xr.open_dataset(f, engine=self.xr_engine, decode_times=False, chunks=self.xr_chunks)
            try:
                tname = self._get_time_name(ds)
                years = self._compute_year_coord(ds, tname)

                da = ds[vname]
                da = self._maybe_select_region_level(da, tname)

                if chosen_units is None:
                    chosen_units = da.attrs.get("units") or self.default_units.get(vname)

                vals = da.values
                if vals.ndim != 1:
                    time_ax = da.get_axis_num(tname) if tname in da.dims else 0
                    if time_ax != 0:
                        vals = np.moveaxis(vals, time_ax, 0)
                    while vals.ndim > 1:
                        vals = np.nanmean(vals, axis=-1)

                year_all.append(np.asarray(years, dtype=float))
                val_all.append(np.asarray(vals, dtype=float))
            finally:
                ds.close()

        years = np.concatenate(year_all)
        vals = np.concatenate(val_all)
        good = np.isfinite(vals) & np.isfinite(years)
        years, vals = years[good], vals[good]
        if years.size == 0:
            raise ValueError("All annual values are NaN after screening.")

        years = years.astype(np.int64, copy=False)

        y0, y1 = int(np.nanmin(years)), int(np.nanmax(years))
        bins = y1 - y0 + 1
        idx = years - y0
        ssum = np.bincount(idx, weights=vals, minlength=bins)
        scount = np.bincount(idx, minlength=bins)

        annual = np.full(bins, np.nan, dtype=float)
        m = scount > 0
        annual[m] = ssum[m] / scount[m]

        year_coord = np.arange(y0, y1 + 1, dtype=int)
        return xr.Dataset(
            data_vars={vname: (("year",), annual, {"units": chosen_units or self.default_units.get(vname, "")})},
            coords={"year": ("year", year_coord)},
            attrs={
                "region_short": str(self.region_short) if self.region_short is not None else "",
                "region_index": int(self.region_index) if self.region_index is not None else -1,
                "level_index": int(self.level_index) if self.level_index is not None else -1,
            }
        )

    def _aggregate_annual_salt_mass_stream(
        self,
        files,
        salinity_name="salinity",
        layer_thickness_name="layerThickness",
        area_name="areaCell",
        rho0=1026.0,
        out_var="ocean_salt_mass_kg",
    ):
        """
        Compute annual-mean GLOBAL ocean salt mass [kg] from MPAS-O cell-wise fields:
          M_salt(t) = rho0 * sum_{cells,k} S_kgkg * areaCell * layerThickness

        IMPORTANT:
          This assumes salinity is in kg/kg (MPAS often labels as "1.e-3").
          Do NOT divide by 1000 here. Only convert by *1000 when you want psu.
        """
        if not files:
            raise FileNotFoundError("No MPAS time series files found.")

        year_all, salt_all = [], []

        for f in files:
            ds = xr.open_dataset(f, engine=self.xr_engine, decode_times=False, chunks=self.xr_chunks)
            try:
                tname = self._get_time_name(ds)
                years = np.asarray(self._compute_year_coord(ds, tname), dtype=float)

                S = ds[salinity_name]
                H = ds[layer_thickness_name]
                A = ds[area_name]

                vol = H * A
                salt_ts = rho0 * (S * vol).sum(dim=[d for d in S.dims if d != tname])

                salt_vals = np.asarray(salt_ts.values, dtype=float)
                while salt_vals.ndim > 1:
                    salt_vals = np.nanmean(salt_vals, axis=-1)

                year_all.append(years)
                salt_all.append(salt_vals)
            finally:
                ds.close()

        years = np.concatenate(year_all)
        vals = np.concatenate(salt_all)

        good = np.isfinite(vals) & np.isfinite(years)
        years, vals = years[good], vals[good]
        if years.size == 0:
            raise ValueError("All salt-mass values are NaN after screening.")

        years = years.astype(np.int64, copy=False)

        y0, y1 = int(np.nanmin(years)), int(np.nanmax(years))
        bins = y1 - y0 + 1
        idx = years - y0

        ssum = np.bincount(idx, weights=vals, minlength=bins)
        scount = np.bincount(idx, minlength=bins)

        annual = np.full(bins, np.nan, dtype=float)
        m = scount > 0
        annual[m] = ssum[m] / scount[m]

        year_coord = np.arange(y0, y1 + 1, dtype=int)
        return xr.Dataset(
            data_vars={out_var: (("year",), annual, {"units": "kg", "long_name": "total ocean salt mass (rho0 const)"})},
            coords={"year": ("year", year_coord)},
            attrs={
                "rho0_kg_m3": float(rho0),
                "salinity_units_assumed": "kg/kg (MPAS '1.e-3' convention)",
                "mass_definition": "rho0 * sum(S_kgkg * areaCell * layerThickness)",
            },
        )

    def _aggregate_annual_ocean_layerregion_primitives_stream(
        self,
        files,
        vol_global_name="timeMonthly_avg_volumeCellGlobal",
        sum_mask_name="timeMonthly_avg_avgValueWithinOceanLayerRegion_sumLayerMaskValue",
        area_name="timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerArea",
        thick_name="timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerThickness",
        sal_name="timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerSalinity",
        region_dim="nOceanRegionsTmp",
        z_dim="nVertLevels",

        choose_region="auto",   # "auto" | "all" | int | str short
        z_sel=None,             # None | int | slice | list/np.ndarray of indices

        area_global_name="timeMonthly_avg_areaCellGlobal",
        use_area_global_if_available=True,

        out_volume="ocean_volume",
        out_volume_reconstructed="ocean_volume_reconstructed",
        out_mean_salinity="mean_salinity_from_integral",
        out_salt_proxy_kgkg_m3="salt_content_proxy_kgkg_m3",
        out_ocean_area_ann="ocean_area_m2_ann",
    ):
        if not files:
            raise FileNotFoundError("No MPAS time series files found.")

        year_all = []
        Vg_all = []
        Vrec_all = []
        Sbar_all = []
        saltproxy_all = []
        Ao_all = []

        choose_region_resolved = choose_region
        if isinstance(choose_region, str) and choose_region not in ("auto", "all"):
            key = choose_region.lower()
            ridx = None
            for i, info in self.ocn_region_map.items():
                if str(info.get("short", "")).lower() == key:
                    ridx = int(i)
                    break
            if ridx is None:
                raise ValueError(
                    f"Unknown ocean region short='{choose_region}'. "
                    f"Valid: {[v['short'] for v in self.ocn_region_map.values()]} plus 'auto'/'all'."
                )
            choose_region_resolved = ridx

        for f in files:
            ds = xr.open_dataset(f, engine=self.xr_engine, decode_times=False, chunks=self.xr_chunks)
            try:
                tname = self._get_time_name(ds)
                years = np.asarray(self._compute_year_coord(ds, tname), dtype=float)

                Vg = ds[vol_global_name]   # (Time,)

                Ao_ts = None
                if use_area_global_if_available and (area_global_name in ds.variables):
                    Ao_ts = ds[area_global_name]  # (Time,)

                mask = ds[sum_mask_name]   # (Time, region, z)
                Aavg = ds[area_name]       # (Time, region, z)
                Havg = ds[thick_name]      # (Time, region, z)
                S = ds[sal_name]           # (Time, region, z), "1.e-3" -> kg/kg

                if z_sel is not None:
                    mask = mask.isel({z_dim: z_sel})
                    Aavg = Aavg.isel({z_dim: z_sel})
                    Havg = Havg.isel({z_dim: z_sel})
                    S    = S.isel({z_dim: z_sel})

                V_rk = (Aavg * mask) * Havg  # (Time, region, z)

                if choose_region_resolved == "all":
                    V_rec = V_rk.sum(dim=(region_dim, z_dim))
                    salt_int_kgkg_m3 = (S * V_rk).sum(dim=(region_dim, z_dim))
                    Sbar_kgkg = salt_int_kgkg_m3 / V_rk.sum(dim=(region_dim, z_dim))

                elif isinstance(choose_region_resolved, int):
                    Vr = V_rk.isel({region_dim: int(choose_region_resolved)})
                    Sr = S.isel({region_dim: int(choose_region_resolved)})
                    V_rec = Vr.sum(dim=z_dim)
                    salt_int_kgkg_m3 = (Sr * Vr).sum(dim=z_dim)
                    Sbar_kgkg = salt_int_kgkg_m3 / Vr.sum(dim=z_dim)

                else:
                    V_reg = V_rk.sum(dim=z_dim)  # (Time, region)
                    diff = abs(V_reg.mean(dim=tname) - Vg.mean(dim=tname))
                    r0 = int(diff.argmin().values)

                    Vr = V_rk.isel({region_dim: r0})
                    Sr = S.isel({region_dim: r0})
                    V_rec = Vr.sum(dim=z_dim)
                    salt_int_kgkg_m3 = (Sr * Vr).sum(dim=z_dim)
                    Sbar_kgkg = salt_int_kgkg_m3 / Vr.sum(dim=z_dim)

                Sbar_psu = 1000.0 * Sbar_kgkg

                year_all.append(years)
                Vg_all.append(np.asarray(Vg.values, dtype=float))
                Vrec_all.append(np.asarray(V_rec.values, dtype=float))
                Sbar_all.append(np.asarray(Sbar_psu.values, dtype=float))
                saltproxy_all.append(np.asarray(salt_int_kgkg_m3.values, dtype=float))
                if Ao_ts is not None:
                    Ao_all.append(np.asarray(Ao_ts.values, dtype=float))
            finally:
                ds.close()

        years = np.concatenate(year_all)
        Vg = np.concatenate(Vg_all)
        Vrec = np.concatenate(Vrec_all)
        Sbar = np.concatenate(Sbar_all)
        saltproxy = np.concatenate(saltproxy_all)

        good = np.isfinite(years) & np.isfinite(Vg) & np.isfinite(Vrec) & np.isfinite(Sbar) & np.isfinite(saltproxy)
        years, Vg, Vrec, Sbar, saltproxy = years[good], Vg[good], Vrec[good], Sbar[good], saltproxy[good]
        if years.size == 0:
            raise ValueError("All values are NaN after screening.")

        years = years.astype(np.int64, copy=False)

        y0, y1 = int(np.nanmin(years)), int(np.nanmax(years))
        bins = y1 - y0 + 1
        idx = years - y0

        def _ann_mean(arr):
            ssum = np.bincount(idx, weights=arr, minlength=bins)
            cnt = np.bincount(idx, minlength=bins)
            out = np.full(bins, np.nan, dtype=float)
            m = cnt > 0
            out[m] = ssum[m] / cnt[m]
            return out, cnt

        Vg_ann, cnt = _ann_mean(Vg)
        Vrec_ann, _ = _ann_mean(Vrec)
        Sbar_ann, _ = _ann_mean(Sbar)
        sp_ann, _ = _ann_mean(saltproxy)

        year_coord = np.arange(y0, y1 + 1, dtype=int)

        ds_out = xr.Dataset(
            data_vars={
                out_volume: (("year",), Vg_ann, {"units": "m^3", "long_name": "total ocean volume (from volumeCellGlobal)"}),
                out_volume_reconstructed: (("year",), Vrec_ann, {"units": "m^3", "long_name": "reconstructed volume from layer-region stats (debug)"}),
                out_mean_salinity: (("year",), Sbar_ann, {"units": "psu", "long_name": "global mean salinity (volume-weighted, from layer-region integral)"}),
                out_salt_proxy_kgkg_m3: (("year",), sp_ann, {"units": "kg/kg*m^3", "long_name": "sum(S_kgkg * V) proxy (integral, before psu scaling)"}),
            },
            coords={"year": ("year", year_coord)},
            attrs={
                "choose_region": str(choose_region),
                "choose_region_resolved": str(choose_region_resolved),
                "z_sel": str(z_sel) if z_sel is not None else "all_levels",
                "salinity_units_assumed": "1.e-3 (kg/kg) -> psu via *1000",
                "area_global_name": str(area_global_name),
            },
        )

        if len(Ao_all) > 0:
            Ao = np.concatenate(Ao_all)
            Ao = Ao[good]
            Ao_sum = np.bincount(idx, weights=Ao, minlength=bins)
            Ao_ann = np.full(bins, np.nan, dtype=float)
            m = cnt > 0
            Ao_ann[m] = Ao_sum[m] / cnt[m]
            ds_out[out_ocean_area_ann] = (("year",), Ao_ann, {"units": "m^2", "long_name": "annual mean ocean area from areaCellGlobal"})

        return ds_out

    def build_combined_layerregion_annual(
        self,
        spinup_nc,
        pctl_nc,
        combined_annual_nc,
        sequence="append",
        force=False,

        ocean_area_m2=None,
        area_cell_name="areaCell",
        use_area_ann_if_available=True,

        rho0=1026.0,

        in_volume="ocean_volume",
        in_volume_reconstructed="ocean_volume_reconstructed",
        in_mean_salinity="mean_salinity_from_integral",
        in_salt_proxy="salt_content_proxy_kgkg_m3",
        in_ocean_area_ann="ocean_area_m2_ann",

        # freshwater-equivalent from salt conservation
        S0_kgkg=None,                 # reference salinity mass fraction (kg/kg), e.g. 34.7e-3
        use_S0_from_first_year=True,  # if True and mean salinity exists, use baseline (first finite) Sbar as S0
        out_salt_fw_mm="ocean_salt_freshwater_equiv_mm",
        out_salt_fw_rate_mm_yr="ocean_salt_freshwater_equiv_mm_per_yr",
    ):
        deps = [spinup_nc, pctl_nc]
        if (not force) and self._is_fresh_multi(combined_annual_nc, deps):
            return combined_annual_nc

        ds_s = xr.open_dataset(spinup_nc, decode_times=False, engine=self.xr_engine)
        ds_p = xr.open_dataset(pctl_nc, decode_times=False, engine=self.xr_engine)
        try:
            for k in (in_volume, in_volume_reconstructed, in_mean_salinity, in_salt_proxy):
                if k not in ds_s.data_vars or k not in ds_p.data_vars:
                    raise KeyError(
                        f"Missing '{k}' in annual per-case files. "
                        f"spinup has: {list(ds_s.data_vars)}; piControl has: {list(ds_p.data_vars)}"
                    )

            y_s = np.asarray(ds_s["year"].values, dtype=int)
            y_p = np.asarray(ds_p["year"].values, dtype=int)
            y_comb = self.make_years(sequence, y_s, y_p).astype(int, copy=False)

            def cat(var):
                return np.concatenate([np.asarray(ds_s[var].values, dtype=float),
                                       np.asarray(ds_p[var].values, dtype=float)])

            Vg   = cat(in_volume)
            Vrec = cat(in_volume_reconstructed)
            Sbar = cat(in_mean_salinity)  # psu
            sp   = cat(in_salt_proxy)     # kg/kg*m^3

            # --- Ao_used (constant) ---
            ocean_area_source = "unknown"
            Ao_used = None

            if ocean_area_m2 is not None:
                Ao_used = float(ocean_area_m2)
                ocean_area_source = "user_provided"
            elif use_area_ann_if_available and (in_ocean_area_ann in ds_s.data_vars) and (in_ocean_area_ann in ds_p.data_vars):
                Ao_comb = np.concatenate([np.asarray(ds_s[in_ocean_area_ann].values, dtype=float),
                                          np.asarray(ds_p[in_ocean_area_ann].values, dtype=float)])
                iA = self._first_finite_index(Ao_comb)
                if iA is not None:
                    Ao_used = float(Ao_comb[iA])
                    ocean_area_source = in_ocean_area_ann

            if Ao_used is None:
                spinup_root = os.path.join(self.base_dir, self.spinup_name)
                spinup_files = self._find_ts_files(spinup_root)
                try:
                    Ao_used = float(self._get_ocean_area_m2(spinup_files, area_cell_name=area_cell_name))
                    ocean_area_source = f"sum({area_cell_name})"
                except KeyError as e:
                    raise KeyError(
                        f"Could not determine ocean area. '{area_cell_name}' not found in timeSeries files "
                        f"and no '{in_ocean_area_ann}' present. Provide ocean_area_m2 explicitly."
                    ) from e

            # --- Closure diagnostics (psu*m^3) ---
            salt_int_psu_m3 = 1000.0 * sp
            salt_from_sbar_psu_m3 = Sbar * Vg
            with np.errstate(invalid="ignore", divide="ignore"):
                salt_resid_frac = (salt_int_psu_m3 - salt_from_sbar_psu_m3) / salt_from_sbar_psu_m3

            # --- Volume baseline (first finite) ---
            iV0 = self._first_finite_index(Vg)
            if iV0 is None:
                raise ValueError("No finite ocean_volume values in combined series.")
            V0 = Vg[iV0]
            dV = Vg - V0
            with np.errstate(invalid="ignore", divide="ignore"):
                dV_frac = dV / V0

            dh_mm = dV / Ao_used * 1000.0
            dh_rate_mm_yr = np.full_like(dh_mm, np.nan, dtype=float)
            if len(Vg) >= 2:
                dh_rate_mm_yr[1:] = (Vg[1:] - Vg[:-1]) / Ao_used * 1000.0

            # --- Salt mass and baseline (first finite) ---
            salt_mass_kg = rho0 * sp
            iS0 = self._first_finite_index(salt_mass_kg)
            if iS0 is None:
                raise ValueError("No finite ocean_salt_mass_kg values in combined series.")
            salt_mass_change_kg = salt_mass_kg - salt_mass_kg[iS0]
            with np.errstate(invalid="ignore", divide="ignore"):
                salt_mass_change_frac = salt_mass_change_kg / salt_mass_kg[iS0]

            # -----------------------------
            # Freshwater-equivalent from salt conservation
            # -----------------------------
            # dV_fw = dM_salt / (rho0 * S0)
            # h_fw_mm = dV_fw / Ao * 1000
            if S0_kgkg is not None:
                S0_used = float(S0_kgkg)
                S0_source = "user_provided"
            elif use_S0_from_first_year and np.isfinite(Sbar[iS0]):
                S0_used = float(Sbar[iS0]) / 1000.0  # psu -> kg/kg
                S0_source = "baseline_mean_salinity"
            else:
                S0_used = 34.7e-3
                S0_source = "default_34.7psu"

            with np.errstate(invalid="ignore", divide="ignore"):
                salt_fw_mm = (salt_mass_change_kg / (rho0 * S0_used * Ao_used)) * 1000.0

            salt_fw_rate_mm_yr = np.full_like(salt_fw_mm, np.nan, dtype=float)
            if len(salt_mass_kg) >= 2:
                with np.errstate(invalid="ignore", divide="ignore"):
                    salt_fw_rate_mm_yr[1:] = (salt_mass_kg[1:] - salt_mass_kg[:-1]) / (rho0 * S0_used * Ao_used) * 1000.0

            ds_out = xr.Dataset(
                data_vars={
                    "ocean_volume": (("year",), Vg, {"units": "m^3"}),
                    "ocean_volume_reconstructed": (("year",), Vrec, {"units": "m^3"}),
                    "mean_salinity_from_integral": (("year",), Sbar, {"units": "psu"}),
                    "salt_content_proxy_kgkg_m3": (("year",), sp, {"units": "kg/kg*m^3"}),

                    "salt_content_integral_psu_m3": (("year",), salt_int_psu_m3, {"units": "psu*m^3"}),
                    "salt_content_from_sbar_psu_m3": (("year",), salt_from_sbar_psu_m3, {"units": "psu*m^3"}),
                    "salt_content_residual_frac": (("year",), salt_resid_frac, {"units": "1"}),

                    "ocean_volume_change_m3": (("year",), dV, {"units": "m^3"}),
                    "ocean_volume_change_frac": (("year",), dV_frac, {"units": "1"}),

                    "ocean_volume_change_equiv_mm": (("year",), dh_mm, {"units": "mm"}),
                    "ocean_volume_tendency_equiv_mm_per_yr": (("year",), dh_rate_mm_yr, {"units": "mm/yr"}),

                    "ocean_salt_mass_kg": (("year",), salt_mass_kg, {"units": "kg"}),
                    "ocean_salt_mass_change_kg": (("year",), salt_mass_change_kg, {"units": "kg"}),
                    "ocean_salt_mass_change_frac": (("year",), salt_mass_change_frac, {"units": "1"}),

                    out_salt_fw_mm: (("year",), salt_fw_mm, {
                        "units": "mm",
                        "long_name": "freshwater-equivalent thickness implied by salt drift (vs salt baseline)"
                    }),
                    out_salt_fw_rate_mm_yr: (("year",), salt_fw_rate_mm_yr, {
                        "units": "mm/yr",
                        "long_name": "freshwater-equivalent tendency implied by salt drift (year-to-year)"
                    }),
                },
                coords={"year": ("year", y_comb)},
                attrs={
                    "series": f"{self.spinup_name} + {self.picontrol_name} (sequence={sequence})",
                    "ocean_area_m2": float(Ao_used),
                    "ocean_area_source": str(ocean_area_source),
                    "rho0_kg_m3_for_salt_mass": float(rho0),
                    "baseline_year_index_volume": int(iV0),
                    "baseline_year_index_salt": int(iS0),
                    "S0_kgkg_for_fw_equiv": float(S0_used),
                    "S0_source_for_fw_equiv": str(S0_source),
                    "fw_equiv_definition": "dV_fw = dM_salt/(rho0*S0); h_fw_mm = dV_fw/Ao*1000",
                },
            )

            enc = {k: {"zlib": True, "complevel": 4} for k in ds_out.data_vars}
            ds_out.to_netcdf(combined_annual_nc, engine=self.xr_engine, encoding=enc)
            return combined_annual_nc
        finally:
            ds_s.close()
            ds_p.close()

    def build_sources(self, out_dir, suffix="annual", var_name="diag", mode="mean", **mode_kwargs):
        os.makedirs(out_dir, exist_ok=True)
        spinup_root = os.path.join(self.base_dir, self.spinup_name)
        pctl_root = os.path.join(self.base_dir, self.picontrol_name)

        spinup_files = self._find_ts_files(spinup_root)
        pctl_files = self._find_ts_files(pctl_root)

        if not spinup_files:
            raise FileNotFoundError(
                f"No MPAS files found under {spinup_root}/{self.subdir} with stem '{self.file_stem}'"
            )
        if not pctl_files:
            raise FileNotFoundError(
                f"No MPAS files found under {pctl_root}/{self.subdir} with stem '{self.file_stem}'"
            )

        spinup_nc = os.path.join(out_dir, f"{self.spinup_name}_{var_name}_timeseries_{suffix}.nc")
        pctl_nc = os.path.join(out_dir, f"{self.picontrol_name}_{var_name}_timeseries_{suffix}.nc")

        # --- FIX 1: ensure force_reprocess truly forces rebuild (skip freshness short-circuit) ---
        if self.force_reprocess:
            for path in (spinup_nc, pctl_nc):
                try:
                    os.remove(path)
                except FileNotFoundError:
                    pass

        def _write_if_needed(target_nc, src_files):
            # If forcing, always rebuild (do not rely on mtimes).
            if (not self.force_reprocess) and self._is_fresh(target_nc, src_files):
                return

            # If forcing but file exists (e.g., deletion failed), remove it to avoid append/partial artifacts.
            if self.force_reprocess and os.path.exists(target_nc):
                try:
                    os.remove(target_nc)
                except Exception:
                    pass

            if mode == "mean":
                ds_out = self._aggregate_annual_mean_stream(src_files)
                enc = {self.var_name: {"zlib": True, "complevel": 4}}

            elif mode == "salt_mass":
                out_var = mode_kwargs.get("out_var", "ocean_salt_mass_kg")
                ds_out = self._aggregate_annual_salt_mass_stream(
                    src_files,
                    out_var=out_var,
                    **{k: v for k, v in mode_kwargs.items() if k != "out_var"}
                )
                enc = {out_var: {"zlib": True, "complevel": 4}}

            elif mode == "ocean_vol_sbar_layerregion":
                ds_out = self._aggregate_annual_ocean_layerregion_primitives_stream(src_files, **mode_kwargs)
                enc = {k: {"zlib": True, "complevel": 4} for k in ds_out.data_vars}

            else:
                raise ValueError(f"Unknown mode='{mode}'. Use 'mean', 'salt_mass', or 'ocean_vol_sbar_layerregion'.")

            ds_out.to_netcdf(target_nc, engine=self.xr_engine, encoding=enc)

        _write_if_needed(spinup_nc, spinup_files)
        _write_if_needed(pctl_nc, pctl_files)

        return {"spinup": spinup_nc, "piControl": pctl_nc}

    def _load_year_val(self, ncpath, var_name=None):
        vname = var_name or self.var_name
        ds = xr.open_dataset(ncpath, decode_times=False, engine=self.xr_engine)
        try:
            if "year" in ds.coords:
                years = ds["year"].values
                vals = ds[vname].values
            else:
                tname = self._get_time_name(ds)
                years = self._compute_year_coord(ds, tname)
                vals = ds[vname].values
            vals = np.asarray(vals, dtype=float)
            years = np.asarray(years, dtype=float)
            good = np.isfinite(years) & np.isfinite(vals)
            years = years[good].astype(int, copy=False)
            vals = vals[good]
            n = min(len(years), len(vals))
            return years[:n], vals[:n]
        finally:
            ds.close()

    def build_combined_tidy(
        self,
        spinup_nc,
        pctl_nc,
        combined_tidy_nc,
        filter_window=11,
        trend_windows=(10,),
        add_per_decade=True,
        sequence="auto",
        force=False,
        strict=True,
        var_name=None,
    ):
        deps = [spinup_nc, pctl_nc]
        if (not force) and self._is_fresh_multi(combined_tidy_nc, deps):
            return combined_tidy_nc

        vname = var_name or self.var_name

        y_s, v_s = self._load_year_val(spinup_nc, var_name=vname)
        y_p, v_p = self._load_year_val(pctl_nc, var_name=vname)
        n_s, n_p = len(y_s), len(y_p)

        years_combined = self.make_years(sequence, y_s, y_p).astype(int, copy=False)
        values_combined = np.concatenate([v_s, v_p])

        if years_combined.size > 1:
            keep = np.ones_like(years_combined, dtype=bool)
            for i in range(len(years_combined) - 1):
                if years_combined[i] == years_combined[i + 1]:
                    keep[i] = False
            if not np.all(keep):
                years_combined = years_combined[keep]
                values_combined = values_combined[keep]

        n_comb = len(years_combined)
        expected = n_s + n_p
        if strict and n_comb != expected:
            raise ValueError(
                f"[combine] Expected {expected} samples (spinup={n_s} + piControl={n_p}) "
                f"but got {n_comb}. sequence='{sequence}'. "
                f"Spinup years {int(y_s.min())}-{int(y_s.max()) if len(y_s) else -1}, "
                f"piControl years {int(y_p.min())}-{int(y_p.max()) if len(y_p) else -1}."
            )

        model_year = np.arange(n_comb, dtype=int)
        rgn = np.zeros(n_comb, dtype=int)

        units = self.default_units.get(vname, "")
        try:
            ds_tmp = xr.open_dataset(spinup_nc, engine=self.xr_engine, decode_times=False)
            units = ds_tmp[vname].attrs.get("units", units)
            ds_tmp.close()
        except Exception:
            pass

        data_vars = {vname: (("time",), values_combined, {"units": units})}

        if filter_window and filter_window > 1:
            data_vars[f"filtered_{int(filter_window)}yr"] = (("time",), self._centered_filter(values_combined, filter_window))

        for W in (trend_windows or ()):
            W = int(W)
            tname = f"sliding_trend_{W}yr"
            trend = self._sliding_trend(values_combined, years_combined, W)
            data_vars[tname] = (("time",), trend)
            if add_per_decade:
                data_vars[f"{tname}_per_decade"] = (("time",), trend * 10.0)

        ds_out = xr.Dataset(
            data_vars=data_vars,
            coords={
                "time": ("time", np.arange(n_comb, dtype=int)),
                "year": ("time", years_combined),
                "model_year": ("time", model_year),
                "rgn": ("time", rgn),
            },
            attrs={
                "series": f"{self.spinup_name} + {self.picontrol_name} (sequence={sequence})",
                "filter_window": int(filter_window) if filter_window else 0,
                "trend_windows": ",".join(str(int(w)) for w in (trend_windows or [])),
                "units": units,
                "spinup_len": int(n_s),
                "picontrol_len": int(n_p),
                "combined_len": int(n_comb),
                "spinup_year_min": int(np.nanmin(y_s)) if n_s else -1,
                "spinup_year_max": int(np.nanmax(y_s)) if n_s else -1,
                "picontrol_year_min": int(np.nanmin(y_p)) if n_p else -1,
                "picontrol_year_max": int(np.nanmax(y_p)) if n_p else -1,
                "combined_year_min": int(np.nanmin(years_combined)) if n_comb else -1,
                "combined_year_max": int(np.nanmax(years_combined)) if n_comb else -1,
                "region_short": str(self.region_short) if self.region_short is not None else "",
                "region_index": int(self.region_index) if self.region_index is not None else -1,
                "level_index": int(self.level_index) if self.level_index is not None else -1,
            },
        )
        enc = {k: {"zlib": True, "complevel": 4} for k in ds_out.data_vars}
        ds_out.to_netcdf(combined_tidy_nc, engine=self.xr_engine, encoding=enc)
        return combined_tidy_nc

    def build_combined_tidy_for_vars(
        self,
        spinup_nc,
        pctl_nc,
        out_dir,
        experiment,
        suffix="annual",
        spinup_period=(0, 0),
        vars_to_build=("ocean_volume", "mean_salinity_from_integral"),
        filter_window=11,
        trend_windows=(10, 30),
        add_per_decade=True,
        sequence="append",
        force=False,
        strict=True,
    ):
        os.makedirs(out_dir, exist_ok=True)
        outputs = {}
        for var in vars_to_build:
            combined_tidy = os.path.join(
                out_dir,
                f"{var}_{experiment}_combined_{suffix}_{spinup_period[0]}-{spinup_period[1]}.nc"
            )
            out_path = self.build_combined_tidy(
                spinup_nc=spinup_nc,
                pctl_nc=pctl_nc,
                combined_tidy_nc=combined_tidy,
                filter_window=filter_window,
                trend_windows=trend_windows,
                add_per_decade=add_per_decade,
                sequence=sequence,
                force=force,
                strict=strict,
                var_name=var,
            )
            outputs[var] = out_path
        return outputs

    def summarize_tidy(self, tidy_nc, var_name=None, filter_window=11, trend_windows=(10, 30)):
        v = var_name or self.var_name
        ds = xr.open_dataset(tidy_nc)
        try:
            if v not in ds.data_vars:
                raise KeyError(f"'{v}' not in {tidy_nc}. Found: {list(ds.data_vars)}")
            years = ds["year"].values if "year" in ds.coords else ds["time"].values
            vals = ds[v].values
            print(
                f"[SUMMARY] {v}: years {int(np.nanmin(years))}–{int(np.nanmax(years))} "
                f"(N={len(years)}). NaNs={int(np.isnan(vals).sum())}"
            )
            extras = []
            if filter_window and filter_window > 1:
                extras.append(f"filtered_{int(filter_window)}yr")
            for W in (trend_windows or ()):
                W = int(W)
                extras.extend([f"sliding_trend_{W}yr", f"sliding_trend_{W}yr_per_decade"])
            for k in extras:
                if k in ds:
                    vv = ds[k].values
                    print(f"[SUMMARY] {v}::{k}: shape={vv.shape}, NaNs={int(np.isnan(vv).sum())}")
        finally:
            ds.close()

    def build_tidy_from_combined_annual(
        self,
        annual_nc: str,
        out_nc: str,
        var_name: str,
        filter_window: int = 11,
        trend_windows: tuple = (10, 30),
        add_per_decade: bool = True,
        force: bool = False,
        preserve_year_axis: bool = True,  # --- FIX 2: keep consistent year axis across variables
    ):
        """
        Build a tidy NetCDF (time=0..N-1, year coord, optional filtered/trends)
        from a SINGLE combined annual diagnostics file.

        Key behavior:
          - preserve_year_axis=True (default): keep all years; values may include NaNs.
            This avoids per-variable year dropping that can desynchronize axes.
          - preserve_year_axis=False: drop samples where either year or value is non-finite
            (legacy behavior).
        """
        if (not force) and os.path.exists(out_nc):
            try:
                if os.path.getmtime(out_nc) >= os.path.getmtime(annual_nc):
                    return out_nc
            except FileNotFoundError:
                pass

        ds = xr.open_dataset(annual_nc, decode_times=False, engine=self.xr_engine)
        try:
            if "year" not in ds.coords and "year" not in ds.variables:
                raise KeyError(f"'year' coord not found in {annual_nc}. Found coords: {list(ds.coords)}")
            if var_name not in ds.data_vars:
                raise KeyError(f"'{var_name}' not in {annual_nc}. Found vars: {list(ds.data_vars)}")

            years = np.asarray(ds["year"].values, dtype=float)
            vals  = np.asarray(ds[var_name].values, dtype=float)

            if preserve_year_axis:
                # Keep year axis intact; just require years to be finite.
                good_year = np.isfinite(years)
                years = years[good_year].astype(int, copy=False)
                vals  = vals[good_year]
            else:
                # Legacy: drop years where vals is non-finite.
                good = np.isfinite(years) & np.isfinite(vals)
                years = years[good].astype(int, copy=False)
                vals  = vals[good]

            n = len(years)
            if n == 0:
                raise ValueError(f"No samples for '{var_name}' in {annual_nc} after screening.")

            model_year = np.arange(n, dtype=int)
            rgn = np.zeros(n, dtype=int)

            units = ds[var_name].attrs.get("units", self.default_units.get(var_name, ""))

            data_vars = {var_name: (("time",), vals, {"units": units})}

            if filter_window and filter_window > 1:
                data_vars[f"filtered_{int(filter_window)}yr"] = (("time",), self._centered_filter(vals, filter_window))

            for W in (trend_windows or ()):
                W = int(W)
                tname = f"sliding_trend_{W}yr"
                trend = self._sliding_trend(vals, years, W)
                data_vars[tname] = (("time",), trend)
                if add_per_decade:
                    data_vars[f"{tname}_per_decade"] = (("time",), trend * 10.0)

            ds_out = xr.Dataset(
                data_vars=data_vars,
                coords={
                    "time": ("time", np.arange(n, dtype=int)),
                    "year": ("time", years),
                    "model_year": ("time", model_year),
                    "rgn": ("time", rgn),
                },
                attrs={
                    "source_annual": os.path.basename(annual_nc),
                    "filter_window": int(filter_window) if filter_window else 0,
                    "trend_windows": ",".join(str(int(w)) for w in (trend_windows or ())),
                    "units": units,
                    "combined_len": int(n),
                    "combined_year_min": int(np.nanmin(years)) if n else -1,
                    "combined_year_max": int(np.nanmax(years)) if n else -1,
                    "preserve_year_axis": int(bool(preserve_year_axis)),
                },
            )

            os.makedirs(os.path.dirname(out_nc) or ".", exist_ok=True)
            enc = {k: {"zlib": True, "complevel": 4} for k in ds_out.data_vars}
            ds_out.to_netcdf(out_nc, engine=self.xr_engine, encoding=enc)
            return out_nc
        finally:
            ds.close()


class AMOCDiagnosticsBuilder:
    """
    Fast AMOC diagnostics:
      - stream over mocTimeSeries_*-*.nc (no open_mfdataset)
      - aggregate to annual means
      - write compact per-source NetCDFs + a tidy file for plotting
      - build a combined tidy file: spin-up + piControl (+ optional diagnostics)
    """

    def __init__(
        self,
        base_dir,
        spinup_name,
        picontrol_name,
        subdir="post/time_series/moc",
        moc_file_pattern=r"mocTimeSeries_(\d{6})-(\d{6})\.nc$",
        xr_engine="h5netcdf",
        xr_chunks=None,
        var_name="mocAtlantic26",
        force_reprocess = False, 
    ):
        self.base_dir = base_dir
        self.spinup_name = spinup_name
        self.picontrol_name = picontrol_name
        self.subdir = subdir
        self.xr_engine = xr_engine
        self.xr_chunks = xr_chunks
        self.var_name = var_name

        self._MOC_FILE_RE = re.compile(moc_file_pattern)
        self.default_units = {var_name: "Sv"}
        self.force_reprocess = force_reprocess
        
    # ---------- utils ----------
    def _is_fresh(self, target, sources):
        if not os.path.exists(target):
            return False
        t_mtime = os.path.getmtime(target)
        for s in sources:
            try:
                if os.path.getmtime(s) > t_mtime:
                    return False
            except FileNotFoundError:
                return False
        return True

    def _is_fresh_multi(self, target, sources):
        if not os.path.exists(target):
            return False
        t_mtime = os.path.getmtime(target)
        for s in sources:
            if not os.path.exists(s) or os.path.getmtime(s) > t_mtime:
                return False
        return True

    def _extract_stamps(self, fname):
        m = self._MOC_FILE_RE.search(os.path.basename(fname))
        if not m:
            return None, None
        return int(m.group(1)), int(m.group(2))

    def _find_moc_files(self, case_root):
        files = glob.glob(os.path.join(case_root, self.subdir, "mocTimeSeries_*-*.nc"))
        stamped = [(f, *self._extract_stamps(f)) for f in files]
        stamped = [t for t in stamped if t[1] is not None]
        stamped.sort(key=lambda t: t[1])
        return [t[0] for t in stamped]

    def _get_time_name(self, ds):
        if "time" in ds.dims:
            return "time"
        if "Time" in ds.dims:
            return "Time"

        def looks_like_time(n):
            n = n.lower()
            return ("time" in n) or ("xtime" in n)

        for name in list(ds.coords) + list(ds.variables):
            da = ds[name]
            if da.ndim == 1 and looks_like_time(name):
                return da.dims[0]
        for name in list(ds.coords) + list(ds.variables):
            da = ds[name]
            if da.ndim == 1:
                return da.dims[0]
        raise KeyError("No recognizable time dimension in dataset.")

    def _compute_year_coord(self, ds, tname):
        v = ds[tname]
        # string-like (including char[time, str_len])
        if (v.dtype.kind in ("S", "U")) or ("str" in str(v.dtype).lower()):
            arr = v.values
            if arr.ndim == 2:
                strings = []
                for row in arr:
                    s = "".join([(c.decode("utf-8") if isinstance(c, (bytes, np.bytes_)) else str(c)) for c in row]).strip()
                    strings.append(s)
            else:
                strings = [(s.decode("utf-8") if isinstance(s, (bytes, np.bytes_)) else str(s)) for s in v.values]
            years = [int(s[:4]) if len(s) >= 4 and s[:4].isdigit() else np.nan for s in strings]
            return np.array(years, dtype=float)

        units = (v.attrs.get("units") or "").lower()
        if "months since" in units:
            m = re.search(r"months since\s+(\d{1,4})", units)
            base_year = int(m.group(1)) if m else 1
            months = np.asarray(v.values, dtype=float)
            return np.floor(base_year + months / 12.0).astype(int)
        if "days since" in units:
            m = re.search(r"days since\s+(\d{1,4})", units)
            base_year = int(m.group(1)) if m else 1
            days = np.asarray(v.values, dtype=float)
            return np.floor(base_year + days / 365.0).astype(int)
        if np.issubdtype(v.dtype, np.datetime64):
            return v.dt.year.values

        n = v.sizes[v.dims[0]]
        return np.arange(n)

    # ---------- streaming aggregation ----------
    def _aggregate_annual_mean_stream(self, files):
        if not files:
            raise FileNotFoundError("No MOC time series files found.")

        year_all, val_all = [], []
        vname = self.var_name

        for f in files:
            ds = xr.open_dataset(f, engine=self.xr_engine, decode_times=False)
            try:
                tname = self._get_time_name(ds)
                years = self._compute_year_coord(ds, tname)
                da = ds[vname]
                vals = da.values
                if vals.ndim != 1:
                    # bring time to front then mean over remaining dims
                    time_ax = da.get_axis_num(tname) if tname in da.dims else 0
                    if time_ax != 0:
                        vals = np.moveaxis(vals, time_ax, 0)
                    while vals.ndim > 1:
                        vals = np.nanmean(vals, axis=-1)
                year_all.append(years.astype(np.int64, copy=False))
                val_all.append(np.asarray(vals, dtype=float))
            finally:
                ds.close()

        years = np.concatenate(year_all)
        vals = np.concatenate(val_all)
        good = np.isfinite(vals) & np.isfinite(years)
        years = years[good]
        vals = vals[good]
        if years.size == 0:
            raise ValueError("All annual values are NaN after screening.")

        y0, y1 = int(years.min()), int(years.max())
        bins = y1 - y0 + 1
        idx = years - y0
        ssum = np.bincount(idx, weights=vals, minlength=bins)
        scount = np.bincount(idx, minlength=bins)
        with np.errstate(invalid="ignore", divide="ignore"):
            annual = ssum / scount

        year_coord = np.arange(y0, y1 + 1, dtype=int)
        ds_out = xr.Dataset(
            data_vars={vname: (("year",), annual, {"units": self.default_units.get(vname, "Sv")})},
            coords={"year": ("year", year_coord)},
            attrs={},
        )
        return ds_out

    # ---------- public API: per-source ----------
    def build_sources(self, out_dir, suffix="annual"):
        """
        Returns dict with paths to annual NetCDFs for spinup & piControl.
        """
        os.makedirs(out_dir, exist_ok=True)
        spinup_root = os.path.join(self.base_dir, self.spinup_name)
        
        pctl_root = os.path.join(self.base_dir, self.picontrol_name)
        
        spinup_files = self._find_moc_files(spinup_root)
        
        pctl_files = self._find_moc_files(pctl_root)
        
        spinup_nc = os.path.join(out_dir, f"{self.spinup_name}_amoc_timeseries_{suffix}.nc")
        pctl_nc = os.path.join(out_dir, f"{self.picontrol_name}_amoc_timeseries_{suffix}.nc")
        vname = self.var_name
        
        if self.force_reprocess:
            for _f in (spinup_nc, pctl_nc):
                if os.path.exists(_f):
                    os.remove(_f)
            
        if not self._is_fresh(spinup_nc, spinup_files):
            ds_s = self._aggregate_annual_mean_stream(spinup_files)
            ds_s.to_netcdf(spinup_nc, engine=self.xr_engine, encoding={vname: {"zlib": True, "complevel": 4}})
        if not self._is_fresh(pctl_nc, pctl_files):
            ds_p = self._aggregate_annual_mean_stream(pctl_files)
            ds_p.to_netcdf(pctl_nc, engine=self.xr_engine, encoding={vname: {"zlib": True, "complevel": 4}})

        return {"spinup": spinup_nc, "piControl": pctl_nc}

    # ---------- diagnostics helpers ----------
    def _centered_filter(self, y, window=11):
        w = max(1, int(window))
        if w % 2 == 0:
            w += 1
        pad = w // 2
        y = np.asarray(y, dtype=float)
        isfin = np.isfinite(y).astype(float)
        y_filled = np.where(np.isfinite(y), y, 0.0)
        k = np.ones(w, dtype=float)
        num = np.convolve(y_filled, k, mode="same")
        den = np.convolve(isfin, k, mode="same")
        out = np.full_like(y, np.nan, dtype=float)
        good = den > 0
        out[good] = num[good] / den[good]
        out[:pad] = np.nan
        out[-pad:] = np.nan
        return out

    def _sliding_trend(self, y, x, window):
        n = len(x)
        w = max(2, int(window))
        if w % 2 == 0:
            w += 1
        half = w // 2
        out = np.full(n, np.nan, dtype=float)
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        for i in range(n):
            i0 = max(0, i - half)
            i1 = min(n, i + half + 1)
            xx = x[i0:i1]
            yy = y[i0:i1]
            m = np.isfinite(xx) & np.isfinite(yy)
            if m.sum() >= 3:
                coef = np.polyfit(xx[m], yy[m], 1)  # units per year
                out[i] = coef[0]
        return out

    def _load_year_val(self, ncpath):
        """Load (years, values) for self.var_name from a per-source annual file."""
        vname = self.var_name  # <-- fix: ensure correct variable used
        ds = xr.open_dataset(ncpath, decode_times=False, engine=self.xr_engine)
        try:
            if "year" in ds.coords:
                years = ds["year"].values
                vals = ds[vname].values
            else:
                tname = self._get_time_name(ds)
                years = self._compute_year_coord(ds, tname)
                vals = ds[vname].values
            vals = np.asarray(vals, dtype=float)
            years = np.asarray(years, dtype=int)
            n = min(len(years), len(vals))
            return years[:n], vals[:n]
        finally:
            ds.close()
            
    # Decide how to construct the 'year' coordinate for the combined series
    def make_years(self, sequence, y_spin, y_pi):
        if sequence == "calendar":
            return np.concatenate([y_spin, y_pi])

        if sequence == "append":
            # Always continue sequentially after spin-up,  e.g., 1..2000 + 2001..2500
            offset = int(np.nanmax(y_spin))
            y_pi_seq = (np.arange(len(y_pi)) + offset + 1).astype(int)
            return np.concatenate([y_spin.astype(int, copy=False), y_pi_seq])

        # sequence == "auto": detect overlap/reset
        y_spin = np.asarray(y_spin, dtype=int)
        y_pi   = np.asarray(y_pi, dtype=int)
        # Heuristics: if piControl doesn't clearly continue after spin-up,
        # or if there's large overlap / reset-to-1, we "append".
        overlap = (y_pi.min() <= y_spin.max())
        reset_like = (np.median(y_pi) <= np.median(y_spin)) or (y_pi.min() <= 5)
        continues = (y_pi.min() > y_spin.max()) and (y_pi.min() - y_spin.max() <= 5)

        if overlap or reset_like:
            offset = int(np.nanmax(y_spin))
            y_pi_seq = (np.arange(len(y_pi)) + offset + 1).astype(int)
            return np.concatenate([y_spin, y_pi_seq])
        elif continues:
            return np.concatenate([y_spin, y_pi])
        else:
            # Fallback: if unsure, append to guarantee contiguity
            offset = int(np.nanmax(y_spin))
            y_pi_seq = (np.arange(len(y_pi)) + offset + 1).astype(int)
            return np.concatenate([y_spin, y_pi_seq])

    # ---------- combined tidy (spin-up + piControl) ----------
    def build_combined_tidy(
        self,
        spinup_nc,
        pctl_nc,
        combined_tidy_nc,
        filter_window=11,
        trend_windows=(10,),
        add_per_decade=True,
        sequence="auto",   # "auto" | "calendar" | "append"
        force=False,       # NEW: bypass freshness and rebuild
        strict=True,       # NEW: raise if combined length != spinup+piControl after dedup
    ):
        """
        Combine spin-up and piControl into one tidy file.

        sequence:
          - "auto": detect overlap/reset and append piControl after spin-up by offsetting years if needed
          - "calendar": keep 'year' exactly as encoded in files
          - "append": force piControl to follow spin-up contiguously in 'year'

        If strict=True, raise an error when the final sample count != sum of inputs (after intended rules).
        """
        print(spinup_nc)
        print(pctl_nc)
        deps = [spinup_nc, pctl_nc]
        if (not force) and self._is_fresh_multi(combined_tidy_nc, deps):
            return combined_tidy_nc

        vname = self.var_name

        # Load sources
        y_s, v_s = self._load_year_val(spinup_nc)
        y_p, v_p = self._load_year_val(pctl_nc)
        n_s, n_p = len(y_s), len(y_p)

        # Construct combined 'year'
        years_combined = self.make_years(sequence, y_s, y_p).astype(int, copy=False)
        values_combined = np.concatenate([v_s, v_p])

        # Defensive: drop only *adjacent* duplicate years (shouldn't happen in "append")
        if years_combined.size > 1:
            keep = np.ones_like(years_combined, dtype=bool)
            for i in range(len(years_combined) - 1):
                if years_combined[i] == years_combined[i + 1]:
                    keep[i] = False    # keep the last duplicate
            if not np.all(keep):
                years_combined  = years_combined[keep]
                values_combined = values_combined[keep]

        # Sanity check expected size
        n_comb = len(years_combined)
        expected = n_s + n_p
        if strict and n_comb != expected:
            raise ValueError(
                f"[combine] Expected {expected} samples (spinup={n_s} + piControl={n_p}) "
                f"but got {n_comb}. sequence='{sequence}'. "
                f"Spinup years {int(y_s.min())}-{int(y_s.max())}, "
                f"piControl years {int(y_p.min())}-{int(y_p.max())}."
            )

        # Build dataset
        model_year = np.arange(n_comb, dtype=int)
        rgn = np.zeros(n_comb, dtype=int)

        data_vars = {vname: (("time",), values_combined, {"units": self.default_units.get(vname, "Sv")})}

        if filter_window and filter_window > 1:
            data_vars[f"filtered_{int(filter_window)}yr"] = (("time",), self._centered_filter(values_combined, filter_window))
        for W in trend_windows or ():
            tname = f"sliding_trend_{int(W)}yr"
            trend = self._sliding_trend(values_combined, years_combined, W)
            data_vars[tname] = (("time",), trend)
            if add_per_decade:
                data_vars[f"{tname}_per_decade"] = (("time",), trend * 10.0)

        ds_out = xr.Dataset(
            data_vars=data_vars,
            coords={
                "time": ("time", np.arange(n_comb, dtype=int)),
                "year": ("time", years_combined),
                "model_year": ("time", model_year),
                "rgn": ("time", rgn),
            },
            attrs={
                "source_series": f"spinup + piControl (sequence={sequence})",
                "filter_window": int(filter_window) if filter_window else 0,
                "trend_windows": ",".join(str(int(w)) for w in (trend_windows or [])),
                "units": self.default_units.get(vname, "Sv"),
                # helpful provenance:
                "spinup_len": int(n_s),
                "picontrol_len": int(n_p),
                "combined_len": int(n_comb),
                "spinup_year_min": int(y_s.min()) if n_s else -1,
                "spinup_year_max": int(y_s.max()) if n_s else -1,
                "picontrol_year_min": int(y_p.min()) if n_p else -1,
                "picontrol_year_max": int(y_p.max()) if n_p else -1,
                "combined_year_min": int(years_combined.min()) if n_comb else -1,
                "combined_year_max": int(years_combined.max()) if n_comb else -1,
            },
        )
        enc = {k: {"zlib": True, "complevel": 4} for k in ds_out.data_vars}
        ds_out.to_netcdf(combined_tidy_nc, engine=self.xr_engine, encoding=enc)
        return combined_tidy_nc


class AtmosphereBalanceAnalysis:
    def __init__(
            self, 
            base_dir, 
            spinup_name, 
            picontrol_name, 
            component='atm', 
            region_map=None,
            var_dict=None,
            subdir="post/time_series/atm", 
            xr_engine='h5netcdf', 
            xr_chunks=None
        ):
        self.base_dir = base_dir
        self._xr_engine = xr_engine
        self._xr_chunks = xr_chunks
        self.spinup_name = spinup_name
        self.picontrol_name = picontrol_name
        self.component = component
        self.subdir = subdir
        self.dataset = None
        self.region_map = self._get_region(region_map) 
        self.var_dict =  self._get_var(var_dict)
        self.seconds_per_year = 365 * 24 * 60 * 60  # unit conversion from second to year 

    def _get_region(self, region_map):
        if region_map is None:
            region_map = {0: "Global", 1: "Northern Hemisphere", 2: "Southern Hemisphere"}
        return region_map
        
    def _get_var(self, var_dict):
        if  var_dict is None:
            var_dict = {
                "net_toa_flux_restom":             {"name": "RESTOM",   'label': "RESTOM",  'unit': 'W m$^{-2}$'},
                "net_sfc_flux_ressurf":            {"name": "RESSURF",  'label': "RESSURF", 'unit': 'W m$^{-2}$'},
                "global_surface_skin_temperature": {"name": "TS",       'label': "TS",      'unit': '$^{o}$C'},
                "global_surface_air_temperature":  {"name": "TAS",      'label': "TAS",     'unit': '$^{o}$C'},
                "toa_radiation":                   {"name": "RESTOA",   'label': "RESTOA",  'unit': 'W m$^{-2}$'},
                "net_atm_energy_imbalance":        {"name": "RESNET",   'label': "RESNET",  'unit': 'W m$^{-2}$'},
                "net_atm_water_imbalance":         {"name": "MASSFLUX", 'label': "E - P",   'unit': 'mm yr$^{-1}$'},
            }
        return var_dict
        
    def _load_lhflx(self, path):
        Lv = 2.501e6
        Lf = 3.337e5
        needed_vars = ['QFLX', 'PRECC', 'PRECL', 'PRECSC', 'PRECSL']
        ds_dict = {}
        for var in needed_vars:
            files = sorted(glob.glob(os.path.join(path, f"{var}_*.nc")))
            if not files:
                continue
            with xr.open_mfdataset(
                files,
                engine=self._xr_engine,
                chunks=self._xr_chunks,
                combine='by_coords',
                decode_times=True
            ) as ds_var:
                ds_dict[var] = ds_var[var].load().copy()
        if all(v in ds_dict for v in needed_vars):
            LHFLX = (Lv + Lf) * ds_dict['QFLX'] - Lf * 1.0e3 * (
                ds_dict['PRECC'] + ds_dict['PRECL'] - ds_dict['PRECSC'] - ds_dict['PRECSL']
            )
        elif 'QFLX' in ds_dict:
            LHFLX = Lv * ds_dict['QFLX']
        else:
            raise ValueError("Missing QFLX for LHFLX derivation.")
        LHFLX.name = "LHFLX"
        LHFLX.attrs["units"] = "W m-2"
        LHFLX.attrs["long_name"] = "Latent Heat Fluxes Inferred from QFLX"
        return LHFLX

    def _load_prec(self, path):
        needed_vars = ['PRECC', 'PRECL']
        ds_dict = {}
        for var in needed_vars:
            files = sorted(glob.glob(os.path.join(path, f"{var}_*.nc")))
            if not files:
                continue
            with xr.open_mfdataset(
                files,
                engine=self._xr_engine,
                chunks=self._xr_chunks,
                combine='by_coords',
                decode_times=True
            ) as ds_var:
                ds_dict[var] = ds_var[var].load().copy()
        if all(v in ds_dict for v in needed_vars):
            PRECT = 1.0e3 * (ds_dict['PRECC'] + ds_dict['PRECL'])
        else:
            raise ValueError("Missing PRECC/PRECL for PRECT derivation.")
        PRECT.name = "PRECT"
        PRECT.attrs["units"] = "kg m-2 s-1"
        PRECT.attrs["long_name"] = "Total Precipitation Rate"
        return PRECT

    def _load_precs(self, path):
        needed_vars = ['PRECSC', 'PRECSL']
        ds_dict = {}
        for var in needed_vars:
            files = sorted(glob.glob(os.path.join(path, f"{var}_*.nc")))
            if not files:
                continue
            with xr.open_mfdataset(
                files,
                engine=self._xr_engine,
                chunks=self._xr_chunks,
                combine='by_coords',
                decode_times=True
            ) as ds_var:
                ds_dict[var] = ds_var[var].load().copy()
        if all(v in ds_dict for v in needed_vars):
            PRECST = 1.0e3 * (ds_dict['PRECSC'] + ds_dict['PRECSL'])
        else:
            raise ValueError("Missing PRECSC/PRECSL for PRECST derivation.")
        PRECST.name = "PRECST"
        PRECST.attrs["units"] = "kg m-2 s-1"
        PRECST.attrs["long_name"] = "Total Snow Precipitation Rate"
        return PRECST

    def _load_flux_dataset(self, run_name):
        path = os.path.join(self.base_dir,run_name,self.subdir)
        needed_vars = ['TREFHT', 'TS', 'FSNTOA', 'FLUT', 'FSNT', 'FLNT', 'FSNS', 'FLNS', 'SHFLX', 'QFLX']
        data_vars = {}
        for var in needed_vars:
            files = sorted(glob.glob(os.path.join(path, f"{var}_*.nc")))
            if not files:
                raise FileNotFoundError(f"[{var}] not found for {run_name}")
            # Use context manager to avoid NetCDF close errors
            with xr.open_mfdataset(
                files,
                engine=self._xr_engine,
                chunks=self._xr_chunks,
                combine='by_coords',
                decode_times=True
            ) as ds_var:
                data_vars[var] = ds_var[var].load()

        data_vars["LHFLX"] = self._load_lhflx(path)
        data_vars["PRECT"] = self._load_prec(path)
        data_vars["PRECST"] = self._load_precs(path)
        return xr.Dataset(data_vars)

    def _load_combined_dataset(self):
        ds_spinup = self._load_flux_dataset(self.spinup_name)
        ds_picontrol = self._load_flux_dataset(self.picontrol_name)
        n_spinup = ds_spinup.sizes["time"]
        n_picontrol = ds_picontrol.sizes["time"]
        spinup_time = xr.cftime_range(
            start=cftime.DatetimeNoLeap(1, 1, 1), 
            periods=n_spinup, 
            freq="MS"
        )
        start_year = 1 + (n_spinup // 12)
        picontrol_time = xr.cftime_range(
            start=cftime.DatetimeNoLeap(start_year, 1, 1), 
            periods=n_picontrol, 
            freq="MS"
        )
        ds_spinup = ds_spinup.assign_coords(time=("time", spinup_time))
        ds_picontrol = ds_picontrol.assign_coords(time=("time", picontrol_time))
        ds_combined = xr.concat([ds_spinup, ds_picontrol], dim="time")
        self.dataset = ds_combined
        return ds_combined

    def compute_fluxes(self):
        if self.dataset is None:
            self._load_combined_dataset()
        ds = self.dataset
        
        restom = ds['FSNT'] - ds['FLNT']
        restoa = ds['FSNTOA'] - ds['FLUT']
        ressurf = ds['FSNS'] - ds['FLNS'] - ds['SHFLX'] - ds['LHFLX']
        radnet = restoa - ressurf
        for var in [restom, restoa, ressurf, radnet]:
            var.attrs["units"] = "W m-2"
            
        massflux = (ds['QFLX'] - ds['PRECT']) * self.seconds_per_year
        massflux.attrs['units'] = "mm yr-1"
        massflux.attrs["long_name"] = "Water Mass Fluxes(E-P)"
        
        tas = ds['TREFHT'] - 273.15 
        tas.attrs["units"] = "degC"
        tas.attrs["long_name"] = "Surface Air Temperature (2m)"
        
        sst = ds['TS'] - 273.15 
        sst.attrs["units"] = "degC"
        sst.attrs["long_name"] = "Surface Skin Temperature (TS)"
        
        self.dataset = xr.merge([
            restom.rename("RESTOM"),
            restoa.rename("RESTOA"),
            ressurf.rename("RESSURF"),
            radnet.rename("RESNET"),
            massflux.rename("MASSFLUX"),
            tas.rename("TAS"),
            sst.rename("SST"),
            ds['FSNT'],
            ds['FLNT'],
            ds['FSNTOA'],
            ds['FLUT'],
            ds['FSNS'],
            ds['FLNS'],
            ds['SHFLX'],
            ds['LHFLX'],
            ds['PRECT'],
            ds['PRECST'],
            ds['QFLX'].rename("EVAP"),
        ])
        return self.dataset

    def to_timeseries_df(
            self, experiment, 
            output_file, 
            annual=True, 
            force_process=False
        ):
        
        if os.path.exists(output_file) and not force_process:
            return 
            
        if force_process and os.path.exists(output_file):
            os.remove(output_file)
            
        flux_ds = self.compute_fluxes()
        if "time" not in flux_ds.coords or flux_ds["time"].isnull().all():
            raise ValueError("Invalid or missing time coordinate.")
        if annual:
            years = xr.DataArray([t.year for t in flux_ds["time"].values], coords={"time": flux_ds["time"]}, dims="time")
            flux_ds.coords["year"] = years
            ds_annual = flux_ds.groupby("year").mean(dim="time", keep_attrs=True)
            ds_annual.to_netcdf(output_file) 
        else:
            flux_ds.to_netcdf(output_file)
        return
