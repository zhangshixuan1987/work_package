from __future__ import annotations
"""Regional / depth-resolved ocean time-series analyses comparing experiments
(MPAS-Analysis mpasTimeSeriesOcean output): ENSO (Nino3.4) and OHC Hovmoller / lag correlation.

  OCNTimeSeriesPlotter (shared base):
    NinoSSTPlotter       (1_enso_compare_analysis)
    NinoIndexPlotter     (1_enso_compare_analysis)
    OHCHovmollerPlotter  (17_ohc_corr_compare_analysis)
  OCNEnsoSpectrum        (1_enso_compare_analysis)
  LaggedOHCCorrelation   (17_ohc_corr_compare_analysis, 18_ohc_corr_full_period_analysis)
"""
from typing import Optional, Tuple, List
import os, re, glob, json
import numpy as np
import xarray as xr
import cftime
from scipy import signal
import matplotlib.pyplot as plt
from matplotlib.ticker import (
    FixedLocator, FixedFormatter, AutoMinorLocator,
    MultipleLocator, NullFormatter
)
from matplotlib.lines import Line2D
import statsmodels.api as sm
from scipy import stats
import os, glob, re
from matplotlib.ticker import FixedLocator, FixedFormatter, MultipleLocator, NullFormatter
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as mticker
from typing import Union, Optional
from datetime import datetime, timezone
import os, glob, re, json
import pandas as pd
from scipy import signal, stats
from scipy.signal import windows
from typing import Dict, Sequence, Optional
from datetime import datetime


class OCNTimeSeriesPlotter:
    """Shared reader for MPAS-Analysis regional time series (mpasTimeSeriesOcean / moc). Use a subclass."""

    def __init__(
        self,
        root_dir,
        sub_dir,
        file_key="mpasTimeSeriesOcean",
        var_name=None,
        compute_anomaly=False,
        annual_mean=False,
        running_mean=False,
        running_mean_months=0,         # 0/None disables; e.g., 12 for 12-mo
        running_mean_center="center",  # 'center' | 'trailing' | 'leading'
        running_mean_min_frac=0.5,     # require ≥ this fraction of valid samples
        region_dim=None,          # e.g., "nOceanRegions" (auto-detect if None)
        region_index=None,        # select by integer index
        region_name=None,         # or select by human-readable name
        region_reduce=None,       # None | "mean" | "median"
        region_names_var="regionNames",  # common mpas-analysis string var
        scale=1.0,
        baseline_years=1,
        baseline_months=12,
        filt_cutoff_months=24,
        filt_order=4,
        trend_years=30,
        trend_unit='per_year',  # or 'per_decade'
        year_min=1,
        year_max=350,
        use_filter=False,
        filt_method='pad',
        filt_padlen=None,
        filt_padtype='odd',
        edge_fix_points=0,
        use_mfdataset=True,
        chunks=None,
        engine=None,
        time_origin=None,
        calendar="noleap"
    ):
        # I/O
        self.time_origin = time_origin   # e.g., "0001-01-01 00:00:00"
        self.calendar = calendar         # "noleap" typical for MPAS/E3SM
        self.root_dir = root_dir
        self.sub_dir = sub_dir
        self.file_key = file_key
        self.var_name = var_name
        if self.var_name is None:
            raise ValueError("var_name must be provided.")

        self.baseline_months = int(baseline_months)
        self.baseline_years = (int(baseline_years)
                               if baseline_years is not None
                               else max(1, int(np.ceil(self.baseline_months / 12.0))))

        self.region_dim       = region_dim
        self.region_index     = None if region_index is None else int(region_index)
        self.region_name      = None if region_name  is None else str(region_name)
        self.region_reduce    = None if region_reduce is None else str(region_reduce).lower()
        self.region_names_var = None if region_names_var in (None, "None") else str(region_names_var)

        # Running mean flags
        self.annual_mean = bool(annual_mean)
        self.running_mean = bool(running_mean)
        self.compute_anomaly = bool(compute_anomaly)
        self.running_mean_months = 0 if running_mean_months in (None, 0) else int(running_mean_months)
        self.running_mean_center = str(running_mean_center).lower()
        self.running_mean_min_frac = float(running_mean_min_frac)

        self.scale = float(scale)
        self.filt_cutoff_months = float(filt_cutoff_months)
        self.filt_order = int(filt_order)
        self.trend_years = int(trend_years)
        self.trend_unit = str(trend_unit)
        self.year_min = int(year_min)
        self.year_max = int(year_max)
        self.use_filter = bool(use_filter)
        self.filt_method   = str(filt_method)
        self.filt_padlen   = None if filt_padlen is None else int(filt_padlen)
        self.filt_padtype  = str(filt_padtype)
        self.edge_fix_points = int(edge_fix_points)

        # xarray options
        self._xr_use_mfdataset = bool(use_mfdataset)
        self._xr_chunks = chunks
        self._xr_engine = engine

        # Filter design if needed
        if self.use_filter:
            fc = 1.0 / self.filt_cutoff_months
            wn = max(min(fc / 0.5, 0.999999), 1e-6)
            self._butter_b, self._butter_a = signal.butter(self.filt_order, wn, btype="low")
        else:
            self._butter_b = self._butter_a = None

        # Storage
        self._meta, self._raw, self._filt, self._stats = [], [], [], []
        self._monthly = []   # NEW: raw monthly series (pre-annual/running/anomaly), used for Niño index overlay

    def _infer_region_dim(self, da, tdim):
        """Return the non-time region-like dim, or None if 1-D."""
        if da.ndim <= 1:
            return None
        # prefer user-specified
        if self.region_dim and self.region_dim in da.dims:
            return self.region_dim
        # otherwise pick any non-time dim that looks like a region
        for d in da.dims:
            if d != tdim and ("region" in d.lower() or "nregion" in d.lower()):
                return d
        # fallback: the first non-time dim
        for d in da.dims:
            if d != tdim:
                return d
        return None

    def _decode_region_names(self, ds, region_dim):
        """
        Try to pull human-readable region names.
        - Prefer a 2D char array like ds['regionNames'] (nOceanRegions, StrLen).
        - Otherwise, look for a string/bytes coordinate on region_dim.
        Returns list[str] or None.
        """
        # option 1: explicit regionNames var (mpas-analysis)
        if self.region_names_var in ds:
            arr = ds[self.region_names_var].values
            # char array -> join along last axis
            if arr.ndim == 2:
                names = ["".join(row.astype(str)).strip().strip("\x00").replace("\x00", "") for row in arr]
                return [n if n != "" else f"{region_dim}[{i}]" for i, n in enumerate(names)]
            # already string or object
            if arr.ndim == 1:
                out = []
                for i, v in enumerate(arr):
                    if isinstance(v, bytes):
                        s = v.decode("utf-8", "ignore").strip()
                    else:
                        s = str(v).strip()
                    out.append(s if s else f"{region_dim}[{i}]")
                return out

        # option 2: coordinate on region_dim
        if region_dim in ds.coords:
            arr = ds.coords[region_dim].values
            if arr.ndim == 1 and arr.size > 0:
                out = []
                for i, v in enumerate(arr):
                    if isinstance(v, bytes):
                        s = v.decode("utf-8", "ignore").strip()
                    else:
                        s = str(v).strip()
                    out.append(s if s else f"{region_dim}[{i}]")
                # guard against numeric-only coords
                if any(ch.isalpha() for ch in "".join(out)):
                    return out
        return None

    def _region_index_from_name(self, ds, region_dim, target_name):
        """Case-insensitive match for region name; returns index or raises."""
        names = self._decode_region_names(ds, region_dim)
        if not names:
            raise KeyError(
                f"Region names not found (looked for '{self.region_names_var or 'regionNames'}' "
                f"or a string coord on '{region_dim}')."
            )
        tgt = str(target_name).strip().lower()
        # exact case-insensitive first
        for i, n in enumerate(names):
            if n.lower() == tgt:
                return i

        # relaxed match: remove spaces/underscores
        def norm(s):
            return "".join(ch for ch in s.lower() if ch.isalnum())

        tgt2 = norm(tgt)
        for i, n in enumerate(names):
            if norm(n) == tgt2:
                return i
        raise KeyError(f"Region name '{target_name}' not found. Available: {names}")

    def _nan_running_mean(self, x, win, center="center", min_valid_frac=0.5):
        """
        NaN-aware running average over window `win`.
        center: 'center' (symmetric), 'trailing' (up to t), 'leading' (from t).
        Keeps array length; returns NaN where valid-count < min_valid_frac*win.
        """
        x = np.asarray(x, float)
        n = len(x)
        if win <= 1 or n == 0:
            return x.copy()

        valid = np.isfinite(x).astype(float)
        x0 = np.where(np.isfinite(x), x, 0.0)
        k = np.ones(int(win), float)

        if center == "trailing":
            s = np.convolve(x0, k, mode="full")[win - 1:win - 1 + n]
            c = np.convolve(valid, k, mode="full")[win - 1:win - 1 + n]
        elif center == "leading":
            s = np.convolve(x0, k, mode="full")[:n]
            c = np.convolve(valid, k, mode="full")[:n]
        else:  # 'center'
            s = np.convolve(x0, k, mode="same")
            c = np.convolve(valid, k, mode="same")

        with np.errstate(invalid="ignore", divide="ignore"):
            y = s / c
        min_count = max(1, int(np.ceil(min_valid_frac * win)))
        y[c < min_count] = np.nan
        return y

    def _rolling_nanstd(self, x, win):
        """
        Rolling nanstd with a centered window length `win` (odd enforced).
        Shrinks at edges. Returns same-length array.
        """
        x = np.asarray(x, float)
        n = len(x)
        if n == 0:
            return x.copy()
        win = max(1, int(win))
        if win % 2 == 0:
            win += 1
        half = win // 2
        out = np.empty(n, dtype=float)
        for i in range(n):
            a = max(0, i - half)
            b = min(n, i + half + 1)
            out[i] = np.nanstd(x[a:b], ddof=1) if (b - a) > 1 else np.nanstd(x[a:b])
        return out

    def to_annual_mean(self, series):
        """Convert monthly 1-D series to annual mean."""
        n = len(series)
        n_years = n // 12
        arr = series[:n_years * 12].reshape(n_years, 12)
        return arr.mean(axis=1)

    def _to_anomaly_annual(self, x_annual, baseline_years=None):
        """Annual anomaly using mean of the first N years."""
        if baseline_years is None:
            baseline_years = self.baseline_years
        if not (0 < baseline_years <= len(x_annual)):
            raise ValueError(f"baseline_years={baseline_years} must be within [1, {len(x_annual)}]")
        clim = float(np.nanmean(np.asarray(x_annual, float)[:baseline_years]))
        return np.asarray(x_annual, float) - clim

    def _to_anomaly(self, x):
        """Monthly anomaly using mean of the first baseline_months."""
        if not (0 < self.baseline_months <= len(x)):
            raise ValueError(f"baseline_months={self.baseline_months} must be within [1, {len(x)}]")
        clim = float(np.nanmean(x[: self.baseline_months]))
        return x - clim

    def _lowpass(self, x):
        if not self.use_filter:
            return x
        x = np.asarray(x, float)
        n = len(x)

        # (re)design filter if needed
        if self._butter_b is None or self._butter_a is None:
            fc = 1.0 / self.filt_cutoff_months
            wn = max(min(fc / 0.5, 0.999999), 1e-6)
            self._butter_b, self._butter_a = signal.butter(self.filt_order, wn, btype="low")
        b, a = self._butter_b, self._butter_a

        if self.filt_method == "gust":
            y = signal.filtfilt(b, a, x, method="gust")
        else:
            default_pad = 3 * (max(len(a), len(b)) - 1)
            cutoff_pad = int(max(1, round(2 * self.filt_cutoff_months)))
            padlen = self.filt_padlen if self.filt_padlen is not None else max(default_pad, cutoff_pad)
            padlen = min(padlen, max(1, n - 1))
            if padlen >= n:
                y = signal.filtfilt(b, a, x, method="gust")
            else:
                y = signal.filtfilt(b, a, x, method="pad", padtype=self.filt_padtype, padlen=padlen)

        k = self.edge_fix_points
        if k > 0 and n >= 2 * k + 1:
            y[:k] = y[k]
            y[-k:] = y[-k - 1]
        return y

    @staticmethod
    def _parse_time_origin(s, calendar="noleap"):
        if s is None:
            return None
        s = str(s).strip()
        try:
            y = int(s[0:4]); m = int(s[5:7]); d = int(s[8:10])
            hh = int(s[11:13]) if len(s) >= 13 else 0
            mm = int(s[14:16]) if len(s) >= 16 else 0
            ss = int(s[17:19]) if len(s) >= 19 else 0
        except Exception:
            return None
        if calendar.lower() in ("noleap", "365_day", "365-day"):
            return cftime.DatetimeNoLeap(y, m, d, hh, mm, ss)
        return cftime.DatetimeNoLeap(y, m, d, hh, mm, ss)

    def _construct_monthly_time(self, n):
        if self.time_origin is None:
            return None
        t0 = self._parse_time_origin(self.time_origin, self.calendar)
        if t0 is None:
            return None
        return xr.cftime_range(start=t0, periods=int(n), freq="MS")

    def _guess_time_dim_da(self, da):
        for d in da.dims:
            if d.lower().startswith("time"):
                return d
        return da.dims[0]

    def _find_mpas_timeseries_files(self, root_dir, exp, sub_dir, file_key):
        base = os.path.join(root_dir, exp, sub_dir)
        print(base)
        single = glob.glob(os.path.join(base, f"{file_key}.nc"))
        if single:
            return single
        segs = glob.glob(os.path.join(base, f"{file_key}_*.nc"))
        if not segs:
            return []
        def start_stamp(path):
            m = re.search(rf"{re.escape(file_key)}_(\d{{6}})-(\d{{6}})\.nc$", os.path.basename(path))
            return int(m.group(1)) if m else 0
        segs.sort(key=start_stamp)
        return segs

    @staticmethod
    def _find_time_slice(year_token):
        a, b = [int(t.strip()) for t in year_token.split("-")]
        if a > b:
            a, b = b, a
        m_start = (a - 1) * 12
        m_end = b * 12
        return m_start, m_end, (b - a + 1) * 12

    def _open_dataset(self, files):
        xr_opts = dict(engine=self._xr_engine, chunks=self._xr_chunks, decode_times=True)
        if len(files) == 1 or not self._xr_use_mfdataset:
            return xr.open_dataset(files[0], **xr_opts)
        try:
            return xr.open_mfdataset(files, combine="by_coords", parallel=True, **xr_opts)
        except Exception:
            dsets = [xr.open_dataset(fp, **xr_opts) for fp in files]
            probe = dsets[0][self.var_name] if self.var_name in dsets[0] else None
            if probe is None:
                probe = list(dsets[0].data_vars.values())[0]
            tdim = self._guess_time_dim_da(probe)
            return xr.concat(dsets, dim=tdim)

    @staticmethod
    def _decode_time_from_Time_var(ds):
        if "Time" not in ds:
            return None
        da = ds["Time"]
        units = da.attrs.get("units", None)
        cal = da.attrs.get("calendar", "noleap")
        if units is None:
            return None
        times = xr.coding.times.decode_cf_datetime(da.values, units, calendar=cal, use_cftime=True)
        return np.asarray(times)

    @staticmethod
    def _decode_time_from_xtime_chars(ds, name="xtime_startMonthly"):
        if name not in ds:
            return None
        a = ds[name].values  # (Time, StrLen)
        s = ["".join(row.astype(str)).strip().strip("\x00").replace("\x00", "") for row in a]
        out = []
        for txt in s:
            txt = txt.replace("T", "_").replace(" ", "_")
            try:
                y = int(txt[0:4]); m = int(txt[5:7]); d = int(txt[8:10])
                hh = int(txt[11:13]); mm = int(txt[14:16]); ss = int(txt[17:19])
                out.append(cftime.DatetimeNoLeap(y, m, d, hh, mm, ss))
            except Exception:
                try:
                    y = int(txt[0:4]); m = int(txt[5:7]); d = int(txt[8:10])
                    out.append(cftime.DatetimeNoLeap(y, m, d))
                except Exception:
                    return None
        return np.asarray(out)

    @staticmethod
    def _slice_indices_from_year_token(time_index, year_token):
        a, b = [int(t.strip()) for t in year_token.split("-")]
        years = np.array([t.year for t in time_index], dtype=int)
        mask = (years >= a) & (years <= b)
        if not mask.any():
            return None
        idx = np.nonzero(mask)[0]
        return int(idx[0]), int(idx[-1] + 1)

    def _derive_ssh(
        self,
        ds,
        var="ssh_anom",
        keep_baseline=True,
        base_var_name="ssh_proxy",
        time_dim=None,
        baseline_mode="calendar_year",
        first_N=12
    ):
        v_area = "timeMonthly_avg_areaCellGlobal"
        v_vol = "timeMonthly_avg_volumeCellGlobal"
        for vname in (v_area, v_vol):
            if vname not in ds:
                raise KeyError(f"Variable '{vname}' not found in dataset.")

        if time_dim is None:
            candidates = ["time", "Time", "nMonths", "nTime", "TimeMonthly", "timeMonthly", "nt"]
            time_dim = next((d for d in candidates if d in ds.dims), None)

        if time_dim is None:
            raise ValueError(
                "Could not find a time-like dimension. Pass time_dim=... explicitly (e.g., 'Time') or check ds.dims."
            )

        ssh = ds[v_vol] / ds[v_area]
        if time_dim not in ssh.dims and time_dim in ds.dims:
            ssh = ssh.expand_dims({time_dim: ds.sizes[time_dim]}) if time_dim in ds.dims else ssh

        if keep_baseline:
            ssh = ssh.rename(base_var_name)

        has_dt_year = False
        year_coord = None
        if time_dim in ds.coords:
            try:
                year_coord = ds[time_dim].dt.year
                has_dt_year = True
            except Exception:
                has_dt_year = False

        if baseline_mode == "calendar_year" and has_dt_year:
            first_year = int(year_coord.min().values)
            baseline = ssh.where(year_coord == first_year).mean(time_dim, keep_attrs=True)
            ref_desc = f"calendar year {first_year}"
        else:
            if ds.sizes.get(time_dim, 0) < first_N:
                first_N = ds.sizes.get(time_dim, 0)
            baseline = ssh.isel({time_dim: slice(0, first_N)}).mean(time_dim, keep_attrs=True)
            ref_desc = f"first {first_N} {time_dim} steps"

        ssh_anom = ssh - baseline
        ssh_anom.attrs["units"] = "m"
        ssh_anom.attrs["long_name"] = "SSH Proxy Anomaly"
        ssh_anom.attrs["description"] = (
            "Sea Surface Height proxy computed as volume/area (column thickness). "
            f"Anomaly computed relative to the {ref_desc} per grid cell."
        )

        if keep_baseline:
            ssh.attrs["units"] = "m"
            ssh.attrs["long_name"] = "Sea Surface Height Proxy (column thickness)"
            ssh.attrs["description"] = "Computed as volume/area."

            baseline = baseline.rename(f"{base_var_name}_baseline")
            baseline.attrs["units"] = "m"
            baseline.attrs["long_name"] = f"Baseline of {base_var_name}"
            baseline.attrs["description"] = f"Mean over the {ref_desc} per grid cell."

            ds[base_var_name] = ssh
            ds[f"{base_var_name}_baseline"] = baseline

        ds[var] = ssh_anom
        return ds

    def _load_series(self, exp, year_token, sub_dir_override=None):
        subdir = sub_dir_override if sub_dir_override is not None else self.sub_dir
        files = self._find_mpas_timeseries_files(self.root_dir, exp, subdir, self.file_key)
        print(files)

        if not files:
            search_root = os.path.join(self.root_dir, exp)
            candidates = glob.glob(os.path.join(
                search_root, "post", "analysis", "mpas_analysis", "ts_*_climo_*", "**", f"{self.file_key}*.nc"
            ), recursive=True)
            if candidates:
                hit = candidates[0]
                subdir = os.path.relpath(os.path.dirname(hit), start=os.path.join(self.root_dir, exp))
                files = self._find_mpas_timeseries_files(self.root_dir, exp, subdir, self.file_key)

        if not files:
            raise FileNotFoundError(
                f"No '{self.file_key}.nc' (or segmented) found for '{exp}' under subdir='{subdir}'."
            )

        ds = self._open_dataset(files)
        if self.var_name == "SSH":
            ds = self._derive_ssh(ds, var="SSH", base_var_name="ssh_proxy")
        elif self.var_name not in ds:
            raise KeyError(f"Variable '{self.var_name}' not found in {files[0]}")

        tdim = self._guess_time_dim_da(ds[self.var_name])
        n = ds[self.var_name].sizes[tdim]

        time_idx = self._decode_time_from_Time_var(ds)
        if time_idx is None:
            ti_start = None; ti_end = None
            for nm in ("xtime_startMonthly", "startTime"):
                ti_start = self._decode_time_from_xtime_chars(ds, nm)
                if ti_start is not None and len(ti_start) == n:
                    break
            for nm in ("xtime_endMonthly", "endTime"):
                ti_end = self._decode_time_from_xtime_chars(ds, nm)
                if ti_end is not None and len(ti_end) == n:
                    break
            time_idx = ti_start if (ti_start is not None and len(ti_start) == n) \
                else (ti_end if (ti_end is not None and len(ti_end) == n) else None)

        if time_idx is None:
            time_idx = self._construct_monthly_time(n)

        use_start = use_end = None
        if time_idx is not None and len(time_idx) == n:
            se = self._slice_indices_from_year_token(time_idx, year_token)
            if se is not None:
                use_start, use_end = se

        if use_start is None or use_end is None:
            m_start, m_end, _ = self._find_time_slice(year_token)
            if m_start >= n or m_end > n or (m_end - m_start) <= 0:
                use_start, use_end = 0, n
            else:
                use_start, use_end = m_start, m_end

        da = ds[self.var_name].isel({tdim: slice(use_start, use_end)}).astype("float64")

        region_dim = self._infer_region_dim(da, tdim)
        if region_dim is not None and region_dim in da.dims:
            if self.region_reduce in ("mean", "median"):
                if self.region_reduce == "mean":
                    da = da.mean(dim=region_dim, skipna=True)
                else:
                    da = da.median(dim=region_dim, skipna=True)
            else:
                if self.region_name is not None:
                    idx = self._region_index_from_name(ds, region_dim, self.region_name)
                    da = da.isel({region_dim: idx})
                elif self.region_index is not None:
                    da = da.isel({region_dim: int(self.region_index)})
                else:
                    da = da.isel({region_dim: 0})

        while da.ndim > 1:
            for d in list(da.dims):
                if d != tdim and da.sizes[d] == 1:
                    da = da.squeeze(d)
                    break
            else:
                for d in list(da.dims):
                    if d != tdim:
                        da = da.mean(dim=d, skipna=True)
                        break

        series = np.asarray(da.values, dtype=float).reshape(-1)
        if not np.isfinite(series).any():
            raise ValueError(f"All values are NaN for exp '{exp}', var '{self.var_name}'")
        return series * self.scale

    def _last_n_year_stats(
        self,
        raw_series,
        filt_series,
        n_years,
        trend_unit=None,
        add_trend_ci=True,
        ci_method="hac",
        hac_lag_years=1,
        glsar_max_iter=10
    ):
        import numpy as np

        if trend_unit is None:
            trend_unit = getattr(self, "trend_unit", "per_year")

        samples_per_year = 1 if getattr(self, "annual_mean", False) else 12
        n_last = max(1, int(n_years * samples_per_year))

        if len(raw_series) < n_last or len(filt_series) < n_last:
            raise ValueError(f"Series too short to compute {n_years}-year stats")

        recent_raw = np.asarray(raw_series[-n_last:], dtype=float)
        mean_val = float(np.nanmean(recent_raw))
        min_val  = float(np.nanmin(recent_raw))
        max_val  = float(np.nanmax(recent_raw))
        std_val  = float(np.nanstd(recent_raw))

        y = np.asarray(filt_series[-n_last:], dtype=float)
        x = np.arange(len(y), dtype=float)

        m, _b = np.polyfit(x, y, 1)
        trend_per_year = m * samples_per_year
        scale = 10.0 if trend_unit == 'per_decade' else 1.0
        trend_out = trend_per_year * scale

        result = dict(mean=mean_val, min=min_val, max=max_val, std=std_val, trend=trend_out)

        if not add_trend_ci:
            return result

        try:
            t_years = x / samples_per_year
            X = sm.add_constant(t_years)

            if ci_method == "hac":
                maxlags = int(max(1, hac_lag_years * samples_per_year))
                ols = sm.OLS(y, X, missing='drop').fit(cov_type='HAC', cov_kwds={'maxlags': maxlags})
                slope_per_year = float(ols.params[1])
                se_per_year    = float(ols.bse[1])
                dof = int(ols.df_resid)
                tcrit = stats.t.ppf(0.975, dof)
                lo, hi = slope_per_year - tcrit * se_per_year, slope_per_year + tcrit * se_per_year
                method = f"HAC(Newey–West, maxlags={maxlags})"

            elif ci_method == "glsar":
                glsar = sm.GLSAR(y, X, rho=1)
                res = glsar.iterative_fit(glsar_max_iter)
                slope_per_year = float(res.params[1])
                se_per_year    = float(res.bse[1])
                dof = int(res.df_resid)
                tcrit = stats.t.ppf(0.975, dof)
                lo, hi = slope_per_year - tcrit * se_per_year, slope_per_year + tcrit * se_per_year
                method = "GLSAR(AR1)"
            else:
                raise ImportError("Force fallback to Neff")

            result.update({
                "trend_se": se_per_year * scale,
                "trend_ci95": (lo * scale, hi * scale),
                "trend_dof": dof,
                "trend_method": method
            })
            return result

        except Exception:
            try:
                t_years = x / samples_per_year
                X = sm.add_constant(t_years)
                ols = sm.OLS(y, X, missing='drop').fit()
                e = ols.resid - np.nanmean(ols.resid)
                if len(e) < 3:
                    raise RuntimeError("Not enough residuals for AR(1) estimate")
                rho = float(np.corrcoef(e[1:], e[:-1])[0, 1])
                N = int(np.isfinite(y).sum())
                Neff = max(3.0, N * (1 - rho) / (1 + rho))
                slope_per_year = float(ols.params[1])
                se_ols = float(ols.bse[1])
                se_adj = se_ols * np.sqrt(N / Neff)
                dof_eff = int(max(3, Neff - 2))
                from scipy import stats as _stats
                tcrit = _stats.t.ppf(0.975, dof_eff)
                lo, hi = slope_per_year - tcrit * se_adj, slope_per_year + tcrit * se_adj

                result.update({
                    "trend_se": se_adj * scale,
                    "trend_ci95": (lo * scale, hi * scale),
                    "trend_dof": dof_eff,
                    "trend_method": f"OLS + AR(1) Neff (rho={rho:.2f}, Neff={Neff:.1f})"
                })
            except Exception:
                result.update({
                    "trend_se": None,
                    "trend_ci95": None,
                    "trend_dof": None,
                    "trend_method": "unavailable"
                })
            return result

    def compute_nino_index_from_monthly(self, sst_monthly, base_years=30, smooth_months=3):
        """
        Compute a Niño-like index from monthly SST series:
          index(t) = SST(t) - monthly_climatology(month(t))
        where monthly_climatology is computed from the first `base_years` years.
        Optionally smooth with `smooth_months` running mean (centered).
        """
        x = np.asarray(sst_monthly, float)
        nb = int(base_years * 12)
        if len(x) < nb:
            raise ValueError(f"Need at least {nb} months to form a {base_years}-yr baseline")

        base = x[:nb]
        clim12 = np.array([np.nanmean(base[m::12]) for m in range(12)], dtype=float)
        month_idx = np.arange(len(x)) % 12
        anom = x - clim12[month_idx]

        if smooth_months and smooth_months > 1:
            anom = self._nan_running_mean(
                anom, int(smooth_months),
                center="center",
                min_valid_frac=self.running_mean_min_frac
            )
        return anom

    def overlay_nino_index(
        self,
        ax,
        colors=None,
        base_years=30,
        smooth_months=3,
        ylabel="Niño3.4 index (°C)",
        linestyle="--",
        linewidth=2.5,
        alpha=0.9,
        add_zero_line=True,
        fontsize=25
    ):
        """
        Overlay the Niño index on a twin y-axis, using stored raw monthly series.
        Call this AFTER plotter.plot(...).
        """
        ax2 = ax.twinx()
        ax2.set_ylabel(ylabel, fontsize=fontsize * 0.9)

        palette = colors or ["black", "#cc0000", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]
        for i, _meta in enumerate(self._meta):
            col = palette[i % len(palette)]
            if i >= len(self._monthly):
                continue
            idx = self.compute_nino_index_from_monthly(
                self._monthly[i],
                base_years=base_years,
                smooth_months=smooth_months
            )
            ax2.plot(idx, color=col, linewidth=linewidth, alpha=alpha, linestyle=linestyle)

        if add_zero_line:
            ax2.axhline(0.0, ls=":", lw=1.2, color="k", alpha=0.6)

        ax2.tick_params(labelsize=fontsize * 0.9, length=6, width=1)
        return ax2

    def add_experiment(self, exp, label=None,
                       year_token="0001-0350",
                       color=None, linewidth=4.0, alpha=1.0,
                       sub_dir=None):
        """Register an experiment and precompute processed series + stats."""
        series = self._load_series(exp, year_token, sub_dir_override=sub_dir)

        # NEW: store raw monthly series for Niño index overlay
        self._monthly.append(np.asarray(series, float).copy())

        if self.annual_mean:
            series_ann = self.to_annual_mean(series)
            series_use = self._to_anomaly_annual(series_ann) if self.compute_anomaly else series_ann
        elif self.running_mean and (self.running_mean_months and self.running_mean_months > 0):
            base = self._nan_running_mean(
                series, self.running_mean_months,
                center=self.running_mean_center,
                min_valid_frac=self.running_mean_min_frac
            )
            series_use = self._to_anomaly(base) if self.compute_anomaly else base
        else:
            base = series
            series_use = self._to_anomaly(base) if self.compute_anomaly else base

        if self.use_filter:
            series_filt = self._lowpass(series_use)
            stats = self._last_n_year_stats(series_use, series_filt,
                                            self.trend_years, self.trend_unit)
            self._filt.append(series_filt)
        else:
            stats = self._last_n_year_stats(series_use, series_use,
                                            self.trend_years, self.trend_unit)
            self._filt.append(None)

        self._meta.append({
            "exp": exp,
            "label": label or exp,
            "year_token": year_token,
            "color": color,
            "linewidth": linewidth,
            "alpha": alpha,
        })
        self._raw.append(series_use)
        self._stats.append(stats)

    def _set_time_axis(
        self, ax, is_annual: bool, step_yrs: int,
        xlim_years: Optional[Tuple[int, int]],
        ref_yrs: Optional[List[int]],
        fontsize: int,
        *, labels_zero_based: bool = False
    ):
        if is_annual:
            ny = self.year_max - self.year_min + 1
            ax.set_xlim(0, ny - 1)
            step = max(1, int(step_yrs))
            xt = np.arange(0, ny, step)
            if labels_zero_based:
                labels = [str(int(t)) for t in xt]
            else:
                labels = [str(int(self.year_min + t)) for t in xt]
            ax.xaxis.set_major_locator(FixedLocator(xt))
            ax.xaxis.set_major_formatter(FixedFormatter(labels))
            ax.xaxis.set_minor_locator(AutoMinorLocator(n=4))
            if ref_yrs:
                for y in ref_yrs:
                    xpos = int(y) - self.year_min
                    if 0 <= xpos < ny:
                        ax.axvline(xpos, color="k", ls="--", lw=1.2, alpha=0.6)
            if xlim_years:
                a, b = xlim_years
                ax.set_xlim(a - self.year_min, b - self.year_min)
        else:
            nm_cfg = (self.year_max - self.year_min + 1) * 12
            ax.set_xlim(0, nm_cfg)
            step = max(1, int(step_yrs))
            tick_step_m = 12 * step
            xt = np.arange(0, nm_cfg + 1, tick_step_m)
            if labels_zero_based:
                labels = [str(int(t / 12)) for t in xt]
            else:
                labels = [str(int(self.year_min + t / 12.0)) for t in xt]
            ax.xaxis.set_major_locator(FixedLocator(xt))
            ax.xaxis.set_major_formatter(FixedFormatter(labels))
            ax.xaxis.set_minor_locator(AutoMinorLocator(n=4))
            if ref_yrs:
                for y in ref_yrs:
                    xpos = (int(y) - self.year_min) * 12
                    if 0 <= xpos <= nm_cfg:
                        ax.axvline(xpos, color="k", ls="--", lw=1.2, alpha=0.6)
            if xlim_years:
                a, b = xlim_years
                ax.set_xlim((a - self.year_min) * 12, (b - self.year_min) * 12)

    def _build_legend(
        self,
        stats,
        ax,
        vunit: str,
        fontsize: int,
        legend_loc: str,
        colors: Optional[List[str]] = None,
        note: Optional[str] = None,
        *,
        compact: bool = False,
        ncol: int = 3,
        bbox: Optional[Tuple[float, float]] = None,
        fontsize_scale: float = 0.9,
        framealpha: float = 0.6,
        borderaxespad: float = 0.0,
        ci_format: str = "pm",
    ):
        from matplotlib.lines import Line2D

        palette = colors or ["black", "#cc0000", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]
        lines, labels = [], []

        for i, meta in enumerate(self._meta):
            col = palette[i % len(palette)]
            lines.append(Line2D([0], [0], color=col, linewidth=meta.get("linewidth", 2.0)))
            labels.append(meta["label"] if compact else f"{meta['label']}")

        if note:
            lines.append(Line2D([0], [0], color="grey", linewidth=3.0))
            labels.append(note)

        return ax.legend(
            lines, labels, fontsize=int(fontsize * fontsize_scale),
            loc=legend_loc, bbox_to_anchor=bbox, ncol=ncol,
            frameon=False, framealpha=framealpha,
            handlelength=2.0, handletextpad=0.4, labelspacing=0.3,
            borderaxespad=borderaxespad, borderpad=0.8,
        )

    def save_nino_plot_data(
        self,
        out_nc: str,
        *,
        base_years: int = 30,
        smooth_months: int = 3,
        window: int = 60,
        ref_yrs=None,
        xlim_years=None,
        labels_zero_based: bool = False,
    ):
        """
        Save plot-ready Niño index + rolling std envelope for each experiment to NetCDF.
    
        Variables:
          - nino_index(exp, time)
          - env_std(exp, time)    # rolling std (>=0)
          - model_year(time)      # year_min + time/12
        Coords:
          - exp (id), label
          - time (month index since start of slice)
        """
        n = len(self._meta)
        if n == 0 or len(self._monthly) == 0:
            raise RuntimeError("No experiments loaded. Call add_experiment(...) before saving.")
    
        idx_list = []
        Lmax = 0
    
        # compute indices for all experiments
        for i in range(n):
            if i >= len(self._monthly) or self._monthly[i] is None:
                idx = np.array([], dtype=float)
            else:
                idx = self.compute_nino_index_from_monthly(
                    self._monthly[i],
                    base_years=base_years,
                    smooth_months=smooth_months
                )
            idx = np.asarray(idx, float)
            idx_list.append(idx)
            Lmax = max(Lmax, len(idx))
    
        if Lmax == 0:
            raise RuntimeError("Computed Niño index length is 0 for all experiments.")
    
        def _pad(x, L, fill=np.nan):
            y = np.full((L,), fill, dtype=float)
            if len(x) > 0:
                y[:len(x)] = x
            return y
    
        idx_pad = np.vstack([_pad(x, Lmax) for x in idx_list])
    
        # rolling envelope (std) for each experiment (same window you use for plotting)
        env_list = []
        for i in range(n):
            idx = idx_list[i]
            env = self._rolling_nanstd(idx, win=window) if len(idx) else np.array([], dtype=float)
            env_list.append(np.asarray(env, float))
        env_pad = np.vstack([_pad(x, Lmax) for x in env_list])
    
        exp_ids    = [m.get("exp", f"exp{i}") for i, m in enumerate(self._meta)]
        exp_labels = [m.get("label", m.get("exp", f"exp{i}")) for i, m in enumerate(self._meta)]
        t = np.arange(Lmax, dtype=int)
    
        ds = xr.Dataset(
            data_vars=dict(
                nino_index=(("exp", "time"), idx_pad),
                env_std=(("exp", "time"), env_pad),
                model_year=("time", self.year_min + t / 12.0),
            ),
            coords=dict(
                exp=("exp", exp_ids),
                label=("exp", exp_labels),
                time=("time", t),
            ),
        )
    
        ds.attrs["plot_settings_json"] = json.dumps(dict(
            base_years=int(base_years),
            smooth_months=int(smooth_months),
            window=int(window),
            ref_yrs=ref_yrs,
            xlim_years=xlim_years,
            labels_zero_based=bool(labels_zero_based),
            year_min=int(self.year_min),
            year_max=int(self.year_max),
            calendar=str(self.calendar),
            time_origin=str(self.time_origin),
        ))
    
        ds.to_netcdf(out_nc)
        return out_nc

    @staticmethod
    def load_nino_plot_data(nc_path: str) -> xr.Dataset:
        return xr.open_dataset(nc_path)

    def _shade_segments(self, ax, segments, color, hatch="///", alpha=0.4, zorder=-50):
        """Hatch (start_year, end_year) windows (model years, inclusive) on ax."""
        for y0, y1 in segments:
            y0, y1 = sorted((int(y0), int(y1)))
            if self.annual_mean:
                x0, x1 = y0 - self.year_min, y1 - self.year_min + 1
            else:
                x0, x1 = (y0 - self.year_min) * 12, (y1 - self.year_min + 1) * 12
            # rasterized: PDF hatch patterns have a fixed 1-inch tile that does not shrink
            # when the figure is scaled in LaTeX (coarse, thick hatching); a bitmap at the
            # savefig dpi scales with the figure like the PNG preview.
            ax.axvspan(x0, x1, facecolor="none", edgecolor=color, hatch=hatch,
                       lw=0.0, alpha=alpha, zorder=zorder, rasterized=True)

    def _sha_for(self, exp_shas, i):
        """Hatching spec for experiment i from exp_shas keyed by index, exp name or label."""
        if not exp_shas:
            return None
        meta = self._meta[i]
        for key in (i, meta.get("exp"), meta.get("label")):
            if key in exp_shas:
                return exp_shas[key]
        return None


class NinoIndexPlotter(OCNTimeSeriesPlotter):
    """Nino3.4 index panels, ENSO event counts, PCC."""

    def _pcc(self, a, b):
        a = np.asarray(a, float)
        b = np.asarray(b, float)
        m = np.isfinite(a) & np.isfinite(b)
        if m.sum() < 3:
            return np.nan
        return float(np.corrcoef(a[m], b[m])[0, 1])

    def _split_month_index(self, ref_year):
        return int((int(ref_year) - int(self.year_min)) * 12)

    def _lowfreq(self, x, win=60):
        return self._nan_running_mean(x, win, center="center")

    def _count_enso_events(self, idx, thr=0.5, min_duration=5, gap=1):
        """
        Count El Niño / La Niña events in a 1D monthly index series.
    
        Definition:
          - El Niño: idx >= +thr for >= min_duration consecutive months
          - La Niña: idx <= -thr for >= min_duration consecutive months
        'gap' merges events separated by <= gap months below threshold.
        Returns: (nino_count, nina_count)
        """
        x = np.asarray(idx, float)
        m = np.isfinite(x)
        x = np.where(m, x, np.nan)
    
        def _count_for_mask(mask):
            # mask: boolean array where condition is met (above/below threshold)
            mask = np.asarray(mask, bool)
    
            # merge short gaps (<= gap) between true segments
            if gap and gap > 0:
                g = int(gap)
                n = len(mask)
                i = 0
                while i < n:
                    if mask[i]:
                        i += 1
                        continue
                    # start of a False run
                    j = i
                    while j < n and (not mask[j]):
                        j += 1
                    # False run is [i, j)
                    if (j - i) <= g:
                        # if it is sandwiched between True segments, fill it
                        left_true = (i - 1) >= 0 and mask[i - 1]
                        right_true = j < n and mask[j]
                        if left_true and right_true:
                            mask[i:j] = True
                    i = j
    
            # count segments with length >= min_duration
            cnt = 0
            n = len(mask)
            i = 0
            while i < n:
                if not mask[i]:
                    i += 1
                    continue
                j = i
                while j < n and mask[j]:
                    j += 1
                if (j - i) >= int(min_duration):
                    cnt += 1
                i = j
            return cnt
    
        nino = _count_for_mask(x >= +thr)
        nina = _count_for_mask(x <= -thr)
        return nino, nina

    def plot_nino_index(
        self,
        title: str = "Niño3.4 index",
        y_label: str = "Niño3.4 index (°C)",
        ylim: Optional[Tuple[float, float]] = None,
        x_label: str = "Model Time (year)",
        colors: Optional[List[str]] = None,
        figsize: Tuple[int, int] = (24, 10),
        fontsize: int = 25,
        step_yrs: int = 10,
        xlim_years: Optional[Tuple[int, int]] = None,
        ref_yrs: Optional[List[int]] = None,
        labels_zero_based: bool = False,
        outname: Optional[str] = None,
        legend_loc: str = "upper left",
        legend_note: Optional[str] = None,
        # index definition
        base_years: int = 30,
        smooth_months: int = 3,
        linestyle: str = "-",
        linewidth: float = 3.0,
        alpha: float = 1.0,
        add_zero_line: bool = True,
        # NEW
        separate_panels: bool = True,
        wspace=0.25,
        hspace=0.35,
        window=60,   
        thr = 0.5,
        min_dur = 5,   # months (classic ONI-ish persistence)
        gap = 1,       # merge events separated by 1 month
        exp_shas: Optional[dict] = None,
        sha_alpha: float = 0.4,
        axes=None,
    ):
        """
        Plot Niño3.4 index.

        If separate_panels=True: one subplot per experiment (stacked, sharex).
        If False: all experiments on one axis (previous behavior).
        exp_shas: hatched year windows per experiment, same format as NinoSSTPlotter.plot
                  (keyed by experiment index, exp name or label; color defaults to the
                  experiment's line color). With separate panels, each experiment's windows
                  are hatched in its own panel only.
        axes: with separate_panels=True, existing axes (one per experiment, top to bottom) to
              draw into instead of creating a figure, e.g. to combine with other panels.
        """
        palette = colors or ["black", "#cc0000", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]
        n = len(self._meta)
        # choose a distinct but neutral color
        sigma_color = "#555555"   # dark gray (recommended)
        # alternatives:
        # sigma_color = "#7f7f7f"  # lighter gray
        # sigma_color = "#4c72b0"  # muted blue

        if not separate_panels:
            # ---- single-axis (old behavior) ----
            fig, ax = plt.subplots(figsize=figsize)
            ymins, ymaxs = [], []

            for i, meta in enumerate(self._meta):
                col = (colors[i] if colors is not None else meta.get("color")) or palette[i % len(palette)]
                if i >= len(self._monthly):
                    continue

                idx = self.compute_nino_index_from_monthly(
                    self._monthly[i], base_years=base_years, smooth_months=smooth_months
                )

                sha = self._sha_for(exp_shas, i)
                if sha and sha.get("segments"):
                    self._shade_segments(ax, sha["segments"], color=sha.get("color", col),
                                         hatch=sha.get("hatch", "///"), alpha=sha.get("alpha", sha_alpha))

                ax.plot(idx, color=col, linewidth=linewidth, alpha=alpha, linestyle=linestyle)
                if idx.size:
                    ymins.append(np.nanmin(idx)); ymaxs.append(np.nanmax(idx))

            self._set_time_axis(
                ax, is_annual=False, step_yrs=step_yrs,
                xlim_years=xlim_years, ref_yrs=ref_yrs,
                fontsize=fontsize, labels_zero_based=labels_zero_based
            )

            ax.set_title(title, fontsize=fontsize * 0.9, loc="left", pad=12)
            ax.set_xlabel(x_label, fontsize=fontsize * 0.9)
            ax.set_ylabel(y_label, fontsize=fontsize * 0.9)

            if add_zero_line:
                ax.axhline(0.0, ls="--", lw=1.2, color="k", alpha=0.8)

            if ylim:
                ax.set_ylim(*ylim)
            elif ymins and ymaxs:
                lo, hi = float(np.nanmin(ymins)), float(np.nanmax(ymaxs))
                pad = 0.08 * max(1e-12, hi - lo)
                ax.set_ylim(lo - pad, hi + pad)

            ax.grid(True, alpha=0.25)
            ax.tick_params(labelsize=fontsize * 0.9, length=6, width=1)

            self._build_legend(
                self._stats, ax, vunit="°C", fontsize=fontsize,
                legend_loc=legend_loc, colors=(colors or palette),
                note=legend_note, compact=True
            )

            if outname:
                fig.savefig(outname, dpi=600, bbox_inches="tight")
            return fig, ax

        # ---- separate panels: one experiment per axis ----
        if axes is None:
            fig, axes = plt.subplots(
                nrows=n, ncols=1,
                figsize=figsize, sharex=True,
                gridspec_kw={"hspace": hspace}
            )
            if n == 1:
                axes = [axes]
        else:
            axes = list(np.ravel(axes))
            if len(axes) != n:
                raise ValueError(f"got {len(axes)} axes for {n} experiments")
            fig = axes[0].figure

        # ---- compute indices for all experiments first ----
        global_min, global_max = np.inf, -np.inf
        idx_all = []

        for i, _meta in enumerate(self._meta):
            if i >= len(self._monthly):
                idx_all.append(None)
                continue
            idx = self.compute_nino_index_from_monthly(
                self._monthly[i], base_years=base_years, smooth_months=smooth_months
            )
            idx_all.append(idx)
            if idx.size and np.isfinite(idx).any():
                global_min = min(global_min, float(np.nanmin(idx)))
                global_max = max(global_max, float(np.nanmax(idx)))

        # ---------------- similarity metrics vs Full-CPL (pre/post ref year) ----------------
        pcc_lf_pre  = [np.nan] * n
        pcc_lf_post = [np.nan] * n
        var_ratio_pre  = [np.nan] * n
        var_ratio_post = [np.nan] * n
        
        nino_pre  = [np.nan] * n
        nino_post = [np.nan] * n
        nina_pre  = [np.nan] * n
        nina_post = [np.nan] * n
        
        # Full-CPL assumed to be the first experiment
        idx_ref = idx_all[0] if (len(idx_all) > 0) else None

        # split point: first ref year
        split_m = self._split_month_index(ref_yrs[0]+1) if (ref_yrs and len(ref_yrs) > 0) else None
        
        # low-frequency reference
        lf_ref = self._lowfreq(idx_ref, win=window) if (idx_ref is not None) else None

        if (idx_ref is not None) and (lf_ref is not None) and (split_m is not None):
            for j in range(1, n):
                idx_j = idx_all[j]
                if idx_j is None:
                    continue

                lf_j = self._lowfreq(idx_j, win=window)

                L = min(len(idx_ref), len(idx_j), len(lf_ref), len(lf_j))
                if L < 3:
                    continue
                s = min(split_m, L)

                # Low-frequency PCC
                pcc_lf_pre[j]  = self._pcc(lf_ref[:s], lf_j[:s])
                pcc_lf_post[j] = self._pcc(lf_ref[s:L], lf_j[s:L])

                # Variance ratio (raw index variability)
                sd_ref_pre  = np.nanstd(idx_ref[:s])
                sd_ref_post = np.nanstd(idx_ref[s:L])
                sd_j_pre    = np.nanstd(idx_j[:s])
                sd_j_post   = np.nanstd(idx_j[s:L])

                var_ratio_pre[j]  = (sd_j_pre  / sd_ref_pre)  if (sd_ref_pre  > 0) else np.nan
                var_ratio_post[j] = (sd_j_post / sd_ref_post) if (sd_ref_post > 0) else np.nan
                
                # ENSO event counts (El Niño / La Niña) in each segment
                nino_pre[0], nina_pre[0]   = self._count_enso_events(idx_ref[:s], thr=thr, min_duration=min_dur, gap=gap)
                nino_post[0], nina_post[0] = self._count_enso_events(idx_ref[s:L], thr=thr, min_duration=min_dur, gap=gap)

                nino_pre[j], nina_pre[j] = self._count_enso_events(idx_j[:s], thr=thr, min_duration=min_dur, gap=gap)
                nino_post[j], nina_post[j] = self._count_enso_events(idx_j[s:L], thr=thr, min_duration=min_dur, gap=gap)
                
        # pad for global auto ylim
        if np.isfinite(global_min) and np.isfinite(global_max):
            pad = 0.08 * max(1e-12, global_max - global_min)
            auto_ylim = (global_min - pad, global_max + pad)
        else:
            auto_ylim = None

        # safe right-title
        title_right = title if (title is not None) else ""

        # ---- plot each panel ----
        for i, ax in enumerate(axes):
            meta = self._meta[i]
            idx = idx_all[i]

            # hatched windows of this experiment, behind everything in its panel
            sha = self._sha_for(exp_shas, i)
            if sha and sha.get("segments"):
                col = (colors[i] if colors is not None else meta.get("color")) or palette[i % len(palette)]
                self._shade_segments(ax, sha["segments"], color=sha.get("color", col),
                                     hatch=sha.get("hatch", "///"), alpha=sha.get("alpha", sha_alpha))

            if idx is not None:
                t = np.arange(len(idx))

                # ENSO-style filled anomalies
                pos = np.where(idx >= 0.0, idx, 0.0)
                neg = np.where(idx <  0.0, idx, 0.0)

                ax.fill_between(t, 0.0, pos, color="red",  alpha=0.5, linewidth=0)
                ax.fill_between(t, 0.0, neg, color="blue", alpha=0.5, linewidth=0)
                
                ax.plot(idx, color="k", linewidth=0.5, alpha=0.95, zorder=4)

                # rolling ENSO amplitude envelope (± std) on RIGHT axis
                env = self._rolling_nanstd(idx, win=window)
                
                ax_r = ax.twinx()
                
                ax_r.plot(+env, color=sigma_color, lw=2.0, alpha=0.75, zorder=6)
                ax_r.plot(-env, color=sigma_color, lw=2.0, alpha=0.75, zorder=6)
                
                # symmetric limits so ±env is centered
                ax_r.set_ylim(*ylim)
                
                # symmetric limits
                #if np.isfinite(env).any():
                #    hi = float(np.nanmax(env))
                #    ax_r.set_ylim(-1.15 * hi, 1.15 * hi)
                
                # axis styling to match the line
                ax_r.set_ylabel(
                    f"Rolling σ ({window/12:.0f}yr)",
                    fontsize=fontsize * 0.75,
                    color=sigma_color
                )
                ax_r.tick_params(
                    axis="y",
                    labelsize=fontsize * 0.75,
                    colors=sigma_color,
                    length=5,
                    width=1
                )
                
                # color the right spine explicitly
                ax_r.spines["right"].set_color(sigma_color)
                ax_r.spines["right"].set_linewidth(1.2)
                
                # keep background transparent
                ax_r.patch.set_alpha(0.0)

            # strong zero line (single)
            if add_zero_line:
                ax.axhline(0.0, color="k", lw=1.6, alpha=0.9, zorder=5)

            # ENSO thresholds
            thr = 0.5
            ax.axhline(+thr, color="k", lw=1.0, ls=":", alpha=0.5)
            ax.axhline(-thr, color="k", lw=1.0, ls=":", alpha=0.5)

            ax.grid(True, axis="y", alpha=0.20)

            # titles
            ax.set_title(f"({chr(97+i)}) {meta['label']}", fontsize=fontsize * 0.95, loc="left", pad=6)
            ax.set_title(f"{title_right}", fontsize=fontsize * 0.95, loc="right", pad=6)
            ax.set_ylabel(y_label, fontsize=fontsize * 0.85)

            ax.tick_params(labelsize=fontsize * 0.85, length=6, width=1)

            if ylim:
                ax.set_ylim(*ylim)
            elif auto_ylim:
                ax.set_ylim(*auto_ylim)

            # ref lines
            if ref_yrs:
                nm_cfg = (self.year_max - self.year_min + 1) * 12
                for y in ref_yrs:
                    xpos = (int(y) - self.year_min) * 12
                    if 0 <= xpos <= nm_cfg:
                        ax.axvline(xpos, color="k", ls="--", lw=1.1, alpha=0.5)

            # annotation
            if i != 0:
                txt = (
                    f"El Niño / La Niña: ({int(nino_pre[i])}/{int(nina_pre[i])}, {int(nino_post[i])}/{int(nina_post[i])})\n"
                    f"σ ratio vs {self._meta[0]['label']}: ({var_ratio_pre[i]:.2f}, {var_ratio_post[i]:.2f})"
                )
            else:
                txt = (
                    f"El Niño / La Niña: ({int(nino_pre[i])}/{int(nina_pre[i])}, {int(nino_post[i])}/{int(nina_post[i])})"
                )
                
            ax.text(
                0.02, 0.95, txt,
                transform=ax.transAxes,
                ha="left", va="top",
                fontsize=fontsize * 0.68,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.7),
                zorder=10
            )

        # shared x axis formatting on the bottom axis only
        self._set_time_axis(
            axes[-1],
            is_annual=False,
            step_yrs=step_yrs,
            xlim_years=xlim_years,
            ref_yrs=None,  # already drawn above
            fontsize=fontsize,
            labels_zero_based=labels_zero_based,
        )
        axes[-1].set_xlabel(x_label, fontsize=fontsize * 0.9)

        if outname:
            fig.savefig(outname, dpi=600, bbox_inches="tight")
        return fig, axes


class NinoSSTPlotter(OCNTimeSeriesPlotter):
    """Nino3.4 SST time series with running mean and std band."""

    def overlay_nino_index(
        self,
        ax,
        colors=None,
        base_years=30,
        smooth_months=3,
        ylabel="Niño3.4 index (°C)",
        linestyle="--",
        linewidth=2.5,
        alpha=0.9,
        add_zero_line=True,
        fontsize=25
    ):
        """
        Overlay the Niño index on a twin y-axis, using stored raw monthly series.
        Call this AFTER plotter.plot(...).
        """
        ax2 = ax.twinx()
        ax2.set_ylabel(ylabel, fontsize=fontsize * 0.9)

        palette = colors or ["black", "#cc0000", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]
        for i, meta in enumerate(self._meta):
            col = palette[i % len(palette)]
            if i >= len(self._monthly):
                continue
            idx = self.compute_nino_index_from_monthly(
                self._monthly[i],
                base_years=base_years,
                smooth_months=smooth_months
            )
            ax2.plot(idx, color=col, linewidth=linewidth, alpha=alpha, linestyle=linestyle)

        if add_zero_line:
            ax2.axhline(0.0, ls=":", lw=1.2, color="k", alpha=0.6)

        ax2.tick_params(labelsize=fontsize * 0.9, length=6, width=1)
        return ax2

    def _build_legend(
        self,
        stats,
        ax,
        vunit: str,
        fontsize: int,
        legend_loc: str,
        colors: Optional[List[str]] = None,
        note: Optional[str] = None,
        *,
        compact: bool = False,
        ncol: int = 1,
        bbox: Optional[Tuple[float, float]] = None,
        fontsize_scale: float = 0.9,
        framealpha: float = 0.6,
        borderaxespad: float = 0.0,
        ci_format: str = "pm",
        indices: Optional[List[int]] = None,
    ):
        palette = colors or ["black", "#cc0000", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]
        lines, labels = [], []

        # indices: subset of experiments to list (default: all)
        for i in (range(len(self._meta)) if indices is None else indices):
            meta = self._meta[i]
            col = palette[i % len(palette)]
            lines.append(Line2D([0], [0], color=col, linewidth=meta.get("linewidth", 2.0)))

            if compact:
                labels.append(meta["label"])
                continue
                
            s = stats[i]
            mean_str = f"{s['mean']:.2f} ± {s['std']:.2f} {vunit}"
            vtrend_unit = f"{vunit}/10yr" if getattr(self, "trend_unit", "per_year") == "per_decade" else f"{vunit}/yr"

            # trend text + CI
            tr_txt = "n/a"
            if np.isfinite(s.get("trend", np.nan)):
                if s.get("trend_ci95") is not None:
                    lo, hi = s["trend_ci95"]
                    if ci_format == "range":
                        ci_txt = f"[{lo:.2f}, {hi:.2f}] (95% CI)"
                    else:
                        half = 0.5 * (hi - lo)
                        ci_txt = f"± {half:.2f} (95% CI)"
                    tr_txt = f"{s['trend']:.2f} {ci_txt} {vtrend_unit}"
                elif s.get("trend_se") is not None:
                    # fallback: approximate 95% CI from SE
                    try:
                        from scipy import stats as _st
                        dof = int(s.get("trend_dof") or 30)
                        tcrit = float(_st.t.ppf(0.975, dof))
                    except Exception:
                        tcrit = 1.96
                    half = tcrit * float(s["trend_se"])
                    tr_txt = f"{s['trend']:.2f} ± {half:.2f} (95% CI) {vtrend_unit}"
                else:
                    tr_txt = f"{s['trend']:.2f} {vtrend_unit}"
                    
            labels.append(f"{meta['label']} (mean: {mean_str}, trend: {tr_txt})")
            #labels.append(f"{meta['label']}")

        if note:
            lines.append(Line2D([0], [0], color="grey", linewidth=3.0))
            labels.append(note)

        return ax.legend(
            lines, labels, fontsize=int(fontsize * fontsize_scale),
            loc=legend_loc, bbox_to_anchor=bbox, ncol=ncol,
            frameon=False, framealpha=framealpha,
            handlelength=2.0, handletextpad=0.4, labelspacing=0.3,
            borderaxespad=borderaxespad, borderpad=0.8,
        )

    def plot(
        self,
        vname: str = "SST",
        vunit: str = "°C",
        ylim: Optional[Tuple[float, float]] = None,
        title: Optional[str] = None,
        y_label: Optional[str] = None,
        x_label: Optional[str] = None,
        show_raw: bool = False,
        show_raw_faint: bool = True,
        colors: Optional[List[str]] = None,
        figsize: Tuple[int, int] = (12, 8),
        fontsize: int = 25,
        step_yrs: int = 10,
        xlim_years: Optional[Tuple[int, int]] = None,
        ref_yrs: Optional[List[int]] = None,
        labels_zero_based: bool = False,
        outname: Optional[str] = None,
        legend_note: Optional[str] = None,
        legend_loc: str = "upper left",
        framealpha: float = 0.6,
        std_band: bool = False,
        std_center: str = "filtered",
        std_source: str = "filtered",
        std_mode: str = "global",
        std_window_months: int = 12,
        std_alpha: float = 0.15,
        panels: Optional[List[List[int]]] = None,
        panel_titles: Optional[List[str]] = None,
        exp_shas: Optional[dict] = None,
        sha_alpha: float = 0.4,
        legend_alpha: float = 0.5,
    ):
        """
        Plot the SST series. By default all experiments share one panel.
        panels: list of experiment-index lists, one per panel stacked vertically
                (e.g. [[0, 1], [0, 2]] -> top: exp0 vs exp1, bottom: exp0 vs exp2).
                Colors stay tied to the experiment index across panels.
        panel_titles: per-panel titles (default: `title` on the top panel only).
        exp_shas: hatched year windows per experiment, keyed by experiment index,
                  exp name or label, e.g. {1: {"segments": [(11, 20)], "hatch": "///",
                  "alpha": 0.4, "color": "#cc0000"}} (color defaults to the line color).
                  Drawn only in panels that contain that experiment.
        legend_alpha: opacity of the legend's white background box (drawn above hatching).
        """
        panels = panels or [list(range(len(self._meta)))]
        npan = len(panels)
        fig, axes = plt.subplots(npan, 1, figsize=figsize, sharex=True, squeeze=False)
        axes = axes[:, 0]
        palette = colors or ["black", "#cc0000", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]

        if self.compute_anomaly:
            title_default = f"{vname} Anomaly (Δ{vname})"
            y_label_default = f"Δ{vname} ({vunit})"
        else:
            title_default = vname
            y_label_default = f"{vname} ({vunit})"

        for p, (ax, idxs) in enumerate(zip(axes, panels)):
            ymins, ymaxs = [], []

            for i in idxs:
                meta = self._meta[i]
                col = (colors[i] if colors is not None else meta.get("color")) or palette[i % len(palette)]
                lw, a = float(meta["linewidth"]), float(meta["alpha"])
                yraw = np.asarray(self._raw[i], float)
                ysrc = np.asarray(self._filt[i], float) if (self.use_filter and (self._filt[i] is not None)) else yraw

                # hatched windows behind everything
                sha = self._sha_for(exp_shas, i)
                if sha and sha.get("segments"):
                    self._shade_segments(
                        ax, sha["segments"], color=sha.get("color", col),
                        hatch=sha.get("hatch", "///"), alpha=sha.get("alpha", sha_alpha),
                    )

                if std_band and yraw.size:
                    center = ysrc if ((std_center.lower() == "filtered") and self.use_filter and (self._filt[i] is not None)) else yraw
                    if (std_source.lower() == "residual") and self.use_filter and (self._filt[i] is not None):
                        base = yraw - ysrc
                    elif (std_source.lower() == "filtered") and self.use_filter and (self._filt[i] is not None):
                        base = ysrc
                    else:
                        base = yraw

                    if std_mode.lower() == "rolling":
                        win = std_window_months if not self.annual_mean else max(1, std_window_months // 12)
                        sigma = self._rolling_nanstd(base, win)
                    else:
                        s = float(np.nanstd(base, ddof=1)) if base.size > 1 else float(np.nanstd(base))
                        sigma = np.full_like(center, s, dtype=float)

                    ax.fill_between(np.arange(center.size), center - sigma, center + sigma,
                                    color=col, alpha=std_alpha, linewidth=0, zorder=0)
                    ymins.append(np.nanmin(center - sigma)); ymaxs.append(np.nanmax(center + sigma))

                if self.use_filter:
                    if show_raw:
                        ax.plot(yraw, color=col, linewidth=max(1.0, lw * 0.7), alpha=a * 0.6, linestyle="--")
                    elif show_raw_faint:
                        ax.plot(yraw, color=col, linewidth=max(1.0, lw * 0.3), alpha=a * 0.2)
                    ax.plot(ysrc, color=col, linewidth=lw, alpha=a)
                else:
                    ax.plot(ysrc, color=col, linewidth=lw, alpha=a)

                if ysrc.size:
                    ymins.append(np.nanmin(ysrc)); ymaxs.append(np.nanmax(ysrc))

            # FIX: honor labels_zero_based argument
            self._set_time_axis(
                ax,
                is_annual=self.annual_mean,
                step_yrs=step_yrs,
                xlim_years=xlim_years,
                ref_yrs=ref_yrs,
                fontsize=fontsize,
                labels_zero_based=labels_zero_based
            )

            if panel_titles is not None:
                ptitle = panel_titles[p]
            else:
                ptitle = (title or title_default) if p == 0 else None
            if ptitle:
                ax.set_title(ptitle, fontsize=fontsize * 0.9, loc='left', pad=15)
            if p == npan - 1:
                ax.set_xlabel(x_label or "Model Time (year)", fontsize=fontsize * 0.9)
            ax.set_ylabel(y_label or y_label_default, fontsize=fontsize * 0.9)

            if ylim:
                ax.set_ylim(*ylim)
            elif ymins and ymaxs:
                lo, hi = float(np.nanmin(ymins)), float(np.nanmax(ymaxs))
                pad = 0.05 * max(1e-12, hi - lo)
                ax.set_ylim(lo - pad, hi + pad)

            # zero line only when it lies inside the axes (absolute SST is far from 0)
            y0, y1 = ax.get_ylim()
            if min(y0, y1) <= 0.0 <= max(y0, y1):
                ax.axhline(0.0, ls="--", lw=1.2, color="k", alpha=0.8, zorder=1000, clip_on=False)
            ax.grid(True, alpha=0.25)
            ax.tick_params(labelsize=fontsize * 0.9, length=6, width=1)

            default_palette = ["black", "#cc0000", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]
            leg = self._build_legend(
                self._stats, ax, vunit=vunit, fontsize=fontsize,
                legend_loc=legend_loc, colors=(colors or default_palette),
                note=legend_note, ci_format="pm", indices=idxs
            )
            # opaque, borderless box so hatching/lines never show through the legend text
            leg.set_frame_on(True)
            leg.get_frame().set_facecolor("white")
            leg.get_frame().set_edgecolor("none")
            leg.get_frame().set_alpha(legend_alpha)
            leg.set_zorder(2000)
            for sp in ax.spines.values():   # keep axis frame above the legend box
                sp.set_zorder(2500)

        if outname:
            plt.savefig(outname, dpi=600, bbox_inches="tight")
        return fig, (axes[0] if npan == 1 else axes)


class OHCHovmollerPlotter(OCNTimeSeriesPlotter):
    """Global OHC depth-time Hovmoller comparison."""

    def __init__(
        self,
        root_dir,
        sub_dir,
        file_key="mpasTimeSeriesOcean",
        mpas_depth="depth.nc",
        layer_info=None,
        var_name=None,
        compute_anomaly=False,
        running_mean=False,
        running_mean_months=0,
        running_mean_center="center",
        running_mean_min_frac=0.5,
        region_dim=None,
        region_index=None,
        region_name=None,
        region_reduce=None,
        region_names_var="regionNames",
        scale=1.0,
        baseline_years=1,
        baseline_months=12,
        filt_cutoff_months=24,
        filt_order=4,
        trend_years=30,
        trend_unit='per_year',  # or 'per_decade'
        year_min=1,
        year_max=350,
        use_filter=False,
        filt_method='pad',
        filt_padlen=None,
        filt_padtype='odd',
        edge_fix_points=0,
        use_mfdataset=True,
        chunks=None,
        engine=None,
        time_origin=None,
        calendar="noleap"
    ):
        # I/O
        self.time_origin = time_origin
        self.calendar = calendar
        self.root_dir = root_dir
        self.sub_dir = sub_dir
        self.file_key = file_key
        self.mpas_depth = mpas_depth
        self.layer_info = layer_info
        self.var_name = var_name
        if self.var_name is None:
            raise ValueError("var_name must be provided.")

        self.baseline_months = int(baseline_months)
        self.baseline_years = (int(baseline_years)
                               if baseline_years is not None
                               else max(1, int(np.ceil(self.baseline_months / 12.0))))

        self.region_dim       = region_dim
        self.region_index     = None if region_index is None else int(region_index)
        self.region_name      = None if region_name  is None else str(region_name)
        self.region_reduce    = None if region_reduce is None else str(region_reduce).lower()
        self.region_names_var = None if region_names_var in (None, "None") else str(region_names_var)

        # Running mean flags
        self.running_mean = bool(running_mean)
        self.compute_anomaly = bool(compute_anomaly)
        self.running_mean_months = 0 if running_mean_months in (None, 0) else int(running_mean_months)
        self.running_mean_center = str(running_mean_center).lower()
        self.running_mean_min_frac = float(running_mean_min_frac)

        self.scale = float(scale)
        self.filt_cutoff_months = float(filt_cutoff_months)
        self.filt_order = int(filt_order)
        self.trend_years = int(trend_years)
        self.trend_unit = str(trend_unit)
        self.year_min = int(year_min)
        self.year_max = int(year_max)
        self.use_filter = bool(use_filter)
        self.filt_method   = str(filt_method)
        self.filt_padlen   = None if filt_padlen is None else int(filt_padlen)
        self.filt_padtype  = str(filt_padtype)
        self.edge_fix_points = int(edge_fix_points)

        # xarray options
        self._xr_use_mfdataset = bool(use_mfdataset)
        self._xr_chunks = chunks
        self._xr_engine = engine

        # Filter design if needed
        if self.use_filter:
            fc = 1.0 / self.filt_cutoff_months
            wn = max(min(fc / 0.5, 0.999999), 1e-6)
            self._butter_b, self._butter_a = signal.butter(self.filt_order, wn, btype="low")
        else:
            self._butter_b = self._butter_a = None

        # Storage
        self._meta, self._raw, self._filt, self._stats = [], [], [], []
        
        self.panel_tags = "lower"      # "lower" | "upper" | None
        self.panel_tag_offset = 0      # start index offset

    def _find_moc_timeseries_files(self, root_dir, exp, sub_dir, file_key, year_token):
        if isinstance(year_token, int):
            yr0, yr1 = year_token, year_token
        elif isinstance(year_token, str):
            yr0, yr1 = map(int, year_token.split("-"))
        else:
            yr0, yr1 = year_token
        base = os.path.join(root_dir, exp, sub_dir)
        selected = []
        for yy in np.arange(yr0, yr1+1):
            pattern = os.path.join(base, f"{file_key}_{yy:04d}-{yy:04d}.nc")
            matches = glob.glob(pattern)
            for fname in matches:
                m = re.search(r"_(\d{4})-\d{4}\.nc$", fname)
                if m:
                    selected.append((int(m.group(1)), fname))
        selected.sort(key=lambda x: x[0])
        return [f for _, f in selected]

    def _find_mpas_timeseries_files(self, root_dir, exp, sub_dir, file_key):
        base = os.path.join(root_dir, exp, sub_dir)
        single = glob.glob(os.path.join(base, f"{file_key}.nc"))
        if single:
            return single
        segs = glob.glob(os.path.join(base, f"{file_key}_*.nc"))
        if not segs:
            return []
        def start_stamp(path):
            m = re.search(rf"{re.escape(file_key)}_(\d{{6}})-(\d{{6}})\.nc$", os.path.basename(path))
            return int(m.group(1)) if m else 0
        segs.sort(key=start_stamp)
        return segs

    def _open_dataset(self, files):
        xr_opts = dict(engine=self._xr_engine, chunks=self._xr_chunks, decode_times=True)
        if len(files) == 1 or not self._xr_use_mfdataset:
            return xr.open_dataset(files[0], **xr_opts)
        try:
            return xr.open_mfdataset(files, combine="by_coords", parallel=True, **xr_opts)
        except Exception:
            dsets = [xr.open_dataset(fp, **xr_opts) for fp in files]
            probe = dsets[0][self.var_name] if self.var_name in dsets[0] else list(dsets[0].data_vars.values())[0]
            tdim = self._guess_time_dim_da(probe)
            return xr.concat(dsets, dim=tdim)

    def _load_thickness(self, ds, lev_dim, fdepth):
        fds = xr.open_dataset(fdepth)
        cand = None
        for k in ("refLayerThickness","layerThickness","dz","thickness"):
            if k in fds: cand = fds[k]; break
        if cand is None and "refBottomDepth" in fds:
            b = fds["refBottomDepth"]; cand = b
        if cand is None:
            raise KeyError("Could not find layer thickness in the fdepth file.")
        fz = cand.squeeze()
        if fz.ndim != 1 or fz.size != ds.sizes[lev_dim]:
            if fz.ndim == 1 and fz.size == ds.sizes[lev_dim]:
                fz = xr.DataArray(fz.values, dims=[lev_dim],
                                  coords={lev_dim: ds[lev_dim] if lev_dim in ds.coords else np.arange(ds.sizes[lev_dim])})
            else:
                raise ValueError("Thickness from file must be 1D and match model levels.")
        return fz

    def _open_timeseries_ds(self, exp, year_token, sub_dir_override=None):
        """Open MPAS time-series dataset(s), slice by year token, and return (ds, tdim, time_idx)."""
        subdir = sub_dir_override if sub_dir_override is not None else self.sub_dir
        files = self._find_mpas_timeseries_files(self.root_dir, exp, subdir, self.file_key)
        if not files:
            search_root = os.path.join(self.root_dir, exp)
            candidates = glob.glob(os.path.join(
                search_root, "post", "analysis", "mpas_analysis", "ts_*_climo_*", "**", f"{self.file_key}*.nc"
            ), recursive=True)
            if candidates:
                hit = candidates[0]
                subdir = os.path.relpath(os.path.dirname(hit), start=os.path.join(self.root_dir, exp))
                files = self._find_mpas_timeseries_files(self.root_dir, exp, subdir, self.file_key)
        if not files:
            raise FileNotFoundError(
                f"No '{self.file_key}.nc' (or segmented) found for '{exp}' under subdir='{subdir}'."
            )
        ds = self._open_dataset(files)

        # pick a representative var to guess time dim
        probe = ds[self.var_name] if self.var_name in ds else next(iter(ds.data_vars.values()))
        tdim = self._guess_time_dim_da(probe)
        n = ds[tdim].sizes[tdim] if tdim in ds.dims else None

        # build time index (CF, xtime, or constructed)
        time_idx = self._decode_time_from_Time_var(ds)
        if time_idx is None:
            ti_start = None; ti_end = None
            for nm in ("xtime_startMonthly", "startTime"):
                ti_start = self._decode_time_from_xtime_chars(ds, nm)
                if ti_start is not None and (n is None or len(ti_start) == n):
                    break
            for nm in ("xtime_endMonthly", "endTime"):
                ti_end = self._decode_time_from_xtime_chars(ds, nm)
                if ti_end is not None and (n is None or len(ti_end) == n):
                    break
            time_idx = ti_start if (ti_start is not None) else ti_end
        if time_idx is None and n is not None:
            time_idx = self._construct_monthly_time(n)

        # slice by years if possible, else by months
        if time_idx is not None:
            se = self._slice_indices_from_year_token(time_idx, year_token)
            if se is not None:
                a, b = se
                ds = ds.isel({tdim: slice(a, b)})
                time_idx = time_idx[a:b]
        else:
            m_start, m_end, _ = self._find_time_slice(year_token)
            ds = ds.isel({tdim: slice(m_start, m_end)})

        return ds, tdim, time_idx

    def _sel_region(self, da, ds=None):
        """Select a single region if a region dimension exists; otherwise return da unchanged."""
        reg_dim = next((d for d in ("nOceanRegions","region","nRegions") if d in da.dims), None)
        if reg_dim is None:
            return da
        # prefer a name if provided & we have ds to resolve it
        if (self.region_name is not None) and (ds is not None):
            idx = self._region_index_from_name(ds, reg_dim, self.region_name)
        else:
            idx = int(self.region_index) if self.region_index is not None else 0
        return da.isel({reg_dim: idx})

    def _force_2d_time_lev(self, da, time_dim=None, lev_dim=None):
        """
        Return DataArray as 2D (time, lev). Squeeze singleton dims; mean over others.
        """
        time_dim, lev_dim = self._infer_time_lev_dims(da, time_dim=time_dim, lev_dim=lev_dim)

        # Squeeze or average leftover dims
        for d in list(da.dims):
            if d not in (time_dim, lev_dim):
                da = da.squeeze(d) if da.sizes[d] == 1 else da.mean(d)

        # Reorder to (time, lev)
        if da.dims != (time_dim, lev_dim):
            da = da.transpose(time_dim, lev_dim)

        return da

    def _infer_time_lev_dims(self, obj, time_dim=None, lev_dim=None):
        """
        Infer names of time and vertical (level/depth) dimensions from an xarray
        Dataset/DataArray. If explicit names are provided, use them.
        """
        import numpy as np
        dims = tuple(obj.dims)

        # --- time dim ---
        if time_dim is None:
            # Look for a datetime-like coordinate first
            for d in dims:
                coord = obj.coords.get(d, None)
                if coord is not None and np.issubdtype(coord.dtype, np.datetime64):
                    time_dim = d
                    break
            # Common fallbacks
            if time_dim is None:
                for cand in ("time", "Time", "xtime", "t"):
                    if cand in dims:
                        time_dim = cand
                        break
            # Last resort: any dim whose name contains 'time'
            if time_dim is None:
                for d in dims:
                    if "time" in d.lower():
                        time_dim = d
                        break
        if time_dim is None:
            raise KeyError(f"Could not infer time dimension from dims={dims}")

        # --- vertical/level dim ---
        if lev_dim is None:
            for cand in ("nVertLevels", "lev", "level", "depth", "z_t", "z", "layer", "nz"):
                if cand in dims:
                    lev_dim = cand
                    break
            if lev_dim is None:
                # pick any other non-time dim with size > 1
                for d in dims:
                    if d != time_dim and obj.sizes.get(d, 0) > 1:
                        lev_dim = d
                        break
        if lev_dim is None:
            raise KeyError(f"Could not infer vertical dimension from dims={dims}")

        return time_dim, lev_dim

    def _parse_start_year_from_units(units: Optional[str]) -> Optional[int]:
        """
        Try to extract a start year from a CF-style units string like:
          'days since 0025-01-15 00:00:00'
        Returns an int year (e.g., 25) or None if not found.
        """
        if not units:
            return None
        m = re.search(r"since\s+(-?\d{1,4})", units)
        if not m:
            # broader pattern to catch 'YYYY-MM-DD'
            m = re.search(r"(-?\d{1,4})-\d{1,2}-\d{1,2}", units)
        if not m:
            return None
        try:
            return int(m.group(1))
        except Exception:
            return None

    def _to_annual_mean_da(self, da: xr.DataArray, ymin: int, ymax: int) -> xr.DataArray:
        """
        Convert monthly (time, depth) data to calendar-year means and
        return only years in [ymin, ymax]. Requires a datetime/CF-time
        coordinate on the time dimension; otherwise raises.
        """
        time_dim, lev_dim = self._infer_time_lev_dims(da)
        coord = da.coords.get(time_dim, None)

        # Must have datetime-like or cftime-like time
        if coord is None or not hasattr(coord, "dt"):
            raise ValueError(
                f"_to_annual_mean_da requires a datetime/CF-time coordinate on '{time_dim}'. "
                "Open with decode_times=True or attach a CFTimeIndex first."
            )

        # Optional sanity: warn if any year has <12 monthly samples (could be incomplete)
        # (We don't try to 'fix' this — just warn.)
        try:
            counts_per_year = xr.ones_like(da.isel({time_dim: 0})).groupby(f"{time_dim}.year").count(dim=...)
        except Exception:
            # safer: compute counts from coord
            years = coord.dt.year.values
            _, counts = np.unique(years, return_counts=True)
            if np.any(counts < 12):
                print(f"[warn] Some years have <12 samples. Annual mean will average available months only.", flush=True)

        # Annual mean via groupby
        out = da.groupby(f"{time_dim}.year").mean(dim=time_dim, skipna=True)

        # Rename 'year' -> time_dim and attach integer-year coord
        if "year" in out.dims:
            out = out.rename({"year": time_dim})
        years = xr.DataArray(np.unique(coord.dt.year.values), dims=(time_dim,))
        out = out.assign_coords({time_dim: years})

        # Subset to requested range
        out = out.sel({time_dim: slice(int(ymin), int(ymax))})

        # Diagnostic string
        yrs = out[time_dim].values
        if yrs.size > 0:
            yr_min, yr_max = int(yrs[0]), int(yrs[-1])
            diag = f"[diag] Annual mean years subset: {yr_min}–{yr_max} (requested {ymin}–{ymax})"
        else:
            diag = f"[diag] Annual mean returned EMPTY for requested range {ymin}–{ymax}"
        print(diag, flush=True)

        # Metadata
        out.attrs = {**da.attrs, **{
            "cell_methods": f"{time_dim}: mean (annual from monthly; groupby calendar year)",
            "notes": f"Annual mean; subset years {ymin}–{ymax} (inclusive).",
            "diagnostic_years": diag,
        }}
        return out

    def _ohc_hovmoller_matrix(
        self,
        ds,
        *,
        rho=1026.0,
        cp=3996.0,
        to_units="GJ/m^2",
        anomaly=True,
        baseline_mode="calendar_year",
        first_N=12,
        mode="per_area",                # "per_area" (J/m^2) or "total" (J)
        normalize_by_surface_area=False # if mode="total": divide by surface area to get J/m^2
    ):
        time_dim, lev_dim = self._infer_time_lev_dims(ds)

        # required MPAS regional-layer means
        t_name = "timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature"
        h_name = "timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerThickness"
        a_name = "timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerArea"
        m_name = "timeMonthly_avg_avgValueWithinOceanLayerRegion_sumLayerMaskValue"
        for k in (t_name, h_name, a_name, m_name):
            if k not in ds: raise KeyError(f"Missing required variable: {k}")

        # region select and force (time, lev)
        T = self._force_2d_time_lev(self._sel_region(ds[t_name], ds))
        H = self._force_2d_time_lev(self._sel_region(ds[h_name], ds))
        Aavg  = self._force_2d_time_lev(self._sel_region(ds[a_name], ds))
        Ncell = self._force_2d_time_lev(self._sel_region(ds[m_name], ds))

        # --- depth coordinate ---
        with xr.open_dataset(self.mpas_depth) as fds:
            if "refBottomDepth" in fds:
                # use reference bottom depth
                z_bot_ref = fds["refBottomDepth"].squeeze()  # (lev_dim,)
                # convert to km if needed
                z_bot_ref = z_bot_ref * 1e-3  
                z_bot_ref = z_bot_ref.rename("depth")
                z_bot_ref.attrs.update(units="km", positive="down")
                depth = z_bot_ref
            else:
                # derive from mean thickness if refBottomDepth not available
                Hmean = H.mean(dim=time_dim, skipna=True)   # mean thickness per layer (m)
                # bottom depth: cumulative sum of layer thickness
                z_bot = Hmean.cumsum(dim=lev_dim)
                # top depth: bottom depth shifted by one layer (top of first layer is 0)
                z_top = z_bot.shift({lev_dim: 1}, fill_value=0.0)
                # midpoint: average of top and bottom
                z_mid = (0.5 * (z_top + z_bot)).rename("depth")
                z_mid = z_mid * 1e-3  # m to km 
                z_mid.attrs.update(units="km", positive="down")
                depth = z_mid
        
        # --- per-layer OHC per unit area ---
        ohc_pa = rho * cp * T * H                 # J/m^2 per layer

        # --- per-layer TOTAL OHC (J) using effective layer area ---
        Aeff_layer = Aavg * Ncell                 # total wet area in that layer
        ohc_total  = ohc_pa * Aeff_layer          # J per layer
        
        # --- per-layer OHC per unit area ---
        ohc = ohc_pa if mode == "per_area" else ohc_total
        intended_units = "J m^-2" if mode == "per_area" else "J"
        ohc.name = "OHC_layer"
        ohc.attrs["units"] = intended_units

        # anomalies can drop attrs -> capture then restore
        if anomaly:
            if baseline_mode == "calendar_year" and hasattr(ohc[time_dim], "dt"):
                yr0  = int(ohc[time_dim].dt.year.min())
                base = ohc.where(ohc[time_dim].dt.year == yr0, drop=True).mean(time_dim, skipna=True)
            else:
                base = ohc.isel({time_dim: slice(0, min(first_N, ohc.sizes.get(time_dim, first_N)))}) \
                       .mean(time_dim, skipna=True)
            ohc = ohc - base
            # restore units after arithmetic
            ohc.attrs["units"] = intended_units

        if (mode == "total") and normalize_by_surface_area:
            A_region = (Aeff_layer.isel({lev_dim: 0}))
            A_region = xr.where(A_region > 0, A_region, np.nan)
            ohc = ohc / A_region
            intended_units = "J m^-2"   # update intended units
            ohc.attrs["units"] = intended_units

        # unit scaling — use .get(...) to avoid KeyError
        u = ohc.attrs.get("units", intended_units)
        if str(to_units).lower() in ("gj/m^2", "gj m^-2") and u == "J m^-2":
            ohc = ohc / 1e9
            ohc.attrs["units"] = "GJ m^-2"
        elif str(to_units).lower() in ("1e22 j", "zj") and u == "J":
            ohc = ohc / 1e22
            ohc.attrs["units"] = "10$^{22}$ J"
            
        # make sure depth uses the lev_dim dimension
        if depth.dims[0] != lev_dim:
            depth = depth.rename({depth.dims[0]: lev_dim})

        # turn the lev_dim → depth mapping into the new dimension index
        ohc = ohc.assign_coords({lev_dim: depth})
        ohc = ohc.swap_dims({lev_dim: "depth"}).sortby("depth")

        # Now explicitly set the coordinate values = depth.values
        ohc = ohc.assign_coords(depth=depth.values)
        
        if lev_dim in ohc.coords:
            ohc = ohc.drop_vars(lev_dim)
    
        # re-assert units just in case
        ohc.attrs["units"] = ohc.attrs.get("units", intended_units)

        return ohc, ohc["depth"]

    def alpha_tag(self, idx):
        # 0 -> a/A, 25 -> z/Z, 26 -> aa/AA, etc.
        s, n = "", idx + (self.panel_tag_offset or 0)
        while True:
            s = chr(97 + (n % 26)) + s
            n = n // 26 - 1
            if n < 0: break
        if self.panel_tags == "upper": s = s.upper()
        return f"({s})"

    def _draw_row(
        self,
        ax_list,
        start_year,
        depth_ref,
        data_list,
        labels,
        vmin,
        vmax,
        *,
        row_title=None,
        annual=True,
        fig_idx=0,
        xtick_interval=None,
        vlines_per_panel=None,      # list of lists, one per panel (years or indices)
        vline_mode="year",          # "year" (default) or "index"
        vline_kwargs=None,
        vline_spin=221,
        cmap="RdBu_r",
        fontsize=18,
    ):
        """
        Draw one Hovmöller row across panels.

        annual=True  -> the time axis is in years (one sample per year)
        annual=False -> the time axis is in months (12 samples per year)
        xtick_interval: spacing in YEARS for labeled ticks (both modes)
        vlines_per_panel: [[years...], [years...], ...] or indices if vline_mode="index"
        """
        if vline_kwargs is None:
            vline_kwargs = dict(color="k", linestyle="--", linewidth=1.0, alpha=0.7)
        if vlines_per_panel is None:
            vlines_per_panel = [[] for _ in range(len(ax_list))]

        # --- Build y-edges once from the shared depth coord
        z = depth_ref.values.astype(float)
        depth_units = (getattr(depth_ref, "attrs", {}) or {}).get("units", "")
        if str(depth_units).lower().startswith("m"):
            z_plot = z / 1000.0
            y_label_unit = "km"
        else:
            z_plot = z
            y_label_unit = "km" if depth_units == "" else depth_units

        if z_plot.size < 2:
            dz = 1.0
            y_edges = np.array([z_plot[0] - 0.5 * dz, z_plot[0] + 0.5 * dz])
        else:
            y_edges = np.empty(z_plot.size + 1, dtype=float)
            y_edges[1:-1] = 0.5 * (z_plot[1:] + z_plot[:-1])
            y_edges[0]    = z_plot[0]  - (y_edges[1] - z_plot[0])
            y_edges[-1]   = z_plot[-1] + (z_plot[-1] - y_edges[-2])

        mappable = None

        for pi, (ax, label, da) in enumerate(zip(ax_list, labels, data_list)):
            # time edges in index space
            nt = da.sizes[da.dims[0]]
            x_edges = np.arange(nt + 1)

            # image
            Z = np.asarray(da.values, float).T
            im = ax.pcolormesh(x_edges, y_edges, Z, cmap=cmap, shading="auto",
                               vmin=vmin, vmax=vmax)
            mappable = im

            # cosmetics
            ax.invert_yaxis()
            tag = self.alpha_tag(pi+fig_idx) if self.panel_tags else "" 
            panel_label = f"{tag} {label}".strip()
            ax.set_title(panel_label, fontsize=fontsize, loc="left") 
            ax.set_title(row_title, fontsize=fontsize, loc="right") 
            ax.tick_params(labelsize=fontsize*0.95, length=5, width=1)

            # --- X ticks (index-space -> year labels)
            nx = nt
            if annual:
                total_years = nx
                tick_interval = (int(xtick_interval) if xtick_interval is not None
                                 else max(1, int(np.ceil(total_years / max(1, min(8, total_years))))))
                xticks = np.arange(0, total_years + 1, tick_interval)
                xlabels = start_year + xticks - 1
            else:
                total_years = nx // 12
                tick_interval = (int(xtick_interval) if xtick_interval is not None
                                 else max(1, int(np.ceil(total_years / max(1, min(8, total_years))))))
                xticks = np.arange(0, nx + 1, 12 * tick_interval)
                xlabels = start_year + (xticks // 12) -1 

            ax.set_xticks(xticks)
            ax.set_xticklabels([f"{int(y)}" for y in xlabels])
            ax.xaxis.set_major_locator(mticker.FixedLocator(xticks))
            ax.set_xlabel("Model Year", fontsize=fontsize * 0.95)

            # --- Vertical reference lines for THIS panel
            lines = vlines_per_panel[pi] if pi < len(vlines_per_panel) else []
            for val in (lines or []):
                if vline_mode == "year":
                    xpos = (val - start_year) if annual else (val - start_year) * 12
                else:  # "index"
                    xpos = int(val)
                if 0 <= xpos <= nx:
                    ax.axvline(x=xpos, **vline_kwargs)
                    
            if vline_spin is not None:
                if vline_mode == "year":
                    xpos = (vline_spin - start_year) if annual else (vline_spin - start_year) * 12
                else:  # "index"
                    xpos = int(vline_spin)
                if 0 <= xpos <= nx:
                    ax.axvline(x=xpos, color="tab:red", linestyle="--", linewidth=2.0, alpha=1.0)

        # y label only on the left
        ax_list[0].set_ylabel(f"Depth ({y_label_unit})", fontsize=fontsize * 0.95)
        for ax in ax_list[1:]:
            ax.set_yticklabels([])
            ax.set_ylabel("")

        return mappable

    def _auto_symmetric_minmax(self, arr_list):
        dmin = min(np.nanmin(o.values) for o in arr_list)
        dmax = max(np.nanmax(o.values) for o in arr_list)
        a = max(abs(dmin), abs(dmax))
        return -a, a

    def _cumulative_along_depth(self, da, *, from_top=True):
        """
        Cumulate along 'depth'.
        from_top=True  -> surface→bottom
        from_top=False -> bottom→surface
        """
        # ensure depth is increasing downward
        if da["depth"].values[0] > da["depth"].values[-1]:
            da = da.sortby("depth")

        if from_top:
            out = da.cumsum(dim="depth")
        else:
            # reverse -> cumsum -> reverse back
            out = da.isel(depth=slice(None, None, -1)).cumsum(dim="depth").isel(depth=slice(None, None, -1))

        return out

    def save_ohc_depth_netcdf(
        self,
        ohc_da: xr.DataArray,
        depth_da: Union[xr.DataArray, np.ndarray],
        out_path: str,
        *,
        var_name: str = "ohc",
        z_name: str = "z",
        time_name: str = "time",         # <--- new: desired canonical time dim name
        title: Optional[str] = None,
        overwrite: bool = True
    ) -> str:
        """
        Save OHC(time,z) and depth(z) to a NetCDF file with CF-friendly metadata.
        Auto-detects and renames the time dim (e.g., 'Time' -> 'time').
        """

        # --------- identify/normalize time dim (case-insensitive) ----------
        # First try exact/CI match
        time_dim = None
        for d in ohc_da.dims:
            if d.lower() == "time":
                time_dim = d
                break
        # If still not found, try a datetime-like heuristic
        if time_dim is None:
            for d in ohc_da.dims:
                c = ohc_da.coords.get(d, None)
                if c is not None and (
                    np.issubdtype(getattr(c, "dtype", object), np.datetime64) or
                    hasattr(c, "dt") or
                    "since" in str(getattr(c, "units", "")).lower()
                ):
                    time_dim = d
                    break
        if time_dim is None:
            raise ValueError(
                f"Could not identify a time dimension in ohc_da dims {ohc_da.dims}. "
                "Ensure there is a coordinate named 'time'/'Time' or with datetime-like values."
            )

        # Rename time dim to canonical 'time_name' if needed
        if time_dim != time_name:
            ohc_da = ohc_da.rename({time_dim: time_name})
            time_dim = time_name  # update local var

        # --------- identify/normalize vertical dim ----------
        z_dim_candidates = [d for d in ohc_da.dims if d != time_dim]
        if len(z_dim_candidates) != 1:
            raise ValueError(
                f"ohc_da must be 2-D (time,z). Found dims: {tuple(ohc_da.dims)}"
            )
        z_dim = z_dim_candidates[0]

        # Wrap depth into a DataArray if it isn't already
        if not isinstance(depth_da, xr.DataArray):
            depth_da = xr.DataArray(depth_da, dims=(z_dim,))

        # Sanity check lengths
        if depth_da.sizes[z_dim] != ohc_da.sizes[z_dim]:
            raise ValueError(
                f"Depth length ({depth_da.sizes[z_dim]}) != ohc vertical size ({ohc_da.sizes[z_dim]})"
            )

        # Standardize vertical dim name to z_name
        if z_dim != z_name:
            ohc_da   = ohc_da.rename({z_dim: z_name})
            depth_da = depth_da.rename({z_dim: z_name})
            z_dim = z_name

        # --------- metadata defaults ----------
        depth_da = depth_da.copy()
        depth_da.attrs.setdefault("standard_name", "depth")
        depth_da.attrs.setdefault("long_name", "depth below sea surface")
        depth_da.attrs.setdefault("units", "m")
        depth_da.attrs.setdefault("positive", "down")

        ohc_to_save = ohc_da.copy()
        ohc_to_save.attrs.setdefault("long_name", "Ocean heat content (layered)")
        ohc_to_save.attrs.setdefault("units", ohc_to_save.attrs.get("units", "J"))

        # Build dataset
        time_coord = ohc_to_save[time_dim]
        ds_out = xr.Dataset({var_name: ohc_to_save},
                            coords={time_dim: time_coord, z_dim: depth_da})

        # --------- globals ----------
        ds_out.attrs.update({
            "Conventions": "CF-1.8",
            "title": title or "OHC Hovmöller output",
            "history": f"{datetime.now(timezone.utc):%Y-%m-%d %H:%M:%S}Z: created by save_ohc_depth_netcdf",
            "source": ohc_da.attrs.get("source", "derived from model diagnostics"),
            "references": "Griffies (2016) for ocean heat budget; model-specific docs."
        })

        # --------- encoding (preserve calendar if present) ----------
        time_enc = {}
        if hasattr(time_coord, "encoding"):
            for k in ("units", "calendar"):
                if k in time_coord.encoding:
                    time_enc[k] = time_coord.encoding[k]

        comp = dict(zlib=True, complevel=4, dtype="float32", _FillValue=np.float32(np.nan))
        enc = {
            var_name: comp,
            z_dim: dict(zlib=True, complevel=4, dtype="float32"),
            time_dim: time_enc or {}
        }
        
        # --------- write ----------
        if (not overwrite) and os.path.exists(out_path):
            raise FileExistsError(f"{out_path} exists and overwrite=False")
        ds_out.to_netcdf(out_path, encoding=enc)
        return out_path

    def plot_ohc_hovmoller(
        self,
        exps,
        labels=None,
        year_token="0001-0350",
        subdirs=None,
        outpath=None,
        figsize_per_panel=(6, 4),
        fontsize=18,
        cmap="RdBu_r",
        vmin_row1=None, vmax_row1=None,     # per-layer OHC color scale
        vmin_row2=None, vmax_row2=None,     # cumulative OHC color scale
        to_units="GJ/m^2",
        mode="total",
        normalize_by_surface_area=False,
        anomaly=True,
        baseline_mode="calendar_year",
        first_N=12,
        title_top="Per-layer ΔOHC",
        title_bottom="Cumulative ΔOHC",
        draw_cumulative=True,
        super_title=None,
        annual_first=False,
        xtick_interval=None,
        vlines_per_exp=None,         # <-- NEW: list of lists (per experiment)
        vline_mode="year",           # <-- "year" or "index"
        vline_kwargs=None,           # <-- style dict
        vline_spin=221,
        from_top=False, 
        outname=None,
    ):
        """
        Two-row Hovmöller: row 1 = per-layer ΔOHC, row 2 = cumulative ΔOHC (surface→bottom).
        Each row has its own (shared-across-panels) colorbar.
        """
        if labels is None: labels = exps
        if subdirs is None: subdirs = [self.sub_dir] * len(exps)
        assert len(labels) == len(exps) == len(subdirs)

        # Resolve start year for x-labels
        try:
            y0, y1 = [int(t.strip()) for t in str(year_token).split("-")]
            start_year = min(y0, y1)
        except Exception:
            start_year = self.year_min

        # ---- Load panels once; build both per-layer and cumulative lists on a common depth grid
        per_list, cum_list = [], []
        depth_ref = None

        for exp, subdir in zip(exps, subdirs):
            ds, _tdim, _time_idx = self._open_timeseries_ds(exp, year_token, sub_dir_override=subdir)

            ohc, depth = self._ohc_hovmoller_matrix(
                ds,
                mode=mode,
                normalize_by_surface_area=normalize_by_surface_area,
                to_units=to_units,
                anomaly=anomaly,
                baseline_mode=baseline_mode,
                first_N=first_N,
            )

            if annual_first:
                # annualize per-layer first
                ohc = self._to_annual_mean_da(ohc,y0, y1)
            
            # Ensure output dir exists
            os.makedirs(outpath, exist_ok=True)

            # (Optional) sanitize pieces for filenames
            safe_exp  = str(exp).replace("/", "-")
            safe_mode = str(mode).replace("/", "-")
            safe_year = str(year_token)

            out_nc = os.path.join(
                outpath,  # use outpath (the directory you just created)
                f"ohc_hov_{safe_exp}_{safe_mode}_annual_{safe_year}.nc"
            )

            # Save (ohc: DataArray with dims ('time', z); depth: 1-D aligned with z)
            self.save_ohc_depth_netcdf(
                ohc_da=ohc,
                depth_da=depth,
                out_path=out_nc,
                var_name="ohc",
                z_name="z",  # set to "depth" if you prefer
                title=f"OHC Hovmöller ({mode}, anomaly={bool(anomaly)})"
            )
            print("Wrote:", out_nc)

            # cumulative (choose direction)
            if draw_cumulative:
                # per-layer field 'ohc' already aligned on depth
                ohc_c = self._cumulative_along_depth(ohc, from_top=from_top)  # False = bottom→surface
                ohc_c = ohc_c.rename("OHC_cumulative")
                ohc_c.attrs["units"] = ohc.attrs.get("units", to_units or "")

            if depth_ref is None:
                depth_ref = depth
                # If annualized, depth_ref is unchanged (depth-only coord)
                per_list.append(ohc)
                if draw_cumulative:
                    cum_list.append(ohc_c)
            else:
                if not np.allclose(depth.values, depth_ref.values, rtol=0, atol=1e-6):
                    per_list.append(ohc.interp(depth=depth_ref))
                    if draw_cumulative:
                        cum_list.append(ohc_c.interp(depth=depth_ref))
                else:
                    per_list.append(ohc)
                    if draw_cumulative:
                        cum_list.append(ohc_c)

        # ---- Compute color scales (independent, symmetric if not provided)
        if vmin_row1 is None or vmax_row1 is None:
            a, b = self._auto_symmetric_minmax(per_list)
            if vmin_row1 is None and vmax_row1 is None:
                vmin_row1, vmax_row1 = a, b
            elif vmin_row1 is None:
                vmin_row1 = -max(abs(a), abs(vmax_row1))
            else:  # vmax_row1 is None
                vmax_row1 =  max(abs(b), abs(vmin_row1))
                
        if draw_cumulative:
            if vmin_row2 is None or vmax_row2 is None:
                a, b = self._auto_symmetric_minmax(cum_list)
                if vmin_row2 is None and vmax_row2 is None:
                    vmin_row2, vmax_row2 = a, b
                elif vmin_row2 is None:
                    vmin_row2 = -max(abs(a), abs(vmax_row2))
                else:  # vmax_row2 is None
                    vmax_row2 =  max(abs(b), abs(vmin_row2))

        # ---- Figure layout: 2 rows, panels + 1 colorbar column
        W, H = figsize_per_panel
        fig = plt.figure(figsize=(W * len(exps) + 0.8*W, 2 * H), constrained_layout=True)
        
        ncol = len(exps)
        if ncol == 0:
            raise ValueError("No experiments to plot.")

        gs = GridSpec(
            nrows=2, ncols=ncol + 1,  # <-- +1 for colorbar column
            width_ratios=[1] * ncol + [0.05],  # last narrow column for cbar
            wspace=0.10, hspace=0.00, figure=fig
        )
        
        axes_row1 = [fig.add_subplot(gs[0, i]) for i in range(len(exps))]
        cax1 = fig.add_subplot(gs[0, -1])
        if draw_cumulative:
            axes_row2 = [fig.add_subplot(gs[1, i]) for i in range(len(exps))]
            cax2 = fig.add_subplot(gs[1, -1])

        # --- Depth edges
        z = depth_ref.values.astype(float)
        depth_units = (getattr(depth_ref, "attrs", {}) or {}).get("units", "")
        if str(depth_units).lower().startswith("m"):
            z = z / 1000.0
            y_label_unit = "km"
        else:
            y_label_unit = "km" if depth_units == "" else depth_units

        if z.size < 2:
            dz = 1.0
            y_edges = np.array([z[0] - 0.5*dz, z[0] + 0.5*dz])
        else:
            z_edges = np.empty(z.size + 1, dtype=float)
            z_edges[1:-1] = 0.5 * (z[1:] + z[:-1])
            z_edges[0]    = z[0]  - (z_edges[1] - z[0])
            z_edges[-1]   = z[-1] + (z[-1] - z_edges[-2])
            y_edges = z_edges

        # Draw rows (tell the drawer whether we’re annual or monthly)
        mappable1 = self._draw_row(
            axes_row1, start_year, depth_ref, per_list, labels,
            vmin_row1, vmax_row1,
            row_title=title_top,
            annual=annual_first,
            xtick_interval=xtick_interval,
            vlines_per_panel=vlines_per_exp,   # e.g., [[25,100],[60],[],[10,50,90]]
            vline_mode="year",                 # or "index"
            vline_kwargs=vline_kwargs,
            vline_spin=vline_spin,
            fig_idx=0,
            cmap=cmap,
            fontsize=fontsize,
        )
        
        if draw_cumulative:
            mappable2 = self._draw_row(
                axes_row2, start_year, depth_ref, cum_list, labels,
                vmin_row2, vmax_row2,
                row_title=title_bottom,
                annual=annual_first,
                xtick_interval=xtick_interval,
                vlines_per_panel=vlines_per_exp,
                vline_mode="year",
                vline_kwargs=vline_kwargs,
                vline_spin=vline_spin,
                fig_idx=3,
                cmap=cmap,
                fontsize=fontsize,
            )

        # --- Colorbars (one per row)
        units1 = per_list[0].attrs.get("units", to_units or "")
        cb1 = fig.colorbar(mappable1, cax=cax1, orientation="vertical")
        cb1.ax.set_ylabel(f"{'Δ' if anomaly else ''}OHC per layer [{units1}]", fontsize=fontsize*0.9)
        cb1.ax.tick_params(labelsize=fontsize*0.95)
        #pos = cax1.get_position()
        #cax1.set_position([
        #    pos.x0 +0.1,  # move left (try 0.01–0.03)
        #    pos.y0 +0.035,
        #    pos.width,
        #    pos.height*1.1
        #])
        
        if draw_cumulative:
            units2 = cum_list[0].attrs.get("units", to_units or "")
            cb2 = fig.colorbar(mappable2, cax=cax2, orientation="vertical")
            cb2.ax.set_ylabel(f"{'Δ' if anomaly else ''}OHC cumulative [{units2}]", fontsize=fontsize*0.9)
            cb2.ax.tick_params(labelsize=fontsize*0.95)
            #pos = cax2.get_position()
            #cax2.set_position([
            #    pos.x0 +0.1,  # move left (try 0.01–0.03)
            #    pos.y0 +0.035,
            #    pos.width,
            #    pos.height*1.1
            #])
        
        # --- Super title & layout
        if super_title:
            fig.suptitle(super_title, fontsize=fontsize)

        #fig.subplots_adjust(left=0.06, right=0.98, top=0.93, bottom=0.10, wspace=0.08, hspace=0.18)

        if outname:
            plt.savefig(outname, bbox_inches="tight")
        if draw_cumulative:
            return fig, (axes_row1, axes_row2)
        else:
            return fig, (axes_row1, axes_row1)


class OCNEnsoSpectrum:
    """
    ENSO (Niño-3.4) power spectrum for MPAS/E3SM timeseries:
      - Customer default: Welch PSD (cycles/month), AR(1) red-noise, 95/99% χ².
      - Optional: MPAS-style periodogram + Daniell smoothing.

    Returns/plots spectra in units of °C² per (cycles·month⁻¹).
    """

    # ---------- Construction ----------
    def __init__(
        self,
        root_dir,
        file_key="mpasTimeSeriesOcean",
        var_name="timeMonthly_avg_avgValueWithinOceanRegion_avgSurfaceTemperature",
        region_dim="nOceanRegions",
        region_index=5,                 # 5 is common for Niño3.4 in mpas-analysis outputs
        region_name=None,
        region_names_var="regionNames",
        engine="netcdf4",
        chunks=None,
        calendar="noleap",
        time_origin="0001-01-01 00:00:00"
    ):
        self.root_dir = root_dir
        self.file_key = file_key
        self.var_name = var_name
        self.region_dim = region_dim
        self.region_index = None if region_index is None else int(region_index)
        self.region_name = None if region_name is None else str(region_name)
        self.region_names_var = None if region_names_var in (None, "None") else str(region_names_var)
        self.engine = engine
        self.chunks = chunks
        self.calendar = calendar
        self.time_origin = time_origin
        if self.region_names_var is None:
            raise KeyError("Region names requested but `region_names_var` is None. "
                           "Pass e.g. region_names_var='regionNames' or select by index.")

    # ---------- File & time helpers ----------
    def _find_mpas_timeseries_files(self, exp, sub_dir):
        base = os.path.join(self.root_dir, exp, sub_dir)
        single = glob.glob(os.path.join(base, f"{self.file_key}.nc"))
        if single:
            return single
        segs = glob.glob(os.path.join(base, f"{self.file_key}_*.nc"))
        if not segs:
            return []
        def start_stamp(path):
            m = re.search(rf"{re.escape(self.file_key)}_(\d{{6}})-(\d{{6}})\.nc$", os.path.basename(path))
            return int(m.group(1)) if m else 0
        segs.sort(key=start_stamp)
        return segs

    def _open_dataset(self, files):
        xr_opts = dict(engine=self.engine, chunks=self.chunks, decode_times=True)
        if len(files) == 1:
            return xr.open_dataset(files[0], **xr_opts)

        try:
            return xr.open_mfdataset(
                files, combine="by_coords", parallel=True,
                data_vars="minimal", coords="minimal", compat="override", **xr_opts
            )
        except Exception:
            dsets = [xr.open_dataset(fp, **xr_opts) for fp in files]
            # remove global attrs that sometimes collide
            dsets = [ds_.assign_attrs({}) for ds_ in dsets]
            # pick a time dim from a probe variable
            probe = dsets[0][self.var_name] if self.var_name in dsets[0] else list(dsets[0].data_vars.values())[0]
            tdim = self._guess_time_dim_da(probe)
            out = xr.concat(dsets, dim=tdim, data_vars="minimal", coords="minimal", compat="override")
            for ds_ in dsets:
                try: ds_.close()
                except Exception: pass
            return out

    @staticmethod
    def _guess_time_dim_da(da):
        for d in da.dims:
            if d.lower().startswith("time"):
                return d
        return da.dims[0]

    @staticmethod
    def _decode_time_from_Time_var(ds):
        if "Time" not in ds: return None
        da = ds["Time"]
        units = da.attrs.get("units", None)
        cal   = da.attrs.get("calendar", "noleap")
        if units is None: return None
        return np.asarray(xr.coding.times.decode_cf_datetime(da.values, units, calendar=cal, use_cftime=True))

    @staticmethod
    def _decode_time_from_xtime_chars(ds, name="xtime_startMonthly"):
        if name not in ds: return None
        a = ds[name].values
        s = ["".join(row.astype(str)).strip().strip("\x00").replace("\x00", "") for row in a]
        out = []
        for txt in s:
            txt = txt.replace("T", "_").replace(" ", "_")
            try:
                y = int(txt[0:4]); m = int(txt[5:7]); d = int(txt[8:10])
                hh = int(txt[11:13]); mm = int(txt[14:16]); ss = int(txt[17:19])
                out.append(cftime.DatetimeNoLeap(y, m, d, hh, mm, ss))
            except Exception:
                try:
                    y = int(txt[0:4]); m = int(txt[5:7]); d = int(txt[8:10])
                    out.append(cftime.DatetimeNoLeap(y, m, d))
                except Exception:
                    return None
        return np.asarray(out)

    def load_nino34_series_from_file(
        self,
        path,
        var_candidates=("nino34", "nino3_4", "sst", "value"),
        already_anomalies=False,
        decode_times=True,
        chunks=None,
        return_time=True,
        nino_lat_bounds=(-5.0, 5.0),
        nino_lon_bounds=(190.0, 240.0),   # 170W–120W in 0–360
    ):
        """
        Load a precomputed Niño-3.4 series from a NetCDF file.
        Handles both:
          - 1x1 grid (lat=1, lon=1): select [0,0] and squeeze to 1-D over time
          - full 2D fields: area-weight average over Niño-3.4 box

        Returns (series, time) by default.
        """

        ds = xr.open_dataset(path, decode_times=decode_times, chunks=chunks)
        print(ds)  # optional, for your current debugging

        # --- detect time / lat / lon names
        tdim = next((d for d in ds.dims if d.lower().startswith("time") or d.lower() in ("t","month")), None)
        lat_dim = next((nm for nm in ("lat","latitude","Latitude","y") if nm in ds.coords or nm in ds.dims), None)
        lon_dim = next((nm for nm in ("lon","longitude","Longitude","x") if nm in ds.coords or nm in ds.dims), None)
        if tdim is None:
            raise KeyError("Could not find a time dimension")
        if lat_dim is None or lon_dim is None:
            lat_dim = lat_dim or next((d for d in ds.dims if "lat" in d.lower()), None)
            lon_dim = lon_dim or next((d for d in ds.dims if "lon" in d.lower()), None)
        if lat_dim is None or lon_dim is None:
            # If still missing, we’ll try to fall back to pure 1-D variable detection below
            pass

        # --- pick the variable
        varname = None
        for v in var_candidates:
            if v in ds:
                varname = v
                break
        if varname is None:
            # choose any data_var that depends on time (and possibly lat/lon)
            candidates = [v for v in ds.data_vars if tdim in ds[v].dims]
            if not candidates:
                raise KeyError("No variable with a time dimension found.")
            # Prefer SST-like names if multiple
            def score(v):
                s = v.lower()
                return int("nino" in s) + int("sst" in s) + int(s.endswith("value"))
            candidates.sort(key=score, reverse=True)
            varname = candidates[0]

        da = ds[varname].astype("float64")

        # --- if lat/lon exist and are size-1, select and drop them
        if lat_dim in da.dims and lon_dim in da.dims:
            nlat = int(da.sizes.get(lat_dim, 0))
            nlon = int(da.sizes.get(lon_dim, 0))
            if nlat == 1 and nlon == 1:
                da = da.isel({lat_dim: 0, lon_dim: 0}).drop_vars(
                    [lat_dim, lon_dim], errors="ignore"
                )
            else:
                # if grid bigger than 1x1, do Niño-3.4 area-mean (cos(lat) weights)
                # ensure lon in 0..360
                if lon_dim in ds.coords or lon_dim in ds:
                    ds = ds.assign_coords({lon_dim: (ds[lon_dim] % 360)})
                da = ds[varname].astype("float64")
                lat = ds[lat_dim]
                lon = ds[lon_dim]
                box = da.where(
                    (lat >= nino_lat_bounds[0]) & (lat <= nino_lat_bounds[1]) &
                    ((lon % 360) >= nino_lon_bounds[0]) & ((lon % 360) <= nino_lon_bounds[1]),
                    drop=True
                )
                if box.sizes.get(lat_dim, 0) == 0 or box.sizes.get(lon_dim, 0) == 0:
                    raise ValueError("Niño-3.4 selection produced empty domain; check bounds/coords.")
                w = np.cos(np.deg2rad(box[lat_dim]))
                # broadcast weights to (lat,lon)
                while w.ndim < box[lat_dim].ndim + box[lon_dim].ndim:
                    w = w.broadcast_like(box.isel({tdim: 0}, drop=True))
                w = w / w.mean()
                da = (box * w).mean(dim=[lat_dim, lon_dim], skipna=True)

        # --- ensure 1-D over time
        if da.ndim != 1 or tdim not in da.dims:
            # squeeze any leftover singleton dims, then require 1-D over time
            for d in list(da.dims):
                if d != tdim and da.sizes.get(d, 1) == 1:
                    da = da.isel({d: 0})
            if da.ndim != 1 or tdim not in da.dims:
                raise ValueError(f"Expected 1-D over time after reduction, got dims: {da.dims}")

        # --- values and anomalies
        series = np.asarray(da.values, float).reshape(-1)
        if not already_anomalies:
            series = self.remove_monthly_seasonal_cycle(series)

        # --- time vector (if available)
        tvals = None
        if return_time:
            try:
                tvals = np.asarray(ds[tdim].values)
            except Exception:
                tvals = None

        # log
        if tvals is not None and len(tvals) == len(series):
            print(f"[INFO] HadSST series window: {tvals[0]} → {tvals[-1]}  |  {len(series)} months (~{len(series)/12:.1f} yrs)")
        else:
            print(f"[INFO] HadSST series length: {len(series)} months (~{len(series)/12:.1f} yrs)")

        return (series, tvals) if return_time else series


    def compute_reference_psd_from_file(
        self,
        path,
        year_token=None,       # e.g., "1870-1969" or "1950-2020"
        last_years=None,       # e.g., 50 (overrides year_token if set)
        already_anomalies=False,
        method="Welch",        # "Welch" | "mpas"
        nperseg=256,
        noverlap=128,
        detrend="linear",
        window="hann",
        normalize=False,
        cvdp_kwargs=None,
        label="HadSST",
    ):
        """
        Compute a reference Niño-3.4 spectrum from a precomputed series file.

        Returns
        -------
        dict : {label: {...}}  # schema matches compute_psd_with_rednoise
        """
        import numpy as np
        from scipy import stats

        series, t = self.load_nino34_series_from_file(
            path,
            already_anomalies=already_anomalies,
            return_time=True
        )

        # Choose window
        if last_years is not None and last_years > 0:
            n_last = int(last_years * 12)
            series = series[-n_last:]
            t = t[-n_last:] if t is not None else None
        elif year_token is not None and t is not None:
            try:
                a, b = [int(s) for s in str(year_token).split("-")]
                years = np.array([getattr(tt, "year", np.nan) for tt in t], dtype=float)
                mask = np.isfinite(years) & (years >= a) & (years <= b)
                if mask.any():
                    series = series[mask]
                    t = t[mask]
            except Exception:
                pass  # fall back to full series

        cvdp_kwargs = cvdp_kwargs or {}
        # remove linear drift so rho/var (red noise) match the detrended PSD
        series = signal.detrend(series, type="linear")

        if method.lower() == "mpas":
            d = self._compute_nino34_spectra_mpas(series, **cvdp_kwargs)
            out = {
                "period_years": d["period_years"],
                "psd": d["psd"],
                "red": d["red"],
                "thr95": d["thr95"],
                "thr99": d["thr99"],
                "rho": d["rho"],
                "segments": None,
                "dof": d["dof"],
            }
        else:
            f_cpm, P, nseg, nov, N = self.welch_psd_monthly(
                series, fs=1.0, nperseg=nperseg, noverlap=noverlap,
                detrend=detrend, window=window, normalize=normalize
            )
            rho = self._lag1_autocorr(series)
            var = float(np.nanvar(series, ddof=1)) if np.isfinite(series).any() else 0.0
            Sred = self.ar1_red_noise_spectrum_cpm(f_cpm, rho=rho, var=var)

            step = max(1, nseg - nov)
            m = 1 + max(0, (N - nseg) // step)
            dof = 2 * m
            q95 = stats.chi2.interval(0.95, dof)[1] / dof
            q99 = stats.chi2.interval(0.99, dof)[1] / dof
            thr95 = Sred * q95
            thr99 = Sred * q99

            period_years = 1.0 / (f_cpm * 12.0)
            idx = np.argsort(period_years)
            out = {
                "period_years": period_years[idx],
                "psd": P[idx],
                "red": Sred[idx],
                "thr95": thr95[idx],
                "thr99": thr99[idx],
                "rho": rho,
                "segments": m,
                "dof": dof,
            }

        return {label: out}

    def _construct_monthly_time(self, n):
        if self.time_origin is None: return None
        y = int(self.time_origin[0:4]); m = int(self.time_origin[5:7]); d = int(self.time_origin[8:10])
        hh = int(self.time_origin[11:13]) if len(self.time_origin) >= 13 else 0
        mm = int(self.time_origin[14:16]) if len(self.time_origin) >= 16 else 0
        ss = int(self.time_origin[17:19]) if len(self.time_origin) >= 19 else 0
        t0 = cftime.DatetimeNoLeap(y, m, d, hh, mm, ss)
        return xr.cftime_range(start=t0, periods=int(n), freq="MS")

    @staticmethod
    def _slice_indices_from_year_token(time_index, year_token):
        a, b = [int(t.strip()) for t in year_token.split("-")]
        years = np.array([t.year for t in time_index], dtype=int)
        mask = (years >= a) & (years <= b)
        if not mask.any(): return None
        idx = np.nonzero(mask)[0]
        return int(idx[0]), int(idx[-1] + 1)

    @staticmethod
    def _find_time_slice(year_token):
        a, b = [int(t.strip()) for t in year_token.split("-")]
        if a > b: a, b = b, a
        m_start = (a - 1) * 12
        m_end   = b * 12
        return m_start, m_end

    def _infer_region_dim(self, da, tdim):
        if da.ndim <= 1: return None
        if self.region_dim and self.region_dim in da.dims: return self.region_dim
        for d in da.dims:
            if d != tdim and ("region" in d.lower() or "nregion" in d.lower()):
                return d
        for d in da.dims:
            if d != tdim: return d
        return None

    def _decode_region_names(self, ds, region_dim):
        if self.region_names_var in ds:
            arr = ds[self.region_names_var].values
            if arr.ndim == 2:
                names = ["".join(row.astype(str)).strip().strip("\x00").replace("\x00","") for row in arr]
                names = [n if n else f"{region_dim}[{i}]" for i,n in enumerate(names)]
                if any(ch.isalpha() for ch in "".join(names)): return names
            if arr.ndim == 1:
                out = []
                for i, v in enumerate(arr):
                    s = v.decode("utf-8","ignore").strip() if isinstance(v, bytes) else str(v).strip()
                    out.append(s if s else f"{region_dim}[{i}]")
                if any(ch.isalpha() for ch in "".join(out)): return out
        if region_dim in ds.coords:
            arr = ds.coords[region_dim].values
            if arr.ndim == 1 and arr.size > 0:
                out = []
                for i, v in enumerate(arr):
                    s = v.decode("utf-8","ignore").strip() if isinstance(v, bytes) else str(v).strip()
                    out.append(s if s else f"{region_dim}[{i}]")
                if any(ch.isalpha() for ch in "".join(out)): return out
        return None

    def _region_index_from_name(self, ds, region_dim, target_name):
        names = self._decode_region_names(ds, region_dim)
        if not names:
            raise KeyError(f"Region names not found (looked for '{self.region_names_var or 'regionNames'}' or coord on '{region_dim}').")
        tgt = str(target_name).strip().lower()
        for i, n in enumerate(names):
            if n.lower() == tgt: return i
        def norm(s): return "".join(ch for ch in s.lower() if ch.isalnum())
        tgt2 = norm(tgt)
        for i, n in enumerate(names):
            if norm(n) == tgt2: return i
        raise KeyError(f"Region name '{target_name}' not found. Available: {names}")

    # ---------- Load monthly series (with diagnostics) ----------
    def load_monthly_series(self, exp, sub_dir, year_token, return_time: bool = False):
        files = self._find_mpas_timeseries_files(exp, sub_dir)
        if not files:
            search_root = os.path.join(self.root_dir, exp)
            candidates = glob.glob(os.path.join(
                search_root, "post", "analysis", "mpas_analysis", "ts_*_climo_*",
                "**", f"{self.file_key}*.nc"
            ), recursive=True)
            if candidates:
                sub_dir = os.path.relpath(os.path.dirname(candidates[0]), start=os.path.join(self.root_dir, exp))
                files = self._find_mpas_timeseries_files(exp, sub_dir)
        if not files:
            raise FileNotFoundError(f"No '{self.file_key}*.nc' for exp='{exp}' under subdir='{sub_dir}'")

        ds = self._open_dataset(files)
        if self.var_name not in ds:
            raise KeyError(f"Variable '{self.var_name}' not found in {files[0]}")
        da = ds[self.var_name]
        tdim = self._guess_time_dim_da(da)
        n = da.sizes[tdim]

        time_idx = self._decode_time_from_Time_var(ds)
        if time_idx is None:
            ti_start = self._decode_time_from_xtime_chars(ds, "xtime_startMonthly")
            ti_end   = self._decode_time_from_xtime_chars(ds, "xtime_endMonthly")
            time_idx = ti_start if (ti_start is not None and len(ti_start) == n) else (
                       ti_end   if (ti_end   is not None and len(ti_end)   == n) else None)
        if time_idx is None:
            time_idx = self._construct_monthly_time(n)

        if time_idx is not None and len(time_idx) == n:
            se = self._slice_indices_from_year_token(time_idx, year_token)
            use_start, use_end = (se if se is not None else (0, n))
        else:
            m_start, m_end = self._find_time_slice(year_token)
            if m_start >= n or m_end > n or (m_end - m_start) <= 0:
                use_start, use_end = 0, n
            else:
                use_start, use_end = m_start, m_end

        use_start = max(0, min(use_start, n))
        use_end   = max(use_start, min(use_end, n))  # never exceed n and keep non-negative length

        da = da.isel({tdim: slice(use_start, use_end)}).astype("float64")

        # Diagnostic: selected file-slice time span
        if time_idx is not None:
            t0 = str(time_idx[use_start])
            t1 = str(time_idx[use_end-1]) if use_end-1 < len(time_idx) else str(time_idx[-1])
            print(f"[INFO] file slice used: idx {use_start}:{use_end}  ({t0} → {t1}),  {use_end-use_start} months")
        else:
            print(f"[INFO] file slice used: idx {use_start}:{use_end}  ({use_end-use_start} months; no calendar index)")

        region_dim = self._infer_region_dim(da, tdim)
        if region_dim is not None and region_dim in da.dims:
            if self.region_name is not None:
                idx = self._region_index_from_name(ds, region_dim, self.region_name)
                da = da.isel({region_dim: idx})
            elif self.region_index is not None:
                da = da.isel({region_dim: int(self.region_index)})
            else:
                da = da.isel({region_dim: 0})

        # squeeze to 1-D time
        while da.ndim > 1:
            for d in list(da.dims):
                if d != tdim and da.sizes[d] == 1:
                    da = da.squeeze(d)
                    break
            else:
                for d in list(da.dims):
                    if d != tdim:
                        da = da.mean(dim=d, skipna=True)
                        break

        series = np.asarray(da.values, dtype=float).reshape(-1)
        if not np.isfinite(series).any():
            raise ValueError(f"All values are NaN for exp '{exp}', var '{self.var_name}'")

        if return_time:
            if time_idx is not None and len(time_idx) == n:
                return series, np.asarray(time_idx[use_start:use_end])
            else:
                return series, None
        return series

    # ---------- ENSO helpers ----------
    @staticmethod
    def last_n_years_monthly(series, n_years, samples_per_year=12):
        n_last = int(n_years * samples_per_year)
        return series.copy() if len(series) <= n_last else series[-n_last:]

    @staticmethod
    def remove_monthly_seasonal_cycle(x):
        x = np.asarray(x, float)
        n = len(x)
        clim = np.full(12, np.nan)
        for m in range(12):
            vals = x[m::12]
            clim[m] = np.nanmean(vals) if np.isfinite(vals).any() else np.nan
        anom = x.copy()
        for i in range(n):
            m = i % 12
            if np.isfinite(anom[i]) and np.isfinite(clim[m]):
                anom[i] -= clim[m]
            else:
                anom[i] = np.nan
        if np.isnan(anom).any():
            idx = np.arange(n)
            good = np.isfinite(anom)
            if good.sum() >= 2:
                anom = np.interp(idx, idx[good], anom[good])
        return anom

    # --- Lag-1 autocorrelation ---
    @staticmethod
    def _lag1_autocorr(x):
        x = np.asarray(x, float)
        x = x[np.isfinite(x)]
        if len(x) < 3: return 0.0
        x = x - np.nanmean(x)
        return float(np.corrcoef(x[1:], x[:-1])[0,1])

    # --- AR(1) red noise in cycles/month space ---
    @staticmethod
    def ar1_red_noise_spectrum_cpm(f_cpm, rho, var):
        # f_cpm in cycles/month; returns ONE-SIDED °C² per (cycles/month),
        # consistent with signal.welch(return_onesided=True): ∫_0^0.5 S df = var
        return 2.0 * var * (1.0 - rho**2) / (1.0 + rho**2 - 2.0 * rho * np.cos(2.0 * np.pi * f_cpm))

    # --- Welch PSD, cycles/month ---
    @staticmethod
    def welch_psd_monthly(series_anom, fs=1.0, nperseg=256, noverlap=128,
                          detrend='linear', window='hann', normalize=False):
        """
        Welch PSD for monthly anomalies with fs=1.0 sample/month.
        Returns f_cpm (cycles/month) and PSD in °C² / (cycles·month⁻1).
        """
        x = np.asarray(series_anom, float)
        N = len(x)
        nperseg = min(int(nperseg), N)
        noverlap = min(int(noverlap), nperseg // 2)
        f, Pxx = signal.welch(
            x, fs=fs, window=window, nperseg=nperseg, noverlap=noverlap,
            detrend=detrend, return_onesided=True, scaling='density'
        )
        keep = f > 0
        f_cpm = f[keep]        # cycles/month
        P = Pxx[keep]          # °C²/(cycles/month)
        if normalize and np.trapezoid(P, f_cpm) > 0:
            P = P / np.trapezoid(P, f_cpm)
        return f_cpm, P, nperseg, noverlap, N

    # --- MPAS-style (periodogram + Daniell smoothing) ---
    @staticmethod
    def _running_mean_weights(x, wgts):
        """Centered weighted running mean with edge reflection (Daniell-style)."""
        x = np.asarray(x, float)
        w = np.asarray(wgts, float)
        half = len(w) // 2
        pad_left  = x[1:half+1][::-1] if half > 0 else np.array([])
        pad_right = x[-half-1:-1][::-1] if half > 0 else np.array([])
        xp = np.concatenate([pad_left, x, pad_right])
        y = np.convolve(xp, w, mode="valid")   # length == len(x)
        return y

    def _compute_nino34_spectra_mpas(self, nino34Index,
                                     sec_per_month=30.0*24*3600.0,
                                     tapcoef=2.661, tukey_alpha=0.10):
        x = nino34Index.values if hasattr(nino34Index, "values") else np.asarray(nino34Index, float)
        x = np.asarray(x, float)
        N = x.size
        if N < 16:
            raise ValueError("Time series is too short for spectral analysis (<16 samples).")

        # Tukey window (robust)
        try:
            w = windows.tukey(N, alpha=tukey_alpha)
        except Exception:
            try:
                w = signal.windows.tukey(N, alpha=tukey_alpha)
            except Exception:
                w = np.ones(N, float)

        # Periodogram (Hz)
        f_hz, Pxx_hz = signal.periodogram(w * x, fs=1.0/sec_per_month, detrend=False, scaling="density")
        
        # Modified Daniell smoothing
        nwts = max(1, int(round(7 * N / 1200)))
        if nwts % 2 == 0:
            nwts += 1
        wgts = np.ones(nwts, float); wgts[0] = wgts[-1] = 0.5; wgts /= wgts.sum()

        PxxSmooth_full = self._running_mean_weights(Pxx_hz, wgts) / sec_per_month  # °C²/(cycles/month)
        half = nwts // 2
        
        if len(f_hz) <= nwts:
            # Degenerate case: skip smoothing
            f_hz_s = f_hz
            PxxSmooth = Pxx_hz / sec_per_month

        # --- FIX: trim smoothed spectrum and frequency *consistently* ---
        if half > 0 and 2*half < len(f_hz):
            f_hz_s = f_hz[half:-half]
            PxxSmooth = PxxSmooth_full[half:-half]     # <-- added
        else:
            f_hz_s = f_hz
            PxxSmooth = PxxSmooth_full                 # <-- added

        # Convert to cycles/month and period in years
        f_cpm = f_hz_s * sec_per_month
        keep = f_cpm > 0
        f_cpm     = f_cpm[keep]
        PxxSmooth = PxxSmooth[keep]                    # same mask length as f_cpm
        period_years = 1.0 / (f_cpm * 12.0)

        # AR(1) red-noise (cycles/month) scaled to smoothed spectrum
        r  = self._lag1_autocorr(x)
        r2 = 2.0 * r; rsq = r * r
        temp2 = r2 * np.cos(2.0 * np.pi * f_cpm)
        mkov = 1.0 / (1.0 + rsq - temp2)
        scale = (PxxSmooth.sum()) / (mkov.sum()) if mkov.sum() != 0 else 1.0
        red = mkov * scale

        # DOF & thresholds (MPAS)
        dof = 2.0 / (tapcoef * np.sum(wgts**2))
        q95  = stats.chi2.interval(0.95, dof)[1] / dof
        q99  = stats.chi2.interval(0.99, dof)[1] / dof
        thr95 = red * q95
        thr99 = red * q99

        return {
            "period_years": period_years,
            "psd": PxxSmooth,
            "red": red,
            "thr95": thr95,
            "thr99": thr99,
            "rho": r,
            "dof": dof
        }

    def summarize_enso_metrics(self, psd_dict, band=(2.0, 7.0), print_table=False, return_df=True):
        """
        Summarize ENSO metrics from PSDs.

        Parameters
        ----------
        psd_dict : dict
            Output from compute_psd_with_rednoise(). Each entry must contain:
            - "period_years": array-like, period in years (ascending or descending)
            - "psd": array-like, power spectral density [°C² / (cycles·month⁻¹)]
            - optionally "rho" (lag-1 autocorr) and "segments" (Welch m) or "dof" (MPAS)
        band : tuple(float,float)
            ENSO period band in years, e.g. (2, 7).
        print_table : bool
            If True, prints a neat text table.
        return_df : bool
            If True, returns a pandas.DataFrame; else returns a list of dict rows.

        Returns
        -------
        pandas.DataFrame or list[dict]
            Columns: experiment, ENSO_power, peak_period_yr, peak_power, rho, segments
            - ENSO_power is the integrated variance (°C²) within the period band.
        """
        import numpy as np
        try:
            import pandas as pd
        except Exception:
            pd = None

        lo_yr, hi_yr = float(band[0]), float(band[1])
        if lo_yr > hi_yr:
            lo_yr, hi_yr = hi_yr, lo_yr

        rows = []
        for exp, d in psd_dict.items():
            period = np.asarray(d["period_years"], float)
            psd    = np.asarray(d["psd"], float)

            # Convert period→frequency in cycles/month for proper integration of density
            # f [cycles/month] = 1 / (period[years] * 12)
            f_cpm = 1.0 / (period * 12.0)

            # Band edges in frequency space (note reversal: shorter period => higher freq)
            f_lo = 1.0 / (hi_yr * 12.0)
            f_hi = 1.0 / (lo_yr * 12.0)

            mask = np.isfinite(f_cpm) & np.isfinite(psd) & (f_cpm >= f_lo) & (f_cpm <= f_hi)

            if np.count_nonzero(mask) >= 2:
                # Integrate PSD over frequency to get variance contribution (°C²)
                # sort by ascending frequency (Welch output is sorted by period)
                o = np.argsort(f_cpm[mask])
                enso_power = float(np.trapezoid(psd[mask][o], f_cpm[mask][o]))
                # Peak within band
                jloc = np.argmax(psd[mask])
                f_peak = float(f_cpm[mask][jloc])
                peak_power = float(psd[mask][jloc])
                peak_period_yr = float(1.0 / (f_peak * 12.0))
            else:
                enso_power = np.nan
                peak_power = np.nan
                peak_period_yr = np.nan

            rho = float(d.get("rho", np.nan))
            # Prefer explicit Welch segment count if present; otherwise, map DOF≈2m
            if "segments" in d and d["segments"] is not None:
                segments = int(d["segments"])
            elif "dof" in d and d["dof"] is not None:
                segments = int(round(float(d["dof"]) / 2.0))
            else:
                segments = np.nan

            rows.append({
                "experiment": exp,
                "ENSO_power": enso_power,      # °C²
                "peak_period_yr": peak_period_yr,
                "peak_power": peak_power,      # °C² / (cycles·month⁻1)
                "rho": rho,
                "segments": segments
            })

        if print_table:
            if pd is not None:
                df_print = pd.DataFrame(rows)
                # Nice ordering
                df_print = df_print[["experiment", "ENSO_power", "peak_period_yr", "peak_power", "rho", "segments"]]
                # Pretty print with limited precision
                with pd.option_context('display.max_rows', None, 'display.max_colwidth', None,
                                       'display.float_format', lambda v: f"{v:.3f}"):
                    print(df_print.to_string(index=False))
            else:
                # Fallback plain print
                for r in rows:
                    print(f"{r['experiment']:>40s}  ENSO_power={r['ENSO_power']:.3f}  "
                          f"peak_period_yr={r['peak_period_yr']:.3f}  peak_power={r['peak_power']:.3f}  "
                          f"rho={r['rho']:.3f}  segments={r['segments']}")

        return pd.DataFrame(rows) if (return_df and pd is not None) else rows

    # ---------- Batch compute (Welch default; MPAS optional) ----------
    def compute_psd_with_rednoise(self, exp_to_subdir, year_token, last_years=50,
                                  nperseg=256, noverlap=128,
                                  detrend='linear', window='hann', normalize=False,
                                  method="Welch",          # "Welch" | "mpas"
                                  cvdp_kwargs=None):
        """
        Returns dict per experiment:
          {"period_years","psd","red","thr95","thr99","rho","segments","dof"}
        Units: °C² per (cycles·month⁻1). Period in years.
        """
        out = {}
        cvdp_kwargs = cvdp_kwargs or {}

        for exp, sub_dir in exp_to_subdir.items():
            try:
                series, tsel = self.load_monthly_series(exp, sub_dir, year_token, return_time=True)

                # last-N years (calendar-aware diagnostic)
                n_last = int(last_years * 12)
                if len(series) <= n_last:
                    series_last = series.copy()
                    tsel_last = tsel.copy() if tsel is not None else None
                else:
                    series_last = series[-n_last:]
                    tsel_last = tsel[-n_last:] if tsel is not None else None

                n_months = len(series_last); yrs_actual = n_months / 12.0
                if not hasattr(self, "_last_window_years"):
                    self._last_window_years = {}
                self._last_window_years[exp] = (
                    (int(tsel_last[0].year), int(tsel_last[-1].year))
                    if tsel_last is not None and len(tsel_last) > 0 else None
                )
                if tsel_last is not None and len(tsel_last) > 0:
                    print(f"[INFO] PSD window ({exp}): {tsel_last[0]} → {tsel_last[-1]}  |  {n_months} months (~{yrs_actual:.1f} yrs)")
                else:
                    print(f"[INFO] PSD window ({exp}): last {n_months} samples (~{yrs_actual:.1f} yrs); no calendar index")

                anom = self.remove_monthly_seasonal_cycle(series_last)
                # remove linear drift so rho/var (red noise) match the detrended PSD
                anom = signal.detrend(anom, type="linear")
                # cache last-used anomaly series (optional)
                if not hasattr(self, "_last_series_anom"):
                    self._last_series_anom = {}
                self._last_series_anom[exp] = anom
                
                if method.lower() == "mpas":  # MPAS-style
                    d = self._compute_nino34_spectra_mpas(anom, **cvdp_kwargs)
                    out[exp] = {
                        "period_years": d["period_years"],
                        "psd": d["psd"],
                        "red": d["red"],
                        "thr95": d["thr95"],
                        "thr99": d["thr99"],
                        "rho": d["rho"],
                        "segments": None,
                        "dof": d["dof"],
                    }
                else:  # Welch (default)
                    f_cpm, P, nseg, nov, N = self.welch_psd_monthly(
                        anom, fs=1.0, nperseg=nperseg, noverlap=noverlap,
                        detrend=detrend, window=window, normalize=normalize
                    )
                    rho = self._lag1_autocorr(anom)
                    var = float(np.nanvar(anom, ddof=1)) if np.isfinite(anom).any() else 0.0
                    Sred = self.ar1_red_noise_spectrum_cpm(f_cpm, rho=rho, var=var)

                    # Welch DOF via segments (≈2m)
                    step = max(1, nseg - nov)
                    m = 1 + max(0, (N - nseg) // step)
                    dof = 2 * m
                    q95 = stats.chi2.interval(0.95, dof)[1] / dof
                    q99 = stats.chi2.interval(0.99, dof)[1] / dof
                    thr95 = Sred * q95
                    thr99 = Sred * q99

                    period_years = 1.0 / (f_cpm * 12.0)
                    idx = np.argsort(period_years)
                    out[exp] = {
                        "period_years": period_years[idx],
                        "psd": P[idx],
                        "red": Sred[idx],
                        "thr95": thr95[idx],
                        "thr99": thr99[idx],
                        "rho": rho,
                        "segments": m,
                        "dof": dof,
                    }

            except Exception as e:
                print(f"[WARN] Skipping '{exp}': {type(e).__name__}: {e}")

        skipped = [e for e in exp_to_subdir if e not in out]
        if skipped:
            print(f"[WARN] {len(skipped)}/{len(exp_to_subdir)} experiments skipped: {skipped}")
        return out

    # ---------- Plot panels (1 row, shared legend bottom) ----------
    @staticmethod
    def plot_psd_panels(psd_dict, labels=None, colors=None,
                        xlim=(10, 1), ylim=None,
                        enso_band=(2.0, 7.0),
                        ylabel=r"$^\circ$C$^2$ / cycles mo$^{-1}$",
                        title="Niño3.4 SST Power Spectrum",
                        figsize=(18, 5), fontsize=13,
                        fig_idx=0,
                        out_path=None,
                        legend_at="bottom",   # "top" | "bottom" | None
                        legend_ncol=4,
                        axes=None):
        """
        Plot spectra side-by-side (1 × n) with shared Y and one shared legend.

        axes: existing axes (one per spectrum) to draw into instead of creating a figure; the
              figure layout is then left to the caller (no tight_layout), the shared legend is
              still added at legend_at and every panel gets its own y label.
        """
        n = len(psd_dict)
        if n == 0:
            raise ValueError("psd_dict is empty (no spectra to plot).")

        if isinstance(labels, (list, tuple)) and len(labels) != n:
            raise ValueError(f"labels length ({len(labels)}) does not match panels ({n}).")

        own_fig = axes is None
        if own_fig:
            fig, axes = plt.subplots(1, n, figsize=figsize, sharey=True)
            if n == 1:
                axes = [axes]
        else:
            axes = list(np.ravel(axes))
            if len(axes) != n:
                raise ValueError(f"got {len(axes)} axes for {n} spectra")
            fig = axes[0].figure
            
        shared_handles, shared_labels = None, None
        items = list(psd_dict.items())

        for i, (exp, d) in enumerate(items):
            ax = axes[i]
            lbl = (
                labels[i] if isinstance(labels, (list, tuple))
                else (labels.get(exp, exp) if isinstance(labels, dict) else exp)
            )

            h_spec, = ax.plot(d["period_years"], d["psd"],   color="black", lw=4.0, label="Spectrum")
            h_red,  = ax.plot(d["period_years"], d["red"],   color="red",   lw=4.0, label="Red noise")
            h_95,   = ax.plot(d["period_years"], d["thr95"], color="green", lw=4.0, label="95% thresh")
            h_99,   = ax.plot(d["period_years"], d["thr99"], color="blue",  lw=4.0, label="99% thresh")

            if shared_handles is None:
                shared_handles = [h_spec, h_red, h_95, h_99]
                shared_labels  = [h.get_label() for h in shared_handles]

            ax.set_xlim(*xlim)
            if ylim is not None:
                ax.set_ylim(*ylim)
                
            panel_label = f"({chr(97+i+fig_idx)}) "  # 97 is ASCII 'a'
            ax.set_title(f"{panel_label} {title}", fontsize=fontsize*1.2,loc="left", pad=15)
            ax.set_title(lbl, fontsize=fontsize*1.2,loc="right", pad=15)
            
            ax.grid(True, alpha=0.25)

            # ENSO band shading
            if enso_band is not None:
                lo, hi = enso_band
                xl0, xl1 = (xlim[0], xlim[1])
                lo_shade = max(min(lo, xl0), min(lo, xl1))
                hi_shade = min(max(hi, xl0), max(hi, xl1))
                if hi_shade > lo_shade:
                    ax.axvspan(lo_shade, hi_shade, color="gray", linewidth = 1.5, alpha=0.15, zorder=0)

            if i == 0 or not own_fig:
                ax.set_ylabel(ylabel, fontsize=fontsize)
            ax.set_xlabel("Period (years)", fontsize=fontsize)
            ax.tick_params(axis="x", labelsize=fontsize*0.95)
            ax.tick_params(axis="y", labelsize=fontsize*0.95)

        if legend_at in ("top", "bottom"):
            if legend_at == "top":
                leg = fig.legend(shared_handles, shared_labels,
                                 loc="upper center", bbox_to_anchor=(0.5, 1.02),
                                 ncol=legend_ncol, fontsize=fontsize*0.95,
                                 framealpha=0.6)
                if own_fig:
                    plt.tight_layout(rect=[0.02, 0.02, 1.0, 0.92])
            elif own_fig:  # bottom
                leg = fig.legend(shared_handles, shared_labels,
                                 loc="lower center", bbox_to_anchor=(0.5, -0.06),
                                 ncol=legend_ncol, fontsize=fontsize,
                                 framealpha=0.6)
                plt.tight_layout(rect=[0.02, 0.06, 1.0, 0.93])
            else:  # bottom, just below the given axes (clear of their x labels)
                x0 = min(ax.get_position().x0 for ax in axes)
                x1 = max(ax.get_position().x1 for ax in axes)
                y0 = min(ax.get_position().y0 for ax in axes)
                pad = 2.6 * fontsize / 72.0 / fig.get_figheight()   # ~ tick labels + x label
                leg = fig.legend(shared_handles, shared_labels,
                                 loc="upper center", bbox_to_anchor=(0.5 * (x0 + x1), y0 - pad),
                                 ncol=legend_ncol, fontsize=fontsize,
                                 framealpha=0.6)
                
            for legline in leg.get_lines():
                legline.set_linewidth(4.0)
                
        elif own_fig:
            plt.tight_layout(rect=[0.02, 0.04, 1.0, 0.95])

        if out_path:
            fig.savefig(out_path, dpi=200, bbox_inches="tight")
        return axes
        
    def save_psd_cache(
        self,
        psd_dict: dict,
        out_nc: str,
        *,
        labels=None,
        meta: dict = None,
        series_anom_dict: dict = None,   # optional: exp -> 1D monthly anomalies used for PSD
    ):
        """
        Save PSD outputs (and optionally the anomaly series used) to a NetCDF cache.
    
        psd_dict[exp] must contain:
          period_years, psd, red, thr95, thr99, rho, segments, dof
        """
        if not psd_dict:
            raise ValueError("psd_dict is empty")
    
        exps = list(psd_dict.keys())
    
        # ---- pad spectra to common length (Kmax) ----
        Kmax = max(len(np.asarray(psd_dict[e]["period_years"])) for e in exps)
    
        def _pad(x, L, fill=np.nan):
            x = np.asarray(x, float)
            y = np.full((L,), fill, dtype=float)
            n = min(L, len(x))
            if n > 0:
                y[:n] = x[:n]
            return y
    
        period = np.vstack([_pad(psd_dict[e]["period_years"], Kmax) for e in exps])
        psd    = np.vstack([_pad(psd_dict[e]["psd"],         Kmax) for e in exps])
        red    = np.vstack([_pad(psd_dict[e]["red"],         Kmax) for e in exps])
        thr95  = np.vstack([_pad(psd_dict[e]["thr95"],       Kmax) for e in exps])
        thr99  = np.vstack([_pad(psd_dict[e]["thr99"],       Kmax) for e in exps])
    
        rho      = np.array([float(psd_dict[e].get("rho", np.nan)) for e in exps], dtype=float)
        segments = np.array([psd_dict[e].get("segments", np.nan) for e in exps], dtype=float)
        dof      = np.array([float(psd_dict[e].get("dof", np.nan)) for e in exps], dtype=float)
    
        # labels handling
        if labels is None:
            lbl = [e for e in exps]
        elif isinstance(labels, dict):
            lbl = [labels.get(e, e) for e in exps]
        else:
            lbl = list(labels)
            if len(lbl) != len(exps):
                raise ValueError("labels length does not match number of experiments")
        lbl = np.asarray(lbl, dtype=object)
    
        ds = xr.Dataset(
            data_vars=dict(
                period_years=(("exp", "k"), period),
                psd=(("exp", "k"), psd),
                red=(("exp", "k"), red),
                thr95=(("exp", "k"), thr95),
                thr99=(("exp", "k"), thr99),
                rho=("exp", rho),
                segments=("exp", segments),
                dof=("exp", dof),
                label=("exp", lbl),
            ),
            coords=dict(
                exp=("exp", np.asarray(exps, dtype=object)),
                k=("k", np.arange(Kmax, dtype=int)),
            ),
        )
    
        # ---- optional: also store anomaly series used for PSD ----
        if series_anom_dict is not None and len(series_anom_dict) > 0:
            Tmax = max(len(np.asarray(series_anom_dict.get(e, []))) for e in exps)
            if Tmax > 0:
                series_pad = np.vstack([_pad(series_anom_dict.get(e, np.array([])), Tmax) for e in exps])
                ds["series_anom"] = (("exp", "time"), series_pad)
                ds = ds.assign_coords(time=("time", np.arange(Tmax, dtype=int)))
    
        ds.attrs["cache_meta_json"] = json.dumps(meta or {}, default=str)
        ds.to_netcdf(out_nc)
        return out_nc
    
    
    @staticmethod
    def load_psd_cache(nc_path: str) -> xr.Dataset:
        return xr.open_dataset(nc_path)


def write_annual_flux(src_nc: str, var: str, out_nc: str, out_var: str,
                      region: str = "Global") -> str:
    """
    Day-weighted (noleap) annual means of a monthly flux series from an
    OHUFigureDataCollector cache (e.g. ``Qnet``/``Qsource`` in
    ``data/ohu/*_ts_<exp>.nc``), written as a 1-D (year) file
    that LaggedOHCCorrelation reads. Years are numbered 1..N.
    """
    dpm = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], float)
    with xr.open_dataset(src_nc) as ds:
        monthly = np.asarray(ds[var].values, float)
        units = ds[var].attrs.get("units", "W m-2")
    nyr = monthly.size // 12
    annual = np.average(monthly[:nyr * 12].reshape(nyr, 12), axis=1, weights=dpm)
    da = xr.DataArray(annual, dims=("year",),
                      coords={"year": np.arange(1, nyr + 1, dtype=np.int32)},
                      name=out_var, attrs={"units": units, "long_name": f"Annual mean {var}"})
    out = xr.Dataset({out_var: da}, attrs={
        "Conventions": "CF-1.8",
        "region": region,
        "source": os.path.basename(src_nc),
        "history": f"{datetime.now(timezone.utc):%Y-%m-%d %H:%M:%S}Z: created by write_annual_flux",
    })
    os.makedirs(os.path.dirname(out_nc) or ".", exist_ok=True)
    out.to_netcdf(out_nc, encoding={out_var: dict(dtype="float32", _FillValue=np.float32(np.nan))})
    return out_nc


class LaggedOHCCorrelation:
    """
    Compute lagged correlations between diagnosed net surface heat flux (Qnet)
    and ocean heat content (OHC) anomalies by depth.

    Correlation (Qnet LEADS OHC):
        corr_lag(L) = corr( Qnet(t), OHC(z, t + L) )

    Assumptions:
    - Annual means for both Qnet and OHC.
    - OHC already anomaly (per your files).
    - Uses xarray default engine (change _try_open if you need a specific one).
    """

    def __init__(self,
                 ohc_files: Dict[str, str],
                 qnet_files: Dict[str, str],
                 *,
                 ohc_var: Optional[str] = None,
                 qnet_var: Optional[str] = None,
                 time_name: Optional[str] = None,
                 depth_name: Optional[str] = None,
                 detrend: bool = False,
                 decode_times: bool = False):
        self.ohc_files = ohc_files
        self.qnet_files = qnet_files
        self.ohc_var = ohc_var
        self.qnet_var = qnet_var
        self.time_name = time_name
        self.depth_name = depth_name
        self.detrend = detrend
        self.decode_times = decode_times

        self._time_candidates = ["year", "time", "Time"]
        self._depth_candidates = ["layer", "nVertLevels", "depth", "z"]
        self._qnet_candidates = ["qnet", "netheat", "netflux", "net_heat", "net_flux"]
        self._ohc_candidates = ["ohc"]
        
        self.panel_tags = "lower"      # "lower" | "upper" | None
        self.panel_tag_offset = 0      # start index offset

    # ========================= public API =========================
    def run(self,
            lags: Sequence[int] = (0, 10, 20, 50),
            *,
            save_nc: Optional[str] = None,
            ts_window: tuple[int, int] = None,   # <--- add this
            make_plots: bool = True,
            fig_dir: str = "./lagcorr_figs",
            fig_name: str= "./lag_corr.png",
            vmax: float = 1.0,
            fig_idx: int = 0, 
            also_return_beta: bool = True,
            xlabel: str= None,
            ylabel: str= None,
            xtick_interval: int = 5, 
            row_title: str= None,
            fontsize: int = 18, 
            figsize: tuple=(None,None)
        ):
        
        exps = sorted(set(self.ohc_files).intersection(self.qnet_files))
        if not exps:
            raise ValueError("No overlapping experiment keys between OHC and Qnet file maps.")

        corr_list = []
        beta_list = []
        depth_ref = None

        for exp in exps:
            ds_ohc, ohc_da = self._open_ohc(self.ohc_files[exp])
            time_name = self._get_time_name(ohc_da)
            depth_name = self._get_depth_name(ohc_da)

            # Standardize time to integer years (1..N)
            years_ohc = self._to_year_index(ohc_da[time_name])
            ohc_da = ohc_da.assign_coords({time_name: years_ohc})

            # Qnet
            ds_qnet, qnet = self._open_qnet(self.qnet_files[exp])
            t_q = self._get_time_name(qnet)
            years_qnet = self._to_year_index(qnet[t_q])
            qnet = qnet.assign_coords({t_q: years_qnet})

            # Intersect time
            common_years = np.intersect1d(ohc_da[time_name].values, qnet[t_q].values)
            if common_years.size < 10:
                raise ValueError(f"{exp}: Too few overlapping annual points ({common_years.size}).")
            ohc_da = ohc_da.sel({time_name: common_years})
            qnet   = qnet.sel({t_q: common_years})
            
            # Restrict to ts_window if provided
            if ts_window is not None:
                ymin, ymax = ts_window
                ohc_da = ohc_da.sel({time_name: slice(ymin, ymax)})
                qnet   = qnet.sel({t_q: slice(ymin, ymax)})

            # Optional detrend
            if self.detrend:
                qnet   = self._detrend_1d(qnet, t_q)
                ohc_da = self._detrend_2d_time(ohc_da, time_dim=time_name)

            # Lagged stats
            corr_exp, beta_exp = self._corr_by_depth(
                ohc_da, qnet, time_name, t_q, depth_name, lags, also_beta=also_return_beta
            )
            corr_exp = corr_exp.assign_coords(exp=exp).expand_dims("exp")
            if also_return_beta:
                beta_exp = beta_exp.assign_coords(exp=exp).expand_dims("exp")

            # Unify depth coordinate
            if depth_ref is None:
                depth_ref = self._depth_mid_from(ds_ohc, depth_name)
                corr_exp = corr_exp.assign_coords({depth_name: depth_ref})
                if also_return_beta:
                    beta_exp = beta_exp.assign_coords({depth_name: depth_ref})
            else:
                if corr_exp.sizes[depth_name] != depth_ref.size:
                    corr_exp = corr_exp.interp({depth_name: depth_ref})
                    if also_return_beta:
                        beta_exp = beta_exp.interp({depth_name: depth_ref})

            corr_list.append(corr_exp)
            if also_return_beta:
                beta_list.append(beta_exp)

            # Close datasets
            for ds in (ds_ohc, ds_qnet):
                try:
                    ds.close()
                except Exception:
                    pass

        ds_out = xr.Dataset(
            data_vars=dict(corr=xr.concat(corr_list, dim="exp")),
            coords=dict(
                exp=("exp", np.array(exps, dtype=object)),
                lag=("lag", np.array(list(lags), dtype=int)),
            ),
            attrs=dict(
                description="Lagged correlation between Qnet(t) and OHC(z, t+lag)",
                note="Positive lag means Qnet leads OHC."
            )
        )
        if also_return_beta:
            ds_out["beta"] = xr.concat(beta_list, dim="exp")
            ds_out["beta"].attrs["description"] = "Standardized regression slope of OHC(z,t+lag) on Qnet(t)"

        if save_nc:
            os.makedirs(os.path.dirname(os.path.abspath(save_nc)), exist_ok=True)
            ds_out.to_netcdf(save_nc)

        if make_plots:
            os.makedirs(fig_dir, exist_ok=True)
            # resolve full path (respect absolute or already-joined names)
            fig_path = fig_name if os.path.isabs(fig_name) else os.path.join(fig_dir, fig_name)
            # enforce left→right sequence
            desired_order = ["FC", "Full-CPL", "FC-FOSI-S1", "FC-FOSI-S2"]
            exps_in_ds = [e for e in desired_order if e in list(ds_out.exp.values)]
            self._plot_multi_exp(
                ds_out,
                exps=exps_in_ds,
                depth_name=self.depth_name or "depth_mid",
                vmax=vmax,
                fig_idx=fig_idx, 
                fig_path=fig_path,
                xlabel=xlabel,
                ylabel=ylabel,
                row_title=row_title,
                fontsize=fontsize,
                xtick_interval=xtick_interval,
                figsize_per_panel=figsize, 
            )
        return ds_out

    # ========================= internals =========================
    def _try_open(self, path: str) -> xr.Dataset:
        return xr.open_dataset(path, decode_times=self.decode_times)

    def _open_ohc(self, path: str) -> tuple[xr.Dataset, xr.DataArray]:
        ds = self._try_open(path)
        ohc_var = self.ohc_var or self._find_var(ds, self._ohc_candidates)
        if ohc_var is None:
            raise ValueError(f"Could not auto-detect OHC variable in {path}. Set ohc_var=.")
        da = ds[ohc_var]
        tname = self._get_time_name(da)
        dname = self._get_depth_name(da)
        if (tname, dname) != tuple(da.dims):
            da = da.transpose(tname, dname, ...)
        da = da.load()
        return ds, da

    def _open_qnet(self, path: str) -> tuple[xr.Dataset, xr.DataArray]:
        ds = self._try_open(path)
        qvar = self.qnet_var or self._find_var(ds, self._qnet_candidates)
        if qvar is None:
            candidates = []
            for v in ds.data_vars:
                da_v = ds[v]
                attrs = " ".join(str(da_v.attrs.get(a, "")).lower()
                                 for a in ("long_name", "standard_name", "description", "units"))
                if any(k in attrs for k in ("qnet", "net heat", "net surface", "heat flux", "surface heat")):
                    candidates.append(v)
            if not candidates:
                raise ValueError(f"No Qnet-like variable in {path}. Vars: {list(ds.data_vars)}")
            qvar = candidates[0]
        da = ds[qvar]

        tname = None
        for t in self._time_candidates:
            if t in da.dims or t in ds.dims or t in ds.coords:
                tname = t
                break
        if tname is None:
            tname = da.dims[0]

        if da.dims[0] != tname:
            dims_new = (tname,) + tuple(d for d in da.dims if d != tname)
            da = da.transpose(*dims_new)

        da = da.load()

        # Squeeze length-1 non-time dims, collapse bounds/meta dims
        other_dims = [d for d in da.dims if d != tname]
        for d in list(other_dims):
            if da.sizes[d] == 1:
                da = da.squeeze(d, drop=True)
        bounds_like = {"nv", "nbnd", "bnds", "bounds"}
        for d in [d for d in da.dims if d != tname]:
            if d.lower() in bounds_like:
                da = da.mean(d, keep_attrs=True)

        other_dims = [d for d in da.dims if d != tname]
        if len(other_dims) == 1 and da.sizes[other_dims[0]] <= 4:
            da = da.mean(other_dims[0], keep_attrs=True)
        if da.ndim != 1 or any(d != tname for d in da.dims):
            raise RuntimeError(f"Qnet still not 1-D after reductions; dims={da.dims}")
        return ds, da

    # ---- utilities ----
    def _find_var(self, ds: xr.Dataset, keys: Sequence[str]) -> Optional[str]:
        name_hits = []
        for v in ds.data_vars:
            lv = v.lower()
            if any(k in lv for k in keys):
                name_hits.append(v)

        def score(varname):
            da = ds[varname]
            time_like = any(t in da.dims for t in self._time_candidates)
            return (int(time_like), -da.ndim)

        if name_hits:
            name_hits.sort(key=score, reverse=True)
            return name_hits[0]

        attr_hits = []
        for v in ds.data_vars:
            da = ds[v]
            attrs = " ".join(str(da.attrs.get(a, "")).lower()
                             for a in ("long_name", "standard_name", "description"))
            if any(k in attrs for k in keys):
                attr_hits.append(v)
        if attr_hits:
            attr_hits.sort(key=score, reverse=True)
            return attr_hits[0]
        return None

    def _get_time_name(self, da: xr.DataArray) -> str:
        if self.time_name and self.time_name in da.dims:
            return self.time_name
        for t in self._time_candidates:
            if t in da.dims:
                return t
        return da.dims[0]

    def _get_depth_name(self, da: xr.DataArray) -> str:
        if self.depth_name and self.depth_name in da.dims:
            return self.depth_name
        for z in self._depth_candidates:
            if z in da.dims:
                return z
        return da.dims[1] if da.ndim >= 2 else "depth"

    def _to_year_index(self, tcoord: xr.DataArray) -> xr.DataArray:
        vals = tcoord.values
        if np.issubdtype(vals.dtype, np.integer):
            return tcoord.astype(int)
        return xr.DataArray(np.arange(1, vals.shape[0] + 1, dtype=int), dims=tcoord.dims)

    def _detrend_1d(self, da: xr.DataArray, time_dim: str) -> xr.DataArray:
        t = np.arange(da.sizes[time_dim], dtype=float)
        x = da.values.astype(float)
        m, b = np.polyfit(t, x, 1) if np.isfinite(x).sum() >= 3 else (0.0, 0.0)
        return da - (m * t + b)

    def _detrend_2d_time(self, da: xr.DataArray, time_dim: str) -> xr.DataArray:
        t = np.arange(da.sizes[time_dim], dtype=float)
        vals = da.values.astype(float)
        vals_dt = np.empty_like(vals)
        for j in range(vals.shape[1]):
            y = vals[:, j]
            mask = np.isfinite(y)
            if mask.sum() >= 3:
                m, b = np.polyfit(t[mask], y[mask], 1)
                vals_dt[:, j] = y - (m * t + b)
            else:
                vals_dt[:, j] = np.nan
        return xr.DataArray(vals_dt, dims=da.dims, coords=da.coords, attrs=da.attrs)

    # ---- fixed: take arrays explicitly ----
    def corr_at_lag(self, Q: np.ndarray, O: np.ndarray, L: int) -> np.ndarray:
        """
        Q: [time]  (Qnet)
        O: [time, depth] (OHC anomalies)
        L: non-negative lag, years that Q leads O
        returns: r[depth]
        """
        if L < 0:
            raise ValueError("Only non-negative lags supported (Qnet leading OHC).")
        nt = min(Q.shape[0], O.shape[0])
        if L >= nt - 3:
            return np.full(O.shape[1], np.nan)

        x = Q[: nt - L].astype(float)
        y = O[L: nt, :].astype(float)

        xm = np.nanmean(x)
        xs = np.nanstd(x)
        if not np.isfinite(xs) or xs == 0.0:
            return np.full(O.shape[1], np.nan)
        x0 = (x - xm) / xs

        ym = np.nanmean(y, axis=0)
        ys = np.nanstd(y, axis=0)
        valid = (ys > 0) & np.isfinite(ys)
        y0 = np.full_like(y, np.nan)
        y0[:, valid] = (y[:, valid] - ym[valid]) / ys[valid]

        r = np.nanmean(x0[:, None] * y0, axis=0)
        return r

    def beta_at_lag(self, Q: np.ndarray, O: np.ndarray, L: int) -> np.ndarray:
        """
        Standardized regression slope of O(z,t+L) on Q(t):
        beta = cov(zscore(Q), zscore(O)) = corr when both standardized.
        Here we return slope in original units of O relative to std(Q):
            beta = cov(Q,O)/var(Q) = r * std(O)/std(Q)
        """
        if L < 0:
            raise ValueError("Only non-negative lags supported.")
        nt = min(Q.shape[0], O.shape[0])
        if L >= nt - 3:
            return np.full(O.shape[1], np.nan)

        x = Q[: nt - L].astype(float)
        y = O[L: nt, :].astype(float)

        xm = np.nanmean(x); xs = np.nanstd(x)
        ym = np.nanmean(y, axis=0); ys = np.nanstd(y, axis=0)
        if not np.isfinite(xs) or xs == 0.0:
            return np.full(O.shape[1], np.nan)

        # covariance across time
        cov = np.nanmean((x - xm)[:, None] * (y - ym), axis=0)
        beta = cov / (xs ** 2)                  # units: O per unit(Q) ← keep & fix docstring
        beta_std = (cov / (xs ** 2)) * xs       # == r * ys / xs  (units: O per std(Q))
        return beta

    def _corr_by_depth(self,
                       ohc_da: xr.DataArray,
                       qnet: xr.DataArray,
                       t_ohc: str,
                       t_q: str,
                       depth_name: str,
                       lags: Sequence[int],
                       also_beta: bool = True) -> tuple[xr.DataArray, xr.DataArray | None]:
        Q = qnet.astype(float).values                      # [time]
        O = ohc_da.astype(float).values                    # [time, depth]
        # Trim to common length (already done above, but keep guard)
        nt = min(Q.shape[0], O.shape[0])
        Q = Q[:nt]
        O = O[:nt, :]

        lag_vals = list(lags)
        R = np.stack([self.corr_at_lag(Q, O, L) for L in lag_vals], axis=0)
        daR = xr.DataArray(
            R, dims=("lag", depth_name),
            coords=dict(lag=("lag", np.array(lag_vals, dtype=int)),
                        **{depth_name: ohc_da[depth_name]})
        )
        if not also_beta:
            return daR, None

        B = np.stack([self.beta_at_lag(Q, O, L) for L in lag_vals], axis=0)
        daB = xr.DataArray(
            B, dims=("lag", depth_name),
            coords=dict(lag=("lag", np.array(lag_vals, dtype=int)),
                        **{depth_name: ohc_da[depth_name]})
        )
        return daR, daB

    # ---- depth coordinate helper ----
    def _depth_mid_from(self, ds: xr.Dataset, depth_name: str) -> xr.DataArray:
        if "refBottomDepth" in ds:
            edges = ds["refBottomDepth"].values
            top_edges = np.concatenate(([0.0], edges[:-1]))
            mids = 0.5 * (top_edges + edges)
            return xr.DataArray(mids, dims=(depth_name,),
                                name="depth_mid",
                                attrs=dict(units="m", long_name="mid-layer depth"))
        coord = ds.coords.get(depth_name, None)
        if coord is not None:
            return coord
        n = ds.dims.get(depth_name, 0)
        return xr.DataArray(np.arange(1, n + 1), dims=(depth_name,), name="depth_mid")
    
    def _as_lag_depth_array(self, item, x_name: str, x_len: int, depth_name: str | None = None) -> np.ndarray:
        """
        Return a float array shaped (x_len, nz) with dims (lag, depth).
        If item is a DataArray, we reorder its dims to (x_name, depth_name_inferred).
        """
        depth_aliases = {"depth", "z", "nvertlevels", "layer"}

        if isinstance(item, xr.DataArray):
            if item.ndim != 2:
                raise ValueError(f"Expect 2-D DataArray, got dims={item.dims}")

            dims = list(item.dims)

            # Infer depth dim if not provided
            if depth_name is None:
                # prefer a non-x_name dim that looks depth-like
                cand = [d for d in dims if d != x_name and d.lower() in depth_aliases]
                depth_name = cand[0] if cand else [d for d in dims if d != x_name][0]

            # Reorder to (x_name, depth_name)
            if tuple(dims) != (x_name, depth_name):
                if x_name not in dims or depth_name not in dims:
                    raise ValueError(f"Cannot find dims ({x_name}, {depth_name}) in {dims}")
                item = item.transpose(x_name, depth_name)

            Z = np.asarray(item.values, float)

        else:
            # plain ndarray path
            Z = np.asarray(item, float)
            if Z.ndim != 2:
                raise ValueError(f"Expect 2-D array, got shape={Z.shape}")

        # Final shape checks
        if Z.shape[0] != x_len:
            raise ValueError(f"Leading dimension must match x_values length ({x_len}), got {Z.shape[0]}")
        return Z
    
    def alpha_tag(self, idx):
        # 0 -> a/A, 25 -> z/Z, 26 -> aa/AA, etc.
        s, n = "", idx + (self.panel_tag_offset or 0)
        while True:
            s = chr(97 + (n % 26)) + s
            n = n // 26 - 1
            if n < 0: break
        if self.panel_tags == "upper": s = s.upper()
        return f"({s})"
    
    # ---- plotting ----
    def _draw_row(
        self,
        ax_list,
        *,
        x_values,                 # 1-D array of x locations (e.g., lag years)
        depth_ref,                # 1-D depth coordinate (m or already km)
        data_list,                # list of 2-D arrays or DataArrays with dims (x, z) or (lag, depth)
        labels,                   # panel titles (left-aligned)
        vmin,
        vmax,
        row_title=None,           # right-aligned row title
        cmap="YlOrBr",
        fontsize=18,
        fig_idx=0, 
        ylabel="Depth",           # base label text
        xlabel="Lag (years)",     # x-axis label
        y_units_prefer_km=True,   # convert meters→km if depth in meters
        xtick_interval=None,      # spacing in x UNITS (e.g., every 5 years)
        figsize_per_panel=(6,4) 
    ):
        """
        Draw one Hovmöller row across panels using pcolormesh with TRUE axis edges.

        x_values:
            numeric axis (e.g., lag years). Can be non-uniform. Length = nx
        depth_ref:
            1-D depth coordinate aligned with columns of Z (after transpose).
            If in meters and y_units_prefer_km=True, the plot uses km.
        data_list:
            list of DataArrays/ndarrays; each item has shape (nx, nz) when using
            dims (x, depth). If your DataArray is (lag, depth), pass da.transpose("lag", depth_name).
        """
        # --- Build X edges from x_values (supports irregular spacing)
        x = np.asarray(x_values, float)
        if x.ndim != 1 or x.size < 1:
            raise ValueError("x_values must be 1-D with length >= 1")
        if x.size == 1:
            dx = 1.0
            x_edges = np.array([x[0] - 0.5 * dx, x[0] + 0.5 * dx])
        else:
            x_edges = np.empty(x.size + 1, dtype=float)
            x_edges[1:-1] = 0.5 * (x[1:] + x[:-1])
            x_edges[0]    = x[0]  - (x_edges[1] - x[0])
            x_edges[-1]   = x[-1] + (x[-1] - x_edges[-2])

        # --- Build Y edges from depth_ref, convert to km if desired
        z = np.asarray(depth_ref.values if hasattr(depth_ref, "values") else depth_ref, float)
        attr_units = (getattr(depth_ref, "attrs", {}) or {}).get("units", "")
        use_km = y_units_prefer_km and str(attr_units).lower().startswith("m")
        z_plot = z
        y_unit = "km" if (use_km or (attr_units == "")) else attr_units

        if z_plot.size == 1:
            dz = 1.0
            y_edges = np.array([z_plot[0] - 0.5 * dz, z_plot[0] + 0.5 * dz])
        else:
            y_edges = np.empty(z_plot.size + 1, dtype=float)
            y_edges[1:-1] = 0.5 * (z_plot[1:] + z_plot[:-1])
            y_edges[0]    = z_plot[0]  - (y_edges[1] - z_plot[0])
            y_edges[-1]   = z_plot[-1] + (z_plot[-1] - y_edges[-2])

        mappable = None
        for pi, (ax, label, da) in enumerate(zip(ax_list, labels, data_list)):
            # Ensure (nx, nz)
            Z = self._as_lag_depth_array(
                da, x_name="lag", 
                x_len=x.size, 
                depth_name=None
            )
            print(Z.shape)
            #print(xxxx)
            
            # pcolormesh expects (len(y_edges)-1, len(x_edges)-1) grid for Z.T
            im = ax.pcolormesh(
                x_edges, 
                y_edges, 
                Z.T, 
                cmap=cmap, 
                shading="auto",               
                vmin=vmin, vmax=vmax
            )
            mappable = im

            # Cosmetics
            tag = self.alpha_tag(pi+fig_idx) if self.panel_tags else "" 
            panel_label = f"{tag} {label}".strip()
            ax.set_title(panel_label, fontsize=fontsize, loc="left") 
            ax.set_title(row_title, fontsize=fontsize, loc="right")             
            ax.tick_params(labelsize=fontsize * 0.95, length=5, width=1)

            # X ticks
            if xtick_interval is not None and xtick_interval > 0:
                xticks = np.arange(np.nanmin(x), np.nanmax(x) + 1e-9, xtick_interval)
                ax.set_xticks(xticks)
            ax.set_xlabel(xlabel, fontsize=fontsize * 0.95)

        # Left-most y label; depth positive-down → show increasing downward
        for ax in ax_list[1:]:
            ax.tick_params(labelleft=False)   # instead of set_yticklabels([])
            ax.set_ylabel("")

        # Ensure left axis shows ticks and label
        left = ax_list[0]
        left.set_ylabel(f"{ylabel} ({y_unit})", fontsize=fontsize * 0.95)
        left.tick_params(left=True, labelleft=True)  # force-visible on the leftmost

        # Depth increases downward → invert y (z_plot is positive down)
        for ax in ax_list:
            ax.invert_yaxis()

        return mappable

    def _plot_multi_exp(
        self,
        ds,
        exps=None,
        depth_name="depth_mid",
        vmax=1.0,
        fig_path="lagcorr_all.pdf",
        vmin=0.0,
        cmap="YlOrBr",             # or "RdBu_r" if you prefer
        fontsize=18,
        fig_idx=0, 
        xtick_interval=5,           # label every 5 years of lag
        xlabel=None,
        ylabel=None,
        row_title=None,
        figsize_per_panel=(6, 4)
    ):
        # fixed order
        default_order = ["FC", "Full-CPL", "FC-FOSI-S1", "FC-FOSI-S2"]
        if exps is None:
            exps = [e for e in default_order if e in ds.exp.values]
            if not exps:
                exps = list(ds.exp.values)
        else:
            exps = [e for e in default_order if e in exps]

        ncol = len(exps)
        if ncol == 0:
            raise ValueError("No experiments to plot.")

        # Collect data arrays, depth coord, and x (lag) axis
        da0 = ds["corr"].sel(exp=exps[0])
        # Find the actual depth coordinate name
        depth_dim = depth_name if depth_name in da0.coords else [d for d in da0.dims if d != "lag"][0]
        depth_ref = da0[depth_dim]

        x_vals = np.asarray(ds["lag"].values, float)   # arbitrary lags OK (0..50 or sparse list)
        data_list = []
        for exp in exps:
            da = ds["corr"].sel(exp=exp).transpose("lag", depth_dim)  # (lag, depth)
            data_list.append(da)

        # Figure & axes using GridSpec; last (ncol+1) column is the colorbar
        W, H = figsize_per_panel
        fig = plt.figure(figsize=(W * ncol + 0.8 * W, H), constrained_layout=True)
        gs = GridSpec(
            nrows=1, ncols=ncol + 1,  # <-- +1 for colorbar column
            width_ratios=[1] * ncol + [0.05],  # last narrow column for cbar
            wspace=0.10, hspace=0.00, figure=fig
        )

        axes = [fig.add_subplot(gs[0, i]) for i in range(ncol)]
        cax   = fig.add_subplot(gs[0, -1])   # colorbar axis

        mappable = self._draw_row(
            axes,
            x_values=x_vals,
            depth_ref=depth_ref,
            data_list=data_list,
            labels=exps,
            vmin=vmin, vmax=vmax,
            cmap=cmap,
            fontsize=fontsize,
            fig_idx=fig_idx,
            ylabel=(ylabel or "Depth"),
            xlabel=(xlabel or "Lag (years)"),
            row_title=row_title,
            y_units_prefer_km=True,
            xtick_interval=xtick_interval,
            figsize_per_panel=figsize_per_panel,
        )


        # one vertical colorbar
        cb = fig.colorbar(mappable, cax=cax, orientation="vertical")
        cb.ax.set_ylabel("Corr Coef. (unitless)", fontsize=fontsize * 0.9)
        cb.ax.tick_params(labelsize=fontsize * 0.95)
        #pos = cax.get_position()
        #cax.set_position([
        #    pos.x0 +0.1,  # move left (try 0.01–0.03)
        #    pos.y0 +0.035,
        #    pos.width,
        #    pos.height
        #])
        
        fig.savefig(fig_path, dpi=300)
        plt.show()
        plt.close(fig)
        print(f"Saved combined plot to {fig_path}")
