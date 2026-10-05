"""Global ocean time-series extraction from MPAS-Analysis mpasTimeSeriesOcean output (used by alt_global_timeseries)."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple, Sequence, Union, Any, List
import os, re, glob

import numpy as np
import xarray as xr
import cftime

YearToken = Union[str, Tuple[int, int], int]


# ======================================================================
# Low-level request / options dataclasses
# ======================================================================

@dataclass
class SeriesRequest:
    """
    What to extract for one variable/series.

    var_name: variable in file, OR a derived key like "SSH" (implemented below).
    file_key: mpasTimeSeriesOcean, seaIceAreaVolNH, etc.
    sub_dir: experiment-specific subdir (relative to root_dir/exp/)
    region_dim/index/name/reduce: optional selection for regional vars.
    scale: unit conversion factor applied to returned 1D numpy array.
    """
    var_name: str
    file_key: str
    sub_dir: str

    scale: float = 1.0

    # region selection (optional)
    region_dim: Optional[str] = None
    region_index: Optional[int] = None
    region_name: Optional[str] = None
    region_reduce: Optional[str] = None  # None | "mean" | "median"
    region_names_var: Optional[str] = "regionNames"


@dataclass
class ExtractOptions:
    """
    General extraction behavior.
    """
    use_mfdataset: bool = True
    chunks: Optional[dict] = None
    engine: Optional[str] = None

    # time handling
    time_origin: Optional[str] = None  # e.g. "0001-01-01 00:00:00"
    calendar: str = "noleap"


# ======================================================================
# Core extractor
# ======================================================================

class OCNTimeSeriesExtractor:
    """
    Minimal extractor:
      - find files
      - open dataset (single or multi)
      - decode time if possible
      - slice by year_token
      - select region if needed
      - derive special vars (SSH proxy anomaly)
      - return 1D numpy series (NaN allowed)
    """

    def __init__(self, root_dir: str, options: Optional[ExtractOptions] = None):
        self.root_dir = root_dir
        self.opt = options or ExtractOptions()

    # ------------------------- file discovery -------------------------
    def _find_mpas_timeseries_files(self, exp: str, sub_dir: str, file_key: str):
        base = os.path.join(self.root_dir, exp, sub_dir)
        single = glob.glob(os.path.join(base, f"{file_key}.nc"))
        if single:
            return single
        segs = glob.glob(os.path.join(base, f"{file_key}_*.nc"))
        if not segs:
            return []
        def start_stamp(path):
            m = re.search(
                rf"{re.escape(file_key)}_(\d{{6}})-(\d{{6}})\.nc$",
                os.path.basename(path),
            )
            return int(m.group(1)) if m else 0
        segs.sort(key=start_stamp)
        return segs

    def _find_moc_timeseries_files(self, exp: str, sub_dir: str, file_key: str, year_token: YearToken):
        if isinstance(year_token, int):
            yr0, yr1 = year_token, year_token
        elif isinstance(year_token, str):
            yr0, yr1 = map(int, year_token.split("-"))
        else:
            yr0, yr1 = year_token

        base = os.path.join(self.root_dir, exp, sub_dir,"moc")
        selected = []
        for yy in range(yr0, yr1 + 1):
            pattern = os.path.join(base, f"{file_key}_{yy:04d}-{yy:04d}.nc")
            for fname in glob.glob(pattern):
                m = re.search(r"_(\d{4})-\d{4}\.nc$", fname)
                if m:
                    selected.append((int(m.group(1)), fname))
        selected.sort(key=lambda x: x[0])
        return [f for _, f in selected]

    # ------------------------- open dataset -------------------------
    def _open_dataset(self, files: Sequence[str]) -> xr.Dataset:
        time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
        xr_opts = dict(
            engine=self.opt.engine,
            chunks=self.opt.chunks,
            decode_times=time_coder,
        )

        if len(files) == 1 or not self.opt.use_mfdataset:
            return xr.open_dataset(files[0], **xr_opts)

        try:
            return xr.open_mfdataset(
                list(files),
                combine="by_coords",
                parallel=False,
                **xr_opts,
            )
        except Exception:
            # safe concat fallback, detached from file handles
            dsets = []
            try:
                dsets = [xr.open_dataset(fp, **xr_opts) for fp in files]
                probe = next(iter(dsets[0].data_vars.values()))
                tdim = self._guess_time_dim_da(probe)
                out = xr.concat(dsets, dim=tdim)
                out.load()
                return out
            finally:
                for ds in dsets:
                    try:
                        ds.close()
                    except Exception:
                        pass

    # ------------------------- time decoding -------------------------
    @staticmethod
    def _guess_time_dim_da(da: xr.DataArray) -> str:
        for d in da.dims:
            if d.lower().startswith("time"):
                return d
        return da.dims[0]

    @staticmethod
    def _decode_time_from_Time_var(ds: xr.Dataset):
        if "Time" not in ds:
            return None
        da = ds["Time"]
        units = da.attrs.get("units", None)
        cal   = da.attrs.get("calendar", "noleap")
        if units is None:
            return None

        coder = xr.coders.CFDatetimeCoder(use_cftime=True)
        times = xr.coding.times.decode_cf_datetime(
            da.values, units, calendar=cal, time_coder=coder
        )
        return np.asarray(times)

    @staticmethod
    def _decode_time_from_xtime_chars(ds: xr.Dataset, name: str):
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
    def _parse_time_origin(s: Optional[str], calendar: str = "noleap"):
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
        # keep it simple: noleap only
        return cftime.DatetimeNoLeap(y, m, d, hh, mm, ss)

    def _construct_monthly_time(self, n: int):
        t0 = self._parse_time_origin(self.opt.time_origin, self.opt.calendar)
        if t0 is None or n <= 0:
            return None
        return xr.date_range(
            start=t0,
            periods=int(n),
            freq="MS",
            use_cftime=True,
            calendar=self.opt.calendar,
        )

    @staticmethod
    def _slice_indices_from_year_token(time_index, year_token: str):
        a, b = [int(t.strip()) for t in year_token.split("-")]
        years = np.array([t.year for t in time_index], dtype=int)
        mask = (years >= a) & (years <= b)
        if not mask.any():
            return None
        idx = np.nonzero(mask)[0]
        return int(idx[0]), int(idx[-1] + 1)

    @staticmethod
    def _find_time_slice_numeric(year_token: str):
        a, b = [int(t.strip()) for t in year_token.split("-")]
        if a > b:
            a, b = b, a
        m_start = (a - 1) * 12
        m_end   = b * 12
        return m_start, m_end

    # ------------------------- region handling -------------------------
    @staticmethod
    def _infer_region_dim(da: xr.DataArray, tdim: str, prefer: Optional[str] = None):
        if da.ndim <= 1:
            return None
        if prefer and prefer in da.dims:
            return prefer
        for d in da.dims:
            if d != tdim and ("region" in d.lower() or "nregion" in d.lower()):
                return d
        for d in da.dims:
            if d != tdim:
                return d
        return None

    def _decode_region_names(self, ds: xr.Dataset, region_dim: str, region_names_var: Optional[str]):
        if region_names_var and region_names_var in ds:
            arr = ds[region_names_var].values
            if arr.ndim == 2:
                names = ["".join(row.astype(str)).strip().strip("\x00").replace("\x00", "") for row in arr]
                return [n if n else f"{region_dim}[{i}]" for i, n in enumerate(names)]
            if arr.ndim == 1:
                out = []
                for i, v in enumerate(arr):
                    out.append(str(v).strip() or f"{region_dim}[{i}]")
                return out
        return None

    def _region_index_from_name(
        self,
        ds: xr.Dataset,
        region_dim: str,
        target_name: str,
        region_names_var: Optional[str],
    ):
        names = self._decode_region_names(ds, region_dim, region_names_var)
        if not names:
            raise KeyError(f"Region names not found for dim '{region_dim}'.")
        tgt = str(target_name).strip().lower()
        for i, n in enumerate(names):
            if n.lower() == tgt:
                return i
        def norm(s): return "".join(ch for ch in s.lower() if ch.isalnum())
        tgt2 = norm(tgt)
        for i, n in enumerate(names):
            if norm(n) == tgt2:
                return i
        raise KeyError(f"Region name '{target_name}' not found. Available: {names}")

    # ------------------------- derived vars -------------------------
    def _derive_ssh_anom(self, ds: xr.Dataset, out_name: str = "SSH"):
        """
        SSH proxy anomaly = (V/A) - baseline (first 12 steps).
        """
        v_area = "timeMonthly_avg_areaCellGlobal"
        v_vol  = "timeMonthly_avg_volumeCellGlobal"
        if v_area not in ds or v_vol not in ds:
            raise KeyError(f"Need '{v_area}' and '{v_vol}' to derive SSH proxy anomaly.")

        # Identify time dim
        candidates = ["time", "Time", "nMonths", "nTime", "TimeMonthly", "timeMonthly", "nt"]
        tdim = next((d for d in candidates if d in ds.dims), None)
        if tdim is None:
            tdim = self._guess_time_dim_da(ds[v_vol])

        ssh = ds[v_vol] / ds[v_area]

        n0 = min(12, ds.sizes.get(tdim, 0))
        baseline = ssh.isel({tdim: slice(0, n0)}).mean(tdim, keep_attrs=True) if n0 > 0 else 0.0

        ds[out_name] = (ssh - baseline).astype("float64")
        ds[out_name].attrs["units"] = "m"
        ds[out_name].attrs["long_name"] = "SSH Proxy Anomaly (V/A, baseline=first 12 steps)"
        return ds

    # ------------------------- single series API -------------------------
    def extract_series(
        self,
        exp: str,
        req: SeriesRequest,
        year_token: str = "0001-0350",
        *,
        allow_fallback_search: bool = False,
    ) -> np.ndarray:
        """
        Returns 1D numpy array (NaN allowed), already sliced and region-selected.
        """
        # files
        if req.file_key.lower().startswith("moc"):
            files = self._find_moc_timeseries_files(exp, req.sub_dir, req.file_key, year_token)
        else:
            files = self._find_mpas_timeseries_files(exp, req.sub_dir, req.file_key)
            
        if (not files) and allow_fallback_search:
            search_root = os.path.join(self.root_dir, exp)
            candidates = glob.glob(
                os.path.join(
                    search_root,
                    "post",
                    "analysis",
                    "mpas_analysis",
                    "ts_*_climo_*",
                    "**",
                    f"{req.file_key}*.nc",
                ),
                recursive=True,
            )
            if candidates:
                hit = candidates[0]
                subdir = os.path.relpath(os.path.dirname(hit), start=os.path.join(self.root_dir, exp))
                files = self._find_mpas_timeseries_files(exp, subdir, req.file_key)

        if not files:
            raise FileNotFoundError(f"No files found for exp={exp}, subdir={req.sub_dir}, key={req.file_key}")

        ds = self._open_dataset(files)

        # derived vars
        if req.var_name == "SSH":
            ds = self._derive_ssh_anom(ds, out_name="SSH")

        if req.var_name not in ds:
            raise KeyError(f"Variable '{req.var_name}' not found in dataset for exp '{exp}'")

        da = ds[req.var_name]

        # time dim + time index
        tdim = self._guess_time_dim_da(da)
        n = da.sizes[tdim]

        time_idx = self._decode_time_from_Time_var(ds)
        if time_idx is None:
            for nm in ("xtime_startMonthly", "startTime", "xtime_endMonthly", "endTime"):
                ti = self._decode_time_from_xtime_chars(ds, nm)
                if ti is not None and len(ti) == n:
                    time_idx = ti
                    break
        if time_idx is None:
            time_idx = self._construct_monthly_time(n)

        # slice
        use_start = use_end = None
        if (time_idx is not None) and (len(time_idx) == n):
            se = self._slice_indices_from_year_token(time_idx, year_token)
            if se is not None:
                use_start, use_end = se

        if use_start is None:
            m0, m1 = self._find_time_slice_numeric(year_token)
            use_start = max(0, min(n, m0))
            use_end   = max(0, min(n, m1))
            if use_end <= use_start:
                use_start, use_end = 0, n

        da = da.isel({tdim: slice(use_start, use_end)}).astype("float64")

        # region select/reduce if needed
        region_dim = self._infer_region_dim(da, tdim, prefer=req.region_dim)
        if region_dim and region_dim in da.dims:
            rr = (req.region_reduce or None)
            if rr in ("mean", "median"):
                da = da.mean(region_dim, skipna=True) if rr == "mean" else da.median(region_dim, skipna=True)
            else:
                if req.region_index is not None:
                    da = da.isel({region_dim: int(req.region_index)})
                elif req.region_name is not None:
                    idx = self._region_index_from_name(ds, region_dim, req.region_name, req.region_names_var)
                    da = da.isel({region_dim: idx})
                else:
                    da = da.isel({region_dim: 0})

        # ensure 1D
        while da.ndim > 1:
            for d in list(da.dims):
                if d != tdim and da.sizes[d] == 1:
                    da = da.squeeze(d)
                    break
            else:
                for d in list(da.dims):
                    if d != tdim:
                        da = da.mean(d, skipna=True)
                        break

        series = np.asarray(da.values, dtype=float).reshape(-1)
        if not np.isfinite(series).any():
            raise ValueError(f"All values NaN for exp='{exp}', var='{req.var_name}'")
        return series * float(req.scale)

    # ------------------------- group API -------------------------
    def extract_group(
        self,
        exp: str,
        requests: Dict[str, SeriesRequest],
        year_token: str,
    ) -> Dict[str, np.ndarray]:
        """
        Extract multiple variables for one experiment.
        Returns dict: {alias: series_1d}
        """
        out: Dict[str, np.ndarray] = {}
        for alias, req in requests.items():
            out[alias] = self.extract_series(exp, req, year_token=year_token)
        return out

    def extract_many(
        self,
        exps: Sequence[str],
        requests: Dict[str, SeriesRequest],
        year_token: str,
    ) -> Dict[str, Dict[str, np.ndarray]]:
        """
        Extract multiple variables for multiple experiments.
        Returns: {exp: {alias: series}}
        """
        out: Dict[str, Dict[str, np.ndarray]] = {}
        for exp in exps:
            out[exp] = self.extract_group(exp, requests, year_token=year_token)
        return out


# ======================================================================
# High-level convenience wrappers (merge of ocn_extract_simple.py)
# ======================================================================

@dataclass
class ExtractConfig:
    """
    Configuration for batch extraction driven by var_dict + exp_subdirs.
    """
    year_token: str

    # region selection (used only for region-style MPAS vars)
    region_dim: Optional[str] = None
    region_index: Optional[int] = None
    region_name: Optional[str] = None
    region_reduce: Optional[str] = None
    region_names_var: Optional[str] = None


def build_requests_from_var_dict(
    var_dict: Dict[str, Dict[str, Any]],
    exp_subdir: str,
    cfg: ExtractConfig,
) -> Dict[str, SeriesRequest]:
    """
    Convert your var_dict entry schema into SeriesRequest objects.

    Heuristic:
      - Apply region selection only if the variable name suggests an ocean region series.
      - Ice diagnostics typically skip region selection.
      - SSH is handled as a derived var by extractor (var_name="SSH").
    """
    reqs: Dict[str, SeriesRequest] = {}

    for key, v in var_dict.items():
        var_name = v["name"]
        file_key = v["key"]
        scale    = float(v.get("fact", 1.0))

        # region selection only for region-like MPAS time series variables
        use_region = (
            cfg.region_dim is not None
            and isinstance(var_name, str)
            and ("WithinOceanRegion" in var_name or "OceanRegion" in var_name or "Region" in var_name)
        )

        reqs[key] = SeriesRequest(
            var_name=var_name,
            file_key=file_key,
            sub_dir=exp_subdir,
            scale=scale,
            region_dim=cfg.region_dim if use_region else None,
            region_index=cfg.region_index if use_region else None,
            region_name=cfg.region_name if use_region else None,
            region_reduce=cfg.region_reduce if use_region else None,
            region_names_var=cfg.region_names_var if use_region else None,
        )

    return reqs


def extract_all(
    root_dir: str,
    exps: List[str],
    exp_subdirs: Dict[str, str],
    var_dict: Dict[str, Dict[str, Any]],
    cfg: ExtractConfig,
    *,
    time_origin: Optional[str] = None,
    calendar: str = "noleap",
    engine: Optional[str] = None,
    chunks: Optional[dict] = None,
    use_mfdataset: bool = True,
) -> Dict[str, Dict[str, np.ndarray]]:
    """
    Batch extract multiple vars for multiple experiments.

    Returns:
      data[exp][var_key] = 1D numpy array (scaled), NaNs allowed.
    """
    opt = ExtractOptions(
        use_mfdataset=use_mfdataset,
        chunks=chunks,
        engine=engine,
        time_origin=time_origin,
        calendar=calendar,
    )
    X = OCNTimeSeriesExtractor(root_dir=root_dir, options=opt)

    out: Dict[str, Dict[str, np.ndarray]] = {}
    for exp in exps:
        subdir = exp_subdirs[exp]
        reqs = build_requests_from_var_dict(var_dict, subdir, cfg)
        out[exp] = X.extract_group(exp, reqs, year_token=cfg.year_token)

    return out

# ======================================================================
# Registry / configuration container (regions + variables)
# ======================================================================

class OCNRegistry:
    """
    Container for region and variable definitions used by time-series extraction.

    This class is intentionally lightweight:
      - stores canonical region_dict and var_dict
      - provides helpers to build ExtractConfig
    """

    def __init__(self):
        self.region_dict = {
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

        self.var_dict = {
            "SST": {
                "name":  "timeMonthly_avg_avgValueWithinOceanRegion_avgSurfaceTemperature",
                "short": "SST",
                "unit":  r"$^o$C",
                "fact":  1.0,
                "min":   17.2,
                "max":   18.6,
                "label": "Sea Surface Temperature (SST)",
                "key":   "mpasTimeSeriesOcean",
            },
            "SSS": {
                "name":  "timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerSalinity",
                "short": "SSS",
                "unit":  r"PSU",
                "fact":  1.0,
                "min":   34.60,
                "max":   34.85,
                "label": "Sea Surface Salinity (SSS)",
                "key":   "mpasTimeSeriesOcean",
            },
            "SSH": {
                "name":  "SSH",
                "short": "SSH",
                "unit":  r"cm",
                "fact":  100.0,
                "min":   -8,
                "max":   15,
                "label": r"Sea Surface Height Anomaly ($\Delta$SSH)",
                "key":   "mpasTimeSeriesOcean",
            },
            "iceVolumeNH": {
                "name":  "iceVolume",
                "short": "Ice Volume",
                "unit":  r"x10$^3$ km$^3$",
                "fact":  1e-12,
                "min":   10,
                "max":   35,
                "label": "Sea Ice Volume (NH)",
                "key":   "seaIceAreaVolNH",
            },
            "iceVolumeSH": {
                "name":  "iceVolume",
                "short": "Ice Volume",
                "unit":  r"x10$^3$ km$^3$",
                "fact":  1e-12,
                "min":   5,
                "max":   30,
                "label": "Sea Ice Volume (SH)",
                "key":   "seaIceAreaVolSH",
            },
            "iceAreaNH": {
                "name":  "iceArea",
                "short": "Ice Area",
                "unit":  r"x10$^6$km$^2$",
                "fact":  1e-12,
                "min":   9,
                "max":   13,
                "label": "Sea Ice Area (NH)",
                "key":   "seaIceAreaVolNH",
            },
            "iceAreaSH": {
                "name":  "iceArea",
                "short": "Ice Area",
                "unit":  r"x10$^6$km$^2$",
                "fact":  1e-12,
                "min":   6,
                "max":   16,
                "label": "Sea Ice Area (SH)",
                "key":   "seaIceAreaVolSH",
            },
            "AMOC26": {
                "name":  "mocAtlantic26",   # <-- actual var in file
                "short": "AMOC (26.5°N)",
                "unit":  "Sv",
                "fact":  1.0,
                "min":   10,
                "max":   25,
                "label": "AMOC (26.5°N)",
                "key":   "mocTimeSeries",   # file stem: mocTimeSeries_YYYY-YYYY.nc
            },
            "AMOC": {
                "name":  "mocAtlantic",   # <-- actual var in file
                "short": "AMOC",
                "unit":  "Sv",
                "fact":  1.0,
                "min":   10,
                "max":   25,
                "label": "AMOC",
                "key":   "mocTimeSeries",   # file stem: mocTimeSeries_YYYY-YYYY.nc
            },
        }

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------
    def get_region(self, region_key: str):
        if region_key not in self.region_dict:
            raise KeyError(f"Unknown region '{region_key}'. Available: {list(self.region_dict)}")
        return self.region_dict[region_key]

    def get_var_dict(self):
        return self.var_dict

    def build_extract_config(
        self,
        *,
        year_token: str,
        region_key: str,
        region_dim: Optional[str] = "nOceanRegions",
        region_names_var: Optional[str] = None,
    ) -> ExtractConfig:
        r = self.get_region(region_key)
        return ExtractConfig(
            year_token=year_token,
            region_dim=region_dim,
            region_index=r["index"],
            region_name=r.get("name"),
            region_reduce=None,
            region_names_var=region_names_var,
        )