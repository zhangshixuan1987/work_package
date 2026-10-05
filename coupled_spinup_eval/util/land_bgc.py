"""Land BGC spin-up stability analysis (ELM time series).

Shared base class SpinupStabilityBGC + one subclass per analysis:
  LandCtrlStability    (10_land_ctrl_analysis)
  LandCompareStability (11_land_compare_analysis)
"""
import os
import re
import string
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Tuple, List, Optional
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator
try:
    import cftime
except ImportError:
    cftime = None
import re, string
from typing import Dict, List, Optional, Sequence, Tuple
try:
    import cftime  # optional for no-leap calendars
except Exception:
    cftime = None


@dataclass
class SpinupStabilityBGC:
    """Shared land BGC time-series reader/processor (ELM post/time_series). Use a subclass."""

    root_dir: str = "/lcrc/group/e3sm/ac.szhang/acme_scratch/e3sm_project/E3SMv3_testings"

    sub_dir: str = "time_series/lnd"

    experiments: Tuple[str, ...] = ("20231209.v3.LR.piControl-spinup.chrysalis",)

    case_ids: Tuple[str, ...] = ("v3.LR.piControl-spinup",)

    exp_labels: Tuple[str, ...] = ("FC",)

    time_strs: Tuple[str, ...] = ("0001-0300",)  # per-experiment span (inclusive)

    exp_subdirs: Optional[Dict[str, str]] = None

    region_index: int = 0

    land_file: str = (
        "/lcrc/group/acme/ac.szhang/acme_scratch/e3sm_project/E3SMv3_testings/"
        "paper_material/data/landmask.nc"
    )

    annual_hist: bool = False

    calendar: str = "noleap"

    time_anchor: str = "start"            # "start", "middle", "end"

    force_time_from_span: bool = True     # construct from time_strs if time is missing

    month_index_correction: str = "auto"  # "auto" | "as_is" | "shift_back_1"

    replace_monthly_time_from_span: bool = True

    verbose: bool = True

    var_setup: Dict[str, Dict[str, object]] = field(
        default_factory=lambda: {
            "TOTECOSYSC":     {"scale": 1e-15, "units": "Pg C",          "title": "Total Ecosystem Carbon (TOTECOSYSC)"},
            "TOTSOMC":        {"scale": 1e-15, "units": "Pg C",          "title": "Total Soil Organic Matter Carbon (TOTSOMC)"},
            "TOTVEGC":        {"scale": 1e-15, "units": "Pg C",          "title": "Total Vegetation Carbon (TOTVEGC)"},
            "TLAI":           {"scale": 1.0,   "units": "m² m⁻²",        "title": "Total Leaf Area Index (TLAI)"},
            "GPP":            {"scale": 1e-15, "units": "Pg C yr⁻¹",     "title": "Gross Primary Production (GPP)"},
            "TWS":            {"scale": 1e-3,  "units": "m",             "title": "Total Water Storage (TWS)"},
            "H2OSNO":         {"scale": 1.0,   "units": "mm",            "title": "Snow Water Equivalent (H2OSNO)"},
            "HTOP":           {"scale": 1.0,   "units": "m",             "title": "Canopy Top (HTOP)"},
            "TSOI_10CM":      {"scale": 1.0,   "units": "°C",            "title": "Soil Temp Top 10 cm (TSOI_10CM)"},
            "SOILWATER_10CM": {"scale": 1.0,   "units": "kg m⁻²",        "title": "Soil Water Top 10 cm"},
            "FSH":            {"scale": 1.0,   "units": "W m⁻²",         "title": "Sensible Heat (FSH)"},
            "HCSOI":          {"scale": 1e-3,  "units": "MJ m⁻²",        "title": "Soil Heat Content (HCSOI)"},
        }
    )

    seconds_in_year: float = 60.0 * 60.0 * 24.0 * 365.0

    def _parse_range_from_fname(self, fname: str) -> tuple[int, int, int, int]:
        _FN_RANGE = re.compile(r'_(\d{6})_(\d{6})\.nc$')
        m = _FN_RANGE.search(Path(fname).name)
        if not m:
            return (0, 0, 0, 0)
        y0m0, y1m1 = m.group(1), m.group(2)
        y0, m0 = int(y0m0[:4]), int(y0m0[4:])
        y1, m1 = int(y1m1[:4]), int(y1m1[4:])
        return (y0, m0, y1, m1)

    def _var_files(self, exp: str, var: str) -> List[str]:
        if self.exp_subdirs and exp in self.exp_subdirs:
            d = Path(self.root_dir) / exp / self.exp_subdirs[exp]
        else:
            d = Path(self.root_dir) / exp / "post" / self.sub_dir
        files = [str(p) for p in d.glob(f"{var}*.nc")]
        if not files:
            if var == "H2OSNO":
                return []  # allow missing snow
            raise FileNotFoundError(f"No files found for var={var} in {d}")
        files = sorted(files, key=lambda f: self._parse_range_from_fname(f))
        return files

    @staticmethod
    def _days_in_month_noleap(month: int) -> int:
        return [31,28,31,30,31,30,31,31,30,31,30,31][month-1]

    def _anchor_monthly_time(self, tvals) -> xr.DataArray:
        """Snap monthly timestamps to month start/middle/end consistently."""
        if len(tvals) == 0:
            return xr.DataArray(tvals, dims=["time"], name="time")
        is_cftime = hasattr(tvals[0], "__class__") and "cftime" in type(tvals[0]).__module__

        def _snap_dt(dt):
            if is_cftime:
                y, m = dt.year, dt.month
                if self.time_anchor == "start":
                    return cftime.DatetimeNoLeap(y, m, 1)
                elif self.time_anchor == "middle":
                    dim = self._days_in_month_noleap(m); d = (dim + 1)//2
                    return cftime.DatetimeNoLeap(y, m, d)
                else:
                    dim = self._days_in_month_noleap(m)
                    return cftime.DatetimeNoLeap(y, m, dim)
            else:
                ts = pd.to_datetime(dt)
                if self.time_anchor == "start":
                    d = 1
                elif self.time_anchor == "middle":
                    d = (self._days_in_month_noleap(ts.month) + 1)//2
                else:
                    d = self._days_in_month_noleap(ts.month)
                return ts.replace(day=d, hour=0, minute=0, second=0, microsecond=0)

        snapped = [_snap_dt(x) for x in tvals]
        return xr.DataArray(snapped, dims=["time"], name="time")

    def _shift_one_month(self, dt, minus=True):
        """Shift a timestamp by ±1 month and re-anchor to self.time_anchor."""
        sign = -1 if minus else +1

        def _anchor_day(year, month):
            if self.time_anchor == "start":
                return 1
            dim = self._days_in_month_noleap(month)
            if self.time_anchor == "middle":
                return (dim + 1)//2
            return dim

        # cftime (no-leap)
        if hasattr(dt, "__class__") and "cftime" in type(dt).__module__:
            y, m = dt.year, dt.month
            m += sign
            if m == 0:
                m = 12; y -= 1
            elif m == 13:
                m = 1; y += 1
            d = _anchor_day(y, m)
            return cftime.DatetimeNoLeap(y, m, d)

        # numpy/pandas
        ts = pd.to_datetime(dt)
        y, m = ts.year, ts.month
        m += sign
        if m == 0:
            m = 12; y -= 1
        elif m == 13:
            m = 1; y += 1
        d = _anchor_day(y, m)
        return ts.replace(year=y, month=m, day=d, hour=0, minute=0, second=0, microsecond=0)

    def _maybe_fix_month_index(self, tvals):
        """
        Fix “next-month stamp” if requested.
        - "as_is": no change
        - "shift_back_1": always subtract one month
        - "auto": detect (start=Feb, end=Jan, length % 12 == 0) then subtract one month
        """
        if self.month_index_correction == "as_is" or len(tvals) < 24:
            return tvals

        try:
            months = [int(getattr(t, "month", pd.to_datetime(t).month)) for t in tvals]
        except Exception:
            return tvals

        do_shift = False
        if self.month_index_correction == "shift_back_1":
            do_shift = True
        elif self.month_index_correction == "auto":
            if months[0] == 2 and months[-1] == 1 and (len(months) % 12 == 0):
                do_shift = True

        if not do_shift:
            return tvals

        return [self._shift_one_month(t, minus=True) for t in tvals]

    def _build_monthly_time_from_span(self, n: int, exp_index: int) -> xr.DataArray:
        """Create a monthly CFTime axis of length n starting at <time_strs[exp_index]>/01."""
        span = self.time_strs[exp_index] if len(self.time_strs) else "0001-0001"
        y0 = int(span.split("-")[0]); m0 = 1
        if cftime is not None and self.calendar.lower() == "noleap":
            months = []
            y, m = y0, m0
            for _ in range(n):
                months.append(cftime.DatetimeNoLeap(y, m, 1))
                m += 1
                if m == 13:
                    m = 1; y += 1
            return xr.DataArray(months, dims=["time"], name="time")
        # Fallback: numeric year.fraction
        vals = []
        y, m = y0, m0
        for _ in range(n):
            vals.append(y + (m - 1) / 12.0)
            m += 1
            if m == 13:
                m = 1; y += 1
        return xr.DataArray(vals, dims=["time"], name="time")

    def _ensure_or_construct_time(self, ds: xr.Dataset, files: List[str], exp_index: int) -> xr.Dataset:
        """Return ds with a consistent monthly/annual 'time' coordinate."""
        # Use existing time if present: possibly replace (monthly) or anchor+shift
        if "time" in ds or "Time" in ds:
            tname = "time" if "time" in ds else "Time"
            if tname != "time":
                ds = ds.rename({tname: "time"})
            if "time" in ds.coords and ds.sizes.get("time", 0) > 0:
                try:
                    if ds.sizes["time"] >= 24:
                        if self.replace_monthly_time_from_span:
                            new_time = self._build_monthly_time_from_span(ds.sizes["time"], exp_index)
                            ds = ds.assign_coords(time=new_time)
                        else:
                            anchored = self._anchor_monthly_time(list(ds["time"].values))
                            fixed = self._maybe_fix_month_index(list(anchored.values))
                            ds = ds.assign_coords(time=xr.DataArray(fixed, dims=["time"], name="time"))
                except Exception:
                    pass

                if self.verbose:
                    n_t = int(ds.sizes.get("time", 0))
                    try:
                        t0 = ds["time"].values[0] if n_t else "N/A"
                        t1 = ds["time"].values[-1] if n_t else "N/A"
                    except Exception as e:
                        t0, t1 = f"ERR:{e}", f"ERR:{e}"
                    print(f"[DEBUG] _ensure_or_construct_time: existing time | "
                          f"n={n_t} | span={t0} → {t1} | anchor={self.time_anchor} | "
                          f"month_index_correction={self.month_index_correction} | "
                          f"replace_from_span={self.replace_monthly_time_from_span}")
                return ds

        # Otherwise construct monthly time from span or filenames
        tdim = None
        for cand in ("time", "Time", "nTime", "t", "TimeMonthly", "Time_counter"):
            if cand in ds.dims:
                tdim = cand; break
        if tdim is None:
            return ds

        n = ds.sizes[tdim]
        if self.force_time_from_span:
            time = self._build_monthly_time_from_span(n, exp_index)
        else:
            y0, m0, *_ = self._parse_range_from_fname(files[0])
            if y0 == 0:
                time = self._build_monthly_time_from_span(n, exp_index)
            else:
                # build from filename start
                if cftime is not None and self.calendar.lower() == "noleap":
                    months = []
                    y, m = y0, m0
                    for _ in range(n):
                        months.append(cftime.DatetimeNoLeap(y, m, 1))
                        m += 1
                        if m == 13:
                            m = 1; y += 1
                    time = xr.DataArray(months, dims=["time"], name="time")
                else:
                    vals, y, m = [], y0, m0
                    for _ in range(n):
                        vals.append(y + (m - 1) / 12.0)
                        m += 1
                        if m == 13:
                            m = 1; y += 1
                    time = xr.DataArray(vals, dims=["time"], name="time")

        if tdim != "time":
            ds = ds.rename({tdim: "time"})
        ds = ds.assign_coords(time=time)

        if self.verbose and n:
            print(f"[DEBUG] _ensure_or_construct_time: constructed monthly time | "
                  f"start={str(time.values[0])} len={n} | replace_from_span={self.force_time_from_span}")
        return ds

    @staticmethod
    def _annualize(da: xr.DataArray, already_annual: bool) -> xr.DataArray:
        if already_annual:
            return da
        if np.issubdtype(da["time"].dtype, np.number):
            years = np.floor(da["time"].values).astype(int)
            df = pd.DataFrame({"v": da.values, "y": years})
            grouped = df.groupby("y")["v"].mean()
            return xr.DataArray(grouped.values, coords={"time": grouped.index}, dims=["time"])
        return da.groupby("time.year").mean("time").rename({"year": "time"})

    @staticmethod
    def _year_centers(years: np.ndarray) -> np.ndarray:
        return 0.5 * (years[:-1] + years[1:])

    def _load_land_mask(self) -> Tuple[xr.DataArray, xr.DataArray, float]:
        fr = xr.open_dataset(self.land_file)
        if "landfrac" not in fr or "area" not in fr:
            raise KeyError("land_file must contain 'landfrac' and 'area'")
        landfrac = fr["landfrac"].load()
        area = fr["area"].load()
        units = str(area.attrs.get("units", "")).lower()
        if any(k in units for k in ("km2", "km^2", "square kilometer", "square kilometres", "square kilometers")):
            area = area * 1e6
        elif not units:
            if float(area.max()) < 1e7:
                area = area * 1e6
        landarea = landfrac * area
        landareaC = float(landarea.sum())
        fr.close()
        return landfrac, area, landareaC

    def _read_series(self, var: str, exp: str, exp_index: int, combine: str, engine: str) -> xr.DataArray:
        files = self._var_files(exp, var)
        if not files:
            raise KeyError(f"{var} not available")

        ds = xr.open_mfdataset(files, combine=combine, parallel=False, chunks=None, engine=engine)
        ds = self._ensure_or_construct_time(ds, files, exp_index)
        if "time" not in ds:
            raise KeyError(f"Unable to find or construct 'time' for {var} in {exp}")

        da = ds[var]

        # Known region dims
        for cand in ["region", "nRegions", "nOceanRegions", "nreg", "ntrsc"]:
            if cand in da.dims:
                if self.region_index >= da.sizes[cand]:
                    raise IndexError(f"region_index={self.region_index} out of range for dim {cand} size={da.sizes[cand]}")
                da = da.isel({cand: self.region_index})
                break

        # Drop bounds-like dims if present
        for bnd in ["bounds", "bnds", "nbnd", "nbounds"]:
            if bnd in da.dims and da.sizes[bnd] >= 1:
                da = da.isel({bnd: 0})

        # Reduce any remaining non-time dims
        non_time = [d for d in da.dims if d != "time"]
        if non_time:
            if len(non_time) == 1:
                d = non_time[0]
                if da.sizes[d] == 1:
                    da = da.squeeze(d)
                elif da.sizes[d] == 3 and self.region_index < 3:
                    da = da.isel({d: self.region_index})
                else:
                    da = da.mean(d)
            else:
                for d in list(non_time):
                    if da.sizes[d] == 1:
                        da = da.squeeze(d)
                non_time = [d for d in da.dims if d != "time"]
                if non_time:
                    da = da.mean(non_time)

        # Final: ensure time sorted & loaded
        da = da.sortby("time").load()

        if self.verbose:
            n_t = int(da.sizes.get("time", 0))
            try:
                t0 = da["time"].values[0] if n_t else "N/A"
                t1 = da["time"].values[-1] if n_t else "N/A"
            except Exception as e:
                t0, t1 = f"ERR:{e}", f"ERR:{e}"
            print(f"[DEBUG] {var} ({exp}): loaded series (monthly or annual) | "
                  f"n={n_t} | span={t0} → {t1} | anchor={self.time_anchor} | "
                  f"month_index_correction={self.month_index_correction} | "
                  f"replace_from_span={self.replace_monthly_time_from_span}")
        return da

    def _processed_path(self, exp: str, var: str, processed_dir: Optional[str]) -> Path:
        base = Path(processed_dir or ".")
        return base / f"{exp}_{var}_processed.nc"

    def save_processed_series(
        self,
        exp: str,
        exp_index: int,
        out_dir: str,
        overwrite: bool = False,
        combine: str = "by_coords",
        engine: str = "netcdf4",
    ) -> None:
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        if exp is None:
            exp = self.experiments[exp_index]
        _, _, landareaC = self._load_land_mask()

        for v in self.var_setup.keys():
            fname = self._processed_path(exp, v, out_dir)
            if fname.exists() and not overwrite:
                print(f"Skip existing {fname}")
                continue

            years, vals, dels, units = self._prep_series(
                v, exp, exp_index, landareaC, combine, engine
            )

            # Optional: persist trend alongside annual/delta
            tmid = np.array([], float); trend = np.array([], float)
            if self.use_trend:
                tmid, trend = self._rolling_trend(
                    years, vals, window=self.trend_window, min_frac=self.trend_min_frac
                )

            data_vars = {
                f"{v}_ann":   (["time"],     vals, {"units": units}),
                f"{v}_delta": (["time_mid"], dels),
            }
            coords = {"time": years, "time_mid": self._year_centers(years)}
            if tmid.size and trend.size:
                data_vars[f"{v}_trend"] = (["trend_mid"], trend)
                coords["trend_mid"] = tmid

            ds_out = xr.Dataset(
                data_vars, coords=coords,
                attrs={
                    "experiment": exp,
                    "trend_window_years": self.trend_window,
                    "trend_min_fraction": self.trend_min_frac,
                },
            )
            ds_out.to_netcdf(fname)
            print(f"Wrote {fname}")

    def _split_units_for_display(self, var: str, stored_units: str) -> tuple[str, str]:
        base = {
            "TOTECOSYSC": "Pg C",
            "TOTSOMC": "Pg C",
            "TOTVEGC": "Pg C",
            "TLAI": "m² m⁻²",
            "GPP": "Pg C yr⁻¹",
            "TWS": "m",
            "H2OSNO": "mm",
            "HTOP": "m",
            "TSOI_10CM": "°C",
            "SOILWATER_10CM": "kg m⁻²",
            "FSH": "W m⁻²",
            "HCSOI": "MJ m⁻²",
        }
        if var in base:
            series_units = base[var]
        else:
            u = stored_units
            u = u.replace("$^{-1}$", "⁻¹").replace("^−1", "⁻¹").replace("^−¹", "⁻¹")
            u = re.sub(r'\s*[/ ]?yr(?:\s*[\-−]?\s*1|⁻¹)?\)?', '', u, flags=re.IGNORECASE)
            u = u.replace("()", "").strip()
            u = re.sub(r'\s+', ' ', u).strip()
            series_units = u if u else stored_units

        # Default right-axis label = per-year version of the series units
        if re.search(r'(?:/| )yr', series_units, flags=re.IGNORECASE) or "yr⁻¹" in series_units:
            delta_units = series_units
        else:
            delta_units = f"{series_units} yr⁻¹"

        return series_units, delta_units


@dataclass
class LandCtrlStability(SpinupStabilityBGC):
    """Full-CPL spin-up + piControl land BGC trends, e-folding and Liao thresholds."""

    use_trend: bool = True

    trend_window: int = 30

    trend_min_frac: float = 0.8

    bgc_focus: bool = True

    bgc_thresholds: Dict[str, float] = field(
        default_factory=lambda: {
            # robust Theil–Sen slope tolerance (PgC/yr)
            "TOTECOSYSC_slope": 0.01,   # ~10 TgC/yr
            "TOTSOMC_slope":   0.006,   # soil slower; a bit stricter
            # relative deviation (unitless) of last-20% window to its mean
            "tot_stock_rel_err": 0.005, # 0.5%
            # acceptable e-folding time ceiling (years)
            "tau_max_yrs": 300.0,
        }
    )

    use_liao_criteria: bool = True                  # turn on globalized thresholds

    liao_flux_thresh_gpm2yr: float = 1.0            # 1 gC m-2 yr-1 for TEC drift

    liao_passive_thresh_gpm2yr: float = 0.5         # 0.5 gC m-2 yr-1 for passive/slow pool

    liao_disequilibrium_frac: float = 0.03          # 3% (defined for future grid-based use)

    thresholds: Dict[str, float] = field(
        default_factory=lambda: {
            "TOTECOSYSC": 0.02,
            "TOTSOMC": 0.02,
            "TOTVEGC": 0.02,
            "TLAI": 0.02,
            "GPP": 0.02,
            "TWS": 0.001,
            "H2OSNO": 0.2,
            "HTOP": 0.001,
            "TSOI_10CM": 0.02,
            "SOILWATER_10CM": 0.1,
            "FSH": 0.2,
            "HCSOI": 0.01,
        }
    )

    def _apply_liao_thresholds(self, landarea_m2: float) -> None:
        """
        Convert per-area drift thresholds (gC m^-2 yr^-1) into global PgC yr^-1
        and update:
          - self.thresholds[...]           (used by right-axis threshold lines)
          - self.bgc_thresholds[...]       (used by spin-up metrics; PgC/yr)
        Also prints the computed values when self.verbose is True.
        """
        if not self.use_liao_criteria:
            return

        # gC m^-2 yr^-1  ->  PgC yr^-1
        tec_pgCyr  = (self.liao_flux_thresh_gpm2yr    * landarea_m2) / 1e15
        slow_pgCyr = (self.liao_passive_thresh_gpm2yr * landarea_m2) / 1e15

        # Update panel thresholds (Δ stock per year) for stocks
        self.thresholds.update({
            "TOTECOSYSC": float(tec_pgCyr),
            "TOTSOMC":    float(slow_pgCyr),
            "TOTVEGC":    float(tec_pgCyr),  # keep same as TEC; tweak if desired
            "GPP":        float(tec_pgCyr),  # keep same as TEC; tweak if desired
        })

        # Update BGC slope tolerances (PgC/yr) for metrics
        self.bgc_thresholds.update({
            "TOTECOSYSC_slope": float(tec_pgCyr),
            "TOTVEGC_slope": float(tec_pgCyr),
            "GPP_slope": float(tec_pgCyr),
            "TOTSOMC_slope":    float(slow_pgCyr),
        })

        # Debug print
        if self.verbose:
            tec_gCyr  = tec_pgCyr  * 1e15  # for readability vs the right axis
            slow_gCyr = slow_pgCyr * 1e15
            print("[LIAO] Using land area = {:.3e} m^2".format(landarea_m2))
            print("[LIAO] TEC threshold:     {:.6f} PgC/yr  ({:.3e} gC/yr)".format(tec_pgCyr,  tec_gCyr))
            print("[LIAO] TOTSOMC threshold: {:.6f} PgC/yr  ({:.3e} gC/yr)".format(slow_pgCyr, slow_gCyr))
            print("[LIAO] TOTVEGC threshold: {:.6f} PgC/yr  ({:.3e} gC/yr)".format(tec_pgCyr,  tec_gCyr))
            print("[LIAO] GPP threshold:     {:.6f} PgC/yr  ({:.3e} gC/yr)".format(tec_pgCyr,  tec_gCyr))

    def _prep_series(
        self,
        var: str,
        exp: str,
        exp_index: int,
        landareaC: float,
        combine: str,
        engine: str,
        processed_dir: Optional[str] = None,
        subper: int = 1
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, str]:

        try:
            da = self._read_series(var, exp, exp_index, combine, engine)
        except KeyError:
            return np.array([], float), np.array([], float), np.array([], float), self.var_setup[var]["units"]

        scale = float(self.var_setup[var]["scale"])
        units = str(self.var_setup[var]["units"])

        # Annualize (after consistent monthly handling)
        ann = self._annualize(da, self.annual_hist)

        if self.verbose:
            n_before = int(da.sizes.get("time", 0))
            n_after  = int(ann.sizes.get("time", 0))
            try:
                t0 = ann["time"].values[0] if n_after else "N/A"
                t1 = ann["time"].values[-1] if n_after else "N/A"
            except Exception as e:
                t0, t1 = f"ERR:{e}", f"ERR:{e}"
            print(f"[DEBUG] {var} ({exp}): annual_hist={self.annual_hist} | "
                  f"annualized n={n_before} → {n_after} | span={t0} → {t1}")

        # Clip to time_strs window
        span = self.time_strs[exp_index]
        y0, y1 = map(int, span.split("-"))
        if np.issubdtype(ann["time"].dtype, np.number):
            ann = ann.where((ann["time"] >= y0) & (ann["time"] <= y1), drop=True)
        else:
            ann = ann.sel(time=slice(f"{y0:04d}", f"{y1:04d}"))
        ann = ann.dropna("time", how="all") if "time" in ann.dims else ann

        if self.verbose:
            n_clip = int(ann.sizes.get("time", 0))
            t0c = ann["time"].values[0] if n_clip else "N/A"
            t1c = ann["time"].values[-1] if n_clip else "N/A"
            print(f"[DEBUG] {var} ({exp}): clipped to {y0}-{y1} | n={n_clip} | span={t0c} → {t1c}")

        if ann.sizes.get("time", 0) < 2:
            years_num = ann["time"].values.astype(float) if "time" in ann.coords else np.array([], float)
            return years_num, ann.to_numpy().astype(float), np.array([], float), units

        # Convert time to numeric years
        if np.issubdtype(ann["time"].dtype, np.number):
            years_num = np.array(ann["time"].values, dtype=float)
        else:
            def _to_year(t):
                return getattr(t, "year", pd.Timestamp(t.astype("datetime64[ns]")).year)
            years_num = np.array([_to_year(t) for t in ann["time"].values], dtype=float)

        vals = ann.load().to_numpy().astype(float)

        # Conversions like NCL
        if var in ("TOTECOSYSC", "TOTSOMC", "TOTVEGC"):
            vals *= landareaC
        if var == "GPP":
            vals *= landareaC
            vals *= self.seconds_in_year
        if var == "TSOI_10CM":
            # convert from K to degC
            vals -= 273.15
        vals *= scale

        dels = np.diff(vals) if vals.size >= 2 else np.array([], float)

        return years_num, vals, dels, units

    @staticmethod
    def _rolling_trend(
        years: np.ndarray,
        values: np.ndarray,
        window: int,
        min_frac: float = 0.8
    ) -> tuple[np.ndarray, np.ndarray]:
        if years.size < 2 or values.size < 2 or years.size != values.size or window < 2:
            return np.array([], float), np.array([], float)
        y = years.astype(float); v = values.astype(float)
        n = len(y); W = int(window)
        mids, slopes = [], []
        need = max(2, int(np.ceil(min_frac * W)))
        for i in range(0, n - W + 1):
            yy = y[i:i+W]; vv = v[i:i+W]
            msk = np.isfinite(yy) & np.isfinite(vv)
            mids.append(0.5 * (yy[0] + yy[-1]))
            if msk.sum() < need:
                slopes.append(np.nan); continue
            yyv = yy[msk]; vvv = vv[msk]
            x = yyv - yyv.mean(); denom = np.dot(x, x)
            s = np.nan if denom == 0.0 else np.dot(x, vvv - vvv.mean()) / denom
            slopes.append(s)
        return np.array(mids, float), np.array(slopes, float)

    def analyze_experiment(
        self,
        var_dict: dict = None,
        exp: str = None,
        exp_index: int = 0,
        combine: str = "by_coords",
        engine: str = "netcdf4",
        subper: int = 1,
        trd_ymin1: int = 10,
        trd_ymax1: int = 160,
        trd_ymin2: int = 350,
        trd_ymax2: int = 500,
        fontz: float = 12.0,
        figsize: tuple = (16, 10.5),
        legend_loc: str = "lower left",
        output_pdf: Optional[str] = None
    ) -> Dict[str, dict]:

        if exp is None:
            exp = self.experiments[exp_index]
        caseid = self.case_ids[min(exp_index, len(self.case_ids)-1)]
        label = self.exp_labels[min(exp_index, len(self.exp_labels)-1)]

        _, _, landareaC = self._load_land_mask()

        self._apply_liao_thresholds(landareaC)

        if var_dict is None:
            var_dict = [
                ("TOTECOSYSC", (2850, 2960), (-8, 8)),
                ("TOTSOMC", (2000, 2062), (-0.6, 0.6)),
                ("TOTVEGC", (710, 750), (-6, 6)),
                ("TLAI", (1.20, 1.50), (-0.2, 0.2)),
                ("GPP", (114, 128), (-6, 6)),
                ("TWS", (9.0, 9.5), (-0.02, 0.02)),
                ("H2OSNO", (112.5, 126.5), (-3.0, 3.0)),
                ("HTOP", (6.0, 7.0), (-0.02, 0.02)),
                ("TSOI_10CM", (282, 285), (-0.8, 0.8)),
                ("SOILWATER_10CM", (30.0, 35.0), (-1.0, 1.0)),
                ("FSH", (32.0, 36.0), (-2.0, 2.0)),
                ("HCSOI", (19.0, 19.5), (-0.02, 0.02)),
            ]

        ordered_vars = [v[0] for v in var_dict]

        results: Dict[str, dict] = {}
        for v in ordered_vars:
            years, vals, dels, units = self._prep_series(
                v, exp, exp_index, landareaC, combine, engine,
                subper=subper,
            )
            results[v] = {
                "years": years,
                "values": vals,
                "deltas": dels,
                "year_mid": self._year_centers(years),
                "units": units,
                "title": self.var_setup[v]["title"],
                "thresh": self.thresholds[v],
                "trend_mid": np.array([], float),
                "trend": np.array([], float),
            }
            # Add BGC spin-up diagnostics for carbon stocks (series is PgC already)
            if self.bgc_focus and v in ("TOTECOSYSC", "TOTSOMC"):
                results[v]["bgc_metrics"] = self._bgc_spinup_metrics(
                    results[v]["years"], results[v]["values"], v
                )
            else:
                results[v]["bgc_metrics"] = {}

        # Align years across variables and compute right-axis series
        years_lists = [r["years"] for r in results.values() if len(r["years"]) > 0]
        if years_lists:
            y_min = max(arr.min() for arr in years_lists)
            y_max = min(arr.max() for arr in years_lists)
            for r in results.values():
                mask = (r["years"] >= y_min) & (r["years"] <= y_max)
                r["years"] = r["years"][mask]
                r["values"] = r["values"][mask]

                if self.use_trend:
                    mid, tr = self._rolling_trend(
                        r["years"], r["values"], window=self.trend_window, min_frac=self.trend_min_frac
                    )
                    r["trend_mid"] = mid; r["trend"] = tr
                    r["deltas"] = np.array([], float); r["year_mid"] = np.array([], float)
                else:
                    r["deltas"] = np.diff(r["values"])
                    r["year_mid"] = 0.5 * (r["years"][:-1] + r["years"][1:])
                    r["trend_mid"] = np.array([], float); r["trend"] = np.array([], float)

        self._plot_panels(
            order=var_dict,
            res=results,
            caseid=caseid,
            label=label,
            fontz=fontz,
            figsize=figsize,
            legend_loc=legend_loc,
            seg_ymin1=trd_ymin1,
            seg_ymax1=trd_ymax1,
            seg_ymin2=trd_ymin2,
            seg_ymax2=trd_ymax2,
            output_pdf=output_pdf
        )

        return results

    def _fmt_slope(self, x):
        if not np.isfinite(x):
            return "—"
        absx = abs(x)
        if absx >= 1e2 or absx < 1e-2:
            return f"{x:.2e}"
        elif absx < 10:
            return f"{x:.2f}".rstrip("0").rstrip(".")
        else:
            return f"{x:.1f}".rstrip("0").rstrip(".")

    def _add_grid(self, ax):
        ax.grid(which="major", linestyle="--", linewidth=0.5, alpha=0.5)

    def _plot_panels(
        self, 
        order: dict,
        res: Dict[str, dict], 
        caseid: str, 
        label: str, 
        output_pdf: str,
        fontz: float = 12.0, 
        figsize: tuple = (24, 10.5), 
        legend_loc: str = "lower left",
        seg_ymin1: int = 10,
        seg_ymax1: int = 500,
        seg_ymin2: int = 2200,
        seg_ymax2: int = 2500,
        x1_min: float = 0.0,
        x1_max: float = 500.0,
        x2_min: float = 2001.0,
        x2_max: float = 2500.0,
        x1_step: float = 100.0,
        x2_step: float = 100.0,
        add_grid: bool = False,
        add_xgrid: bool = False,
    ) -> None:
            
        nvars = len(order)
        ncols = 3
        nrows = (nvars + ncols - 1) // ncols
        letters = list(string.ascii_lowercase)
    
        series_color = "black"
        right_color  = "tab:blue"
    
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, constrained_layout=True)
        axes = np.atleast_2d(axes)
        fig.set_constrained_layout_pads(
            wspace=0.15,  # relative spacing (column gap)
            hspace=0.1,  # relative spacing (row gap)
        )
        
        # small gap between the two broken pieces
        gap = 0.04
        left_w  = 0.5 - gap / 2.0
        right_w = left_w
    
        for i, (v, ylims, dlims) in enumerate(order):
            row, col = divmod(i, ncols)
            container = axes[row, col]
    
            # container is just a holder (no frame/ticks)
            container.set_frame_on(False)
            container.tick_params(
                left=False, right=False, bottom=False, top=False,
                labelleft=False, labelbottom=False
            )
    
            # two equal-width inset axes inside the container
            ax_left  = container.inset_axes([0.0, 0.0, left_w, 1.0])
            ax_right = container.inset_axes([0.5 + gap / 2.0, 0.0, right_w, 1.0])
    
            if add_grid:
                self._add_grid(ax_left)
    
            ax_left2  = ax_left.twinx()
            ax_right2 = ax_right.twinx()
            # matplotlib >= 3.11 gives twins of inset axes an opaque patch that hides the series
            ax_left2.patch.set_visible(False)
            ax_right2.patch.set_visible(False)

            # --- turn OFF the right-side ticks/labels on the LEFT broken part ---
            ax_left2.tick_params(
                axis="y",
                right=False,
                labelright=False,
                length=8,          # <-- removes tick marks completely
                which="both"
            )
            ax_left2.spines["right"].set_visible(False)
            
            # --- turn OFF the left-side ticks/labels on the RIGHT broken part ---
            ax_right2.tick_params(
                axis="y",
                left=False,
                labelleft=False,
                length=8,          # <-- removes tick marks completely
                which="both"
            )
            ax_right2.spines["left"].set_visible(False)
            
            # left broken panel → hide ticks on the RIGHT edge
            ax_left.tick_params(
                axis="y",
                which="both",   # major + minor
                right=False,
                length=8
            )
            
            # right broken panel → hide ticks on the LEFT edge
            ax_right.tick_params(
                axis="y",
                which="both",
                left=False,
                length=8
            )
            
            # --- make right axes (trend) red on BOTH parts ---
            for a2 in (ax_left2, ax_right2):
                a2.tick_params(axis="y", colors=right_color)
                a2.spines["right"].set_color(right_color)
    
            r = res[v]
            years = r["years"]
            vals  = r["values"]
    
            panel_str = f"({letters[i]})"
            ax_left.set_title(f"{panel_str} {v}", fontsize=fontz, loc="left")
    
            # ===== main series (annual means), clipped =====
            ln1_left = ln1_right = None
    
            if years.size and vals.size:
                mask_left  = (years >= x1_min) & (years <= x1_max)
                mask_right = (years >= x2_min) & (years <= x2_max)
    
                if np.any(mask_left):
                    ln1_left = ax_left.plot(
                        years[mask_left], vals[mask_left],
                        color=series_color, linewidth=2.0, label="Annual Mean"
                    )[0]
                    ax_left.set_ylim(*ylims)
    
                if np.any(mask_right):
                    ln1_right = ax_right.plot(
                        years[mask_right], vals[mask_right],
                        color=series_color, linewidth=2.0
                    )[0]
                    ax_right.set_ylim(*ylims)
            else:
                ax_left.text(
                    0.5, 0.5, "No series data",
                    ha="center", va="center",
                    fontsize=fontz * 0.95,
                    transform=ax_left.transAxes,
                )
    
            # ===== trend / delta (right axis), clipped =====
            right_x = r.get("trend_mid", np.array([]))
            right_y = r.get("trend",     np.array([]))
            have_trend = right_x.size and right_y.size
            if not have_trend:
                right_x = r.get("year_mid", np.array([]))
                right_y = r.get("deltas",   np.array([]))
    
            ln2_left = ln2_right = None
            if right_x.size and right_y.size:
                mask_left_r  = (right_x >= x1_min) & (right_x <= x1_max)
                mask_right_r = (right_x >= x2_min) & (right_x <= x2_max)
    
                if np.any(mask_left_r):
                    ln2_left = ax_left2.plot(
                        right_x[mask_left_r], right_y[mask_left_r],
                        color=right_color, linestyle="dotted", linewidth=2.0,
                        label=f"Trend({self.trend_window}yr)" if have_trend else f"Δ{v}"
                    )[0]
    
                if np.any(mask_right_r):
                    ln2_right = ax_right2.plot(
                        right_x[mask_right_r], right_y[mask_right_r],
                        color=right_color, linestyle="dotted", linewidth=2.0,
                    )[0]
    
                for a2 in (ax_left2, ax_right2):
                    a2.axhline(0.0, color="gray", linewidth=0.5)
                    th = r.get("thresh", np.nan)
                    if np.isfinite(th):
                        a2.axhline(th,  linestyle="--", linewidth=1.0, color="blue", alpha=0.7)
                        a2.axhline(-th, linestyle="--", linewidth=1.0, color="blue", alpha=0.7)
                    a2.set_ylim(*dlims)

            # ===== x-limits, ticks, broken cosmetics =====
            # lock the ranges exactly
            ax_left.set_xlim(x1_min, x1_max)    # e.g. 0   → 500
            ax_right.set_xlim(x2_min, x2_max)   # e.g. 2001 → 2500
            ax_left.set_autoscale_on(False)
            ax_right.set_autoscale_on(False)
            ax_left.set_xmargin(0.0)
            ax_right.set_xmargin(0.0)

            # explicit tick locations so we *start* at 0 and *end* at 2500
            if x1_step is not None and x1_step > 0:
                left_ticks = np.arange(x1_min, x1_max + 0.5 * x1_step, x1_step)
                # e.g. [0, 100, 200, 300, 400, 500]
                ax_left.set_xticks(left_ticks)

            if x2_step is not None and x2_step > 0:
                right_ticks = np.arange(x2_min, x2_max + 0.5 * x2_step, x2_step)
                # e.g. [2000/2001, 2100, 2200, 2300, 2400, 2500]
                ax_right.set_xticks(right_ticks)

            # optional: keep grid aligned with ticks if requested
            if add_xgrid:
                ax_left.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.5)
                ax_right.grid(axis="x", linestyle="--", linewidth=0.5, alpha=0.5)

            # y-labels: attach series units to the container (one per panel),
            # and trend units to the RIGHT broken part only
            series_units, delta_units = self._split_units_for_display(v, r["units"])
            ax_left.set_ylabel(
                series_units, fontsize=fontz * 0.95, 
                color=series_color, rotation=90, labelpad=15
            )
    
            ax_right2.set_ylabel(
                delta_units, fontsize=fontz * 0.95,
                color=right_color, rotation=270, labelpad=15
            )
            
            ax_left.tick_params(labelsize=fontz * 0.95)
            ax_right.tick_params(labelsize=fontz * 0.95)
            ax_left2.tick_params(labelsize=fontz * 0.95)
            ax_right2.tick_params(labelsize=fontz * 0.95)
            
            ax_left.spines["right"].set_visible(False)
            ax_right.spines["left"].set_visible(False)
            ax_left.tick_params(labelright=False)
            ax_right.tick_params(labelleft=False)
    
            # diagonal break markers
            d = 0.015
            kwargs = dict(color="k", clip_on=False)
            ax_left.plot((1 - d, 1 + d), (-d, +d), transform=ax_left.transAxes, **kwargs)
            ax_left.plot((1 - d, 1 + d), (1 - d, 1 + d), transform=ax_left.transAxes, **kwargs)
            ax_right.plot((-d, +d), (-d, +d), transform=ax_right.transAxes, **kwargs)
            ax_right.plot((-d, +d), (1 - d, 1 + d), transform=ax_right.transAxes, **kwargs)
    
            # ===== trend annotation text =====
            try:
                seg_year1 = seg_ymax1
                seg_year2 = seg_ymin2
    
                trend_mid  = r.get("trend_mid", np.array([]))
                trend_vals = r.get("trend",     np.array([]))
    
                if trend_mid.size and trend_vals.size:
                    pre_mask  = (trend_mid < seg_ymax1) & (trend_mid > seg_ymin1)
                    post_mask = (trend_mid >= seg_ymin2) & (trend_mid < seg_ymax2)
                    pre_valid  = pre_mask  & np.isfinite(trend_vals)
                    post_valid = post_mask & np.isfinite(trend_vals)
                    slope_pre  = np.nanmean(trend_vals[pre_valid])  if np.any(pre_valid)  else np.nan
                    slope_post = np.nanmean(trend_vals[post_valid]) if np.any(post_valid) else np.nan
                else:
                    yrs  = r["years"]
                    vals = r["values"]
                    pre_idx  = (yrs < seg_year1)  & np.isfinite(yrs) & np.isfinite(vals)
                    post_idx = (yrs >= seg_year2) & np.isfinite(yrs) & np.isfinite(vals)
                    slope_pre  = self._theil_sen_slope(yrs[pre_idx],  vals[pre_idx])  if np.count_nonzero(pre_idx)  >= 2 else np.nan
                    slope_post = self._theil_sen_slope(yrs[post_idx], vals[post_idx]) if np.count_nonzero(post_idx) >= 2 else np.nan
    
                s_pre  = self._fmt_slope(slope_pre)
                s_post = self._fmt_slope(slope_post)
    
                #note = (
                #    f"trend<{int(seg_year1)} = {s_pre} {delta_units}\n"
                #    f"trend≥{int(seg_year2)} = {s_post} {delta_units}"
                #)
                note = f"Long-term Trend: ({s_pre}, {s_post}) {delta_units}"
                
                # put text inside the left data panel, just under the top edge
                container.text(
                    0.05, 0.95, note,
                    transform=container.transAxes,
                    ha="left", va="top",
                    fontsize=fontz * 0.95,
                    color="red",
                    linespacing=1.0,
                    fontstyle="italic",
                    clip_on=False,
                    zorder=20,
                    bbox=dict(
                        boxstyle="round",   # "square", "round", "round4", etc.
                        fc="white",         # facecolor
                        ec="black",         # edgecolor (frame color)
                        lw=1.0,             # line width
                        alpha=0.9,          # transparency
                        pad=0.25            # padding inside the box
                    ),
                )
        
            except Exception:
                pass
                
            # legend across the full panel (container-level)
            handles, labels_list = [], []
            if ln1_left is not None:
                handles.append(ln1_left)
                labels_list.append(ln1_left.get_label())
            if ln2_left is not None:
                handles.append(ln2_left)
                labels_list.append(ln2_left.get_label())
            
            if handles:
                container.legend(
                    handles, labels_list,
                    loc="lower left",
                    bbox_to_anchor=(0.02, -0.02),   # inside but below data area
                    bbox_transform=container.transAxes,
                    ncol=len(handles),              # <-- forces a single row
                    fontsize=fontz * 0.9,
                    frameon=False,
                )
                
        # hide unused containers
        for j in range(nrows * ncols):
            if j >= nvars:
                row, col = divmod(j, ncols)
                axes[row, col].axis("off")
    
        with PdfPages(output_pdf) as pdf:
            fig.canvas.draw()
            pdf.savefig(fig, bbox_inches="tight", pad_inches=0.1)
        plt.show()
        plt.close(fig)

    def _theil_sen_slope(self, years: np.ndarray, values: np.ndarray) -> float:
        y = years.astype(float); v = values.astype(float)
        m = np.isfinite(y) & np.isfinite(v)
        y, v = y[m], v[m]
        n = y.size
        if n < 3:
            return np.nan
        Yi = y[:, None]
        Yj = y[None, :]
        Vi = v[:, None]
        Vj = v[None, :]
        mask = np.triu(np.ones((n, n), dtype=bool), k=1)
        dy = (Yj - Yi)[mask]
        dv = (Vj - Vi)[mask]
        ok = np.isfinite(dy) & np.isfinite(dv) & (dy != 0)
        if not np.any(ok):
            return np.nan
        return float(np.median(dv[ok] / dy[ok]))

    def _running_mean(self, a: np.ndarray, w: int) -> np.ndarray:
        if a.size == 0 or w < 1 or w > a.size:
            return np.full_like(a, np.nan, dtype=float)
        c = np.cumsum(np.insert(a, 0, 0.0))
        rm = (c[w:] - c[:-w]) / float(w)
        pad = (w - 1) // 2
        out = np.full_like(a, np.nan, dtype=float)
        out[pad:pad+rm.size] = rm
        return out

    def _exp_e_fold_tau(self, years: np.ndarray, values: np.ndarray) -> float:
        """Fit v(t) ~ A + B * exp(-(t - t0)/tau) by OLS on log|v-A|; returns tau (yrs)."""
        y = years.astype(float); v = values.astype(float)
        m = np.isfinite(y) & np.isfinite(v)
        y, v = y[m], v[m]
        if y.size < 10:
            return np.nan
        k = max(10, int(0.2 * y.size))  # last 20% window (>=10)
        A = float(np.nanmean(v[-k:]))   # asymptote
        r = np.abs(v - A)
        idx = np.where(r > 0)[0]
        if idx.size < 6:
            return np.nan
        t, rr = y[idx], r[idx]
        q20 = np.quantile(rr, 0.2)      # drop tiny tail near asymptote
        use = rr > q20
        t, rr = t[use], rr[use]
        if t.size < 6:
            return np.nan
        X = np.vstack([np.ones_like(t), -t]).T  # log(rr)=log|B| - t/tau
        beta, *_ = np.linalg.lstsq(X, np.log(rr), rcond=None)
        slope = beta[1]
        tau = (1.0 / slope) if slope != 0 else np.nan
        return float(tau) if tau > 0 else np.nan

    def _bgc_spinup_metrics(self, years: np.ndarray, series_PgC: np.ndarray, var: str) -> dict:
        """
        Spin-up diagnostics for stocks (assumes series is PgC).
        Returns: sen slope (PgC/yr), 30y slope (PgC/yr), late rel error, tau, meets_criteria.
        """
        out = {"sen_slope_PgC_per_yr": np.nan,
               "slope30yr_PgC_per_yr": np.nan,
               "rel_err_late_mean": np.nan,
               "tau_years": np.nan,
               "meets_criteria": False}
        if years.size < 5 or series_PgC.size < 5:
            return out

        out["sen_slope_PgC_per_yr"] = self._theil_sen_slope(years, series_PgC)

        W = min(self.trend_window, max(5, years.size // 4))
        rm = self._running_mean(series_PgC, W)
        valid = np.isfinite(rm)
        if valid.sum() >= 2:
            yv, vv = years[valid], rm[valid]
            x = yv - yv.mean(); den = np.dot(x, x)
            out["slope30yr_PgC_per_yr"] = float(np.dot(x, vv - vv.mean()) / den) if den > 0 else np.nan

        k = max(10, int(0.2 * series_PgC.size))
        late_mean = float(np.nanmean(series_PgC[-k:]))
        if np.isfinite(late_mean) and late_mean != 0:
            out["rel_err_late_mean"] = float(np.nanmean(np.abs(series_PgC[-k:] - late_mean)) / abs(late_mean))

        out["tau_years"] = self._exp_e_fold_tau(years, series_PgC)

        th = self.bgc_thresholds
        sen_ok = (abs(out["sen_slope_PgC_per_yr"]) <= th["TOTSOMC_slope"] if var == "TOTSOMC"
                  else abs(out["sen_slope_PgC_per_yr"]) <= th["TOTECOSYSC_slope"])
        rel_ok = out["rel_err_late_mean"] <= th["tot_stock_rel_err"]
        tau_ok = (np.isfinite(out["tau_years"]) and out["tau_years"] <= th["tau_max_yrs"])
        out["meets_criteria"] = bool(sen_ok and rel_ok and tau_ok)
        return out

    def _load_processed_series(
        self,
        exp: str,
        var: str,
        processed_dir: str,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, str]:
        """
        Load preprocessed annual series and deltas from disk
        (written by save_processed_series).
        """
        fname = self._processed_path(exp, var, processed_dir)
        if not fname.exists():
            raise FileNotFoundError(f"Processed file not found: {fname}")
        else:
            print(f"load processed data: {fname}")

        ds = xr.open_dataset(fname)
        ann_name   = f"{var}_ann"
        delta_name = f"{var}_delta"

        years = ds["time"].values.astype(float)
        vals  = ds[ann_name].values.astype(float)

        if getattr(self, "verbose", False):
            print(f"[LOAD PROCESSED] {exp} / {var}: "
                  f"years(size={years.size}), vals(size={vals.size})")

        if delta_name in ds:
            dels = ds[delta_name].values.astype(float)
        else:
            dels = np.array([], float)

        units = str(ds[ann_name].attrs.get("units", "")) or str(
            self.var_setup.get(var, {}).get("units", "")
        )

        ds.close()
        return years, vals, dels, units

    def analyze_combined_fullcpl_picontrol(
        self,
        var_dict: list,
        full_index: int = 0,
        pi_index: int = 1,
        combine: str = "by_coords",
        engine: str = "netcdf4",
        subper: int = 1,
        trd_ymin1: int = 10,
        trd_ymax1: int = 160,
        trd_ymin2: int = 350,
        trd_ymax2: int = 500,
        ax_ymin1: int = 0,
        ax_ymax1: int = 500,
        ax_ymin2: int = 2001,
        ax_ymax2: int = 2500,
        x1_step: float = 100.0,
        x2_step: float = 100.0,
        add_grid: bool = False,
        add_xgrid: bool = False,
        fontz: float = 12.0,
        figsize: tuple = (16, 10.5),
        legend_loc: str = "lower left",
        output_pdf: Optional[str] = None,
        processed_dir: Optional[str] = None,
    ) -> Dict[str, dict]:
        """
        Build a combined time series for each variable:
          - Full-CPL segment treated as model years 0 .. (N_full-1)
          - piControl segment treated as model years 2001 .. (2001+N_pi-1)
        """

        full_exp = self.experiments[full_index]
        pi_exp   = self.experiments[pi_index]

        full_label = self.exp_labels[min(full_index, len(self.exp_labels) - 1)]
        pi_label   = self.exp_labels[min(pi_index,   len(self.exp_labels) - 1)]

        caseid = f"{full_label} + {pi_label}"
        label  = caseid

        # Land mask + Liao thresholds (landarea only needed for threshold scaling)
        _, _, landareaC = self._load_land_mask()
        self._apply_liao_thresholds(landareaC)

        results: Dict[str, dict] = {}

        for (v, ylims, dlims) in var_dict:
            # --- Load Full-CPL & piControl series ---
            if processed_dir:
                years_f, vals_f, _, units_f = self._load_processed_series(
                    full_exp, v, processed_dir
                )
                years_p, vals_p, _, units_p = self._load_processed_series(
                    pi_exp, v, processed_dir
                )
            else:
                years_f, vals_f, _, units_f = self._prep_series(
                    v, full_exp, full_index, landareaC, combine, engine,
                    subper=subper,
                )
                years_p, vals_p, _, units_p = self._prep_series(
                    v, pi_exp, pi_index, landareaC, combine, engine,
                    subper=subper,
                )

            # If neither has data, store empty record
            if years_f.size == 0 and years_p.size == 0:
                results[v] = {
                    "years": np.array([], float),
                    "values": np.array([], float),
                    "deltas": np.array([], float),
                    "year_mid": np.array([], float),
                    "units": units_f or units_p,
                    "title": self.var_setup[v]["title"],
                    "thresh": self.thresholds.get(v, np.nan),
                    "trend_mid": np.array([], float),
                    "trend": np.array([], float),
                    "bgc_metrics": {},
                }
                continue

            # --- Synthetic model-year axes using actual year indices ---
            # Full-CPL: map (e.g.) 1..501 → 0..500
            if years_f.size:
                year0_f = float(years_f.min())
                new_years_f = years_f.astype(float) - year0_f
            else:
                new_years_f = np.array([], float)

            # piControl: map (e.g.) 1..500 → 2001..2500
            if years_p.size:
                new_years_p = years_p.astype(float) + 2000.0
            else:
                new_years_p = np.array([], float)

            if self.verbose:
                print(f"[COMBINE] {v}: Full-CPL years {years_f.min()}–{years_f.max()} "
                      f"→ {new_years_f.min() if new_years_f.size else 'NA'}–"
                      f"{new_years_f.max() if new_years_f.size else 'NA'}")
                print(f"[COMBINE] {v}: piControl years {years_p.min()}–{years_p.max()} "
                      f"→ {new_years_p.min() if new_years_p.size else 'NA'}–"
                      f"{new_years_p.max() if new_years_p.size else 'NA'}")

            # --- Combine segments ---
            combined_years = np.concatenate([new_years_f, new_years_p])
            combined_vals  = np.concatenate([vals_f,      vals_p])

            if combined_years.size:
                order_idx = np.argsort(combined_years)
                combined_years = combined_years[order_idx]
                combined_vals  = combined_vals[order_idx]

            # --- Right-axis: trend or deltas ---
            if self.use_trend and combined_years.size >= 2:
                trend_mid, trend_vals = self._rolling_trend(
                    combined_years, combined_vals,
                    window=self.trend_window,
                    min_frac=self.trend_min_frac,
                )
                deltas   = np.array([], float)
                year_mid = np.array([], float)
            else:
                trend_mid  = np.array([], float)
                trend_vals = np.array([], float)
                if combined_years.size >= 2:
                    deltas   = np.diff(combined_vals)
                    year_mid = 0.5 * (combined_years[:-1] + combined_years[1:])
                else:
                    deltas   = np.array([], float)
                    year_mid = np.array([], float)

            series_units = units_f or units_p or self.var_setup[v]["units"]

            # Store in results
            results[v] = {
                "years": combined_years,
                "values": combined_vals,
                "deltas": deltas,
                "year_mid": year_mid,
                "units": series_units,
                "title": self.var_setup[v]["title"],
                "thresh": self.thresholds.get(v, np.nan),
                "trend_mid": trend_mid,
                "trend": trend_vals,
            }

            # BGC metrics on combined series
            if self.bgc_focus and v in ("TOTECOSYSC", "TOTSOMC"):
                results[v]["bgc_metrics"] = self._bgc_spinup_metrics(
                    results[v]["years"], results[v]["values"], v
                )
            else:
                results[v]["bgc_metrics"] = {}

        # Plot using existing panel routine
        self._plot_panels(
            order=var_dict,
            res=results,
            caseid=caseid,
            label=label,
            fontz=fontz,
            figsize=figsize,
            legend_loc=legend_loc,
            seg_ymin1=trd_ymin1,
            seg_ymax1=trd_ymax1,
            seg_ymin2=trd_ymin2,
            seg_ymax2=trd_ymax2,
            x1_min=ax_ymin1,
            x1_max=ax_ymax1,
            x2_min=ax_ymin2,
            x2_max=ax_ymax2,
            x1_step=x1_step,
            x2_step=x2_step,
            add_xgrid=add_xgrid,
            add_grid=add_grid,
            output_pdf=output_pdf,
        )

        return results


@dataclass
class LandCompareStability(SpinupStabilityBGC):
    """Land BGC comparison of Full-CPL vs alternate spin-ups."""

    show_monthly: bool = True

    monthly_alpha: float = 0.30

    monthly_lw: float = 0.6

    show_filtered: bool = True

    filt_window_months: int = 12    # rolling window in months

    filtered_lw: float = 2.0

    show_variability: bool = True

    band_alpha: float = 0.18

    sigma_mult: float = 1.0         # shade = mean ± sigma_mult * std

    filter_mode: str = "deseasoned" # "deseasoned" | "raw"

    remove_season: bool = True

    var_method: str = "rolling"

    var_window: int = 30

    anomalies: bool = False

    def __post_init__(self):
        assert len(self.experiments) == len(self.time_strs), "Each experiment needs a matching time_strs span"

    def _parse_range_from_fname(self, fname: str) -> tuple[int, int, int, int]:
        _FN_RANGE = re.compile(r'_(\d{6})_(\d{6})\.nc(?:4|\.gz|\.bz2|\.zst)?$')
        m = _FN_RANGE.search(Path(fname).name)
        if not m: return (0, 0, 0, 0)
        y0m0, y1m1 = m.group(1), m.group(2)
        return (int(y0m0[:4]), int(y0m0[4:]), int(y1m1[:4]), int(y1m1[4:]))

    def _var_files(self, exp: str, var: str) -> List[str]:
        if self.exp_subdirs and exp in self.exp_subdirs:
            d = Path(self.root_dir) / exp / self.exp_subdirs[exp]
        else:
            d = Path(self.root_dir) / exp / "post" / self.sub_dir
        files = sorted([str(p) for p in d.glob(f"{var}*.nc*")], key=lambda f: self._parse_range_from_fname(f))
        if not files:
            if var == "H2OSNO": return []  # allow missing snow
            raise FileNotFoundError(f"No files found for var={var} in {d}")
        return files

    def _anchor_monthly_time(self, tvals) -> xr.DataArray:
        if len(tvals) == 0: return xr.DataArray(tvals, dims=["time"], name="time")
        is_cf = (cftime is not None) and hasattr(tvals[0], "__class__") and ("cftime" in type(tvals[0]).__module__)
        def _snap_dt(dt):
            if is_cf:
                y, m = dt.year, dt.month
                if self.time_anchor == "start":  d = 1
                elif self.time_anchor == "middle": d = (self._days_in_month_noleap(m) + 1)//2
                else: d = self._days_in_month_noleap(m)
                return cftime.DatetimeNoLeap(y, m, d)
            ts = pd.to_datetime(dt); m = ts.month
            d = 1 if self.time_anchor=="start" else ((self._days_in_month_noleap(m)+1)//2 if self.time_anchor=="middle" else self._days_in_month_noleap(m))
            return ts.replace(day=d, hour=0, minute=0, second=0, microsecond=0)
        return xr.DataArray([_snap_dt(x) for x in tvals], dims=["time"], name="time")

    def _shift_one_month(self, dt, minus=True):
        sign = -1 if minus else +1
        def _anchor_day(year, month):
            if self.time_anchor == "start": return 1
            dim = self._days_in_month_noleap(month)
            return (dim + 1)//2 if self.time_anchor == "middle" else dim
        if (cftime is not None) and hasattr(dt, "__class__") and ("cftime" in type(dt).__module__):
            y, m = dt.year, dt.month
            m += sign
            if m == 0: m = 12; y -= 1
            elif m == 13: m = 1; y += 1
            return cftime.DatetimeNoLeap(y, m, _anchor_day(y, m))
        ts = pd.to_datetime(dt); y, m = ts.year, ts.month
        m += sign
        if m == 0: m = 12; y -= 1
        elif m == 13: m = 1; y += 1
        return ts.replace(year=y, month=m, day=_anchor_day(y, m), hour=0, minute=0, second=0, microsecond=0)

    def _maybe_fix_month_index(self, tvals):
        if self.month_index_correction == "as_is" or len(tvals) < 24: return tvals
        try:
            months = [int(getattr(t, "month", pd.to_datetime(t).month)) for t in tvals]
        except Exception:
            return tvals
        do_shift = (self.month_index_correction == "shift_back_1") or (
            self.month_index_correction == "auto" and months[0] == 2 and months[-1] == 1 and (len(months) % 12 == 0)
        )
        return [self._shift_one_month(t, minus=True) for t in tvals] if do_shift else tvals

    def _build_monthly_time_from_span(self, n: int, exp_index: int) -> xr.DataArray:
        span = self.time_strs[exp_index] if len(self.time_strs) else "0001-0001"
        y0 = int(span.split("-")[0]); m0 = 1
        if cftime is not None and self.calendar.lower() == "noleap":
            months = []; y, m = y0, m0
            for _ in range(n):
                months.append(cftime.DatetimeNoLeap(y, m, 1))
                m += 1;  y += m//13; m = 1 if m==13 else m
            return xr.DataArray(months, dims=["time"], name="time")
        vals = []; y, m = y0, m0
        for _ in range(n):
            vals.append(y + (m - 1) / 12.0)
            m += 1;  y += m//13; m = 1 if m==13 else m
        return xr.DataArray(vals, dims=["time"], name="time")

    def _ensure_or_construct_time(self, ds: xr.Dataset, files: List[str], exp_index: int) -> xr.Dataset:
        if "time" in ds or "Time" in ds:
            tname = "time" if "time" in ds else "Time"
            if tname != "time": ds = ds.rename({tname: "time"})
            if "time" in ds.coords and ds.sizes.get("time", 0) > 0:
                try:
                    if ds.sizes["time"] >= 24:
                        if self.replace_monthly_time_from_span:
                            ds = ds.assign_coords(time=self._build_monthly_time_from_span(ds.sizes["time"], exp_index))
                        else:
                            anchored = self._anchor_monthly_time(list(ds["time"].values))
                            fixed = self._maybe_fix_month_index(list(anchored.values))
                            ds = ds.assign_coords(time=xr.DataArray(fixed, dims=["time"], name="time"))
                except Exception:
                    pass
                if self.verbose:
                    n_t = int(ds.sizes.get("time", 0))
                    t0 = ds["time"].values[0] if n_t else "N/A"
                    t1 = ds["time"].values[-1] if n_t else "N/A"
                    print(f"[DEBUG] ensure_time: n={n_t} | span={t0} → {t1}")
                return ds

        tdim = next((cand for cand in ("time", "Time", "nTime", "t", "TimeMonthly", "Time_counter") if cand in ds.dims), None)
        if tdim is None: return ds
        n = ds.sizes[tdim]
        if self.force_time_from_span:
            time = self._build_monthly_time_from_span(n, exp_index)
        else:
            y0, m0, *_ = self._parse_range_from_fname(files[0])
            if y0 == 0:
                time = self._build_monthly_time_from_span(n, exp_index)
            else:
                if cftime is not None and self.calendar.lower() == "noleap":
                    months=[]; y, m = y0, m0
                    for _ in range(n):
                        months.append(cftime.DatetimeNoLeap(y, m, 1))
                        m += 1; y += m//13; m = 1 if m==13 else m
                    time = xr.DataArray(months, dims=["time"], name="time")
                else:
                    vals=[]; y, m = y0, m0
                    for _ in range(n):
                        vals.append(y + (m - 1) / 12.0)
                        m += 1; y += m//13; m = 1 if m==13 else m
                    time = xr.DataArray(vals, dims=["time"], name="time")
        if tdim != "time": ds = ds.rename({tdim: "time"})
        ds = ds.assign_coords(time=time)
        if self.verbose and n:
            print(f"[DEBUG] ensure_time: constructed monthly | start={str(time.values[0])} len={n}")
        return ds

    def _to_numeric_years(self, time_vals) -> np.ndarray:
        out = []
        for t in time_vals:
            if hasattr(t, "year") and hasattr(t, "month"):
                y, m = int(t.year), int(t.month)
                out.append(y + (m - 0.5) / 12.0)
            else:
                ts = pd.to_datetime(t)
                out.append(ts.year + (ts.month - 0.5) / 12.0)
        return np.array(out, dtype=float)

    def _load_land_mask(self) -> Tuple[xr.DataArray, xr.DataArray, float]:
        fr = xr.open_dataset(self.land_file)
        try:
            if "landfrac" not in fr or "area" not in fr:
                raise KeyError("land_file must contain 'landfrac' and 'area'")
            landfrac = fr["landfrac"].load()
            area = fr["area"].load()
            units = str(area.attrs.get("units", "")).lower()
            if any(k in units for k in ("km2","km^2","square kilometer","square kilometres","square kilometers")):
                area = area * 1e6
            elif not units and float(area.max()) < 1e7:
                area = area * 1e6
            landareaC = float((landfrac * area).sum())
            return landfrac, area, landareaC
        finally:
            fr.close()

    def _shared_period(self, exp_indices: Optional[Sequence[int]] = None) -> Tuple[int, int]:
        if exp_indices is None: exp_indices = list(range(len(self.experiments)))
        y0s, y1s = [], []
        for ei in exp_indices:
            a, b = self.time_strs[ei].split("-")
            y0s.append(int(a)); y1s.append(int(b))
        return max(y0s), min(y1s)

    def _read_series(self, var: str, exp: str, exp_index: int, combine: str, engine: str) -> xr.DataArray:
        files = self._var_files(exp, var)
        if not files: raise KeyError(f"{var} not available")
        engines = [engine] + [e for e in ("netcdf4","h5netcdf","scipy") if e != engine]
        last_err = None; ds = None
        for eng in engines:
            try:
                ds = xr.open_mfdataset(files, combine=combine, engine=eng, decode_times=False, parallel=False, chunks="auto")
                break
            except Exception as e: last_err = e
        if ds is None:
            raise RuntimeError(f"Failed to open {var} for {exp} with engines {engines}: {last_err}")
        ds = self._ensure_or_construct_time(ds, files, exp_index)
        if "time" not in ds: raise KeyError(f"Unable to find or construct 'time' for {var} in {exp}")
        da = ds[var]

        for cand in ["region","nRegions","nOceanRegions","nreg","ntrsc"]:
            if cand in da.dims:
                if self.region_index >= da.sizes[cand]:
                    raise IndexError(f"region_index={self.region_index} out of range for dim {cand} size={da.sizes[cand]}")
                da = da.isel({cand: self.region_index})
                break

        for bnd in ["bounds","bnds","nbnd","nbounds"]:
            if bnd in da.dims and da.sizes[bnd] >= 1: da = da.isel({bnd:0})

        non_time = [d for d in da.dims if d != "time"]
        if non_time:
            if len(non_time) == 1:
                d = non_time[0]
                da = da.squeeze(d) if da.sizes[d] == 1 else da.isel({d: self.region_index}) if (da.sizes[d]==3 and self.region_index<3) else da.mean(d)
            else:
                for d in list(non_time):
                    if da.sizes[d] == 1: da = da.squeeze(d)
                non_time = [d for d in da.dims if d != "time"]
                if non_time: da = da.mean(non_time)

        da = da.sortby("time").load()
        if self.verbose:
            n_t = int(da.sizes.get("time", 0))
            t0 = da["time"].values[0] if n_t else "N/A"
            t1 = da["time"].values[-1] if n_t else "N/A"
            print(f"[DEBUG] {var} ({exp}): loaded series | n={n_t} | span={t0} → {t1}")
        return da

    def _processed_path(self, exp_index: int, var: str, processed_dir: Optional[str]) -> Path:
        base = Path(processed_dir or "."); exp = self.experiments[exp_index]
        return base / f"{exp}_{var}_processed.nc"

    def _prep_monthly_series(self, var: str, exp: str, exp_index: int, landareaC: float,
                             combine: str, engine: str, subper: int = 1) -> Tuple[np.ndarray, np.ndarray, str]:
        try:
            da = self._read_series(var, exp, exp_index, combine, engine)
        except KeyError:
            return np.array([], float), np.array([], float), self.var_setup[var]["units"]
        
        scale = float(self.var_setup[var]["scale"])
        units = str(self.var_setup[var]["units"])
        
        y0, y1 = map(int, self.time_strs[exp_index].split("-"))
        t = da["time"]
        da = (da.where((t >= y0) & (t <= y1), drop=True) if np.issubdtype(t.dtype, np.number)
              else da.sel(time=slice(f"{y0:04d}-01-01", f"{y1:04d}-12-31")))
        if da.sizes.get("time", 0) == 0:
            return np.array([], float), np.array([], float), self.var_setup[var]["units"]

        years_num_mo = self._to_numeric_years(da["time"].values)
        vals = da.load().to_numpy().astype(float)
        
        # Conversions like NCL
        if var in ("TOTECOSYSC", "TOTSOMC", "TOTVEGC"):
            vals *= landareaC
        if var == "GPP":
            vals *= landareaC
            vals *= self.seconds_in_year
        if var == "TSOI_10CM":
            # convert from K to degC
            vals -= 273.15
        vals *= scale
        
        return years_num_mo, vals, units

    def _monthly_climatology(self, months: np.ndarray, vals: np.ndarray) -> np.ndarray:
        clim = np.full(12, np.nan)
        for m in range(1, 13):
            w = (months == m) & np.isfinite(vals)
            if np.any(w): clim[m-1] = float(np.nanmean(vals[w]))
        return clim

    def _rolling_mean_nan(self, x: np.ndarray, win: int, min_frac: float = 0.8) -> np.ndarray:
        n = x.size; win = max(1, min(int(win), n)); k = (win - 1)//2
        out = np.full(n, np.nan, dtype=float)
        for i in range(n):
            lo = max(0, i - k); hi = min(n, i + k + 1)
            seg = x[lo:hi]; ok = np.isfinite(seg)
            need = int(np.ceil(min_frac * (hi - lo)))
            if ok.sum() >= need and need > 0: out[i] = float(np.nanmean(seg[ok]))
        return out

    def _rolling_std_nan(self, x: np.ndarray, win: int, min_frac: float = 0.8) -> np.ndarray:
        n = x.size; win = max(1, min(int(win), n)); k = (win - 1)//2
        out = np.full(n, np.nan, dtype=float)
        for i in range(n):
            lo = max(0, i - k); hi = min(n, i + k + 1)
            seg = x[lo:hi]; ok = np.isfinite(seg)
            need = int(np.ceil(min_frac * (hi - lo)))
            if ok.sum() >= need and need > 0: out[i] = float(np.nanstd(seg[ok], ddof=1))
        return out

    def _running_mean_std(self, time_vals, vals, mode: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        months = np.array([self._month_of_safely(t) for t in time_vals])
        years_num = self._to_numeric_years(time_vals)
        x = np.array(vals, dtype=float, copy=False)
        if mode.lower() == "deseasoned":
            clim = self._monthly_climatology(months, x)
            x = x - clim[months - 1]
        mean12 = self._rolling_mean_nan(x, self.filt_window_months, min_frac=0.8)
        std12  = self._rolling_std_nan(x,  self.filt_window_months, min_frac=0.8)
        if mode.lower() == "deseasoned":
            mu = float(np.nanmean(vals)) if np.isfinite(vals).any() else 0.0
            mean12 = mean12 + mu
        return years_num, mean12, std12

    def _month_of_safely(self, t) -> int:
        try:
            return int(t.month)
        except Exception:
            return int(pd.to_datetime(t).month)

    def _load_figdata_one(self, var: str, figdata_dir: str) -> dict:
        # Match what _save_figure_data writes
        f_nc  = Path(figdata_dir) / f"land_ts_compare_figdata_{var}.nc"
        f_csv = Path(figdata_dir) / f"land_ts_compare_figdata_{var}.csv"

        if f_nc.exists():
            ds = xr.open_dataset(f_nc)
            name_mo = next((k for k in ds.data_vars if k.endswith("_monthly") and k.startswith(var)), None)
            name_rm = next((k for k in ds.data_vars if (k.endswith("_rollmean") or k.endswith("_filtered")) and k.startswith(var)), None)
            name_rs = next((k for k in ds.data_vars if (k.endswith("_rolling_std") or k.endswith("_rollstd")) and k.startswith(var)), None)
            if not (name_mo and name_rm):
                raise KeyError(f"{f_nc} missing monthly/rollmean for {var}")

            mo = ds[name_mo].values
            rm = ds[name_rm].values
            rs = ds[name_rs].values if name_rs else np.full_like(rm, np.nan)

            time = ds["time"].values.astype(float)
            exps = [str(x) for x in ds["experiment"].values.tolist()]
            units = str(ds[name_mo].attrs.get("units", ""))
            title = str(ds.attrs.get("var_title", var))
            ds.close()

            # Normalize keys for the rest of the pipeline
            return {"time": time, "experiments": exps, "monthly": mo,
                    "filtered": rm, "rollstd": rs, "units": units, "title": title}

        if f_csv.exists():
            df = pd.read_csv(f_csv)
            df["time"] = df["time"].astype(float)
            exps = df["experiment"].dropna().unique().tolist()
            time = np.sort(df["time"].dropna().unique())
            nexp, nt = len(exps), len(time)
            mo = np.full((nexp, nt), np.nan); rm = np.full((nexp, nt), np.nan); rs = np.full((nexp, nt), np.nan)
            exp_to_i = {e: i for i, e in enumerate(exps)}
            t_to_j = {t: j for j, t in enumerate(time)}

            for _, r in df.iterrows():
                i = exp_to_i.get(r["experiment"]); j = t_to_j.get(float(r["time"]))
                if i is None or j is None: continue
                if "monthly" in r:      mo[i, j] = r["monthly"]
                if "rollmean" in df:    rm[i, j] = r.get("rollmean", np.nan)
                elif "filtered" in df:  rm[i, j] = r.get("filtered", np.nan)
                if "rolling_std" in df: rs[i, j] = r.get("rolling_std", np.nan)
                elif "rollstd" in df:   rs[i, j] = r.get("rollstd", np.nan)

            units = self.var_setup.get(var, {}).get("units", "")
            title = self.var_setup.get(var, {}).get("title", var)
            return {"time": time, "experiments": exps, "monthly": mo,
                    "filtered": rm, "rollstd": rs, "units": units, "title": title}

        raise FileNotFoundError(f"No figdata found for {var} in {figdata_dir}")

    def load_figdata(self, figdata_dir: str, variables: list[str] | None = None) -> dict:
        if variables is None: variables = [v for v in self.var_setup.keys()]
        out = {}
        for v in variables:
            try:
                out[v] = self._load_figdata_one(v, figdata_dir)
            except FileNotFoundError:
                continue
        if not out: raise FileNotFoundError(f"No figdata files in {figdata_dir}")
        return out

    def _figdata_to_results(self, figdata: dict, exp_indices: Sequence[int],
                            order_for_plot: Sequence[Tuple[str, Optional[Tuple[float, float]]]]
                            ) -> Dict[int, Dict[str, dict]]:
        idx_to_label = {ei: self.exp_labels[min(ei, len(self.exp_labels)-1)] for ei in exp_indices}
        out: Dict[int, Dict[str, dict]] = {ei: {} for ei in exp_indices}
        for (var, _ylims) in order_for_plot:
            if var not in figdata:
                for ei in exp_indices:
                    out[ei][var] = {
                        "years_mo": np.array([], float), "values_mo": np.array([], float),
                        "values_filtered": np.array([], float), "values_std": np.array([], float),
                        "units": self.var_setup.get(var, {}).get("units", ""), "label": idx_to_label[ei],
                    }
                continue
            payload = figdata[var]
            time_union = np.asarray(payload["time"], dtype=float)
            exps_in_file = list(payload["experiments"])
            pos_by_label = {lab: i for i, lab in enumerate(exps_in_file)}
            monthly = np.asarray(payload["monthly"])
            filtered = np.asarray(payload["filtered"])
            rollstd = np.asarray(payload.get("rollstd", np.full_like(filtered, np.nan)))
            units = str(payload.get("units", self.var_setup.get(var, {}).get("units", "")))
            # cached under different legend labels (e.g. "Full-CPL" before it became "FC"):
            # fall back to plotting order when the experiment count matches
            by_position = (not all(idx_to_label[ei] in pos_by_label for ei in exp_indices)
                           and len(exps_in_file) == len(exp_indices))
            if by_position and self.verbose:
                print(f"[figdata] {var}: labels {exps_in_file} in file → matched by position")
            for epos, ei in enumerate(exp_indices):
                label = idx_to_label[ei]
                if by_position or label in pos_by_label:
                    pos = epos if by_position else pos_by_label[label]
                    ymo = time_union
                    vmo = monthly[pos, :].astype(float, copy=False)
                    vf  = filtered[pos, :].astype(float, copy=False)
                    vs  = rollstd[pos, :].astype(float, copy=False) if rollstd.ndim == 2 else np.full_like(vf, np.nan)
                else:
                    ymo = vmo = vf = vs = np.array([], float)
                out[ei][var] = {"years_mo": ymo, "values_mo": vmo, "values_filtered": vf,
                                "values_std": vs, "units": units, "label": label}
        return out

    def analyze_comparison(
        self,
        exp_indices: Optional[Sequence[int]] = None,
        var_dict: dict = None, 
        combine: str = "by_coords",
        engine: str = "netcdf4",
        subper: int = 1,
        fontz: float = 12.0,
        figsize: tuple = (16, 10.5),
        legend_loc: str = "lower left",
        output_pdf: Optional[str] = None,
        var_order: Optional[List[Tuple[str, Tuple[float, float]]]] = None,
        save_data: bool = False,
        save_data_dir: Optional[str] = None,
        save_fmt: str = "nc",
        force_recompute: bool = False,
        xlim: tuple = (220,350),
        xstep: int = 10, 
    ) -> Dict[int, Dict[str, dict]]:

        if exp_indices is None:
            exp_indices = list(range(len(self.experiments)))
        else:
            exp_indices = list(exp_indices)

        if var_dict is None:
            ordered_vars = [
                "TOTECOSYSC", "TOTSOMC", "TOTVEGC", "TLAI", "GPP", "TWS",
                "H2OSNO", "HTOP", "TSOI_10CM", "SOILWATER_10CM", "HCSOI", "FSH"
            ]
            order_for_plot = ([(v, None) for v in ordered_vars]
                              if var_order is None else var_order)
        else:
            ordered_vars = [v[0] for v in var_dict]

            order_for_plot = var_dict 
            

        # ---------- FAST PATH: try figdata ----------
        results_by_exp: Dict[int, Dict[str, dict]] = {}
        tried_figdata = False
        if save_data_dir and not force_recompute:
            try:
                var_list = [v for v, _ in order_for_plot]
                figdata = self.load_figdata(save_data_dir, variables=var_list)
                results_by_exp = self._figdata_to_results(figdata, exp_indices, order_for_plot)
                tried_figdata = True
                if self.verbose: print(f"[figdata] Loaded re-plot data from {save_data_dir}")
            except Exception as e:
                if self.verbose: print(f"[figdata] Could not load figdata ({e}); falling back to raw series.")

        # ---------- SLOW PATH: compute ----------
        if not results_by_exp:
            _, _, landareaC = self._load_land_mask()

            results_by_exp = {}
            for ei in exp_indices:
                exp = self.experiments[ei]
                label = self.exp_labels[min(ei, len(self.exp_labels)-1)]
                res_one: Dict[str, dict] = {}
                for v, _yl in order_for_plot:
                    ymo, vmo, units = self._prep_monthly_series(v, exp, ei, landareaC, combine, engine, subper=subper)
                    res_one[v] = {"years_mo": ymo, "values_mo": vmo, "units": units, "label": label}
                results_by_exp[ei] = res_one

            y_shared0, y_shared1 = self._shared_period(exp_indices)
            if y_shared0 > y_shared1:
                raise ValueError(f"No shared period across experiments {exp_indices}: intersection ({y_shared0}, {y_shared1}).")

            for v, _yl in order_for_plot:
                for ei in exp_indices:
                    exp = self.experiments[ei]
                    try:
                        da = self._read_series(v, exp, ei, combine, engine)
                    except KeyError:
                        r = results_by_exp[ei].get(v, {})
                        if r:
                            r["years_mo"] = np.array([], dtype=float)
                            r["values_mo"] = np.array([], dtype=float)
                            r["values_filtered"] = np.array([], dtype=float)
                            r["values_std"] = np.array([], dtype=float)
                        continue
                        
                    scale = float(self.var_setup[v]["scale"])
                    units = str(self.var_setup[v]["units"])
                    
                    if np.issubdtype(da["time"].dtype, np.number):
                        da = da.where((da["time"] >= y_shared0) & (da["time"] <= y_shared1), drop=True)
                    else:
                        da = da.sel(time=slice(f"{y_shared0:04d}-01-01", f"{y_shared1:04d}-12-31"))

                    if da.sizes.get("time", 0) == 0:
                        r = results_by_exp[ei][v]
                        r["years_mo"] = np.array([], dtype=float)
                        r["values_mo"] = np.array([], dtype=float)
                        r["values_filtered"] = np.array([], dtype=float)
                        r["values_std"] = np.array([], dtype=float)
                        continue

                    vals = da.load().to_numpy().astype(float)
                    if v in ("TOTECOSYSC", "TOTSOMC", "TOTVEGC"):
                        vals *= landareaC
                    if v == "GPP":
                        vals *= landareaC
                        vals *= self.seconds_in_year
                    if v == "TSOI_10CM" and np.nanmedian(vals) > 150:
                        vals -= 273.15
                    vals *= scale

                    t_native = da["time"].values
                    y_mo, mean12, std12 = self._running_mean_std(t_native, vals, mode=self.filter_mode)

                    r = results_by_exp[ei][v]
                    r["years_mo"] = y_mo
                    r["values_mo"] = vals
                    r["values_filtered"] = mean12
                    r["values_std"] = std12

            if save_data_dir is not None and not tried_figdata:
                self._save_figure_data(
                    res_by_exp=results_by_exp,
                    exp_indices=exp_indices,
                    order_for_plot=order_for_plot,
                    out_dir=save_data_dir,
                    fmt=save_fmt,
                )

        if output_pdf is None:
            output_pdf = "spinup_stability_runningmean_sigma.pdf"

        self._plot_panels_multi(
            res_by_exp=results_by_exp,
            exp_indices=exp_indices,
            order_for_plot=order_for_plot,
            fontz=fontz,
            figsize=figsize,
            legend_loc=legend_loc,
            output_pdf=output_pdf,
            xlim=xlim,
            xstep=xstep, 
        )

        return results_by_exp

    def _save_figure_data(
        self,
        res_by_exp: Dict[int, Dict[str, dict]],
        exp_indices: Sequence[int],
        order_for_plot: Sequence[Tuple[str, Optional[Tuple[float, float]]]],
        out_dir: str,
        fmt: str = "nc",
    ) -> None:
        outp = Path(out_dir); outp.mkdir(parents=True, exist_ok=True)
        exp_labels = {ei: res_by_exp[ei][next(iter(res_by_exp[ei]))]["label"]
                      for ei in exp_indices if res_by_exp.get(ei)}

        for var, _yl in order_for_plot:
            all_times = []
            for ei in exp_indices:
                r = res_by_exp[ei].get(var, {})
                if r.get("years_mo", np.array([])).size:
                    all_times.append(r["years_mo"])
            if not all_times:
                continue
            time_union = np.unique(np.concatenate(all_times)); time_union.sort()

            exps = [exp_labels[ei] for ei in exp_indices]
            monthly = np.full((len(exps), time_union.size), np.nan, dtype=float)
            rollmean = np.full((len(exps), time_union.size), np.nan, dtype=float)
            rollstd  = np.full((len(exps), time_union.size), np.nan, dtype=float)
            units   = ""

            for epos, ei in enumerate(exp_indices):
                r = res_by_exp[ei].get(var, {})
                if not r or r["years_mo"].size == 0:
                    continue
                t = r["years_mo"]; vmo = r["values_mo"]
                vf = r.get("values_filtered", np.full_like(vmo, np.nan))
                vs = r.get("values_std",      np.full_like(vmo, np.nan))
                units = r.get("units", units)

                idx = np.searchsorted(time_union, t)
                ok = (idx >= 0) & (idx < time_union.size) & np.isclose(time_union[idx], t, atol=1e-9)
                monthly[epos, idx[ok]] = vmo[ok]
                rollmean[epos, idx[ok]] = vf[ok]
                rollstd [epos, idx[ok]] = vs[ok]

            if fmt in ("nc", "both"):
                ds = xr.Dataset(
                    data_vars={
                        f"{var}_monthly":     (("experiment", "time"), monthly, {"long_name": f"{var} monthly", "units": units}),
                        f"{var}_rollmean":    (("experiment", "time"), rollmean, {"long_name": f"{var} 12-mo running mean", "units": units}),
                        f"{var}_rolling_std": (("experiment", "time"), rollstd,  {"long_name": f"{var} 12-mo running std",  "units": units}),
                    },
                    coords={"experiment": np.array(exps, dtype=object), "time": time_union},
                    attrs={
                        "note": "Figure-ready data: monthly, 12-mo running mean, 12-mo running std.",
                        "var_title": self.var_setup.get(var, {}).get("title", var),
                        "filter_mode": self.filter_mode,
                        "window_months": self.filt_window_months,
                    },
                )
                f_nc = outp / f"land_ts_compare_figdata_{var}.nc"
                ds.to_netcdf(f_nc); print(f"[figdata] wrote {f_nc}")

            if fmt in ("csv", "both"):
                rows = []
                for epos, label in enumerate(exps):
                    rows.append(pd.DataFrame({
                        "experiment": label,
                        "time": time_union,
                        "monthly": monthly[epos, :],
                        "rollmean": rollmean[epos, :],
                        "rolling_std": rollstd[epos, :],
                    }))
                df = pd.concat(rows, ignore_index=True)
                f_csv = outp / f"land_ts_compare_figdata_{var}.csv"
                df.to_csv(f_csv, index=False); print(f"[figdata] wrote {f_csv}")

    def save_processed_series(self, exp_index: int, out_dir: str, overwrite: bool = False,
                              combine: str = "by_coords", engine: str = "netcdf4") -> None:
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        exp = self.experiments[exp_index]
        _, _, landareaC = self._load_land_mask()

        for v in self.var_setup.keys():
            fname = self._processed_path(exp_index, v, out_dir)
            if fname.exists() and not overwrite:
                print(f"Skip existing {fname}"); continue

            da = self._read_series(v, exp, exp_index, combine, engine)
            ann = self._annualize(da, self.annual_hist)

            y0, y1 = map(int, self.time_strs[exp_index].split("-"))
            ann = (ann.where((ann["time"] >= y0) & (ann["time"] <= y1), drop=True)
                   if np.issubdtype(ann["time"].dtype, np.number)
                   else ann.sel(time=slice(f"{y0:04d}", f"{y1:04d}")))

            years = (ann["time"].values.astype(float) if np.issubdtype(ann["time"].dtype, np.number)
                     else np.array([getattr(t, "year", pd.Timestamp(t.astype("datetime64[ns]")).year)
                                    for t in ann["time"].values], dtype=float))

            vals = ann.load().to_numpy().astype(float)
            
            units = str(self.var_setup[v]["units"])
            ds_out = xr.Dataset({f"{v}_ann": (["time"], vals, {"units": units})},
                                coords={"time": years}, attrs={"experiment": exp})
            ds_out.to_netcdf(fname); print(f"Wrote {fname}")

    def _split_units_for_display(self, var: str, stored_units: str) -> tuple[str, str]:
        if var in ("TOTECOSYSC","TOTSOMC","TOTVEGC"): return "Pg C", "Pg C yr⁻¹"
        base = {"TLAI":"m² m⁻²","GPP":"Pg C yr⁻¹","TWS":"m","H2OSNO":"mm","HTOP":"m",
                "TSOI_10CM":"°C","SOILWATER_10CM":"kg m⁻²","FSH":"W m⁻²","HCSOI":"MJ m⁻²"}
        if var in base: series_units = base[var]
        else:
            u = stored_units
            u = u.replace("$^{-1}$","⁻¹").replace("^−1","⁻¹").replace("^−¹","⁻¹")
            u = re.sub(r'\s*[/ ]?yr(?:\s*[\-−]?\s*1|⁻¹)?\)?', '', u, flags=re.IGNORECASE)
            u = u.replace("()","").strip(); u = re.sub(r'\s+',' ',u).strip()
            series_units = u if u else stored_units
        return series_units, f"{series_units} yr⁻¹"

    def _plot_panels_multi(
        self,
        res_by_exp: Dict[int, Dict[str, dict]],
        exp_indices: Sequence[int],
        order_for_plot: Sequence[Tuple[str, Optional[Tuple[float, float]]]],
        output_pdf: str,
        fontz: float = 12.0,
        figsize: tuple = (16, 10.5),
        legend_loc: str = "lower left",
        xlim: tuple=(220,350), 
        xstep: int = 10 
    ) -> None:

        nvars = len(order_for_plot)
        ncols = 3
        nrows = (nvars + ncols - 1) // ncols

        base_colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd",
                       "#ff7f0e", "#17becf", "#8c564b"]
        exp_colors = {ei: base_colors[i % len(base_colors)] for i, ei in enumerate(exp_indices)}
        letters = list(string.ascii_lowercase)

        fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
        axes = np.atleast_2d(axes)

        for i, (v, ylims) in enumerate(order_for_plot):
            row, col = divmod(i, ncols)
            ax = axes[row, col]

            # --- titles / labels sized off fontz ---
            units = None
            for ei in exp_indices:
                units = res_by_exp[ei][v]["units"] or units
            series_units, _ = self._split_units_for_display(v, units or "")

            ax.set_title(f"({letters[i]}) {v}", fontsize=fontz, loc="left")
            ax.set_ylabel(series_units if not self.anomalies else f"{series_units} (anomaly)",
                          fontsize=fontz * 0.9)
            ax.tick_params(axis="both", labelsize=fontz * 0.9)

            if ylims is not None:
                try:
                    ax.set_ylim(*ylims)
                except Exception:
                    pass
                    
            if xlim is not None:
                try:
                    ax.set_xlim(*xlim)
                except Exception:
                    pass

            if xstep is not None and xstep > 0:
                x_min, x_max = xlim
                x_ticks = np.arange(x_min+2, x_max, xstep)
                ax.set_xticks(x_ticks)

            
            ax.grid(True, ls="--", lw=0.5, alpha=0.5)

            # --- plot data ---
            for ei in exp_indices:
                r = res_by_exp[ei][v]
                ymo = r.get("years_mo", np.array([]))
                vmo = r.get("values_mo", np.array([]))
                if ymo.size == 0:
                    continue

                color = exp_colors[ei]
                label = r.get("label", f"Exp{ei}")

                # monthly (faint; hidden from legend)
                if self.show_monthly and vmo.size:
                    ax.plot(ymo, vmo, color=color, lw=self.monthly_lw,
                            alpha=self.monthly_alpha, label="_monthly")

                # running mean + ±σ band
                if self.show_filtered and "values_filtered" in r:
                    m = r["values_filtered"]
                    ax.plot(ymo, m, color=color, lw=self.filtered_lw, label=label)
                    if self.show_variability and "values_std" in r and np.isfinite(r["values_std"]).any():
                        s = r["values_std"]
                        ax.fill_between(
                            ymo, m - self.sigma_mult * s, m + self.sigma_mult * s,
                            facecolor=color, alpha=self.band_alpha, linewidth=0, label="_band"
                        )

            # --- legend: one entry per experiment, hide _monthly/_band ---
            hs, ls = ax.get_legend_handles_labels()
            keep = [(h, l) for h, l in zip(hs, ls) if not l.startswith("_")]
            seen, uniq_h, uniq_l = set(), [], []
            for h, l in keep:
                if l in seen:
                    continue
                seen.add(l)
                uniq_h.append(h)
                uniq_l.append(l)

            # show legend only in first row of each column
            if row == 0 and uniq_l:
                old = ax.get_legend()
                if old is not None:
                    old.remove()
                ax.legend(
                    uniq_h,
                    uniq_l,
                    loc=legend_loc,
                    ncol=len(uniq_l),      # put all items in one row
                    fontsize=fontz * 0.9,
                    frameon=False
                )
                
        # turn off any empty panels
        for j in range(nrows * ncols - nvars):
            r, c = divmod(nvars + j, ncols)
            axes[r, c].axis("off")

        # --- adjust layout spacing ---
        plt.subplots_adjust(
            left=0.06,   # padding from left
            right=0.98,  # padding from right
            top=0.96,    # padding from top
            bottom=0.06, # padding from bottom
            wspace=0.2, # horizontal space between columns
            hspace=0.30  # vertical space between rows
        )

        with PdfPages(output_pdf) as pdf:
            pdf.savefig(fig, bbox_inches="tight", pad_inches=0.03)

        plt.show()
        plt.close(fig)
