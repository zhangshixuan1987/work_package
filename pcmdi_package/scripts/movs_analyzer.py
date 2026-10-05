"""
Mode-of-variability analysis pipeline.

Classes
-------
ModeConfigManager
    Loads and queries a JSON configuration describing datasets and modes.
BaseModeAnalysis
    Low-level helpers: I/O, pre-processing, EOF computation, area index,
    and projection onto an external EOF (common-base approach).
ModeAnalyzer
    High-level API: analyse obs / a single model / all models.
    For EOF modes each saved model Dataset contains both sets of variables:
      Independent EOF   : eof, pc, frac, slope, slope_pval, corr
      Common-base proj  : pc_proj, eof_lr_proj, frac_proj, slope_proj, slope_pval_proj, corr_proj
    The common basis is always the obs EOF solver (passed automatically via
    analyze_or_load_all → ref_solver from the observational EOF analysis).
"""

import os
import glob
import json
import sys
from importlib import import_module
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import xarray as xr
from scipy import stats

try:
    from eofs.xarray import Eof
except ModuleNotFoundError:
    Eof = None


class ModeConfigManager:
    """Load and query mode / dataset configuration from a JSON file or dict."""

    def __init__(self, config_source):
        if isinstance(config_source, dict):
            self.config = config_source
        else:
            with open(config_source, "r") as f:
                self.config = json.load(f)
        self.datasets: Dict = self.config["datasets"]
        self.modes: Dict    = self.config["modes"]

    def get_mode_info(self, mode_name: str) -> dict:
        if mode_name not in self.modes:
            raise KeyError(f"Mode '{mode_name}' not found in config.")
        return self.modes[mode_name]

    def get_dataset_info(self, case_name: str) -> dict:
        if case_name not in self.datasets:
            raise KeyError(f"Dataset '{case_name}' not found in config.")
        return self.datasets[case_name]

    def get_default_obs(self, mode_name: str) -> Optional[dict]:
        return self.get_mode_info(mode_name).get("obs")


class BaseModeAnalysis:
    """Low-level helpers shared by all mode analyses."""

    LOCAL_PCMDI_METRICS_ROOT = os.environ.get(
        "PCMDI_METRICS_SOURCE_ROOT",
        "/home/ac.szhang/code/pcmdi_metrics",
    )

    # ------------------------------------------------------------------
    # Coordinate normalisation
    # ------------------------------------------------------------------
    @staticmethod
    def normalize_lon(da: xr.DataArray) -> xr.DataArray:
        """Shift longitudes to [-180, 180) and sort."""
        if "lon" in da.coords:
            da = da.assign_coords(lon=((da.lon + 180) % 360) - 180).sortby("lon")
        return da

    @staticmethod
    def subset_time(da: xr.DataArray, period: Optional[Tuple[int, int]] = None) -> xr.DataArray:
        if period is None:
            return da
        y0, y1 = period
        return da.sel(time=slice(f"{y0}-01-01", f"{y1}-12-31"))

    @staticmethod
    def subset_latlon(
        da: xr.DataArray,
        lat_bnds: Optional[Tuple[float, float]] = None,
        lon_bnds: Optional[Tuple[float, float]] = None,
    ) -> xr.DataArray:
        if lat_bnds is not None:
            lat0, lat1 = lat_bnds
            da = da.sel(lat=slice(lat0, lat1))
        if lon_bnds is not None:
            lon0, lon1 = lon_bnds
            if lon0 <= lon1:
                da = da.sel(lon=slice(lon0, lon1))
            else:
                # dateline-crossing region
                da = xr.concat(
                    [da.sel(lon=slice(lon0, 180)), da.sel(lon=slice(-180, lon1))],
                    dim="lon",
                ).sortby("lon")
        return da

    # ------------------------------------------------------------------
    # Temporal aggregation
    # ------------------------------------------------------------------
    @staticmethod
    def compute_monthly_anomaly(da: xr.DataArray) -> xr.DataArray:
        clim = da.groupby("time.month").mean("time")
        return da.groupby("time.month") - clim

    @staticmethod
    def compute_time_aggregation(
        da: xr.DataArray,
        season: Optional[str],
        period: Optional[Tuple[int, int]] = None,
    ) -> xr.DataArray:
        if season is None:
            return da

        season = season.upper()

        if season == "MONTHLY":
            return da

        if season == "ANNUAL":
            return da.groupby("time.year").mean("time").rename({"year": "time"})

        if season == "DJF":
            winter_year = xr.where(
                da["time"].dt.month == 12,
                da["time"].dt.year + 1,
                da["time"].dt.year,
            )
            djf = da.where(da["time"].dt.month.isin([12, 1, 2]), drop=True)
            djf = djf.assign_coords(
                winter_year=("time", winter_year.sel(time=djf.time).data)
            )
            grouped = djf.groupby("winter_year")
            counts = grouped.count("time")
            out = grouped.mean("time").where(counts >= 3)
            out = out.dropna("winter_year", how="all").rename({"winter_year": "time"})
            if period is not None:
                y0, y1 = period
                out = out.sel(time=slice(y0 + 1, y1 + 1))
            return out

        if season in ("MAM", "JJA", "SON"):
            sub = da.where(da["time"].dt.season == season, drop=True)
            grouped = sub.groupby("time.year")
            counts = grouped.count("time")
            out = grouped.mean("time").where(counts >= 3)
            out = out.dropna("year", how="all").rename({"year": "time"})
            if period is not None:
                y0, y1 = period
                out = out.sel(time=slice(y0, y1))
            return out

        raise ValueError(f"Unsupported season: {season!r}")

    # ------------------------------------------------------------------
    # Weights
    # ------------------------------------------------------------------
    @staticmethod
    def compute_sqrtcoslat_weights(da: xr.DataArray) -> xr.DataArray:
        return np.sqrt(np.cos(np.deg2rad(da.lat)))

    @staticmethod
    def compute_coslat_weights_2d(da: xr.DataArray) -> xr.DataArray:
        w_lat = np.cos(np.deg2rad(da.lat))
        return w_lat.broadcast_like(da.isel(time=0))

    @classmethod
    def remove_mode_domain_mean(
        cls,
        da: xr.DataArray,
        mode_info: dict,
    ) -> xr.DataArray:
        """
        Subtract the area-weighted mode-domain mean at each time step.

        PCMDI Metrics applies this residual step after anomaly/seasonal
        aggregation and before EOF/regression diagnostics.
        """
        domain = cls.subset_latlon(
            da,
            lat_bnds=mode_info.get("lat_bnds"),
            lon_bnds=mode_info.get("lon_bnds"),
        )
        weights = cls.compute_coslat_weights_2d(domain)
        mask = np.isfinite(domain)
        weights = weights.where(mask)
        domain_mean = (domain.where(mask) * weights).sum(("lat", "lon")) / weights.sum(("lat", "lon"))
        return da - domain_mean

    @staticmethod
    def _ensure_local_pcmdi_on_path() -> None:
        pcmdi_root = BaseModeAnalysis.LOCAL_PCMDI_METRICS_ROOT
        if pcmdi_root and os.path.isdir(pcmdi_root) and pcmdi_root not in sys.path:
            sys.path.insert(0, pcmdi_root)

    @classmethod
    def preprocess_field_pcmdi(
        cls,
        da: xr.DataArray,
        mode_info: dict,
        period: Optional[Tuple[int, int]] = None,
        season_override: Optional[str] = None,
    ) -> Optional[xr.DataArray]:
        """Run PCMDI Metrics' own anomaly/season/residual preprocessing."""
        cls._ensure_local_pcmdi_on_path()
        try:
            adjust_mod = import_module("pcmdi_metrics.variability_mode.lib.adjust_timeseries")
            io_mod = import_module("pcmdi_metrics.io")
            adjust_timeseries = adjust_mod.adjust_timeseries
            load_regions_specs = io_mod.load_regions_specs
        except (ImportError, ModuleNotFoundError, AttributeError):
            return None

        var_name = da.name or mode_info.get("var", "field")
        mode_name = mode_info.get("_name", "")
        season = season_override if season_override is not None else mode_info.get("season", "monthly")
        parent_ds = da.attrs.get("_pcmdi_parent_dataset")
        if isinstance(parent_ds, xr.Dataset) and var_name in parent_ds:
            ds = parent_ds.copy()
        else:
            ds = da.to_dataset(name=var_name)

        if "lon" in ds.coords and float(ds.lon.min()) < 0.0:
            ds = ds.assign_coords(lon=(ds.lon % 360)).sortby("lon")

        if period is not None:
            y0, y1 = period
            ds = ds.sel(time=slice(f"{y0:04d}-01-01", f"{y1:04d}-12-31"))

        if "time" in ds.coords and "calendar" not in ds["time"].encoding:
            ds["time"].encoding["calendar"] = "standard"

        ds = cls._add_latlon_bounds(ds)
        try:
            ds = ds.bounds.add_missing_bounds()
        except Exception:
            pass

        try:
            regions_specs = load_regions_specs()
            out = adjust_timeseries(
                ds,
                var_name,
                mode_name,
                season,
                regions_specs,
                mode_info.get("remove_domain_mean", True),
            )
        except Exception as exc:
            print(f"  PCMDI preprocessing unavailable for this field ({exc}); using local preprocessing.")
            return None

        out_da = out[var_name]
        # Copy attributes from original DataArray to preserve units
        for k, v in da.attrs.items():
            if k not in out_da.attrs:
                out_da.attrs[k] = v
        out_da.attrs["_pcmdi_parent_dataset"] = out
        out_da.attrs["preprocessing"] = "pcmdi_metrics.adjust_timeseries"
        return out_da

    @classmethod
    def subset_mode_domain_pcmdi(
        cls,
        da: xr.DataArray,
        mode_info: dict,
    ) -> Optional[xr.DataArray]:
        """Subset a field to the EOF mode domain using PCMDI region specs."""
        cls._ensure_local_pcmdi_on_path()
        try:
            io_mod = import_module("pcmdi_metrics.io")
            load_regions_specs = io_mod.load_regions_specs
            region_subset = io_mod.region_subset
        except (ImportError, ModuleNotFoundError, AttributeError):
            return None

        mode_name = mode_info.get("_name", "")
        var_name = da.name or mode_info.get("var", "field")
        parent_ds = da.attrs.get("_pcmdi_parent_dataset")
        if isinstance(parent_ds, xr.Dataset) and var_name in parent_ds:
            ds = parent_ds.copy()
        else:
            ds = da.to_dataset(name=var_name)
        ds = cls._add_latlon_bounds(ds)
        try:
            out = region_subset(ds, mode_name, data_var=var_name, regions_specs=load_regions_specs())
        except Exception as exc:
            print(f"  PCMDI region_subset unavailable for this field ({exc}); using local subset.")
            return None
        out_da = out[var_name]
        # Copy attributes from original DataArray to preserve units
        for k, v in da.attrs.items():
            if k not in out_da.attrs:
                out_da.attrs[k] = v
        out_da.attrs["_pcmdi_parent_dataset"] = out
        return out_da

    @classmethod
    def subset_mode_domain(
        cls,
        da: xr.DataArray,
        mode_info: dict,
    ) -> xr.DataArray:
        out = cls.subset_mode_domain_pcmdi(da, mode_info)
        if out is not None:
            return out
        return cls.subset_latlon(
            da,
            lat_bnds=mode_info.get("lat_bnds"),
            lon_bnds=mode_info.get("lon_bnds"),
        )

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------
    @staticmethod
    def read_model_data(dataset_info: dict, var_name: str) -> xr.DataArray:
        pattern = os.path.join(
            dataset_info["dir"],
            dataset_info["name"],
            dataset_info["subdir"],
            f"{var_name}_Amon_*.nc",
        )
        files = sorted(glob.glob(pattern))
        if not files:
            raise FileNotFoundError(f"No model files found for pattern:\n{pattern}")
        return xr.open_mfdataset(files, combine="by_coords", data_vars="minimal", coords="minimal", compat="override")[var_name]

    @staticmethod
    def read_obs_data(obs_path: str, var_name: str) -> xr.DataArray:
        return xr.open_mfdataset(obs_path, combine="by_coords", data_vars="minimal", coords="minimal", compat="override")[var_name]

    @classmethod
    def read_pcmdi_data(
        cls,
        path: Sequence[str] | str,
        var_name: str,
        period: Tuple[int, int],
    ) -> Optional[xr.DataArray]:
        """Read monthly input through PCMDI Metrics' own xcdat reader."""
        cls._ensure_local_pcmdi_on_path()
        try:
            lib_mod = import_module("pcmdi_metrics.variability_mode.lib")
            read_data_in = lib_mod.read_data_in
        except (ImportError, ModuleNotFoundError, AttributeError):
            return None

        units_adjust = (True, "multiply", 0.01) if var_name == "psl" else None
        try:
            ds = read_data_in(
                path,
                var_name,
                var_name,
                period[0],
                period[1],
                UnitsAdjust=units_adjust,
                LandMask=False,
                debug=False,
            )
        except Exception as exc:
            print(f"  PCMDI read_data_in unavailable for this input ({exc}); using xarray reader.")
            return None

        da = ds[var_name]
        # Keep the exact PCMDI/xcdat Dataset available at runtime.  The
        # lat_bnds/lon_bnds variables cannot be attached directly to a 3-D
        # DataArray as coordinates because they carry a separate ``bnds`` dim.
        # Reusing the parent Dataset avoids rebuilding synthetic bounds before
        # PCMDI's area-weighted EOF calculation.
        da.attrs["_pcmdi_parent_dataset"] = ds
        if var_name == "psl":
            da.attrs["units"] = "hPa"
            da.attrs["pcmdi_units_adjustment"] = "pcmdi read_data_in UnitsAdjust x0.01"
        da.attrs["analysis_reader"] = "pcmdi_metrics.read_data_in"
        return da

    @staticmethod
    def apply_pcmdi_unit_adjustment(da: xr.DataArray, var_name: str) -> xr.DataArray:
        """
        Apply the unit convention used by the PCMDI variability-mode runs.

        The PCMDI diagnostics for sea-level pressure are written in hPa.  The
        upstream E3SM/obs monthly files commonly store ``psl`` in Pa, so convert
        those inputs before anomaly/EOF analysis.
        """
        if var_name != "psl":
            return da

        units = str(da.attrs.get("units", "")).lower()
        needs_pa_to_hpa = units in {"pa", "pascal", "pascals"} or not units
        if not needs_pa_to_hpa:
            return da

        out = da * 0.01
        out.attrs.update(da.attrs)
        out.attrs["units"] = "hPa"
        out.attrs["pcmdi_units_adjustment"] = "psl Pa to hPa (x0.01)"
        return out

    # ------------------------------------------------------------------
    # Pre-processing pipeline
    # ------------------------------------------------------------------
    def preprocess_field(
        self,
        da: xr.DataArray,
        mode_info: dict,
        period: Optional[Tuple[int, int]] = None,
        season_override: Optional[str] = None,
        subset_space: bool = True,
    ) -> xr.DataArray:
        pcmdi_da = self.preprocess_field_pcmdi(
            da,
            mode_info,
            period=period,
            season_override=season_override,
        )
        if pcmdi_da is not None:
            if subset_space:
                pcmdi_da = self.subset_mode_domain(pcmdi_da, mode_info)
            return pcmdi_da

        da = self.normalize_lon(da)
        season = season_override if season_override is not None else mode_info.get("season", "monthly")
        if period is not None and str(season).upper() == "DJF":
            y0, y1 = period
            da = da.sel(time=slice(f"{y0}-01-01", f"{y1 + 1}-02-28"))
        else:
            da = self.subset_time(da, period=period)
        if subset_space:
            da = self.subset_mode_domain(da, mode_info)
        da = self.compute_monthly_anomaly(da)
        da = self.compute_time_aggregation(da, season, period=period)
        if mode_info.get("remove_domain_mean", True):
            da = self.remove_mode_domain_mean(da, mode_info)
        return da

    # ------------------------------------------------------------------
    # EOF analysis
    # ------------------------------------------------------------------
    @staticmethod
    def _linear_regression_maps(
        da: xr.DataArray, pc: xr.DataArray
    ) -> Tuple[xr.DataArray, xr.DataArray]:
        """Return OLS slope and intercept maps for field = slope * pc + intercept."""
        da_mean = da.mean("time")
        pc_mean = pc.mean("time")
        da_anom = da - da_mean
        pc_anom = pc - pc_mean
        pc_var = (pc_anom ** 2).mean("time")
        slope = (da_anom * pc_anom).mean("time") / pc_var
        intercept = da_mean - slope * pc_mean
        return slope, intercept

    @staticmethod
    def _slope_corr_pval(
        da: xr.DataArray, pc_std: xr.DataArray
    ) -> Tuple[xr.DataArray, xr.DataArray, xr.DataArray, xr.DataArray]:
        """
        Compute the OLS regression slope of *da* onto *pc_std* at each grid point,
        the Pearson correlation, and the two-tailed p-value of that correlation.

        Returns
        -------
        slope : xr.DataArray  (lat, lon)
        corr  : xr.DataArray  (lat, lon)
        pval  : xr.DataArray  (lat, lon)  – two-tailed p-value via t-test
        """
        n = da.sizes["time"]

        # OLS slope/intercept for field = slope * standardized_pc + intercept.
        slope, intercept = BaseModeAnalysis._linear_regression_maps(da, pc_std)

        # Pearson r
        corr = xr.corr(da, pc_std, dim="time")

        # Two-tailed p-value: t = r*sqrt(n-2)/sqrt(1-r^2), df = n-2
        r_np   = corr.values
        t_stat = r_np * np.sqrt(n - 2) / np.sqrt(np.maximum(1 - r_np ** 2, 1e-15))
        pval   = xr.DataArray(
            2.0 * stats.t.sf(np.abs(t_stat), df=n - 2),
            dims=corr.dims,
            coords=corr.coords,
        )
        return slope, intercept, corr, pval

    @classmethod
    def _pcmdi_regression_maps(
        cls,
        pc_raw: xr.DataArray,
        reg_field: xr.DataArray,
        *,
        var_name: Optional[str] = None,
        rm_domain_mean: bool = True,
        eof_scaling: bool = False,
    ) -> Optional[Tuple[xr.DataArray, xr.DataArray, xr.DataArray]]:
        """Use PCMDI Metrics' regression/eof reconstruction routine exactly."""
        cls._ensure_local_pcmdi_on_path()
        try:
            eof_mod = import_module("pcmdi_metrics.variability_mode.lib.eof_analysis")
            stat_mod = import_module("pcmdi_metrics.variability_mode.lib.calc_stat")
            linear_regression_on_globe_for_teleconnection = (
                eof_mod.linear_regression_on_globe_for_teleconnection
            )
            calcSTD = stat_mod.calcSTD
        except (ImportError, ModuleNotFoundError, AttributeError):
            return None

        data_var = var_name or reg_field.name or "field"
        parent_ds = reg_field.attrs.get("_pcmdi_parent_dataset")
        if isinstance(parent_ds, xr.Dataset) and data_var in parent_ds:
            reg_ds = parent_ds.copy()
        else:
            reg_ds = reg_field.to_dataset(name=data_var)
        reg_ds = cls._add_latlon_bounds(reg_ds)
        stdv_pc = calcSTD(pc_raw)
        eof_lr, slope, intercept = linear_regression_on_globe_for_teleconnection(
            pc_raw,
            reg_ds,
            data_var,
            stdv_pc,
            rm_domain_mean,
            eof_scaling,
            debug=False,
        )
        return eof_lr, slope, intercept

    @classmethod
    def _pcmdi_regrid_to_reference(
        cls,
        ds: xr.Dataset,
        data_var: str,
        ref_grid_global: Optional[xr.DataArray],
    ) -> xr.Dataset:
        """Mirror PCMDI CBF regridding: model adjusted field -> obs global grid."""
        if ref_grid_global is None:
            return ds

        cls._ensure_local_pcmdi_on_path()
        try:
            utils_mod = import_module("pcmdi_metrics.utils")
            regrid = utils_mod.regrid
        except (ImportError, ModuleNotFoundError, AttributeError):
            return ds

        target = cls._add_latlon_bounds(ref_grid_global.to_dataset(name=data_var))
        try:
            return regrid(ds, data_var, target, regrid_tool="regrid2", fill_zero=True)
        except Exception as exc:
            print(f"  PCMDI reference-grid regrid unavailable ({exc}); using native projection grid.")
            return ds

    def project_onto_eof(
        self, da: xr.DataArray, ref_eof: xr.DataArray,
        da_global: Optional[xr.DataArray] = None,
        ref_solver: Optional[object] = None,
        ref_reverse_sign: bool = False,
        ref_grid_global: Optional[xr.DataArray] = None,
        eof_num: int = 1,
        mode_name: str = "",
        mode_info: Optional[dict] = None,
    ) -> xr.Dataset:
        """
        Project *da* onto an external EOF pattern (common-base approach).

        The reference EOF *ref_eof* (typically from observations) defines the
        spatial basis.  The field *da* is projected onto it via a cosine-latitude-
        weighted inner product to produce a standardised PC, from which regression
        diagnostics are derived.  Results are stored with the ``_proj`` suffix so
        they can coexist with the independent-EOF variables in the same Dataset.

        Parameters
        ----------
        da : xr.DataArray  (time, lat, lon)
            Pre-processed anomaly field (regional, used only for PC computation).
        ref_eof : xr.DataArray  (lat, lon)
            Reference EOF spatial pattern in physical units (e.g. obs EOF1).
        da_global : xr.DataArray, optional
            Full-domain anomaly field (time, lat, lon) used for regression /
            correlation diagnostics.  If None, *da* is used instead.

        Returns
        -------
        xr.Dataset
            Variables: ``pc_proj``, ``slope_proj``, ``slope_pval_proj``, ``corr_proj``.
        """
        coslat  = np.cos(np.deg2rad(da.lat))
        pc_raw  = (da * coslat * ref_eof).sum(("lat", "lon"))
        reg_field = da_global if da_global is not None else da
        common_base_method = "local_weighted_dot_product"
        frac_proj = xr.DataArray(np.nan, name="frac_proj")

        if ref_solver is not None:
            self._ensure_local_pcmdi_on_path()
            try:
                eof_mod = import_module("pcmdi_metrics.variability_mode.lib.eof_analysis")
                gain_pseudo_pcs = eof_mod.gain_pseudo_pcs
                gain_pcs_fraction = eof_mod.gain_pcs_fraction
                linear_regression_on_globe_for_teleconnection = eof_mod.linear_regression_on_globe_for_teleconnection

                var_name = da.name or "field"
                reg_var = reg_field.name or var_name
                reg_ds = self._add_latlon_bounds(reg_field.to_dataset(name=reg_var))
                regrid_ds = self._pcmdi_regrid_to_reference(reg_ds, reg_var, ref_grid_global)
                regrid_da = regrid_ds[reg_var]
                if mode_info is not None:
                    da_for_projection = self.subset_mode_domain(regrid_da, mode_info)
                else:
                    da_for_projection = self.subset_latlon(
                        regrid_da,
                        lat_bnds=(float(da.lat.min()), float(da.lat.max())),
                        lon_bnds=(float(da.lon.min()), float(da.lon.max())),
                    )
                pc_raw = gain_pseudo_pcs(
                    ref_solver,
                    da_for_projection,
                    eofn=eof_num,
                    reverse_sign=ref_reverse_sign,
                    EofScaling=False,
                ).rename("pc_proj_raw")
                pc_raw = pc_raw.assign_coords(time=da.time)
                stdv_pc = float(pc_raw.std("time"))
                pc_proj = ((pc_raw - pc_raw.mean("time")) / pc_raw.std("time")).rename("pc_proj")

                eof_lr_proj, slope_proj, intercept_proj = linear_regression_on_globe_for_teleconnection(
                    pc_raw,
                    reg_ds,
                    reg_var,
                    stdv_pc=stdv_pc,
                    RmDomainMean=True,
                    EofScaling=False,
                    debug=False,
                )
                eof_lr_proj = eof_lr_proj.rename("eof_lr_proj")
                slope_proj = slope_proj.rename("slope_proj")
                intercept_proj = intercept_proj.rename("intercept_proj")
                corr_proj = xr.corr(reg_field, pc_proj, dim="time")

                r_np = corr_proj.values
                n = reg_field.sizes["time"]
                t_stat = r_np * np.sqrt(n - 2) / np.sqrt(np.maximum(1 - r_np ** 2, 1e-15))
                pval_proj = xr.DataArray(
                    2.0 * stats.t.sf(np.abs(t_stat), df=n - 2),
                    dims=corr_proj.dims,
                    coords=corr_proj.coords,
                )

                if mode_info is not None:
                    frac_field = self.subset_mode_domain(reg_field, mode_info)
                    eof_lr_region = self.subset_mode_domain(eof_lr_proj, mode_info)
                else:
                    frac_field = da
                    eof_lr_region = self.subset_latlon(
                        eof_lr_proj,
                        lat_bnds=(float(da.lat.min()), float(da.lat.max())),
                        lon_bnds=(float(da.lon.min()), float(da.lon.max())),
                    )
                frac_val = gain_pcs_fraction(
                    self._add_latlon_bounds(frac_field.to_dataset(name=reg_var)),
                    reg_var,
                    eof_lr_region.rename("eof_lr_proj").to_dataset(name="eof_lr_proj"),
                    "eof_lr_proj",
                    pc_raw / stdv_pc,
                    debug=False,
                )
                frac_proj = xr.DataArray(frac_val, name="frac_proj")
                common_base_method = "pcmdi_metrics.gain_pseudo_pcs"
            except Exception as exc:
                print(f"  PCMDI common-base projection unavailable ({exc}); using local projection.")
                pc_proj = (pc_raw - pc_raw.mean("time")) / pc_raw.std("time")
                pc_raw = pc_raw.rename("pc_proj")
                slope_proj, intercept_proj, corr_proj, pval_proj = self._slope_corr_pval(reg_field, pc_proj)
                eof_lr_proj = slope_proj + intercept_proj
        else:
            pc_proj = (pc_raw - pc_raw.mean("time")) / pc_raw.std("time")
            pc_raw = pc_raw.rename("pc_proj")
            slope_proj, intercept_proj, corr_proj, pval_proj = self._slope_corr_pval(reg_field, pc_proj)
            eof_lr_proj = slope_proj + intercept_proj

        out = xr.Dataset({
            "pc_proj":         pc_raw.rename("pc_proj"),
            "pc_proj_std":     pc_proj.rename("pc_proj_std"),
            "eof_lr_proj":     eof_lr_proj,
            "frac_proj":       frac_proj,
            "slope_proj":      slope_proj,
            "intercept_proj":  intercept_proj,
            "slope_pval_proj": pval_proj,
            "corr_proj":       corr_proj,
        })
        out["pc_proj"        ].attrs["long_name"] = "Raw PC from projection onto obs EOF"
        out["pc_proj_std"    ].attrs["long_name"] = "Standardized PC from projection onto obs EOF"
        out["eof_lr_proj"    ].attrs["long_name"] = "Reconstructed common-base EOF map from obs-projected PC"
        out["slope_proj"     ].attrs["long_name"] = "Regression slope onto raw obs-projected PC"
        out["intercept_proj" ].attrs["long_name"] = "Regression intercept for obs-projected PC"
        out["slope_pval_proj"].attrs["long_name"] = "Two-tailed p-value of obs-projected regression slope"
        out["slope_pval_proj"].attrs["description"] = "Derived from t-test on Pearson r; df = n_time - 2"
        out["corr_proj"      ].attrs["long_name"] = "Pearson r with obs-projected standardized PC"
        out["frac_proj"      ].attrs["long_name"] = "Explained variance fraction from common-base pseudo-PC"
        out.attrs["common_base_projection_method"] = common_base_method
        out.attrs["pc_proj_saved_as"] = "raw_pcmdi_pseudo_pc"
        out.attrs["common_base_regrid"] = (
            "pcmdi_metrics.regrid_to_reference_grid"
            if ref_grid_global is not None else "native_grid"
        )
        return out

    @staticmethod
    def _weighted_spatial_corr(a: xr.DataArray, b: xr.DataArray) -> float:
        """Area-weighted spatial correlation over finite overlapping grid cells."""
        a, b = xr.align(a, b, join="inner")
        if "lat" not in a.coords:
            return np.nan

        mask = np.isfinite(a) & np.isfinite(b)
        if int(mask.sum()) < 2:
            return np.nan

        weights = np.cos(np.deg2rad(a["lat"]))
        weights = xr.broadcast(weights, a)[0].where(mask)
        a = a.where(mask)
        b = b.where(mask)

        wsum = weights.sum()
        if not np.isfinite(float(wsum)) or float(wsum) == 0.0:
            return np.nan

        amean = (a * weights).sum() / wsum
        bmean = (b * weights).sum() / wsum
        aanom = a - amean
        banom = b - bmean
        denom = np.sqrt(float((weights * aanom ** 2).sum()) * float((weights * banom ** 2).sum()))
        if denom == 0.0 or not np.isfinite(denom):
            return np.nan
        return float((weights * aanom * banom).sum()) / denom

    @classmethod
    def orient_eof_to_reference(cls, ds: xr.Dataset, ref_eof: xr.DataArray) -> xr.Dataset:
        """
        Flip independent EOF diagnostics so model EOFs have positive spatial
        correlation with the observation EOF.

        EOF signs are arbitrary.  This keeps model pattern panels comparable to
        the obs panel while preserving variance fraction and p-values.
        """
        if "eof" not in ds:
            return ds

        corr = cls._weighted_spatial_corr(ds["eof"], ref_eof)
        sign = -1.0 if np.isfinite(corr) and corr < 0.0 else 1.0
        if sign < 0.0:
            ds = ds.copy()
            for name in ("eof", "pc", "slope", "corr"):
                if name in ds:
                    attrs = ds[name].attrs.copy()
                    ds[name] = -ds[name]
                    ds[name].attrs.update(attrs)
            if "eof_lr" in ds and "slope" in ds and "intercept" in ds:
                attrs = ds["eof_lr"].attrs.copy()
                ds["eof_lr"] = ds["slope"] + ds["intercept"]
                ds["eof_lr"].attrs.update(attrs)
        ds.attrs["eof_oriented_to_ref"] = "true"
        ds.attrs["eof_ref_spatial_corr_before_orientation"] = corr
        return ds

    def run_eof(
        self,
        da: xr.DataArray,
        eof_num: int = 1,
        da_global: Optional[xr.DataArray] = None,
    ) -> xr.Dataset:
        """
        Parameters
        ----------
        da : xr.DataArray  (time, lat, lon)
            Regional anomaly field used for EOF computation.
        eof_num : int
            Which EOF to extract (1-based).
        da_global : xr.DataArray, optional
            Full-domain anomaly field used for regression / correlation
            diagnostics (``slope``, ``slope_pval``, ``corr``).  If None,
            the regional *da* is used — slope will only cover the mode domain.
        """
        weights = self.compute_sqrtcoslat_weights(da)
        pcmdi_result = self._run_eof_pcmdi(
            da,
            eof_num=eof_num,
            mode_name=da.attrs.get("mode_name", ""),
        )
        solver = None
        if pcmdi_result is not None:
            eof, pc, frac, _reverse_sign, solver = pcmdi_result
        elif Eof is not None:
            solver = Eof(da * weights)

            eof_w = solver.eofs(neofs=eof_num)[eof_num - 1]
            eof = eof_w / weights
            pc = solver.pcs(npcs=eof_num, pcscaling=1)[:, eof_num - 1]
            frac = solver.varianceFraction(neigs=eof_num)[eof_num - 1]
        else:
            eof, pc, frac = self._run_eof_svd(da, weights, eof_num=eof_num)

        pc_raw = pc.rename("pc")
        pc_std = ((pc_raw - pc_raw.mean("time")) / pc_raw.std("time")).rename("pc_std")
        reg_field = da_global if da_global is not None else da
        corr = xr.corr(reg_field, pc_std, dim="time")
        r_np = corr.values
        n = reg_field.sizes["time"]
        t_stat = r_np * np.sqrt(n - 2) / np.sqrt(np.maximum(1 - r_np ** 2, 1e-15))
        pval = xr.DataArray(
            2.0 * stats.t.sf(np.abs(t_stat), df=n - 2),
            dims=corr.dims,
            coords=corr.coords,
        )

        pcmdi_reg = self._pcmdi_regression_maps(
            pc_raw,
            reg_field,
            var_name=reg_field.name or da.name or "field",
        )
        if pcmdi_reg is not None:
            eof_lr, slope, intercept = pcmdi_reg
        else:
            slope, intercept, _corr_unused, _pval_unused = self._slope_corr_pval(reg_field, pc_std)
            eof_lr = slope + intercept

        # PCMDI's saved obs/model "eof" pattern is the reconstructed regression
        # map: slope(original PC) * stdv_pc + intercept(original PC).  We use
        # PCMDI's own regression routine above when available.

        # Expand the regional EOF to the global grid (NaN outside the mode
        # domain) so that all variables share the same lat/lon coordinates and
        # can be stored together in a single NetCDF file without coordinate
        # conflicts between the regional eof and the global slope/corr maps.
        if da_global is not None:
            ref_spatial = da_global.isel(time=0, drop=True)
            eof = eof.reindex_like(ref_spatial, fill_value=np.nan)

        out = xr.Dataset({
            "pc":         pc_raw,
            "pc_std":     pc_std,
            "eof":        eof,
            "eof_lr":     eof_lr,
            "slope":      slope,
            "intercept":  intercept,
            "slope_pval": pval,
            "corr":       corr,
            "frac":       frac,
        })
        out["pc"        ].attrs["long_name"] = f"Raw PC of EOF{eof_num}"
        out["pc_std"    ].attrs["long_name"] = f"Standardized PC of EOF{eof_num}"
        out["eof"       ].attrs["long_name"] = f"EOF{eof_num} spatial pattern"
        out["eof_lr"    ].attrs["long_name"] = f"Reconstructed EOF{eof_num} regression map"
        out["slope"     ].attrs["long_name"] = f"Regression slope onto raw EOF{eof_num} PC"
        out["intercept" ].attrs["long_name"] = f"Regression intercept for EOF{eof_num} PC"
        out["slope_pval"].attrs["long_name"] = f"Two-tailed p-value of regression slope (EOF{eof_num})"
        out["slope_pval"].attrs["description"] = "Derived from t-test on Pearson r; df = n_time - 2"
        out["corr"      ].attrs["long_name"] = f"Pearson r with standardized EOF{eof_num} PC"
        out["frac"      ].attrs["long_name"] = f"Explained variance fraction of EOF{eof_num}"
        for attr_name in ("eof_solver", "pcmdi_metrics_module"):
            if attr_name in eof.attrs:
                out.attrs[attr_name] = eof.attrs[attr_name]
        out.attrs["pcmdi_regression_method"] = (
            "pcmdi_metrics.linear_regression_on_globe_for_teleconnection"
            if pcmdi_reg is not None else "local_standardized_pc_regression"
        )
        out.attrs["pc_saved_as"] = "raw_pcmdi_pc"
        if "preprocessing" in da.attrs:
            out.attrs["preprocessing"] = da.attrs["preprocessing"]
        if solver is not None:
            out.attrs["_pcmdi_solver_available_runtime_only"] = "true"
            self._last_pcmdi_solver = solver
            self._last_pcmdi_reverse_sign = bool(_reverse_sign)
        return out

    @staticmethod
    def _run_eof_pcmdi(
        da: xr.DataArray,
        eof_num: int = 1,
        mode_name: str = "",
    ) -> Optional[Tuple[xr.DataArray, xr.DataArray, xr.DataArray, bool, object]]:
        """Run EOF analysis with PCMDI Metrics when available."""
        pcmdi_root = BaseModeAnalysis.LOCAL_PCMDI_METRICS_ROOT
        if pcmdi_root and os.path.isdir(pcmdi_root) and pcmdi_root not in sys.path:
            sys.path.insert(0, pcmdi_root)

        try:
            eof_mod = import_module("pcmdi_metrics.variability_mode.lib.eof_analysis")
            eof_func = eof_mod.eof_analysis_get_variance_mode
        except (ImportError, ModuleNotFoundError, AttributeError):
            return None

        var_name = da.name or "field"
        parent_ds = da.attrs.get("_pcmdi_parent_dataset")
        if isinstance(parent_ds, xr.Dataset) and var_name in parent_ds:
            ds = parent_ds.copy()
        else:
            ds = da.to_dataset(name=var_name)
        ds = BaseModeAnalysis._add_latlon_bounds(ds)
        try:
            eof, pc, frac, reverse_sign, solver = eof_func(
                mode=mode_name,
                ds=ds,
                data_var=var_name,
                eofn=eof_num,
                eofn_max=eof_num,
                debug=False,
                EofScaling=False,
                save_multiple_eofs=False,
            )
        except Exception as exc:
            print(f"  PCMDI EOF analysis unavailable for this field ({exc}); using fallback EOF solver.")
            return None

        eof = eof.rename("eof")
        pc = pc.rename("pc")
        frac = frac.rename("frac")
        eof.attrs["eof_solver"] = "pcmdi_metrics.eof_analysis_get_variance_mode"
        eof.attrs["pcmdi_metrics_module"] = getattr(eof_mod, "__file__", "")
        return eof, pc, frac, reverse_sign, solver

    @staticmethod
    def _add_latlon_bounds(ds: xr.Dataset) -> xr.Dataset:
        """Add simple CF-style lat/lon bounds required by PCMDI grid-area utilities."""
        ds = ds.copy()

        def _bounds_1d(coord: xr.DataArray, *, clip: Optional[Tuple[float, float]] = None) -> xr.DataArray:
            values = np.asarray(coord.values, dtype=float)
            if values.ndim != 1 or values.size < 2:
                raise ValueError(f"Cannot build bounds for coordinate {coord.name!r}.")
            mids = 0.5 * (values[:-1] + values[1:])
            first = values[0] - 0.5 * (values[1] - values[0])
            last = values[-1] + 0.5 * (values[-1] - values[-2])
            edges = np.concatenate([[first], mids, [last]])
            if clip is not None:
                edges = np.clip(edges, clip[0], clip[1])
            bnds = np.column_stack([edges[:-1], edges[1:]])
            return xr.DataArray(
                bnds,
                dims=(coord.name, "bnds"),
                coords={coord.name: coord, "bnds": [0, 1]},
            )

        if "lat" in ds.coords and "lat_bnds" not in ds:
            ds["lat_bnds"] = _bounds_1d(ds["lat"], clip=(-90.0, 90.0))
            ds["lat"].attrs["bounds"] = "lat_bnds"
        if "lon" in ds.coords and "lon_bnds" not in ds:
            ds["lon_bnds"] = _bounds_1d(ds["lon"])
            ds["lon"].attrs["bounds"] = "lon_bnds"
        return ds

    @staticmethod
    def _run_eof_svd(
        da: xr.DataArray,
        weights: xr.DataArray,
        eof_num: int = 1,
    ) -> Tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
        """
        Compute EOFs using NumPy SVD when the optional ``eofs`` package is absent.

        The input field is weighted by sqrt(cos(lat)), centered in time, reshaped
        to time x space, and decomposed with SVD.  The returned EOF is converted
        back to physical units by dividing by the same latitude weights.
        """
        weighted = (da * weights).transpose("time", "lat", "lon")
        centered = weighted - weighted.mean("time")

        arr = np.asarray(centered.values)
        n_time = arr.shape[0]
        flat = arr.reshape(n_time, -1)
        valid = np.isfinite(flat).all(axis=0)
        if valid.sum() == 0:
            raise ValueError("No finite grid points available for EOF analysis.")
        if eof_num > min(n_time, valid.sum()):
            raise ValueError(
                f"Requested EOF{eof_num}, but only {min(n_time, valid.sum())} modes are available."
            )

        u, s, vt = np.linalg.svd(flat[:, valid], full_matrices=False)
        mode_idx = eof_num - 1

        eof_flat = np.full(flat.shape[1], np.nan, dtype=float)
        eof_flat[valid] = vt[mode_idx]
        eof_w = xr.DataArray(
            eof_flat.reshape(weighted.sizes["lat"], weighted.sizes["lon"]),
            dims=("lat", "lon"),
            coords={"lat": weighted.lat, "lon": weighted.lon},
        )
        eof = eof_w / weights

        pc = xr.DataArray(
            u[:, mode_idx] * s[mode_idx],
            dims=("time",),
            coords={"time": weighted.time},
        )
        eigvals = (s ** 2) / max(n_time - 1, 1)
        frac = xr.DataArray(eigvals[mode_idx] / eigvals.sum())
        return eof, pc, frac

    # ------------------------------------------------------------------
    # Area-mean index
    # ------------------------------------------------------------------
    def run_index(
        self,
        da: xr.DataArray,
        standardize: bool = False,
        running_mean: Optional[int] = None,
    ) -> xr.Dataset:
        w2d = self.compute_coslat_weights_2d(da)
        idx = da.weighted(w2d).mean(("lat", "lon"))

        if running_mean is not None and running_mean > 1:
            idx = idx.rolling(time=running_mean, center=True).mean()
        if standardize:
            idx = (idx - idx.mean("time")) / idx.std("time")

        out = xr.Dataset({"index": idx})
        out["index"].attrs["long_name"] = "Area-mean anomaly index"
        return out


class ModeAnalyzer(BaseModeAnalysis):
    """
    High-level analysis driver.

    Parameters
    ----------
    config_path : str or dict
        Path to the JSON configuration file or an in-memory dict with
        ``datasets`` and ``modes`` entries.
    """

    def __init__(self, config_path: str):
        super().__init__()
        self.cfg = ModeConfigManager(config_path)
        self._last_pcmdi_solver = None
        self._last_pcmdi_reverse_sign = False
        self._last_ref_grid_global = None
        self._last_runtime_key = None

    # ------------------------------------------------------------------
    # Observation helpers
    # ------------------------------------------------------------------
    def _resolve_observation(
        self,
        mode_name: str,
        custom_obs_path: Optional[str] = None,
        custom_obs_name: Optional[str] = None,
    ) -> dict:
        if custom_obs_path is not None:
            return {"name": custom_obs_name or "custom_obs", "data": custom_obs_path}
        obs = self.cfg.get_default_obs(mode_name)
        if obs is None:
            raise ValueError(f"No default observation configured for mode '{mode_name}'.")
        return obs

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def analyze_observation(
        self,
        mode_name: str,
        custom_obs_path: Optional[str] = None,
        custom_obs_name: Optional[str] = None,
        period: Optional[Tuple[int, int]] = None,
        season: Optional[str] = None,
    ) -> xr.Dataset:
        """Load, pre-process and analyse an observation dataset."""
        mode_info  = dict(self.cfg.get_mode_info(mode_name))
        mode_info["_name"] = mode_name
        obs_info   = self._resolve_observation(mode_name, custom_obs_path, custom_obs_name)
        raw        = self.read_pcmdi_data(obs_info["data"], mode_info["var"], period) if period is not None else None
        if raw is None:
            raw    = self.read_obs_data(obs_info["data"], mode_info["var"])
            raw    = self.apply_pcmdi_unit_adjustment(raw, mode_info["var"])
        # Pre-process globally first (no spatial subsetting), then subset.
        # This guarantees da_global is the full-domain field and avoids calling
        # preprocess_field twice on the same lazy DataArray (which can silently
        # produce a regional result the second time due to lazy-graph sharing).
        da_global  = self.preprocess_field(raw, mode_info, period=period,
                                           season_override=season, subset_space=False)
        self._last_ref_grid_global = da_global
        self._last_runtime_key = (mode_name, period, season, obs_info.get("name"))
        obs_pre    = self.subset_mode_domain(da_global, mode_info)
        return self._analyze_preprocessed(obs_pre, mode_info, da_global=da_global)

    def compute_observation_pcmdi_solver(
        self,
        mode_name: str,
        custom_obs_path: Optional[str] = None,
        custom_obs_name: Optional[str] = None,
        period: Optional[Tuple[int, int]] = None,
        season: Optional[str] = None,
    ) -> Tuple[Optional[object], bool]:
        """Rebuild the observational PCMDI EOF solver for common-base projection."""
        mode_info = dict(self.cfg.get_mode_info(mode_name))
        if mode_info.get("type", "eof") != "eof":
            return None, False

        mode_info["_name"] = mode_name
        obs_info = self._resolve_observation(mode_name, custom_obs_path, custom_obs_name)
        raw = self.read_pcmdi_data(obs_info["data"], mode_info["var"], period) if period is not None else None
        if raw is None:
            raw = self.read_obs_data(obs_info["data"], mode_info["var"])
            raw = self.apply_pcmdi_unit_adjustment(raw, mode_info["var"])
        da_global = self.preprocess_field(
            raw,
            mode_info,
            period=period,
            season_override=season,
            subset_space=False,
        )
        obs_pre = self.subset_mode_domain(da_global, mode_info)
        result = self._run_eof_pcmdi(
            obs_pre,
            eof_num=mode_info["eof_num"],
            mode_name=mode_name,
        )
        if result is None:
            return None, False

        _eof, _pc, _frac, reverse_sign, solver = result
        return solver, bool(reverse_sign)

    def analyze_model(
        self,
        mode_name: str,
        case_name: str,
        period: Optional[Tuple[int, int]] = None,
        season: Optional[str] = None,
        ref_eof: Optional[xr.DataArray] = None,
        ref_solver: Optional[object] = None,
        ref_reverse_sign: bool = False,
        ref_grid_global: Optional[xr.DataArray] = None,
    ) -> xr.Dataset:
        """Load, pre-process and analyse a single model case.

        Parameters
        ----------
        ref_eof : xr.DataArray, optional
            If provided (typically ``obs_ds["eof"]``), the field is also
            projected onto this common basis and the results are stored as
            ``slope_proj``, ``slope_pval_proj``, ``corr_proj``, ``pc_proj``.
        """
        mode_info    = dict(self.cfg.get_mode_info(mode_name))
        mode_info["_name"] = mode_name
        dataset_info = self.cfg.get_dataset_info(case_name)
        pattern = os.path.join(
            dataset_info["dir"],
            dataset_info["name"],
            dataset_info["subdir"],
            f"{mode_info['var']}_Amon_*.nc",
        )
        files = sorted(glob.glob(pattern))
        raw = self.read_pcmdi_data(files, mode_info["var"], period) if period is not None and files else None
        if raw is None:
            raw = self.read_model_data(dataset_info, mode_info["var"])
            raw = self.apply_pcmdi_unit_adjustment(raw, mode_info["var"])
        raw          = self.apply_pcmdi_unit_adjustment(raw, mode_info["var"])
        # Pre-process globally first (no spatial subsetting), then subset.
        # This guarantees da_global is the full-domain field and avoids calling
        # preprocess_field twice on the same lazy DataArray (which can silently
        # produce a regional result the second time due to lazy-graph sharing).
        da_global    = self.preprocess_field(raw, mode_info, period=period,
                                             season_override=season, subset_space=False)
        model_pre    = self.subset_mode_domain(da_global, mode_info)
        return self._analyze_preprocessed(
            model_pre,
            mode_info,
            ref_eof=ref_eof,
            ref_solver=ref_solver,
            ref_reverse_sign=ref_reverse_sign,
            ref_grid_global=ref_grid_global,
            da_global=da_global,
        )

    def analyze_all_models(
        self,
        mode_name: str,
        model_list: Optional[List[str]] = None,
        period: Optional[Tuple[int, int]] = None,
        season: Optional[str] = None,
    ) -> Dict[str, xr.Dataset]:
        """Analyse all models (or a specified subset) for a given mode."""
        if model_list is None:
            model_list = list(self.cfg.datasets.keys())
        results = {}
        for case_name in model_list:
            print(f"  Processing: {case_name}")
            results[case_name] = self.analyze_model(
                mode_name=mode_name,
                case_name=case_name,
                period=period,
                season=season,
            )
        return results

    # ------------------------------------------------------------------
    # Save / load helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _nc_path(out_dir: str, tag: str, prefix: str = "") -> str:
        """Return the NetCDF path.

        The filename is ``{prefix}{tag}.nc`` where *prefix* encodes
        mode / season / period / obs so that files are unique even when
        stored in a shared flat directory.
        """
        return os.path.join(out_dir, f"{prefix}{tag}.nc")

    @staticmethod
    def save_result(ds: xr.Dataset, path: str) -> None:
        """Write an analysis Dataset to NetCDF, creating parent dirs as needed."""
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        # Remove existing file first so xarray can write without a lock conflict.
        if os.path.exists(path):
            os.remove(path)
        ds = ds.copy()
        ds.attrs.pop("_pcmdi_parent_dataset", None)
        for name in list(ds.variables):
            ds[name].attrs.pop("_pcmdi_parent_dataset", None)
        ds.to_netcdf(path)
        print(f"  Saved → {path}")

    @staticmethod
    def load_result(path: str) -> xr.Dataset:
        """Read an analysis Dataset from NetCDF into memory (closes file handle)."""
        with xr.open_dataset(path) as ds:
            return ds.load()
        return xr.open_dataset(path)

    def _requires_saved_anomaly(self, mode_name: str, ds: xr.Dataset) -> bool:
        """Return True when a cached result is missing the saved anomaly field."""
        mode_info = self.cfg.get_mode_info(mode_name)
        var_name = mode_info.get("var")
        if var_name != "psl":
            return False
        return f"{var_name}_anom" not in ds.data_vars

    def _requires_local_pcmdi_eof(self, mode_name: str, ds: xr.Dataset) -> bool:
        """Return True when an EOF cache was not made with the local PCMDI helper."""
        mode_info = self.cfg.get_mode_info(mode_name)
        if mode_info.get("type", "eof") != "eof":
            return False

        expected = os.path.realpath(os.path.join(
            self.LOCAL_PCMDI_METRICS_ROOT,
            "pcmdi_metrics",
            "variability_mode",
            "lib",
            "eof_analysis.py",
        ))
        used = os.path.realpath(ds.attrs.get("pcmdi_metrics_module", ""))
        return used != expected

    # ------------------------------------------------------------------
    # analyse-or-load wrappers
    # ------------------------------------------------------------------
    def analyze_or_load_observation(
        self,
        mode_name: str,
        out_dir: str,
        custom_obs_path: Optional[str] = None,
        custom_obs_name: Optional[str] = None,
        period: Optional[Tuple[int, int]] = None,
        season: Optional[str] = None,
        file_prefix: str = "",
        overwrite: bool = False,
    ) -> xr.Dataset:
        """
        Load obs result from *out_dir* if the file already exists,
        otherwise run the analysis and save to *out_dir*.

        Parameters
        ----------
        file_prefix : str
            Prefix prepended to the filename, e.g.
            ``"NAO_DJF_1985-2014_HadISST2_"``.
        overwrite : bool
            If True, ignore any cached file and recompute from scratch.
        """
        obs_info = self._resolve_observation(mode_name, custom_obs_path, custom_obs_name)
        tag  = obs_info.get("name", "obs")
        path = self._nc_path(out_dir, tag, prefix=file_prefix)

        if os.path.exists(path) and not overwrite:
            ds = self.load_result(path)
            needs_anom = self._requires_saved_anomaly(mode_name, ds)
            needs_pcmdi = self._requires_local_pcmdi_eof(mode_name, ds)
            needs_eof_lr = self.cfg.get_mode_info(mode_name).get("type", "eof") == "eof" and not {"eof_lr", "intercept"}.issubset(ds.data_vars)
            needs_raw_pc = self.cfg.get_mode_info(mode_name).get("type", "eof") == "eof" and ds.attrs.get("pc_saved_as") != "raw_pcmdi_pc"
            needs_pcmdi_reg = (
                self.cfg.get_mode_info(mode_name).get("type", "eof") == "eof"
                and ds.attrs.get("pcmdi_regression_method")
                != "pcmdi_metrics.linear_regression_on_globe_for_teleconnection"
            )
            needs_units = (
                self.cfg.get_mode_info(mode_name).get("var") == "psl"
                and ds.attrs.get("analysis_units") != "hPa"
            )
            needs_reader = ds.attrs.get("analysis_reader") != "pcmdi_metrics.read_data_in"
            if (
                not needs_anom and not needs_pcmdi and not needs_eof_lr
                and not needs_raw_pc and not needs_pcmdi_reg and not needs_units
                and not needs_reader
            ):
                print(f"  Loading obs from cache: {path}")
                return ds
            missing = []
            if needs_anom:
                missing.append("psl_anom")
            if needs_pcmdi:
                missing.append("local PCMDI EOF stamp")
            if needs_eof_lr:
                missing.append("eof_lr/intercept variables")
            if needs_raw_pc:
                missing.append("raw PCMDI pc storage")
            if needs_pcmdi_reg:
                missing.append("PCMDI regression method")
            if needs_units:
                missing.append("PCMDI hPa unit normalization")
            if needs_reader:
                missing.append("PCMDI read_data_in reader")
            print(f"  Cache for obs ({tag}) missing {', '.join(missing)}, recomputing …")
        if os.path.exists(path) and overwrite:
            print(f"  Overwrite=True — recomputing obs ({tag}) …")

        print(f"  Computing obs ({tag}) …")
        ds = self.analyze_observation(
            mode_name=mode_name,
            custom_obs_path=custom_obs_path,
            custom_obs_name=custom_obs_name,
            period=period,
            season=season,
        )
        self.save_result(ds, path)
        return ds

    def analyze_or_load_model(
        self,
        mode_name: str,
        case_name: str,
        out_dir: str,
        period: Optional[Tuple[int, int]] = None,
        season: Optional[str] = None,
        ref_eof: Optional[xr.DataArray] = None,
        ref_solver: Optional[object] = None,
        ref_reverse_sign: bool = False,
        ref_grid_global: Optional[xr.DataArray] = None,
        file_prefix: str = "",
        overwrite: bool = False,
    ) -> xr.Dataset:
        """
        Load model result from *out_dir* if the file already exists,
        otherwise run the analysis and save to *out_dir*.

        If *ref_eof* is provided, the model data is recomputed unless the cache
        already contains common-base variables and a ref-oriented EOF.

        Parameters
        ----------
        file_prefix : str
            Prefix prepended to the filename, e.g.
            ``"NAO_DJF_1985-2014_HadISST2_"``.
        overwrite : bool
            If True, ignore any cached file and recompute from scratch.
        """
        path = self._nc_path(out_dir, case_name, prefix=file_prefix)

        if os.path.exists(path) and not overwrite:
            ds = self.load_result(path)
            # Cache hit — only reuse if proj/anomaly variables are already present.
            has_proj = ref_eof is None or {"slope_proj", "eof_lr_proj", "intercept_proj", "frac_proj"}.issubset(ds.data_vars)
            has_pcmdi_common_base = (
                ref_eof is None
                or ds.attrs.get("common_base_projection_method") == "pcmdi_metrics.gain_pseudo_pcs"
            )
            has_oriented_eof = ref_eof is None or ds.attrs.get("eof_oriented_to_ref") in {
                "true",
                "pcmdi_sign_rule_only",
            }
            has_anom = not self._requires_saved_anomaly(mode_name, ds)
            has_local_pcmdi = not self._requires_local_pcmdi_eof(mode_name, ds)
            has_eof_lr = {"eof_lr", "intercept"}.issubset(ds.data_vars)
            has_raw_pc = ds.attrs.get("pc_saved_as") == "raw_pcmdi_pc"
            has_raw_proj_pc = ref_eof is None or ds.attrs.get("pc_proj_saved_as") == "raw_pcmdi_pseudo_pc"
            has_pcmdi_reg = (
                ds.attrs.get("pcmdi_regression_method")
                == "pcmdi_metrics.linear_regression_on_globe_for_teleconnection"
            )
            has_units = (
                self.cfg.get_mode_info(mode_name).get("var") != "psl"
                or ds.attrs.get("analysis_units") == "hPa"
            )
            has_reader = ds.attrs.get("analysis_reader") == "pcmdi_metrics.read_data_in"
            has_ref_regrid = (
                ref_eof is None
                or ds.attrs.get("common_base_regrid") == "pcmdi_metrics.regrid_to_reference_grid"
            )
            if (
                has_proj and has_pcmdi_common_base and has_oriented_eof and has_anom
                and has_local_pcmdi and has_eof_lr and has_raw_pc and has_raw_proj_pc
                and has_pcmdi_reg and has_units and has_reader and has_ref_regrid
            ):
                print(f"  Loading {case_name} from cache: {path}")
                return ds
            missing = []
            if not has_proj:
                missing.append("proj eof_lr/intercept/frac variables")
            if not has_pcmdi_common_base:
                missing.append("PCMDI common-base projection stamp")
            if not has_oriented_eof:
                missing.append("ref-oriented EOF")
            if not has_anom:
                missing.append("psl_anom")
            if not has_local_pcmdi:
                missing.append("local PCMDI EOF stamp")
            if not has_eof_lr:
                missing.append("eof_lr/intercept variables")
            if not has_raw_pc:
                missing.append("raw PCMDI pc storage")
            if not has_raw_proj_pc:
                missing.append("raw PCMDI pseudo-pc storage")
            if not has_pcmdi_reg:
                missing.append("PCMDI regression method")
            if not has_units:
                missing.append("PCMDI hPa unit normalization")
            if not has_reader:
                missing.append("PCMDI read_data_in reader")
            if not has_ref_regrid:
                missing.append("PCMDI reference-grid CBF regrid")
            print(f"  Cache for {case_name} missing {', '.join(missing)}, recomputing …")
        elif os.path.exists(path) and overwrite:
            print(f"  Overwrite=True — recomputing {case_name} …")

        print(f"  Computing {case_name} …")
        ds = self.analyze_model(
            mode_name=mode_name,
            case_name=case_name,
            period=period,
            season=season,
            ref_eof=ref_eof,
            ref_solver=ref_solver,
            ref_reverse_sign=ref_reverse_sign,
            ref_grid_global=ref_grid_global,
        )
        self.save_result(ds, path)
        return ds

    def analyze_or_load_all(
        self,
        mode_name: str,
        model_list: Optional[List[str]],
        out_dir: str,
        custom_obs_path: Optional[str] = None,
        custom_obs_name: Optional[str] = None,
        period: Optional[Tuple[int, int]] = None,
        season: Optional[str] = None,
        overwrite: bool = False,
    ) -> Tuple[xr.Dataset, Dict[str, xr.Dataset]]:
        """
        Convenience wrapper: process (or load) obs + all models.

        Parameters
        ----------
        overwrite : bool
            If True, ignore cached files and recompute everything from scratch.

        Returns
        -------
        obs_ds : xr.Dataset
        model_results : dict[case_name -> xr.Dataset]
        """
        os.makedirs(out_dir, exist_ok=True)
        if model_list is None:
            model_list = list(self.cfg.datasets.keys())

        # Build a unique file prefix that encodes all analysis parameters.
        # Format: {mode}_{season}_{y0}-{y1}_{obs_name}_
        obs_info  = self._resolve_observation(mode_name, custom_obs_path, custom_obs_name)
        obs_name  = obs_info.get("name", "obs")
        _season   = season if season is not None else "default"
        _period   = f"{period[0]}-{period[1]}" if period is not None else "full"
        file_prefix = f"{mode_name}_{_season}_{_period}_{obs_name}_"

        self._last_pcmdi_solver = None
        self._last_pcmdi_reverse_sign = False
        self._last_ref_grid_global = None
        self._last_runtime_key = None

        obs_ds = self.analyze_or_load_observation(
            mode_name=mode_name,
            out_dir=out_dir,
            custom_obs_path=custom_obs_path,
            custom_obs_name=custom_obs_name,
            period=period,
            season=season,
            file_prefix=file_prefix,
            overwrite=overwrite,
        )

        # Use the obs EOF/EOF-solver as the common spatial basis for all models.
        # For index-type modes there is no EOF, so ref_eof/ref_solver stay None.
        ref_eof: Optional[xr.DataArray] = obs_ds.get("eof")
        var_name = self.cfg.get_mode_info(mode_name).get("var")
        runtime_key = (mode_name, period, season, obs_name)
        use_runtime_obs_state = (
            self._last_runtime_key == runtime_key
            and self._last_pcmdi_solver is not None
            and self._last_ref_grid_global is not None
        )
        ref_grid_global: Optional[xr.DataArray] = (
            self._last_ref_grid_global
            if use_runtime_obs_state
            else (obs_ds.get(f"{var_name}_anom") if var_name else None)
        )
        if ref_eof is not None:
            if use_runtime_obs_state:
                print("  Using in-memory obs EOF solver/reference grid for PCMDI CBF projection.")
                ref_solver = self._last_pcmdi_solver
                ref_reverse_sign = self._last_pcmdi_reverse_sign
            else:
                ref_solver, ref_reverse_sign = self.compute_observation_pcmdi_solver(
                    mode_name=mode_name,
                    custom_obs_path=custom_obs_path,
                    custom_obs_name=custom_obs_name,
                    period=period,
                    season=season,
                )
        else:
            ref_solver, ref_reverse_sign = None, False

        model_results: Dict[str, xr.Dataset] = {}
        for case_name in model_list:
            model_results[case_name] = self.analyze_or_load_model(
                mode_name=mode_name,
                case_name=case_name,
                out_dir=out_dir,
                period=period,
                season=season,
                ref_eof=ref_eof,
                ref_solver=ref_solver,
                ref_reverse_sign=ref_reverse_sign,
                ref_grid_global=ref_grid_global,
                file_prefix=file_prefix,
                overwrite=overwrite,
            )

        return obs_ds, model_results

    # ------------------------------------------------------------------
    # Internal dispatch
    # ------------------------------------------------------------------
    def _analyze_preprocessed(
        self,
        da: xr.DataArray,
        mode_info: dict,
        ref_eof: Optional[xr.DataArray] = None,
        ref_solver: Optional[object] = None,
        ref_reverse_sign: bool = False,
        ref_grid_global: Optional[xr.DataArray] = None,
        da_global: Optional[xr.DataArray] = None,
    ) -> xr.Dataset:
        def _with_saved_anomaly(ds: xr.Dataset) -> xr.Dataset:
            var_name = mode_info.get("var")
            if var_name != "psl":
                return ds
            anom = da_global if da_global is not None else da
            anom_name = f"{var_name}_anom"
            anom = anom.rename(anom_name)
            anom.attrs["long_name"] = "Full-domain PSL anomaly used for regression diagnostics"
            anom.attrs["description"] = (
                "Monthly climatology removed, then aggregated using the configured "
                "season/frequency and period before EOF analysis."
            )
            merged = xr.merge([ds, anom.to_dataset()])
            merged.attrs["analysis_units"] = str(anom.attrs.get("units", ""))
            merged.attrs["pcmdi_units_adjustment"] = str(anom.attrs.get("pcmdi_units_adjustment", "none"))
            merged.attrs["analysis_reader"] = str(anom.attrs.get("analysis_reader", "xarray.open_mfdataset"))
            return merged

        mode_type = mode_info.get("type", "eof")
        if mode_type == "eof":
            da = da.copy()
            da.attrs["mode_name"] = mode_info.get("_name", "")
            ds = self.run_eof(da, eof_num=mode_info["eof_num"], da_global=da_global)
            if ref_eof is not None:
                if mode_info.get("orient_independent_eof_to_reference", False):
                    ds = self.orient_eof_to_reference(ds, ref_eof)
                else:
                    ds.attrs["eof_oriented_to_ref"] = "pcmdi_sign_rule_only"
                ds_proj = self.project_onto_eof(
                    da,
                    ref_eof,
                    da_global=da_global,
                    ref_solver=ref_solver,
                    ref_reverse_sign=ref_reverse_sign,
                    ref_grid_global=ref_grid_global,
                    eof_num=mode_info["eof_num"],
                    mode_name=mode_info.get("_name", ""),
                    mode_info=mode_info,
                )
                ds = xr.merge([ds, ds_proj])
                for attr_name in (
                    "common_base_projection_method",
                    "common_base_regrid",
                    "pc_proj_saved_as",
                ):
                    if attr_name in ds_proj.attrs:
                        ds.attrs[attr_name] = ds_proj.attrs[attr_name]
            return _with_saved_anomaly(ds)
        if mode_type == "index":
            ds = self.run_index(
                da,
                standardize=mode_info.get("standardize", False),
                running_mean=mode_info.get("running_mean"),
            )
            return _with_saved_anomaly(ds)
        raise ValueError(f"Unsupported mode type: {mode_type!r}")
