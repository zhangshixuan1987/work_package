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
      Common-base proj  : pc_proj, slope_proj, slope_pval_proj, corr_proj
    The common basis is always the obs EOF (passed automatically via
    analyze_or_load_all → ref_eof = obs_ds["eof"]).
"""

import os
import glob
import json
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import xarray as xr
from scipy import stats
from eofs.xarray import Eof


class ModeConfigManager:
    """Load and query mode / dataset configuration from a JSON file."""

    def __init__(self, config_path: str):
        with open(config_path, "r") as f:
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
    def compute_time_aggregation(da: xr.DataArray, season: Optional[str]) -> xr.DataArray:
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
            return djf.groupby("winter_year").mean("time").rename({"winter_year": "time"})

        if season in ("MAM", "JJA", "SON"):
            sub = da.where(da["time"].dt.season == season, drop=True)
            return sub.groupby("time.year").mean("time").rename({"year": "time"})

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
        da = self.normalize_lon(da)
        da = self.subset_time(da, period=period)
        if subset_space:
            da = self.subset_latlon(
                da,
                lat_bnds=mode_info.get("lat_bnds"),
                lon_bnds=mode_info.get("lon_bnds"),
            )
        da = self.compute_monthly_anomaly(da)
        season = season_override if season_override is not None else mode_info.get("season", "monthly")
        da = self.compute_time_aggregation(da, season)
        return da

    # ------------------------------------------------------------------
    # EOF analysis
    # ------------------------------------------------------------------
    @staticmethod
    def _slope_and_pval(
        da: xr.DataArray, pc_std: xr.DataArray
    ) -> Tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
        """
        Compute the OLS regression slope of *da* onto *pc_std* at each grid point,
        the Pearson correlation, and the two-tailed p-value of that correlation.

        Returns
        -------
        slope : xr.DataArray  (lat, lon)
        corr  : xr.DataArray  (lat, lon)
        pval  : xr.DataArray  (lat, lon)  – two-tailed p-value via t-test
        """
        n       = da.sizes["time"]
        da_anom = da - da.mean("time")
        pc_var  = (pc_std ** 2).mean("time")

        # OLS slope: cov(field, pc) / var(pc)  [pc is standardised so var≈1]
        slope = (da_anom * pc_std).mean("time") / pc_var

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
        return slope, corr, pval

    def project_onto_eof(
        self, da: xr.DataArray, ref_eof: xr.DataArray,
        da_global: Optional[xr.DataArray] = None,
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
        pc_proj = (pc_raw - pc_raw.mean("time")) / pc_raw.std("time")

        reg_field = da_global if da_global is not None else da
        slope_proj, corr_proj, pval_proj = self._slope_and_pval(reg_field, pc_proj)

        out = xr.Dataset({
            "pc_proj":         pc_proj,
            "slope_proj":      slope_proj,
            "slope_pval_proj": pval_proj,
            "corr_proj":       corr_proj,
        })
        out["pc_proj"        ].attrs["long_name"] = "Standardized PC from projection onto obs EOF"
        out["slope_proj"     ].attrs["long_name"] = "Regression slope onto obs-projected PC"
        out["slope_pval_proj"].attrs["long_name"] = "Two-tailed p-value of obs-projected regression slope"
        out["slope_pval_proj"].attrs["description"] = "Derived from t-test on Pearson r; df = n_time - 2"
        out["corr_proj"      ].attrs["long_name"] = "Pearson r with obs-projected standardized PC"
        return out

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
        solver  = Eof(da * weights)

        eof_w = solver.eofs(neofs=eof_num)[eof_num - 1]
        eof   = eof_w / weights
        pc    = solver.pcs(npcs=eof_num, pcscaling=1)[:, eof_num - 1]
        frac  = solver.varianceFraction(neigs=eof_num)[eof_num - 1]

        pc_std = (pc - pc.mean("time")) / pc.std("time")
        reg_field = da_global if da_global is not None else da
        slope, corr, pval = self._slope_and_pval(reg_field, pc_std)

        # Expand the regional EOF to the global grid (NaN outside the mode
        # domain) so that all variables share the same lat/lon coordinates and
        # can be stored together in a single NetCDF file without coordinate
        # conflicts between the regional eof and the global slope/corr maps.
        if da_global is not None:
            ref_spatial = da_global.isel(time=0, drop=True)
            eof = eof.reindex_like(ref_spatial, fill_value=np.nan)

        out = xr.Dataset({
            "pc":         pc_std,
            "eof":        eof,
            "slope":      slope,
            "slope_pval": pval,
            "corr":       corr,
            "frac":       frac,
        })
        out["pc"        ].attrs["long_name"] = f"Standardized PC of EOF{eof_num}"
        out["eof"       ].attrs["long_name"] = f"EOF{eof_num} spatial pattern"
        out["slope"     ].attrs["long_name"] = f"Regression slope onto standardized EOF{eof_num} PC"
        out["slope_pval"].attrs["long_name"] = f"Two-tailed p-value of regression slope (EOF{eof_num})"
        out["slope_pval"].attrs["description"] = "Derived from t-test on Pearson r; df = n_time - 2"
        out["corr"      ].attrs["long_name"] = f"Pearson r with standardized EOF{eof_num} PC"
        out["frac"      ].attrs["long_name"] = f"Explained variance fraction of EOF{eof_num}"
        return out

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
    config_path : str
        Path to the JSON configuration file (datasets + modes).
    """

    def __init__(self, config_path: str):
        super().__init__()
        self.cfg = ModeConfigManager(config_path)

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
        mode_info  = self.cfg.get_mode_info(mode_name)
        obs_info   = self._resolve_observation(mode_name, custom_obs_path, custom_obs_name)
        raw        = self.read_obs_data(obs_info["data"], mode_info["var"])
        # Pre-process globally first (no spatial subsetting), then subset.
        # This guarantees da_global is the full-domain field and avoids calling
        # preprocess_field twice on the same lazy DataArray (which can silently
        # produce a regional result the second time due to lazy-graph sharing).
        da_global  = self.preprocess_field(raw, mode_info, period=period,
                                           season_override=season, subset_space=False)
        obs_pre    = self.subset_latlon(
            da_global,
            lat_bnds=mode_info.get("lat_bnds"),
            lon_bnds=mode_info.get("lon_bnds"),
        )
        return self._analyze_preprocessed(obs_pre, mode_info, da_global=da_global)

    def analyze_model(
        self,
        mode_name: str,
        case_name: str,
        period: Optional[Tuple[int, int]] = None,
        season: Optional[str] = None,
        ref_eof: Optional[xr.DataArray] = None,
    ) -> xr.Dataset:
        """Load, pre-process and analyse a single model case.

        Parameters
        ----------
        ref_eof : xr.DataArray, optional
            If provided (typically ``obs_ds["eof"]``), the field is also
            projected onto this common basis and the results are stored as
            ``slope_proj``, ``slope_pval_proj``, ``corr_proj``, ``pc_proj``.
        """
        mode_info    = self.cfg.get_mode_info(mode_name)
        dataset_info = self.cfg.get_dataset_info(case_name)
        raw          = self.read_model_data(dataset_info, mode_info["var"])
        # Pre-process globally first (no spatial subsetting), then subset.
        # This guarantees da_global is the full-domain field and avoids calling
        # preprocess_field twice on the same lazy DataArray (which can silently
        # produce a regional result the second time due to lazy-graph sharing).
        da_global    = self.preprocess_field(raw, mode_info, period=period,
                                             season_override=season, subset_space=False)
        model_pre    = self.subset_latlon(
            da_global,
            lat_bnds=mode_info.get("lat_bnds"),
            lon_bnds=mode_info.get("lon_bnds"),
        )
        return self._analyze_preprocessed(model_pre, mode_info, ref_eof=ref_eof, da_global=da_global)

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
        ds.to_netcdf(path)
        print(f"  Saved → {path}")

    @staticmethod
    def load_result(path: str) -> xr.Dataset:
        """Read an analysis Dataset from NetCDF into memory (closes file handle)."""
        with xr.open_dataset(path) as ds:
            return ds.load()
        return xr.open_dataset(path)

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
            print(f"  Loading obs from cache: {path}")
            return self.load_result(path)
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
        file_prefix: str = "",
        overwrite: bool = False,
    ) -> xr.Dataset:
        """
        Load model result from *out_dir* if the file already exists,
        otherwise run the analysis and save to *out_dir*.

        If *ref_eof* is provided and the cached file does not yet contain
        ``slope_proj``, the model data is recomputed from scratch so that both
        the independent-EOF and common-base variables are saved together.

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
            # Cache hit — only reuse if proj variables are already present (or not needed)
            if ref_eof is None or "slope_proj" in ds.data_vars:
                print(f"  Loading {case_name} from cache: {path}")
                return ds
            print(f"  Cache for {case_name} missing proj variables, recomputing …")
        elif os.path.exists(path) and overwrite:
            print(f"  Overwrite=True — recomputing {case_name} …")

        print(f"  Computing {case_name} …")
        ds = self.analyze_model(
            mode_name=mode_name,
            case_name=case_name,
            period=period,
            season=season,
            ref_eof=ref_eof,
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

        # Use the obs EOF as the common spatial basis for all models.
        # For index-type modes there is no EOF, so ref_eof stays None.
        ref_eof: Optional[xr.DataArray] = obs_ds.get("eof")

        model_results: Dict[str, xr.Dataset] = {}
        for case_name in model_list:
            model_results[case_name] = self.analyze_or_load_model(
                mode_name=mode_name,
                case_name=case_name,
                out_dir=out_dir,
                period=period,
                season=season,
                ref_eof=ref_eof,
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
        da_global: Optional[xr.DataArray] = None,
    ) -> xr.Dataset:
        mode_type = mode_info.get("type", "eof")
        if mode_type == "eof":
            ds = self.run_eof(da, eof_num=mode_info["eof_num"], da_global=da_global)
            if ref_eof is not None:
                ds_proj = self.project_onto_eof(da, ref_eof, da_global=da_global)
                ds = xr.merge([ds, ds_proj])
            return ds
        if mode_type == "index":
            return self.run_index(
                da,
                standardize=mode_info.get("standardize", False),
                running_mean=mode_info.get("running_mean"),
            )
        raise ValueError(f"Unsupported mode type: {mode_type!r}")
