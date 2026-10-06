"""Cached ensemble atmosphere-land cross-correlation diagnostics."""
from __future__ import annotations

import glob
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence, Dict, Tuple

import numpy as np
import pandas as pd
import xarray as xr
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature

@dataclass
class VariableSpec:
    """Specification for one analysis variable."""
    name: str
    level_index: Optional[int] = None   # None for 2D variables like PS
    dim_name: Optional[str] = None      # e.g. "lev", "levtot"


class AtmosLandCrossCorrelationAnalyzer:
    """
    Analyze ensemble-based atmosphere-land cross correlations using:
      - atmosphere on native EAM ne30 ncol
      - land aggregated from ELM column -> gridcell

    Main workflow:
      1) discover paired eam.i / elm.r files
      2) aggregate land columns to land gridcells
      3) build nearest-neighbor mapping land_gridcell -> atm_ncol
      4) compute ensemble cross correlations
    """

    def __init__(
        self,
        base_dir: str,
        time_tag: str,
        ensemble_ids: Optional[Sequence[str]] = None,
    ):
        self.base_dir = base_dir
        self.time_tag = time_tag
        self.ensemble_ids = ensemble_ids

        self.eam_files: Dict[str, str] = {}
        self.elm_files: Dict[str, str] = {}

        self.land_gridcell_ids: Optional[np.ndarray] = None
        self.land_lat: Optional[np.ndarray] = None
        self.land_lon: Optional[np.ndarray] = None

        self.atm_lat: Optional[np.ndarray] = None
        self.atm_lon: Optional[np.ndarray] = None

        self.land_to_atm_index: Optional[np.ndarray] = None
        self.land_to_atm_distance_km: Optional[np.ndarray] = None

        self._discover_files()

    # ------------------------------------------------------------------
    # File discovery
    # ------------------------------------------------------------------
    def _discover_files(self):
        patterns = (
            (
                os.path.join(self.base_dir, "EN*", "archive", "post", "init_rgd", self.time_tag, f"*.eam.i.{self.time_tag}.nc"),
                os.path.join(self.base_dir, "EN*", "archive", "post", "init_rgd", self.time_tag, f"*.elm.r.{self.time_tag}.nc"),
            ),
            (
                os.path.join(self.base_dir, f"*.eam.i.{self.time_tag}.nc"),
                os.path.join(self.base_dir, f"*.elm.r.{self.time_tag}.nc"),
            ),
        )
        eam_all, elm_all = [], []
        for eam_pattern, elm_pattern in patterns:
            eam_all = sorted(glob.glob(eam_pattern))
            elm_all = sorted(glob.glob(elm_pattern))
            if eam_all or elm_all:
                break

        def extract_ens_id(path: str) -> str:
            fname = os.path.basename(path)
            parts = fname.split(".")
            for p in parts:
                if p.startswith("EN"):
                    return p
            raise ValueError(f"Could not parse ensemble ID from {fname}")

        eam_map = {extract_ens_id(f): f for f in eam_all}
        elm_map = {extract_ens_id(f): f for f in elm_all}

        common_ids = sorted(set(eam_map) & set(elm_map))
        if self.ensemble_ids is not None:
            common_ids = [eid for eid in self.ensemble_ids if eid in common_ids]

        if not common_ids:
            raise ValueError("No matched eam.i / elm.r ensemble pairs found.")

        self.eam_files = {eid: eam_map[eid] for eid in common_ids}
        self.elm_files = {eid: elm_map[eid] for eid in common_ids}

        print(f"Found {len(common_ids)} matched ensemble pairs.")

    # ------------------------------------------------------------------
    # Coordinate helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _normalize_lon(lon: np.ndarray) -> np.ndarray:
        """Convert longitude to [-180, 180)."""
        lon = np.asarray(lon, dtype=float)
        return ((lon + 180.0) % 360.0) - 180.0

    @staticmethod
    def _latlon_to_xyz(lat_deg: np.ndarray, lon_deg: np.ndarray) -> np.ndarray:
        """Convert lat/lon to unit-sphere Cartesian coordinates."""
        lat = np.deg2rad(lat_deg)
        lon = np.deg2rad(lon_deg)
        x = np.cos(lat) * np.cos(lon)
        y = np.cos(lat) * np.sin(lon)
        z = np.sin(lat)
        return np.column_stack((x, y, z))

    @staticmethod
    def _great_circle_distance_km(
        lat1: np.ndarray, lon1: np.ndarray, lat2: np.ndarray, lon2: np.ndarray
    ) -> np.ndarray:
        """Vectorized great-circle distance in km."""
        r_earth = 6371.0
        lat1 = np.deg2rad(lat1)
        lon1 = np.deg2rad(lon1)
        lat2 = np.deg2rad(lat2)
        lon2 = np.deg2rad(lon2)

        dlat = lat2 - lat1
        dlon = lon2 - lon1

        a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0) ** 2
        c = 2.0 * np.arcsin(np.sqrt(a))
        return r_earth * c

    @staticmethod
    def _infer_atm_latlon_names(ds: xr.Dataset) -> Tuple[str, str]:
        candidates = [
            ("lat", "lon"),
            ("lats", "lons"),
            ("grid_center_lat", "grid_center_lon"),
        ]
        for lat_name, lon_name in candidates:
            if lat_name in ds.variables and lon_name in ds.variables:
                return lat_name, lon_name
        raise ValueError("Could not infer atmosphere lat/lon variable names.")

    # ------------------------------------------------------------------
    # Low-level extraction helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _extract_1d_column_values(ds: xr.Dataset, var_spec: VariableSpec) -> np.ndarray:
        """
        Extract a land variable as a 1D numpy array over 'column',
        avoiding xarray coordinate slicing issues.
        """
        da = ds[var_spec.name]
        dims = da.dims
        arr = np.asarray(da.values)

        if "time" in dims:
            tidx = dims.index("time")
            if arr.shape[tidx] < 1:
                raise ValueError(f"{var_spec.name}: time dimension has length 0.")
            arr = np.take(arr, indices=0, axis=tidx)
            dims = tuple(d for d in dims if d != "time")

        if var_spec.level_index is not None:
            if var_spec.dim_name is None:
                raise ValueError(f"{var_spec.name}: dim_name is required when level_index is given.")
            if var_spec.dim_name not in dims:
                raise ValueError(
                    f"{var_spec.name}: dim '{var_spec.dim_name}' not found after time slicing. "
                    f"Remaining dims: {dims}"
                )
            lidx = dims.index(var_spec.dim_name)
            if arr.shape[lidx] <= var_spec.level_index:
                raise IndexError(
                    f"{var_spec.name}: level_index={var_spec.level_index} out of bounds "
                    f"for dim '{var_spec.dim_name}' with size {arr.shape[lidx]}."
                )
            arr = np.take(arr, indices=var_spec.level_index, axis=lidx)
            dims = tuple(d for d in dims if d != var_spec.dim_name)

        if dims != ("column",):
            raise ValueError(
                f"{var_spec.name} did not reduce to 1D column field. Remaining dims: {dims}"
            )

        return np.asarray(arr)

    @staticmethod
    def _extract_1d_ncol_values(ds: xr.Dataset, var_spec: VariableSpec, source_name: str = "") -> np.ndarray:
        """
        Extract an atmosphere variable as a 1D numpy array over 'ncol',
        avoiding xarray coordinate slicing issues.
        """
        da = ds[var_spec.name]
        dims = da.dims
        arr = np.asarray(da.values)

        if "time" in dims:
            tidx = dims.index("time")
            if arr.shape[tidx] < 1:
                raise ValueError(
                    f"{source_name}: variable {var_spec.name} has time dimension of length 0. "
                    f"dims={dims}, shape={arr.shape}"
                )
            arr = np.take(arr, indices=0, axis=tidx)
            dims = tuple(d for d in dims if d != "time")

        if var_spec.level_index is not None:
            if var_spec.dim_name is None:
                raise ValueError(f"{source_name}: {var_spec.name}: dim_name is required when level_index is given.")
            if var_spec.dim_name not in dims:
                raise ValueError(
                    f"{source_name}: {var_spec.name}: dim '{var_spec.dim_name}' not found after time slicing. "
                    f"Remaining dims: {dims}"
                )
            lidx = dims.index(var_spec.dim_name)
            if arr.shape[lidx] <= var_spec.level_index:
                raise IndexError(
                    f"{source_name}: {var_spec.name}: level_index={var_spec.level_index} out of bounds "
                    f"for dim '{var_spec.dim_name}' with size {arr.shape[lidx]}."
                )
            arr = np.take(arr, indices=var_spec.level_index, axis=lidx)
            dims = tuple(d for d in dims if d != var_spec.dim_name)

        if dims != ("ncol",):
            raise ValueError(
                f"{source_name}: {var_spec.name} did not reduce to 1D ncol field. "
                f"Remaining dims: {dims}, shape={arr.shape}"
            )

        return np.asarray(arr)

    # ------------------------------------------------------------------
    # Land aggregation
    # ------------------------------------------------------------------
    def _aggregate_elm_to_gridcell(
        self,
        ds_elm: xr.Dataset,
        var_spec: VariableSpec,
        active_only: bool = True,
    ) -> pd.DataFrame:
        """
        Aggregate one ELM variable from column -> gridcell.
        Returns a DataFrame with columns: gc, lat, lon, value
        """
        values = self._extract_1d_column_values(ds_elm, var_spec)

        active = (
            ds_elm["cols1d_active"].values == 1
            if active_only
            else np.ones(ds_elm.sizes["column"], dtype=bool)
        )
        gc = ds_elm["cols1d_gridcell_index"].values.astype(int)

        lat = np.asarray(ds_elm["cols1d_lat"].values)
        lon = self._normalize_lon(np.asarray(ds_elm["cols1d_lon"].values))

        df = pd.DataFrame(
            {
                "gc": gc[active],
                "lat": lat[active],
                "lon": lon[active],
                "value": values[active],
            }
        )

        out = (
            df.groupby("gc", as_index=False)
            .agg(
                {
                    "lat": "mean",
                    "lon": "mean",
                    "value": "mean",
                }
            )
            .sort_values("gc")
            .reset_index(drop=True)
        )
        return out

    def build_land_grid(self, land_var: VariableSpec, active_only: bool = True):
        """Build static land gridcell coordinates from the first ensemble file."""
        first_eid = next(iter(self.elm_files))
        with xr.open_dataset(self.elm_files[first_eid]) as ds_elm:
            land_df = self._aggregate_elm_to_gridcell(ds_elm, land_var, active_only=active_only)

        self.land_gridcell_ids = land_df["gc"].to_numpy()
        self.land_lat = land_df["lat"].to_numpy()
        self.land_lon = land_df["lon"].to_numpy()

        print(f"Built land grid with {len(self.land_gridcell_ids)} aggregated gridcells.")

    # ------------------------------------------------------------------
    # Atmosphere grid + nearest-neighbor map
    # ------------------------------------------------------------------
    def build_land_to_atm_mapping(self):
        """Build nearest-neighbor map from aggregated land gridcells to atmosphere ncol."""
        if self.land_lat is None or self.land_lon is None:
            raise RuntimeError("Call build_land_grid() first.")

        first_eid = next(iter(self.eam_files))
        with xr.open_dataset(self.eam_files[first_eid]) as ds_eam:
            lat_name, lon_name = self._infer_atm_latlon_names(ds_eam)
            self.atm_lat = np.asarray(ds_eam[lat_name].values)
            self.atm_lon = self._normalize_lon(np.asarray(ds_eam[lon_name].values))

        atm_xyz = self._latlon_to_xyz(self.atm_lat, self.atm_lon)
        land_xyz = self._latlon_to_xyz(self.land_lat, self.land_lon)

        tree = cKDTree(atm_xyz)
        _, idx = tree.query(land_xyz, k=1)

        self.land_to_atm_index = idx.astype(int)
        self.land_to_atm_distance_km = self._great_circle_distance_km(
            self.land_lat,
            self.land_lon,
            self.atm_lat[self.land_to_atm_index],
            self.atm_lon[self.land_to_atm_index],
        )

        print("Built land -> nearest atmosphere mapping.")
        print(f"Mean distance: {np.nanmean(self.land_to_atm_distance_km):.2f} km")
        print(f"Max distance : {np.nanmax(self.land_to_atm_distance_km):.2f} km")

    # ------------------------------------------------------------------
    # Ensemble extraction
    # ------------------------------------------------------------------
    def get_land_ensemble_matrix(
        self,
        land_var: VariableSpec,
        active_only: bool = True,
    ) -> xr.DataArray:
        """
        Return land values as DataArray with dims:
            ensemble, land_point
        """
        if self.land_gridcell_ids is None:
            raise RuntimeError("Call build_land_grid() first.")

        rows = []
        ensemble_ids = []

        for eid in self.elm_files:
            with xr.open_dataset(self.elm_files[eid]) as ds_elm:
                land_df = self._aggregate_elm_to_gridcell(ds_elm, land_var, active_only=active_only)

            land_df = land_df.set_index("gc").reindex(self.land_gridcell_ids)
            rows.append(land_df["value"].to_numpy())
            ensemble_ids.append(eid)

        arr = np.vstack(rows)
        return xr.DataArray(
            arr,
            dims=("ensemble", "land_point"),
            coords={
                "ensemble": ensemble_ids,
                "land_point": np.arange(arr.shape[1]),
                "gridcell": ("land_point", self.land_gridcell_ids),
                "lat": ("land_point", self.land_lat),
                "lon": ("land_point", self.land_lon),
            },
            name=land_var.name,
        )

    def get_atm_ensemble_matrix(
        self,
        atm_var: VariableSpec,
    ) -> xr.DataArray:
        """
        Return atmosphere values as DataArray with dims:
            ensemble, ncol
        """
        if self.atm_lat is None or self.atm_lon is None:
            raise RuntimeError("Call build_land_to_atm_mapping() first.")

        rows = []
        ensemble_ids = []

        for eid in self.eam_files:
            fpath = self.eam_files[eid]
            try:
                with xr.open_dataset(fpath, decode_times=False) as ds_eam:
                    vals = self._extract_1d_ncol_values(ds_eam, atm_var, source_name=f"{eid} | {fpath}")
            except Exception as e:
                raise RuntimeError(f"Failed reading atmosphere variable {atm_var.name} from {eid}: {fpath}\n{e}") from e

            rows.append(vals)
            ensemble_ids.append(eid)

        arr = np.vstack(rows)
        return xr.DataArray(
            arr,
            dims=("ensemble", "ncol"),
            coords={
                "ensemble": ensemble_ids,
                "ncol": np.arange(arr.shape[1]),
                "lat": ("ncol", self.atm_lat),
                "lon": ("ncol", self.atm_lon),
            },
            name=atm_var.name,
        )

    def filter_valid_ensembles(
        self,
        atm_var: VariableSpec,
        land_var: Optional[VariableSpec] = None,
        active_only: bool = True,
        verbose: bool = True,
    ):
        """
        Remove ensemble members whose atmosphere or land files are invalid
        for the requested variables.
        """
        valid_eam = {}
        valid_elm = {}

        for eid in list(self.eam_files.keys()):
            eam_ok = True
            elm_ok = True

            # Check atmosphere file
            try:
                with xr.open_dataset(self.eam_files[eid], decode_times=False) as ds_eam:
                    _ = self._extract_1d_ncol_values(
                        ds_eam, atm_var, source_name=f"{eid} | {self.eam_files[eid]}"
                    )
            except Exception as e:
                eam_ok = False
                if verbose:
                    print(f"[DROP ATM] {eid}: {e}")

            # Check land file if requested
            if land_var is not None:
                try:
                    with xr.open_dataset(self.elm_files[eid], decode_times=False) as ds_elm:
                        _ = self._aggregate_elm_to_gridcell(ds_elm, land_var, active_only=active_only)
                except Exception as e:
                    elm_ok = False
                    if verbose:
                        print(f"[DROP LND] {eid}: {e}")

            if eam_ok and elm_ok:
                valid_eam[eid] = self.eam_files[eid]
                valid_elm[eid] = self.elm_files[eid]

        self.eam_files = valid_eam
        self.elm_files = valid_elm

        print(f"Retained {len(self.eam_files)} valid ensemble pairs.")

    # ------------------------------------------------------------------
    # Correlation methods
    # ------------------------------------------------------------------
    @staticmethod
    def _corr_2d(x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """
        Correlate x and y across ensemble axis (axis=0).
        x, y shape: (nens, npoint)
        Returns r(point), shape (npoint,)
        """
        x = x - np.nanmean(x, axis=0, keepdims=True)
        y = y - np.nanmean(y, axis=0, keepdims=True)

        num = np.nansum(x * y, axis=0)
        den = np.sqrt(np.nansum(x * x, axis=0) * np.nansum(y * y, axis=0))

        r = np.full(num.shape, np.nan, dtype=float)
        mask = den > 0
        r[mask] = num[mask] / den[mask]
        return r

    def compute_collocated_correlation(
        self,
        atm_var: VariableSpec,
        land_var: VariableSpec,
        active_only: bool = True,
    ) -> xr.Dataset:
        """
        Compute correlation across ensemble members between:
          - land variable at each aggregated land gridcell
          - nearest collocated atmosphere variable at matched atm ncol
        """
        if self.land_to_atm_index is None:
            raise RuntimeError("Call build_land_to_atm_mapping() first.")

        land_da = self.get_land_ensemble_matrix(land_var, active_only=active_only)
        atm_da = self.get_atm_ensemble_matrix(atm_var)

        atm_collocated = atm_da.isel(ncol=xr.DataArray(self.land_to_atm_index, dims="land_point"))
        atm_collocated = atm_collocated.assign_coords(
            land_point=("land_point", np.arange(len(self.land_to_atm_index)))
        )
        atm_collocated = atm_collocated.transpose("ensemble", "land_point")

        r = self._corr_2d(atm_collocated.values, land_da.values)

        ds_out = xr.Dataset(
            {
                "correlation": ("land_point", r),
                "atm_ncol": ("land_point", self.land_to_atm_index),
                "distance_km": ("land_point", self.land_to_atm_distance_km),
                "land_value_mean": ("land_point", np.nanmean(land_da.values, axis=0)),
                "atm_value_mean": ("land_point", np.nanmean(atm_collocated.values, axis=0)),
            },
            coords={
                "land_point": np.arange(len(r)),
                "gridcell": ("land_point", self.land_gridcell_ids),
                "lat": ("land_point", self.land_lat),
                "lon": ("land_point", self.land_lon),
            },
            attrs={
                "atm_var": atm_var.name,
                "land_var": land_var.name,
                "atm_level_index": atm_var.level_index if atm_var.level_index is not None else -1,
                "land_level_index": land_var.level_index if land_var.level_index is not None else -1,
            },
        )
        return ds_out

    def compute_single_landpoint_to_all_atm_correlation(
        self,
        atm_var: VariableSpec,
        land_var: VariableSpec,
        land_point_index: int,
        active_only: bool = True,
    ) -> xr.Dataset:
        """
        For one land point, correlate its ensemble perturbation with ALL atmosphere ncol points.
        Useful for localization-distance diagnostics.
        """
        land_da = self.get_land_ensemble_matrix(land_var, active_only=active_only)
        atm_da = self.get_atm_ensemble_matrix(atm_var)

        y = land_da.isel(land_point=land_point_index).values[:, None]
        x = atm_da.values
        r = self._corr_2d(x, np.repeat(y, x.shape[1], axis=1))

        dist = self._great_circle_distance_km(
            np.full_like(self.atm_lat, self.land_lat[land_point_index], dtype=float),
            np.full_like(self.atm_lon, self.land_lon[land_point_index], dtype=float),
            self.atm_lat,
            self.atm_lon,
        )

        return xr.Dataset(
            {
                "correlation": ("ncol", r),
                "distance_km": ("ncol", dist),
            },
            coords={
                "ncol": np.arange(len(r)),
                "lat": ("ncol", self.atm_lat),
                "lon": ("ncol", self.atm_lon),
            },
            attrs={
                "land_point_index": int(land_point_index),
                "land_gridcell": int(self.land_gridcell_ids[land_point_index]),
                "land_lat": float(self.land_lat[land_point_index]),
                "land_lon": float(self.land_lon[land_point_index]),
                "atm_var": atm_var.name,
                "land_var": land_var.name,
            },
        )

class CrossCorrelationPlotter:
    """
    Visualization utilities for atmosphere-land cross-correlation diagnostics.

    Assumes 1-D variables indexed by land_point (or a flattened 1-D record axis),
    including:
      - correlation
      - distance_km
      - lat
      - lon

    Dataset attributes optionally used:
      - atm_var
      - land_var
      - atm_level_index
      - land_level_index
    """

    def __init__(self, nc_file):
        self.ds = xr.open_dataset(nc_file)

    # -------------------------------------------------
    # Internal helpers
    # -------------------------------------------------
    def _get_meta_title(self, include_levels=False):
        atm_var = self.ds.attrs.get("atm_var", "atm_var")
        land_var = self.ds.attrs.get("land_var", "land_var")

        if include_levels:
            atm_lev = self.ds.attrs.get("atm_level_index", "NA")
            land_lev = self.ds.attrs.get("land_level_index", "NA")
            return f"{atm_var}(lev={atm_lev}) vs {land_var}(lev={land_lev})"

        return f"{atm_var} vs {land_var}"

    def get_tag(self, time_tag=None):
        atm_var = self.ds.attrs.get("atm_var", "atm")
        land_var = self.ds.attrs.get("land_var", "land")
        atm_lev = self.ds.attrs.get("atm_level_index", "NA")
        land_lev = self.ds.attrs.get("land_level_index", "NA")

        tag = f"{atm_var}_lev{atm_lev}_vs_{land_var}_lev{land_lev}"
        if time_tag is not None:
            tag = f"{tag}_{time_tag}"
        return tag

    def _get_valid_arrays(self):
        lon = self.ds["lon"].values
        lat = self.ds["lat"].values
        corr = self.ds["correlation"].values

        mask = np.isfinite(lon) & np.isfinite(lat) & np.isfinite(corr)
        return lon[mask], lat[mask], corr[mask]

    def _get_valid_distance_corr(self, use_abs=False):
        dist = self.ds["distance_km"].values
        corr = self.ds["correlation"].values

        mask = np.isfinite(dist) & np.isfinite(corr)
        dist = dist[mask]
        corr = corr[mask]

        if use_abs:
            corr = np.abs(corr)

        return dist, corr

    def _subsample_arrays(self, *arrays, subsample=None, seed=42):
        if len(arrays) == 0:
            return arrays

        n = len(arrays[0])
        if not all(len(arr) == n for arr in arrays):
            raise ValueError("All input arrays must have the same length.")

        if subsample is None or subsample >= n:
            return arrays

        rng = np.random.default_rng(seed)
        idx = rng.choice(n, size=subsample, replace=False)
        return tuple(arr[idx] for arr in arrays)

    def _make_distance_bins(self, dist, bins):
        bin_edges = np.linspace(dist.min(), dist.max(), bins + 1)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        bin_ids = np.digitize(dist, bin_edges, right=False) - 1
        bin_ids[bin_ids == bins] = bins - 1
        return bin_edges, bin_centers, bin_ids

    def _finalize_figure(self, outfile=None, dpi=200, show=True):
        plt.tight_layout()
        if outfile is not None:
            plt.savefig(outfile, dpi=dpi, bbox_inches="tight")
            print(f"Saved figure: {outfile}")
        if show:
            plt.show()
        else:
            plt.close()

    # -------------------------------------------------
    # Basic scatter map with geographic context
    # -------------------------------------------------
    def plot_global_scatter(
        self,
        vmin=None,
        vmax=None,
        cmap="RdBu_r",
        point_size=1,
        subsample=None,
        seed=42,
        title=None,
        outfile=None,
        dpi=200,
        show=True,
        alpha=0.7,
        add_borders=True,
        add_land=False,
    ):
        lon, lat, corr = self._get_valid_arrays()
        lon, lat, corr = self._subsample_arrays(
            lon, lat, corr, subsample=subsample, seed=seed
        )

        # enforce symmetric color range around zero unless provided
        if vmin is None or vmax is None:
            vmax_auto = np.nanmax(np.abs(corr))
            vmin = -vmax_auto
            vmax = vmax_auto

        fig = plt.figure(figsize=(14, 5))
        ax = plt.axes(projection=ccrs.PlateCarree())

        sc = ax.scatter(
            lon,
            lat,
            c=corr,
            s=point_size,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            alpha=alpha,
            transform=ccrs.PlateCarree(),
        )

        ax.set_global()
        ax.coastlines(linewidth=0.8)

        if add_borders:
            ax.add_feature(cfeature.BORDERS, linewidth=0.4)

        if add_land:
            ax.add_feature(cfeature.LAND, facecolor="lightgray", alpha=0.15)

        gl = ax.gridlines(
            draw_labels=True,
            linewidth=0.5,
            color="gray",
            alpha=0.5,
            linestyle="--",
        )
        gl.top_labels = False
        gl.right_labels = False
        gl.xlocator = mticker.FixedLocator(np.arange(-180, 181, 60))
        gl.ylocator = mticker.FixedLocator(np.arange(-90, 91, 30))

        cb = plt.colorbar(
            sc,
            ax=ax,
            orientation="horizontal",
            pad=0.06,
            fraction=0.05,
        )
        cb.set_label("Correlation")

        ax.set_title(title or self._get_meta_title())
        self._finalize_figure(outfile=outfile, dpi=dpi, show=show)

    # -------------------------------------------------
    # Histogram
    # -------------------------------------------------
    def plot_histogram(self, bins=50, outfile=None, dpi=200, show=True):
        corr = self.ds["correlation"].values
        corr = corr[np.isfinite(corr)]

        plt.figure(figsize=(6, 4))
        plt.hist(corr, bins=bins)
        plt.xlabel("Correlation")
        plt.ylabel("Count")
        plt.title("Distribution of cross-correlations")
        self._finalize_figure(outfile=outfile, dpi=dpi, show=show)

    # -------------------------------------------------
    # Correlation vs latitude
    # -------------------------------------------------
    def plot_latitudinal_mean(
        self,
        lat_bins=np.linspace(-90, 90, 37),
        outfile=None,
        dpi=200,
        show=True,
    ):
        corr = self.ds["correlation"]
        lat = self.ds["lat"]

        ds_lat = corr.groupby_bins(lat, bins=lat_bins).mean()
        lat_mid = [b.mid for b in ds_lat["lat_bins"].values]

        plt.figure(figsize=(6, 4))
        plt.plot(lat_mid, ds_lat.values)
        plt.xlabel("Latitude")
        plt.ylabel("Mean correlation")
        plt.title("Latitudinal mean correlation")
        plt.grid(True)
        self._finalize_figure(outfile=outfile, dpi=dpi, show=show)

    # -------------------------------------------------
    # Raw distance vs correlation scatter
    # -------------------------------------------------
    def plot_distance_relationship(
        self,
        subsample=None,
        seed=42,
        outfile=None,
        dpi=200,
        show=True,
    ):
        dist, corr = self._get_valid_distance_corr(use_abs=False)
        dist, corr = self._subsample_arrays(
            dist, corr, subsample=subsample, seed=seed
        )

        plt.figure(figsize=(6, 4))
        plt.scatter(dist, corr, s=1, alpha=0.3)
        plt.xlabel("Distance (km)")
        plt.ylabel("Correlation")
        plt.title("Correlation vs land-atmosphere separation distance")
        self._finalize_figure(outfile=outfile, dpi=dpi, show=show)

    # -------------------------------------------------
    # Distance-binned mean
    # -------------------------------------------------
    def plot_binned_distance_mean(
        self,
        bins=30,
        use_abs=False,
        show_counts=False,
        outfile=None,
        count_outfile=None,
        dpi=200,
        show=True,
    ):
        dist, corr = self._get_valid_distance_corr(use_abs=use_abs)
        bin_edges, bin_centers, bin_ids = self._make_distance_bins(dist, bins)

        mean_vals = np.full(bins, np.nan)
        counts = np.zeros(bins, dtype=int)

        for i in range(bins):
            m = bin_ids == i
            counts[i] = m.sum()
            if counts[i] > 0:
                mean_vals[i] = corr[m].mean()

        plt.figure(figsize=(6, 4))
        plt.plot(bin_centers, mean_vals, marker="o")
        plt.xlabel("Distance (km)")
        plt.ylabel("Mean |correlation|" if use_abs else "Mean correlation")
        plt.title("Distance-binned correlation")
        plt.grid(True)
        self._finalize_figure(outfile=outfile, dpi=dpi, show=show)

        if show_counts:
            plt.figure(figsize=(6, 3))
            plt.bar(bin_centers, counts, width=np.diff(bin_edges))
            plt.xlabel("Distance (km)")
            plt.ylabel("Count")
            plt.title("Counts per distance bin")
            self._finalize_figure(outfile=count_outfile, dpi=dpi, show=show)

    # -------------------------------------------------
    # Distance-binned percentiles
    # -------------------------------------------------
    def plot_binned_distance_percentiles(
        self,
        bins=30,
        use_abs=False,
        show_counts=False,
        outfile=None,
        count_outfile=None,
        dpi=200,
        show=True,
    ):
        dist, corr = self._get_valid_distance_corr(use_abs=use_abs)
        bin_edges, bin_centers, bin_ids = self._make_distance_bins(dist, bins)

        p10 = np.full(bins, np.nan)
        p50 = np.full(bins, np.nan)
        p90 = np.full(bins, np.nan)
        counts = np.zeros(bins, dtype=int)

        for i in range(bins):
            m = bin_ids == i
            counts[i] = m.sum()
            if counts[i] > 0:
                vals = corr[m]
                p10[i] = np.percentile(vals, 10)
                p50[i] = np.percentile(vals, 50)
                p90[i] = np.percentile(vals, 90)

        plt.figure(figsize=(6, 4))
        plt.fill_between(bin_centers, p10, p90, alpha=0.3, label="10–90%")
        plt.plot(bin_centers, p50, linewidth=2, label="Median")
        plt.xlabel("Distance (km)")
        plt.ylabel("|Correlation|" if use_abs else "Correlation")
        plt.title("Distance-binned correlation distribution")
        plt.grid(True)
        plt.legend()
        self._finalize_figure(outfile=outfile, dpi=dpi, show=show)

        if show_counts:
            plt.figure(figsize=(6, 3))
            plt.bar(bin_centers, counts, width=np.diff(bin_edges))
            plt.xlabel("Distance (km)")
            plt.ylabel("Count")
            plt.title("Counts per distance bin")
            self._finalize_figure(outfile=count_outfile, dpi=dpi, show=show)

    # -------------------------------------------------
    # Fraction exceeding |r| threshold vs distance
    # -------------------------------------------------
    def plot_distance_exceedance(
        self,
        bins=30,
        threshold=0.2,
        show_counts=False,
        outfile=None,
        count_outfile=None,
        dpi=200,
        show=True,
    ):
        dist, corr = self._get_valid_distance_corr(use_abs=True)
        bin_edges, bin_centers, bin_ids = self._make_distance_bins(dist, bins)

        frac = np.full(bins, np.nan)
        counts = np.zeros(bins, dtype=int)

        for i in range(bins):
            m = bin_ids == i
            counts[i] = m.sum()
            if counts[i] > 0:
                frac[i] = np.mean(corr[m] > threshold)

        plt.figure(figsize=(6, 4))
        plt.plot(bin_centers, frac, marker="o")
        plt.xlabel("Distance (km)")
        plt.ylabel(f"Fraction with |r| > {threshold}")
        plt.title("Correlation exceedance vs distance")
        plt.grid(True)
        self._finalize_figure(outfile=outfile, dpi=dpi, show=show)

        if show_counts:
            plt.figure(figsize=(6, 3))
            plt.bar(bin_centers, counts, width=np.diff(bin_edges))
            plt.xlabel("Distance (km)")
            plt.ylabel("Count")
            plt.title("Counts per distance bin")
            self._finalize_figure(outfile=count_outfile, dpi=dpi, show=show)

    # -------------------------------------------------
    # Convenience wrapper for old usage
    # -------------------------------------------------
    def save_scatter(
        self,
        outfile,
        subsample=30000,
        vmin=None,
        vmax=None,
        cmap="RdBu_r",
        seed=42,
        dpi=200,
    ):
        self.plot_global_scatter(
            vmin=vmin,
            vmax=vmax,
            cmap=cmap,
            point_size=2,
            subsample=subsample,
            seed=seed,
            title=self._get_meta_title(include_levels=True),
            outfile=outfile,
            dpi=dpi,
            show=False,
        )

    # -------------------------------------------------
    # Optional close method
    # -------------------------------------------------
    def close(self):
        self.ds.close()

def correlation_cache_path(cache_dir, time_tag, atm_var, atm_level, land_var, land_level, mode="collocated"):
    def level_tag(value):
        return "surface" if value is None else f"lev{int(value):02d}"
    name = f"{mode}_{atm_var}_{level_tag(atm_level)}_vs_{land_var}_{level_tag(land_level)}_{time_tag}.nc"
    return Path(cache_dir) / name


def cache_cross_correlation(path, analyzer=None, atm_var=None, land_var=None, mode="collocated", land_point_index=0, active_only=True, force_compute=False):
    path = Path(path)
    if path.is_file() and not force_compute:
        return xr.load_dataset(path)
    if analyzer is None or atm_var is None or land_var is None:
        raise ValueError("analyzer, atm_var, and land_var are required for a missing cache")
    if mode == "collocated":
        ds = analyzer.compute_collocated_correlation(atm_var, land_var, active_only=active_only)
    elif mode == "landpoint_all_atm":
        ds = analyzer.compute_single_landpoint_to_all_atm_correlation(atm_var, land_var, land_point_index=land_point_index, active_only=active_only)
        ds.attrs.update(atm_level_index=atm_var.level_index if atm_var.level_index is not None else -1, land_level_index=land_var.level_index if land_var.level_index is not None else -1)
    else:
        raise ValueError("mode must be collocated or landpoint_all_atm")
    ds.attrs["diagnostic_mode"] = mode
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    ds.to_netcdf(tmp)
    os.replace(tmp, path)
    return ds
