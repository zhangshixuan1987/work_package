import os
import glob
from dataclasses import dataclass
from typing import Dict, Sequence, Optional, Tuple, Literal, List

import numpy as np
import xarray as xr


@dataclass(frozen=True)
class ModeFileSpec:
    mode: str                  # "NAM", "AMO", ...
    var: str                   # "psl" or "ts"
    eof: str                   # "eof1"/"eof2"/"eof3"
    period: Tuple[int, int]    # (1985, 2014)
    season_or_freq: str        # "DJF"/"MAM"/... OR "monthly"/"yearly"


class EMOVDiagReader:
    """
    Read PCMDI PMP 'variability_modes' NetCDF outputs with a fixed schema:

      dims: time, lat, lon
      vars:
        - eof(lat, lon)        : reconstructed EOF / pattern map (teleconnection-considered)
        - slope(lat, lon)      : regression slope map
        - intercept(lat, lon)  : regression intercept map
        - pc(time)             : principal component time series (index)
        - frac                 : variance fraction (scalar, per eof)
        - lat(lat), lon(lon), time(time)

    Files live under:
      {data_dir}/{model_tag}/pcmdi_diags/model_vs_obs/metrics_data/variability_modes/{MODE}/{OBS_DATASET}/

    File name patterns (examples):
      NAM_psl_EOF1_DJF_obs_1985-2014.nc
      NAM_psl_EOF1_DJF_e3sm_*_1985-2014(_cbf).nc

      AMO_ts_EOF1_monthly_obs_1985-2014.nc
      AMO_ts_EOF1_monthly_e3sm_*_1985-2014(_cbf).nc
    """

    # fixed schema names in PMP output
    MAP_VARS = ("eof", "slope", "intercept")
    TS_VARS  = ("pc",)
    SCALARS  = ("frac",)

    def __init__(
        self,
        data_dir: str,
        *,
        prefer_cbf: bool = True,
        obs_key: str = "obs",
        model_key: str = "e3sm",
        target_lat_name: str = "latitude",
        target_lon_name: str = "longitude_a",
        obs_dataset_map: Optional[Dict[str, str]] = None,  # mode -> obs folder name (NOAA-20C, HadSST2, etc.)
        normalize_lon_360: bool = True,
        sort_latlon: bool = True,
        open_engine: Optional[str] = None,  # e.g. "netcdf4"
    ):
        self.data_dir = data_dir
        self.prefer_cbf = prefer_cbf
        self.obs_key = obs_key
        self.model_key = model_key
        self.target_lat_name = target_lat_name
        self.target_lon_name = target_lon_name
        self.obs_dataset_map = {} if obs_dataset_map is None else dict(obs_dataset_map)
        self.normalize_lon_360 = normalize_lon_360
        self.sort_latlon = sort_latlon
        self.open_engine = open_engine

    # -------------------------
    # path helpers
    # -------------------------
    def _mode_root(self, model_tag: str, mode: str) -> str:
        return os.path.join(
            self.data_dir,
            model_tag,
            "pcmdi_diags",
            "model_vs_obs",
            "metrics_data",
            "variability_modes",
            mode,
        )

    def _find_obs_dataset_dir(self, model_tag: str, mode: str) -> str:
        root = self._mode_root(model_tag, mode)
        if not os.path.isdir(root):
            raise FileNotFoundError(f"Missing mode directory: {root}")

        obs_name = self.obs_dataset_map.get(mode, None)
        if obs_name is not None:
            d = os.path.join(root, obs_name)
            if not os.path.isdir(d):
                raise FileNotFoundError(f"Obs dataset dir not found for mode={mode}: {d}")
            return d

        subdirs = sorted([d for d in glob.glob(os.path.join(root, "*")) if os.path.isdir(d)])
        if not subdirs:
            raise FileNotFoundError(f"No obs dataset subdir found under: {root}")
        return subdirs[0]

    @staticmethod
    def _EOFn(eof: str) -> str:
        eof = eof.strip().lower()
        if not eof.startswith("eof"):
            raise ValueError(f"Expected eof like 'eof1', got {eof}")
        return "EOF" + eof.replace("eof", "")

    def _obs_filename(self, spec: ModeFileSpec) -> str:
        y0, y1 = spec.period
        EOFn = self._EOFn(spec.eof)
        return f"{spec.mode}_{spec.var}_{EOFn}_{spec.season_or_freq}_{self.obs_key}_{y0}-{y1}.nc"

    def _model_glob(self, spec: ModeFileSpec) -> str:
        y0, y1 = spec.period
        EOFn = self._EOFn(spec.eof)
        return f"{spec.mode}_{spec.var}_{EOFn}_{spec.season_or_freq}_{self.model_key}_*_{y0}-{y1}*.nc"

    def _pick_best_match(self, matches: List[str]) -> str:
        matches = sorted(matches)
        if self.prefer_cbf:
            cbf = [m for m in matches if m.endswith("_cbf.nc")]
            if cbf:
                return cbf[0]
        return matches[0]

    # -------------------------
    # IO helpers
    # -------------------------
    def _open_dataset(self, path: str) -> xr.Dataset:
        if self.open_engine is None:
            return xr.open_dataset(path)
        return xr.open_dataset(path, engine=self.open_engine)

    def _rename_and_postprocess_latlon(self, da: xr.DataArray) -> xr.DataArray:
        # rename lat/lon -> target names
        rename = {}
        if "lat" in da.dims or "lat" in da.coords:
            rename["lat"] = self.target_lat_name
        if "lon" in da.dims or "lon" in da.coords:
            rename["lon"] = self.target_lon_name
        if rename:
            da = da.rename(rename)

        # normalize lon to [0, 360) if needed
        if self.normalize_lon_360 and (self.target_lon_name in da.coords):
            lon = da[self.target_lon_name]
            if float(lon.min()) < 0.0:
                da = da.assign_coords({self.target_lon_name: (lon + 360.0) % 360.0})

        # sort for stable plotting
        if self.sort_latlon:
            if self.target_lat_name in da.coords:
                da = da.sortby(self.target_lat_name)
            if self.target_lon_name in da.coords:
                da = da.sortby(self.target_lon_name)

        return da

    # -------------------------
    # read "known variables"
    # -------------------------
    def read_map(
        self,
        path: str,
        *,
        which: Literal["eof", "slope", "intercept"] = "eof",
    ) -> xr.DataArray:
        """Read a (lat,lon) map variable from a single PMP file."""
        with self._open_dataset(path) as ds:
            if which not in ds:
                raise KeyError(f"'{which}' not found in {path}. Vars={list(ds.data_vars)}")
            da = ds[which].load()
        return self._rename_and_postprocess_latlon(da)

    def read_pc(self, path: str) -> xr.DataArray:
        """Read the principal component time series pc(time) from a single PMP file."""
        with self._open_dataset(path) as ds:
            if "pc" not in ds:
                raise KeyError(f"'pc' not found in {path}. Vars={list(ds.data_vars)}")
            da = ds["pc"].load()
        return da  # pc is 1D time series; keep original time coord

    def read_frac(self, path: str) -> float:
        """Read variance fraction scalar 'frac' from a single PMP file."""
        with self._open_dataset(path) as ds:
            if "frac" not in ds:
                raise KeyError(f"'frac' not found in {path}. Vars={list(ds.data_vars)}")
            v = ds["frac"].load().values
        return float(np.asarray(v).squeeze())

    # -------------------------
    # public path-resolvers
    # -------------------------
    def obs_path(self, model_tag: str, spec: ModeFileSpec) -> str:
        """Resolve the exact obs file path for (model_tag, spec)."""
        obs_dir = self._find_obs_dataset_dir(model_tag, spec.mode)
        f = os.path.join(obs_dir, self._obs_filename(spec))
        if not os.path.exists(f):
            raise FileNotFoundError(f"Obs file not found: {f}")
        return f

    def model_path(self, model_tag: str, spec: ModeFileSpec) -> str:
        """Resolve the best-matching model file path for (model_tag, spec)."""
        obs_dir = self._find_obs_dataset_dir(model_tag, spec.mode)
        patt = os.path.join(obs_dir, self._model_glob(spec))
        matches = glob.glob(patt)
        if not matches:
            raise FileNotFoundError(f"No model file matches: {patt}")
        return self._pick_best_match(matches)

    # -------------------------
    # high-level helpers for plotting
    # -------------------------
    def read_reference_map(
        self,
        model_tag: str,
        spec: ModeFileSpec,
        *,
        which: Literal["eof", "slope", "intercept"] = "eof",
    ) -> xr.DataArray:
        return self.read_map(self.obs_path(model_tag, spec), which=which)

    def read_model_map(
        self,
        model_tag: str,
        spec: ModeFileSpec,
        *,
        which: Literal["eof", "slope", "intercept"] = "eof",
    ) -> xr.DataArray:
        return self.read_map(self.model_path(model_tag, spec), which=which)

    def build_multimodel_stack(
        self,
        model_tags: Sequence[str],
        spec: ModeFileSpec,
        *,
        which: Literal["eof", "slope", "intercept"] = "eof",
        member_dim: str = "member",
    ) -> Dict[str, xr.DataArray]:
        """
        Returns:
          {"reference": (lat,lon), "hist": (member,lat,lon)}
        where hist stacks across model_tags.
        """
        if not model_tags:
            raise ValueError("model_tags is empty")

        ref = self.read_reference_map(model_tags[0], spec, which=which)
        mods = [self.read_model_map(mt, spec, which=which) for mt in model_tags]
        ens = xr.concat(mods, dim=member_dim).assign_coords({member_dim: list(model_tags)})
        return {"reference": ref, "hist": ens}