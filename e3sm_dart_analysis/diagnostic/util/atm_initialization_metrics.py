"""Reusable DART analysis-initialization readers and diagnostics."""

from __future__ import annotations

import os, re
import numpy as np
import warnings
import pandas as pd
import math

from typing import Dict, List, Optional

from datetime import datetime
import xarray as xr
import xskillscore as xs

from xcdat.dataset import open_dataset
from xcdat.bounds import create_bounds
from xcdat.dataset import open_mfdataset
from pathlib import Path
from typing import Sequence, Optional

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

import cartopy.crs as ccrs
import cartopy.feature as cfeature

from xarray.conventions import SerializationWarning
#warnings.filterwarnings("ignore", category=SerializationWarning)
warnings.filterwarnings("once", category=SerializationWarning)
def extract_exp_info(
    data_path: str,
    *,
    resolution: str = "ne30pg2_r05_IcoswISC30E3r5",
    machine: str = "compy",
    atm_subdir: str = "archive/post/atm/180x360_aave",
    lnd_subdir: str = "archive/post/lnd/180x360_aave",
) -> Dict[str, dict]:
    """
    Build a standardized experiment metadata dictionary for E3SM ensemble and DA runs.

    Returns a dict keyed by experiment name with:
      - nens, season, group_key
      - runs: {da|fc|wc} -> sub-run dict or None
      - default_run, key, period
    """

    exps = {
        "CTRL": {
            "nens": 1,
            "key": "ctrl",
            "da_run": None,
            "fc_run": {"compset": "F20TR", "name": "CTRL", "period": "201201-201212"},
            "wc_run": {"compset": "WCYCL20TR", "name": "CTRL", "period": "201201-201212"},
        },
        "CTRL10-S0": {
            "nens": 10,
            "key": "dart_en10",
            "da_run": {
                "compset": "F20TR",
                "name": "CTRLEN10",
                "period":
                "201112-201112",
                "obs_diag": ["obs_seq", "obs_diag", "obs_common", "closest_member"],
            },
            "fc_run": {"compset": "F20TR", "name": "CTRLEN10_15day", "period": "201201-201202"},
            "wc_run": {"compset": "WCYCL20TR", "name": "CTRLEN10_15day", "period": "201201-201202"},
        },
        "CAPT10-S0": {
            "nens": 10,
            "key": "dart_en10",
            "da_run": None,
            "fc_run": {"compset": "F20TR", "name": "CAPTEN10_15day", "period": "201201-201202"},
            "wc_run": {"compset": "WCYCL20TR", "name": "CAPTEN10_15day", "period": "201201-201202"},
        },
        "DART10-S0": {
            "nens": 10,
            "key": "dart_en10",
            "da_run": {
                "compset": "F20TR", "name": "DARTEN10", "period": "201112-201112",
                "obs_diag": ["obs_seq", "obs_diag", "obs_common", "closest_member"],
            },
            "fc_run": None,
            "wc_run": None,
        },
        "DART20-S0": {
            "nens": 20,
            "key": "dart_en20",
            "da_run": {
                "compset": "F20TR", "name": "DARTEN20", "period": "201112-201112",
                "obs_diag": ["obs_seq", "obs_diag", "obs_common", "closest_member"],
            },
            "fc_run": {"compset": "F20TR", "name": "DARTEN20_15day", "period": "201201-201202"},
            "wc_run": {"compset": "WCYCL20TR", "name": "DARTEN20_15day", "period": "201201-201202"},
        },
        "DART40-S0": {
            "nens": 40,
            "key": "dart_en40",
            "da_run": {
                "compset": "F20TR", "name": "DARTEN40", "period": "201112-201112",
                "obs_diag": ["obs_seq", "obs_diag", "obs_common", "closest_member"],
            },
            "fc_run": {"compset": "F20TR", "name": "DARTEN40_15day", "period": "201201-201202"},
            "wc_run": {"compset": "WCYCL20TR", "name": "DARTEN40_15day", "period": "201201-201202"},
        },
        "CAM80-S0": {
            "nens": 80,
            "key": "dart_en80",
            "da_run": {
                "run_id": "f.e21.FHIST_BGC.f09_025.CAM6assim.011", "name": "DARTEN40", "period": "201112-201112",
                "compset": "F20TR", "obs_diag": ["obs_seq", "obs_diag", "obs_common", "closest_member"],
            },
            "fc_run": None,
            "wc_run": None,
        },
        "CTRL10-S1": {
            "nens": 10,
            "key": "ctrl_en10",
            "da_run": None,
            "fc_run": {"compset": "F20TR", "name": "CTRLEN10s1_15day", "period": "201206-201207"},
            "wc_run": {"compset": "WCYCL20TR", "name": "CTRLEN10s1_15day", "period": "201206-201207"},
        },
        "CAPT10-S1": {
            "nens": 10,
            "key": "capt_en10",
            "da_run": None,
            "fc_run": {"compset": "F20TR", "name": "CAPTEN10S1_15day", "period": "201206-201207"},
            "wc_run": {"compset": "WCYCL20TR", "name": "CAPTEN10S1_15day", "period": "201206-201207"},
        },
        "DART40-S1": {
            "nens": 40,
            "key": "dart_en40",
            "da_run": {
                "compset": "F20TR", "name": "DARTEN40S1", "period": "201205-201205",
                "obs_diag": ["obs_seq", "obs_diag", "obs_common", "closest_member"],
            },
            "fc_run": {"compset": "F20TR", "name": "DARTEN40S1_15day", "period": "201206-201207"},
            "wc_run": {"compset": "WCYCL20TR", "name": "DARTEN40S1_15day", "period": "201206-201207"},
        },
        "CAM80-S1": {
            "nens": 80,
            "key": "dart_en80",
            "da_run": {
                "run_id": "f.e22.FHIST_BGC.f09_025.CAM6assim.011", "name": "DARTEN40", "period": "201205-201205",
                "compset": "F20TR", "obs_diag": ["obs_seq", "obs_diag", "obs_common", "closest_member"],
            },
            "fc_run": None,
            "wc_run": None,
        },
    }

    def _season_from_name(name: str) -> Optional[str]:
        _SEASON_RE = re.compile(r"(?:-|_)(S\d+)\b")
        m = _SEASON_RE.search(name)
        return m.group(1) if m else None

    def _require_valid_period(p: str) -> str:
        _PERIOD_RE = re.compile(r"^\d{6}-\d{6}$")
        if not _PERIOD_RE.match(p):
            raise ValueError(f"Invalid period '{p}' for experiment; expected 'YYYYMM-YYYYMM'.")
        return p

    def _build_run(exp_name: str, spec: Optional[dict]) -> Optional[dict]:
        if spec is None:
            return None
        period = _require_valid_period(spec["period"])
        run_id = spec.get(
            "run_id",
            f"{spec['name']}_{spec['compset']}_{resolution}_{machine}",
        )
        atm_path = os.path.join(data_path, run_id, atm_subdir)
        lnd_path = os.path.join(data_path, run_id, lnd_subdir)
        out = {
            "run_id": run_id,
            "name": spec["name"],
            "compset": spec["compset"],
            "period": period,
            "atm": atm_subdir,
            "lnd": lnd_subdir,
            "atm_path": atm_path,
            "lnd_path": lnd_path,
        }
        if "obs_diag" in spec:
            out["obs_diag"] = list(spec["obs_diag"])
        return out

    exp_dict: Dict[str, dict] = {}
    for exp_name, meta in sorted(exps.items()):
        runs = {
            "da": _build_run(exp_name, meta.get("da_run")),
            "fc": _build_run(exp_name, meta.get("fc_run")),
            "wc": _build_run(exp_name, meta.get("wc_run")),
        }
        default_run = runs["fc"] or runs["wc"] or runs["da"]
        exp_dict[exp_name] = {
            "nens": meta["nens"],
            "season": _season_from_name(exp_name),
            "group_key": meta.get("key"),
            "runs": runs,
            "default_run": default_run,
            "key": meta.get("key"),
            "period": (default_run or {}).get("period"),
        }

    return exp_dict
class InvalidSpreadFile(ValueError):
    """A remapped *_sd file that was not built with the corrected spread remap."""


def check_spread_file(path):
    """Reject remapped *_sd files built before the spread-remap fix.

    Corrected files carry PS_ENSMEAN (the ensemble-mean surface pressure used
    for the vertical interpolation). Older files interpolated with the PS
    spread as surface pressure and are invalid below ~100 hPa; one CAM6 file
    (CAM6-S0 2011-12-31 18Z, no native source) even holds ensemble-mean values.
    """
    with xr.open_dataset(path) as dataset:
        if "PS_ENSMEAN" not in dataset.variables:
            raise InvalidSpreadFile(
                f"{os.path.basename(path)} lacks PS_ENSMEAN: not built with the corrected "
                "spread remap (regenerate it, or skip this cycle)"
            )




class DartAssimOutput:
    """
    Read, analyze, and plot DART data assimilation outputs for one experiment
    and one analysis time.

    This version is tailored to files named like:

        <root>/<exp>/<exp>.forecast_mean.<timestamp>.nc
        <root>/<exp>/<exp>.forecast_sd.<timestamp>.nc
        <root>/<exp>/<exp>.output_mean.<timestamp>.nc
        <root>/<exp>/<exp>.output_sd.<timestamp>.nc

    where <exp> could be 'DART40-S0', 'CTRL10-S0', 'CAM6-S0', etc., and
    <timestamp> is e.g. '2011-12-26-00000', '2011-12-26-21600', ...
    """

    def __init__(self, name, analysis_time,
                 prior_mean, prior_sd, post_mean, post_sd):
        # Basic metadata
        self.name = name
        self.analysis_time = analysis_time

        # Core datasets
        self.prior_mean = prior_mean   # forecast_mean
        self.prior_sd = prior_sd       # forecast_sd
        self.post_mean = post_mean     # output_mean
        self.post_sd = post_sd         # output_sd

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------
    @classmethod
    def from_experiment(cls, root_dir, exp_name, timestamp):
        """
        Build a DartAssimOutput from your processed files.

        root_dir:  base directory that contains the experiment folders.
        exp_name:  e.g. 'DART40-S0', 'CTRL10-S0', 'CAM6-S0'
        timestamp: e.g. '2011-12-26-00000'
        """
        exp_dir = os.path.join(root_dir, exp_name)
        prefix = os.path.join(exp_dir, exp_name)

        prior_mean_path = f"{prefix}.forecast_mean.{timestamp}.nc"
        prior_sd_path   = f"{prefix}.forecast_sd.{timestamp}.nc"
        post_mean_path  = f"{prefix}.output_mean.{timestamp}.nc"
        post_sd_path    = f"{prefix}.output_sd.{timestamp}.nc"

        if not os.path.exists(prior_mean_path):
            raise FileNotFoundError(prior_mean_path)
        if not os.path.exists(prior_sd_path):
            raise FileNotFoundError(prior_sd_path)
        if not os.path.exists(post_mean_path):
            raise FileNotFoundError(post_mean_path)
        if not os.path.exists(post_sd_path):
            raise FileNotFoundError(post_sd_path)

        prior_mean = xr.open_dataset(prior_mean_path)
        prior_sd   = cls._open_spread(prior_sd_path, prior_mean)
        post_mean  = xr.open_dataset(post_mean_path)
        post_sd    = cls._open_spread(post_sd_path, post_mean)

        return cls(exp_name, timestamp, prior_mean, prior_sd, post_mean, post_sd)

    @staticmethod
    def _open_spread(path, mean_dataset):
        """Open a remapped *_sd file; an invalid one is treated as missing.

        The returned stand-in has the mean file's structure with every
        floating-point field set to NaN, so spread-based metrics for this
        cycle are missing while mean-based metrics are still derived.
        """
        try:
            check_spread_file(path)
        except InvalidSpreadFile as error:
            print(f"[MISSING SPREAD] {error}; spread-based metrics for this cycle are NaN.")
            missing = mean_dataset.copy()
            for name, field in missing.data_vars.items():
                if np.issubdtype(field.dtype, np.floating):
                    missing[name] = field * np.nan
            missing.attrs["missing_spread"] = f"invalid source: {os.path.basename(path)}"
            return missing
        return xr.open_dataset(path)

    # ------------------------------------------------------------------
    # Small helpers for coordinates and region
    # ------------------------------------------------------------------
    @staticmethod
    def _guess_lat_name(ds_or_da):
        for d in ds_or_da.dims:
            if d.lower().startswith("lat"):
                return d
        raise KeyError("No latitude-like dimension found.")

    @staticmethod
    def _guess_lon_name(ds_or_da):
        for d in ds_or_da.dims:
            if d.lower().startswith("lon"):
                return d
        raise KeyError("No longitude-like dimension found.")

    @staticmethod
    def _guess_lev_name(ds_or_da):
        for d in ds_or_da.dims:
            if d.lower() in ("lev", "plev", "level"):
                return d
        return None

    @staticmethod
    def normalize_longitude(ds_or_da, lon_name):
        """
        Normalize longitudes to [-180, 180] and sort by longitude.
        Works for both Dataset and DataArray.
        """
        lon = ds_or_da[lon_name]
        lon_new = ((lon + 180.0) % 360.0) - 180.0  # 0–360 -> -180–180

        ds_or_da = ds_or_da.assign_coords({lon_name: lon_new})
        ds_or_da = ds_or_da.sortby(lon_name)
        return ds_or_da

    def subset_region(self, ds, lat_range=None, lon_range=None):
        """
        Subset a dataset to [lat_range, lon_range] if provided.

        If lon_range includes negative values but the dataset longitudes
        are in 0–360, normalize to [-180, 180] first so the slice behaves
        as expected.
        """
        out = ds

        # Latitude
        if lat_range is not None:
            lat_name = self._guess_lat_name(out)
            out = out.sel({lat_name: slice(lat_range[0], lat_range[1])})

        # Longitude
        if lon_range is not None:
            lon_name = self._guess_lon_name(out)

            # If user asks for negative longitudes but data are 0–360,
            # normalize to [-180, 180] first.
            if (lon_range[0] < 0 or lon_range[1] < 0) and float(out[lon_name].max()) > 180.0:
                out = self.normalize_longitude(out, lon_name)

            out = out.sel({lon_name: slice(lon_range[0], lon_range[1])})

        return out

    # ------------------------------------------------------------------
    # Bias / RMSE / spread vs reference
    # ------------------------------------------------------------------
    def compute_stats(self, ref, varname,
                      lat_range=None, lon_range=None,
                      mean_over_space=True):
        """
        Compute bias, RMSE and spread for PRIOR (forecast) and POSTERIOR (analysis)
        versus a reference dataset.

        Returns a dict of DataArray:
          bias_prior, rmse_prior, spread_prior,
          bias_post,  rmse_post,  spread_post
        """
        pri = self.subset_region(self.prior_mean[[varname]], lat_range, lon_range)
        pos = self.subset_region(self.post_mean[[varname]], lat_range, lon_range)
        spr_pri = self.subset_region(self.prior_sd[[varname]], lat_range, lon_range)
        spr_pos = self.subset_region(self.post_sd[[varname]], lat_range, lon_range)
        ref_sub = self.subset_region(ref[[varname]], lat_range, lon_range)

        pri, pos, spr_pri, spr_pos, ref_sub = xr.align(
            pri, pos, spr_pri, spr_pos, ref_sub, join="inner"
        )

        diff_prior = pri[varname] - ref_sub[varname]
        diff_post  = pos[varname] - ref_sub[varname]

        if mean_over_space:
            # Average over all dims except vertical (lev/plev/level) and time
            dims = list(diff_prior.dims)
            for cand in ("lev", "plev", "level", "time"):
                if cand in dims:
                    dims.remove(cand)
            bias_prior   = diff_prior.mean(dims)
            bias_post    = diff_post.mean(dims)
            rmse_prior   = np.sqrt((diff_prior ** 2).mean(dims))
            rmse_post    = np.sqrt((diff_post ** 2).mean(dims))
            spread_prior = spr_pri[varname].mean(dims)
            spread_post  = spr_pos[varname].mean(dims)
        else:
            bias_prior   = diff_prior
            bias_post    = diff_post
            rmse_prior   = np.sqrt(diff_prior ** 2)
            rmse_post    = np.sqrt(diff_post ** 2)
            spread_prior = spr_pri[varname]
            spread_post  = spr_pos[varname]

        return {
            "bias_prior":   bias_prior,
            "rmse_prior":   rmse_prior,
            "spread_prior": spread_prior,
            "bias_post":    bias_post,
            "rmse_post":    rmse_post,
            "spread_post":  spread_post,
        }

    # ------------------------------------------------------------------
    # Vertical profile plot (multi-panel)
    # ------------------------------------------------------------------
    @staticmethod
    def plot_profile(
        daos,
        ref_ds,
        varname="T",
        metric="rmse",
        lat_range=(-30, 30),
        lon_range=(0, 360),
        title_prefix="",
        label_map=None,
    ):
        """
        Make ONE figure with multiple panels (subplots), each showing a PRIOR vs POST
        vertical profile for one experiment.

        daos      : dict {exp_name: DartAssimOutput}
        ref_ds    : reference dataset (e.g., CAM6-S0 analysis or ERA5)
        varname   : variable to plot, e.g. 'T'
        metric    : 'rmse', 'bias', or 'spread'
        label_map : optional dict mapping exp_name -> pretty label
        """

        exp_names = list(daos.keys())
        n_exp = len(exp_names)

        # layout: up to 3 columns
        ncols = min(4, n_exp)
        nrows = math.ceil(n_exp / ncols)

        fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 6 * nrows), squeeze=False)

        for idx, exp_name in enumerate(exp_names):
            row = idx // ncols
            col = idx % ncols
            ax = axes[row, col]

            dao = daos[exp_name]
            stats = dao.compute_stats(
                ref_ds, varname,
                lat_range=lat_range,
                lon_range=lon_range
            )

            # find vertical coord
            lev_name = None
            for d in stats["rmse_prior"].dims:
                if d.lower() in ("lev", "plev", "level"):
                    lev_name = d
            if lev_name is None:
                raise ValueError(f"No vertical coordinate found for {exp_name}.")

            lev = stats["rmse_prior"][lev_name]

            # helper: reduce to 1D(lev)
            def _to_1d(arr):
                dims = [d for d in arr.dims if d != lev_name]
                if dims:
                    arr = arr.mean(dims)
                return arr

            # prior/post curves (same metric)
            if metric == "rmse":
                y_pri = _to_1d(stats["rmse_prior"])
                y_pos = _to_1d(stats["rmse_post"])
                xlabel = f"{varname} RMSE"
            elif metric == "bias":
                y_pri = _to_1d(stats["bias_prior"])
                y_pos = _to_1d(stats["bias_post"])
                xlabel = f"{varname} Bias"
            elif metric == "spread":
                y_pri = _to_1d(stats["spread_prior"])
                y_pos = _to_1d(stats["spread_post"])
                xlabel = f"{varname} Spread"
            else:
                raise ValueError("metric must be 'rmse', 'bias', or 'spread'")

            ax.plot(y_pri, lev, label="Prior", marker="o")
            ax.plot(y_pos, lev, label="Post", marker="s")

            ax.invert_yaxis()
            ax.set_xlabel(xlabel)
            ax.set_ylabel(lev_name)

            label = label_map.get(exp_name, exp_name) if label_map else exp_name
            ax.set_title(label)
            ax.grid(True, linestyle="--", alpha=0.4)
            ax.legend(fontsize=8)

        # turn off unused axes if any
        for k in range(n_exp, nrows * ncols):
            row = k // ncols
            col = k % ncols
            axes[row, col].axis("off")

        if title_prefix:
            fig.suptitle(title_prefix, fontsize=14)
            fig.tight_layout(rect=[0, 0, 1, 0.95])
        else:
            fig.tight_layout()

        return fig

    # ------------------------------------------------------------------
    # Lat–lon map (multi-panel, with Cartopy + longitude normalization)
    # ------------------------------------------------------------------
    @staticmethod
    def plot_map(
        daos,
        ref_ds,
        varname="PS",
        metric="bias",      # 'bias', 'rmse', or 'increment'
        level=None,
        lat_range=(-30, 30),
        lon_range=(-180, 180),
        figsize=(8,24),
        title_prefix="",
        label_map=None,
    ):
        """
        ONE figure with multiple panels, each showing a lat–lon map for one experiment.

        metric:
          'increment' : post - prior
          'bias'      : post - ref
          'rmse'      : RMSE(post vs ref)
        """

        exp_names = list(daos.keys())
        n_exp = len(exp_names)

        # Layout: force 2x2 for 4 exps, otherwise behave sensibly
        if n_exp <= 2:
            ncols, nrows = n_exp, 1
        elif n_exp <= 4:
            ncols, nrows = 2, 2
        else:
            ncols = 3
            nrows = math.ceil(n_exp / ncols)

        # First pass: compute fields for all experiments to get shared vmin/vmax
        fields = {}
        vmin = None
        vmax = None

        for exp_name in exp_names:
            dao = daos[exp_name]
            ds_pri = dao.subset_region(dao.prior_mean[[varname]], lat_range, lon_range)
            ds_pos = dao.subset_region(dao.post_mean[[varname]], lat_range, lon_range)

            lev_name = dao._guess_lev_name(ds_pos)
            idx = None
            if level is not None and lev_name is not None and lev_name in ds_pos.dims:
                idx = int(np.argmin(np.abs(ds_pos[lev_name] - level)))
                ds_pos = ds_pos.isel({lev_name: idx})
                ds_pri = ds_pri.isel({lev_name: idx})

            # normalize longitude to [-180, 180] if needed
            try:
                lon_name = dao._guess_lon_name(ds_pos)
                if float(ds_pos[lon_name].max()) > 180.0:
                    ds_pos = dao.normalize_longitude(ds_pos, lon_name)
                    ds_pri = dao.normalize_longitude(ds_pri, lon_name)
            except KeyError:
                lon_name = None  # no lon dimension

            if metric == "increment":
                field = ds_pos[varname] - ds_pri[varname]
            else:
                ref_sub = dao.subset_region(ref_ds[[varname]], lat_range, lon_range)
                if level is not None and lev_name is not None and lev_name in ref_sub.dims:
                    ref_sub = ref_sub.isel({lev_name: idx})

                # normalize reference longitudes the same way
                if lon_name is not None and float(ref_sub[lon_name].max()) > 180.0:
                    ref_sub = dao.normalize_longitude(ref_sub, lon_name)

                ds_pos_aligned, ref_sub = xr.align(ds_pos, ref_sub, join="inner")
                if metric == "bias":
                    field = ds_pos_aligned[varname] - ref_sub[varname]
                elif metric == "rmse":
                    field = np.sqrt((ds_pos_aligned[varname] - ref_sub[varname]) ** 2)
                else:
                    raise ValueError("metric must be 'bias', 'rmse', or 'increment'")

            # ensure 2D (lat, lon)
            lat_name = dao._guess_lat_name(field)
            lon_name = dao._guess_lon_name(field)
            extras = [d for d in field.dims if d not in (lat_name, lon_name)]
            for d in extras:
                field = field.isel({d: 0})
            field = field.squeeze(drop=True)

            fields[exp_name] = (field, lat_name, lon_name)

            fmin = float(field.min())
            fmax = float(field.max())
            vmin = fmin if vmin is None else min(vmin, fmin)
            vmax = fmax if vmax is None else max(vmax, fmax)

        # symmetric color scale often nicer for bias/increment
        if metric in ("bias", "increment"):
            m = max(abs(vmin), abs(vmax))
            vmin, vmax = -m, m

        proj = ccrs.PlateCarree(central_longitude=0)

        # --- Second pass: plot with Cartopy ---
        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=figsize,
            subplot_kw=dict(projection=proj),
            constrained_layout=True,   # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            squeeze=False,
        )

        im_last = None  # for shared colorbar

        for idx, exp_name in enumerate(exp_names):
            row = idx // ncols
            col = idx % ncols
            ax = axes[row, col]

            field, lat_name, lon_name = fields[exp_name]
            lon_vals = field[lon_name].values
            lat_vals = field[lat_name].values
            lon2d, lat2d = np.meshgrid(lon_vals, lat_vals)

            # Plot
            im = ax.pcolormesh(
                lon2d,
                lat2d,
                field.values,
                shading="auto",
                transform=proj,
                vmin=vmin,
                vmax=vmax,
            )
            im_last = im

            # gridlines & labels: only on left column / bottom row to reduce clutter
            gl = ax.gridlines(
                draw_labels=True,
                linewidth=0.5,
                color="gray",
                alpha=0.5,
                linestyle="--",
            )
            gl.top_labels = False
            gl.right_labels = False
            if row > 0:
                gl.bottom_labels = False
            if col > 0:
                gl.left_labels = False
            gl.xlabel_style = {"size": 8}
            gl.ylabel_style = {"size": 8}

            ax.coastlines(color="black", linewidth=0.7)

            # focus on requested region
            ax.set_extent(
                [lon_range[0], lon_range[1],
                 lat_range[0], lat_range[1]],
                crs=proj,
            )

            label = label_map.get(exp_name, exp_name) if label_map else exp_name
            ax.set_title(label, fontsize=10)

        # turn off unused axes
        for k in range(len(exp_names), nrows * ncols):
            r = k // ncols
            c = k % ncols
            axes[r, c].axis("off")

        # --- shared colorbar on the right ---
        fig.subplots_adjust(right=0.88)  # a bit more room
        cax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
        cbar = fig.colorbar(im_last, cax=cax)
        cbar.ax.tick_params(labelsize=9)

        if title_prefix:
            fig.suptitle(title_prefix, fontsize=14)
            fig.tight_layout(rect=[0, 0, 0.88, 0.93])
        else:
            fig.tight_layout(rect=[0, 0, 0.88, 1])

        return fig



# ---- Cached diagnostic products ----

def profile_cache_path(cache_dir, experiment, timestamp, variable, metric):
    """Return the request-keyed profile cache path."""
    cache_dir = Path(cache_dir)
    return cache_dir / f"{experiment}_{timestamp}_profile_{variable}_{metric}.nc"


def map_cache_path(cache_dir, experiment, timestamp, variable, metric, level=None):
    """Return the request-keyed map cache path."""
    cache_dir = Path(cache_dir)
    level_tag = "surface" if level is None else f"level-{level:g}"
    return cache_dir / f"{experiment}_{timestamp}_map_{variable}_{metric}_{level_tag}.nc"


def _atomic_to_netcdf(dataset, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    dataset.to_netcdf(temporary)
    os.replace(temporary, path)


def derive_profile_metric(dao, reference, variable, metric, lat_range, lon_range):
    """Derive prior/posterior profile data for one variable and metric."""
    if metric not in {"bias", "rmse", "spread"}:
        raise ValueError("Profile metric must be bias, rmse, or spread.")
    stats = dao.compute_stats(
        reference,
        variable,
        lat_range=lat_range,
        lon_range=lon_range,
        mean_over_space=True,
    )
    dataset = xr.Dataset(
        {
            "prior": stats[f"{metric}_prior"],
            "posterior": stats[f"{metric}_post"],
        }
    )
    dataset.attrs.update(
        {
            "experiment": dao.name,
            "analysis_time": dao.analysis_time,
            "variable": variable,
            "metric": metric,
            "product_type": "profile",
        }
    )
    return dataset


def derive_map_metric(
    dao,
    reference,
    variable,
    metric,
    lat_range,
    lon_range,
    level=None,
):
    """Derive a two-dimensional initialization metric map."""
    if metric not in {"bias", "rmse", "spread", "increment"}:
        raise ValueError("Map metric must be bias, rmse, spread, or increment.")

    prior = dao.subset_region(
        dao.prior_mean[[variable]], lat_range, lon_range
    )
    posterior = dao.subset_region(
        dao.post_mean[[variable]], lat_range, lon_range
    )
    posterior_spread = dao.subset_region(
        dao.post_sd[[variable]], lat_range, lon_range
    )
    level_name = dao._guess_lev_name(posterior)
    level_index = None
    if level is not None and level_name is not None and level_name in posterior.dims:
        level_values = np.asarray(posterior[level_name].values, dtype=float)
        level_units = str(posterior[level_name].attrs.get("units", "")).lower()
        pressure_in_pa = (
            level_units in {"pa", "pascal", "pascals"}
            or (level_name.lower() == "plev" and np.nanmax(np.abs(level_values)) > 2000)
        )
        selection_values = level_values / 100.0 if pressure_in_pa else level_values
        level_index = int(np.abs(selection_values - level).argmin())
        prior = prior.isel({level_name: level_index})
        posterior = posterior.isel({level_name: level_index})
        posterior_spread = posterior_spread.isel({level_name: level_index})

    longitude_name = dao._guess_lon_name(posterior)
    if float(posterior[longitude_name].max()) > 180.0:
        prior = dao.normalize_longitude(prior, longitude_name)
        posterior = dao.normalize_longitude(posterior, longitude_name)
        posterior_spread = dao.normalize_longitude(
            posterior_spread, longitude_name
        )

    if metric == "increment":
        field = posterior[variable] - prior[variable]
    elif metric == "spread":
        field = posterior_spread[variable]
    else:
        ref = dao.subset_region(reference[[variable]], lat_range, lon_range)
        if level_index is not None and level_name in ref.dims:
            ref = ref.isel({level_name: level_index})
        if float(ref[longitude_name].max()) > 180.0:
            ref = dao.normalize_longitude(ref, longitude_name)
        posterior, ref = xr.align(posterior, ref, join="inner")
        difference = posterior[variable] - ref[variable]
        field = difference if metric == "bias" else np.sqrt(difference ** 2)

    latitude_name = dao._guess_lat_name(field)
    longitude_name = dao._guess_lon_name(field)
    for dimension in list(field.dims):
        if dimension not in (latitude_name, longitude_name):
            field = field.isel({dimension: 0})
    field = field.squeeze(drop=True)
    field.name = "metric"
    dataset = field.to_dataset()
    dataset.attrs.update(
        {
            "experiment": dao.name,
            "analysis_time": dao.analysis_time,
            "variable": variable,
            "metric": metric,
            "product_type": "map",
            "level": "surface" if level is None else float(level),
        }
    )
    return dataset


def cache_profile_metric(
    path,
    dao,
    reference,
    variable,
    metric,
    lat_range,
    lon_range,
    force_compute=False,
):
    """Load a profile cache, deriving it only when missing or forced."""
    path = Path(path)
    if path.exists() and not force_compute:
        return xr.load_dataset(path)
    dataset = derive_profile_metric(
        dao, reference, variable, metric, lat_range, lon_range
    )
    _atomic_to_netcdf(dataset, path)
    return dataset


def cache_map_metric(
    path,
    dao,
    reference,
    variable,
    metric,
    lat_range,
    lon_range,
    level=None,
    force_compute=False,
):
    """Load a map cache, deriving it only when missing or forced."""
    path = Path(path)
    if path.exists() and not force_compute:
        return xr.load_dataset(path)
    dataset = derive_map_metric(
        dao,
        reference,
        variable,
        metric,
        lat_range,
        lon_range,
        level=level,
    )
    _atomic_to_netcdf(dataset, path)
    return dataset


def plot_cached_profiles(
    datasets,
    variable,
    metric,
    label_map=None,
    title="",
):
    """Plot cached profile products for multiple experiments."""
    experiment_names = list(datasets)
    columns = min(4, len(experiment_names))
    rows = math.ceil(len(experiment_names) / columns)
    fig, axes = plt.subplots(
        rows,
        columns,
        figsize=(4 * columns, 6 * rows),
        squeeze=False,
    )
    for index, experiment in enumerate(experiment_names):
        axis = axes[index // columns, index % columns]
        dataset = datasets[experiment]
        prior = dataset["prior"]
        posterior = dataset["posterior"]
        level_name = next(
            (
                dimension
                for dimension in prior.dims
                if dimension.lower() in ("lev", "plev", "level")
            ),
            None,
        )
        if level_name is None:
            raise ValueError(f"No vertical coordinate found for {experiment}.")
        other_dimensions = [
            dimension for dimension in prior.dims if dimension != level_name
        ]
        if other_dimensions:
            prior = prior.mean(other_dimensions)
            posterior = posterior.mean(other_dimensions)
        level_values = prior[level_name]
        axis.plot(prior, level_values, label="Prior", marker="o")
        axis.plot(posterior, level_values, label="Posterior", marker="s")
        axis.invert_yaxis()
        axis.set_xlabel(f"{variable} {metric.upper()}")
        axis.set_ylabel(level_name)
        axis.set_title(
            label_map.get(experiment, experiment) if label_map else experiment
        )
        axis.grid(True, linestyle="--", alpha=0.4)
        axis.legend(fontsize=8)

    for index in range(len(experiment_names), rows * columns):
        axes[index // columns, index % columns].axis("off")
    if title:
        fig.suptitle(title, fontsize=14)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
    else:
        fig.tight_layout()
    return fig


def plot_cached_maps(
    datasets,
    variable,
    metric,
    lat_range,
    lon_range,
    label_map=None,
    title="",
    figsize=(8, 24),
):
    """Plot cached map products with a shared color scale."""
    experiment_names = list(datasets)
    if len(experiment_names) <= 2:
        columns, rows = len(experiment_names), 1
    elif len(experiment_names) <= 4:
        columns, rows = 2, 2
    else:
        columns, rows = 3, math.ceil(len(experiment_names) / 3)

    fields = {name: dataset["metric"] for name, dataset in datasets.items()}
    minimum = min(float(field.min()) for field in fields.values())
    maximum = max(float(field.max()) for field in fields.values())
    if metric in {"bias", "increment"}:
        magnitude = max(abs(minimum), abs(maximum))
        minimum, maximum = -magnitude, magnitude

    projection = ccrs.PlateCarree()
    fig, axes = plt.subplots(
        rows,
        columns,
        figsize=figsize,
        subplot_kw={"projection": projection},
        squeeze=False,
    )
    image = None
    for index, experiment in enumerate(experiment_names):
        axis = axes[index // columns, index % columns]
        field = fields[experiment]
        latitude_name = next(
            dim for dim in field.dims if dim.lower().startswith("lat")
        )
        longitude_name = next(
            dim for dim in field.dims if dim.lower().startswith("lon")
        )
        longitude, latitude = np.meshgrid(
            field[longitude_name].values,
            field[latitude_name].values,
        )
        image = axis.pcolormesh(
            longitude,
            latitude,
            field.values,
            shading="auto",
            transform=projection,
            vmin=minimum,
            vmax=maximum,
        )
        axis.coastlines(linewidth=0.7)
        axis.set_extent(
            [lon_range[0], lon_range[1], lat_range[0], lat_range[1]],
            crs=projection,
        )
        axis.set_title(
            label_map.get(experiment, experiment) if label_map else experiment
        )

    for index in range(len(experiment_names), rows * columns):
        axes[index // columns, index % columns].axis("off")
    if image is not None:
        fig.colorbar(image, ax=axes.ravel().tolist(), shrink=0.75)
    if title:
        fig.suptitle(title, fontsize=14)
    return fig


def _area_weighted_mean(field, latitude_name):
    """Cosine-latitude-weighted mean over finite grid cells."""
    weights = np.cos(np.deg2rad(field[latitude_name])).broadcast_like(field)
    valid = np.isfinite(field)
    total = float(weights.where(valid).sum())
    if total <= 0.0:
        return np.nan
    return float((field * weights).where(valid).sum()) / total


def plot_combined_cached_maps(
    datasets,
    requests,
    experiments,
    lat_range,
    lon_range,
    label_map=None,
    title="",
    figsize=(16, 15),
    font_size=18,
    vmax_percentile=99.0,
    cmap="viridis",
    missing_color="0.85",
    show_area_mean=True,
):
    """Plot variable rows by experiment columns with one color scale per row.

    ``datasets`` is keyed by ``(variable, level)`` and then experiment. Each
    dataset must contain a two-dimensional ``metric`` field. Plot-only scale
    factors and units are supplied by the corresponding request dictionaries.

    Each row's upper color limit is the ``vmax_percentile`` of all finite
    values in that row (``None`` uses the maximum); larger values saturate and
    the colorbar is extended. Missing cells, such as pressure surfaces below
    terrain, are drawn in ``missing_color``. When ``show_area_mean`` is true,
    each panel is annotated with its cosine-latitude-weighted mean.
    """
    if not requests:
        raise ValueError("At least one combined-map request is required.")
    if not experiments:
        raise ValueError("At least one experiment is required.")

    projection = ccrs.PlateCarree()
    figure, axes = plt.subplots(
        len(requests),
        len(experiments),
        figsize=figsize,
        subplot_kw={"projection": projection},
        squeeze=False,
        layout="constrained",
    )
    colormap = plt.get_cmap(cmap)

    panel_index = 0
    for row, request in enumerate(requests):
        variable = request["variable"]
        level = request.get("level")
        request_key = (variable, level)
        if request_key not in datasets:
            raise KeyError(f"Missing combined-map request {request_key}.")
        scale = float(request.get("scale", 1.0))
        fields = {
            experiment: datasets[request_key][experiment]["metric"] * scale
            for experiment in experiments
        }
        row_values = np.concatenate([
            np.asarray(field.values, dtype=float).ravel()
            for field in fields.values()
        ])
        row_values = row_values[np.isfinite(row_values)]
        if row_values.size:
            row_maximum = float(row_values.max())
            maximum = (
                row_maximum if vmax_percentile is None
                else float(np.percentile(row_values, vmax_percentile))
            )
        else:
            row_maximum = maximum = 1.0
        if maximum <= 0.0:
            maximum = row_maximum if row_maximum > 0.0 else 1.0

        row_image = None
        for column, experiment in enumerate(experiments):
            axis = axes[row, column]
            axis.set_facecolor(missing_color)
            field = fields[experiment].squeeze(drop=True)
            latitude_name = next(
                dim for dim in field.dims if dim.lower().startswith("lat")
            )
            longitude_name = next(
                dim for dim in field.dims if dim.lower().startswith("lon")
            )
            longitude, latitude = np.meshgrid(
                field[longitude_name].values,
                field[latitude_name].values,
            )
            row_image = axis.pcolormesh(
                longitude,
                latitude,
                field.values,
                shading="auto",
                transform=projection,
                cmap=colormap,
                vmin=0.0,
                vmax=maximum,
                rasterized=True,
            )
            axis.coastlines(linewidth=0.7)
            axis.set_extent(
                [lon_range[0], lon_range[1], lat_range[0], lat_range[1]],
                crs=projection,
            )
            gridlines = axis.gridlines(
                draw_labels=True,
                linewidth=0.4,
                color="0.5",
                alpha=0.4,
                linestyle="--",
                x_inline=False,
                y_inline=False,
            )
            gridlines.top_labels = False
            gridlines.right_labels = False
            gridlines.bottom_labels = row == len(requests) - 1
            gridlines.left_labels = column == 0
            gridlines.xlabel_style = {"size": 0.65 * font_size}
            gridlines.ylabel_style = {"size": 0.65 * font_size}

            panel_letter = chr(ord("a") + panel_index)
            panel_index += 1
            axis.set_title(
                f"({panel_letter})",
                loc="left",
                fontweight="normal",
                fontsize=0.78 * font_size,
            )
            if show_area_mean:
                area_mean = _area_weighted_mean(field, latitude_name)
                axis.set_title(
                    f"mean {area_mean:.2f}",
                    loc="right",
                    fontweight="normal",
                    fontsize=0.70 * font_size,
                )
            if row == 0:
                experiment_label = (
                    label_map.get(experiment, experiment)
                    if label_map else experiment
                )
                axis.annotate(
                    experiment_label,
                    xy=(0.5, 1.0),
                    xycoords="axes fraction",
                    xytext=(0, 1.6 * font_size),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=0.90 * font_size,
                    fontweight="bold",
                )
            if column == 0:
                axis.annotate(
                    request.get("label", variable),
                    xy=(0.0, 0.5),
                    xycoords="axes fraction",
                    xytext=(-3.2 * font_size, 0),
                    textcoords="offset points",
                    ha="right",
                    va="center",
                    rotation=90,
                    fontsize=0.90 * font_size,
                    fontweight="bold",
                )

        if row_image is not None:
            unit = request.get("unit", "")
            colorbar = figure.colorbar(
                row_image,
                ax=axes[row, :].tolist(),
                orientation="vertical",
                shrink=0.85,
                pad=0.015,
                extend="max" if row_maximum > maximum else "neither",
            )
            if unit:
                colorbar.set_label(unit, fontsize=0.72 * font_size)
            colorbar.ax.tick_params(labelsize=0.65 * font_size)

    if title:
        figure.suptitle(title, fontsize=font_size)
    return figure



def cross_section_cache_path(
    cache_dir, experiment, timestamp, variable, metric, horizontal_axis
):
    """Return the request-keyed vertical cross-section cache path."""
    if horizontal_axis not in {"latitude", "longitude"}:
        raise ValueError("horizontal_axis must be 'latitude' or 'longitude'.")
    cache_dir = Path(cache_dir)
    return cache_dir / (
        f"{experiment}_{timestamp}_cross-section_{horizontal_axis}_"
        f"{variable}_{metric}.nc"
    )


def derive_cross_section_metric(
    dao,
    reference,
    variable,
    metric,
    lat_range,
    lon_range,
    horizontal_axis="latitude",
):
    """Derive a level-by-latitude or level-by-longitude state metric."""
    if metric not in {"bias", "rmse", "spread"}:
        raise ValueError("Cross-section metric must be bias, rmse, or spread.")
    if horizontal_axis not in {"latitude", "longitude"}:
        raise ValueError("horizontal_axis must be 'latitude' or 'longitude'.")

    prior = dao.subset_region(dao.prior_mean[[variable]], lat_range, lon_range)
    posterior = dao.subset_region(dao.post_mean[[variable]], lat_range, lon_range)
    prior_spread = dao.subset_region(dao.prior_sd[[variable]], lat_range, lon_range)
    posterior_spread = dao.subset_region(
        dao.post_sd[[variable]], lat_range, lon_range
    )
    reference = dao.subset_region(reference[[variable]], lat_range, lon_range)
    prior, posterior, prior_spread, posterior_spread, reference = xr.align(
        prior, posterior, prior_spread, posterior_spread, reference, join="inner"
    )

    level_name = dao._guess_lev_name(posterior)
    if level_name is None:
        raise ValueError(f"{variable} has no vertical coordinate.")
    latitude_name = dao._guess_lat_name(posterior)
    longitude_name = dao._guess_lon_name(posterior)
    retained_name = latitude_name if horizontal_axis == "latitude" else longitude_name
    reduced_name = longitude_name if horizontal_axis == "latitude" else latitude_name

    if horizontal_axis == "longitude" and float(posterior[longitude_name].max()) > 180.0:
        prior = dao.normalize_longitude(prior, longitude_name)
        posterior = dao.normalize_longitude(posterior, longitude_name)
        prior_spread = dao.normalize_longitude(prior_spread, longitude_name)
        posterior_spread = dao.normalize_longitude(posterior_spread, longitude_name)
        reference = dao.normalize_longitude(reference, longitude_name)

    if metric == "spread":
        prior_field = prior_spread[variable].mean(reduced_name)
        posterior_field = posterior_spread[variable].mean(reduced_name)
    else:
        prior_difference = prior[variable] - reference[variable]
        posterior_difference = posterior[variable] - reference[variable]
        if metric == "bias":
            prior_field = prior_difference.mean(reduced_name)
            posterior_field = posterior_difference.mean(reduced_name)
        else:
            prior_field = np.sqrt((prior_difference ** 2).mean(reduced_name))
            posterior_field = np.sqrt((posterior_difference ** 2).mean(reduced_name))

    def _to_section(field):
        for dimension in list(field.dims):
            if dimension not in (level_name, retained_name):
                field = field.isel({dimension: 0})
        return field.squeeze(drop=True).transpose(level_name, retained_name)

    dataset = xr.Dataset(
        {
            "prior": _to_section(prior_field),
            "posterior": _to_section(posterior_field),
        }
    )
    dataset.attrs.update(
        {
            "experiment": dao.name,
            "analysis_time": dao.analysis_time,
            "variable": variable,
            "metric": metric,
            "product_type": "vertical_cross_section",
            "horizontal_axis": horizontal_axis,
        }
    )
    return dataset


def cache_cross_section_metric(
    path,
    dao,
    reference,
    variable,
    metric,
    lat_range,
    lon_range,
    horizontal_axis="latitude",
    force_compute=False,
):
    """Load a cross-section cache, deriving it only when missing or forced."""
    path = Path(path)
    if path.exists() and not force_compute:
        return xr.load_dataset(path)
    dataset = derive_cross_section_metric(
        dao,
        reference,
        variable,
        metric,
        lat_range,
        lon_range,
        horizontal_axis=horizontal_axis,
    )
    _atomic_to_netcdf(dataset, path)
    return dataset


def plot_cached_cross_sections(
    datasets,
    variable,
    metric,
    label_map=None,
    title="",
):
    """Plot cached prior/posterior vertical sections for each experiment."""
    experiment_names = list(datasets)
    fields = [
        datasets[name][state]
        for name in experiment_names
        for state in ("prior", "posterior")
    ]
    minimum = min(float(field.min()) for field in fields)
    maximum = max(float(field.max()) for field in fields)
    if metric == "bias":
        magnitude = max(abs(minimum), abs(maximum))
        minimum, maximum = -magnitude, magnitude

    fig, axes = plt.subplots(
        len(experiment_names),
        2,
        figsize=(12, 4 * len(experiment_names)),
        squeeze=False,
    )
    image = None
    for row, experiment in enumerate(experiment_names):
        dataset = datasets[experiment]
        for column, state in enumerate(("prior", "posterior")):
            axis = axes[row, column]
            field = dataset[state]
            level_name = next(
                dimension
                for dimension in field.dims
                if dimension.lower() in ("lev", "plev", "level")
            )
            horizontal_name = next(
                dimension for dimension in field.dims if dimension != level_name
            )
            image = axis.pcolormesh(
                field[horizontal_name],
                field[level_name],
                field,
                shading="auto",
                vmin=minimum,
                vmax=maximum,
            )
            axis.invert_yaxis()
            axis.set_xlabel(horizontal_name)
            axis.set_ylabel(level_name)
            label = label_map.get(experiment, experiment) if label_map else experiment
            axis.set_title(f"{label} - {state.title()}")
    if image is not None:
        fig.colorbar(image, ax=axes.ravel().tolist(), shrink=0.8)
    if title:
        fig.suptitle(title, fontsize=14)
        fig.subplots_adjust(top=0.94)
    return fig


def increment_cross_section_cache_path(
    cache_dir, experiment, timestamp, variable, horizontal_axis, increment_mode
):
    """Return the request-keyed increment/spread-ratio cross-section cache path."""
    if horizontal_axis not in {"latitude", "longitude"}:
        raise ValueError("horizontal_axis must be latitude or longitude.")
    if increment_mode not in {"signed", "absolute"}:
        raise ValueError("increment_mode must be signed or absolute.")
    return Path(cache_dir) / (
        f"{experiment}_{timestamp}_cross-section_{horizontal_axis}_"
        f"{variable}_{increment_mode}-increment-spread-ratio-hpa.nc"
    )


def derive_increment_cross_section(
    dao, variable, lat_range, lon_range, horizontal_axis="latitude",
    increment_mode="absolute",
):
    """Derive spatial sections of increment and posterior/prior spread ratio."""
    if horizontal_axis not in {"latitude", "longitude"}:
        raise ValueError("horizontal_axis must be latitude or longitude.")
    if increment_mode not in {"signed", "absolute"}:
        raise ValueError("increment_mode must be signed or absolute.")

    prior = dao.subset_region(dao.prior_mean[[variable]], lat_range, lon_range)
    posterior = dao.subset_region(dao.post_mean[[variable]], lat_range, lon_range)
    prior_spread = dao.subset_region(
        dao.prior_sd[[variable]], lat_range, lon_range
    )
    posterior_spread = dao.subset_region(
        dao.post_sd[[variable]], lat_range, lon_range
    )
    prior, posterior, prior_spread, posterior_spread = xr.align(
        prior, posterior, prior_spread, posterior_spread, join="inner"
    )

    level_name = dao._guess_lev_name(posterior)
    if level_name is None:
        raise ValueError(f"{variable} has no vertical coordinate.")
    latitude_name = dao._guess_lat_name(posterior)
    longitude_name = dao._guess_lon_name(posterior)

    if float(posterior[longitude_name].max()) > 180.0:
        prior = dao.normalize_longitude(prior, longitude_name)
        posterior = dao.normalize_longitude(posterior, longitude_name)
        prior_spread = dao.normalize_longitude(prior_spread, longitude_name)
        posterior_spread = dao.normalize_longitude(
            posterior_spread, longitude_name
        )

    increment = posterior[variable] - prior[variable]
    if increment_mode == "absolute":
        increment = np.abs(increment)
    spread_ratio = xr.where(
        np.isfinite(prior_spread[variable]) & (prior_spread[variable] > 0),
        posterior_spread[variable] / prior_spread[variable],
        np.nan,
    )

    retained_name = (
        latitude_name if horizontal_axis == "latitude" else longitude_name
    )
    reduced_name = (
        longitude_name if horizontal_axis == "latitude" else latitude_name
    )

    def _horizontal_mean(field, dimension):
        if dimension == latitude_name:
            latitude_weights = np.cos(np.deg2rad(field[latitude_name])).clip(min=0.0)
            return field.weighted(latitude_weights).mean(dimension, skipna=True)
        return field.mean(dimension, skipna=True)

    def _to_section(field):
        field = _horizontal_mean(field, reduced_name)
        extra_dimensions = [
            dimension for dimension in field.dims
            if dimension not in (level_name, retained_name)
        ]
        if extra_dimensions:
            field = field.mean(extra_dimensions, skipna=True)
        return field.transpose(level_name, retained_name)

    dataset = xr.Dataset(
        {
            "increment": _to_section(increment),
            "prior_spread": _to_section(prior_spread[variable]),
            "posterior_spread": _to_section(posterior_spread[variable]),
            "spread_ratio": _to_section(spread_ratio),
        }
    )
    level_values = np.asarray(dataset[level_name].values, dtype=float)
    level_units = str(dataset[level_name].attrs.get("units", "")).lower()
    pressure_in_pa = (
        level_units in {"pa", "pascal", "pascals"}
        or (level_name.lower() == "plev" and np.nanmax(np.abs(level_values)) > 2000)
    )
    if pressure_in_pa:
        dataset = dataset.assign_coords({level_name: dataset[level_name] / 100.0})
    if level_name.lower() == "plev" or pressure_in_pa:
        dataset[level_name].attrs.update({
            "long_name": "Pressure",
            "units": "hPa",
        })
    dataset.attrs.update(
        {
            "experiment": dao.name,
            "analysis_time": dao.analysis_time,
            "variable": variable,
            "metric": "increment",
            "increment_mode": increment_mode,
            "product_type": "vertical_increment_cross_section",
            "horizontal_axis": horizontal_axis,
            "horizontal_mean": (
                "cosine-latitude weighted" if reduced_name == latitude_name
                else "arithmetic longitude mean"
            ),
        }
    )
    return dataset


def cache_increment_cross_section(
    path, dao, variable, lat_range, lon_range, horizontal_axis="latitude",
    increment_mode="absolute", force_compute=False,
):
    """Load an increment section cache, deriving it only when required."""
    path = Path(path)
    if path.exists() and not force_compute:
        return xr.load_dataset(path)
    dataset = derive_increment_cross_section(
        dao, variable, lat_range, lon_range,
        horizontal_axis=horizontal_axis, increment_mode=increment_mode,
    )
    _atomic_to_netcdf(dataset, path)
    return dataset


def increment_zonal_mean_square_cache_path(
    cache_dir, experiment, timestamp, variable
):
    """Return the per-cycle zonal-mean squared-increment cache path."""
    return Path(cache_dir) / (
        f"{experiment}_{timestamp}_zonal-section_{variable}_"
        "increment-mean-square-hpa.nc"
    )


def derive_increment_zonal_mean_square(dao, variable, lat_range, lon_range):
    """Derive the zonal mean of the squared analysis increment.

    Returns a ``metric`` field on (pressure in hPa, latitude). Cells where the
    pressure surface lies below terrain are missing and excluded from the
    zonal mean; ``valid_fraction`` records the fraction of longitudes used.
    The square is kept (not its root) so cycles can be averaged before the
    window RMS is taken.
    """
    prior = dao.subset_region(dao.prior_mean[[variable]], lat_range, lon_range)
    posterior = dao.subset_region(dao.post_mean[[variable]], lat_range, lon_range)
    prior, posterior = xr.align(prior, posterior, join="inner")

    level_name = dao._guess_lev_name(posterior)
    if level_name is None:
        raise ValueError(f"{variable} has no vertical coordinate.")
    latitude_name = dao._guess_lat_name(posterior)
    longitude_name = dao._guess_lon_name(posterior)

    squared = (posterior[variable] - prior[variable]) ** 2
    extra_dimensions = [
        dimension for dimension in squared.dims
        if dimension not in (level_name, latitude_name, longitude_name)
    ]
    if extra_dimensions:
        squared = squared.mean(extra_dimensions, skipna=True)
    valid_fraction = np.isfinite(squared).mean(longitude_name)
    zonal_mean = squared.mean(longitude_name, skipna=True)

    dataset = xr.Dataset({
        "metric": zonal_mean.transpose(level_name, latitude_name),
        "valid_fraction": valid_fraction.transpose(level_name, latitude_name),
    })
    level_values = np.asarray(dataset[level_name].values, dtype=float)
    level_units = str(dataset[level_name].attrs.get("units", "")).lower()
    if level_units in {"pa", "pascal", "pascals"} or np.nanmax(level_values) > 2000:
        dataset = dataset.assign_coords({level_name: dataset[level_name] / 100.0})
    dataset[level_name].attrs.update({"long_name": "Pressure", "units": "hPa"})
    dataset.attrs.update({
        "experiment": dao.name,
        "analysis_time": dao.analysis_time,
        "variable": variable,
        "metric": "increment_zonal_mean_square",
        "product_type": "zonal_mean_section",
    })
    return dataset


def cache_increment_zonal_mean_square(
    path, dao, variable, lat_range, lon_range, force_compute=False,
):
    """Load a zonal-mean squared-increment cache, deriving it when required."""
    path = Path(path)
    if path.exists() and not force_compute:
        return xr.load_dataset(path)
    dataset = derive_increment_zonal_mean_square(
        dao, variable, lat_range, lon_range
    )
    _atomic_to_netcdf(dataset, path)
    return dataset


def _global_rms_from_map(rms_field, latitude_name):
    """Cosine-latitude-weighted global RMS of a gridpoint RMS map."""
    return float(np.sqrt(_area_weighted_mean(rms_field ** 2, latitude_name)))


def _global_rms_from_section(mean_square, level_name, latitude_name):
    """Mass-weighted (cos-latitude x layer thickness) RMS of a zonal section."""
    pressure = np.asarray(mean_square[level_name].values, dtype=float)
    thickness = np.abs(np.gradient(pressure))
    weights = (
        xr.DataArray(thickness, dims=[level_name], coords={level_name: pressure})
        * np.cos(np.deg2rad(mean_square[latitude_name]))
    ).broadcast_like(mean_square)
    valid = np.isfinite(mean_square)
    total = float(weights.where(valid).sum())
    if total <= 0.0:
        return np.nan
    return float(np.sqrt(float((mean_square * weights).where(valid).sum()) / total))


SECTION_PRESSURE_TICKS_HPA = (1000, 700, 500,300, 200, 100, 50, 30, 20, 10, 5, 2, 1)


def _format_pressure_axis(axis, bottom, top, log_pressure=True):
    """Orient a pressure (hPa) axis from ``bottom`` (surface) up to ``top``."""
    if log_pressure:
        axis.set_yscale("log")
        axis.set_yticks([
            tick for tick in SECTION_PRESSURE_TICKS_HPA
            if top <= tick <= bottom
        ])
        axis.yaxis.set_major_formatter(plt.FuncFormatter(lambda value, _: f"{value:g}"))
        axis.yaxis.set_minor_locator(plt.NullLocator())
    axis.set_ylim(bottom, top)


def plot_map_and_zonal_sections(
    map_datasets,
    map_request,
    section_datasets,
    section_requests,
    experiments,
    lat_range,
    lon_range,
    label_map=None,
    title="",
    figsize=(16, 17),
    font_size=18,
    pressure_top=100.0,
    log_pressure=False,
    vmax_percentile=99.0,
    cmap="YlGnBu",
    missing_color="0.85",
    color_bins=10,
    show_global_rms=True,
    gap_before=None,
    gap_width=0.12,
):
    """Plot an RMS map row above zonal-mean RMS section rows.

    ``map_datasets`` is keyed by experiment and holds a gridpoint RMS
    ``metric`` map. ``section_datasets`` is keyed by variable and then
    experiment and holds a window-mean squared increment ``metric`` on
    (pressure, latitude); its square root is plotted. Each row shares one
    discrete color scale of about ``color_bins`` round-valued intervals
    whose top covers the ``vmax_percentile`` of that row; larger values
    saturate and the colorbar is extended.
    With ``show_global_rms``, panels are annotated with the area-weighted
    (map) or mass-weighted (section, below ``pressure_top``) global RMS.
    With ``log_pressure``, section rows use a logarithmic pressure axis.
    ``gap_before`` (a column index) inserts an empty spacer column of relative
    width ``gap_width`` before that column, e.g. to set apart another period.
    """
    if not experiments:
        raise ValueError("At least one experiment is required.")
    colormap = plt.get_cmap(cmap)
    projection = ccrs.PlateCarree()
    n_rows = 1 + len(section_requests)
    n_columns = len(experiments)
    figure = plt.figure(figsize=figsize, layout="constrained")
    width_ratios = [1.0] * n_columns
    if gap_before is not None:
        if not 0 < gap_before < n_columns:
            raise ValueError("gap_before must index an interior column.")
        width_ratios.insert(gap_before, gap_width)
    grid = figure.add_gridspec(
        n_rows, len(width_ratios),
        height_ratios=[1.0] + [1.05] * len(section_requests),
        width_ratios=width_ratios,
    )

    def _grid_column(column):
        """Gridspec column for a panel column, skipping the spacer."""
        return column + (gap_before is not None and column >= gap_before)

    def _row_scale(values):
        """Return round color levels, a norm, and the colorbar extend mode."""
        values = np.concatenate([np.asarray(v, dtype=float).ravel() for v in values])
        values = values[np.isfinite(values)]
        row_maximum = float(values.max()) if values.size else 1.0
        maximum = (
            row_maximum if vmax_percentile is None or not values.size
            else float(np.percentile(values, vmax_percentile))
        )
        if maximum <= 0.0:
            maximum = row_maximum if row_maximum > 0.0 else 1.0
        levels = MaxNLocator(nbins=color_bins).tick_values(0.0, maximum)
        levels = levels[levels >= 0.0]
        if levels[-1] < maximum:
            levels = np.append(levels, 2 * levels[-1] - levels[-2])
        extend = "max" if row_maximum > levels[-1] else "neither"
        norm = BoundaryNorm(levels, colormap.N, extend=extend)
        return levels, norm, extend

    def _annotate(axis, row, column, row_label, value, unit=""):
        """Titles "(a) EXPERIMENT" / "(VAR)"; global RMS boxed inside the panel."""
        panel_letter = chr(ord("a") + row * n_columns + column)
        experiment = experiments[column]
        experiment_label = (
            label_map.get(experiment, experiment) if label_map else experiment
        )
        axis.set_title(
            f"({panel_letter}) {experiment_label}",
            loc="left", fontsize=0.72 * font_size,
        )
        axis.set_title(f"({row_label})", loc="right", fontsize=0.72 * font_size)
        if show_global_rms and np.isfinite(value):
            # Map: bottom-left corner; sections: upper left (quiet
            # stratosphere, away from the near-surface Q signal).
            x, y, ha, va = (
                (0.02, 0.06, "left", "bottom") if row == 0
                else (0.02, 0.96, "left", "top")
            )
            axis.text(
                x, y,
                f"Mean: {value:.2f}" + (f" {unit}" if unit else ""),
                transform=axis.transAxes, ha=ha, va=va,
                fontsize=0.62 * font_size, zorder=10,
                bbox={
                    "boxstyle": "round,pad=0.25", "facecolor": "none",
                    "edgecolor": "0.6", "linewidth": 0.6,
                },
            )

    # Row 0: horizontal RMS map.
    scale = float(map_request.get("scale", 1.0))
    map_fields = {
        experiment: (map_datasets[experiment]["metric"] * scale).squeeze(drop=True)
        for experiment in experiments
    }
    levels, norm, extend = _row_scale(f.values for f in map_fields.values())
    map_axes = []
    image = None
    for column, experiment in enumerate(experiments):
        axis = figure.add_subplot(grid[0, _grid_column(column)], projection=projection)
        map_axes.append(axis)
        axis.set_facecolor(missing_color)
        field = map_fields[experiment]
        latitude_name = next(d for d in field.dims if d.lower().startswith("lat"))
        longitude_name = next(d for d in field.dims if d.lower().startswith("lon"))
        image = axis.pcolormesh(
            field[longitude_name].values, field[latitude_name].values,
            field.values, shading="auto", transform=projection,
            cmap=colormap, norm=norm, rasterized=True,
        )
        axis.coastlines(linewidth=0.7)
        axis.set_extent(
            [lon_range[0], lon_range[1], lat_range[0], lat_range[1]],
            crs=projection,
        )
        gridlines = axis.gridlines(
            draw_labels=True, linewidth=0.4, color="0.5", alpha=0.4,
            linestyle="--", x_inline=False, y_inline=False,
        )
        gridlines.top_labels = False
        gridlines.right_labels = False
        gridlines.left_labels = column == 0
        gridlines.xlabel_style = {"size": 0.65 * font_size}
        gridlines.ylabel_style = {"size": 0.65 * font_size}
        _annotate(
            axis, 0, column, map_request.get("label", map_request["variable"]),
            _global_rms_from_map(field, latitude_name),
            map_request.get("unit", ""),
        )
    colorbar = figure.colorbar(
        image, ax=map_axes, orientation="vertical", shrink=0.85, pad=0.015,
        extend=extend, ticks=levels[::2],
    )
    colorbar.set_label(map_request.get("unit", ""), fontsize=0.72 * font_size)
    colorbar.ax.tick_params(labelsize=0.65 * font_size)

    # Rows 1..: zonal-mean RMS sections.
    shared_axis = None
    for row, request in enumerate(section_requests, start=1):
        variable = request["variable"]
        if variable not in section_datasets:
            raise KeyError(f"Missing section request {variable}.")
        scale = float(request.get("scale", 1.0))
        fields = {}
        for experiment in experiments:
            mean_square = section_datasets[variable][experiment]["metric"]
            level_name = next(
                d for d in mean_square.dims if d.lower() in ("lev", "plev", "level")
            )
            mean_square = mean_square.where(
                mean_square[level_name] >= pressure_top, drop=True
            ) * scale ** 2
            fields[experiment] = (mean_square, level_name)
        levels, norm, extend = _row_scale(
            np.sqrt(ms.values) for ms, _ in fields.values()
        )
        row_axes = []
        image = None
        for column, experiment in enumerate(experiments):
            axis = figure.add_subplot(
                grid[row, _grid_column(column)],
                sharex=shared_axis, sharey=shared_axis,
            )
            shared_axis = shared_axis or axis
            row_axes.append(axis)
            axis.set_facecolor(missing_color)
            mean_square, level_name = fields[experiment]
            latitude_name = next(
                d for d in mean_square.dims if d.lower().startswith("lat")
            )
            rms = np.sqrt(mean_square)
            image = axis.contourf(
                rms[latitude_name].values, rms[level_name].values, rms.values,
                levels=levels, cmap=colormap, norm=norm, extend=extend,
            )
            axis.contour(
                rms[latitude_name].values, rms[level_name].values, rms.values,
                levels=levels[2::2], colors="k", linewidths=0.4, alpha=0.5,
            )
            axis.set_xlim(lat_range[0], lat_range[1])
            _format_pressure_axis(
                axis, float(rms[level_name].max()), pressure_top, log_pressure
            )
            axis.set_xticks(np.arange(-90, 91, 30))
            axis.xaxis.set_major_formatter(plt.FuncFormatter(
                lambda value, _: (
                    "EQ" if value == 0
                    else f"{abs(value):.0f}°{'N' if value > 0 else 'S'}"
                )
            ))
            axis.grid(color="0.5", alpha=0.3, linestyle="--", linewidth=0.4)
            axis.tick_params(labelsize=0.65 * font_size)
            if column == 0:
                axis.set_ylabel("Pressure (hPa)", fontsize=0.72 * font_size)
            else:
                axis.tick_params(labelleft=False)
            if row != n_rows - 1:
                axis.tick_params(labelbottom=False)
            _annotate(
                axis, row, column, request.get("label", variable),
                _global_rms_from_section(mean_square, level_name, latitude_name),
                request.get("unit", ""),
            )
        colorbar = figure.colorbar(
            image, ax=row_axes, orientation="vertical", shrink=0.85, pad=0.015,
            ticks=levels[::2],
        )
        colorbar.set_label(request.get("unit", ""), fontsize=0.72 * font_size)
        colorbar.ax.tick_params(labelsize=0.65 * font_size)

    if title:
        figure.suptitle(title, fontsize=font_size)
    return figure


def plot_cached_increment_cross_sections(
    datasets, variable, increment_mode="absolute", label_map=None, title="",
    pressure_top=None, log_pressure=False,
):
    """Plot spatial sections of increment and posterior/prior spread ratio.

    ``pressure_top`` (hPa) drops levels above it before color limits are set;
    ``log_pressure`` uses a logarithmic pressure axis.
    """
    if increment_mode not in {"signed", "absolute"}:
        raise ValueError("increment_mode must be signed or absolute.")
    experiment_names = list(datasets)
    if not experiment_names:
        raise ValueError("datasets is empty; nothing to plot.")
    if pressure_top is not None:
        trimmed = {}
        for name in experiment_names:
            dataset = datasets[name]
            level_name = next(
                dimension for dimension in dataset["increment"].dims
                if dimension.lower() in ("lev", "plev", "level")
            )
            trimmed[name] = dataset.where(
                dataset[level_name] >= pressure_top, drop=True
            )
        datasets = trimmed

    fields_by_column = {
        key: [datasets[name][key] for name in experiment_names]
        for key in ("increment", "spread_ratio")
    }
    limits = {}
    for key, fields in fields_by_column.items():
        minimum = min(float(field.min(skipna=True)) for field in fields)
        maximum = max(float(field.max(skipna=True)) for field in fields)
        if key == "increment" and increment_mode == "signed":
            magnitude = max(abs(minimum), abs(maximum))
            minimum, maximum = -magnitude, magnitude
        elif key == "spread_ratio":
            magnitude = max(abs(minimum - 1.0), abs(maximum - 1.0))
            if np.isclose(magnitude, 0.0):
                magnitude = 0.05
            minimum, maximum = max(0.0, 1.0 - magnitude), 1.0 + magnitude
        elif minimum >= 0:
            minimum = 0.0
        if np.isclose(minimum, maximum):
            maximum = minimum + 1.0
        limits[key] = (minimum, maximum)

    fig, axes = plt.subplots(
        len(experiment_names), 2,
        figsize=(12, 4 * len(experiment_names)), squeeze=False,
    )
    images = [None, None]
    for row, experiment in enumerate(experiment_names):
        dataset = datasets[experiment]
        for column, key in enumerate(("increment", "spread_ratio")):
            axis = axes[row, column]
            field = dataset[key]
            level_name = next(
                dimension for dimension in field.dims
                if dimension.lower() in ("lev", "plev", "level")
            )
            horizontal_name = next(
                dimension for dimension in field.dims if dimension != level_name
            )
            images[column] = axis.pcolormesh(
                field[horizontal_name], field[level_name], field, shading="auto",
                vmin=limits[key][0], vmax=limits[key][1],
                cmap=(
                    "RdBu_r"
                    if key == "spread_ratio"
                    or (key == "increment" and increment_mode == "signed")
                    else None
                ),
            )
            levels_hpa = field[level_name].values
            _format_pressure_axis(
                axis, float(np.max(levels_hpa)), float(np.min(levels_hpa)),
                log_pressure,
            )
            axis.set_xlabel(horizontal_name)
            level_units = field[level_name].attrs.get("units", "")
            axis.set_ylabel(
                "Pressure (hPa)" if level_units == "hPa" else level_name
            )
            label = label_map.get(experiment, experiment) if label_map else experiment
            heading = (
                "Increment" if key == "increment"
                else "Posterior/prior spread ratio"
            )
            axis.set_title(f"{label} - {heading}")

    for column, image in enumerate(images):
        fig.colorbar(image, ax=axes[:, column].tolist(), shrink=0.8)
    if title:
        fig.suptitle(title, fontsize=14)
        fig.subplots_adjust(top=0.94)
    return fig
