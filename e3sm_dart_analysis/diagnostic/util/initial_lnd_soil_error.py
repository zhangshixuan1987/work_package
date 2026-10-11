"""Initial-condition soil-moisture error against observations.

``SoilBiasAssessorFromCollectedFiles`` computes, for one target date, the
ensemble mean, spread, bias, t-test significance, pattern correlation and
RMSE of each experiment against an observed soil-moisture product, from the
file lists of :mod:`util.initial_lnd_file_collection`.
``SoilErrorMapPlotter`` draws two experiments' bias and spread maps (2x2)
from the saved product.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.colors import BoundaryNorm
from scipy.ndimage import gaussian_filter
from scipy.stats import ttest_1samp


@dataclass(frozen=True)
class SoilBiasFromFilesConfig:
    out_dir: str | Path

    # Variable names
    model_var_name: str = "SOILWATER_10CM"  # model variable (already kg/m2)
    obs_var_name: str = "H2OSOI"            # obs variable (m3/m3)
    out_var_name: str = "SOILWATER_10CM"    # for metadata / file naming

    # Unit conversion (obs m3/m3 -> kg/m2 via rho*depth)
    depth_cm: float = 10.0
    water_density: float = 1000.0

    # If model is already an integrated product, no DZ integration needed
    model_is_integrated: bool = True

    # Only used if model_is_integrated=False
    dz_var: str = "DZSOI"
    lev_dim: str = "levgrnd"
    soilayer_file: str | None = None

    # Masking
    landmask_file: str | None = None
    landmask_var: str = "landfrac"
    landfrac_threshold: float = 0.1
    mask_land: bool = True

    # Significance
    alpha: float = 0.05
    significance_threshold: float | None = None

    # NEW: t-test robustness
    sig_min_n: int = 5          # require at least this many finite ensemble samples per grid cell
    sig_eps_std: float = 1e-12  # treat near-zero variance cells as not-testable

    # Coords
    lat_name: str = "lat"
    lon_name: str = "lon"

    # NEW: time selection robustness
    nearest_tolerance: str = "12H"  # tolerance used for method="nearest" time selection


class SoilBiasAssessorFromCollectedFiles:
    """
    Input expected (for one var):
      var_files[exp]["obs"]   -> list[str] (obs files)
      var_files[exp]["model"] -> dict[ens]-> list[str] (model files)
    """

    def __init__(self, cfg: SoilBiasFromFilesConfig):
        self.cfg = cfg
        self.out_dir = Path(cfg.out_dir)
        self.out_dir.mkdir(parents=True, exist_ok=True)

        self.landmask = None
        if cfg.landmask_file and os.path.exists(cfg.landmask_file):
            ds_mask = xr.open_dataset(cfg.landmask_file)
            self.landmask = ds_mask[cfg.landmask_var]
            self.landmask = self._to_0_360(self.landmask)

        self.fixed_dz = None
        if (not cfg.model_is_integrated) and cfg.soilayer_file and os.path.exists(cfg.soilayer_file):
            self.fixed_dz = xr.open_dataset(cfg.soilayer_file)[cfg.dz_var]

    def compute_one_date_from_varfiles(
        self,
        var_files: dict,
        *,
        target_date: str | pd.Timestamp,
        obs_key: str,
    ) -> xr.Dataset:
        t = pd.to_datetime(target_date)

        dsets = []
        for exp, info in var_files.items():
            obs_files = info["obs"]
            model_by_ens = info["model"]

            obs_field = self._load_obs_snapshot_from_files(obs_files, t)
            ens_stack = self._load_model_ensemble_snapshot_from_files(model_by_ens, t)

            ens_stack, obs_aligned = self._align_ens_to_obs(ens_stack, obs_field)

            model_mean = ens_stack.mean("ensemble")
            model_std = ens_stack.std("ensemble")
            bias = model_mean - obs_aligned

            sig_mask = self._compute_significance_mask(
                ens_stack=ens_stack,
                obs_field=obs_aligned,
                alpha=self.cfg.alpha,
                threshold=self.cfg.significance_threshold,
            )

            pcc, rmse = self._pcc_rmse(model_mean, obs_aligned)

            ds_exp = xr.Dataset(
                data_vars=dict(
                    obs=obs_aligned,
                    model_mean=model_mean,
                    model_std=model_std,
                    bias=bias,
                    sig_mask=sig_mask,
                    pcc=xr.DataArray(pcc, dims=()),
                    rmse=xr.DataArray(rmse, dims=()),
                ),
                coords={
                    self.cfg.lat_name: obs_aligned[self.cfg.lat_name],
                    self.cfg.lon_name: obs_aligned[self.cfg.lon_name],
                },
                attrs=dict(
                    target_date=str(t.date()),
                    exp=exp,
                    obs=obs_key,
                    model_var=self.cfg.model_var_name,
                    obs_var=self.cfg.obs_var_name,
                    out_var=self.cfg.out_var_name,
                    depth_cm=float(self.cfg.depth_cm),
                    units="kg/m2",
                    alpha=float(self.cfg.alpha),
                    sig_min_n=int(self.cfg.sig_min_n),
                ),
            ).expand_dims(exp=[exp])

            dsets.append(ds_exp)

        ds = xr.concat(dsets, dim="exp").expand_dims(time=[np.datetime64(t.to_datetime64())])
        return ds

    def save(self, ds: xr.Dataset, out_nc: str | Path, overwrite: bool = False) -> Path:
        out_nc = Path(out_nc)
        if out_nc.exists() and (not overwrite):
            raise FileExistsError(f"Output exists: {out_nc}")
        ds.to_netcdf(out_nc)
        return out_nc

    # ---------- helpers ----------

    def _open_mfd(self, files: list[str]) -> xr.Dataset:
        missing = [f for f in files if not os.path.exists(f)]
        if missing:
            raise FileNotFoundError(f"Missing files (showing up to 5): {missing[:5]}")
        return xr.open_mfdataset(files, combine="by_coords")

    def _to_0_360(self, da_or_ds):
        lon_name = self.cfg.lon_name
        if lon_name not in da_or_ds.coords:
            return da_or_ds
        lon = da_or_ds[lon_name]
        if lon.min() < 0 or lon.max() > 360:
            new_lon = (lon % 360).sortby(lon)
            da_or_ds = da_or_ds.assign_coords({lon_name: new_lon}).sortby(lon_name)
        return da_or_ds

    def _sel_time_robust(self, ds: xr.Dataset, t: pd.Timestamp) -> xr.Dataset:
        """
        Prefer exact match; fall back to nearest within tolerance.
        Helps when times are 00:00 vs 12:00 or slightly shifted.
        """
        if "time" not in ds.coords and "time" not in ds.dims:
            return ds
        try:
            return ds.sel(time=t)
        except Exception:
            tol = np.timedelta64(pd.to_timedelta(self.cfg.nearest_tolerance))
            return ds.sel(time=t, method="nearest", tolerance=tol)

    def _integrate_top_layers(self, ds: xr.Dataset) -> xr.DataArray:
        # Only used if model_is_integrated=False
        sm = ds[self.cfg.model_var_name]

        if self.cfg.dz_var in ds:
            dz = ds[self.cfg.dz_var]
            if "time" in dz.dims:
                dz = dz.isel(time=0)
        elif self.fixed_dz is not None:
            dz = self.fixed_dz
        else:
            raise ValueError("No DZSOI in dataset and no external soilayer_file provided.")

        lev = self.cfg.lev_dim
        dz_cumsum = dz.cumsum(dim=lev)
        top_mask = dz_cumsum <= (self.cfg.depth_cm / 100.0)

        sm_top = sm.where(top_mask)
        dz_masked = dz.where(top_mask)

        sm_mass = sm_top * dz_masked * self.cfg.water_density
        integrated = sm_mass.sum(dim=lev)
        integrated.attrs["units"] = "kg/m2"

        if self.cfg.mask_land and (self.landmask is not None):
            integrated = self._to_0_360(integrated)
            integrated = integrated.where(self.landmask > self.cfg.landfrac_threshold)

        return integrated

    def _load_model_ensemble_snapshot_from_files(
        self,
        model_by_ens: dict[str, list[str]],
        target_date: pd.Timestamp,
    ) -> xr.DataArray:
        ens_fields = []
        ens_labels = []
        for ens, files in model_by_ens.items():
            ds = self._open_mfd(files)
            ds = self._to_0_360(ds)
            ds_day = self._sel_time_robust(ds, target_date)

            if self.cfg.model_is_integrated:
                mod = ds_day[self.cfg.model_var_name]
            else:
                mod = self._integrate_top_layers(ds_day)

            mod = self._to_0_360(mod).squeeze()

            if self.cfg.mask_land and (self.landmask is not None):
                mod = mod.where(self.landmask > self.cfg.landfrac_threshold)

            ens_fields.append(mod)
            ens_labels.append(ens)

        out = xr.concat(ens_fields, dim="ensemble")
        out = out.assign_coords(ensemble=np.array(ens_labels, dtype=str))
        return out

    def _load_obs_snapshot_from_files(self, obs_files: list[str], target_date: pd.Timestamp) -> xr.DataArray:
        ds = self._open_mfd(obs_files)
        ds_day = self._sel_time_robust(ds, target_date)

        obs = ds_day[self.cfg.obs_var_name]

        # mask missing
        fill = obs.attrs.get("_FillValue", None)
        if fill is None:
            fill = obs.attrs.get("missing_value", -9999)
        obs = obs.where(obs != fill)

        obs = self._to_0_360(obs)

        # Convert volumetric (m3/m3) -> kg/m2 using effective depth
        obs = obs * self.cfg.water_density * (self.cfg.depth_cm / 100.0)
        obs.attrs["units"] = "kg/m2"

        if self.cfg.mask_land and (self.landmask is not None):
            obs = obs.where(self.landmask > self.cfg.landfrac_threshold)

        return obs

    def _align_ens_to_obs(self, ens_stack: xr.DataArray, obs: xr.DataArray):
        if (ens_stack.sizes.get(self.cfg.lat_name) != obs.sizes.get(self.cfg.lat_name)) or \
           (ens_stack.sizes.get(self.cfg.lon_name) != obs.sizes.get(self.cfg.lon_name)):
            ens_stack = ens_stack.interp_like(obs)
        return ens_stack, obs

    def _compute_significance_mask(
        self,
        ens_stack: xr.DataArray,
        obs_field: xr.DataArray,
        alpha: float,
        threshold: float | None,
    ) -> xr.DataArray:
        """
        Robust 1-sample t-test on ensemble differences (ens - obs):
          - require at least cfg.sig_min_n finite ensemble samples per grid cell
          - skip cells with ~zero variance (cfg.sig_eps_std)
          - mark not-testable cells as NaN in sig_mask
        """
        diff = ens_stack - obs_field
        arr = diff.transpose("ensemble", self.cfg.lat_name, self.cfg.lon_name).values

        n_valid = np.sum(np.isfinite(arr), axis=0)
        std = np.nanstd(arr, axis=0)

        pvals = np.ones_like(std, dtype=float)  # default: not significant
        ok = (n_valid >= self.cfg.sig_min_n) & (std > self.cfg.sig_eps_std)

        if np.any(ok):
            # compute once, then apply ok mask
            _, p_all = ttest_1samp(arr, popmean=0.0, axis=0, nan_policy="omit")
            pvals[ok] = p_all[ok]

        sig = xr.DataArray(
            (pvals < alpha),
            coords={
                self.cfg.lat_name: obs_field[self.cfg.lat_name],
                self.cfg.lon_name: obs_field[self.cfg.lon_name],
            },
            dims=(self.cfg.lat_name, self.cfg.lon_name),
        )

        # not-testable cells -> NaN (cleaner than warnings/garbage)
        sig = sig.where(ok)

        if threshold is not None:
            mean_bias = (ens_stack.mean("ensemble") - obs_field)
            sig = sig.where(np.abs(mean_bias) > threshold)

        return sig

    def _pcc_rmse(self, model_mean: xr.DataArray, obs: xr.DataArray) -> tuple[float, float]:
        mv = model_mean.values.ravel()
        ov = obs.values.ravel()
        valid = np.isfinite(mv) & np.isfinite(ov)
        if valid.sum() < 2:
            return np.nan, np.nan
        pcc = float(np.corrcoef(mv[valid], ov[valid])[0, 1])
        rmse = float(np.sqrt(((mv[valid] - ov[valid]) ** 2).mean()))
        return pcc, rmse


@dataclass(frozen=True)
class SoilErrorMapPlotConfig:
    # data coords
    lat_name: str = "lat"
    lon_name: str = "lon"
    exp_dim: str = "exp"
    time_dim: str = "time"

    # variable names in file
    bias_var: str = "bias"
    spread_var: str = "model_std"   # "Spread" panel
    sigmask_var: str = "sig_mask"   # optional
    
    smooth: bool = False
    smooth_sigma_bias: float = 0.8     # ~1 grid cell
    smooth_sigma_spread: float = 0.6

    # map projection
    projection: ccrs.CRS = ccrs.Robinson()
    data_crs: ccrs.CRS = ccrs.PlateCarree()

    # figure layout
    figsize: tuple[float, float] = (11.5, 7.8)
    wspace: float = 0.08
    hspace: float = 0.10
    title_fontsize: int = 18

    # coastlines
    coast_lw: float = 1.2
    borders_lw: float = 0.6

    # Bias colorbar (diverging)
    bias_levels: tuple[float, ...] = (-8, -6, -4, -2, -1, 0, 1, 2, 4, 6, 8)
    bias_cmap: str = "BrBG_r"

    # Spread colorbar (sequential)
    spread_levels: tuple[float, ...] = (0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50)
    spread_cmap: str = "YlGnBu"

    # Caption control (default OFF)
    add_caption: bool = False
    caption_fontsize: int = 16

    # Stippling
    stipple: bool = False
    stipple_stride: int = 3
    stipple_size: float = 2.0
    stipple_color: str = "0.3"
    stipple_alpha: float = 0.6

    # Rendering tweaks
    shading: str = "nearest"     # "nearest" helps avoid mesh artifacts
    rasterized: bool = True      # good for pdf + reduces screen-door look


class SoilErrorMapPlotter:
    def __init__(self, cfg: SoilErrorMapPlotConfig):
        self.cfg = cfg

    # -------------------------
    # public API
    # -------------------------
    def plot_2x2(
        self,
        *,
        nc_file: str | Path,
        exp_left: str,
        exp_right: str,
        time_index: int = 0,
        obs_label: str | None = None,   # if None -> ds.attrs["obs"] if present
        caption: str | None = None,     # only used if cfg.add_caption=True
        out_png: str | Path | None = None,
        show: bool = True,
        dpi: int = 200,
    ):
        ds = xr.open_dataset(nc_file)

        # normalize lon to 0..360 for clean global plotting
        ds = self._to_0_360(ds, lon_name=self.cfg.lon_name)

        lat = ds[self.cfg.lat_name]
        lon = ds[self.cfg.lon_name]

        # pick time
        if self.cfg.time_dim in ds.dims:
            ds_t = ds.isel({self.cfg.time_dim: time_index})
        else:
            ds_t = ds

        # obs label
        if obs_label is None:
            obs_label = ds.attrs.get("obs", "OBS")

        # extract fields
        bL, sL, pccL, rmseL, sigL = self._get_one(ds_t, exp_left)
        bR, sR, pccR, rmseR, sigR = self._get_one(ds_t, exp_right)

        # norms
        bias_levels = np.asarray(self.cfg.bias_levels, dtype=float)
        bias_norm = BoundaryNorm(bias_levels, ncolors=plt.get_cmap(self.cfg.bias_cmap).N, clip=False)

        spread_levels = np.asarray(self.cfg.spread_levels, dtype=float)
        spread_norm = BoundaryNorm(spread_levels, ncolors=plt.get_cmap(self.cfg.spread_cmap).N, clip=False)

        # ---- figure ----
        fig = plt.figure(figsize=self.cfg.figsize)
        gs = fig.add_gridspec(
            nrows=2, ncols=3,
            width_ratios=[1.0, 1.0, 0.06],
            wspace=self.cfg.wspace,
            hspace=self.cfg.hspace,
        )

        ax_bL = fig.add_subplot(gs[0, 0], projection=self.cfg.projection)
        ax_bR = fig.add_subplot(gs[0, 1], projection=self.cfg.projection)
        ax_sL = fig.add_subplot(gs[1, 0], projection=self.cfg.projection)
        ax_sR = fig.add_subplot(gs[1, 1], projection=self.cfg.projection)

        cax_bias = fig.add_subplot(gs[0, 2])
        cax_spre = fig.add_subplot(gs[1, 2])

        # ---- plot panels ----
        m_bL = self._plot_map(ax_bL, lon, lat, bL,
                              cmap=self.cfg.bias_cmap, norm=bias_norm, is_bias=True)
        m_bR = self._plot_map(ax_bR, lon, lat, bR,
                              cmap=self.cfg.bias_cmap, norm=bias_norm, is_bias=True)
        
        m_sL = self._plot_map(ax_sL, lon, lat, sL,
                              cmap=self.cfg.spread_cmap, norm=spread_norm, is_bias=False)
        m_sR = self._plot_map(ax_sR, lon, lat, sR,
                              cmap=self.cfg.spread_cmap, norm=spread_norm, is_bias=False)

        # titles
        ax_bL.set_title(f"{exp_left} - {obs_label} (Bias)", fontsize=self.cfg.title_fontsize, pad=10)
        ax_bR.set_title(f"{exp_right} - {obs_label} (Bias)", fontsize=self.cfg.title_fontsize, pad=10)
        ax_sL.set_title(f"{exp_left} (Spread)", fontsize=self.cfg.title_fontsize, pad=10)
        ax_sR.set_title(f"{exp_right} (Spread)", fontsize=self.cfg.title_fontsize, pad=10)

        # PCC/RMSE boxes on bias panels
        self._stats_box(ax_bL, pccL, rmseL)
        self._stats_box(ax_bR, pccR, rmseR)

        # stippling (ONLY over valid-data land points)
        if self.cfg.stipple:
            if sigL is not None:
                self._stipple(ax_bL, lon, lat, sigL, bL)
            if sigR is not None:
                self._stipple(ax_bR, lon, lat, sigR, bR)

        # ---- colorbars ----
        cb1 = plt.colorbar(m_bR, cax=cax_bias, orientation="vertical", extend="both")
        cb1.set_label("Bias (kg m$^{-2}$)", fontsize=self.cfg.title_fontsize)
        cb1.set_ticks(list(self.cfg.bias_levels))
        cb1.ax.tick_params(labelsize=self.cfg.title_fontsize - 2)

        cb2 = plt.colorbar(m_sR, cax=cax_spre, orientation="vertical", extend="max")
        cb2.set_label("Spread (kg m$^{-2}$)", fontsize=self.cfg.title_fontsize)
        cb2.set_ticks(list(self.cfg.spread_levels))
        cb2.ax.tick_params(labelsize=self.cfg.title_fontsize - 2)

        # caption: OFF by default; only if enabled explicitly
        if self.cfg.add_caption and caption:
            fig.text(0.5, 0.02, caption, ha="center", va="bottom", fontsize=self.cfg.caption_fontsize)

        if out_png is not None:
            out_png = Path(out_png)
            out_png.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
        if show:
            plt.show()
        else:
            plt.close(fig)

        return fig

    # -------------------------
    # internals
    # -------------------------
    def _get_one(self, ds_t: xr.Dataset, exp: str):
        exp_dim = self.cfg.exp_dim
        if exp_dim in ds_t.dims:
            one = ds_t.sel({exp_dim: exp})
        else:
            one = ds_t

        bias = one[self.cfg.bias_var]
        spread = one[self.cfg.spread_var]

        pcc = float(one["pcc"].values) if "pcc" in one else np.nan
        rmse = float(one["rmse"].values) if "rmse" in one else np.nan

        sig = None
        if self.cfg.sigmask_var in one.data_vars:
            sig = one[self.cfg.sigmask_var]

        return bias, spread, pcc, rmse, sig

    def _plot_map(self, ax, lon, lat, field, *, cmap, norm, is_bias: bool):
        ax.set_global()
        ax.add_feature(cfeature.COASTLINE, linewidth=self.cfg.coast_lw)
        ax.add_feature(cfeature.BORDERS, linewidth=self.cfg.borders_lw)
    
        data = np.asarray(field.values)
    
        # ---- optional smoothing (visualization only) ----
        if self.cfg.smooth:
            sigma = (
                self.cfg.smooth_sigma_bias
                if is_bias
                else self.cfg.smooth_sigma_spread
            )
    
            # preserve NaNs during smoothing
            mask = np.isfinite(data)
            data_filled = np.where(mask, data, 0.0)
    
            smooth_data = gaussian_filter(data_filled, sigma=sigma, mode="nearest")
            smooth_mask = gaussian_filter(mask.astype(float), sigma=sigma, mode="nearest")
    
            with np.errstate(invalid="ignore"):
                data = np.where(smooth_mask > 0.5, smooth_data / smooth_mask, np.nan)
    
        m = ax.pcolormesh(
            lon, lat, data,
            transform=self.cfg.data_crs,
            cmap=plt.get_cmap(cmap),
            norm=norm,
            shading=self.cfg.shading,
            rasterized=self.cfg.rasterized,
        )
        return m


    def _stats_box(self, ax, pcc: float, rmse: float):
        txt = f"PCC: {pcc:.2f}\nRMSE: {rmse:.2f}"
        ax.text(
            0.98, 0.05, txt,
            transform=ax.transAxes,
            ha="right", va="bottom",
            fontsize=self.cfg.title_fontsize - 2,
            bbox=dict(boxstyle="square,pad=0.25", facecolor="white", edgecolor="black", linewidth=1.0),
        )

    def _stipple(self, ax, lon, lat, sigmask: xr.DataArray, field: xr.DataArray):
        """
        Stipple only where:
          - sigmask is True
          - field is finite (i.e., where you have valid land data)
        """
        sig = np.asarray(sigmask.values).astype(bool)
        val = np.isfinite(np.asarray(field.values))

        mask = sig & val
        if mask.ndim != 2 or not np.any(mask):
            return

        jj, ii = np.where(mask)

        stride = max(1, int(self.cfg.stipple_stride))
        jj = jj[::stride]
        ii = ii[::stride]

        LON, LAT = np.meshgrid(lon.values, lat.values)

        ax.scatter(
            LON[jj, ii], LAT[jj, ii],
            s=self.cfg.stipple_size,
            c=self.cfg.stipple_color,
            alpha=self.cfg.stipple_alpha,
            transform=self.cfg.data_crs,
            linewidths=0,
            zorder=5,
        )

    @staticmethod
    def _to_0_360(ds: xr.Dataset, lon_name: str = "lon") -> xr.Dataset:
        if lon_name not in ds.coords:
            return ds
        lon = ds[lon_name]
        if lon.min() < 0 or lon.max() > 360:
            ds = ds.assign_coords({lon_name: (lon % 360)}).sortby(lon_name)
        return ds
