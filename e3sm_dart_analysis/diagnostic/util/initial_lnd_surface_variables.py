"""Ensemble bias and spread of surface variables (PRECT, TREFHT) on target dates.

``SurfaceVarDiffProcessor`` computes, for one experiment and date, the
ensemble mean, spread, bias, t-test significance and global PCC/RMSE against
an observation and writes a compact NetCDF; ``SurfaceVarDiffPlotter`` plots
bias (top) and spread (bottom) maps for several experiments.
"""

from __future__ import annotations

import os
from typing import Dict, List, Optional, Sequence

import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from util.model_daily import open_model_daily
from scipy.stats import ttest_1samp


class SurfaceVarDiffProcessor:
    """
    Compute ensemble mean, spread, bias, and significance vs obs for a single date
    and save results into a compact NetCDF file for later plotting.
    """

    def __init__(self, data_dict: Dict, *,
                 variable: str = "PRECT",
                 obs_variable: Optional[str] = None,
                 fscale: Optional[float] = None,
                 obs_fscale: Optional[float] = None,
                 landmask_file: Optional[str] = None,
                 landmask_var: str = "landfrac",
                 mask_land: bool = True,
                 confidence_level: float = 0.05):
        self.data_dict = data_dict
        self.variable = variable
        self.obs_variable = obs_variable or variable
        self.fscale = 1.0 if fscale is None else fscale
        self.obs_fscale = 1.0 if obs_fscale is None else obs_fscale
        self.landmask_var = landmask_var
        self.mask_land = mask_land
        self.confidence_level = float(confidence_level)

        self.landmask = None
        if landmask_file and os.path.exists(landmask_file):
            ds_mask = xr.open_dataset(landmask_file)
            lm = ds_mask[self.landmask_var]
            lon_name = next((k for k in lm.coords if k.lower() in ["lon","longitude"]), None)
            if lon_name:
                lon = lm[lon_name]
                if lon.min() < 0 or lon.max() > 360:
                    new_lon = (lon % 360).sortby(lon)
                    lm = lm.assign_coords({lon_name: new_lon}).sortby(lon_name)
            self.landmask = lm

    # ----------------- helpers -----------------
    def _get_lat_lon_names(self, da: xr.DataArray):
        lon_name = next((k for k in da.coords if k.lower() in ["lon", "longitude"]), None)
        lat_name = next((k for k in da.coords if k.lower() in ["lat", "latitude"]), None)
        return lon_name, lat_name

    def _apply_landmask(self, da: xr.DataArray):
        if self.mask_land and (self.landmask is not None):
            lm = self.landmask
            # interpolate nearest to avoid fill; works for 2D lat/lon grid
            lm_i = lm.interp_like(da, method="nearest")
            return da.where(lm_i >= 0.5)   # keep land
        return da

    def _wrap_longitude_for_plot(self, da: xr.DataArray):
        lon_name, _ = self._get_lat_lon_names(da)
        if lon_name is None:
            return da
        lon = da[lon_name]
        if float(lon.max() - lon.min()) < 359:  # pad to avoid seam
            pad = da.isel({lon_name: 0}).assign_coords({lon_name: lon.isel({lon_name: 0}) + 360})
            da = xr.concat([da, pad], dim=lon_name)
        return da

    def _weighted_stats(self, model: xr.DataArray, obs: xr.DataArray):
        # area weights by latitude
        lat_name = next((k for k in model.coords if k.lower() in ["lat","latitude"]), None)
        w1d = np.cos(np.deg2rad(model[lat_name]))
        w1d = w1d / w1d.mean()
        w2d = xr.ones_like(model) * w1d

        diff = (model - obs)
        rmse = np.sqrt(((diff**2) * w2d).where(np.isfinite(diff)).mean()).item()

        m = (model * w2d).where(np.isfinite(model))
        o = (obs   * w2d).where(np.isfinite(obs))
        mmean = m.mean()
        omean = o.mean()
        cov = ((m - mmean) * (o - omean)).mean()
        stdm = np.sqrt(((m - mmean)**2).mean())
        stdo = np.sqrt(((o - omean)**2).mean())
        pcc = (cov / (stdm * stdo + 1e-12)).item()
        return float(pcc), float(rmse)

    def _load_dataset(self, info: Dict, year: int, ensemble: Optional[str] = None) -> xr.Dataset:
        template = info["template"]
        path = info["path"]
        if "%(" in template:  # dict style
            vals = {"year": year}
            if ensemble is not None:
                vals["ensemble"] = ensemble
            filename = template % vals
        else:
            try:
                filename = template % year
            except Exception:
                filename = template % {"year": year}
        fp = os.path.join(path, filename)
        if not os.path.exists(fp):
            raise FileNotFoundError(f"[ERROR] File not found: {fp}")
        return xr.open_dataset(fp)

    def _load_model_ensemble_snapshot(self, model_info: Dict, target_date: pd.Timestamp) -> xr.DataArray:
        nens = int(model_info["nens"])
        ens_list: List[xr.DataArray] = []
        for n in range(1, nens + 1):
            en = f"EN{n:02d}"
            if "template" in model_info:
                ds = self._load_dataset(model_info, target_date.year, ensemble=en)
                da = ds[self.variable].sel(time=target_date, method="nearest", tolerance="12h") * self.fscale
            else:
                # Case root: daily values per util.model_daily, labelled by the day they cover.
                ts_root = os.path.join(model_info["path"], en, "archive/post/atm/180x360_aave/ts")
                daily = open_model_daily(ts_root, self.variable, en, [target_date.year])
                da = daily.sel(time=target_date.floor("D")) * self.fscale
            ens_list.append(da.squeeze())
        ens = xr.concat(ens_list, dim="ensemble").assign_coords(ensemble=np.arange(1, nens + 1))
        # NOTE: no auto-chunk here; leave chunking choice to caller or to_disk
        return ens

    def _load_obs_snapshot(self, obs_info: Dict, target_date: pd.Timestamp) -> xr.DataArray:
        ds = self._load_dataset(obs_info, target_date.year)
        obs = ds[self.obs_variable].sel(time=target_date, method="nearest", tolerance="12h") * self.obs_fscale
        # normalize longitudes
        lon_name = next((k for k in obs.coords if k.lower() in ["lon","longitude"]), None)
        if lon_name:
            lon = obs[lon_name]
            if lon.min() < 0 or lon.max() > 360:
                new_lon = (lon % 360).sortby(lon)
                obs = obs.assign_coords({lon_name: new_lon}).sortby(lon_name)
        return self._apply_landmask(obs)

    def _compute_significance_mask(self, ens: xr.DataArray, obs: xr.DataArray, alpha=0.05) -> xr.DataArray:
        if ens.shape[1:] != obs.shape:
            ens = ens.interp_like(obs)
        diff = ens - obs
        # t-test against 0 (population mean 0) along ensemble axis
        tstat, pvals = ttest_1samp(
            diff.transpose("ensemble", ...).values, popmean=0.0, axis=0, nan_policy="omit"
        )
        sig = xr.DataArray((pvals < alpha), coords=obs.coords, dims=obs.dims)
        return sig

    def process_and_save(self, *,
                         model_key: str,
                         obs_key: str = "OBS",
                         target_date: str = "2012-01-01",
                         out_nc: str):
        """
        Compute fields for one model vs one obs on one date; write NetCDF.
        """
        td = pd.to_datetime(target_date)

        model_info = self.data_dict[model_key]
        obs_info = self.data_dict[obs_key]

        ens = self._load_model_ensemble_snapshot(model_info, td)
        obs = self._load_obs_snapshot(obs_info, td)

        if ens.shape[1:] != obs.shape:
            ens = ens.interp_like(obs)

        mean_da = ens.mean("ensemble")
        std_da  = ens.std("ensemble", ddof=1)
        bias    = mean_da - obs

        # apply land mask (obs already masked in loader; but re-apply to keep consistency)
        mean_da = self._apply_landmask(mean_da)
        std_da  = self._apply_landmask(std_da)
        bias    = self._apply_landmask(bias)

        # global stats
        pcc, rmse = self._weighted_stats(mean_da, obs.interp_like(mean_da))

        # significance
        sig = self._compute_significance_mask(ens, obs, alpha=self.confidence_level)
        sig = self._apply_landmask(sig)

        # (Optional) wrap longitudes for plotting convenience; store unwrapped too if desired.
        mean_wr = self._wrap_longitude_for_plot(mean_da)
        std_wr  = self._wrap_longitude_for_plot(std_da)
        bias_wr = self._wrap_longitude_for_plot(bias)
        sig_wr  = self._wrap_longitude_for_plot(sig)
        obs_wr  = self._wrap_longitude_for_plot(obs)

        # Build dataset to save
        ds_out = xr.Dataset(
            {
                "model_mean": mean_wr,
                "model_std":  std_wr,
                "obs_field":  obs_wr,
                "bias":       bias_wr,
                "sig_mask":   sig_wr.astype("i1"),  # 0/1 to save space
            }
        )
        ds_out.attrs.update({
            "variable": self.variable,
            "obs_variable": self.obs_variable,
            "model_key": model_key,
            "obs_key": obs_key,
            "target_date": str(pd.to_datetime(target_date).date()),
            "fscale": float(self.fscale),
            "obs_fscale": float(self.obs_fscale),
            "mask_land": int(self.mask_land),   # NetCDF attributes cannot be bool
            "confidence_level": float(self.confidence_level),
            "pcc_global": pcc,
            "rmse_global": rmse,
        })

        # Lightweight encodings (set zlib if you like)
        enc = {vn: {"zlib": True, "complevel": 4} for vn in ds_out.data_vars}
        ds_out.to_netcdf(out_nc, encoding=enc)
        print(f"[SAVED] {out_nc}  (PCC={pcc:.3f}, RMSE={rmse:.3f})")


class SurfaceVarDiffPlotter:
    """
    Read preprocessed NetCDF files (from SurfaceVarDiffProcessor) and
    generate a 2xN multi-panel plot: Bias (top) and Spread (bottom).
    """

    def __init__(self, *, confidence_level: float = 0.05):
        self.confidence_level = confidence_level

    def _get_lat_lon_names(self, da: xr.DataArray):
        lon_name = next((k for k in da.coords if k.lower() in ["lon","longitude"]), None)
        lat_name = next((k for k in da.coords if k.lower() in ["lat","latitude"]), None)
        return lon_name, lat_name

    def plot_bias_and_stddev_multi_panel(self,
                                         processed_files: Sequence[str],
                                         model_labels: Optional[Sequence[str]] = None,
                                         bias_levels=None,
                                         spread_levels=None,
                                         cmap_bias="RdBu_r", cmap_std="viridis",
                                         figsize=(12,5), savepath=None, fontz=8, dpi=200):
        """
        processed_files: list of .nc paths (one per model case for the SAME variable/date)
        """
        ncols = len(processed_files)
        if model_labels is None:
            # derive from attrs if possible
            model_labels = []
            for fp in processed_files:
                attrs = xr.open_dataset(fp).attrs
                label = attrs.get("model_key", os.path.basename(fp).replace(".nc",""))
                model_labels.append(label)

        if bias_levels is None:
            bias_levels = [-5, -3, -1, -0.5, -0.1, 0.1, 0.5, 1, 3, 5]
        if spread_levels is None:
            spread_levels = [0.1, 0.5, 1, 1.5, 2, 2.5]

        bias_cmap   = plt.get_cmap(cmap_bias)
        spread_cmap = plt.get_cmap(cmap_std)
        bias_norm   = mcolors.BoundaryNorm(bias_levels,   bias_cmap.N)
        spread_norm = mcolors.BoundaryNorm(spread_levels, spread_cmap.N)

        fig, axes = plt.subplots(2, ncols, figsize=figsize,
                                 subplot_kw={"projection": ccrs.PlateCarree()})

        # plot panels
        pcc_rmse_texts = []
        for j, (fp, label) in enumerate(zip(processed_files, model_labels)):
            ds = xr.open_dataset(fp)
            bias = ds["bias"]
            std  = ds["model_std"]
            obs  = ds["obs_field"]
            sig  = ds["sig_mask"]

            # coordinates
            lon_name_b, lat_name_b = self._get_lat_lon_names(bias)
            lon_name_s, lat_name_s = self._get_lat_lon_names(std)

            # bias panel
            ax0 = axes[0, j]
            cf0 = ax0.contourf(bias[lon_name_b], bias[lat_name_b], bias,
                               levels=bias_levels, cmap=bias_cmap, norm=bias_norm,
                               transform=ccrs.PlateCarree(), extend="both")
            ax0.coastlines(linewidth=0.4)
            ax0.set_title(f"{label} - {ds.attrs.get('obs_key','OBS')} (Bias)", fontsize=fontz)
            ax0.set_xticks(np.arange(-180, 181, 60), crs=ccrs.PlateCarree())
            ax0.set_yticks(np.arange(-90, 91, 30),  crs=ccrs.PlateCarree())
            ax0.tick_params(labelsize=fontz-2)
            ax0.set_xlabel("Longitude", fontsize=fontz-1)
            ax0.set_ylabel("Latitude",  fontsize=fontz-1)

            # hatching for significant bias (p<alpha)
            sig.plot.contourf(ax=ax0, transform=ccrs.PlateCarree(),
                              colors='none', hatches=['..'], add_colorbar=False, add_labels=False)

            # annotate global PCC/RMSE
            pcc = float(ds.attrs.get("pcc_global", np.nan))
            rmse = float(ds.attrs.get("rmse_global", np.nan))
            ax0.text(0.98, 0.02, f"PCC: {pcc:.2f}\nRMSE: {rmse:.2f}",
                     transform=ax0.transAxes, fontsize=fontz-2, va='bottom', ha='right',
                     bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.1'))

            # spread panel
            ax1 = axes[1, j]
            cf1 = ax1.contourf(std[lon_name_s], std[lat_name_s], std,
                               levels=spread_levels, cmap=spread_cmap, norm=spread_norm,
                               transform=ccrs.PlateCarree(), extend="both")
            ax1.coastlines(linewidth=0.4)
            ax1.set_title(f"{label} Spread", fontsize=fontz)
            ax1.set_xticks(np.arange(-180, 181, 60), crs=ccrs.PlateCarree())
            ax1.set_yticks(np.arange(-90, 91, 30),  crs=ccrs.PlateCarree())
            ax1.tick_params(labelsize=fontz-2)
            ax1.set_xlabel("Longitude", fontsize=fontz-1)
            ax1.set_ylabel("Latitude",  fontsize=fontz-1)

        # shared colorbars on the right
        plt.subplots_adjust(wspace=0.1, hspace=0.2, right=0.88)

        cbar_ax1 = fig.add_axes([0.90, 0.60, 0.012, 0.30])
        sm1 = mpl.cm.ScalarMappable(norm=bias_norm, cmap=bias_cmap)
        sm1.set_array([])
        fig.colorbar(sm1, cax=cbar_ax1, label="Bias", extend='both')

        cbar_ax2 = fig.add_axes([0.90, 0.18, 0.012, 0.30])
        sm2 = mpl.cm.ScalarMappable(norm=spread_norm, cmap=spread_cmap)
        sm2.set_array([])
        fig.colorbar(sm2, cax=cbar_ax2, label="Spread", extend='both')

        if savepath:
            fig.savefig(savepath, dpi=dpi, bbox_inches='tight')
            print(f"[SAVED] Figure saved to: {savepath}")
        return fig
