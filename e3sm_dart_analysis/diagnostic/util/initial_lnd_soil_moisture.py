"""Soil-moisture bias/spread diagnostics against a gridded observation for one date."""

import os

import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import ttest_1samp

from util.dask_helpers import configure_dask, open_dataset_lazy


class SoilMoistureDifferencePlotter:
    def __init__(self, data_dict, landmask_file=None, soilayer_file=None,
                 variable='H2OSOI', dz_var='DZSOI', landmask_var='landfrac',
                 depth_cm=10, mask_land=True, confidence_level=0.05,
                 use_mass_units=True, use_dask=True, dask_chunks=None):
        self.data_dict = data_dict
        self.variable = variable
        self.dz_var = dz_var
        self.landmask_var = landmask_var
        self.depth_cm = depth_cm
        self.mask_land = mask_land
        self.significance_threshold = confidence_level
        self.confidence_level = confidence_level
        self.use_mass_units = use_mass_units
        self.water_density = 1000.0
        self.use_dask = use_dask
        self.dask_chunks = dask_chunks or {"time": 1, "levgrnd": -1, "lat": 90, "lon": 180}
        self.dask_client = configure_dask(use_distributed=False) if use_dask else None
        
        self.landmask = None
        if landmask_file and os.path.exists(landmask_file):
            ds_mask = open_dataset_lazy(landmask_file, chunks=self.dask_chunks) if use_dask else xr.open_dataset(landmask_file)
            self.landmask = ds_mask['landmask']
        
            # Detect longitude coordinate name (case-insensitive)
            lon_name = None
            lower_coords = {k.lower(): k for k in self.landmask.coords}
            for test in ['lon', 'longitude']:
                if test in lower_coords:
                    lon_name = lower_coords[test]
                    break
        
            # Convert to 0–360 if needed
            if lon_name:
                lon = self.landmask[lon_name]
                if lon.min() < 0 or lon.max() > 360:
                    new_lon = (lon % 360).sortby(lon)
                    self.landmask = self.landmask.assign_coords({lon_name: new_lon})
                    self.landmask = self.landmask.sortby(lon_name)

    
        self.fixed_dz = None
        if soilayer_file and os.path.exists(soilayer_file):
            ds_dz = open_dataset_lazy(soilayer_file, chunks=self.dask_chunks) if self.use_dask else xr.open_dataset(soilayer_file)
            self.fixed_dz = ds_dz[dz_var]

    def _load_dataset(self, info, year, ensemble=None):
        template = info['template']
        path = info['path']
        filename = template % {'year': year, 'ensemble': ensemble} if ensemble else template % {'year': year}
        full_path = os.path.join(path, filename)
        if not os.path.exists(full_path):
            raise FileNotFoundError(f"[ERROR] File not found: {full_path}")
        ds = open_dataset_lazy(full_path, chunks=self.dask_chunks) if self.use_dask else xr.open_dataset(full_path)

        return ds 

    def normalize_longitude_and_sort(self, ds, lon_name='lon'):
        """
        Normalize longitude in the dataset to [-180, 180] and sort by longitude.
    
        Parameters
        ----------
        ds : xr.DataArray or xr.Dataset
            Input data with a longitude coordinate.
        lon_name : str, default='lon'
            Name of the longitude coordinate.
    
        Returns
        -------
        xr.DataArray or xr.Dataset
            Dataset with normalized and sorted longitude.
        """
        if lon_name in ds.coords:
            lon = ds[lon_name]
            if lon.max() > 180:
                print("[INFO] Normalizing longitude to [-180, 180].")
                ds[lon_name] = ((lon + 180) % 360) - 180
                ds = ds.sortby(lon_name)
        return ds
    
    
        
    def _integrate_top_layers(self, ds):
        sm = ds[self.variable]

        if self.dz_var in ds:
            dz = ds[self.dz_var]
            if "time" in dz.dims:
                dz = dz.isel(time=0)
        elif self.fixed_dz is not None:
            dz = self.fixed_dz
        else:
            raise ValueError("No valid DZSOI found in dataset or external source.")

        dz_cumsum = dz.cumsum(dim='levgrnd')
        top_mask = dz_cumsum <= (self.depth_cm / 100.0)
        sm_top = sm.where(top_mask)
        dz_masked = dz.where(top_mask)
        if self.use_mass_units:
            sm_mass = sm_top * dz_masked * self.water_density
            integrated = sm_mass.sum(dim='levgrnd')
            integrated.attrs['units'] = 'kg/m2'
        else:
            integrated = (sm_top * dz_masked).sum(dim='levgrnd') / dz_masked.sum(dim='levgrnd')
            integrated.attrs['units'] = 'm3/m3'

        if self.mask_land and self.landmask is not None:
            integrated = integrated.where(self.landmask > 0.1)

        return integrated

    def _load_model_ensemble_snapshot(self, model_info, target_date):
        nens = model_info['nens']
        ens_data = []
        for n in range(1, nens + 1):
            ens = f'EN{n:02d}'
            ds = self._load_dataset(model_info, target_date.year, ensemble=ens)
            ds_day = ds.sel(time=target_date)
            sm_top = self._integrate_top_layers(ds_day)
            ens_data.append(sm_top.squeeze())
        return xr.concat(ens_data, dim='ensemble')

    def _load_obs_snapshot(self, obs_info, target_date):
        ds = self._load_dataset(obs_info, target_date.year)
        obs = ds[self.variable].sel(time=target_date)
    
        # Detect longitude coordinate name (case-insensitive)
        lon_name = None
        lower_coords = {k.lower(): k for k in obs.coords}
        for test in ['lon', 'longitude']:
            if test in lower_coords:
                lon_name = lower_coords[test]
                break
                
        # Convert longitude from [-180, 180] to [0, 360] if needed
        if lon_name:
            lon = obs[lon_name]
            if lon.min() < 0 or lon.max() > 360:
                new_lon = (lon % 360).sortby(lon)
                obs = obs.assign_coords({lon_name: new_lon})
                obs = obs.sortby(lon_name)
    
        if self.use_mass_units:
            obs = obs * 1000.0 * self.depth_cm / 100.0
    
        if self.mask_land and self.landmask is not None:
            obs = obs.where(self.landmask > 0.1)
    
        return obs

    def _wrap_longitude_for_plot(self, da):
        if 'lon' not in da.coords:
            return da
    
        lon = da.lon
        if (lon[-1] - lon[0]) < 359:  # Avoid double-padding
            da = da.pad(lon=(0, 1), constant_values=np.nan)
            da['lon'][-1] = da['lon'][0] + 360  # Wrap around
        return da
        
    def _compute_significance_mask(self, ensemble_data, obs_mean, alpha=0.05, threshold=None):
        if threshold is None:
            threshold = self.significance_threshold

        diff = ensemble_data - obs_mean
        diff_np = diff.transpose("ensemble", "lat", "lon").compute().values
        tstat, pvals = ttest_1samp(diff_np,
                                   popmean=0.0, axis=0, nan_policy='omit')

        sig_mask = xr.DataArray(pvals < alpha, coords=obs_mean.coords, dims=obs_mean.dims)
        model_mean = ensemble_data.mean(dim='ensemble')
        valid_mask = (abs(model_mean) > threshold) & (abs(obs_mean) > threshold)
        return sig_mask.where(valid_mask)

    def process_bias_and_spread(self, model_keys, target_date="2012-01-01", obs_key="ESA_CCI"):
        # Compute reusable gridded diagnostics before any plotting.
        target_date = pd.to_datetime(target_date)
        obs_field = self._load_obs_snapshot(self.data_dict[obs_key], target_date)
        processed = []
        for model_key in model_keys:
            ensemble = self._load_model_ensemble_snapshot(self.data_dict[model_key], target_date)
            if ensemble.shape[1:] != obs_field.shape:
                ensemble = ensemble.interp_like(obs_field)
            model_mean = ensemble.mean(dim="ensemble")
            spread = ensemble.std(dim="ensemble")
            bias = model_mean - obs_field
            significance = self._compute_significance_mask(
                ensemble, obs_field, alpha=self.confidence_level
            )
            model_values = np.asarray(model_mean.values).ravel()
            obs_values = np.asarray(obs_field.values).ravel()
            valid = np.isfinite(model_values) & np.isfinite(obs_values)
            if valid.any():
                pcc = float(np.corrcoef(model_values[valid], obs_values[valid])[0, 1])
                rmse = float(np.sqrt(np.mean((model_values[valid] - obs_values[valid]) ** 2)))
            else:
                pcc = np.nan
                rmse = np.nan
            processed.append(xr.Dataset({
                "bias": bias,
                "spread": spread,
                "significant": significance,
                "pcc": xr.DataArray(pcc),
                "rmse": xr.DataArray(rmse),
            }).expand_dims(experiment=[model_key]))
        output = xr.concat(processed, dim="experiment")
        output.attrs.update({
            "observation": obs_key,
            "target_date": target_date.strftime("%Y-%m-%d"),
            "depth_cm": self.depth_cm,
            "use_mass_units": int(self.use_mass_units),
            "mask_land": int(self.mask_land),
            "confidence_level": self.confidence_level,
        })
        return output


    def plot_bias_and_stddev_multi_panel(self, model_keys, target_date='2012-01-01',
                                         obs_key='ESA_CCI', bias_levels=None, spread_levels=None,
                                         cmap_bias='RdBu_r', cmap_std='viridis',
                                         figsize=(14.0, 6.6), savepath=None, fontz=9, show=True,
                                         processed_data=None):
        target_date = pd.to_datetime(target_date)
        ncols = len(model_keys)
    
        if bias_levels is None:
            bias_levels = [-0.5, -0.3, -0.1, -0.05, -0.01, 0.01, 0.05, 0.1, 0.3, 0.5]
        if spread_levels is None:
            spread_levels = [0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1]
    
        bias_cmap = mpl.colors.ListedColormap(plt.get_cmap(cmap_bias)(np.linspace(0, 1, len(bias_levels) - 1)))
        bias_norm = mcolors.BoundaryNorm(bias_levels, ncolors=bias_cmap.N)
        spread_cmap = mpl.colors.ListedColormap(plt.get_cmap(cmap_std)(np.linspace(0, 1, len(spread_levels) - 1)))
        spread_norm = mcolors.BoundaryNorm(spread_levels, ncolors=spread_cmap.N)
        bias_label = "Bias (kg m$^{-2}$)" if self.use_mass_units else "Bias (m$^3$ m$^{-3}$)"
        spread_label = "Spread (kg m$^{-2}$)" if self.use_mass_units else "Spread (m$^3$ m$^{-3}$)"
    
        fig, axes = plt.subplots(
            nrows=2, ncols=ncols, figsize=figsize,
            subplot_kw={'projection': ccrs.Robinson()},
            constrained_layout=False
        )
    
        for j, model_key in enumerate(model_keys):
            if processed_data is None:
                current = self.process_bias_and_spread(
                    [model_key], target_date=target_date, obs_key=obs_key
                ).sel(experiment=model_key)
            else:
                current = processed_data.sel(experiment=model_key)
            bias = current["bias"]
            model_std = current["spread"]
            sig_mask = current["significant"]
            pcc = float(current["pcc"].values)
            rmse_val = float(current["rmse"].values)
            bias = self._wrap_longitude_for_plot(bias)
            model_std = self._wrap_longitude_for_plot(model_std)
    
            # Bias panel
            ax0 = axes[0, j]
            bias_mesh = ax0.contourf(
                bias['lon'], bias['lat'], bias.values,
                levels=bias_levels,
                cmap=bias_cmap,
                norm=bias_norm,
                transform=ccrs.PlateCarree(),
                extend='both'
            )
    
            ax0.set_title(f"{model_key} - {obs_key} (Bias)", fontsize=fontz, pad=8)
            ax0.coastlines(linewidth=0.4)
            ax0.tick_params(labelsize=fontz - 1)
            
            ax0.text(
                0.98, 0.02, f"PCC: {pcc:.2f}\nRMSE: {rmse_val:.2f}",
                transform=ax0.transAxes, fontsize=fontz - 2, va='bottom', ha='right',
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.1')
            )
    
            sig_mask.plot.contourf(
                ax=ax0, transform=ccrs.PlateCarree(),
                colors='none', hatches=['..', '..'], 
                add_colorbar=False, 
                add_labels=False,
                linewidths=0
            )
    
            # Spread panel
            ax1 = axes[1, j]
            spread_mesh = model_std.plot.pcolormesh(
                ax=ax1, transform=ccrs.PlateCarree(),
                cmap=spread_cmap,
                norm=spread_norm,
                add_labels=False,
                add_colorbar=False
            )
    
            ax1.set_title(f"{model_key} (Spread)", fontsize=fontz, pad=8)
            ax1.coastlines(linewidth=0.4)
            ax1.tick_params(labelsize=fontz - 1)

            bias_cbar = fig.colorbar(
                bias_mesh, ax=ax0, orientation='horizontal',
                ticks=bias_levels[::2], fraction=0.055, pad=0.08, aspect=24
            )
            bias_cbar.set_label(bias_label, fontsize=fontz - 1, labelpad=2)
            bias_cbar.ax.tick_params(labelsize=fontz - 2, pad=1)

            spread_cbar = fig.colorbar(
                spread_mesh, ax=ax1, orientation='horizontal',
                ticks=spread_levels[::2], fraction=0.055, pad=0.08, aspect=24
            )
            spread_cbar.set_label(spread_label, fontsize=fontz - 1, labelpad=2)
            spread_cbar.ax.tick_params(labelsize=fontz - 2, pad=1)
    
        # Adjust layout after per-panel colorbars are attached.
        plt.subplots_adjust(left=0.035, right=0.985, bottom=0.09, top=0.92, wspace=0.12, hspace=0.46)
    
        if savepath:
            fig.savefig(savepath, dpi=600, bbox_inches='tight')
            print(f"[SAVED] Figure saved to: {savepath}")

        if show:
            plt.show()
        else:
            plt.close(fig)

        return fig, axes
