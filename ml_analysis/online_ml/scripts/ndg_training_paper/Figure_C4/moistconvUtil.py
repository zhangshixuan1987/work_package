import os
import json
import time
import xarray as xr
import numpy as np
import colorcet as cc

import scipy.ndimage
from scipy.ndimage import gaussian_filter1d
from scipy.ndimage import gaussian_filter
from scipy.stats import ttest_ind

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import TwoSlopeNorm, ListedColormap, LinearSegmentedColormap
from matplotlib import ticker
import matplotlib.gridspec as gridspec

import cartopy.crs as ccrs
import cartopy.feature as cfeature

def create_masked_colormap(cmap, clevels):
    base_cmap = cmap(np.linspace(0, 1, len(clevels) - 1))
    alphas = np.ones(len(clevels) - 1)
    mid = len(alphas) // 2
    alphas[mid - 1: mid + 1] = 0.0  # transparent center
    rgba = np.concatenate([base_cmap[:, :3], alphas[:, None]], axis=1)
    return LinearSegmentedColormap.from_list("masked", rgba, N=len(rgba))

# === Color utilities ===
def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def rgb_to_dec(rgb_tuple):
    return [val / 255 for val in rgb_tuple]

def get_continuous_cmap(hex_list, float_list=None):
    if float_list is not None:
        if len(float_list) != len(hex_list):
            raise ValueError("float_list must be the same length as hex_list")
        if float_list[0] != 0.0 or float_list[-1] != 1.0:
            raise ValueError("float_list must start with 0.0 and end with 1.0")
    else:
        float_list = np.linspace(0, 1, len(hex_list))

    rgb_list = [rgb_to_dec(hex_to_rgb(hx)) for hx in hex_list]

    cdict = {"red": [], "green": [], "blue": []}
    for i, (r, g, b) in enumerate(rgb_list):
        f = float_list[i]
        cdict["red"].append([f, r, r])
        cdict["green"].append([f, g, g])
        cdict["blue"].append([f, b, b])

    return mcolors.LinearSegmentedColormap('custom_cmap', segmentdata=cdict, N=256)

class LSMCPlotter:
    def __init__(self, model_names, ref_model, data_path, varname="LSMC",
                 clim_model="CLIM",
                 dif_cmap="RdBu_r", clevels=None, vmin=None, vmax=None,
                 raw_cmap="YlGnBu", raw_clevels=None,
                 sig_level=0.05, smooth_sigma=1.0, export_nc=False):
        self.varname = varname
        self.dif_cmap = dif_cmap
        self.clevels = clevels
        self.vmin = vmin if vmin is not None else clevels[0]
        self.vmax = vmax if vmax is not None else clevels[-1]
        self.sig_level = sig_level
        self.smooth_sigma = smooth_sigma
        self.model_names = model_names
        self.ref_model = ref_model
        self.clim_model = clim_model
        self.data_path = data_path
        self.raw_cmap = raw_cmap
        self.raw_clevels = raw_clevels if raw_clevels is not None else np.linspace(0, 0.005, 11)
        self.raw_vmin = self.raw_clevels[0]
        self.raw_vmax = self.raw_clevels[-1]
        self.diff_clevels = clevels
        self.diff_vmin = self.vmin
        self.diff_vmax = self.vmax
        self.export_nc = export_nc

    def _normalize_longitude(self, ds):
        lon = ds.lon.values
        if np.any(lon > 180):
            lon = np.where(lon > 180, lon - 360, lon)
            ds = ds.assign_coords(lon=lon)
            ds = ds.sortby("lon")
        return ds

    def _open_stack(self, file_list):
        datasets = []
        for i, f in enumerate(file_list):
            da = xr.open_dataset(f)[self.varname]
            da = self._normalize_longitude(da)
            da = da.expand_dims(time=[i])
            datasets.append(da)
        return xr.concat(datasets, dim="time")

    def _get_model_file(self, season, years):
        model_dict = {}
        for model in self.model_names:
            model_files = [
                os.path.join(self.data_path, f"{model}_LSMC_{season}_{year}.nc")
                for year in years
                if os.path.exists(os.path.join(self.data_path, f"{model}_LSMC_{season}_{year}.nc"))
            ]
            if len(model_files) >= 5:
                model_dict[model] = model_files
            else:
                print(f"[SKIP] Not enough data for model {model} in {season}")

        ref_files = [
            os.path.join(self.data_path, f"{self.ref_model}_LSMC_{season}_{year}.nc")
            for year in years
            if os.path.exists(os.path.join(self.data_path, f"{self.ref_model}_LSMC_{season}_{year}.nc"))
        ]

        if len(ref_files) < 5 or not model_dict:
            print(f"[SKIP] Skipping {season}: missing reference or no valid models")
            return None, None

        return model_dict, ref_files

    def _maybe_smooth(self, da):
        if self.smooth_sigma <= 0:
            return da
        valid_mask = ~np.isnan(da)
        filled = da.fillna(0)
        smoothed = gaussian_filter(filled, sigma=(self.smooth_sigma, self.smooth_sigma))
        norm = gaussian_filter(valid_mask.astype(float), sigma=(self.smooth_sigma, self.smooth_sigma))
        result = smoothed / norm
        return xr.DataArray(result, coords=da.coords, dims=da.dims).where(valid_mask)

    def _plot_reference_panel(self, ax, da, title, units):
        da_clipped = da.clip(min=self.raw_vmin, max=self.raw_vmax)
        mappable = da_clipped.plot.contourf(
            ax=ax,
            levels=self.raw_clevels,
            cmap=self.raw_cmap,
            extend="neither",
            transform=ccrs.PlateCarree(),
            add_colorbar=False
        )
        ax.coastlines(linewidth=0.5, color="black")
        ax.set_title(title, fontsize=12)
        gl = ax.gridlines(draw_labels=True, linewidth=0.2, color="gray", alpha=0.5)
        gl.top_labels = gl.right_labels = False
        return mappable

    def _plot_model_diff_panel(self, ax, model, model_data, ref_data, units):
        if ref_data.time.size < 5 or model_data.time.size < 5:
            print(f"[WARNING] Limited samples for t-test in {model}")

        model_mean = model_data.mean(dim="time")
        diff = model_mean - ref_data.mean(dim="time")
        smoothed_diff = self._maybe_smooth(diff).clip(min=self.diff_vmin, max=self.diff_vmax)

        tstat, pval = ttest_ind(model_data, ref_data, axis=0, equal_var=False, nan_policy="omit")
        sig_mask_np = pval < self.sig_level
        sig_mask = xr.DataArray(sig_mask_np, coords=diff.coords, dims=diff.dims)
        sig_mask = sig_mask.where(~np.isnan(diff))

        masked_cmap = create_masked_colormap(self.dif_cmap, self.diff_clevels)
        mappable = smoothed_diff.plot.contourf(
            ax=ax,
            levels=self.diff_clevels,
            cmap=masked_cmap,
            extend="neither",
            transform=ccrs.PlateCarree(),
            add_colorbar=False
        )

        try:
            sig_plot_mask = sig_mask.where(sig_mask, 0).fillna(0).astype(float)
            ax.contourf(
                sig_plot_mask.lon, sig_plot_mask.lat, sig_plot_mask.values,
                levels=[0.5, 1.5], hatches=[".."], colors="none",
                transform=ccrs.PlateCarree()
            )
        except Exception as e:
            print(f"[WARNING] Hatching failed for {model}: {e}")

        ax.coastlines(linewidth=0.5, color="black")
        ax.set_title(f"{model} - {self.ref_model}", fontsize=12)
        gl = ax.gridlines(draw_labels=True, linewidth=0.2, color="gray", alpha=0.5)
        gl.top_labels = gl.right_labels = False
        return mappable

    def plot_multi_model_comparison(
            self,
            season="ANN", 
            years=None, 
            save_path="./",
            fig_name=None
        ):
        model_dict, ref_files = self._get_model_file(season, years)
        if model_dict is None or ref_files is None:
            return

        if self.clim_model not in model_dict:
            print(f"[ERROR] '{self.clim_model}' model must be included in model_names for layout.")
            return

        ref_data = self._open_stack(ref_files)
        ref_mean = ref_data.mean(dim="time")
        clim_data = self._open_stack(model_dict[self.clim_model])
        clim_mean = clim_data.mean(dim="time")
        units = ref_mean.attrs.get("units", "unknown")

        if self.export_nc:
            ref_mean.to_netcdf(f"{save_path}/{self.ref_model}_mean_{season}.nc")
            clim_mean.to_netcdf(f"{save_path}/{self.clim_model}_mean_{season}.nc")

        fig = plt.figure(figsize=(10, 10), constrained_layout=True)
        gs = fig.add_gridspec(5, 2, height_ratios=[1, 0.02, 0.7, 0.7, 0.01])

        axs_row1 = [fig.add_subplot(gs[0, i], projection=ccrs.PlateCarree()) for i in range(2)]
        m1 = self._plot_reference_panel(axs_row1[0], self._maybe_smooth(ref_mean), f"{self.ref_model} (Reference)", units)
        m2 = self._plot_reference_panel(axs_row1[1], self._maybe_smooth(clim_mean), f"{self.clim_model}", units)

        raw_cbar_ax = fig.add_axes([0.15, 0.62, 0.75, 0.015])
        plt.colorbar(
            m1,
            cax=raw_cbar_ax,
            orientation="horizontal",
            boundaries=self.raw_clevels,
            ticks=self.raw_clevels,
            extend="neither"
        ).set_label(f"{self.varname} ({units})")

        diff_models = [m for m in model_dict if m != self.ref_model]
        all_axes = [fig.add_subplot(gs[r, c], projection=ccrs.PlateCarree())
                    for r in [2, 3] for c in [0, 1]]

        if len(diff_models) > len(all_axes):
            print(f"[WARNING] Too many models; only showing first {len(all_axes)}")

        mappable_list = []
        for i, model in enumerate(diff_models[:len(all_axes)]):
            ax = all_axes[i]
            model_data = self._open_stack(model_dict[model])
            if self.export_nc:
                model_data.mean(dim="time").to_netcdf(f"{save_path}/{model}_mean_{season}.nc")
            mappable = self._plot_model_diff_panel(ax, model, model_data, ref_data, units)
            mappable_list.append(mappable)

        diff_cbar_ax = fig.add_axes([0.15, -0.02, 0.75, 0.015])
        if mappable_list:
            plt.colorbar(
                mappable_list[0],
                cax=diff_cbar_ax,
                orientation="horizontal",
                boundaries=self.diff_clevels,
                ticks=self.diff_clevels,
                extend="neither"
            ).set_label(f"Δ{self.varname} ({units})")

        if fig_name is not None:
            fig_name = os.path.join(save_path, fig_name)
        else:
            fig_name = os.path.join(save_path, f"raw_and_diff_LSMC_{season}.pdf")
            
        plt.savefig(fig_name, dpi=600, bbox_inches="tight", pad_inches=0.05)
        plt.show()
        plt.close()
        print(f"[PLOT] 3-row comparison saved: {fig_name}")

class LSMCPlotterROW2:
    def __init__(self, model_names, ref_model, data_path, varname="LSMC",
                 dif_cmap="RdBu_r", clevels=None, vmin=None, vmax=None,
                 raw_cmap="YlGnBu", raw_clevels=None,
                 sig_level=0.05, smooth_sigma=1.0):
        self.varname = varname
        self.dif_cmap = dif_cmap
        self.clevels = clevels
        self.vmin = vmin if vmin is not None else clevels[0]
        self.vmax = vmax if vmax is not None else clevels[-1]
        self.sig_level = sig_level
        self.smooth_sigma = smooth_sigma
        self.model_names = model_names
        self.ref_model = ref_model
        self.data_path = data_path
        self.raw_cmap = raw_cmap
        self.raw_clevels = raw_clevels if raw_clevels is not None else np.linspace(0, 0.005, 11)
        self.raw_vmin = self.raw_clevels[0]
        self.raw_vmax = self.raw_clevels[-1]
        self.diff_clevels = clevels
        self.diff_vmin = self.vmin
        self.diff_vmax = self.vmax

    def _normalize_longitude(self, ds):
        lon = ds.lon.values
        if np.any(lon > 180):
            lon = np.where(lon > 180, lon - 360, lon)
            ds = ds.assign_coords(lon=lon)
            ds = ds.sortby("lon")
        return ds

    def _open_stack(self, file_list):
        datasets = []
        for i, f in enumerate(file_list):
            da = xr.open_dataset(f)[self.varname]
            da = self._normalize_longitude(da)
            da = da.expand_dims(time=[i])
            datasets.append(da)
        return xr.concat(datasets, dim="time")

    def _get_model_file(self, season, years):
        model_dict = {}
        for model in self.model_names:
            model_files = [
                os.path.join(self.data_path, f"{model}_LSMC_{season}_{year}.nc")
                for year in years
                if os.path.exists(os.path.join(self.data_path, f"{model}_LSMC_{season}_{year}.nc"))
            ]
            if len(model_files) >= 5:
                model_dict[model] = model_files
            else:
                print(f"[SKIP] Not enough data for model {model} in {season}")

        ref_files = [
            os.path.join(self.data_path, f"{self.ref_model}_LSMC_{season}_{year}.nc")
            for year in years
            if os.path.exists(os.path.join(self.data_path, f"{self.ref_model}_LSMC_{season}_{year}.nc"))
        ]

        if len(ref_files) < 5 or not model_dict:
            print(f"[SKIP] Skipping {season}: missing reference or no valid models")
            return None, None

        return model_dict, ref_files

    def _latlon_gaussian_filter(self, da, sigma=(1, 1)):
        valid_mask = ~np.isnan(da)
        filled = da.fillna(0)
        smoothed = gaussian_filter(filled, sigma=sigma)
        norm = gaussian_filter(valid_mask.astype(float), sigma=sigma)
        result = smoothed / norm
        result = xr.DataArray(result, coords=da.coords, dims=da.dims)
        return result.where(valid_mask)

    def plot_multi_model_comparison(
            self, 
            season="ANN", 
            years=None, 
            save_path="./",
            fig_name=None,
        ):
        model_dict, ref_files = self._get_model_file(season, years)
        if model_dict is None or ref_files is None:
            return

        if "CLIM" not in model_dict:
            print("[ERROR] 'CLIM' model must be included in model_names for layout.")
            return

        ref_data = self._open_stack(ref_files)
        ref_mean = ref_data.mean(dim="time")

        clim_data = self._open_stack(model_dict["CLIM"])
        clim_mean = clim_data.mean(dim="time")

        diff_models = [m for m in model_dict if m != self.ref_model]
        ncols = len(diff_models)

        fig = plt.figure(figsize=(5 * ncols, 9))
        gs = gridspec.GridSpec(4, ncols, height_ratios=[1, 0.03, 1, 0.03], hspace=0.01, wspace=0.15)

        top_width = 0.24
        left_start = 0.5 - (top_width * 2 + 0.05) / 2
        axs_top = [
            fig.add_axes([left_start + i * (top_width + 0.04), 0.57, top_width, 0.25], projection=ccrs.PlateCarree())
            for i in range(2)
        ]

        axs_bot = [fig.add_subplot(gs[2, j], projection=ccrs.PlateCarree()) for j in range(ncols)]

        row1_data = [
            self._latlon_gaussian_filter(ref_mean, sigma=(self.smooth_sigma, self.smooth_sigma)),
            self._latlon_gaussian_filter(clim_mean, sigma=(self.smooth_sigma, self.smooth_sigma))
        ]
        row1_titles = [f"{self.ref_model} (Reference)", "CLIM (EAMv2)"]
        raw_plot_list = []
        for i, ax in enumerate(axs_top):
            da = row1_data[i].clip(min=self.raw_vmin, max=self.raw_vmax)
            mappable = da.plot.contourf(
                ax=ax,
                levels=self.raw_clevels,
                cmap=self.raw_cmap,
                extend="neither",
                transform=ccrs.PlateCarree(),
                add_colorbar=False
            )
            raw_plot_list.append(mappable)
            ax.coastlines(linewidth=0.5, color="black")
            gl = ax.gridlines(draw_labels=True, linewidth=0.2, color="gray", alpha=0.5)
            gl.top_labels = gl.right_labels = False
            ax.set_title(row1_titles[i], fontsize=fontz)

        raw_cbar_ax = fig.add_axes([0.5 - 0.2, 0.535, 0.4, 0.015])
        raw_cbar = plt.colorbar(
            raw_plot_list[0],
            cax=raw_cbar_ax,
            orientation="horizontal",
            boundaries=self.raw_clevels,
            ticks=self.raw_clevels,
            spacing="uniform",
            extend="neither"
        )
        raw_cbar.outline.set_visible(True)
        raw_cbar.outline.set_edgecolor("black")
        raw_cbar.outline.set_linewidth(0.3)
        raw_cbar.set_label(f"{self.varname} ({ref_mean.attrs.get('units', 'mm/day')})")
        raw_cbar.ax.tick_params(direction='out', color="black", length=5, width=0.1)

        diff_plot_list = []
        for j, model in enumerate(diff_models):
            ax = axs_bot[j]
            model_data = self._open_stack(model_dict[model])
            model_mean = model_data.mean(dim="time")
            diff = model_mean - ref_mean
            smoothed_diff = self._latlon_gaussian_filter(diff, sigma=(self.smooth_sigma, self.smooth_sigma))
            clipped_diff = smoothed_diff.clip(min=self.diff_vmin, max=self.diff_vmax)

            tstat, pval = ttest_ind(model_data, ref_data, axis=0, equal_var=False, nan_policy="omit")
            sig_mask_np = pval < self.sig_level
            sig_mask = xr.DataArray(sig_mask_np, coords=diff.coords, dims=diff.dims)
            sig_mask = sig_mask.where(~np.isnan(diff))

            masked_cmap = create_masked_colormap(self.dif_cmap, self.diff_clevels)
            mappable = clipped_diff.plot.contourf(
                ax=ax,
                levels=self.diff_clevels,
                cmap=masked_cmap,
                extend="neither",
                transform=ccrs.PlateCarree(),
                add_colorbar=False
            )
            diff_plot_list.append(mappable)

            # === FIXED SIGNIFICANCE HATCHING ===
            try:
                sig_plot_mask = sig_mask.where(sig_mask, 0).fillna(0).astype(float)
                lon = sig_plot_mask.lon.values
                lat = sig_plot_mask.lat.values
                hatch_contour = ax.contourf(
                    lon, lat,
                    sig_plot_mask.values,
                    levels=[0.5, 1.5],
                    hatches=["..."],
                    colors="none",
                    transform=ccrs.PlateCarree()
                )
                for c in hatch_contour.collections:
                    c.set_edgecolor((0.3, 0.3, 0.3, 0.6))  # RGBA: grey with transparency
                    c.set_facecolor("none")  # Ensure no shading is shown
                    c.set_linewidth(0)

            except Exception as e:
                print(f"[WARNING] Failed to plot significance contours: {e}")

            ax.coastlines(linewidth=0.5, color="black")
            gl = ax.gridlines(draw_labels=True, linewidth=0.2, color="gray", alpha=0.5)
            gl.top_labels = gl.right_labels = False
            ax.set_title(f"{model} - {self.ref_model}", fontsize=fontz)

        diff_cbar_ax = fig.add_axes([0.5 - 0.25, 0.14, 0.5, 0.015])
        diff_cbar = plt.colorbar(
            diff_plot_list[0],
            cax=diff_cbar_ax,
            orientation="horizontal",
            boundaries=self.diff_clevels,
            ticks=self.diff_clevels,
            spacing="uniform",
            extend="neither"
        )
        diff_cbar.outline.set_visible(True)
        diff_cbar.outline.set_edgecolor("black")
        diff_cbar.outline.set_linewidth(0.3)
        diff_cbar.set_label(f"Δ{self.varname} ({ref_mean.attrs.get('units', 'mm/day')})")
        diff_cbar.ax.tick_params(direction='out', color="black", length=5, width=0.1)

        if fig_name is not None:
            fig_name = os.path.join(save_path, fig_name)
        else:
            fig_name = os.path.join(save_path, f"raw_and_diff_LSMC_{season}.pdf")
        plt.savefig(fig_name, dpi=600, bbox_inches="tight", pad_inches=0.05)
        plt.show()
        plt.close()
        print(f"[PLOT] 2-row comparison saved: {fig_name}")

