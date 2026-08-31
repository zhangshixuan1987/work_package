import os
import time
import xarray as xr
import numpy as np
import scipy.ndimage
from scipy.ndimage import gaussian_filter1d
from scipy.ndimage import gaussian_filter
import colorcet as cc

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import TwoSlopeNorm, ListedColormap, LinearSegmentedColormap
from matplotlib import ticker
from matplotlib.colors import BoundaryNorm

import cartopy.crs as ccrs
import cartopy.feature as cfeature

def get_reference_data(freq="monthly", data_root=None):
    # Validate frequency
    valid_freqs = {"monthly", "daily"}
    if freq not in valid_freqs:
        raise ValueError(f"[ERROR] Unsupported frequency: '{freq}'. Supported options: {valid_freqs}")

    # Define base data directory
    base_root = data_root or "/anvil/scratch/x-szhang3/post_data"

    ref_dict = {
        "TREFHT": {
            "source": "ERA5", "alias": "TREFHT", "fscl": 1.0,
            "run": os.path.join(base_root, "ERA5", freq),
            "template": "ERA5_analysis_monthly_%(year).nc",
            "period": "200801_201712"
        },
        "TMQ": {
            "source": "ERA5", "alias": "TMQ", "fscl": 1.0,
            "run": os.path.join(base_root, "ERA5", freq),
            "template": "ERA5_analysis_monthly_%(year).nc",
            "period": "200801_201712"
        },
        "TS": {
            "source": "ERA5", "alias": "TS", "fscl": 1.0,
            "run": os.path.join(base_root, "ERA5", freq),
            "template": "ERA5_analysis_monthly_%(year).nc",
            "period": "200801_201712"
        },
        "FLUT": {
            "source": "CERES", "alias": "OLR", "fscl": 1.0,
            "run": os.path.join(base_root, "OBS", "CERES-OAFlux", freq),
            "template": "CERES-OAFlux_%(year).nc",
            "period": "200101-201812"
        },
        "SWCF": {
            "source": "CERES", "alias": "SWCRE", "fscl": 1.0,
            "run": os.path.join(base_root, "OBS", "CERES-OAFlux", freq),
            "template": "CERES-OAFlux_%(year).nc",
            "period": "200101-201812"
        },
        "LWCF": {
            "source": "CERES", "alias": "LWCRE", "fscl": 1.0,
            "run": os.path.join(base_root, "OBS", "CERES-OAFlux", freq),
            "template": "CERES-OAFlux_%(year).nc",
            "period": "200101-201812"
        },
        "NETCRE": {
            "source": "CERES", "alias": "CRE", "fscl": 1.0,
            "components": ["SWCF", "LWCF"], "operation": "sum", 
            "run": os.path.join(base_root, "OBS", "CERES-OAFlux", freq),
            "template": "CERES-OAFlux_%(year).nc",
            "period": "200101-201812"
        },
        "CLDTOT": {
            "source": "CERES", "alias": "CLDTOT", "fscl": 1.0,
            "run": os.path.join(base_root, "OBS", "CERES-OAFlux", freq),
            "template": "CERES-OAFlux_%(year).nc",
            "period": "200101-201812"
        },
        "PRECT": {
            "source": "GPCP", "alias": "PRECT", "fscl": 1.0,
            "run": os.path.join(base_root, "OBS", "GPCP", freq),
            "template": "GPCP_%(year).nc",
            "period": "197901-201912"
        }
    }

    return ref_dict
    
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

def create_centered_colormap_with_white(base_cmap, clevels, white_width=1):
    """
    Create a custom diverging colormap that passes through white at the center.

    Parameters:
        base_cmap : a matplotlib or colorcet continuous colormap
        clevels : array-like of contour levels (length = bins + 1)
        white_width : number of color bins to make white around 0

    Returns:
        ListedColormap for use in contourf
    """
    n_bins = len(clevels) + 1
    colors = base_cmap(np.linspace(0, 1, n_bins))

    # Find the index of the bin closest to zero
    bin_centers = clevels[:-1] + np.diff(clevels) / 2
    zero_idx = np.argmin(np.abs(bin_centers))

    half_width = int(white_width // 2)
    white_rgb = np.array([1, 1, 1, 1])  # RGBA white

    for i in range(zero_idx - half_width+2, zero_idx + half_width+2):
        if 0 <= i < n_bins:
            colors[i] = white_rgb

    return ListedColormap(colors)

def plot_surface_metrics_comparison(base_dir, fig_dir, model_dict, var_dict,
                                    sig_level=0.1, cmap=None, fontz=None,
                                    season="ANN",
                                    fig_name=None,
                                    freq="monthly",
                                    output_formats=["pdf", "png"],
                                    panel_label_mode="rowcol"):
    os.makedirs(fig_dir, exist_ok=True)
    ref_dict = get_reference_data(freq=freq)
    if cmap is None:
        import colorcet as cc
        cmap = cc.cm["CET_C9s"]

    nrow, ncol = len(var_dict), len(model_dict)
    if fontz is not None:
        fontz = fontz
        tick_fontz = fontz*0.85
    else:
        fontz = 14
        tick_fontz = 12
    fig = plt.figure(figsize=(4.5 * ncol + 2, 3.0 * nrow))
    scale = nrow / 4.0
    start_time = time.time()
    missing_files = []

    for i, (var_key, vinfo) in enumerate(var_dict.items()):
        clevels = np.unique(np.sort(np.array(vinfo["level"])))
        cmap_custom = create_centered_colormap_with_white(cmap, clevels, white_width=2)
        norm = BoundaryNorm(boundaries=clevels, ncolors=cmap_custom.N, extend='both')
        cf_last = None

        for j, model_key in enumerate(model_dict):
            ax = plt.subplot(nrow, ncol, i * ncol + j + 1, projection=ccrs.PlateCarree())
            fname = os.path.join(base_dir, f"{var_key}_global_clim_ttest_{model_key}.nc")

            if not os.path.exists(fname):
                ax.set_title(f"{model_key}\n[Missing]", fontsize=fontz)
                ax.axis("off")
                missing_files.append(fname)
                continue

            try:
                ds = xr.open_dataset(fname)
            except Exception:
                ax.set_title(f"{model_key}\n[Read Error]", fontsize=fontz)
                ax.axis("off")
                continue

            if season not in ds.season.values:
                ax.set_title(f"{model_key}\n[No {season}]", fontsize=fontz)
                ax.axis("off")
                continue

            diff = ds["diff"].sel(season=season)
            pval = ds["p_value"].sel(season=season)

            smoothed_diff = xr.apply_ufunc(
                gaussian_filter, diff,
                kwargs={"sigma": 1.0},
                input_core_dims=[["lat", "lon"]],
                output_core_dims=[["lat", "lon"]],
                vectorize=True, dask="parallelized",
                output_dtypes=[diff.dtype],
            )

            cf = ax.contourf(diff.lon, diff.lat, smoothed_diff, levels=clevels,
                             cmap=cmap_custom, norm=norm, extend="both", transform=ccrs.PlateCarree())
            cf_last = cf

            masked = xr.where((smoothed_diff < -0.1) | (smoothed_diff > 0.1), smoothed_diff, np.nan)
            contour_levels = [lvl for lvl in clevels if abs(lvl) > 0.1]
            if contour_levels:
                ax.contour(masked.lon, masked.lat, masked,
                           levels=contour_levels, colors='black',
                           linewidths=0.2, alpha=0.3, transform=ccrs.PlateCarree())

            sig_mask = xr.where(pval < sig_level, 1, np.nan)
            ax.contourf(sig_mask.lon, sig_mask.lat, sig_mask,
                        levels=[0, 2], colors='none', hatches=['..'],
                        transform=ccrs.PlateCarree())

            ax.coastlines(linewidth=0.5)
            gl = ax.gridlines(draw_labels=True, linewidth=0.2, color='gray', linestyle='--')
            gl.top_labels = False
            gl.right_labels = False
            gl.xlabel_style = {'size': tick_fontz}
            gl.ylabel_style = {'size': tick_fontz}

            # Panel label logic
            if panel_label_mode == "global":
                panel_label = chr(97 + i * ncol + j)  # a, b, c, ...
            else:
                panel_label = f"{chr(97 + i)}{j + 1}"  # a1, a2, ...

            exp_display = model_dict[model_key]
            ref_name = ref_dict.get(var_key, {}).get("source", "OBS")
            ax.set_title(f"({panel_label}) {exp_display} - {ref_name}", loc="left", fontsize=fontz)
            ax.set_title(vinfo['alias'], loc="right", fontsize=fontz)

        # Colorbar per row
        if cf_last:
            cbar_ax = fig.add_axes([0.93, 0.78 - i * 0.23 / scale - 0.3 * (1 - scale), 0.014, 0.17 / scale])
            cbar = fig.colorbar(cf_last, cax=cbar_ax, ticks=clevels)
            cbar.set_label(f"{vinfo['alias']} ({vinfo['unit']})", fontsize=fontz)
            cbar.ax.tick_params(labelsize=tick_fontz)
            for val in clevels:
                lw = 1.5 if val == 0 else 0.5
                cbar.ax.hlines(val, *cbar.ax.get_xlim(), colors='black', linewidth=lw, alpha=0.6)
            for spine in cbar.ax.spines.values():
                spine.set_visible(True)
                spine.set_linewidth(1.0)
                spine.set_edgecolor("black")

    plt.subplots_adjust(left=0.06, right=0.90, top=0.95, bottom=0.08, hspace=0.2, wspace=0.15)

    # Save all requested formats
    for fmt in output_formats:
        if fig_name is not None:
            out_file = os.path.join(fig_dir, fig_name)
        else:
            out_file = os.path.join(fig_dir, f"surface_mean_climate_{season}_{fig_key}.{fmt}")
        plt.savefig(out_file, dpi=300, bbox_inches="tight")
        print(f"[INFO] Saved: {out_file}")

    plt.show()
    plt.close()

    print(f"[INFO] Completed in {time.time() - start_time:.2f} seconds.")
    if missing_files:
        print("[WARNING] Missing files:")
        for f in missing_files:
            print(f"  - {f}")

