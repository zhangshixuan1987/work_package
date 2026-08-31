import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import ListedColormap, BoundaryNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.ndimage import gaussian_filter
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error
from typing import Tuple, Optional

def define_region(
    regnam: str = 'global',
    custom_lat: Optional[Tuple[float, float]] = None,
    custom_lon: Optional[Tuple[float, float]] = None,
    use_360_lon: bool = False
) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    """
    Return latitude and longitude bounds for a named or custom region.

    Parameters
    ----------
    regnam : str
        Name of the region (e.g., 'global', 'NHMidLat', 'Tropics').
    custom_lat : Optional[Tuple[float, float]]
        Custom latitude bounds if regnam is 'custom'.
    custom_lon : Optional[Tuple[float, float]]
        Custom longitude bounds if regnam is 'custom'.
    use_360_lon : bool
        If True, convert lon range to [0, 360] convention.

    Returns
    -------
    Tuple[Tuple[float, float], Tuple[float, float]]
        A tuple containing (lat_min, lat_max) and (lon_min, lon_max).
    """
    reg_dict = {
        'global': [(-90, 90), (-180, 180)],
        'nhmidlat': [(25, 50), (-60, 150)],
        'tropics': [(-10, 10), (-90, 60)],
        'atlantic': [(5, 55), (-95, -40)],
        'conus': [(25, 50), (-125, -95)],
        'antarctic': [(-90, -50), (-180, 180)],
        'polarn': [(50, 90), (-180, 180)],
        'greenland': [(60, 85), (-75, -10)]
    }

    regnam_lower = regnam.lower()
    if regnam_lower == "custom":
        if custom_lat is None or custom_lon is None:
            raise ValueError("For 'custom' region, you must provide both `custom_lat` and `custom_lon`.")
        lat_bounds = custom_lat
        lon_bounds = custom_lon
    elif regnam_lower in reg_dict:
        lat_bounds, lon_bounds = reg_dict[regnam_lower]
    else:
        raise ValueError(f"Region '{regnam}' not defined. Available: {list(reg_dict.keys()) + ['custom']}")

    if use_360_lon:
        lon_bounds = tuple((lon + 360) if lon < 0 else lon for lon in lon_bounds)

    return lat_bounds, lon_bounds

def format_coord_labels(coords, axis="lat"):
    labels = []
    for c in coords:
        val = int(round(c))
        if axis == "lat":
            if val > 0:
                labels.append(f"{val}°N")
            elif val < 0:
                labels.append(f"{abs(val)}°S")
            else:
                labels.append("0°")
        elif axis == "lon":
            val = val % 360
            if val == 0 or val == 360:
                labels.append("0°")
            elif val <= 180:
                labels.append(f"{val}°E")
            else:
                labels.append(f"{360 - val}°W")
        else:
            labels.append(str(val))
    return labels


def plot_multi_model_hovmoller_from_files(
    data_dir,
    model_dict,
    varname,
    model_to_varname=None,
    reg_dict=None,
    process_year=2011,
    hovmoller_axis="lat",
    season_suffix="ALL",
    base_cmap="YlGnBu",
    norm=None,
    n_skip=2,
    save_path=None,
    smooth_sigma=(1, 1),
    panel_width=5.0,
    panel_height=3.5,
    resample_to_daily=True,
    reverse_x=False,
    reverse_y=False,
    show_plot=False,
    fig_name=None,
):
    regnams = list(reg_dict.keys())
    model_list = list(model_dict.keys())
    nrows, ncols = len(regnams), len(model_list)
    fontz = int(24 * panel_width / 5.0)
    axis_suffix = f"time-{hovmoller_axis}"

    width_ratios = [1.0] * (ncols - 1) + [1.2]
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(panel_width * sum(width_ratios), panel_height * nrows),
        sharex=False, sharey=False,
        gridspec_kw={"width_ratios": width_ratios, "hspace": 0.3, "wspace": 0.15}
    )
    plt.subplots_adjust(left=0.07, right=0.93)
    axes = np.atleast_2d(axes)

    for i_r, regnam in enumerate(regnams):
        reg_info = reg_dict[regnam]
        clevels = reg_info['clevels']
        cmap_name = base_cmap if isinstance(base_cmap, str) else base_cmap.name
        base = plt.get_cmap(base_cmap)  # e.g., cmocean.cm.CET_R1
        nlevs = len(clevels) - 1

        # Create finely resolved color values and select middle portion
        n_total = nlevs + 2 * n_skip
        full_vals = np.linspace(0, 1, n_total)
        main_colors = base(full_vals[n_skip:-n_skip])
        under_color = base(full_vals[n_skip - 1])
        over_color = base(full_vals[-(n_skip)])
    
        # Create and customize colormap
        cmap = ListedColormap(main_colors)
        cmap.set_under(under_color)
        cmap.set_over(over_color)
    
        norm = BoundaryNorm(clevels, ncolors=cmap.N)

        # Retrieve plot info
        reg_label = reg_info['reg_str']
        lat_range, lon_range = reg_info['region']
        xticks = reg_info['xticks']

        ref_name = next((m for m in model_list if model_to_varname and model_to_varname.get(m) == "reference"), None)
        if not ref_name:
            print(f"[ERROR] No reference found for {regnam}.")
            continue

        ref_path = os.path.join(data_dir, f"hovmoller_{ref_name}_{regnam}_{varname}_{process_year}_{axis_suffix}.nc")
        if not os.path.exists(ref_path):
            print(f"[ERROR] Missing reference file: {ref_path}")
            continue

        ref_var = model_to_varname.get(ref_name, 'reference')
        ref_ds = xr.open_dataset(ref_path)
        da_ref = ref_ds[ref_var].resample(time="1D").mean() if resample_to_daily else ref_ds[ref_var]

        if hovmoller_axis == "lon" and (da_ref[hovmoller_axis] < 0).any():
            da_ref = da_ref.assign_coords({hovmoller_axis: da_ref[hovmoller_axis] % 360})
            da_ref = da_ref.sortby(hovmoller_axis)

        for j_m, model in enumerate(model_list):
            ax = axes[i_r, j_m]
            model_path = os.path.join(data_dir, f"hovmoller_{model}_{regnam}_{varname}_{process_year}_{axis_suffix}.nc")
            var_name = model_to_varname.get(model, 'model')

            if not os.path.exists(model_path):
                print(f"[SKIP] Missing file: {model_path}")
                ax.set_visible(False)
                continue

            ds = xr.open_dataset(model_path)
            if var_name not in ds:
                print(f"[SKIP] Variable '{var_name}' missing in {model_path}")
                ax.set_visible(False)
                continue

            da = ds[var_name].resample(time="1D").mean() if resample_to_daily else ds[var_name]
            space, time = da.coords[hovmoller_axis], da.coords["time"]

            if hovmoller_axis == "lon" and (space < 0).any():
                space = space % 360
                da = da.assign_coords({hovmoller_axis: space})
                da = da.sortby(hovmoller_axis)

            mod_data = da.values
            ref_data = da_ref.sel({hovmoller_axis: space, "time": time}).values

            # Compute metrics
            x, y = ref_data.flatten(), mod_data.flatten()
            mask = np.isfinite(x) & np.isfinite(y)
            pcc = pearsonr(x[mask], y[mask])[0] if np.any(mask) else np.nan
            rmse = np.sqrt(mean_squared_error(x[mask], y[mask])) if np.any(mask) else np.nan

            # Smoothing
            smoothed = gaussian_filter(mod_data, sigma=smooth_sigma)
            data = smoothed
            time_plot = mdates.date2num(time.values)

            dx = (space.values[-1] - space.values[0]) / (len(space) - 1)
            dy = (time_plot[-1] - time_plot[0]) / (len(time_plot) - 1)

            x0 = float(space.values[0]) - dx / 2
            x1 = float(space.values[-1]) + dx / 2
            y0 = time_plot[0] - dy / 2
            y1 = time_plot[-1] + dy / 2

            if reverse_x:
                x0, x1 = x1, x0
                space = space[::-1]
                data = data[:, ::-1]

            if reverse_y:
                y0, y1 = y1, y0
                time_plot = time_plot[::-1]
                data = data[::-1, :]

            extent = [x0, x1, y0, y1]

            #ax.imshow(
            #    data,
            #    aspect="auto",
            #    cmap=cmap,
            #    norm=norm,
            #    vmin=clevels[0],
            #    vmax=clevels[-1],
            #    extent=extent,
            #    origin="lower"
            #)

            ax.imshow(
                data,
                aspect="auto",
                cmap=cmap,
                norm=norm,             # ✅ norm already contains vmin/vmax info
                extent=extent,
                origin="lower"
            )

            panel_id = f"{chr(97 + i_r * ncols + j_m)})"
            ax.set_title(f"{panel_id} {model_dict[model]}", fontsize=fontz, loc="left")
            ax.set_title(f"{reg_label}", fontsize=fontz, loc="right")

            if j_m == 0:
                ax.set_ylabel("Month", fontsize=fontz)
                ax.yaxis.set_major_locator(mdates.MonthLocator())
                ax.yaxis.set_major_formatter(mdates.DateFormatter('%b'))
                ax.tick_params(axis='y', labelsize=fontz)
            else:
                ax.set_yticklabels([])

            if hovmoller_axis in ("lat", "lon"):
                if hovmoller_axis == "lon":
                    xticks_vals = [x % 360 for x in xticks]
                    xlim_vals = [lon_range[0] % 360, lon_range[1] % 360]
                else:
                    xticks_vals = xticks
                    xlim_vals = lat_range

                if reverse_x:
                    xlim_vals = xlim_vals[::-1]
                    xticks_vals = xticks_vals[::-1]

                ax.set_xticks(xticks_vals)
                ax.set_xticklabels(format_coord_labels(xticks_vals, hovmoller_axis), fontsize=int(fontz * 0.95))
                ax.set_xlim(xlim_vals)
                ax.set_xlabel(hovmoller_axis.capitalize(), fontsize=fontz)

            if model != ref_name:
                ax.text(
                    0.98, 0.02,
                    f"PCC={pcc:.2f}\nRMSE={rmse:.2f}",
                    transform=ax.transAxes,
                    fontsize=fontz * 0.9,
                    va="bottom",
                    ha="right",
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8)
                )

        # === Colorbar per row ===
        last_ax = axes[i_r, -1]
        divider = make_axes_locatable(last_ax)
        cax = divider.append_axes("right", size="7%", pad=0.1)
        cbar = fig.colorbar(last_ax.images[0], cax=cax, extend="both")
        units = getattr(da, "attrs", {}).get("units", "")
        cbar.set_label(f"{varname} ({units})", fontsize=fontz)
        cbar.ax.tick_params(labelsize=int(fontz * 0.95))

    # === Save and/or show ===
    if save_path:
        os.makedirs(save_path, exist_ok=True)
        if fig_name is not None:
            out_file = os.path.join(save_path, fig_name)
        else:
            out_file = os.path.join(save_path, f"hovmoller_MULTI_REGIONS_{varname}_ALL_{process_year}.pdf")
        plt.savefig(out_file, bbox_inches="tight", dpi=600)
        print(f"[SAVED] {out_file}")

    if show_plot:
        plt.show()

    plt.close()
