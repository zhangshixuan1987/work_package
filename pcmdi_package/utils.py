import glob
import os
import re
from datetime import datetime
from typing import Dict, List, Tuple
from matplotlib.transforms import Bbox
from matplotlib import colorbar as mpl_colorbar

import pandas as pd

from logger import _setup_child_logger

logger = _setup_child_logger(__name__)


def find_latest_file_list(
    path: str,
    file_pattern: str,
    var_pattern=r"\.(\w+)\.\d{8}\.nc$",
    time_pattern=r"\.(\d{8})\.nc$",
) -> List[str]:
    """
    Find the latest NetCDF file for each variable in the directory based on timestamps in filenames.

    Args:
        path (str): Directory to search.
        file_pattern (str): Regex to search file lists.
        var_pattern (str): Regex to extract variable name.
        time_pattern (str): Regex to extract date.

    Returns:
        List[str]: List of file paths, one for each variable (latest by timestamp).
    """
    latest_files: Dict[str, Tuple[datetime, str]] = {}
    files = glob.glob(os.path.join(path, file_pattern))
    if not files:
        """
        FAILURE

        No files found in /lcrc/group/e3sm/public_html/diagnostic_output/ac.forsyth2/zppy_pr719_output/unique_id_21/v3.LR.amip_0101/pcmdi_diags/model_vs_obs/metrics_data/variability_modes/*/* that match pattern: var_mode_*.json

        ls /lcrc/group/e3sm/public_html/diagnostic_output/ac.forsyth2/zppy_pr_719_output/unique_id_21/v3.LR.amip_0101/pcmdi_diags/model_vs_obs/metrics_data/variability_modes/
        AMO  NAM  NAO  NPGO  NPO  PDO  PNA  PSA1  PSA2  SAM

        ls /lcrc/group/e3sm/public_html/diagnostic_output/ac.forsyth2/zppy_pr719_output/unique_id_21/v3.LR.amip_0101/pcmdi_diags/model_vs_obs/metrics_data/variability_modes/AMO/HadISST2/
        AMO_ts_EOF1_monthly_obs_1869-2014.nc  AMO_ts_EOF1_yearly_obs_1869-2014.nc

        SYNTHETIC PLOTS ERROR #2: No json files produced by variability modes, even though those jobs completed successfully!
        """
        logger.error(f"No files found in {path} that match pattern: {file_pattern}")
    for f in files:
        fname = os.path.basename(f)
        var_match = re.search(var_pattern, fname)
        time_match = re.search(time_pattern, fname)

        if var_match and time_match:
            logger.info(f"{fname} matched var and time patterns")
            var = var_match.group(1)
            try:
                timestamp = datetime.strptime(time_match.group(1), "%Y%m%d")
            except ValueError:
                continue

            if var not in latest_files or timestamp > latest_files[var][0]:
                latest_files[var] = (timestamp, f)
        else:
            logger.warning(f"{fname} failed to match both var and time patterns")

    return [file for _, file in latest_files.values()]


def get_highlight_models(all_models, model_name):
    """
    Prioritize models containing 'e3sm' and then any additional specified models.

    Parameters:
        data_dict (dict): Dictionary with a 'model' key containing a list of model names.
        model_name (list): List of models to also highlight (after e3sm models).

    Returns:
        list: Ordered list of unique models to highlight.
    """
    highlight_model1 = []

    # First, collect all models that contain "e3sm" (case-insensitive)
    e3sm_models = [m for m in all_models if "e3sm" in m.lower()]

    # Then collect models in model_name that are not already in e3sm_models
    additional_models = [
        m for m in all_models if m in model_name and m not in e3sm_models
    ]

    # Combine both lists
    highlight_model1 = e3sm_models + additional_models

    return highlight_model1


def shift_row_to_bottom(df, index_to_shift):
    """
    Moves the specified row to the bottom of the DataFrame and resets the index.

    Parameters:
        df (pd.DataFrame): The input DataFrame.
        index_to_shift (int): The index of the row to move to the bottom.

    Returns:
        pd.DataFrame: A new DataFrame with the row moved to the bottom and index reset.
    """
    if index_to_shift not in df.index:
        raise IndexError(f"Index {index_to_shift} not found in DataFrame.")

    df_top = df.drop(index=index_to_shift)
    df_bottom = df.loc[[index_to_shift]]

    new_df = pd.concat([df_top, df_bottom], ignore_index=True)
    return new_df

def drop_vars(data_dict, var_names, var_units=None):
    """
    Drop variables (columns) from data_dict where more than 90% of the values are NaN.

    Parameters:
        data_dict (pd.DataFrame): Data containing variable columns.
        var_names (list): List of variable names matching data_dict columns.
        var_units (list, optional): List of units for variables. Must match var_names in order.

    Returns:
        Tuple of (filtered_data_dict, updated_var_names, updated_var_units)
    """
    protected_columns = {"model", "run", "model_run", "num_runs"}
    columns_to_drop = []

    for column in data_dict.columns:
        if column in protected_columns:
            continue
        nan_ratio = data_dict[column].isna().mean()
        if nan_ratio > 0.9:
            columns_to_drop.append(column)

    # Drop columns from DataFrame
    data_dict = data_dict.drop(columns=columns_to_drop)

    # Update var_names and var_units if applicable
    updated_var_names = [v for v in var_names if v not in columns_to_drop]
    updated_var_units = None
    if var_units is not None:
        # Keep units only for remaining variables
        name_to_unit = dict(zip(var_names, var_units))
        updated_var_units = [
            name_to_unit[v] for v in updated_var_names if v in name_to_unit
        ]

    return data_dict, updated_var_names, updated_var_units

def archive_data(
    region, stat, season, data_dict, model_name, var_names, var_units, outdir
):
    """
    Archive processed data into a CSV file with variable units in column headers if available.

    Parameters:
        region (str): Region name.
        stat (str): Statistic type (e.g., mean, std).
        season (str): Season name.
        data_dict (dict or DataFrame): Data to archive.
        model_name (str): Model identifier.
        var_names (list): List of variable names.
        var_units (list): List of variable units (optional, same order as var_names).
        outdir (str): Directory to save the CSV file.
    """
    df = pd.DataFrame(data_dict)

    # Determine the index of the first variable column (assumes first 3 are metadata)
    metadata_cols = df.columns[:3].tolist()
    variable_cols = df.columns[3:]

    filtered_cols = []
    new_column_names = df.columns.tolist()

    for var in variable_cols:
        if var in var_names:
            filtered_cols.append(var)
            if var_units:
                idx = df.columns.get_loc(var)
                unit_label = var_units[var_names.index(var)]
                new_column_names[idx] = f"{var} ({unit_label})"

    # Subset dataframe and rename columns if units provided
    df = df[metadata_cols + filtered_cols]
    df.columns = new_column_names[: len(df.columns)]

    # Ensure output directory exists
    os.makedirs(outdir, exist_ok=True)

    # Construct and save the output filename
    outfile = f"{stat}_{region}_{season}_{model_name}.csv"
    df.to_csv(os.path.join(outdir, outfile), index=False)

    return

def realign_cbar_and_legend(
    fig, ax, cbar,
    cbar_width_in=0.18,
    gap_in=0.1,
    right_pad_in=0.12,
    y_gap_above_cbar_in=None,
    top_guard=0.98,
    debug=False,
    cbar_side="right",          # "right" or "left"
    min_buffer_frac=0.002,      # tiny buffer so cbar never touches axes
    auto_nudge=True,            # derive position from tight bbox of axes
    label_clearance_in=0.06,    # inches between labels and colorbar
):
    """
    Align colorbar with ax height; optionally rebuild it on the LEFT; then place
    legend/season-diamond tight above the colorbar.

    If auto_nudge is True, position the colorbar just outside the rendered
    x-tick labels (from ax.get_tightbbox) plus `label_clearance_in`.
    """
    fig.canvas.draw()
    fig_w_in, fig_h_in = fig.get_figwidth(), fig.get_figheight()
    ax_box = ax.get_position()  # figure-fraction coords

    # convert lengths to figure fraction
    if y_gap_above_cbar_in is None:
        y_gap_above_cbar_in = 0.02
    y_gap = max(0.0, y_gap_above_cbar_in) / fig_h_in

    cbar_w = cbar_width_in / fig_w_in
    gap    = gap_in / fig_w_in
    rpad   = right_pad_in / fig_w_in
    buf    = max(0.0, float(min_buffer_frac))
    clearance = max(0.0, label_clearance_in) / fig_w_in

    # tight bbox of axes (includes tick labels), in figure coords
    tight_fig = None
    try:
        bb_disp = ax.get_tightbbox(fig.canvas.get_renderer())
        if bb_disp is not None:
            tight_fig = bb_disp.transformed(fig.transFigure.inverted())
    except Exception:
        tight_fig = None

    # --- 1) Place (or rebuild) the colorbar ---
    if cbar_side == "left":
        # rebuild cbar on left
        mappable = cbar.mappable
        label    = cbar.ax.get_ylabel() if getattr(cbar, "ax", None) is not None else ""
        ticks    = [t for t in (cbar.get_ticks() if hasattr(cbar, "get_ticks") else [])]
        try:
            if getattr(cbar, "ax", None) is not None:
                fig.delaxes(cbar.ax)
        except Exception:
            pass

        # desired right edge of cbar
        if auto_nudge and (tight_fig is not None):
            target_x1 = min(ax_box.x0 - buf, tight_fig.x0 - clearance)
        else:
            target_x1 = ax_box.x0 - max(gap, buf)

        # clamp inside figure [0,1]
        cbar_x1 = max(cbar_w + buf, min(1.0 - rpad, target_x1))
        cbar_x0 = cbar_x1 - cbar_w

        new_cax = fig.add_axes([cbar_x0, ax_box.y0, cbar_w, ax_box.height])
        cbar = mpl_colorbar.Colorbar(new_cax, mappable=mappable, orientation="vertical")
        if label:
            cbar.set_label(label)
        if ticks:
            cbar.set_ticks(ticks)
        new_cax.yaxis.set_ticks_position("left")
        new_cax.yaxis.set_label_position("left")

        needed_right = ax_box.x1 + rpad
        needed_left  = cbar_x0 - rpad

    else:  # right
        if auto_nudge and (tight_fig is not None):
            target_x0 = max(ax_box.x1 + buf, tight_fig.x1 + clearance)
        else:
            target_x0 = ax_box.x1 + max(gap, buf)

        # clamp inside figure [0,1]
        cbar_x0 = max(buf, min(1.0 - rpad - cbar_w, target_x0))
        cbar_x1 = cbar_x0 + cbar_w

        cbar.ax.set_position([cbar_x0, ax_box.y0, cbar_w, ax_box.height])
        needed_right = cbar_x1 + rpad
        needed_left  = None

    # --- 2) Legend / seasonal diamond just above the colorbar ---
    leg = ax.get_legend()
    if leg is not None:
        leg.set_loc("lower center")
        leg.set_bbox_to_anchor((0.5, 1.0 + y_gap), transform=cbar.ax.transAxes)
    else:
        # find the small square-ish inset
        candidates = []
        for a in fig.axes:
            if a is ax or a is cbar.ax:
                continue
            b = a.get_position()
            w, h = b.width, b.height
            if 0.0015 <= (w * h) <= 0.04 and 0.6 <= (w / h) <= 1.4:
                candidates.append(a)
        if candidates:
            inset = min(candidates, key=lambda a: a.get_position().width * a.get_position().height)
            ib = inset.get_position()
            w, h = ib.width, ib.height
            cb = cbar.ax.get_position()
            cx = 0.5 * (cb.x0 + cb.x1)
            x0 = max(0.0, min(1.0 - w, cx - w / 2.0))
            y0 = min(0.995 - h, max(0.0, cb.y1 + y_gap))
            inset.set_position([x0, y0, w, h])

    # --- 3) Adjust margins if needed ---
    if cbar_side == "right":
        if needed_right > 1.0:
            new_right = max(0.55, 1.0 - (needed_right - 1.0))
            fig.subplots_adjust(right=new_right)
    else:
        if (needed_left is not None) and (needed_left < 0.0):
            grow = -needed_left
            new_left = min(0.45, fig.subplotpars.left + grow)
            fig.subplots_adjust(left=new_left)

    fig.subplots_adjust(
        top=min(top_guard, fig.subplotpars.top),
        bottom=max(0.04, fig.subplotpars.bottom),
    )
    fig.canvas.draw_idle()

    return cbar  # <-- important when cbar_side="left"

