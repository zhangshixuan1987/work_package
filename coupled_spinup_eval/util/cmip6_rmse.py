"""
Seasonal climatology RMSE of the spin-up runs in the context of the CMIP6 ensemble.

* CMIP6 benchmark: e3sm_diags seasonal RMSE of the CMIP6 historical models
  (``cmip6_historical_seasonal_rmse_<date>.csv``: header rows = field, units, season; one row per model,
  "--" = missing). The last three rows are E3SM-2-0 and its composite base/best versions.
* Runs: RMSE read from the e3sm_diags ``viewer/table-data/<season>_metrics_table.csv`` of a
  ``model_vs_obs_<years>`` run, matched by "<variable> <region> <obs>" row label.

FIELDS lists the nine fields of the benchmark with the e3sm_diags row that matches each CSV column.
"""
import os

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SEASONS = ["ANN", "DJF", "MAM", "JJA", "SON"]

# csv field name -> (e3sm_diags row label, unit label, panel title)
FIELDS = {
    "Net TOA":  ("RESTOM global ceres_ebaf_toa_v4.1", r"W m$^{-2}$",     "Net TOA radiation"),
    "SW CRE":   ("SWCF global ceres_ebaf_toa_v4.1",   r"W m$^{-2}$",     "SW cloud radiative effect"),
    "LW CRE":   ("LWCF global ceres_ebaf_toa_v4.1",   r"W m$^{-2}$",     "LW cloud radiative effect"),
    "prec":     ("PRECT global GPCP_v2.3",            r"mm day$^{-1}$",  "Precipitation"),
    "tas land": ("TREFHT land ERA5",                  "K",               "Surface air temperature (land)"),
    "SLP":      ("PSL global ERA5",                   "hPa",             "Sea-level pressure"),
    "u-200":    ("U-200mb global ERA5",               r"m s$^{-1}$",     "Zonal wind, 200 hPa"),
    "u-850":    ("U-850mb global ERA5",               r"m s$^{-1}$",     "Zonal wind, 850 hPa"),
    "Zg-500":   ("Z3-500mb global ERA5",              "hm",              "Geopotential height, 500 hPa"),
}


def read_cmip6_csv(path, n_e3sm=3):
    """CMIP6 RMSE table -> (cmip6 (model, field, season), e3sm (model, field, season)).

    The last ``n_e3sm`` rows (E3SM-2-0 and its composites) are returned separately."""
    raw = pd.read_csv(path, header=None, dtype=str)
    fields, seasons = raw.iloc[0, 1:].str.strip().values, raw.iloc[2, 1:].str.strip().values
    models = raw.iloc[3:, 0].str.strip().values
    vals = raw.iloc[3:, 1:].replace("--", np.nan).astype(float).values
    names = list(dict.fromkeys(fields))
    arr = np.full((len(models), len(names), len(SEASONS)), np.nan)
    for j, (f, s) in enumerate(zip(fields, seasons)):
        arr[:, names.index(f), SEASONS.index(s)] = vals[:, j]
    da = xr.DataArray(arr, dims=("model", "field", "season"),
                      coords={"model": models, "field": names, "season": SEASONS}, name="rmse")
    return da.isel(model=slice(0, -n_e3sm)), da.isel(model=slice(-n_e3sm, None))


def read_e3sm_diags_rmse(table_dir, fields=FIELDS):
    """RMSE (field, season) from the <season>_metrics_table.csv files in ``table_dir``."""
    arr = np.full((len(fields), len(SEASONS)), np.nan)
    for k, s in enumerate(SEASONS):
        path = os.path.join(table_dir, f"{s}_metrics_table.csv")
        if not os.path.isfile(path):
            continue
        tab = pd.read_csv(path)
        tab["Variables"] = tab["Variables"].str.strip()
        tab = tab.set_index("Variables")
        for i, (row, _, _) in enumerate(fields.values()):
            if row in tab.index:
                arr[i, k] = float(np.atleast_1d(tab.loc[row, "RMSE"])[0])
    return xr.DataArray(arr, dims=("field", "season"),
                        coords={"field": list(fields), "season": SEASONS}, name="rmse")


def plot_rmse_boxes(cmip6, runs, *, reference=None, fields=FIELDS, ncols=3, figsize=None,
                    fontz=12.0, panel_letters=True, save=None):
    """
    One panel per field: CMIP6 spread per season (box = 25-75 %, whiskers = min-max, line = median,
    dot = mean), with one marker per run.

    runs: {label: dict(rmse=DataArray (field, season), color=..., marker=...)}
    reference: optional {label: dict(rmse=..., color=..., marker=...)} drawn like runs (e.g. E3SM-2-0)
    """
    allruns = {**(reference or {}), **runs}
    names = list(fields)
    nrows = int(np.ceil(len(names) / ncols))
    figsize = figsize or (5.0 * ncols, 3.4 * nrows)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False, layout="constrained")
    x = np.arange(len(SEASONS))
    nr = len(allruns)
    offs = (np.arange(nr) - (nr - 1) / 2) * 0.12
    letters = iter("abcdefghijklmnopqrstuvwxyz")
    for ax, f in zip(axes.flat, names):
        data = [cmip6.sel(field=f, season=s).dropna("model").values for s in SEASONS]
        ax.boxplot(data, positions=x, widths=0.5, whis=(0, 100), showfliers=False,
                   medianprops=dict(color="k"), boxprops=dict(color="0.3"), whiskerprops=dict(color="0.3", ls="--"))
        ax.plot(x, [np.mean(d) for d in data], "o", color="k", ms=6, zorder=3)
        for (lab, r), dx in zip(allruns.items(), offs):
            ax.plot(x + dx, r["rmse"].sel(field=f).values, linestyle="none", marker=r.get("marker", "^"),
                    mfc="none", mec=r["color"], mew=1.8, ms=8, zorder=4)
        _, unit, title = fields[f]
        ax.set_title((f"({next(letters)}) " if panel_letters else "") + f"{title} ({unit})",
                     fontsize=fontz, loc="left")
        ax.set_xticks(x, SEASONS)
        ax.tick_params(labelsize=0.9 * fontz)
        ax.grid(axis="y", alpha=0.3)
    for ax in axes.flat[len(names):]:
        ax.set_visible(False)
    handles = [Line2D([], [], color="k", marker="o", ls="none", ms=6, label=f"CMIP6 mean ({cmip6.sizes['model']} models)")]
    handles += [Line2D([], [], ls="none", marker=r.get("marker", "^"), mfc="none", mec=r["color"], mew=1.8,
                       ms=8, label=lab) for lab, r in allruns.items()]
    fig.legend(handles=handles, loc="outside lower center", ncol=min(len(handles), 5), fontsize=fontz, frameon=False)
    if save:
        fig.savefig(save, dpi=300, bbox_inches="tight")
    return fig, axes
