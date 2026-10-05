"""
Ocean surface freshwater-flux budget by source, from the MPAS-Ocean global statistics written to
``post/time_series/ocn/<field>_YYYYMM-YYYYMM.nc``.

The ``accumulated*Flux`` and ``netMassFlux`` series are global integrals in kg s-1 (positive into the ocean;
the ``kg m^-2 s^-1`` units attribute of the accumulated fields is wrong). The net flux is the sum of the
sources: ``netMassFlux`` itself is 0 in the first month after every restart (accumulator reset), which biases
the annual mean of restart years by ~10 mm yr-1; it is kept as ``net_mpas`` for reference.
They are converted to ocean-area means in mm yr-1 with the global ocean area (MPAS-Analysis
``timeMonthly_avg_areaCellGlobal``).
"""
import glob
import os
import re

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from util.area_mean_drift import DPM

# component -> (MPAS field, label)
COMPONENTS = {
    "rain":   ("accumulatedRainFlux",        "Rain"),
    "snow":   ("accumulatedSnowFlux",        "Snow"),
    "evap":   ("accumulatedEvaporationFlux", "Evaporation"),
    "river":  ("accumulatedRiverRunoffFlux", "River runoff"),
    "irun":   ("accumulatedIceRunoffFlux",   "Ice runoff"),
    "seaice": ("accumulatedSeaIceFlux",      "Sea-ice melt/growth"),
    "frazil": ("accumulatedFrazilFlux",      "Frazil"),
    "landice": ("accumulatedLandIceFlux",    "Land-ice melt"),
    "berg":   ("accumulatedIcebergFlux",     "Iceberg melt"),
}
NET_MPAS = ("netMassFlux", "Net (MPAS netMassFlux)")
SEC_PER_YEAR = 86400.0 * 365.0
_FNAME = re.compile(r"_(\d{4})01-(\d{4})12\.nc$")


def ocean_area(mpas_ts_file):
    """Global ocean area (m2) from an MPAS-Analysis mpasTimeSeriesOcean.nc."""
    with xr.open_dataset(mpas_ts_file, decode_times=False) as ds:
        return float(ds["timeMonthly_avg_areaCellGlobal"].isel(Time=0))


def read_fw_annual(ts_dir, years, area_m2, components=COMPONENTS, with_mpas_net=True):
    """Annual (days-in-month weighted) ocean-area-mean freshwater fluxes (mm yr-1), Dataset (year,).

    Years without all 12 months are NaN."""
    y0, y1 = years
    yrs = np.arange(y0, y1 + 1)
    out = {}
    w = DPM / DPM.sum()
    fields = dict(components, **({"net_mpas": NET_MPAS} if with_mpas_net else {}))
    for key, (field, label) in fields.items():
        mon = np.full((yrs.size, 12), np.nan)
        for f in sorted(glob.glob(os.path.join(ts_dir, f"{field}_*.nc"))):
            m = _FNAME.search(f)
            if not m or int(m[2]) < y0 or int(m[1]) > y1:
                continue
            with xr.open_dataset(f, decode_times=False) as ds:
                vals = ds[field].values.astype(float)
            fy = int(m[1])
            for i, v in enumerate(vals):
                yr, mo = fy + i // 12, i % 12
                if y0 <= yr <= y1:
                    mon[yr - y0, mo] = v
        ann = (mon * w).sum(1)
        ann[np.isnan(mon).any(1)] = np.nan
        out[key] = xr.DataArray(ann / area_m2 * SEC_PER_YEAR, dims="year", coords={"year": yrs},
                                attrs=dict(units="mm yr-1", long_name=label, source=field))
    ds = xr.Dataset(out)
    ds["net"] = sum(ds[k] for k in components)          # sources only (not net_mpas)
    ds["net"].attrs = dict(units="mm yr-1", long_name="Net (sum of sources)")
    return ds


def running_mean(da, n):
    """Centred running mean over n years (NaN where the window is incomplete)."""
    return da.rolling(year=n, center=True, min_periods=n).mean() if n and n > 1 else da


def plot_fw_components(series, components, *, smooth=10, segments=None, markers=(), ncols=2, figsize=None,
                       fontz=12.0, panel_letters=True, save=None):
    """
    One panel per component: annual flux (mm yr-1) of each run, ``smooth``-yr running mean.

    series:   {label: dict(ds=read_fw_annual() Dataset, color=...)}
    segments: optional {label: [(y0, y1), ...]} shaded per run (e.g. its fully coupled segments)
    """
    names = list(components)
    nrows = int(np.ceil(len(names) / ncols))
    figsize = figsize or (7.0 * ncols, 2.6 * nrows)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False, sharex=True, layout="constrained")
    letters = iter("abcdefghijklmnopqrstuvwxyz")
    for ax, k in zip(axes.flat, names):
        for lab, r in series.items():
            ax.plot(r["ds"]["year"], running_mean(r["ds"][k], smooth), color=r["color"], lw=1.3, label=lab)
            for y0, y1 in (segments or {}).get(lab, []):
                ax.axvspan(y0, y1, color=r["color"], alpha=0.08, lw=0)
        for yr in markers:
            ax.axvline(yr, color="0.4", ls=":", lw=1)
        ax.set_title((f"({next(letters)}) " if panel_letters else "") + components[k][1], fontsize=fontz, loc="left")
        ax.set_ylabel(r"mm yr$^{-1}$", fontsize=0.9 * fontz)
        ax.tick_params(labelsize=0.85 * fontz)
        ax.grid(alpha=0.3)
    for ax in axes.flat[len(names):]:
        ax.set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel("Model year", fontsize=0.9 * fontz)
    h, l = axes.flat[0].get_legend_handles_labels()
    fig.legend(h, l, loc="outside lower center", ncol=len(l), fontsize=fontz, frameon=False)
    if save:
        fig.savefig(save, dpi=300, bbox_inches="tight")
    return fig, axes


def plot_fw_differences(diffs, components, *, colors, fontz=12.0, figsize=(12, 4), save=None):
    """Grouped bars: mean difference from the reference per component and period.

    diffs: {period label: {run label: Dataset of mean differences}}; colors: {run label: color}"""
    names = list(components)
    periods = list(diffs)
    runs = list(next(iter(diffs.values())))
    fig, axes = plt.subplots(1, len(periods), figsize=figsize, sharey=True, squeeze=False, layout="constrained")
    x = np.arange(len(names))
    wd = 0.8 / len(runs)
    for ax, p in zip(axes[0], periods):
        for j, r in enumerate(runs):
            ax.bar(x + (j - (len(runs) - 1) / 2) * wd, [float(diffs[p][r][k]) for k in names], wd,
                   color=colors[r], label=r)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks(x, [components[k][1] for k in names], rotation=40, ha="right", fontsize=0.85 * fontz)
        ax.set_title(p, fontsize=fontz, loc="left")
        ax.grid(axis="y", alpha=0.3)
    axes[0, 0].set_ylabel(r"difference (mm yr$^{-1}$)", fontsize=0.9 * fontz)
    axes[0, -1].legend(fontsize=0.85 * fontz, frameon=False)
    if save:
        fig.savefig(save, dpi=300, bbox_inches="tight")
    return fig, axes
