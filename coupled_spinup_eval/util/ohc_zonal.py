"""
Zonally integrated ocean heat content as a latitude-depth section, from the remapped (0.5°)
decadal annual climatologies of MPAS-Ocean written by the zppy-style post-processing:

    <MODEL_ROOT>/<case>/post/time_series/ocn_2d/<case>_ANN_YYYY01_YYYY12_climo.nc
    (timeMonthly_avg_activeTracers_temperature and timeMonthly_avg_layerThickness on
     (Time, nVertLevels, lat, lon), cell area in steradian)

For latitude band j and model level k

    H[j, k] = rho0 * cp * sum_lon( T * h * area * R^2 ) / (dz_k * dlat)      [J m^-1 deg^-1]

i.e. the heat content of the band per metre of (reference) depth and per degree of latitude, so
sections with uneven levels and different grids are comparable. T is potential temperature in
°C, so H is relative to 0 °C; differences (between runs or windows) are the physical quantity.
A window is the mean of the decadal climatologies covering it.
"""

import os

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm

RHO0 = 1026.0           # kg m^-3 (MPAS-Analysis OHC constants)
CP = 3996.0             # J kg^-1 K^-1
R_EARTH = 6371229.0     # m, MPAS sphere radius
TVAR = "timeMonthly_avg_activeTracers_temperature"
HVAR = "timeMonthly_avg_layerThickness"


def decadal_climo_files(case_dir, case, years, decade=10):
    """Annual decadal climatology files <case>_ANN_YYYY01_YYYY12_climo.nc covering years (y0, y1)."""
    y0, y1 = years
    if (y1 - y0 + 1) % decade:
        raise ValueError(f"window {years} is not a whole number of {decade}-yr decades")
    return [os.path.join(case_dir, f"{case}_ANN_{y:04d}01_{y + decade - 1:04d}12_climo.nc")
            for y in range(y0, y1 + 1, decade)]


class ZonalOHCSection:
    """
    Latitude-depth sections of zonally integrated heat content (see module docstring).

    sections are cached per (case, window) in cache_dir as
    ohc_zonal_<case>_<y0>-<y1>.nc and reused unless force=True.
    """

    def __init__(self, model_root, depth_file, cache_dir, *,
                 sub_dir="post/time_series/ocn_2d", rho0=RHO0, cp=CP, radius=R_EARTH,
                 engine="netcdf4"):
        self.model_root = model_root
        self.sub_dir = sub_dir
        self.cache_dir = cache_dir
        self.rho0, self.cp, self.radius = float(rho0), float(cp), float(radius)
        self.engine = engine
        with xr.open_dataset(depth_file) as d:
            bottom = d["refBottomDepth"].values.astype(float)
        top = np.concatenate([[0.0], bottom[:-1]])
        self.dz = bottom - top                       # reference layer thickness (m)
        self.depth_mid = 0.5 * (top + bottom)        # level mid-depth (m)
        self.depth_edges = np.concatenate([[0.0], bottom])

    def _one(self, path):
        """H(lat, level) of one climatology file."""
        with xr.open_dataset(path, engine=self.engine, decode_times=False) as ds:
            T = ds[TVAR].isel(Time=0)
            h = ds[HVAR].isel(Time=0)
            area = ds["area"] * self.radius ** 2                      # m^2
            heat = (T * h * area).sum("lon", skipna=True).load()     # J/(rho0 cp), (level, lat)
            dlat = float(np.abs(np.diff(ds["lat_bnds"].isel(lat=0).values))[0])
            lat = ds["lat"].values
        H = self.rho0 * self.cp * heat.values / (self.dz[:, None] * dlat)
        return xr.DataArray(H.T, dims=("lat", "depth"), coords={"lat": lat, "depth": self.depth_mid})

    def section(self, case, years, *, force=False):
        """Window-mean H(lat, depth) of `case` over years (y0, y1), from its decadal climatologies."""
        os.makedirs(self.cache_dir, exist_ok=True)
        cache = os.path.join(self.cache_dir, f"ohc_zonal_{case}_{years[0]:04d}-{years[1]:04d}.nc")
        if os.path.isfile(cache) and not force:
            with xr.open_dataset(cache) as ds:
                return ds["ohc_zonal"].load()
        files = decadal_climo_files(os.path.join(self.model_root, case, self.sub_dir), case, years)
        missing = [f for f in files if not os.path.isfile(f)]
        if missing:
            raise FileNotFoundError(f"{len(missing)} decadal climatologies missing, e.g. {missing[0]}")
        H = sum(self._one(f) for f in files) / len(files)
        H.name = "ohc_zonal"
        H.attrs = dict(units="J m-1 deg-1", case=case, years=f"{years[0]:04d}-{years[1]:04d}",
                       long_name="zonally integrated heat content per metre depth per degree latitude",
                       rho0=self.rho0, cp=self.cp, radius=self.radius,
                       source=f"{self.sub_dir}/{os.path.basename(files[0])} ... {os.path.basename(files[-1])}")
        H["depth"].attrs = dict(units="m", long_name="level mid-depth (refBottomDepth)")
        H.to_dataset().to_netcdf(cache)
        return H

    def plot_compare(self, rows, labels, *, scale=1e18, unit=r"$10^{18}$ J m$^{-1}$ deg$^{-1}$",
                     ref_levels, diff_levels, ref_cmap="RdBu_r", diff_cmap="RdBu_r",
                     ref_title="{label}: {window} $-$ {baseline}", diff_title="{label} $-$ {ref}",
                     ylim_km=(0.0, 5.5), lat_band=(-80.0, 90.0), figsize=None, fontz=14.0,
                     panel_letters=True, save=None):
        """
        One row per window: reference anomaly | target 1 - reference | target 2 - reference ...

        rows: list of dicts with keys
            window   (y0, y1) label of the row
            baseline (y0, y1) label of the reference anomaly's baseline
            ref      H of the reference over `window` minus H of the reference over `baseline`
            diffs    [H_target - H_reference over `window`, ...] (same order as labels[1:])
        labels: [reference label, target labels...]
        """
        ncols = 1 + len(labels) - 1
        nrows = len(rows)
        figsize = figsize or (6.0 * ncols, 3.6 * nrows)
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharey=True, squeeze=False,
                                 layout="constrained")
        depth_km = self.depth_edges / 1e3
        ref_norm = BoundaryNorm(ref_levels, plt.get_cmap(ref_cmap).N, extend="both")
        diff_norm = BoundaryNorm(diff_levels, plt.get_cmap(diff_cmap).N, extend="both")
        tok = lambda w: f"{w[0]:04d}–{w[1]:04d}"
        letters = iter("abcdefghijklmnopqrstuvwxyz")

        for r, row in enumerate(rows):
            panels = [(row["ref"], ref_cmap, ref_norm,
                       ref_title.format(label=labels[0], window=tok(row["window"]),
                                        baseline=tok(row["baseline"])))]
            panels += [(d, diff_cmap, diff_norm, diff_title.format(label=lab, ref=labels[0]))
                       for d, lab in zip(row["diffs"], labels[1:])]
            mesh = []
            for c, (da, cmap, norm, title) in enumerate(panels):
                ax = axes[r, c]
                da = da.sel(lat=slice(*lat_band)) if lat_band else da
                lat = da["lat"].values
                lat_edges = np.concatenate([[lat[0] - 0.25], 0.5 * (lat[1:] + lat[:-1]), [lat[-1] + 0.25]])
                Z = np.ma.masked_invalid(da.transpose("depth", "lat").values / scale)
                Z = np.ma.masked_where(Z == 0, Z)                     # levels below the bathymetry
                mesh.append(ax.pcolormesh(lat_edges, depth_km, Z, cmap=cmap, norm=norm, shading="flat",
                                          rasterized=True))      # no cell seams in the PDF
                ax.set_ylim(ylim_km[1], ylim_km[0])
                ax.set_facecolor("0.75")
                ax.set_title((f"({next(letters)}) " if panel_letters else "") + title,
                             fontsize=fontz, loc="left")
                ax.tick_params(labelsize=0.9 * fontz)
                if r == nrows - 1:
                    ax.set_xlabel("Latitude (°)", fontsize=fontz)
                if c == 0:
                    ax.set_ylabel("Depth (km)", fontsize=fontz)
            cb = fig.colorbar(mesh[0], ax=axes[r, 0], extend="both", pad=0.01)
            cb.set_label(rf"$\Delta$OHC ({unit})", fontsize=0.9 * fontz)
            cb.ax.tick_params(labelsize=0.8 * fontz)
            if len(mesh) > 1:
                cb = fig.colorbar(mesh[1], ax=list(axes[r, 1:]), extend="both", pad=0.01)
                cb.set_label(rf"$\Delta$OHC ({unit})", fontsize=0.9 * fontz)
                cb.ax.tick_params(labelsize=0.8 * fontz)
        if save:
            fig.savefig(save, dpi=300, bbox_inches="tight")
        return fig, axes
