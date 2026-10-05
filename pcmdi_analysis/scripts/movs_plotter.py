import os
import math
from typing import Dict, Sequence, Optional, Tuple, Any, List

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter


def add_sig_dots(
    ax,
    pval: xr.DataArray,
    sig_level: float = 0.05,
    dot_color: str = "k",
    dot_size: float = 0.8,
    dot_density: int = 2,
    transform=None,
) -> None:
    """Overlay stippling dots on *ax* where *pval* < *sig_level*.

    Parameters
    ----------
    ax : cartopy GeoAxes (or any matplotlib Axes)
    pval : xr.DataArray (lat, lon)
        Raw two-tailed p-value field saved by ModeAnalyzer
        (``ds["slope_pval"]`` or ``ds["slope_pval_proj"]``).
    sig_level : float
        Significance threshold (e.g. 0.05 or 0.10).
    dot_color : str
        Marker colour.
    dot_size : float
        Marker size in points.
    dot_density : int
        Stride for sub-sampling the grid (1 = every point, 2 = every other, …).
        Increase to reduce dot density on fine grids.
    transform : cartopy CRS, optional
        Coordinate transform passed to ``ax.scatter``.
        Defaults to ``ccrs.PlateCarree()``.
    """
    if transform is None:
        transform = ccrs.PlateCarree()

    sig_mask = pval.values < sig_level          # (lat, lon) bool
    lat = pval.lat.values
    lon = pval.lon.values

    # Sub-sample for readability
    lat_idx = np.arange(0, len(lat), dot_density)
    lon_idx = np.arange(0, len(lon), dot_density)
    mask_sub = sig_mask[np.ix_(lat_idx, lon_idx)]
    lon2d, lat2d = np.meshgrid(lon[lon_idx], lat[lat_idx])

    ax.scatter(
        lon2d[mask_sub],
        lat2d[mask_sub],
        s=dot_size,
        c=dot_color,
        marker=".",
        linewidths=0,
        transform=transform,
        zorder=5,
    )


class ExtrapropicalModeMapPlotter:
    """
    Plotter for PMP extratropical modes of variability.

    Two plot families:

    1) Multi-product grid figures (one row per product, one column per group):
       - plot_mode_season_maps(...)

    2) Multi-panel "Obs + all models" with per-panel r / RMSE annotations:
       - plot_multimodel_mode_pattern_with_stats(...)   # regional (mode bounds)
       - plot_multimodel_teleconnection_with_stats(...) # global
    """

    def __init__(
        self,
        fig_dir: str,
        plot_dict: Dict[str, dict],
        group_order: Sequence[str] = ("hist", "future"),
        obs_key: str = "reference",
        lat_name: str = "latitude",
        lon_name: str = "longitude_a",
    ):
        self.fig_dir = fig_dir
        self.plot_dict = plot_dict
        self.group_order = tuple(group_order)
        self.obs_key = obs_key
        self.lat_name = lat_name
        self.lon_name = lon_name
        os.makedirs(self.fig_dir, exist_ok=True)

    # ==========================================================
    # Region helpers (wrap-safe lon subset)
    # ==========================================================
    def _normalize_lon_for_bounds(
        self,
        da: xr.DataArray,
        lon_convention: str,
    ) -> xr.DataArray:
        """
        Normalize lon coordinate to a convention and sort to monotonic.
        - "negpos": [-180, 180)
        - "0_360": [0, 360)
        """
        lon = da[self.lon_name]
        if lon_convention == "negpos":
            lon2 = ((lon + 180.0) % 360.0) - 180.0
            da = da.assign_coords({self.lon_name: lon2}).sortby(self.lon_name)
        elif lon_convention == "0_360":
            lon2 = lon % 360.0
            da = da.assign_coords({self.lon_name: lon2}).sortby(self.lon_name)
        else:
            raise ValueError(f"Unknown lon_convention={lon_convention!r}")
        return da

    def _subset_latlon(
        self,
        da: xr.DataArray,
        lat_bnds: Optional[Tuple[float, float]],
        lon_bnds: Optional[Tuple[float, float]],
        lon_convention: str,
    ) -> xr.DataArray:
        """
        Subset da on lat/lon bounds with wrap handling.
        For negative/positive lon bounds (e.g., -80..40), use lon_convention="negpos".
        """
        da2 = da

        if lon_bnds is not None:
            da2 = self._normalize_lon_for_bounds(da2, lon_convention)

        if lat_bnds is not None:
            la0, la1 = lat_bnds
            lo, hi = (la0, la1) if la0 <= la1 else (la1, la0)
            da2 = da2.sel({self.lat_name: slice(lo, hi)})

        if lon_bnds is not None:
            lo0, lo1 = lon_bnds
            if lo0 <= lo1:
                da2 = da2.sel({self.lon_name: slice(lo0, lo1)})
            else:
                # wrap (e.g. 120 -> -120 in negpos)
                if lon_convention == "negpos":
                    a = da2.sel({self.lon_name: slice(lo0, 180)})
                    b = da2.sel({self.lon_name: slice(-180, lo1)})
                    da2 = xr.concat([a, b], dim=self.lon_name)
                else:  # 0_360 wrap
                    a = da2.sel({self.lon_name: slice(lo0, 360)})
                    b = da2.sel({self.lon_name: slice(0, lo1)})
                    b = b.assign_coords({self.lon_name: b[self.lon_name] + 360.0})
                    da2 = xr.concat([a, b], dim=self.lon_name)

        if da2.sizes.get(self.lat_name, 0) == 0 or da2.sizes.get(self.lon_name, 0) == 0:
            raise ValueError("Region bounds produced empty lat/lon selection.")

        return da2

    # ==========================================================
    # Stats helper
    # ==========================================================
    @staticmethod
    def _weighted_corr_rmse(
        obs2d: np.ndarray,
        mod2d: np.ndarray,
        w2d: Optional[np.ndarray] = None,
    ) -> Tuple[float, float]:
        """Return (corr, rmse) between two 2D arrays with optional area weights."""
        m = np.isfinite(obs2d) & np.isfinite(mod2d)
        if not np.any(m):
            return np.nan, np.nan

        x = obs2d[m].astype(float)
        y = mod2d[m].astype(float)

        if w2d is None:
            x0 = x - x.mean()
            y0 = y - y.mean()
            denom = np.sqrt((x0 * x0).mean()) * np.sqrt((y0 * y0).mean())
            corr = np.nan if denom == 0 else float((x0 * y0).mean() / denom)
            rmse = float(np.sqrt(((y - x) ** 2).mean()))
            return corr, rmse

        w = w2d[m].astype(float)
        w = w / np.sum(w)
        mx = np.sum(w * x)
        my = np.sum(w * y)
        x0 = x - mx
        y0 = y - my
        cov = np.sum(w * x0 * y0)
        sx = np.sqrt(np.sum(w * x0 * x0))
        sy = np.sqrt(np.sum(w * y0 * y0))
        corr = np.nan if (sx == 0 or sy == 0) else float(cov / (sx * sy))
        rmse = float(np.sqrt(np.sum(w * (y - x) ** 2)))
        return corr, rmse

    # ==========================================================
    # Internal helpers for multi-row plot API
    # ==========================================================
    def _panel_label(self, key: str) -> str:
        return self.plot_dict.get(key, {}).get("label", key)

    def _build_mean_maps(self, data_dict: Dict[str, xr.DataArray]) -> Dict[str, xr.DataArray]:
        if self.obs_key not in data_dict:
            raise ValueError(f"data_dict must contain obs_key='{self.obs_key}'")
        mean_maps: Dict[str, xr.DataArray] = {}
        ref = data_dict[self.obs_key]
        extra_dims = [d for d in ref.dims if d not in (self.lat_name, self.lon_name)]
        mean_maps[self.obs_key] = ref.mean(dim=extra_dims) if extra_dims else ref
        for key in self.group_order:
            if key not in data_dict:
                continue
            da = data_dict[key]
            extra_dims = [d for d in da.dims if d not in (self.lat_name, self.lon_name)]
            mean_maps[key] = da.mean(dim=extra_dims) if extra_dims else da
        return mean_maps

    def _build_spread_maps(self, data_dict: Dict[str, xr.DataArray]) -> Dict[str, xr.DataArray]:
        spread_maps: Dict[str, xr.DataArray] = {}
        for key in self.group_order:
            if key not in data_dict:
                continue
            da = data_dict[key]
            extra_dims = [d for d in da.dims if d not in (self.lat_name, self.lon_name)]
            if extra_dims:
                spread_maps[key] = da.std(dim=extra_dims)
        return spread_maps

    @staticmethod
    def _auto_symmetric_levels(
        arrs: List[np.ndarray],
        nlevels: int = 17,
        vlim: Optional[float] = None,
    ) -> Tuple[np.ndarray, float, float]:
        vals = np.concatenate([a.ravel() for a in arrs])
        finite = np.isfinite(vals)
        if not np.any(finite):
            vmax = 1.0 if vlim is None else float(vlim)
        else:
            vmax = float(np.nanmax(np.abs(vals[finite]))) if vlim is None else float(vlim)
        vmin = -vmax
        return np.linspace(vmin, vmax, nlevels), vmin, vmax

    @staticmethod
    def _auto_spread_level_from_quantile(
        spreads: List[np.ndarray],
        quantile: float = 0.75,
    ) -> Optional[float]:
        if not spreads:
            return None
        vals = np.concatenate([s.ravel() for s in spreads])
        finite = np.isfinite(vals)
        return float(np.nanpercentile(vals[finite], quantile * 100.0)) if np.any(finite) else None

    # ==========================================================
    # Multi-product grid figure
    # ==========================================================
    def plot_mode_season_maps(
        self,
        *,
        mode: str,
        season: str,
        products: Dict[str, Dict[str, xr.DataArray]],
        product_order: Optional[Sequence[str]] = None,
        product_labels: Optional[Dict[str, str]] = None,
        filename: Optional[str] = None,
        cmap: str = "RdBu_r",
        central_lon: float = 180.0,
        figsize: Tuple[float, float] = (15, 9),
        fontz: int = 13,
        yticks: Optional[np.ndarray] = None,
        xtick_step: float = 30.0,
        extent: Optional[Tuple[float, float, float, float]] = None,
        mlevels_by_product: Optional[Dict[str, Sequence[float]]] = None,
        nlevels: int = 17,
        overlay_spread: bool = True,
        spread_quantile: float = 0.75,
        spread_level_by_product: Optional[Dict[str, float]] = None,
        hatch: str = "....",
        significance_by_product: Optional[Dict[str, Dict[str, xr.DataArray]]] = None,
        sig_level: float = 0.05,
        sig_dot_color: str = "k",
        sig_dot_size: float = 2.0,
        sig_dot_density: int = 2,
        annotate_stats: bool = False,
        stats_area_weight: bool = True,
        stats_font_scale: float = 0.78,
        hide_inner_ylabels: bool = False,
        cb_labels_by_product: Optional[Dict[str, str]] = None,
        one_colorbar_per_row: bool = True,
        fig_format: str = "pdf",
        fig_dpi: int = 300,
        fig_idx_start: int = 0,
        left: float = 0.06,
        right: float = 0.98,
        bottom: float = 0.06,
        top: float = 0.93,
        hspace: float = 0.12,
        wspace: float = 0.12,
        row_cbar_pad: float = 0.028,
        row_cbar_height: float = 0.018,
    ):
        """
        Make a single figure for (mode, season) with multiple rows of map products.

        Parameters
        ----------
        products : dict
            product_key -> data_dict mapping (reference / hist / future DataArrays).
        one_colorbar_per_row : bool
            True  → one horizontal colorbar per product row.
            False → one shared vertical colorbar on the right.
        """
        product_order = list(product_order) if product_order is not None else list(products.keys())
        product_labels = product_labels or {}
        cb_labels_by_product = cb_labels_by_product or {}
        mlevels_by_product = mlevels_by_product or {}
        spread_level_by_product = spread_level_by_product or {}
        significance_by_product = significance_by_product or {}

        panel_keys = [self.obs_key] + list(self.group_order)
        ncols = len(panel_keys)
        nrows = len(product_order)

        mean_maps_by_prod: Dict[str, Dict[str, xr.DataArray]] = {}
        spread_maps_by_prod: Dict[str, Dict[str, xr.DataArray]] = {}
        for pk in product_order:
            dd = products[pk]
            mean_maps_by_prod[pk] = self._build_mean_maps(dd)
            spread_maps_by_prod[pk] = self._build_spread_maps(dd)

        pk0 = product_order[0]
        ref0 = mean_maps_by_prod[pk0][self.obs_key]
        lat = ref0[self.lat_name].values
        lon = ref0[self.lon_name].values

        if yticks is None:
            yticks = np.arange(-90, 91, 30)
        if extent is None:
            extent = (float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max()))
        xticks = np.arange(
            np.ceil(extent[0] / xtick_step) * xtick_step,
            extent[1] + 1e-9,
            xtick_step,
        )

        fig = plt.figure(figsize=figsize)
        proj = ccrs.PlateCarree(central_longitude=central_lon)
        data_crs = ccrs.PlateCarree()
        axes: List[Any] = []
        ims: Dict[str, Any] = {}

        fig.suptitle(f"{mode} — {season}", fontsize=fontz * 1.25, y=0.98)

        if not one_colorbar_per_row:
            right = 0.92

        for r, pk in enumerate(product_order):
            mp = mean_maps_by_prod[pk]
            sp = spread_maps_by_prod[pk]

            if pk in mlevels_by_product and mlevels_by_product[pk] is not None:
                mlevels = np.asarray(list(mlevels_by_product[pk]), dtype=float)
                vmin, vmax = float(np.min(mlevels)), float(np.max(mlevels))
            else:
                arrs = [mp[k].values for k in panel_keys]
                mlevels, vmin, vmax = self._auto_symmetric_levels(arrs, nlevels=nlevels)

            spread_level = spread_level_by_product.get(pk)
            if overlay_spread and spread_level is None and sp:
                spread_arrays = [sp[k].values for k in self.group_order if k in sp]
                spread_level = self._auto_spread_level_from_quantile(spread_arrays, quantile=spread_quantile)

            for c, key in enumerate(panel_keys):
                idx = r * ncols + c + 1
                ax = fig.add_subplot(nrows, ncols, idx, projection=proj)
                ax.set_aspect("auto")
                axes.append(ax)

                da_map = mp[key]
                im = ax.contourf(
                    lon, lat, da_map, levels=mlevels, cmap=cmap,
                    vmin=vmin, vmax=vmax, transform=data_crs, extend="both",
                )
                ims[pk] = im
                ax.contour(
                    lon, lat, da_map, levels=mlevels,
                    colors="k", linewidths=0.3, transform=data_crs,
                )
                ax.coastlines(linewidth=0.5)
                ax.set_extent(extent, crs=data_crs)
                ax.set_yticks(yticks, crs=data_crs)
                ax.set_xticks(xticks, crs=data_crs)
                ax.xaxis.set_major_formatter(LongitudeFormatter(".0f"))
                ax.yaxis.set_major_formatter(LatitudeFormatter(".0f"))

                if overlay_spread and (key in sp) and (spread_level is not None):
                    spread_mask = np.where(sp[key].values > spread_level, 1.0, np.nan)
                    ax.contourf(
                        lon, lat, spread_mask, levels=[0.5, 1.5],
                        hatches=[hatch], colors="none", transform=data_crs,
                    )

                pval = significance_by_product.get(pk, {}).get(key)
                if pval is not None:
                    pval = pval.transpose(self.lat_name, self.lon_name)
                    if not (
                        pval[self.lat_name].size == da_map[self.lat_name].size
                        and pval[self.lon_name].size == da_map[self.lon_name].size
                        and np.allclose(pval[self.lat_name], da_map[self.lat_name])
                        and np.allclose(pval[self.lon_name], da_map[self.lon_name])
                    ):
                        pval = pval.interp({
                            self.lat_name: da_map[self.lat_name],
                            self.lon_name: da_map[self.lon_name],
                        })
                    add_sig_dots(
                        ax, pval.rename({self.lat_name: "lat", self.lon_name: "lon"}),
                        sig_level=sig_level, dot_color=sig_dot_color,
                        dot_size=sig_dot_size, dot_density=sig_dot_density,
                        transform=data_crs,
                    )

                panel_letter = chr(97 + (idx - 1 + fig_idx_start))
                ax.set_title(f"({panel_letter}) {self._panel_label(key)}", loc="left", fontsize=fontz)
                product_title = product_labels.get(pk, pk)
                if c == ncols - 1 and product_title:
                    ax.set_title(product_title, loc="right", fontsize=fontz)

                if annotate_stats and key != self.obs_key:
                    lon_convention = "negpos" if extent[0] < 0 else "0_360"
                    ref_stats = self._normalize_lon_for_bounds(mp[self.obs_key], lon_convention)
                    map_stats = self._normalize_lon_for_bounds(da_map, lon_convention)
                    lat_lo, lat_hi = sorted((extent[2], extent[3]))
                    lon_lo, lon_hi = sorted((extent[0], extent[1]))
                    ref_stats = ref_stats.sel({
                        self.lat_name: slice(lat_lo, lat_hi),
                        self.lon_name: slice(lon_lo, lon_hi),
                    }).transpose(self.lat_name, self.lon_name)
                    map_stats = map_stats.interp({
                        self.lat_name: ref_stats[self.lat_name],
                        self.lon_name: ref_stats[self.lon_name],
                    }).transpose(self.lat_name, self.lon_name)
                    if stats_area_weight:
                        weights = np.cos(np.deg2rad(ref_stats[self.lat_name].values))
                        weights_2d = weights[:, None] * np.ones(ref_stats.shape, dtype=float)
                    else:
                        weights_2d = None
                    pcc, rmsd = self._weighted_corr_rmse(
                        ref_stats.values, map_stats.values, w2d=weights_2d,
                    )
                    ax.text(
                        0.03, 0.04, f"PCC = {pcc:.2f}\nRMSD = {rmsd:.2f}",
                        transform=ax.transAxes, ha="left", va="bottom",
                        fontsize=fontz * stats_font_scale,
                        bbox=dict(
                            facecolor="white", edgecolor="0.7", alpha=0.65,
                            boxstyle="round,pad=0.22",
                        ),
                        zorder=7,
                    )

                ax.tick_params(labelsize=fontz * 0.9)
                if c == 0:
                    ax.set_ylabel("Latitude", fontsize=fontz)
                elif hide_inner_ylabels:
                    ax.tick_params(axis="y", labelleft=False)
                ax.set_xlabel("Longitude", fontsize=fontz)

        fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top, hspace=hspace, wspace=wspace)

        if one_colorbar_per_row:
            for r, pk in enumerate(product_order):
                row_height = (top - bottom - hspace * (nrows - 1)) / nrows
                row_bottom = top - (r + 1) * row_height - r * hspace
                cax = fig.add_axes([left, row_bottom - row_cbar_pad, right - left, row_cbar_height])
                cbar = fig.colorbar(ims[pk], cax=cax, orientation="horizontal")
                cbar.ax.tick_params(labelsize=fontz * 0.85)
                cbar.set_label(cb_labels_by_product.get(pk, "Map value"), fontsize=fontz * 0.9)
        else:
            cax = fig.add_axes([0.94, 0.15, 0.015, 0.70])
            cbar = fig.colorbar(ims[product_order[-1]], cax=cax, orientation="vertical")
            cbar.ax.tick_params(labelsize=fontz * 0.85)
            cbar.set_label("Map value", fontsize=fontz * 0.9)

        if filename is None:
            filename = f"{mode}_{season}_mode_maps.{fig_format}"
        out_path = os.path.join(self.fig_dir, filename)
        fig.savefig(out_path, dpi=fig_dpi, format=fig_format, bbox_inches="tight", pad_inches=0.05)
        print(f"Saved: {out_path}")
        return fig, axes

    # ==========================================================
    # Core engine: multimodel panel plot (regional OR global)
    # ==========================================================
    def plot_multimodel_panel_with_stats(
        self,
        *,
        title: str,
        obs_map: xr.DataArray,
        model_stack: xr.DataArray,     # (member, lat, lon)
        member_labels: Sequence[str],
        obs_pval_map: Optional[xr.DataArray] = None,
        model_pval_stack: Optional[xr.DataArray] = None,
        member_dim: str = "member",
        filename: str = "mode_multimodel.pdf",
        cmap: str = "RdBu_r",
        central_lon: float = 180.0,
        fontz: int = 12,
        fig_dpi: int = 300,
        fig_format: str = "pdf",
        area_weight: bool = True,
        region_lat_bounds: Optional[Tuple[float, float]] = None,
        region_lon_bounds: Optional[Tuple[float, float]] = None,
        lon_convention: str = "negpos",
        ncols: int = 2,
        figsize_per_panel: Tuple[float, float] = (5.0, 3.3),
        wspace: float = 0.10,
        hspace: float = 0.18,
        cb_pad: float = 0.08,
        xtick_step: float = 20.0,
        yticks: Optional[np.ndarray] = None,
        mlevels: Optional[Sequence[float]] = None,
        nlevels: int = 17,
        cbar_ticks: Optional[Sequence[float]] = None,
        cbar_tick_labels: Optional[Sequence[str]] = None,
        cbar_label: str = "EOF pattern (units as provided)",
        extent_override: Optional[Tuple[float, float, float, float]] = None,
        show_significance: bool = False,
        sig_level: float = 0.05,
        dot_color: str = "k",
        dot_size: float = 0.8,
        dot_density: int = 2,
    ):
        """
        Plot Obs + each model in a (nrows x ncols) grid with per-panel r/RMSE annotations.

        region_* provided → crop maps and compute stats over that region.
        region_* None     → global plotting + global stats.
        """
        panels: List[Tuple[str, xr.DataArray]] = [(self._panel_label(self.obs_key), obs_map)]
        for i in range(model_stack.sizes[member_dim]):
            panels.append((member_labels[i], model_stack.isel({member_dim: i})))

        n_panels = len(panels)
        ncols = max(1, int(ncols))
        nrows = int(math.ceil(n_panels / ncols))

        panels_plot: List[Tuple[str, xr.DataArray]] = []
        pvals_plot: List[Optional[xr.DataArray]] = []
        pval_panels: List[Optional[xr.DataArray]] = [obs_pval_map]
        for i in range(model_stack.sizes[member_dim]):
            if model_pval_stack is not None and i < model_pval_stack.sizes.get(member_dim, 0):
                pval_panels.append(model_pval_stack.isel({member_dim: i}))
            else:
                pval_panels.append(None)

        # Align all models globally to obs grid first to avoid edge-interpolation artifacts.
        obs_map_aligned = panels[0][1].transpose(self.lat_name, self.lon_name)
        panels[0] = (panels[0][0], obs_map_aligned)
        for idx in range(1, len(panels)):
            name, da = panels[idx]
            da = da.transpose(self.lat_name, self.lon_name)
            if not (
                da[self.lon_name].size == obs_map_aligned[self.lon_name].size
                and da[self.lat_name].size == obs_map_aligned[self.lat_name].size
                and np.allclose(da[self.lon_name].values, obs_map_aligned[self.lon_name].values, atol=1e-6)
                and np.allclose(da[self.lat_name].values, obs_map_aligned[self.lat_name].values, atol=1e-6)
            ):
                da = da.interp({self.lon_name: obs_map_aligned[self.lon_name], self.lat_name: obs_map_aligned[self.lat_name]})
            else:
                da = da.assign_coords({self.lon_name: obs_map_aligned[self.lon_name], self.lat_name: obs_map_aligned[self.lat_name]})
            panels[idx] = (name, da)

        for idx, pval in enumerate(pval_panels):
            if pval is None:
                continue
            pval = pval.transpose(self.lat_name, self.lon_name)
            if not (
                pval[self.lon_name].size == obs_map_aligned[self.lon_name].size
                and pval[self.lat_name].size == obs_map_aligned[self.lat_name].size
                and np.allclose(pval[self.lon_name].values, obs_map_aligned[self.lon_name].values, atol=1e-6)
                and np.allclose(pval[self.lat_name].values, obs_map_aligned[self.lat_name].values, atol=1e-6)
            ):
                pval = pval.interp({self.lon_name: obs_map_aligned[self.lon_name], self.lat_name: obs_map_aligned[self.lat_name]})
            else:
                pval = pval.assign_coords({self.lon_name: obs_map_aligned[self.lon_name], self.lat_name: obs_map_aligned[self.lat_name]})
            pval_panels[idx] = pval

        if region_lat_bounds is None and region_lon_bounds is None:
            obs0 = self._normalize_lon_for_bounds(obs_map, "0_360")
            panels_plot.append((panels[0][0], obs0))
            pvals_plot.append(
                self._normalize_lon_for_bounds(obs_pval_map, "0_360")
                if obs_pval_map is not None else None
            )
            for i in range(model_stack.sizes[member_dim]):
                da = self._normalize_lon_for_bounds(model_stack.isel({member_dim: i}), "0_360")
                panels_plot.append((member_labels[i], da))
                pval = pval_panels[i + 1]
                pvals_plot.append(
                    self._normalize_lon_for_bounds(pval, "0_360")
                    if pval is not None else None
                )
        else:
            for idx, (name, da) in enumerate(panels):
                panels_plot.append((name, self._subset_latlon(da, region_lat_bounds, region_lon_bounds, lon_convention)))
                pval = pval_panels[idx]
                pvals_plot.append(
                    self._subset_latlon(pval, region_lat_bounds, region_lon_bounds, lon_convention)
                    if pval is not None else None
                )

        obs_plot = panels_plot[0][1]
        obs_plot = obs_plot.transpose(self.lat_name, self.lon_name)
        panels_plot[0] = (panels_plot[0][0], obs_plot)
        for idx in range(1, len(panels_plot)):
            name, da = panels_plot[idx]
            da = da.transpose(self.lat_name, self.lon_name)
            if not (
                da[self.lon_name].size == obs_plot[self.lon_name].size
                and da[self.lat_name].size == obs_plot[self.lat_name].size
                and np.allclose(da[self.lon_name].values, obs_plot[self.lon_name].values, atol=1e-6)
                and np.allclose(da[self.lat_name].values, obs_plot[self.lat_name].values, atol=1e-6)
            ):
                da = da.interp({self.lon_name: obs_plot[self.lon_name], self.lat_name: obs_plot[self.lat_name]})
            else:
                # Snap coordinates exactly to obs grid to avoid downstream mismatches.
                da = da.assign_coords({self.lon_name: obs_plot[self.lon_name], self.lat_name: obs_plot[self.lat_name]})
            panels_plot[idx] = (name, da.transpose(self.lat_name, self.lon_name))

        for idx, pval in enumerate(pvals_plot):
            if pval is None:
                continue
            pval = pval.transpose(self.lat_name, self.lon_name)
            if not (
                pval[self.lon_name].size == obs_plot[self.lon_name].size
                and pval[self.lat_name].size == obs_plot[self.lat_name].size
                and np.allclose(pval[self.lon_name].values, obs_plot[self.lon_name].values, atol=1e-6)
                and np.allclose(pval[self.lat_name].values, obs_plot[self.lat_name].values, atol=1e-6)
            ):
                pval = pval.interp({self.lon_name: obs_plot[self.lon_name], self.lat_name: obs_plot[self.lat_name]})
            else:
                pval = pval.assign_coords({self.lon_name: obs_plot[self.lon_name], self.lat_name: obs_plot[self.lat_name]})
            pvals_plot[idx] = pval.transpose(self.lat_name, self.lon_name)

        lat_plot = obs_plot[self.lat_name].values
        lon_plot = obs_plot[self.lon_name].values

        if yticks is None:
            yticks = np.arange(np.floor(lat_plot.min() / 10) * 10, np.ceil(lat_plot.max() / 10) * 10 + 1e-9, 10)
        xticks = np.arange(
            np.floor(lon_plot.min() / xtick_step) * xtick_step,
            np.ceil(lon_plot.max() / xtick_step) * xtick_step + 1e-9,
            xtick_step,
        )

        if mlevels is None:
            all_vals = np.concatenate([da.values.ravel() for _, da in panels_plot])
            finite = np.isfinite(all_vals)
            vmax = float(np.nanmax(np.abs(all_vals[finite]))) if np.any(finite) else 1.0
            mlevels_arr = np.linspace(-vmax, vmax, nlevels)
            vmin = -vmax
        else:
            mlevels_arr = np.asarray(list(mlevels), dtype=float)
            vmin, vmax = float(np.min(mlevels_arr)), float(np.max(mlevels_arr))

        if area_weight:
            w = np.cos(np.deg2rad(lat_plot)).astype(float)
            w = w / np.nanmean(w)
            w2d = w[:, None] * np.ones((lat_plot.size, lon_plot.size), dtype=float)
        else:
            w2d = None

        fig_w = figsize_per_panel[0] * ncols
        fig_h = figsize_per_panel[1] * nrows + 0.9
        fig = plt.figure(figsize=(fig_w, fig_h))
        proj = ccrs.PlateCarree(central_longitude=central_lon)
        data_crs = ccrs.PlateCarree()
        fig.suptitle(title, fontsize=fontz * 1.25, y=0.98)

        axes: List[Any] = []
        im_last = None
        extent = extent_override or (
            float(lon_plot.min()), float(lon_plot.max()),
            float(lat_plot.min()), float(lat_plot.max()),
        )

        for p in range(n_panels):
            ax = fig.add_subplot(nrows, ncols, p + 1, projection=proj)
            axes.append(ax)
            name, da = panels_plot[p]

            im_last = ax.contourf(
                lon_plot, lat_plot, da.values,
                levels=mlevels_arr, cmap=cmap, vmin=vmin, vmax=vmax,
                transform=data_crs, extend="both",
            )
            ax.contour(
                lon_plot, lat_plot, da.values,
                levels=mlevels_arr, colors="k", linewidths=0.3, transform=data_crs,
            )
            if show_significance and pvals_plot[p] is not None:
                pval = pvals_plot[p]
                if not (
                    np.array_equal(pval[self.lon_name].values, lon_plot)
                    and np.array_equal(pval[self.lat_name].values, lat_plot)
                ):
                    pval = pval.interp({self.lon_name: lon_plot, self.lat_name: lat_plot})
                sig_mask = np.isfinite(pval.values) & (pval.values < sig_level)
                lat_idx = np.arange(0, len(lat_plot), max(1, int(dot_density)))
                lon_idx = np.arange(0, len(lon_plot), max(1, int(dot_density)))
                mask_sub = sig_mask[np.ix_(lat_idx, lon_idx)]
                lon2d, lat2d = np.meshgrid(lon_plot[lon_idx], lat_plot[lat_idx])
                ax.scatter(
                    lon2d[mask_sub],
                    lat2d[mask_sub],
                    s=dot_size,
                    c=dot_color,
                    marker=".",
                    linewidths=0,
                    transform=data_crs,
                    zorder=5,
                )
            ax.coastlines(linewidth=0.6)
            ax.set_extent(extent, crs=data_crs)
            ax.set_xticks(xticks, crs=data_crs)
            ax.set_yticks(yticks, crs=data_crs)
            ax.xaxis.set_major_formatter(LongitudeFormatter(".0f"))
            ax.yaxis.set_major_formatter(LatitudeFormatter(".0f"))
            ax.tick_params(labelsize=fontz * 0.85)
            ax.set_title(f"({chr(97 + p)}) {name}", loc="left", fontsize=fontz)

            if p > 0:
                r0, rmsd0 = self._weighted_corr_rmse(panels_plot[0][1].values, da.values, w2d=w2d)
                ax.text(
                    0.05, 0.05, f"PCC = {r0:.2f}\nRMSD = {rmsd0:.2f}",
                    transform=ax.transAxes, ha="left", va="bottom", fontsize=fontz * 0.9,
                    bbox=dict(facecolor="white", edgecolor="0.7", alpha=0.5, boxstyle="round,pad=0.25"),
                )

            if (p % ncols) == 0:
                ax.set_ylabel("Latitude", fontsize=fontz)
            ax.set_xlabel("Longitude", fontsize=fontz)

        fig.subplots_adjust(left=0.06, right=0.98, top=0.92, bottom=0.12, wspace=wspace, hspace=hspace)

        cbar = fig.colorbar(
            im_last,
            ax=axes,
            orientation="horizontal",
            fraction=0.05,
            pad=cb_pad,
            aspect=45,
            ticks=cbar_ticks,
        )
        if cbar_tick_labels is not None:
            if cbar_ticks is None:
                raise ValueError("cbar_tick_labels requires cbar_ticks.")
            if len(cbar_tick_labels) != len(cbar_ticks):
                raise ValueError("cbar_tick_labels must have the same length as cbar_ticks.")
            cbar.ax.set_xticklabels(cbar_tick_labels)
        cbar.ax.tick_params(labelsize=fontz * 0.85)
        cbar.set_label(cbar_label, fontsize=fontz)

        out_path = os.path.join(self.fig_dir, filename)
        fig.savefig(out_path, dpi=fig_dpi, format=fig_format, bbox_inches="tight", pad_inches=0.05)
        print(f"Saved: {out_path}")
        return fig, axes

    # ==========================================================
    # Convenience wrappers
    # ==========================================================
    def plot_multimodel_mode_pattern_with_stats(
        self,
        *,
        mode: str,
        token: str,
        obs_map: xr.DataArray,
        model_stack: xr.DataArray,
        member_labels: Sequence[str],
        filename: str,
        region_lat_bounds: Tuple[float, float],
        region_lon_bounds: Tuple[float, float],
        lon_convention: str = "negpos",
        central_lon: float = 0.0,
        cbar_label: str = "EOF pattern (units as provided)",
        title: Optional[str] = None,
        **kwargs,
    ):
        return self.plot_multimodel_panel_with_stats(
            title=title or f"{mode} pattern — {token}",
            obs_map=obs_map,
            model_stack=model_stack,
            member_labels=member_labels,
            filename=filename,
            region_lat_bounds=region_lat_bounds,
            region_lon_bounds=region_lon_bounds,
            lon_convention=lon_convention,
            central_lon=central_lon,
            cbar_label=cbar_label,
            **kwargs,
        )

    def plot_multimodel_teleconnection_with_stats(
        self,
        *,
        mode: str,
        token: str,
        obs_map: xr.DataArray,
        model_stack: xr.DataArray,
        member_labels: Sequence[str],
        filename: str,
        central_lon: float = 180.0,
        cbar_label: str = "Teleconnection slope (units as provided)",
        title: Optional[str] = None,
        **kwargs,
    ):
        return self.plot_multimodel_panel_with_stats(
            title=title or f"{mode} teleconnection — {token}",
            obs_map=obs_map,
            model_stack=model_stack,
            member_labels=member_labels,
            filename=filename,
            region_lat_bounds=None,
            region_lon_bounds=None,
            lon_convention="0_360",
            central_lon=central_lon,
            cbar_label=cbar_label,
            **kwargs,
        )

class MultimodelPCTimeSeriesPlotter:
    def __init__(
        self,
        fig_dir: str,
        plot_dict: Dict[str, dict],
        group_order: Sequence[str] = ("hist",),
        obs_key: str = "reference",
    ):
        self.fig_dir = fig_dir
        self.plot_dict = plot_dict
        self.group_order = tuple(group_order)
        self.obs_key = obs_key
        os.makedirs(self.fig_dir, exist_ok=True)

    def _panel_label(self, key: str) -> str:
        return self.plot_dict.get(key, {}).get("label", key)

    def plot_multimodel_pc_timeseries_with_stats(
        self,
        mode,
        token,
        obs_pc,
        model_stack,
        member_labels,
        filename,
        **kwargs,
    ):
        from scipy import stats

        cfg = dict(
            pc_var="pc_proj",
            member_dim="member",

            obs_label=None,
            model_styles={},

            figsize=(14, 6),
            fontz=18,
            fig_dpi=400,
            fig_format="pdf",

            obs_bar_width=0.80,
            obs_bar_alpha=0.18,
            obs_bar_color="black",
            obs_bar_edgecolor="black",
            obs_bar_linewidth=0.8,

            linewidth_model=2.3,
            markersize=6.0,
            markeredgewidth=1.0,
            alpha_model=0.95,

            use_time_coord=True,
            time_coord_offset=0,
            xtick_step=2,
            xtick_start=None,
            xtick_end=None,
            xlim=None,
            sort_by_corr=True,

            zero_line_color="0.25",
            zero_linewidth=1.2,
            zero_linestyle="--",
            zero_alpha=0.85,

            grid=True,
            grid_axis="y",
            grid_linestyle="--",
            grid_linewidth=0.6,
            grid_alpha=0.35,

            legend_below=True,
            legend_frameon=True,
            legend_edgecolor="black",
            legend_facecolor="white",
            legend_framealpha=1.0,
            legend_ncol=3,

            legend_stat="corr",
            show_corr=True,
            show_corr_significance=True,
            show_corr_pvalue=False,

            ylim=None,
            xlabel="Time",
            ylabel="Standardized PC",
            title=None,
            corr_prefix="r",
            corr_fmt="{:.2f}",
            corr_label_style="paren",
        )
        cfg.update(kwargs)

        if cfg["obs_label"] is None:
            cfg["obs_label"] = self._panel_label(self.obs_key)

        obs_raw = self._as_1d(obs_pc)
        obs_vals = self._standardize(obs_raw)
        obs_x = self._time_coord(obs_pc, fallback_size=obs_vals.size) if cfg["use_time_coord"] else np.arange(obs_vals.size)
        obs_x = self._offset_time_coord(obs_x, cfg["time_coord_offset"])

        fig, ax = plt.subplots(figsize=cfg["figsize"])

        ax.bar(
            obs_x,
            obs_vals,
            width=cfg["obs_bar_width"],
            color=cfg["obs_bar_color"],
            edgecolor=cfg["obs_bar_edgecolor"],
            linewidth=cfg["obs_bar_linewidth"],
            alpha=cfg["obs_bar_alpha"],
            label=cfg["obs_label"],
            zorder=1,
        )

        model_records = []
        for i, label in enumerate(member_labels):
            pc = self._select_member(
                model_stack,
                i,
                label,
                member_dim=cfg["member_dim"],
            )
            model_raw = self._as_1d(pc)
            model_vals = self._standardize(model_raw)
            model_x = self._time_coord(pc, fallback_size=model_vals.size) if cfg["use_time_coord"] else np.arange(model_vals.size)
            model_x = self._offset_time_coord(model_x, cfg["time_coord_offset"])
            r, p = self._corr_p_aligned(obs_vals, model_vals, obs_x, model_x)
            std_ratio = self._std_ratio_aligned(obs_raw, model_raw, obs_x, model_x)
            var_ratio = std_ratio**2 if np.isfinite(std_ratio) else np.nan
            model_records.append((label, model_x, model_vals, r, p, std_ratio, var_ratio))

        if cfg["sort_by_corr"]:
            model_records = sorted(
                model_records,
                key=lambda item: -999.0 if not np.isfinite(item[3]) else item[3],
            )

        for label, model_x, model_vals, r, p, std_ratio, var_ratio in model_records:
            style = cfg["model_styles"].get(label, {})
            plot_label = self._build_label(
                label,
                r,
                p,
                std_ratio=std_ratio,
                var_ratio=var_ratio,
                legend_stat=cfg["legend_stat"],
                show_corr=cfg["show_corr"],
                show_corr_significance=cfg["show_corr_significance"],
                show_corr_pvalue=cfg["show_corr_pvalue"],
                corr_prefix=cfg["corr_prefix"],
                corr_fmt=cfg["corr_fmt"],
                corr_label_style=cfg["corr_label_style"],
            )

            ax.plot(
                model_x,
                model_vals,
                color=style.get("color", None),
                linestyle=style.get("linestyle", "-"),
                marker=style.get("marker", "o"),
                linewidth=cfg["linewidth_model"],
                markersize=cfg["markersize"],
                markeredgewidth=cfg["markeredgewidth"],
                alpha=cfg["alpha_model"],
                label=plot_label,
                zorder=3,
            )

        ax.axhline(
            0.0,
            color=cfg["zero_line_color"],
            linewidth=cfg["zero_linewidth"],
            linestyle=cfg["zero_linestyle"],
            alpha=cfg["zero_alpha"],
        )
        ax.set_xlabel(cfg["xlabel"], fontsize=cfg["fontz"])
        ax.set_ylabel(cfg["ylabel"], fontsize=cfg["fontz"])
        ax.tick_params(labelsize=cfg["fontz"] - 2)

        if cfg["ylim"] is not None:
            ax.set_ylim(cfg["ylim"])

        if cfg["grid"]:
            ax.grid(
                True,
                axis=cfg["grid_axis"],
                linestyle=cfg["grid_linestyle"],
                linewidth=cfg["grid_linewidth"],
                alpha=cfg["grid_alpha"],
            )

        if cfg["xlim"] is not None:
            ax.set_xlim(cfg["xlim"])

        self._set_time_ticks(
            ax,
            obs_x,
            cfg["xtick_step"],
            start=cfg["xtick_start"],
            end=cfg["xtick_end"],
            xlim=cfg["xlim"],
        )

        title = cfg["title"] or f"{mode} {token} standardized {cfg['pc_var']} time series"
        ax.set_title(title, fontsize=cfg["fontz"])

        self._add_legend(ax, cfg)

        fig.tight_layout()
        if cfg["legend_below"]:
            fig.subplots_adjust(bottom=0.25)

        outpath = os.path.join(self.fig_dir, filename)
        fig.savefig(
            outpath,
            dpi=cfg["fig_dpi"],
            format=cfg["fig_format"],
            bbox_inches="tight",
        )
        plt.close(fig)

        print(f"Saved: {outpath}")
        return fig, ax

    @staticmethod
    def _as_1d(da):
        vals = np.asarray(da).squeeze()
        return vals.astype(float)

    @staticmethod
    def _standardize(x):
        x = np.asarray(x, dtype=float)
        out = np.full(x.shape, np.nan)
        mask = np.isfinite(x)

        if mask.sum() < 2:
            return out

        std = x[mask].std()
        if not np.isfinite(std) or std == 0:
            return out

        out[mask] = (x[mask] - x[mask].mean()) / std
        return out

    @staticmethod
    def _corr_p(x, y):
        from scipy import stats

        mask = np.isfinite(x) & np.isfinite(y)
        if mask.sum() < 3:
            return np.nan, np.nan

        r, p = stats.pearsonr(x[mask], y[mask])
        return r, p

    @classmethod
    def _corr_p_aligned(cls, ref_vals, test_vals, ref_time, test_time):
        ref_vals = np.asarray(ref_vals, dtype=float).ravel()
        test_vals = np.asarray(test_vals, dtype=float).ravel()
        ref_time = np.asarray(ref_time).ravel()
        test_time = np.asarray(test_time).ravel()

        if ref_time.size == ref_vals.size and test_time.size == test_vals.size:
            common_time = np.intersect1d(ref_time, test_time)
            if common_time.size >= 3:
                ref_by_time = dict(zip(ref_time, ref_vals))
                test_by_time = dict(zip(test_time, test_vals))
                ref_aligned = np.asarray([ref_by_time[t] for t in common_time], dtype=float)
                test_aligned = np.asarray([test_by_time[t] for t in common_time], dtype=float)
                return cls._corr_p(ref_aligned, test_aligned)

        n = min(ref_vals.size, test_vals.size)
        return cls._corr_p(ref_vals[:n], test_vals[:n])

    @classmethod
    def _std_ratio_aligned(cls, ref_vals, test_vals, ref_time, test_time):
        ref_vals, test_vals = cls._align_by_time_or_length(ref_vals, test_vals, ref_time, test_time)
        mask = np.isfinite(ref_vals) & np.isfinite(test_vals)
        if mask.sum() < 2:
            return np.nan

        ref_std = np.nanstd(ref_vals[mask])
        test_std = np.nanstd(test_vals[mask])
        if not np.isfinite(ref_std) or ref_std == 0 or not np.isfinite(test_std):
            return np.nan

        return test_std / ref_std

    @staticmethod
    def _align_by_time_or_length(ref_vals, test_vals, ref_time, test_time):
        ref_vals = np.asarray(ref_vals, dtype=float).ravel()
        test_vals = np.asarray(test_vals, dtype=float).ravel()
        ref_time = np.asarray(ref_time).ravel()
        test_time = np.asarray(test_time).ravel()

        if ref_time.size == ref_vals.size and test_time.size == test_vals.size:
            common_time = np.intersect1d(ref_time, test_time)
            if common_time.size >= 2:
                ref_by_time = dict(zip(ref_time, ref_vals))
                test_by_time = dict(zip(test_time, test_vals))
                return (
                    np.asarray([ref_by_time[t] for t in common_time], dtype=float),
                    np.asarray([test_by_time[t] for t in common_time], dtype=float),
                )

        n = min(ref_vals.size, test_vals.size)
        return ref_vals[:n], test_vals[:n]

    @staticmethod
    def _sig_marker(p):
        if not np.isfinite(p):
            return ""
        if p < 0.05:
            return "**"
        if p < 0.10:
            return "*"
        return ""

    def _build_label(
        self,
        label,
        r,
        p,
        *,
        std_ratio=np.nan,
        var_ratio=np.nan,
        legend_stat="corr",
        show_corr=True,
        show_corr_significance=True,
        show_corr_pvalue=False,
        corr_prefix="r",
        corr_fmt="{:.2f}",
        corr_label_style="paren",
    ):
        plot_label = str(label)
        
        if legend_stat == "std_ratio":
            if not np.isfinite(std_ratio):
                return label
            return rf"{label} ($\sigma/\sigma_{{\mathrm{{obs}}}}={std_ratio:.2f}$)"

        if legend_stat == "var_ratio":
            if not np.isfinite(var_ratio):
                return plot_label
            return f"{label} (Var Ratio={var_ratio:.2f})"

        if legend_stat in ("none", None):
            return plot_label

        if not show_corr:
            return plot_label

        if not np.isfinite(r):
            return plot_label

        corr_text = corr_fmt.format(r)
        if corr_label_style == "raw":
            corr_label = f"{corr_prefix}={corr_text}"
            if show_corr_significance:
                corr_label += self._sig_marker(p)
            if show_corr_pvalue:
                corr_label += f" (p={p:.3f})"
            return f"{label} ({corr_label})"

        if show_corr_pvalue:
            return f"{label} ({corr_prefix}={corr_text}, p={p:.3f})"

        if show_corr_significance:
            return f"{label} ({corr_prefix}={corr_text}{self._sig_marker(p)})"

        return f"{label} ({corr_prefix}={corr_text})"

    @staticmethod
    def _select_member(model_stack, i, label, member_dim="member"):
        if hasattr(model_stack, "dims") and member_dim in model_stack.dims:
            return model_stack.isel({member_dim: i})

        if hasattr(model_stack, "sel"):
            return model_stack.sel({member_dim: label})

        return model_stack[i]

    @staticmethod
    def _time_coord(da, fallback_size):
        if hasattr(da, "coords") and "time" in da.coords:
            time = np.asarray(da.time.values).ravel()
            if time.size == fallback_size:
                if np.issubdtype(time.dtype, np.datetime64):
                    return time.astype("datetime64[Y]").astype(int) + 1970
                if time.size > 0 and hasattr(time[0], "year"):
                    return np.asarray([t.year for t in time])
                return time
        return np.arange(fallback_size)

    @staticmethod
    def _offset_time_coord(x, offset):
        if offset in (None, 0):
            return x
        vals = np.asarray(x)
        if np.issubdtype(vals.dtype, np.number):
            return vals + offset
        return x

    @staticmethod
    def _set_time_ticks(ax, x, step, *, start=None, end=None, xlim=None):
        if step is None:
            return

        try:
            vals = np.asarray(x)
            if vals.size == 0 or not np.issubdtype(vals.dtype, np.number):
                return

            step = int(step)
            if step <= 0:
                return

            if start is not None or end is not None or xlim is not None:
                if start is None:
                    start = xlim[0] if xlim is not None else np.nanmin(vals)
                if end is None:
                    end = xlim[1] if xlim is not None else np.nanmax(vals)
                ticks = np.arange(np.ceil(start), np.floor(end) + 1e-9, step)
            else:
                ticks = vals[::step]

            ax.set_xticks(ticks)
            ax.set_xticklabels([str(int(t)) if float(t).is_integer() else str(t) for t in ticks])
        except Exception:
            return

    @staticmethod
    def _add_legend(ax, cfg):
        if cfg["legend_below"]:
            leg = ax.legend(
                loc="upper center",
                bbox_to_anchor=(0.5, -0.18),
                ncol=cfg["legend_ncol"],
                fontsize=cfg["fontz"] - 4,
                frameon=cfg["legend_frameon"],
            )
        else:
            leg = ax.legend(
                loc="best",
                ncol=cfg["legend_ncol"],
                fontsize=cfg["fontz"] - 4,
                frameon=cfg["legend_frameon"],
            )

        if leg is not None and cfg["legend_frameon"]:
            frame = leg.get_frame()
            frame.set_edgecolor(cfg["legend_edgecolor"])
            frame.set_facecolor(cfg["legend_facecolor"])
            frame.set_alpha(cfg["legend_framealpha"])
