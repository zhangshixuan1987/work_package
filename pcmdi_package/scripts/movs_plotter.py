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
        cb_labels_by_product: Optional[Dict[str, str]] = None,
        one_colorbar_per_row: bool = True,
        fig_format: str = "pdf",
        fig_dpi: int = 300,
        fig_idx_start: int = 0,
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
        xticks = np.arange(np.floor(lon.min()), np.ceil(lon.max()) + 1e-9, xtick_step)

        if extent is None:
            extent = (float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max()))

        fig = plt.figure(figsize=figsize)
        proj = ccrs.PlateCarree(central_longitude=central_lon)
        data_crs = ccrs.PlateCarree()
        axes: List[Any] = []
        ims: Dict[str, Any] = {}

        fig.suptitle(f"{mode} — {season}", fontsize=fontz * 1.25, y=0.98)

        left, right, bottom, top = 0.06, 0.98, 0.06, 0.93
        hspace, wspace = 0.12, 0.12
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

                panel_letter = chr(97 + (idx - 1 + fig_idx_start))
                ax.set_title(f"({panel_letter}) {self._panel_label(key)}", loc="left", fontsize=fontz)
                if c == ncols - 1:
                    ax.set_title(product_labels.get(pk, pk), loc="right", fontsize=fontz)

                ax.tick_params(labelsize=fontz * 0.9)
                if c == 0:
                    ax.set_ylabel("Latitude", fontsize=fontz)
                ax.set_xlabel("Longitude", fontsize=fontz)

        fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top, hspace=hspace, wspace=wspace)

        if one_colorbar_per_row:
            for r, pk in enumerate(product_order):
                row_height = (top - bottom - hspace * (nrows - 1)) / nrows
                row_bottom = top - (r + 1) * row_height - r * hspace
                cax = fig.add_axes([left, row_bottom - 0.028, right - left, 0.018])
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
        cbar_label: str = "EOF pattern (units as provided)",
        extent_override: Optional[Tuple[float, float, float, float]] = None,
    ):
        """
        Plot Obs + each model in a (nrows x ncols) grid with per-panel r/RMSE annotations.

        region_* provided → crop maps and compute stats over that region.
        region_* None     → global plotting + global stats.
        """
        panels: List[Tuple[str, xr.DataArray]] = [("Obs", obs_map)]
        for i in range(model_stack.sizes[member_dim]):
            panels.append((member_labels[i], model_stack.isel({member_dim: i})))

        n_panels = len(panels)
        ncols = max(1, int(ncols))
        nrows = int(math.ceil(n_panels / ncols))

        panels_plot: List[Tuple[str, xr.DataArray]] = []
        if region_lat_bounds is None and region_lon_bounds is None:
            obs0 = self._normalize_lon_for_bounds(obs_map, "0_360")
            panels_plot.append(("Obs", obs0))
            for i in range(model_stack.sizes[member_dim]):
                da = self._normalize_lon_for_bounds(model_stack.isel({member_dim: i}), "0_360")
                panels_plot.append((member_labels[i], da))
        else:
            for name, da in panels:
                panels_plot.append((name, self._subset_latlon(da, region_lat_bounds, region_lon_bounds, lon_convention)))

        obs_plot = panels_plot[0][1]
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
            ax.coastlines(linewidth=0.6)
            ax.set_extent(extent, crs=data_crs)
            ax.set_xticks(xticks, crs=data_crs)
            ax.set_yticks(yticks, crs=data_crs)
            ax.xaxis.set_major_formatter(LongitudeFormatter(".0f"))
            ax.yaxis.set_major_formatter(LatitudeFormatter(".0f"))
            ax.tick_params(labelsize=fontz * 0.85)
            ax.set_title(f"({chr(97 + p)}) {name}", loc="left", fontsize=fontz)

            if p > 0:
                r0, rmse0 = self._weighted_corr_rmse(panels_plot[0][1].values, da.values, w2d=w2d)
                ax.text(
                    0.98, 0.98, f"r = {r0:.2f}\nRMSE = {rmse0:.2f}",
                    transform=ax.transAxes, ha="right", va="top", fontsize=fontz * 0.9,
                    bbox=dict(facecolor="white", edgecolor="0.7", alpha=0.9, boxstyle="round,pad=0.25"),
                )

            if (p % ncols) == 0:
                ax.set_ylabel("Latitude", fontsize=fontz)
            ax.set_xlabel("Longitude", fontsize=fontz)

        fig.subplots_adjust(left=0.06, right=0.98, top=0.92, bottom=0.12, wspace=wspace, hspace=hspace)

        cbar = fig.colorbar(im_last, ax=axes, orientation="horizontal", fraction=0.05, pad=cb_pad, aspect=45)
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
        **kwargs,
    ):
        return self.plot_multimodel_panel_with_stats(
            title=f"{mode} pattern — {token}",
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
        **kwargs,
    ):
        return self.plot_multimodel_panel_with_stats(
            title=f"{mode} teleconnection — {token}",
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
