from __future__ import annotations

from pathlib import Path
from typing import Callable, Mapping, Optional, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import cartopy.crs as ccrs
import cartopy.feature as cfeature


class ClimatologyMapPlotter:
    """Plot seasonal climatology and model-reference bias maps for selected variables."""

    def __init__(
        self,
        exp_dict: Mapping[str, Mapping[str, str | Path]],
        save_dir: Optional[str | Path] = None,
        ref_name: str = "ERA5",
        unit_map: Optional[Mapping[str, str]] = None,
        scale_map: Optional[Mapping[str, float | Callable[[xr.DataArray], xr.DataArray]]] = None,
        season_order: Optional[Sequence[str]] = None,
    ) -> None:
        self.exp_dict = {
            exp: {var: Path(path) for var, path in paths.items()}
            for exp, paths in exp_dict.items()
        }
        self.save_dir = Path(save_dir) if save_dir is not None else None
        self.ref_name = ref_name
        self.season_order = list(season_order or ["ANN", "DJF", "MAM", "JJA", "SON"])
        self.unit_map = {
            "pr": "pr",
            "psl": "psl",
            "tas": "tas",
            "tauu": r"tauu (10$^{-2}$ N m$^{-2}$)",
            "tauv": r"tauv (10$^{-2}$ N m$^{-2}$)",
            "tau_mag": r"tau_mag (10$^{-2}$ N m$^{-2}$)",
        }
        if unit_map:
            self.unit_map.update(unit_map)
        self.scale_map: dict[str, float | Callable[[xr.DataArray], xr.DataArray]] = {
            "tauu": 100.0,
            "tauv": 100.0,
        }
        if scale_map:
            self.scale_map.update(scale_map)

    def _load_variable(self, path: Path, var_name: str) -> xr.DataArray:
        ds = xr.open_dataset(path, decode_times=True)
        data_var = var_name if var_name in ds.data_vars else self._find_data_variable(ds, var_name)
        da = ds[data_var]
        da = self._mask_bad_values(da)
        da = self._normalize_longitude(da)
        da = self._standardize_units(da, var_name)
        scale = self.scale_map.get(var_name)
        if callable(scale):
            da = scale(da)
        elif scale is not None:
            da = da * scale
        return da

    @staticmethod
    def _mask_bad_values(da: xr.DataArray) -> xr.DataArray:
        da = da.where(np.isfinite(da))
        da = da.where(np.abs(da) < 1.0e20)

        valid_range = da.attrs.get("valid_range")
        if valid_range is not None and len(valid_range) == 2:
            da = da.where((da >= valid_range[0]) & (da <= valid_range[1]))
        return da

    @staticmethod
    def _standardize_units(da: xr.DataArray, var_name: str) -> xr.DataArray:
        units = str(da.attrs.get("units", "")).lower().replace(" ", "")
        if var_name == "pr" and units in {
            "kgm-2s-1",
            "kg/m2/s",
            "kgm**-2s**-1",
            "kgm^-2s^-1",
        }:
            da = da * 86400.0
            da.attrs = dict(da.attrs)
            da.attrs["units"] = "mm day-1"
        elif var_name == "psl" and units in {"pa", "pascal", "pascals"}:
            da = da / 100.0
            da.attrs = dict(da.attrs)
            da.attrs["units"] = "hPa"
        return da

    @staticmethod
    def _normalize_longitude(da: xr.DataArray) -> xr.DataArray:
        if "lon" in da.coords:
            da = da.assign_coords(lon=(((da.lon + 180) % 360) - 180))
            da = da.sortby("lon")
        return da

    @staticmethod
    def _reference_on_model_grid(
        reference: xr.DataArray, model: xr.DataArray
    ) -> xr.DataArray:
        if all(dim in reference.coords and dim in model.coords for dim in ("lat", "lon")):
            if reference.sizes.get("lat") != model.sizes.get("lat") or reference.sizes.get(
                "lon"
            ) != model.sizes.get("lon"):
                return reference.interp(lat=model["lat"], lon=model["lon"])
        return reference

    @staticmethod
    def _area_weights_like(da: xr.DataArray) -> xr.DataArray:
        weights = np.cos(np.deg2rad(da["lat"]))
        weights = weights / weights.mean()
        return weights.broadcast_like(da)

    @staticmethod
    def _finite_pair_weights(
        model: xr.DataArray, reference: xr.DataArray, weights: xr.DataArray
    ) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
        valid = np.isfinite(model) & np.isfinite(reference) & np.isfinite(weights)
        return model.where(valid), reference.where(valid), weights.where(valid)

    @classmethod
    def _weighted_rmse(
        cls, model: xr.DataArray, reference: xr.DataArray, weights: xr.DataArray
    ) -> float:
        model, reference, weights = cls._finite_pair_weights(model, reference, weights)
        weight_sum = weights.sum(dim=["lat", "lon"], skipna=True)
        weight_sum_value = float(weight_sum)
        if not np.isfinite(weight_sum_value) or weight_sum_value == 0:
            return float("nan")

        mse = (((model - reference) ** 2) * weights).sum(
            dim=["lat", "lon"], skipna=True
        ) / weight_sum
        return float(np.sqrt(mse))

    @classmethod
    def _weighted_pcor(
        cls, model: xr.DataArray, reference: xr.DataArray, weights: xr.DataArray
    ) -> float:
        model, reference, weights = cls._finite_pair_weights(model, reference, weights)
        weight_sum = weights.sum(dim=["lat", "lon"], skipna=True)
        weight_sum_value = float(weight_sum)
        if not np.isfinite(weight_sum_value) or weight_sum_value == 0:
            return float("nan")

        model_mean = (model * weights).sum(dim=["lat", "lon"], skipna=True) / weight_sum
        ref_mean = (reference * weights).sum(dim=["lat", "lon"], skipna=True) / weight_sum
        model_anom = model - model_mean
        ref_anom = reference - ref_mean

        covariance = (model_anom * ref_anom * weights).sum(
            dim=["lat", "lon"], skipna=True
        ) / weight_sum
        model_variance = (model_anom**2 * weights).sum(
            dim=["lat", "lon"], skipna=True
        ) / weight_sum
        ref_variance = (ref_anom**2 * weights).sum(
            dim=["lat", "lon"], skipna=True
        ) / weight_sum

        denominator = float(np.sqrt(model_variance * ref_variance))
        if not np.isfinite(denominator) or denominator == 0:
            return float("nan")
        return float(covariance / denominator)

    def _load_components(
        self, tauu_path: Path, tauv_path: Path
    ) -> tuple[xr.DataArray, xr.DataArray]:
        tauu = self._load_variable(tauu_path, "tauu")
        tauv = self._load_variable(tauv_path, "tauv")
        return tauu, tauv

    @staticmethod
    def _find_data_variable(ds: xr.Dataset, var_name: str) -> str:
        candidates = [
            name
            for name, da in ds.data_vars.items()
            if name.lower() == var_name.lower() and "time" in da.dims
        ]
        if not candidates:
            candidates = [
                name
                for name, da in ds.data_vars.items()
                if var_name.lower() in name.lower() and "time" in da.dims
            ]
        if not candidates:
            candidates = [
                name for name, da in ds.data_vars.items() if "time" in da.dims
            ]
        if not candidates:
            raise ValueError(f"Could not find a time-dependent variable for {var_name!r}.")
        return candidates[0]

    @staticmethod
    def _find_tau_variable(ds: xr.Dataset) -> str:
        candidates = [
            name
            for name, da in ds.data_vars.items()
            if "tau" in name.lower() and da.dims and da.dims[0] == "time"
        ]
        if not candidates:
            raise ValueError("Could not find a time-dependent tau variable in dataset.")
        return candidates[0]

    def _compute_seasonal_means(self, da: xr.DataArray) -> dict[str, xr.DataArray]:
        da = da.copy()
        da["time"] = xr.cftime_range(start="2000-01-01", periods=da.sizes["time"], freq="MS")
        seasonal = {"ANN": da.mean(dim="time", skipna=True)}
        seasonal.update(
            {
                season: da.sel(time=da["time"].dt.season == season).mean(
                    dim="time", skipna=True
                )
                for season in ["DJF", "MAM", "JJA", "SON"]
            }
        )
        return seasonal

    def _compute_all_fields(self, var_name: str) -> dict[str, dict[str, xr.DataArray]]:
        seasonal_data: dict[str, dict[str, xr.DataArray]] = {}
        for exp, paths in self.exp_dict.items():
            if var_name == "tau_mag":
                tauu, tauv = self._load_components(paths["tauu"], paths["tauv"])
                tauu_season = self._compute_seasonal_means(tauu)
                tauv_season = self._compute_seasonal_means(tauv)
                seasonal_data[exp] = {
                    season: np.sqrt(tauu_season[season] ** 2 + tauv_season[season] ** 2)
                    for season in self.season_order
                }
                continue

            if var_name not in paths:
                raise ValueError(f"Missing path for variable {var_name!r} in experiment {exp!r}.")
            seasonal_data[exp] = self._compute_seasonal_means(
                self._load_variable(paths[var_name], var_name)
            )
        return seasonal_data

    @staticmethod
    def _add_horizontal_colorbar(
        fig: plt.Figure,
        cax: plt.Axes,
        im: mpl.cm.ScalarMappable,
        label: str,
        font_size: float,
        label_position: str = "bottom",
    ) -> None:
        cbar = fig.colorbar(im, cax=cax, orientation="horizontal", extend="both")
        if hasattr(im, "levels"):
            cbar.set_ticks(im.levels)
        labelpad = 8 if label_position == "top" else 4
        cbar.set_label(label, fontsize=font_size, labelpad=labelpad)
        cbar.ax.xaxis.set_label_position(label_position)
        cbar.ax.xaxis.set_ticks_position(label_position)
        cbar.ax.tick_params(labelsize=font_size * 0.95)

    def plot_variable(
        self,
        var_name: str = "tau_mag",
        cmap1: str = "viridis",
        cmap2: str = "RdBu_r",
        levels1: Optional[Sequence[float]] = None,
        levels2: Optional[Sequence[float]] = None,
        center_zero: bool = True,
        font_size: float = 11,
        fig_size: tuple[float, float] = (20, 12),
        dpi: int = 150,
        show_borders: bool = False,
        ref_label: Optional[str] = None,
        metric_box_alpha: float = 0.65,
        panel_layout: str = "model_rows",
    ) -> None:
        ref = self.ref_name
        if ref not in self.exp_dict:
            raise ValueError(f"Reference {ref!r} is not in exp_dict.")
        ref_display = ref_label or ref

        all_fields = self._compute_all_fields(var_name)
        ref_data = all_fields[ref]
        model_exps = [exp for exp in self.exp_dict if exp != ref]
        panel_layout = panel_layout.lower()
        if panel_layout not in {"model_rows", "season_rows"}:
            raise ValueError("panel_layout must be 'model_rows' or 'season_rows'.")

        if panel_layout == "season_rows":
            map_nrow = len(self.season_order)
            map_ncol = len(model_exps) + 1
            grid_height_ratios = [1.0] * map_nrow + [0.10, 0.10]
            map_grid_rows = list(range(map_nrow))
            ref_cbar_spec = (map_nrow, slice(None))
            diff_cbar_spec = (map_nrow + 1, slice(None)) if model_exps else None
            hspace = 0.16
        else:
            map_nrow = len(model_exps) + 1
            map_ncol = len(self.season_order)
            if model_exps:
                grid_height_ratios = [1.0, 0.18, *([1.0] * len(model_exps)), 0.10]
            else:
                grid_height_ratios = [1.0, 0.12]
            map_grid_rows = [0] + [i + 1 for i in range(1, map_nrow)]
            ref_cbar_spec = (1, slice(None))
            diff_cbar_spec = (-1, slice(None)) if model_exps else None
            hspace = 0.10

        fig = plt.figure(figsize=fig_size)
        grid = fig.add_gridspec(
            len(grid_height_ratios),
            map_ncol,
            height_ratios=grid_height_ratios,
            left=0.05,
            right=0.95,
            top=0.94,
            bottom=0.10,
            wspace=0.18,
            hspace=hspace,
        )

        axes = np.empty((map_nrow, map_ncol), dtype=object)
        for i, grid_row in enumerate(map_grid_rows):
            for j in range(map_ncol):
                axes[i, j] = fig.add_subplot(grid[grid_row, j], projection=ccrs.PlateCarree())

        ref_cax = fig.add_subplot(grid[ref_cbar_spec])
        diff_cax = fig.add_subplot(grid[diff_cbar_spec]) if diff_cbar_spec else None

        ref_im = None
        diff_im = None
        for i in range(map_nrow):
            for j in range(map_ncol):
                if panel_layout == "season_rows":
                    season = self.season_order[i]
                    exp = ref if j == 0 else model_exps[j - 1]
                else:
                    exp = ([ref] + model_exps)[i]
                    season = self.season_order[j]

                if exp == ref:
                    plot_data = ref_data[season]
                    cmap = cmap1
                    levels = levels1
                else:
                    model_data = all_fields[exp][season]
                    ref_on_model_grid = self._reference_on_model_grid(ref_data[season], model_data)
                    plot_data = model_data - ref_on_model_grid
                    cmap = cmap2
                    levels = levels2

                norm = None
                kwargs = {}
                if levels is not None:
                    norm = mpl.colors.BoundaryNorm(
                        boundaries=levels, ncolors=plt.get_cmap(cmap).N
                    )
                    kwargs.update({"levels": levels, "vmin": levels[0], "vmax": levels[-1]})
                else:
                    kwargs["levels"] = 20

                ax = axes[i, j]
                ax.set_global()
                ax.coastlines()
                if show_borders:
                    ax.add_feature(cfeature.BORDERS, linewidth=0.3)
                gl = ax.gridlines(draw_labels=True, linewidth=0.2, alpha=0.5)
                gl.top_labels = False
                gl.right_labels = False
                gl.left_labels = True
                gl.bottom_labels = True
                gl.xlabel_style = {"size": font_size * 0.85}
                gl.ylabel_style = {"size": font_size * 0.85}

                im = ax.contourf(
                    plot_data["lon"].values,
                    plot_data["lat"].values,
                    plot_data,
                    transform=ccrs.PlateCarree(),
                    cmap=cmap,
                    norm=norm,
                    extend="both",
                    **kwargs,
                )
                if exp == ref and ref_im is None:
                    ref_im = im
                elif exp != ref and diff_im is None:
                    diff_im = im

                if exp == ref:
                    ax.set_title(f"{ref_display} ({season})", fontsize=font_size)
                else:
                    ax.set_title(f"{exp} - {ref_display} ({season})", fontsize=font_size)
                    model_data = all_fields[exp][season]
                    ref_on_model_grid = self._reference_on_model_grid(ref_data[season], model_data)
                    weights = self._area_weights_like(model_data)
                    rmse = self._weighted_rmse(model_data, ref_on_model_grid, weights)
                    pcor = self._weighted_pcor(model_data, ref_on_model_grid, weights)
                    ax.text(
                        0.98,
                        0.02,
                        f"RMSE={rmse:.2f}\nPCOR={pcor:.2f}",
                        transform=ax.transAxes,
                        fontsize=font_size * 0.90,
                        va="bottom",
                        ha="right",
                        bbox=dict(
                            facecolor=(1.0, 1.0, 1.0, metric_box_alpha),
                            edgecolor="black",
                            boxstyle="round,pad=0.3",
                        ),
                    )

        unit_label = self.unit_map.get(var_name, var_name)
        if ref_im is not None:
            self._add_horizontal_colorbar(
                fig,
                ref_cax,
                ref_im,
                unit_label,
                font_size,
                label_position="top",
            )
        else:
            ref_cax.set_visible(False)

        if diff_im is not None and diff_cax is not None:
            diff_label = f"Bias of {unit_label}" if center_zero else unit_label
            self._add_horizontal_colorbar(
                fig,
                diff_cax,
                diff_im,
                diff_label,
                font_size,
            )
        elif diff_cax is not None:
            diff_cax.set_visible(False)

        if self.save_dir:
            self.save_dir.mkdir(parents=True, exist_ok=True)
            out_file = self.save_dir / f"climatology_bias_grid_{var_name}.png"
            plt.savefig(out_file, dpi=dpi)
            print(f"[INFO] Saved {out_file}")

        plt.show()
