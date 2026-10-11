"""Soil-moisture time series at selected points: model ensemble vs observations.

``SoilPointTimeSeriesFromCollectedFiles`` extracts observed and per-member
model series at named points from the file lists of
:mod:`util.initial_lnd_file_collection`, caches them as NetCDF
(``save_collected``) and plots member spaghetti or mean +/- spread.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


@dataclass(frozen=True)
class SoilPointTSConfig:
    model_var_name: str = "SOILWATER_10CM"  # already kg/m2 in your setup
    obs_var_name: str = "H2OSOI"            # m3/m3 in ESA_CCI
    depth_cm: float = 10.0
    water_density: float = 1000.0

    lat_name: str = "lat"
    lon_name: str = "lon"
    time_name: str = "time"

    # land mask
    landmask_file: str | None = None
    landmask_var: str = "landfrac"
    landfrac_threshold: float = 0.1
    mask_land: bool = True


class SoilPointTimeSeriesFromCollectedFiles:
    """
    Input expected (for one var):
      var_files[exp]["obs"]   -> list[str] (obs files)
      var_files[exp]["model"] -> dict[ens]-> list[str] (model files)
    """

    def __init__(self, cfg: SoilPointTSConfig):
        self.cfg = cfg
        self.landmask = None
        if cfg.landmask_file and os.path.exists(cfg.landmask_file):
            ds_mask = xr.open_dataset(cfg.landmask_file)
            self.landmask = ds_mask[cfg.landmask_var]
            self.landmask = self._to_0_360(self.landmask)

    def save_collected(
        self,
        var_files: dict,
        *,
        points: Dict[str, dict],
        obs_key: str,
        exp_order: List[str],
        out_nc: str | Path,
        overwrite: bool = False,
    ) -> Path:
        """
        Save a compact point dataset to NetCDF:
          obs(point,time)
          model(exp,ensemble,point,time)
        """
        out_nc = Path(out_nc)
        out_nc.parent.mkdir(parents=True, exist_ok=True)
        if out_nc.exists() and not overwrite:
            return out_nc
    
        # obs from the first exp (same obs list for all exps)
        obs_ts = self._load_obs_timeseries_at_points(var_files[exp_order[0]]["obs"], points)
    
        # model spaghetti: keep ensemble dimension
        per_exp = []
        for exp in exp_order:
            ens_ts = self._load_model_ensemble_timeseries_at_points(var_files[exp]["model"], points)
            ens_ts, obs_aligned = xr.align(ens_ts, obs_ts, join="inner")  # align time
    
            ds_exp = xr.Dataset(
                data_vars=dict(
                    model=ens_ts,     # (ensemble, point, time)
                ),
                coords=dict(
                    time=ens_ts["time"],
                    point=ens_ts["point"],
                    ensemble=ens_ts["ensemble"],
                    point_name=obs_aligned["point_name"],
                ),
                attrs=dict(exp=exp, obs=obs_key, units="kg/m2"),
            ).expand_dims(exp=[exp])
            per_exp.append(ds_exp)
    
        ds_model = xr.concat(per_exp, dim="exp")
        ds_out = xr.Dataset(
            data_vars=dict(
                obs=obs_ts,                     # (point,time)
                model=ds_model["model"],         # (exp,ensemble,point,time)
            ),
            coords=dict(
                time=ds_model["time"],
                exp=ds_model["exp"],
                ensemble=ds_model["ensemble"],
                point=ds_model["point"],
                point_name=ds_model["point_name"],
            ),
            attrs=dict(
                obs=obs_key,
                model_var=self.cfg.model_var_name,
                obs_var=self.cfg.obs_var_name,
                depth_cm=float(self.cfg.depth_cm),
                units="kg/m2",
            ),
        )
    
        # compression helps a lot (esp. model exp×ens×time)
        encoding = {}
        for v in ds_out.data_vars:
            encoding[v] = dict(zlib=True, complevel=4)
    
        ds_out.to_netcdf(out_nc, encoding=encoding)
        return out_nc
    
    
    @staticmethod
    def open_collected(nc_file: str | Path) -> xr.Dataset:
        return xr.open_dataset(nc_file)
    
        
    def compute(
        self,
        var_files: dict,
        *,
        points: Dict[str, dict],
        obs_key: str,
        exp_order: List[str] | None = None,
    ) -> xr.Dataset:
        if exp_order is None:
            exp_order = list(var_files.keys())

        obs_ts = self._load_obs_timeseries_at_points(var_files[exp_order[0]]["obs"], points)

        per_exp = []
        for exp in exp_order:
            model_by_ens = var_files[exp]["model"]
            ens_ts = self._load_model_ensemble_timeseries_at_points(model_by_ens, points)

            # align time with obs
            ens_ts, obs_ts_aligned = xr.align(ens_ts, obs_ts, join="inner")

            model_mean = ens_ts.mean("ensemble")
            model_std = ens_ts.std("ensemble")
            bias = model_mean - obs_ts_aligned

            ds_exp = xr.Dataset(
                data_vars=dict(
                    obs=obs_ts_aligned,
                    model_mean=model_mean,
                    model_std=model_std,
                    bias=bias,
                ),
                coords=dict(
                    point=obs_ts_aligned["point"],
                    time=obs_ts_aligned["time"],
                    point_name=obs_ts_aligned["point_name"],
                ),
                attrs=dict(exp=exp, obs=obs_key, units="kg/m2"),
            ).expand_dims(exp=[exp])

            per_exp.append(ds_exp)

        return xr.concat(per_exp, dim="exp")

    def plot_spaghetti(
        self,
        ds: xr.Dataset,
        *,
        exps: list[str] | None = None,
        out_png: str | Path | None = None,
        show: bool = True,
        title: str | None = None,
        max_members: int | None = None,
        dpi: int = 200,
    ):
        if exps is None:
            exps = list(ds["exp"].values)
    
        # split exps for styling
        capt_exps = [e for e in exps if str(e).startswith("CAPT")]
        dart_exps = [e for e in exps if str(e).startswith("DART")]
        other_exps = [e for e in exps if (e not in capt_exps and e not in dart_exps)]
    
        point_keys = list(ds["point"].values)
        fig, axes = plt.subplots(
            nrows=len(point_keys),
            ncols=1,
            figsize=(12, 3.6 * len(point_keys)),
            sharex=True,
        )
        if len(point_keys) == 1:
            axes = [axes]
    
        for ax, pk in zip(axes, point_keys):
            # -----------------------------
            # OBS: RED dashed, thick
            # -----------------------------
            obs0 = ds["obs"].sel(point=pk)
            ax.plot(
                ds["time"].values,
                obs0.values,
                color="red",
                linestyle="--",
                linewidth=2.8,
                label="OBS",
                zorder=5,
            )
    
            # helper to plot spaghetti for a list of experiments with one common style
            def _plot_spaghetti_for(explist, *, color, alpha, lw, legend_label):
                did_legend = False
                for exp in explist:
                    ens = ds["model"].sel(exp=exp, point=pk)  # (ensemble,time)
    
                    if max_members is not None and ens.sizes["ensemble"] > max_members:
                        ens = ens.isel(ensemble=slice(0, max_members))
    
                    # draw each member
                    for i in range(ens.sizes["ensemble"]):
                        ax.plot(
                            ds["time"].values,
                            ens.isel(ensemble=i).values,
                            color=color,
                            alpha=alpha,
                            linewidth=lw,
                            zorder=2,
                        )
    
                    # one legend proxy per category (not per exp)
                    if not did_legend:
                        ax.plot([], [], color=color, linewidth=1.6, label=legend_label)
                        did_legend = True
    
            # -----------------------------
            # CAPT: GREY
            # -----------------------------
            _plot_spaghetti_for(
                capt_exps,
                color="0.6",         # grey
                alpha=0.35,
                lw=1.5,
                legend_label="CAPT members",
            )
    
            # -----------------------------
            # DART: BLUE
            # -----------------------------
            _plot_spaghetti_for(
                dart_exps,
                color="tab:blue",    # blue
                alpha=0.45,
                lw=0.8,
                legend_label="DART members",
            )
    
            # -----------------------------
            # Any other experiments (optional): light black
            # -----------------------------
            _plot_spaghetti_for(
                other_exps,
                color="0.2",
                alpha=0.25,
                lw=0.6,
                legend_label="Other members",
            )
    
            pname = str(ds["point_name"].sel(point=pk).values) if "point_name" in ds.coords else str(pk)
            ax.set_title(pname, fontsize=14)
            ax.set_ylabel("Soil moisture (kg m$^{-2}$)")
            ax.grid(True, alpha=0.3)
    
            # cleaner legend: OBS + CAPT + DART (and Other only if present)
            ax.legend(loc="upper right", frameon=False, ncol=3)
    
        axes[-1].set_xlabel("Time")
    
        if title:
            fig.suptitle(title, y=1.01, fontsize=16)
    
        fig.tight_layout()
    
        if out_png is not None:
            out_png = Path(out_png)
            out_png.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    
        if show:
            plt.show()
        else:
            plt.close(fig)
    
        return fig


    def plot(
        self,
        ds: xr.Dataset,
        *,
        exps: List[str],
        point_keys: List[str] | None = None,
        out_png: str | Path | None = None,
        show: bool = True,
        title: str | None = None,
        dpi: int = 200,
    ):
        if point_keys is None:
            point_keys = list(ds["point"].values)

        fig, axes = plt.subplots(
            nrows=len(point_keys), ncols=1,
            figsize=(12, 3.6 * len(point_keys)),
            sharex=True,
        )
        if len(point_keys) == 1:
            axes = [axes]

        for ax, pk in zip(axes, point_keys):
            # OBS (same for all exp after align)
            obs0 = ds["obs"].sel(point=pk, exp=exps[0])
            ax.plot(obs0["time"].values, obs0.values, linestyle="--", linewidth=2.0, label="OBS")

            # each experiment mean + spread
            for exp in exps:
                mm = ds["model_mean"].sel(point=pk, exp=exp)
                ss = ds["model_std"].sel(point=pk, exp=exp)

                ax.plot(mm["time"].values, mm.values, label=f"{exp} mean")
                ax.fill_between(
                    mm["time"].values,
                    (mm - ss).values,
                    (mm + ss).values,
                    alpha=0.15,
                )

            pname = str(ds["point_name"].sel(point=pk).values)
            ax.set_title(pname, fontsize=14)
            ax.set_ylabel("Soil moisture (kg m$^{-2}$)")
            ax.grid(True, alpha=0.3)

        axes[-1].set_xlabel("Time")
        axes[0].legend(loc="upper right", ncol=2)

        if title:
            fig.suptitle(title, y=1.01, fontsize=16)

        fig.tight_layout()

        if out_png is not None:
            out_png = Path(out_png)
            out_png.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(out_png, dpi=dpi, bbox_inches="tight")

        if show:
            plt.show()
        else:
            plt.close(fig)

        return fig

    # ---------- helpers ----------
    def _open_mfd(self, files: list[str]) -> xr.Dataset:
        missing = [f for f in files if not os.path.exists(f)]
        if missing:
            raise FileNotFoundError(f"Missing files (showing up to 5): {missing[:5]}")
        return xr.open_mfdataset(files, combine="by_coords")

    def _to_0_360(self, da_or_ds):
        lon_name = self.cfg.lon_name
        if lon_name not in da_or_ds.coords:
            return da_or_ds
        lon = da_or_ds[lon_name]
        if lon.min() < 0 or lon.max() > 360:
            da_or_ds = da_or_ds.assign_coords({lon_name: (lon % 360)}).sortby(lon_name)
        return da_or_ds

    def _load_obs_timeseries_at_points(self, obs_files: list[str], points: Dict[str, dict]) -> xr.DataArray:
        ds = self._open_mfd(obs_files)
        ds = self._to_0_360(ds)

        obs = ds[self.cfg.obs_var_name]

        # mask missing
        fill = obs.attrs.get("_FillValue", None)
        if fill is None:
            fill = obs.attrs.get("missing_value", -9999)
        obs = obs.where(obs != fill)

        # convert m3/m3 -> kg/m2
        obs = obs * self.cfg.water_density * (self.cfg.depth_cm / 100.0)
        obs.attrs["units"] = "kg/m2"

        if self.cfg.mask_land and (self.landmask is not None):
            obs = obs.where(self.landmask > self.cfg.landfrac_threshold)

        series = []
        pkeys, pnames = [], []
        for k, p in points.items():
            da = obs.sel(
                {self.cfg.lat_name: float(p["lat"]), self.cfg.lon_name: float(p["lon"])},
                method="nearest",
            )
            series.append(da)
            pkeys.append(k)
            pnames.append(p.get("name", k))

        out = xr.concat(series, dim="point")
        out = out.assign_coords(point=np.array(pkeys, dtype=str))
        out = out.rename({self.cfg.time_name: "time"})
        out = out.assign_coords(point_name=("point", np.array(pnames, dtype=str)))
        return out

    def _load_model_ensemble_timeseries_at_points(self, model_by_ens: dict[str, list[str]], points: Dict[str, dict]) -> xr.DataArray:
        ens_series = []
        ens_labels = []

        for ens, files in model_by_ens.items():
            ds = self._open_mfd(files)
            ds = self._to_0_360(ds)
            mod = ds[self.cfg.model_var_name]

            if self.cfg.mask_land and (self.landmask is not None):
                mod = mod.where(self.landmask > self.cfg.landfrac_threshold)

            series = []
            pkeys = []
            for k, p in points.items():
                da = mod.sel(
                    {self.cfg.lat_name: float(p["lat"]), self.cfg.lon_name: float(p["lon"])},
                    method="nearest",
                )
                series.append(da)
                pkeys.append(k)

            da_pts = xr.concat(series, dim="point").assign_coords(point=np.array(pkeys, dtype=str))
            ens_series.append(da_pts)
            ens_labels.append(ens)

        out = xr.concat(ens_series, dim="ensemble").assign_coords(ensemble=np.array(ens_labels, dtype=str))
        out = out.rename({self.cfg.time_name: "time"})
        return out
