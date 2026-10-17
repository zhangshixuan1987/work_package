"""Monthly stratospheric polar-vortex (SSW) indices for the AMIP hindcasts.

From monthly ``ua``, ``ta``, ``zg`` (model CMIP output and reanalyses on the
model pressure levels):

* ``u_55_65N``: cos-lat weighted zonal-mean U over 55-65N (the u60 index);
* ``u_vortex``: the same over 55-75N;
* ``z_pcap`` / ``t_pcap``: polar-cap (60-90N) zonal-mean Z and T.

``compare_to_reference`` adds, per index, the bias series
``<idx>_bias_<REF>(time, plev)`` and RMSE profiles over all months
(``rmse_all``), winter months (``rmse_winter``, DJF), the SSW season
(``rmse_ssw``, Nov 2012-Mar 2013 by default) and lead-month windows counted
from the first month (``rmse_leadXX_YY``; Jan 2012 = lead 0, inclusive).
Files: ``ssw_monthly_<SOURCE>_vs_<REF>.nc``; a reference's own indices are
``ssw_monthly_<REF>_vs_<REF>.nc``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import xarray as xr

from configs.hind_experiment_config import (
    DEFAULT_INPUT_DATA_ROOT,
    END_DATE,
    MODEL_PLEV_PA,
    START_DATE,
    check_experiments,
)
from util.hind_model_data import open_model_monthly, select_plev
from util.hind_observation_data import (
    DEFAULT_OBSERVATION_ROOT,
    REANALYSES,
    open_reference_monthly,
)

SSW_FIELDS = ("ua", "ta", "zg")
SSW_INDICES = ("u_55_65N", "u_vortex", "z_pcap", "t_pcap")


@dataclass(frozen=True)
class SswConfig:
    plev_pa: tuple = MODEL_PLEV_PA
    polar_cap: tuple[float, float] = (60.0, 90.0)
    u_lat_band: tuple[float, float] = (55.0, 65.0)
    vortex_lat_band: tuple[float, float] = (55.0, 75.0)
    winter_months: tuple[int, ...] = (12, 1, 2)
    ssw_start: tuple[int, int] = (2012, 11)   # (year, month), inclusive
    ssw_end: tuple[int, int] = (2013, 3)
    lead_windows: tuple[tuple[int, int], ...] = field(default_factory=lambda: ((0, 3), (0, 6), (10, 14)))


def lead_key(window: tuple[int, int]) -> str:
    return f"lead{window[0]:02d}_{window[1]:02d}"


def ssw_path(data_dir, source: str, ref: str) -> Path:
    return Path(data_dir) / f"ssw_monthly_{source}_vs_{ref}.nc"


class SswAnalyzer:
    def __init__(self, cfg: SswConfig = SswConfig()):
        self.cfg = cfg

    @staticmethod
    def _cap_mean(zm: xr.DataArray, band: tuple[float, float]) -> xr.DataArray:
        sub = zm.sel(lat=slice(*band))
        return sub.weighted(np.cos(np.deg2rad(sub.lat))).mean("lat")

    def compute_indices(self, ua: xr.DataArray, ta: xr.DataArray, zg: xr.DataArray, name: str,
                        compute: bool = True) -> xr.Dataset:
        """Indices (time, plev[Pa]); fields are reduced lazily (lat >= 50N, zonal mean) then computed.

        ``compute=False`` returns the lazy indices, e.g. to compute many sources in one ``compute_all``.
        """
        c = self.cfg
        zm = {}
        for key, da in (("ua", ua), ("ta", ta), ("zg", zg)):
            da = select_plev(da, c.plev_pa).sel(lat=slice(50.0, 90.0))
            zm[key] = da.mean("lon") if "lon" in da.dims else da
        ds = xr.Dataset({
            f"u_{c.u_lat_band[0]:.0f}_{c.u_lat_band[1]:.0f}N": self._cap_mean(zm["ua"], c.u_lat_band),
            "u_vortex": self._cap_mean(zm["ua"], c.vortex_lat_band),
            "z_pcap": self._cap_mean(zm["zg"], c.polar_cap),
            "t_pcap": self._cap_mean(zm["ta"], c.polar_cap),
        }).transpose("time", "plev")
        if compute:
            ds = ds.compute()
        ds.attrs = {
            "name": name,
            "polar_cap": f"{c.polar_cap[0]:g}-{c.polar_cap[1]:g}N",
            "u_lat_band": f"{c.u_lat_band[0]:g}-{c.u_lat_band[1]:g}N",
            "vortex_lat_band": f"{c.vortex_lat_band[0]:g}-{c.vortex_lat_band[1]:g}N",
        }
        return ds

    def _masks(self, time: xr.DataArray) -> dict[str, np.ndarray]:
        c = self.cfg
        years = np.array([t.year for t in time.values])
        months = np.array([t.month for t in time.values])
        key = years * 100 + months
        masks = {
            "all": np.ones(key.size, bool),
            "winter": np.isin(months, c.winter_months),
            "ssw": (key >= c.ssw_start[0] * 100 + c.ssw_start[1]) & (key <= c.ssw_end[0] * 100 + c.ssw_end[1]),
        }
        lead = np.arange(key.size)
        for w in c.lead_windows:
            masks[lead_key(w)] = (lead >= w[0]) & (lead <= w[1])
        return masks

    def compare_to_reference(self, ds: xr.Dataset, ref: xr.Dataset, label: str) -> xr.Dataset:
        exp, ref = xr.align(ds, ref, join="inner")
        out = exp.copy()
        masks = self._masks(out.time)
        for v in [v for v in exp.data_vars if v in ref.data_vars]:
            err = exp[v] - ref[v]
            out[f"{v}_bias_{label}"] = err
            for wname, m in masks.items():
                rmse = np.sqrt((err.isel(time=np.flatnonzero(m)) ** 2).mean("time")) if m.any() \
                    else xr.full_like(err.isel(time=0, drop=True), np.nan)
                out[f"{v}_rmse_{wname}_{label}"] = rmse
        c = self.cfg
        out.attrs.update(reference=label,
                         ssw_window=f"{c.ssw_start[0]}-{c.ssw_start[1]:02d}..{c.ssw_end[0]}-{c.ssw_end[1]:02d}",
                         lead_windows=", ".join(lead_key(w) for w in c.lead_windows),
                         lead_definition="lead month counted from the first month (Jan 2012 = 0), inclusive")
        return out


def open_ssw_input(
    source: str,
    *,
    data_root=DEFAULT_INPUT_DATA_ROOT,
    observation_root=DEFAULT_OBSERVATION_ROOT,
    start: str = START_DATE,
    end: str = END_DATE,
) -> dict[str, xr.DataArray]:
    """Lazy monthly ``ua``, ``ta``, ``zg`` of an experiment or reanalysis (``ERA5``, ``ERAI``)."""
    if source in REANALYSES:
        return {v: open_reference_monthly(source, v, start=start, end=end, observation_root=observation_root)
                for v in SSW_FIELDS}
    check_experiments([source])
    return {v: open_model_monthly(source, v, data_root=data_root, start=start, end=end) for v in SSW_FIELDS}
