"""NAO index of the AMIP hindcasts by projection on an observed EOF base.

1. ``compute_nao_eof_base``: EOF1 of NDJFM monthly SLP anomalies (NOAA-20C,
   1985-2014 base climatology) over 20-80N, 90W-40E, sqrt(cos lat) weighted,
   unit norm, sign fixed positive near the Azores (35N, 20W). PC1 is the
   projection of the anomalies on the EOF. The file holds the EOF, PC1, the
   explained variance, and the base monthly climatology (all 12 months) used
   for every later anomaly.
2. ``verify_nao_eof_base``: norm/variance checks and re-projection of the
   source data, which must reproduce the saved PC1.
3. ``project_daily``: daily SLP anomalies from the base monthly climatology,
   projected on the EOF and standardized with the base NDJFM PC1 mean/std.
4. ``window_means`` / ``window_correlation``: consecutive lead windows from
   the initial date (day 1 = init + 1); the correlation is across days
   within each window, for the single initialization.
"""
from __future__ import annotations

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from util.hind_model_data import to_datetime64


@dataclass(frozen=True)
class NaoConfig:
    lat_bounds: tuple[float, float] = (20.0, 80.0)
    lon_bounds: tuple[float, float] = (-90.0, 40.0)    # [-180, 180)
    months: tuple[int, ...] = (11, 12, 1, 2, 3)
    base_start: str = "1985-01-01"
    base_end: str = "2014-12-31"
    sign_point: tuple[float, float] = (35.0, -20.0)    # EOF positive here (Azores)
    window_days: int = 14
    max_lead_days: int = 84


def _domain(da: xr.DataArray, cfg: NaoConfig) -> xr.DataArray:
    da = da.rename({k: v for k, v in (("latitude", "lat"), ("longitude", "lon")) if k in da.dims})
    da = da.assign_coords(lon=((da.lon + 180.0) % 360.0) - 180.0).sortby("lon").sortby("lat")
    return da.sel(lat=slice(*cfg.lat_bounds), lon=slice(*cfg.lon_bounds))


def _weights(lat: xr.DataArray) -> xr.DataArray:
    return np.sqrt(np.cos(np.deg2rad(lat))).clip(min=0.0)


def _decoded_months(time) -> np.ndarray:
    return np.array([t.month for t in time.values]) if time.dtype == object else time.dt.month.values


def compute_nao_eof_base(psl_monthly: xr.DataArray, *, cfg: NaoConfig = NaoConfig(),
                         dataset_name: str = "NOAA-20C") -> xr.Dataset:
    """EOF1 base from monthly SLP (lazy input; reduced to the domain and NDJFM before computing)."""
    da = _domain(psl_monthly, cfg)
    clim = da.sel(time=slice(cfg.base_start, cfg.base_end)).groupby("time.month").mean("time").compute()
    da = da.sel(time=da.time.dt.month.isin(cfg.months))
    anom = (da.groupby("time.month") - clim).drop_vars("month", errors="ignore").compute()
    w = _weights(anom.lat)
    X = (anom * w).stack(space=("lat", "lon"))
    valid = np.isfinite(X).all("time").values
    Xv = X.values[:, valid]
    Xv = Xv - Xv.mean(axis=0)
    U, s, VT = np.linalg.svd(Xv, full_matrices=False)
    eof = np.full(X.sizes["space"], np.nan)
    eof[valid] = VT[0]
    eof_map = xr.DataArray(eof, coords={"space": X.space}, dims="space").unstack("space")
    norm = float(np.sqrt(np.nansum(eof_map.values ** 2)))
    eof_map = eof_map / norm
    if float(eof_map.sel(lat=cfg.sign_point[0], lon=cfg.sign_point[1], method="nearest")) < 0:
        eof_map = -eof_map
    # PC1 is the plain projection of the anomalies, the same operation as the daily index
    pc1 = (anom * w * eof_map.fillna(0)).sum(("lat", "lon")).values
    out = xr.Dataset(
        {
            "nao_eof1": eof_map.transpose("lat", "lon").assign_attrs(units="1", long_name="NAO EOF1 (weighted space, unit norm)"),
            "nao_pc1": ("time", pc1, {"units": "Pa", "long_name": "NAO PC1 (projection in weighted space)"}),
            "explained_variance_fraction": ((), float(s[0] ** 2 / np.sum(s ** 2))),
            "base_climatology": clim.assign_attrs(long_name="base-period monthly climatology (all months)"),
            "lat_weights_sqrtcos": w,
        },
        coords={"time": anom.time},
        attrs=dict(nao_definition="EOF1 of monthly NDJFM SLP anomalies over 20-80N, 90W-40E",
                   dataset_name=dataset_name, base_period_start=cfg.base_start, base_period_end=cfg.base_end,
                   months=",".join(map(str, cfg.months)), lat_bounds=str(cfg.lat_bounds),
                   lon_bounds=str(cfg.lon_bounds), weighting="sqrt(cos(lat))", normalized_unit_norm="True"),
    )
    return out


def verify_nao_eof_base(base: xr.Dataset, psl_monthly: xr.DataArray, *, cfg: NaoConfig = NaoConfig()) -> dict:
    """Checks: unit norm, 0 < variance fraction < 1, re-projection reproduces PC1."""
    da = _domain(psl_monthly, cfg)
    da = da.sel(time=da.time.dt.month.isin(cfg.months)).interp(lat=base.lat, lon=base.lon).compute()
    anom = da.groupby("time.month") - base.base_climatology
    pc = (anom * _weights(anom.lat) * base.nao_eof1.fillna(0)).sum(("lat", "lon")).drop_vars("month", errors="ignore")
    pc, saved = xr.align(pc, base.nao_pc1, join="inner")
    x, y = pc.values, saved.values
    corr = float(np.corrcoef(x, y)[0, 1])
    flags = []
    norm = float(np.sqrt(np.nansum(base.nao_eof1.values ** 2)))
    frac = float(base.explained_variance_fraction)
    if abs(norm - 1) > 5e-3:
        flags.append(f"EOF norm {norm:.4f} != 1")
    if not 0 < frac < 1:
        flags.append(f"explained variance fraction {frac}")
    if corr < 0.999:
        flags.append(f"re-projected PC1 correlation {corr:.5f} < 0.999")
    return dict(eof_norm=norm, explained_variance_fraction=frac, n_months=int(x.size), pc_corr=corr,
                pc_rmse=float(np.sqrt(np.mean((x - y) ** 2))), pc_std=float(y.std()), flags=flags,
                reprojected_pc=pc)


def plot_nao_eof_base(base: xr.Dataset, check: dict | None = None) -> plt.Figure:
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.2), gridspec_kw=dict(width_ratios=(1.5, 2, 1)))
    im = axes[0].pcolormesh(base.lon, base.lat, base.nao_eof1, cmap="RdBu_r", shading="auto",
                            vmin=-float(abs(base.nao_eof1).max()), vmax=float(abs(base.nao_eof1).max()))
    fig.colorbar(im, ax=axes[0], label="EOF1 (unit norm)")
    axes[0].set_title(f"(a) NAO EOF1 ({base.attrs.get('dataset_name', '')}, "
                      f"{float(base.explained_variance_fraction):.0%} variance)", loc="left")
    t = to_datetime64(base.time.values)
    axes[1].plot(t, base.nao_pc1, linewidth=0.8, label="saved PC1")
    if check is not None:
        axes[1].plot(to_datetime64(check["reprojected_pc"].time.values), check["reprojected_pc"], "--",
                     linewidth=0.8, label=f"re-projected (r={check['pc_corr']:.4f})")
        axes[1].legend(frameon=False)
    axes[1].set_title("(b) NAO PC1, NDJFM months", loc="left")
    axes[1].set_ylabel("Pa")
    axes[2].hist(base.nao_pc1, bins=40, color="0.6")
    axes[2].set_title("(c) PC1 distribution", loc="left")
    for ax in axes:
        ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig


class NaoProjector:
    """Daily NAO index on a fixed EOF base."""

    def __init__(self, base: xr.Dataset, cfg: NaoConfig = NaoConfig()):
        self.cfg = cfg
        self.base = base
        pc = base.nao_pc1.sel(time=slice(cfg.base_start, cfg.base_end))
        pc = pc.sel(time=pc.time.dt.month.isin(cfg.months))
        self.ref_mean, self.ref_std = float(pc.mean()), float(pc.std())

    def project_daily(self, psl_daily: xr.DataArray, name: str, compute: bool = True) -> xr.Dataset:
        """Daily ``nao_raw`` (Pa) and ``nao_std`` from daily SLP (time on datetime64 days).

        ``compute=False`` returns the lazy index, e.g. to compute many sources in one ``compute_all``.
        """
        da = _domain(psl_daily, self.cfg).interp(lat=self.base.lat, lon=self.base.lon)
        anom = da.groupby("time.month") - self.base.base_climatology
        raw = (anom * _weights(anom.lat) * self.base.nao_eof1.fillna(0)).sum(("lat", "lon"))
        raw = raw.drop_vars("month", errors="ignore")
        if compute:
            raw = raw.compute()
        return xr.Dataset({"nao_raw": raw.assign_attrs(units="Pa"),
                           "nao_std": ((raw - self.ref_mean) / self.ref_std).assign_attrs(units="1")},
                          attrs=dict(name=name, ref_mean=self.ref_mean, ref_std=self.ref_std,
                                     anomaly="daily SLP minus base monthly climatology",
                                     standardization="(raw - mean) / std of base NDJFM PC1"))

    def windows(self) -> list[tuple[int, int]]:
        w = self.cfg.window_days
        return [(1 + i * w, (i + 1) * w) for i in range(self.cfg.max_lead_days // w)]

    @staticmethod
    def _label(w: tuple[int, int]) -> str:
        return f"days {w[0]}-{w[1]}"

    def lead(self, series: xr.DataArray, init_date: str) -> np.ndarray:
        return ((series.time.values.astype("datetime64[D]") - np.datetime64(init_date, "D")) /
                np.timedelta64(1, "D")).astype(int)

    def window_means(self, series: xr.DataArray, init_date: str) -> xr.DataArray:
        """Mean over each lead window (NaN unless every day is present; Feb 29 is not
        required, since the noleap model calendar has none)."""
        lead = self.lead(series, init_date)
        t0 = np.datetime64(init_date, "D")
        vals = []
        for w in self.windows():
            days = t0 + np.arange(w[0], w[1] + 1)
            need = int(np.sum(np.char.find(days.astype(str), "-02-29") < 0))
            m = (lead >= w[0]) & (lead <= w[1]) & np.isfinite(series.values)
            vals.append(float(series.values[m].mean()) if m.sum() >= need else np.nan)
        return xr.DataArray(vals, dims="window", coords={"window": [self._label(w) for w in self.windows()]})

    def window_correlation(self, fcst: xr.DataArray, obs: xr.DataArray, init_date: str,
                           min_samples: int = 5) -> xr.Dataset:
        """Correlation across days in each lead window, on the days both series have."""
        f, o = xr.align(fcst, obs, join="inner")
        lead = self.lead(f, init_date)
        r, n = [], []
        for w in self.windows():
            m = (lead >= w[0]) & (lead <= w[1]) & np.isfinite(f.values) & np.isfinite(o.values)
            n.append(int(m.sum()))
            r.append(float(np.corrcoef(f.values[m], o.values[m])[0, 1]) if m.sum() >= min_samples else np.nan)
        labels = [self._label(w) for w in self.windows()]
        return xr.Dataset({"nao_corr": ("window", r), "n_days": ("window", n)}, coords={"window": labels},
                          attrs=dict(init_date=init_date, definition="corr over days within each lead window "
                                                                     "(single initialization); day 1 = init + 1"))


def nao_tables(daily: dict[str, xr.Dataset], projector: NaoProjector, init_date: str, ref: str = "ERA5"):
    """Window means of every source and window correlations vs ``ref`` as DataFrames."""
    means = pd.DataFrame({k: projector.window_means(v.nao_std, init_date).to_series() for k, v in daily.items()})
    corr = pd.DataFrame({k: projector.window_correlation(v.nao_std, daily[ref].nao_std, init_date)
                         .nao_corr.to_series() for k, v in daily.items() if k != ref})
    return means, corr


def plot_nao_daily(daily: dict[str, xr.Dataset], init_date: str, max_lead_days: int, *, ref: str = "ERA5",
                   style_map=None, corr: pd.DataFrame | None = None) -> plt.Figure:
    from configs.hind_experiment_config import STYLE_MAP
    style_map = style_map or STYLE_MAP
    t0 = np.datetime64(init_date)
    fig, axes = plt.subplots(1, 2 if corr is not None else 1, figsize=(16, 4.8), squeeze=False,
                             gridspec_kw=dict(width_ratios=(2, 1)) if corr is not None else None)
    ax = axes[0, 0]
    for k, ds in daily.items():
        s = ds.nao_std.sel(time=slice(t0, t0 + np.timedelta64(max_lead_days, "D")))
        ax.plot(s.time, s, label=k, **style_map.get(k, {}))
    ax.axhline(0, color="k", linewidth=0.5)
    ax.set_title(f"(a) Daily NAO index from {init_date}", loc="left")
    ax.set_ylabel("NAO (std)")
    ax.grid(True, alpha=0.3)
    ax.legend(ncol=5, fontsize=8, frameon=False)
    if corr is not None:
        ax2 = axes[0, 1]
        x = np.arange(len(corr.index))
        width = 0.8 / max(len(corr.columns), 1)
        for j, k in enumerate(corr.columns):
            ax2.bar(x - 0.4 + width * (j + 0.5), corr[k], width, label=k,
                    color=style_map.get(k, {}).get("color"))
        ax2.set_xticks(x)
        ax2.set_xticklabels(corr.index, rotation=30, ha="right")
        ax2.axhline(0, color="k", linewidth=0.5)
        ax2.set_ylim(-1, 1)
        ax2.set_title(f"(b) Correlation with {ref} within lead windows", loc="left")
        ax2.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    return fig
