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
4. ``window_means``: consecutive lead windows from the initial date
   (day 1 = init + 1).
5. ``slp_acc``: daily cos-lat weighted, centred spatial anomaly correlation of
   the SLP anomalies over the NAO domain with the reference (lead 0 = initial
   day), with ``acc_window_table`` and the e-folding lead. A spatial ACC has many
   more degrees of freedom than a correlation of the one-dimensional index
   across the ~14 strongly autocorrelated days of a window (effective sample
   size ~2-3), which is no longer computed.
"""
from __future__ import annotations

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from util.hind_factorial import plot_group_curves
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
        """Daily ``nao_raw`` (Pa), ``nao_std`` and the SLP anomaly field ``slp_anom`` (Pa, float32,
        NAO domain on the base grid) from daily SLP (time on datetime64 days).

        ``compute=False`` returns the lazy index, e.g. to compute many sources in one ``compute_all``.
        """
        da = _domain(psl_daily, self.cfg).interp(lat=self.base.lat, lon=self.base.lon)
        anom = da.groupby("time.month") - self.base.base_climatology
        raw = (anom * _weights(anom.lat) * self.base.nao_eof1.fillna(0)).sum(("lat", "lon"))
        raw = raw.drop_vars("month", errors="ignore")
        anom = anom.drop_vars("month", errors="ignore").astype("float32")
        if compute:
            raw, anom = raw.compute(), anom.compute()
        return xr.Dataset({"nao_raw": raw.assign_attrs(units="Pa"),
                           "nao_std": ((raw - self.ref_mean) / self.ref_std).assign_attrs(units="1"),
                           "slp_anom": anom.assign_attrs(units="Pa", long_name="SLP minus base monthly climatology")},
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


def slp_acc(fcst: xr.Dataset, ref: xr.Dataset, init_date: str, max_lead_days: int) -> xr.DataArray:
    """Daily centred, cos-lat weighted ACC of ``slp_anom`` with ``ref`` over the NAO domain (lead 0 = init)."""
    t0 = np.datetime64(init_date, "D")
    sl = slice(t0, t0 + np.timedelta64(max_lead_days, "D"))
    f, r = xr.align(fcst.slp_anom.sel(time=sl), ref.slp_anom.sel(time=sl), join="inner")
    w = np.cos(np.deg2rad(f.lat))
    w = (w / w.mean()).broadcast_like(f.isel(time=0, drop=True)).where(np.isfinite(f.isel(time=0, drop=True)))
    wm = lambda x: (x * w).sum(("lat", "lon")) / w.sum()
    a, b = f - wm(f), r - wm(r)
    acc = wm(a * b) / np.sqrt(wm(a * a) * wm(b * b))
    lead = ((f.time.values.astype("datetime64[D]") - t0) / np.timedelta64(1, "D")).astype(int)
    return acc.assign_coords(lead=("time", lead)).swap_dims(time="lead").drop_vars("time").rename("ACC")


def acc_window_table(acc: dict[str, xr.DataArray], windows) -> pd.DataFrame:
    """Window-mean ACC (columns ``ACC dAA-BB``), day-0 ACC and e-folding lead per source."""
    from util.tropical_wind_skill import efold_lead
    rows = {}
    for k, a in acc.items():
        row = {f"ACC d{w0:02d}-{w1:02d}": float(a.sel(lead=slice(w0, w1)).mean()) for w0, w1 in windows}
        row.update(ACC_day0=float(a.sel(lead=0)), efold_days=efold_lead(a.values))
        rows[k] = row
    return pd.DataFrame(rows).T


def nao_tables(daily: dict[str, xr.Dataset], projector: NaoProjector, init_date: str, ref: str = "ERA5", *,
               acc_windows=((0, 4), (5, 9), (10, 14), (15, 28)), acc_max_lead: int = 42):
    """Window means of the index of every source, daily SLP ACC vs ``ref`` and its window table."""
    means = pd.DataFrame({k: projector.window_means(v.nao_std, init_date).to_series() for k, v in daily.items()})
    acc = {k: slp_acc(v, daily[ref], init_date, acc_max_lead) for k, v in daily.items() if k != ref}
    return means, acc, acc_window_table(acc, acc_windows)


def plot_nao_daily(daily: dict[str, xr.Dataset], init_date: str, max_lead_days: int, *, ref: str = "ERA5",
                   style_map=None, acc: dict[str, xr.DataArray] | None = None, acc_max_lead: int = 30) -> plt.Figure:
    """(a) daily NAO index; (b) daily North-Atlantic SLP ACC vs ``ref`` grouped by atmosphere initial state."""
    from configs.hind_experiment_config import STYLE_MAP
    style_map = style_map or STYLE_MAP
    t0 = np.datetime64(init_date)
    fig, axes = plt.subplots(1, 2 if acc is not None else 1, figsize=(17, 4.8), squeeze=False,
                             gridspec_kw=dict(width_ratios=(1.6, 1)) if acc is not None else None)
    ax = axes[0, 0]
    for k, ds in daily.items():
        s = ds.nao_std.sel(time=slice(t0, t0 + np.timedelta64(max_lead_days, "D")))
        ax.plot(s.time, s, label=k, **style_map.get(k, {}))
    ax.axhline(0, color="k", linewidth=0.5)
    ax.set_title(f"(a) Daily NAO index from {init_date}", loc="left")
    ax.set_ylabel("NAO (std)")
    ax.grid(True, alpha=0.3)
    ax.legend(ncol=5, fontsize=8, frameon=False)
    if acc is not None:
        ax2 = axes[0, 1]
        sub = {k: a.sel(lead=slice(0, acc_max_lead)) for k, a in acc.items()}
        plot_group_curves(ax2, next(iter(sub.values())).lead, {k: a.values for k, a in sub.items()})
        ax2.axhline(0, color="k", linewidth=0.5)
        ax2.set_xlim(0, acc_max_lead)
        ax2.set_xlabel("Lead day")
        ax2.set_ylabel("ACC")
        ax2.set_title(f"(b) North-Atlantic SLP ACC vs {ref}", loc="left")
        ax2.legend(fontsize=8, frameon=False)
        ax2.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig
