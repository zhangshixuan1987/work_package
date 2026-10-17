"""MJO Hovmöller diagrams and convection-wind coupling of the AMIP hindcasts.

Hovmöller: daily anomalies (seasonal cycle removed, see below) averaged over
``lat_band`` (cos-lat weighted) for each longitude in ``lon_sector`` and
standardized by their own std over the period; written per source
(``mjo_hov_<SOURCE>_<var>_<period>.nc``), observations as ``OBS``.

Coupling: area means over ``lat_band`` x ``lon_sector`` of OLR and U850
anomalies, band-passed (FFT, ``bandpass_days``) and standardized, then the
lag correlation corr(OLR(t), U850(t + lag)); positive lag = U850 lags OLR
(``mjo_coupling_OLR_U850_<SOURCE>_<period>.nc``). The convective MJO phase
(low OLR) leads westerly U850, giving a negative minimum at small positive lag.

Seasonal cycle: ``"harmonics"`` (default) removes the period mean plus the
first ``n_harmonics`` annual harmonics, a smooth estimate from a two-year
record; ``"daily_climatology"`` removes the day-of-year mean of the period
(noisy with two years, kept to reproduce earlier figures); ``"none"``.
Model fields keep their own sign: model and NOAA OLR are both outgoing
(positive up), so no sign alignment with the reference is applied.
"""
from __future__ import annotations

import string
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from configs.hind_experiment_config import (
    DEFAULT_INPUT_DATA_ROOT,
    END_DATE,
    PERIOD,
    START_DATE,
    STYLE_MAP,
    check_experiments,
)
from util.hind_model_data import open_model_daily
from util.hind_observation_data import (
    DEFAULT_OBSERVATION_ROOT,
    open_reference_daily,
)

# variable -> (reference product, unit scale applied to model and reference)
MJO_FIELDS = {
    "FLUT": dict(product="NOAA_OLR", scale=1.0, units="W m-2", label="OLR"),
    "PRECT": dict(product="IMERG", scale=86400.0 * 1000.0, units="mm day-1", label="PRECT"),
    "U850": dict(product="ERA5_DAILY", scale=1.0, units="m s-1", label="U850"),
    "U200": dict(product="ERA5_DAILY", scale=1.0, units="m s-1", label="U200"),
}


@dataclass(frozen=True)
class MjoConfig:
    lat_band: tuple[float, float] = (-10.0, 10.0)
    hov_lon_sector: tuple[float, float] = (0.0, 360.0)
    coupling_lon_sector: tuple[float, float] = (60.0, 150.0)
    bandpass_days: tuple[float, float] = (20.0, 100.0)
    lags: np.ndarray = field(default_factory=lambda: np.arange(-20, 21))
    seasonal_cycle: str = "harmonics"   # "harmonics" | "daily_climatology" | "none"
    n_harmonics: int = 3


def remove_seasonal_cycle(da: xr.DataArray, method: str = "harmonics", n_harmonics: int = 3) -> xr.DataArray:
    if method == "none":
        return da
    if method == "daily_climatology":
        g = da.groupby("time.dayofyear")
        return (g - g.mean("time")).drop_vars("dayofyear", errors="ignore")
    if method != "harmonics":
        raise ValueError(f"Unknown seasonal-cycle method {method!r}")
    doy = da.time.dt.dayofyear.values.astype(float)
    phase = 2 * np.pi * doy / 365.25
    basis = [np.ones_like(phase)]
    for k in range(1, n_harmonics + 1):
        basis += [np.cos(k * phase), np.sin(k * phase)]
    A = np.stack(basis, axis=1)                          # (time, nbasis)
    x = da.transpose("time", ...).values
    shape = x.shape
    coef, *_ = np.linalg.lstsq(A, x.reshape(shape[0], -1), rcond=None)
    fit = (A @ coef).reshape(shape)
    return da.transpose("time", ...).copy(data=x - fit).transpose(*da.dims)


def _band(da: xr.DataArray, lat_band, lon_sector) -> xr.DataArray:
    sub = da.sel(lat=slice(*lat_band), lon=slice(*lon_sector))
    if sub.sizes["lat"] == 0 or sub.sizes["lon"] == 0:
        raise ValueError(f"Empty selection lat={lat_band} lon={lon_sector}")
    return sub.weighted(np.cos(np.deg2rad(sub.lat))).mean("lat")


def standardize(da: xr.DataArray, dim: str = "time") -> xr.DataArray:
    sd = da.std(dim)
    return ((da - da.mean(dim)) / sd.where(sd > 0)).fillna(0.0)


def bandpass_fft(da: xr.DataArray, low_days: float, high_days: float) -> xr.DataArray:
    """FFT band-pass along ``time`` of an evenly spaced daily series (gaps filled linearly first)."""
    x = da.transpose("time", ...).values
    n = x.shape[0]
    freqs = np.fft.rfftfreq(n, d=1.0)
    mask = (freqs >= 1.0 / high_days) & (freqs <= 1.0 / low_days)
    X = np.fft.rfft(x, axis=0) * mask.reshape((-1,) + (1,) * (x.ndim - 1))
    return da.transpose("time", ...).copy(data=np.fft.irfft(X, n=n, axis=0)).transpose(*da.dims)


def _fill_daily(da: xr.DataArray) -> xr.DataArray:
    """Reindex to every calendar day and fill gaps (e.g. noleap Feb 29) linearly."""
    days = np.arange(da.time.values[0].astype("datetime64[D]"), da.time.values[-1].astype("datetime64[D]") + 1)
    return da.reindex(time=days.astype("datetime64[ns]")).interpolate_na("time")


class MjoDiagnostics:
    def __init__(self, cfg: MjoConfig = MjoConfig()):
        self.cfg = cfg

    def anomalies(self, da: xr.DataArray) -> xr.DataArray:
        return remove_seasonal_cycle(_fill_daily(da), self.cfg.seasonal_cycle, self.cfg.n_harmonics)

    # ``*_reduce``: lazy area means (the expensive part, batched with ``compute_all``);
    # ``*_finish``: anomalies and filtering of the in-memory result.
    def hov_reduce(self, da: xr.DataArray) -> xr.DataArray:
        return _band(da, self.cfg.lat_band, self.cfg.hov_lon_sector)

    def hov_finish(self, band: xr.DataArray) -> xr.DataArray:
        return standardize(self.anomalies(band.compute()))

    def hovmoller(self, da: xr.DataArray) -> xr.DataArray:
        """Standardized lat-band-mean anomalies (time, lon)."""
        return self.hov_finish(self.hov_reduce(da))

    def sector_reduce(self, da: xr.DataArray) -> xr.DataArray:
        return _band(da, self.cfg.lat_band, self.cfg.coupling_lon_sector).mean("lon")

    def sector_finish(self, ts: xr.DataArray) -> xr.DataArray:
        return standardize(bandpass_fft(self.anomalies(ts.compute()), *self.cfg.bandpass_days))

    def sector_series(self, da: xr.DataArray) -> xr.DataArray:
        return self.sector_finish(self.sector_reduce(da))

    def lag_correlation(self, a: xr.DataArray, b: xr.DataArray) -> xr.DataArray:
        """corr(a(t), b(t + lag)) for every lag."""
        a, b = xr.align(a, b, join="inner")
        av, bv = a.values, b.values
        out = []
        for lag in self.cfg.lags:
            x, y = (av[:len(av) - lag], bv[lag:]) if lag >= 0 else (av[-lag:], bv[:lag])
            out.append(np.corrcoef(x, y)[0, 1])
        return xr.DataArray(out, dims="lag", coords={"lag": self.cfg.lags}, name="lagcorr")


def hov_path(data_dir, source, variable, period=PERIOD) -> Path:
    return Path(data_dir) / f"mjo_hov_{source}_{variable}_{period}.nc"


def coupling_path(data_dir, source, period=PERIOD) -> Path:
    return Path(data_dir) / f"mjo_coupling_OLR_U850_{source}_{period}.nc"


def open_mjo_input(
    source: str,
    variable: str,
    *,
    data_root=DEFAULT_INPUT_DATA_ROOT,
    observation_root=DEFAULT_OBSERVATION_ROOT,
    start: str = START_DATE,
    end: str = END_DATE,
    lat_range=None,
) -> xr.DataArray:
    """Lazy daily ``variable`` of an experiment or ``OBS`` (product from ``MJO_FIELDS``), in common units."""
    meta = MJO_FIELDS[variable]
    if source == "OBS":
        da = open_reference_daily(meta["product"], variable, start=start, end=end,
                                  observation_root=observation_root, lat_range=lat_range)
    else:
        check_experiments([source])
        da = open_model_daily(source, variable, data_root=data_root, start=start, end=end)
    return (da * meta["scale"]).assign_attrs(units=meta["units"])


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def _skill(model: xr.DataArray, obs: xr.DataArray) -> tuple[float, float]:
    m, o = xr.align(model, obs, join="inner")
    a, b = m.values.ravel(), o.values.ravel()
    ok = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[ok], b[ok])[0, 1]), float(np.sqrt(np.mean((a[ok] - b[ok]) ** 2)))


def plot_hovmoller_grid(data_dir, variable: str, experiments, *, period=PERIOD, days=(0, 120), ncols=3,
                        robust_q=(2.0, 98.0), cmap="RdBu_r", figsize=None, font_scale=11) -> plt.Figure:
    """OBS and every experiment for ``days`` since the start; corr/RMSE vs OBS in titles."""
    def load(src):
        da = xr.load_dataset(hov_path(data_dir, src, variable, period)).hov
        return da.isel(time=slice(*days))
    obs = load("OBS")
    panels = [("OBS", obs)] + [(e, load(e)) for e in experiments]
    vals = np.concatenate([p.values.ravel() for _, p in panels])
    vmax = float(np.nanmax(np.abs(np.nanpercentile(vals, robust_q))))
    nrows = int(np.ceil(len(panels) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize or (4.6 * ncols, 5.2 * nrows), squeeze=False,
                             sharey=True)
    t0 = obs.time.values[0]
    for i, (name, da) in enumerate(panels):
        ax = axes.flat[i]
        y = (da.time.values - t0) / np.timedelta64(1, "D")
        im = ax.pcolormesh(da.lon, y, da.transpose("time", "lon"), cmap=cmap, vmin=-vmax, vmax=vmax,
                           shading="auto")
        title = name if name == "OBS" else "{} (r={:.2f}, rmse={:.2f})".format(name, *_skill(da, obs))
        ax.set_title(f"({string.ascii_lowercase[i]}) {title}", loc="left", fontsize=font_scale)
        ax.set_xlabel("Longitude (°E)")
        if i % ncols == 0:
            ax.set_ylabel(f"Days since {np.datetime_as_string(t0, unit='D')}")
    for j in range(len(panels), nrows * ncols):
        axes.flat[j].axis("off")
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.8, label="Standardized anomaly")
    fig.suptitle(f"MJO Hovmöller: {MJO_FIELDS[variable]['label']}", fontsize=font_scale * 1.2)
    return fig


def load_coupling(data_dir, sources, period=PERIOD) -> dict[str, xr.DataArray]:
    return {s: xr.load_dataset(coupling_path(data_dir, s, period)).lagcorr for s in sources}


def coupling_summary(curves: dict[str, xr.DataArray]):
    import pandas as pd
    rows = []
    for name, c in curves.items():
        i = int(np.nanargmin(c.values))
        rows.append(dict(source=name, min_corr=float(c.values[i]), lag_at_min=int(c.lag.values[i]),
                         corr_lag0=float(c.sel(lag=0))))
    return pd.DataFrame(rows).set_index("source")


def plot_coupling(curves: dict[str, xr.DataArray], *, style_map=None, ylim=None, ref="OBS",
                  title="OLR–U850 lag correlation (20–100 d, 10°S–10°N, 60–150°E)") -> plt.Figure:
    """Lag-correlation curves and their minimum and its lag; ``ylim=None`` fits every curve."""
    style_map = style_map or STYLE_MAP
    summary = coupling_summary(curves)
    fig, (ax, ax2, ax3) = plt.subplots(1, 3, figsize=(17, 5.4), gridspec_kw=dict(width_ratios=(2, 1, 1)))
    for name, c in curves.items():
        ax.plot(c.lag, c, label=name, **style_map.get(name, {}))
    ax.axhline(0, color="k", linewidth=0.6)
    ax.axvline(0, color="k", linewidth=0.6)
    ax.set_xlabel("Lag (days): positive = U850 lags OLR")
    ax.set_ylabel("Correlation")
    if ylim is None:
        vals = np.concatenate([c.values for c in curves.values()])
        lo, hi = float(np.nanmin(vals)), float(np.nanmax(vals))
        ylim = (lo - 0.05 * (hi - lo), hi + 0.05 * (hi - lo))
    ax.set_ylim(*ylim)
    ax.set_title(f"(a) {title}", loc="left")
    ax.legend(ncol=5, fontsize=8, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.16))
    x = np.arange(len(summary))
    colors = [style_map.get(n, {}).get("color", "0.5") for n in summary.index]
    ax2.bar(x, summary.min_corr, color=colors)
    ax2.set_title("(b) Minimum correlation", loc="left")
    ax3.bar(x, summary.lag_at_min, color=colors)
    ax3.set_title("(c) Lag at minimum (days)", loc="left")
    if ref in summary.index:
        for a, col in ((ax2, "min_corr"), (ax3, "lag_at_min")):
            a.axhline(summary.loc[ref, col], color="k", linestyle="--", linewidth=0.8, label=ref)
        ax3.legend(fontsize=8, frameon=False)
    for a in (ax, ax2, ax3):
        a.grid(True, alpha=0.3)
    for a in (ax2, ax3):
        a.set_xticks(x)
        a.set_xticklabels(summary.index, rotation=45, ha="right")
    fig.tight_layout()
    return fig
