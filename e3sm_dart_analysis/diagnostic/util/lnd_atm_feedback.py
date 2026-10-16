"""Land-atmosphere lag coupling of hindcast ensembles from intra-ensemble departures.

For an experiment, a land-side variable x (e.g. top-5 cm soil moisture) and an
atmosphere-side variable y (e.g. LHFLX, PRECT), each land point gets the lagged
correlation

    r(lag) = corr( x'(m, t), y'(m, t + lag) )   pooled over members m and days t,

where x' and y' are departures of each member from the ensemble mean of that
day. Removing the ensemble mean removes what all members share (the seasonal
drift over the 60-day hindcast, the initial shock, the forced response), so
r measures how member-to-member land differences go with atmosphere
differences. Positive lag: the land variable leads; an asymmetric curve shows
which side drives the other (e.g. SM -> PRECT at lag > 0 is the soil-moisture
feedback, lag < 0 is precipitation wetting the soil).

Daily values come from ``util.model_daily``: ``ts/daily_mean`` 00Z-24Z means
for the atmosphere, the ELM ``ts/daily`` series for land; only daily means are
used. Land points: land fraction >= 0.5 north of 60S.

Significance: two-sided t-test of r with the pooled sample size reduced for
the lag-1 autocorrelation of x' and y' (n_eff = n (1 - a_x a_y) / (1 + a_x a_y)).
Regional curves are cos-lat weighted means of r over the region's land points,
with a confidence interval from a circular block bootstrap over days (all
members kept, so the departures are unchanged; blocks keep the day-to-day
persistence). The interval can be asymmetric when a few high-spread days
carry much of the coupling.
"""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from scipy import stats

from configs.experiment_labels import experiment_label
from util.bias_metrics import _period_years
from util.forecast_ensemble_ts import region_masks
from util.forecast_regional_timeseries import REGIONS as BOX_REGIONS
from util.lnd_sm_daily_skill import top_layer_weights
from util.model_daily import open_model_daily
from util.model_evaluation import experiment_colors

CACHE_VERSION = 1  # bump when the computation changes
REGIONS = ("land",) + tuple(r for r in BOX_REGIONS if r != "global")
# name -> (realm, model variable, units, scale to units)
VARIABLES = {
    "SOILM_5CM": ("lnd", "H2OSOI", "m3 m-3", 1.0),
    "SOILWATER_10CM": ("lnd", "SOILWATER_10CM", "kg m-2", 1.0),
    "TSA": ("lnd", "TSA", "K", 1.0),
    "LHFLX": ("atm", "LHFLX", "W m-2", 1.0),
    "SHFLX": ("atm", "SHFLX", "W m-2", 1.0),
    "PRECT": ("atm", "PRECT", "mm day-1", 86400e3),
    "TREFHT": ("atm", "TREFHT", "K", 1.0),
    "TS": ("atm", "TS", "K", 1.0),
}
_SOUTH_LIMIT = -60.0  # Antarctica left out


def hindcast_days(period, skip_days=0):
    """Days of a ``YYYYMM-YYYYMM`` hindcast period, without the first ``skip_days``."""
    y0, m0, y1, m1, _ = _period_years(period)
    days = pd.date_range(pd.Timestamp(y0, m0, 1), pd.Timestamp(y1, m1, 1) + pd.offsets.MonthEnd(0), freq="D")
    return days[skip_days:]


def feedback_cache_path(cache_dir, experiment, x_variable, y_variable, period, skip_days):
    return Path(cache_dir) / f"{experiment}_{x_variable}_vs_{y_variable}_{period}_skip{skip_days}d_feedback.nc"


def _land_points(lat, lon, landmask_file):
    land = region_masks(lat, lon, ("land",), landmask_file).isel(region=0)
    return land & (land.lat >= _SOUTH_LIMIT)


def _member_field(config, member, variable, days, dz):
    """One member's daily ``variable`` on ``days`` (time, lat, lon), in VARIABLES units."""
    realm, name, _, scale = VARIABLES[variable]
    ts_root = Path(config["path"]) / member / "archive/post" / realm / "180x360_aave/ts"
    field = open_model_daily(ts_root, name, member, sorted({day.year for day in days}))
    if field.attrs.get("sampling") != "daily_mean":
        raise ValueError(f"{variable}: {member} has only 00Z states, not daily means")
    field = field.reindex(time=days)
    if variable == "SOILM_5CM":
        weights = top_layer_weights(dz)
        weights = weights.isel(levgrnd=slice(0, int(weights.notnull().any(("lat", "lon")).sum())))
        field = field.isel(levgrnd=slice(0, weights.sizes["levgrnd"]))
        field = (field * weights.assign_coords(levgrnd=field.levgrnd)).sum("levgrnd", min_count=1)
    return (field * scale).astype("float32")


def _read_ensemble(config, variable, days, land, dz, max_members, log):
    """(member, time, point) float32 values at the land points."""
    count = config["nens"] if max_members is None else min(max_members, config["nens"])
    index = np.flatnonzero(land.values.ravel())
    out = None
    for number in range(1, count + 1):
        field = _member_field(config, f"EN{number:02d}", variable, days, dz).load()
        values = field.values.reshape(len(days), -1)[:, index]
        if out is None:
            out = np.empty((count, len(days), len(index)), dtype="float32")
        out[number - 1] = values
    log(f"  {variable}: {count} members, {len(days)} days, {len(index)} land points")
    return out


def _lag_pairs(x, y, lag):
    """x(t) and y(t + lag) over the days both exist; (member, time, point)."""
    if lag >= 0:
        return x[:, :x.shape[1] - lag], y[:, lag:]
    return x[:, -lag:], y[:, :y.shape[1] + lag]


def _day_sums(a, b):
    """Sums over members per day (stat, time, point) for pooled correlations; missing pairs skipped."""
    valid = np.isfinite(a) & np.isfinite(b)
    a = np.where(valid, a, 0.0).astype("float64")
    b = np.where(valid, b, 0.0).astype("float64")
    return np.stack([valid.sum(0), a.sum(0), b.sum(0), (a * a).sum(0), (b * b).sum(0), (a * b).sum(0)])


def _block_weights(n_days, n_boot, block_days, rng):
    """(boot, day) counts of a circular block bootstrap of ``n_days`` days (every day equally likely)."""
    block = max(1, min(block_days, n_days))
    starts = rng.integers(0, n_days, size=(n_boot, int(np.ceil(n_days / block))))
    picks = ((starts[..., None] + np.arange(block)) % n_days).reshape(n_boot, -1)[:, :n_days]
    weights = np.zeros((n_boot, n_days))
    np.add.at(weights, (np.arange(n_boot)[:, None], picks), 1.0)
    return weights


def _pooled_r(sums):
    n, sa, sb, saa, sbb, sab = sums
    with np.errstate(invalid="ignore", divide="ignore"):
        cov = sab - sa * sb / n
        var = (saa - sa * sa / n) * (sbb - sb * sb / n)
        return np.where((n > 2) & (var > 0), cov / np.sqrt(var), np.nan)


def compute_feedback(config, x_variable, y_variable, days, lags, *, landmask_file, land_support_dir,
                     regions=REGIONS, alpha=0.05, n_boot=500, block_days=5, seed=0, max_members=None, log=print):
    """Lagged correlation maps, p-values and regional curves for one experiment (see module docstring)."""
    with xr.open_dataset(Path(land_support_dir) / "dzsoi_elm.nc") as ds:
        dz = ds["DZSOI"].load()
    probe = _member_field(config, "EN01", y_variable, days[:1], dz)
    lat, lon = probe.lat, probe.lon
    land = _land_points(lat, lon, landmask_file)
    x = _read_ensemble(config, x_variable, days, land, dz, max_members, log)
    y = _read_ensemble(config, y_variable, days, land, dz, max_members, log)
    if x.shape[0] < 2:
        raise ValueError("Departures from the ensemble mean need at least 2 members")
    x -= np.nanmean(x, axis=0, keepdims=True)
    y -= np.nanmean(y, axis=0, keepdims=True)
    autocorr = [_pooled_r(_day_sums(*_lag_pairs(v, v, 1)).sum(1)) for v in (x, y)]
    persistence = np.clip(autocorr[0] * autocorr[1], -0.99, 0.99)

    point_lat = np.broadcast_to(lat.values[:, None], land.shape)[land.values]
    weight = np.cos(np.deg2rad(point_lat))
    masks = region_masks(lat, lon, regions, landmask_file)
    in_region = masks.values.reshape(len(regions), -1)[:, np.flatnonzero(land.values.ravel())]
    rng = np.random.default_rng(seed)
    members = x.shape[0]

    r_points, n_points, curves, low, high = [], [], [], [], []
    for lag in lags:
        sums = _day_sums(*_lag_pairs(x, y, lag))  # (stat, day, point)
        r_points.append(_pooled_r(sums.sum(1)))
        n_points.append(sums[0].sum(0))
        resample = _block_weights(sums.shape[1], n_boot, block_days, rng)
        boot = _pooled_r(np.einsum("bt,stp->sbp", resample, sums))  # (boot, point)
        row, row_low, row_high = [], [], []
        for mask in in_region:
            w = np.where(mask, weight, 0.0)
            full = np.nansum(r_points[-1] * w) / np.sum(w * np.isfinite(r_points[-1]))
            sample = np.nansum(boot * w, axis=1) / np.sum(w * np.isfinite(boot), axis=1)
            row.append(full)
            row_low.append(np.nanpercentile(sample, 100 * alpha / 2))
            row_high.append(np.nanpercentile(sample, 100 * (1 - alpha / 2)))
        curves.append(row)
        low.append(row_low)
        high.append(row_high)

    r_points, n_points = np.array(r_points), np.array(n_points)
    n_eff = np.clip(n_points * (1 - persistence) / (1 + persistence), 3, None)
    with np.errstate(invalid="ignore", divide="ignore"):
        t = r_points * np.sqrt((n_eff - 2) / (1 - r_points ** 2))
    p_points = 2 * stats.t.sf(np.abs(t), n_eff - 2)

    def to_grid(values):
        grid = np.full((len(lags), land.size), np.nan, dtype="float32")
        grid[:, np.flatnonzero(land.values.ravel())] = values
        return xr.DataArray(grid.reshape(len(lags), *land.shape), dims=("lag", "lat", "lon"),
                            coords={"lag": list(lags), "lat": lat.values, "lon": lon.values})

    curve_coords = {"lag": list(lags), "region": list(regions)}
    ds = xr.Dataset({
        "r": to_grid(r_points), "p_value": to_grid(p_points), "n_eff": to_grid(n_eff),
        "curve": xr.DataArray(np.array(curves), dims=("lag", "region"), coords=curve_coords),
        "curve_low": xr.DataArray(np.array(low), dims=("lag", "region"), coords=curve_coords),
        "curve_high": xr.DataArray(np.array(high), dims=("lag", "region"), coords=curve_coords),
    })
    # No "units" attribute: "days" would make xarray decode lag as a timedelta on reload.
    ds["lag"].attrs.update(long_name=f"lag in days ({x_variable} leads for lag > 0)")
    ds["r"].attrs.update(long_name=f"correlation of {x_variable}(t) and {y_variable}(t + lag) ensemble departures")
    ds["curve"].attrs.update(long_name="cos-lat weighted land mean of r")
    ds.attrs.update(cache_version=CACHE_VERSION, x_variable=x_variable, y_variable=y_variable,
                    x_units=VARIABLES[x_variable][2], y_units=VARIABLES[y_variable][2],
                    members=members, first_day=str(days[0].date()), last_day=str(days[-1].date()),
                    alpha=alpha, n_boot=n_boot, block_days=block_days,
                    anomalies="departures of each member from the ensemble mean of the day",
                    significance="two-sided t-test, n_eff = n (1 - a_x a_y) / (1 + a_x a_y), lag-1 autocorrelations",
                    curve_interval="circular block bootstrap over days, percentiles")
    return ds


def cache_feedback(path, *, force_compute=False, lags=None, **kwargs):
    """Reuse ``path`` when it holds this cache version and lags, else compute and write it."""
    path = Path(path)
    if path.is_file() and not force_compute:
        ds = xr.load_dataset(path)
        if ds.attrs.get("cache_version") == CACHE_VERSION and (lags is None or list(ds.lag.values) == list(lags)):
            return ds, False
    ds = compute_feedback(lags=lags, **kwargs)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    ds.to_netcdf(tmp)
    os.replace(tmp, path)
    return ds, True


def plot_lag_curves(results, regions, *, title=""):
    """One panel per region: r(lag) per experiment with its block-bootstrap interval."""
    colors = experiment_colors(results)
    ncol = min(3, len(regions))
    nrow = int(np.ceil(len(regions) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.2 * nrow), sharex=True, squeeze=False)
    for ax, region in zip(axes.flat, regions):
        for experiment, ds in results.items():
            curve = ds.sel(region=region)
            ax.fill_between(curve.lag, curve.curve_low, curve.curve_high, color=colors[experiment], alpha=0.18, lw=0)
            ax.plot(curve.lag, curve.curve, color=colors[experiment], lw=1.8, label=experiment_label(experiment))
        ax.axhline(0, color="gray", lw=0.7)
        ax.axvline(0, color="gray", lw=0.7, ls="--")
        ax.set_title(region)
        ax.grid(alpha=0.3)
    for ax in axes.flat[len(regions):]:
        ax.set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel("Lag (days; land leads for lag > 0)")
    for ax in axes[:, 0]:
        ax.set_ylabel("r of ensemble departures")
    axes.flat[0].legend(fontsize=8, frameon=False)
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_lag_maps(results, lag, *, alpha=0.05, title="", vmax=0.6):
    """Map of r at ``lag`` per experiment; dots where p < alpha."""
    import cartopy.crs as ccrs

    names = list(results)
    ncol = min(2, len(names))
    nrow = int(np.ceil(len(names) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.5 * ncol, 3.4 * nrow), squeeze=False,
                             subplot_kw={"projection": ccrs.PlateCarree()})
    levels = np.linspace(-vmax, vmax, 13)
    for ax, experiment in zip(axes.flat, names):
        ds = results[experiment].sel(lag=lag)
        ds = ds.assign_coords(lon=((ds.lon + 180) % 360) - 180).sortby("lon")
        r, p = ds.r, ds.p_value
        mesh = ax.contourf(r.lon, r.lat, r, levels=levels, cmap="RdBu_r", extend="both", transform=ccrs.PlateCarree())
        lon2d, lat2d = np.meshgrid(p.lon, p.lat)
        dots = (p.values < alpha)
        dots[::2, :] = False  # thin the stipple
        dots[:, ::2] = False
        ax.scatter(lon2d[dots], lat2d[dots], s=0.3, c="k", transform=ccrs.PlateCarree())
        ax.coastlines(lw=0.5)
        ax.set_extent([-180, 180, _SOUTH_LIMIT, 90], crs=ccrs.PlateCarree())
        ax.set_title(f"{experiment_label(experiment)} (lag {lag} d)")
    for ax in axes.flat[len(names):]:
        ax.set_visible(False)
    fig.colorbar(mesh, ax=axes, orientation="horizontal", shrink=0.6, pad=0.06,
                 label=f"r of ensemble departures; dots: p < {alpha}")
    fig.suptitle(title, y=0.98 + 0.04 / nrow)
    return fig
