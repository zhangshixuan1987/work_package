"""WH04 RMM index from observed EOFs, for continuous runs and ensemble hindcasts.

Every model and observed RMM in the ``hind_atm_mjo`` and ``fcst_atm_mjo``
workflows is a projection onto one observed EOF pair, so all of them are the
same index (Wheeler and Hendon 2004; Gottschalck et al. 2010).

1. ``compute_rmm_eofs``: combined EOFs of daily OLR, U850 and U200, 15S-15N
   (NOAA OLR, ERA5; 1986-2017 by default): Feb 29 dropped; anomalies from the
   365-day climatology (optionally smoothed); low-frequency variability
   removed by subtracting the mean of the previous ``lowfreq_rm_days`` days
   (WH04; ``lowfreq_method="trailing"``); sqrt(cos lat) weighting; each field
   divided by its space-time std. PCs are divided by their std, then scaled so
   mean(RMM1² + RMM2²) = 1. With the BoM RMM index (``read_bom_rmm``) the
   standardized PC pair is rotated in the RMM plane (2x2 orthogonal Procrustes
   fit, ``rmm_rotation``) onto it, so RMM1/RMM2 and phases 1-8 have their usual
   meaning; the leading EOF pair has nearly equal variance, so its orientation
   is otherwise arbitrary.
   The file keeps the climatologies and normalisation scalars.
2. ``verify_rmm_eofs``: re-projection of the training data (must reproduce the
   saved RMMs), EOF orthogonality, and correlation with the BoM index.
3. ``RmmProjector.project``: the same preprocessing for any daily fields. The
   trailing mean needs ``lowfreq_rm_days`` earlier days: observations are read
   with that much lead-in; a run that starts at its initial date gets the
   observed anomalies of the preceding days as ``history`` (Gottschalck et al.
   2010), so days from the start are usable.
4. Skill (Rashid et al. 2011; Lin et al. 2008): bivariate correlation
   COR = Σ(f1 o1 + f2 o2) / sqrt(Σ(f1² + f2²) Σ(o1² + o2²)) and bivariate
   RMSE, without removing means (biases count as errors); amplitude error
   A_f − A_o and phase error θ_f − θ_o (degrees; positive = forecast ahead
   in the eastward propagation). ``rmm_running_skill`` / ``rmm_leadbin_skill``
   pool the days of one continuous run; ``rmm_ensemble_scores`` gives
   per-lead-day ensemble-mean errors and spread of a hindcast ensemble, and
   ``rmm_window_skill`` pools lead windows.
"""
from __future__ import annotations

import string
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from configs.hind_experiment_config import STYLE_MAP

RMM_FIELDS = ("FLUT", "U850", "U200")


@dataclass(frozen=True)
class RmmConfig:
    lat_range: tuple[float, float] = (-15.0, 15.0)
    train_start: str = "1986-01-01"
    train_end: str = "2017-12-31"
    doy_smooth_days: int = 31          # <= 1 disables smoothing of the 365-day climatology
    lowfreq_rm_days: int = 120         # <= 1 disables the low-frequency removal
    lowfreq_method: str = "trailing"   # "trailing" (WH04, previous days) | "centred"
    n_modes: int = 2
    seed: int = 0


# ---------------------------------------------------------------------------
# Observed BoM index
# ---------------------------------------------------------------------------
def read_bom_rmm(path) -> xr.Dataset:
    """BoM ``rmm.74toRealtime.txt`` as daily RMM1, RMM2, phase, amplitude (missing -> NaN)."""
    path = Path(path)
    if path.suffix == ".nc":
        ds = xr.load_dataset(path)
        return ds.assign_coords(time=ds.time.values.astype("datetime64[D]").astype("datetime64[ns]"))
    names = ["year", "month", "day", "RMM1", "RMM2", "phase", "amplitude"]
    df = pd.read_csv(path, sep=r"\s+", header=None, usecols=range(7), names=names, comment=None)
    df = df.apply(pd.to_numeric, errors="coerce").dropna(subset=["year", "month", "day"])
    time = pd.to_datetime(dict(year=df.year.astype(int), month=df.month.astype(int), day=df.day.astype(int)))
    out = df[["RMM1", "RMM2", "phase", "amplitude"]].where(df[["RMM1", "RMM2", "phase", "amplitude"]] < 999)
    return xr.Dataset.from_dataframe(out.set_index(pd.DatetimeIndex(time, name="time")))


# ---------------------------------------------------------------------------
# Shared preprocessing
# ---------------------------------------------------------------------------
def tropics_subset(da: xr.DataArray, lat_range) -> xr.DataArray:
    """Lazy ``lat_range`` subset (time, lat, lon); load it as float32 (e.g. with ``compute_all``) before
    ``compute_rmm_eofs`` / ``RmmProjector.project``, which then skip their own per-field reads."""
    return da.sel(lat=slice(*lat_range)).transpose("time", "lat", "lon")


_tropics = tropics_subset


def drop_feb29(da: xr.DataArray) -> xr.DataArray:
    t = da.time
    return da.sel(time=~((t.dt.month == 2) & (t.dt.day == 29)))


def doy365(da: xr.DataArray) -> xr.DataArray:
    """Day of year 1..365 with Feb 29 removed (leap-year days after Feb 28 shifted by -1)."""
    doy = da.time.dt.dayofyear
    return xr.where(da.time.dt.is_leap_year & (doy > 59), doy - 1, doy).rename("doy365")


def smooth_doy(clim: xr.DataArray, win: int) -> xr.DataArray:
    if win is None or win <= 1:
        return clim
    pad = win // 2
    n = clim.sizes["doy365"]
    padded = xr.concat([clim.isel(doy365=slice(n - pad, n)), clim, clim.isel(doy365=slice(0, pad))], "doy365")
    sm = padded.rolling(doy365=win, center=True, min_periods=win).mean().isel(doy365=slice(pad, pad + n))
    return sm.assign_coords(doy365=clim.doy365.values)


def lat_weight(da: xr.DataArray) -> xr.DataArray:
    return da * np.sqrt(np.cos(np.deg2rad(da.lat)))


def remove_lowfreq(anom: xr.DataArray, window: int, method: str = "trailing") -> xr.DataArray:
    """Subtract the mean of the previous ``window`` days (``trailing``, WH04) or a centred mean.

    Days without a full window are NaN. The series is assumed daily and
    contiguous apart from the dropped Feb 29.
    """
    if window is None or window <= 1:
        return anom
    if method == "trailing":
        low = anom.rolling(time=window, min_periods=window).mean().shift(time=1)
    elif method == "centred":
        low = anom.rolling(time=window, center=True, min_periods=window).mean()
    else:
        raise ValueError(f"Unknown low-frequency method {method!r}")
    return anom - low


def _anomaly(da: xr.DataArray, clim: xr.DataArray) -> xr.DataArray:
    da = drop_feb29(da)
    return (da.groupby(doy365(da)) - clim).drop_vars("doy365", errors="ignore")


def _stack(fields: list[xr.DataArray]) -> np.ndarray:
    return np.concatenate([f.stack(feature=("lat", "lon")).values for f in fields], axis=1)


def _eof_vectors(ds: xr.Dataset) -> np.ndarray:
    return np.stack([np.concatenate([ds[f"EOF_{v}"].sel(mode=m).values.ravel() for v in RMM_FIELDS])
                     for m in (1, 2)], axis=1)


def _phase(rmm1, rmm2):
    """WH04 phase 1-8 (phase 1: angle 180-225 degrees, counter-clockwise)."""
    return (np.floor(((np.degrees(np.arctan2(rmm2, rmm1)) + 180.0) % 360.0) / 45.0) + 1).astype(int)


def _finish(ds: xr.Dataset) -> xr.Dataset:
    ds["RMM_amp"] = np.sqrt(ds.RMM1 ** 2 + ds.RMM2 ** 2)
    ds["RMM_phase"] = _phase(ds.RMM1, ds.RMM2)
    return ds


# ---------------------------------------------------------------------------
# EOF build and verification
# ---------------------------------------------------------------------------
def _randomized_svd(X: np.ndarray, k: int, oversample: int = 8, niter: int = 2, seed: int = 0):
    rng = np.random.default_rng(seed)
    Y = X @ rng.standard_normal((X.shape[1], k + oversample)).astype(X.dtype)
    for _ in range(niter):
        Y = X @ (X.T @ Y)
    Q, _ = np.linalg.qr(Y)
    Ub, S, Vt = np.linalg.svd(Q.T @ X, full_matrices=False)
    return (Q @ Ub)[:, :k], S[:k], Vt[:k]


def preprocess_training(fields: dict[str, xr.DataArray], cfg: RmmConfig):
    """Climatologies, normalised training anomalies, and std scalars."""
    clim, proc, stds = {}, {}, {}
    need = np.datetime64(cfg.train_end) - np.datetime64(cfg.train_start)
    for v, da in fields.items():
        da = drop_feb29(_tropics(da, cfg.lat_range).sel(time=slice(cfg.train_start, cfg.train_end)))
        span = da.time.values[-1] - da.time.values[0] if da.sizes["time"] else np.timedelta64(0, "D")
        if span < 0.95 * need:
            raise ValueError(f"{v}: data span {str(da.time.values[0])[:10] if da.sizes['time'] else '-'}.."
                             f"{str(da.time.values[-1])[:10] if da.sizes['time'] else '-'} does not cover the "
                             f"training period {cfg.train_start}..{cfg.train_end}")
        da = da.astype("float32").compute()
        c = smooth_doy(da.groupby(doy365(da)).mean("time"), cfg.doy_smooth_days)
        a = lat_weight(remove_lowfreq(_anomaly(da, c), cfg.lowfreq_rm_days, cfg.lowfreq_method))
        s = float(a.std())
        clim[v], proc[v], stds[v] = c, a / s, s
    common = np.intersect1d(np.intersect1d(proc["FLUT"].time, proc["U850"].time), proc["U200"].time)
    proc = {v: p.sel(time=common) for v, p in proc.items()}
    ok = np.logical_and.reduce([np.isfinite(p).all(("lat", "lon")).values for p in proc.values()])
    proc = {v: p.isel(time=np.flatnonzero(ok)) for v, p in proc.items()}
    return clim, proc, stds


def _rotation_to_bom(pcs: np.ndarray, time, bom: xr.Dataset):
    """2x2 orthogonal matrix R minimizing |pcs @ R - BoM| (Procrustes) and the resulting correlations.

    ``pcs`` must have unit variance per mode; the leading pair of WH04 EOFs has nearly
    equal variance, so its orientation in the RMM plane is set by this rotation, which
    keeps RMM1 and RMM2 uncorrelated (sign flips and swaps are special cases).
    """
    b = bom.reindex(time=time)
    ok = np.isfinite(b.RMM1.values) & np.isfinite(b.RMM2.values)
    obs = np.stack([b.RMM1.values[ok], b.RMM2.values[ok]], axis=1)
    U, _, Vt = np.linalg.svd(pcs[ok].T @ obs)
    R = U @ Vt
    rot = pcs[ok] @ R
    corr = tuple(float(np.corrcoef(rot[:, k], obs[:, k])[0, 1]) for k in range(2))
    return R, corr, int(ok.sum())


def compute_rmm_eofs(fields: dict[str, xr.DataArray], cfg: RmmConfig = RmmConfig(), source: str = "",
                     bom: xr.Dataset | None = None) -> xr.Dataset:
    """Observed RMM EOFs from daily ``fields`` = {"FLUT", "U850", "U200"} (lazy; 1-degree grid).

    ``bom``: the BoM index (``read_bom_rmm``) that sets the orientation of the RMM pair.
    """
    clim, proc, stds = preprocess_training(fields, cfg)
    X = _stack([proc[v] for v in RMM_FIELDS]).astype("float64")
    _, S, Vt = _randomized_svd(X, cfg.n_modes, seed=cfg.seed)
    pcs = X @ Vt.T  # projection on the EOFs: the same operation as every later projection
    pc_std = pcs.std(axis=0)
    pcs = pcs / pc_std
    R, attrs = np.eye(2), dict(bom_aligned=0)
    if bom is not None:
        R, corr, n = _rotation_to_bom(pcs, proc["FLUT"].time.values, bom)
        pcs = pcs @ R
        attrs = dict(bom_aligned=1, bom_corr_rmm1=corr[0], bom_corr_rmm2=corr[1], bom_n_days=n,
                     bom_rotation_deg=float(np.degrees(np.arctan2(R[1, 0], R[0, 0]))))
    scale = float(np.sqrt(np.mean(pcs[:, 0] ** 2 + pcs[:, 1] ** 2)))
    pcs = pcs / scale
    ref = proc["FLUT"].isel(time=0, drop=True)
    ngrid = ref.size
    data = {f"EOF_{v}": (("mode", "lat", "lon"), Vt[:, i * ngrid:(i + 1) * ngrid].reshape(cfg.n_modes, *ref.shape))
            for i, v in enumerate(RMM_FIELDS)}
    data.update({f"clim_{v}": clim[v].astype("float32") for v in RMM_FIELDS})
    data.update({f"norm_std_{v}": ((), stds[v]) for v in RMM_FIELDS})
    ds = xr.Dataset(
        dict(data, RMM1=("time", pcs[:, 0]), RMM2=("time", pcs[:, 1]),
             pc_std_mode=("mode", pc_std), rmm_scale=((), scale), singular_values=("mode", S),
             rmm_rotation=(("mode", "rmm"), R)),
        coords=dict(mode=np.arange(1, cfg.n_modes + 1), rmm=[1, 2], lat=ref.lat, lon=ref.lon,
                    time=proc["FLUT"].time),
        attrs=dict(description="Observed RMM EOFs (WH04)", source=source, train_start=cfg.train_start,
                   train_end=cfg.train_end, lat_range=f"{cfg.lat_range[0]:g},{cfg.lat_range[1]:g}",
                   doy_smooth_days=cfg.doy_smooth_days, lowfreq_rm_days=cfg.lowfreq_rm_days,
                   lowfreq_method=cfg.lowfreq_method,
                   normalization="sqrt(cos lat) weighting; each field / space-time std; "
                                 "RMM = ((X @ EOF) / pc_std_mode) @ rmm_rotation / rmm_scale",
                   **attrs),
    )
    return _finish(ds)


class RmmProjector:
    def __init__(self, eofs: xr.Dataset):
        self.eofs = eofs
        a = eofs.attrs
        self.cfg = RmmConfig(lat_range=tuple(float(x) for x in str(a["lat_range"]).strip("()").split(",")),
                             train_start=a["train_start"], train_end=a["train_end"],
                             doy_smooth_days=int(a["doy_smooth_days"]), lowfreq_rm_days=int(a["lowfreq_rm_days"]),
                             lowfreq_method=a.get("lowfreq_method", "centred"))
        self.E = _eof_vectors(eofs)
        self.R = eofs["rmm_rotation"].values if "rmm_rotation" in eofs else np.eye(2)

    def anomalies(self, da: xr.DataArray, var: str) -> xr.DataArray:
        """Anomaly from the EOF-file climatology on the EOF grid (Feb 29 dropped)."""
        e = self.eofs
        da = _tropics(da, self.cfg.lat_range).astype("float32").compute()  # same precision as the training
        if not (np.array_equal(da.lat, e.lat) and np.array_equal(da.lon, e.lon)):
            da = da.interp(lat=e.lat, lon=e.lon)
        return _anomaly(da, e[f"clim_{var}"])

    def preprocess(self, da: xr.DataArray, var: str, history: xr.DataArray | None = None) -> xr.DataArray:
        anom = self.anomalies(da, var)
        if history is not None:
            past = self.anomalies(history, var)
            past = past.sel(time=past.time < anom.time[0])
            full = xr.concat([past, anom], "time")
        else:
            full = anom
        full = remove_lowfreq(full, self.cfg.lowfreq_rm_days, self.cfg.lowfreq_method).sel(time=anom.time)
        return lat_weight(full) / float(self.eofs[f"norm_std_{var}"])

    def project(self, fields: dict[str, xr.DataArray], name: str,
                history: dict[str, xr.DataArray] | None = None) -> xr.Dataset:
        """RMM of daily ``fields`` (time labelled by day); ``history``: observed fields of the preceding days."""
        proc = {v: self.preprocess(fields[v], v, None if history is None else history[v]) for v in RMM_FIELDS}
        common = np.intersect1d(np.intersect1d(proc["FLUT"].time, proc["U850"].time), proc["U200"].time)
        proc = {v: p.sel(time=common) for v, p in proc.items()}
        X = _stack([proc[v] for v in RMM_FIELDS])
        ok = np.isfinite(X).all(axis=1)
        pcs = ((X[ok] @ self.E) / self.eofs.pc_std_mode.values) @ self.R / float(self.eofs.rmm_scale)
        ds = xr.Dataset({"RMM1": ("time", pcs[:, 0]), "RMM2": ("time", pcs[:, 1])}, coords={"time": common[ok]},
                        attrs=dict(name=name, lowfreq_method=self.cfg.lowfreq_method,
                                   history="observed anomalies before the first day" if history else "none"))
        return _finish(ds)

    def project_members(self, members: dict[int, dict[str, xr.DataArray]], name: str,
                        history: dict[str, xr.DataArray] | None = None) -> xr.Dataset:
        """RMM (member, time) of an ensemble: ``members`` = {member: {var: field}}."""
        out = [self.project(f, f"{name} member {m}", history).expand_dims(member=[m]) for m, f in members.items()]
        return xr.concat(out, "member", join="inner").assign_attrs(name=name)


def open_forecast_member_fields(entry: dict, years, *, max_members: int | None = None,
                                variables=RMM_FIELDS) -> dict[int, dict[str, xr.DataArray]]:
    """Daily fields of each member of one ``util.bias_metrics.experiment_registry`` entry.

    Members live in ``<run>/ENxx/archive`` (DA runs) or share ``<run>/archive``
    (hindcasts, member in the file name). Daily values follow
    ``util.model_daily.open_model_daily_estimate``: daily means where they exist,
    else the mean of the two 00Z states bracketing the day (U850, U200).
    """
    from util.model_daily import open_model_daily_estimate

    nens = entry["nens"] if max_members is None else min(entry["nens"], max_members)
    out = {}
    for k in range(1, nens + 1):
        member = f"EN{k:02d}"
        root = Path(entry["path"]) / member / "archive"
        if not root.is_dir():
            root = Path(entry["path"]) / "archive"
        ts = root / "post/atm/180x360_aave/ts"
        out[k] = {v: open_model_daily_estimate(ts, v, member, list(years)) for v in variables}
    return out


def verify_rmm_eofs(eofs: xr.Dataset, fields: dict[str, xr.DataArray], bom: xr.Dataset | None = None) -> dict:
    """Re-project the training data; compare with the saved RMMs and (optionally) the BoM index."""
    proj = RmmProjector(eofs)
    c = proj.cfg
    train = {v: f.sel(time=slice(c.train_start, c.train_end)) for v, f in fields.items()}
    re_ = proj.project(train, "verify")
    r, s = xr.align(re_, eofs[["RMM1", "RMM2"]], join="inner")
    E = proj.E
    out = dict(n_days=int(r.sizes["time"]),
               corr_rmm1=float(np.corrcoef(r.RMM1, s.RMM1)[0, 1]),
               corr_rmm2=float(np.corrcoef(r.RMM2, s.RMM2)[0, 1]),
               rmse_rmm1=float(np.sqrt(((r.RMM1 - s.RMM1) ** 2).mean())),
               rmse_rmm2=float(np.sqrt(((r.RMM2 - s.RMM2) ** 2).mean())),
               eof_dot=float(E[:, 0] @ E[:, 1]),
               mean_amp2=float((eofs.RMM1 ** 2 + eofs.RMM2 ** 2).mean()))
    flags = [k for k in ("corr_rmm1", "corr_rmm2") if out[k] < 0.999]
    if abs(out["eof_dot"]) > 1e-6:
        flags.append("eof_dot")
    if bom is not None:
        b, m = xr.align(bom[["RMM1", "RMM2"]].dropna("time"), s, join="inner")
        out.update(bom_n_days=int(b.sizes["time"]),
                   bom_corr_rmm1=float(np.corrcoef(b.RMM1, m.RMM1)[0, 1]),
                   bom_corr_rmm2=float(np.corrcoef(b.RMM2, m.RMM2)[0, 1]),
                   bom_bivariate_cor=bivariate_cor(m.RMM1.values, m.RMM2.values, b.RMM1.values, b.RMM2.values))
        flags += [k for k in ("bom_corr_rmm1", "bom_corr_rmm2") if out[k] < 0.8]
    out["flags"] = flags
    return out


def plot_rmm_eofs(eofs: xr.Dataset) -> plt.Figure:
    fig, axes = plt.subplots(2, 3, figsize=(16, 6), constrained_layout=True)
    for j, v in enumerate(RMM_FIELDS):
        vmax = float(np.abs(eofs[f"EOF_{v}"]).max())
        for i, m in enumerate((1, 2)):
            im = axes[i, j].pcolormesh(eofs.lon, eofs.lat, eofs[f"EOF_{v}"].sel(mode=m), cmap="RdBu_r",
                                       vmin=-vmax, vmax=vmax, shading="auto")
            axes[i, j].set_title(f"({string.ascii_lowercase[3 * i + j]}) EOF{m} {v}", loc="left")
        fig.colorbar(im, ax=axes[:, j], shrink=0.8)
    fig.suptitle(f"Observed RMM EOFs ({eofs.train_start}..{eofs.train_end}; climatology smoothed "
                 f"{eofs.doy_smooth_days} d; previous-{eofs.lowfreq_rm_days}-day mean removed)")
    return fig


def plot_rmm_diagnostics(rmm: xr.Dataset, title: str = "", bom: xr.Dataset | None = None) -> plt.Figure:
    """Phase space, amplitude histogram, RMM1 power spectrum (and RMM1/RMM2 vs BoM if given)."""
    n = 4 if bom is not None else 3
    fig, axes = plt.subplots(1, n, figsize=(5.4 * n, 4.6), constrained_layout=True)
    axes[0].plot(rmm.RMM1, rmm.RMM2, linewidth=0.3)
    th = np.linspace(0, 2 * np.pi, 200)
    axes[0].plot(np.cos(th), np.sin(th), "k", linewidth=0.8)
    axes[0].set_aspect("equal")
    axes[0].set_xlabel("RMM1")
    axes[0].set_ylabel("RMM2")
    axes[0].set_title("(a) Phase space", loc="left")
    axes[1].hist(rmm.RMM_amp, bins=50, color="0.6")
    axes[1].set_title("(b) Amplitude", loc="left")
    x = rmm.RMM1.values - rmm.RMM1.values.mean()
    f = np.fft.rfftfreq(x.size)
    axes[2].plot(1 / f[1:], (np.abs(np.fft.rfft(x)) ** 2 / x.size)[1:])
    axes[2].set_xscale("log")
    axes[2].set_xlim(2, 200)
    axes[2].set_xlabel("Period (days)")
    axes[2].set_title("(c) RMM1 power spectrum", loc="left")
    if bom is not None:
        b, m = xr.align(bom[["RMM1", "RMM2"]].dropna("time"), rmm, join="inner")
        for k, c in (("RMM1", "tab:blue"), ("RMM2", "tab:orange")):
            r = np.corrcoef(b[k], m[k])[0, 1]
            axes[3].scatter(b[k], m[k], s=1, alpha=0.3, color=c, label=f"{k} r={r:.2f}")
        axes[3].plot([-4, 4], [-4, 4], "k", linewidth=0.6)
        axes[3].set_xlabel("BoM")
        axes[3].set_ylabel("This index")
        axes[3].legend(markerscale=8, frameon=False)
        axes[3].set_title("(d) vs BoM RMM", loc="left")
    fig.suptitle(title)
    return fig


# ---------------------------------------------------------------------------
# Skill
# ---------------------------------------------------------------------------
DEFAULT_LEAD_BINS = (("Wk1-2", 1, 14), ("Wk3-4", 15, 28), ("Wk5-6", 29, 42), ("Wk7-8", 43, 56))


def bivariate_cor(f1, f2, o1, o2) -> float:
    """Rashid et al. (2011) bivariate correlation (means not removed)."""
    den = np.sqrt(np.sum(f1 ** 2 + f2 ** 2) * np.sum(o1 ** 2 + o2 ** 2))
    return float(np.sum(f1 * o1 + f2 * o2) / den) if den > 0 else np.nan


def bivariate_rmse(f1, f2, o1, o2) -> float:
    return float(np.sqrt(np.mean((f1 - o1) ** 2 + (f2 - o2) ** 2)))


def phase_error(f1, f2, o1, o2):
    """θ_f − θ_o in degrees wrapped to (-180, 180]; positive = forecast ahead (counter-clockwise)."""
    d = np.degrees(np.arctan2(f2, f1) - np.arctan2(o2, o1))
    return (d + 180.0) % 360.0 - 180.0


def lead_days(time, init: str) -> np.ndarray:
    return ((np.asarray(time).astype("datetime64[D]") - np.datetime64(init, "D")) / np.timedelta64(1, "D")).astype(int)


def _paired(fcst, obs, start, amp_thresh):
    f, o = xr.align(fcst, obs, join="inner")
    if amp_thresh is not None:
        keep = o.RMM_amp.values >= amp_thresh
        f, o = f.isel(time=keep), o.isel(time=keep)
    return f, o, lead_days(f.time.values, start)


def rmm_running_skill(fcst, obs, start, *, window_days=60, max_lead_days=None, amp_thresh=1.0,
                      min_samples=15) -> xr.Dataset:
    """COR/RMSE over the forward window [L, L + window_days) of retained days, per lead L (one run)."""
    f, o, lead = _paired(fcst, obs, start, amp_thresh)
    max_lead = int(max_lead_days if max_lead_days is not None else lead.max())
    cor, rmse, n = (np.full(max_lead + 1, np.nan) for _ in range(3))
    for L in range(max_lead + 1):
        m = (lead >= L) & (lead < L + window_days)
        n[L] = m.sum()
        if m.sum() >= min_samples:
            args = (f.RMM1.values[m], f.RMM2.values[m], o.RMM1.values[m], o.RMM2.values[m])
            cor[L], rmse[L] = bivariate_cor(*args), bivariate_rmse(*args)
    return xr.Dataset({"COR": ("lead_day", cor), "RMSE": ("lead_day", rmse), "n_samples": ("lead_day", n)},
                      coords={"lead_day": np.arange(max_lead + 1)},
                      attrs=dict(window_days=window_days, amp_thresh=np.nan if amp_thresh is None else amp_thresh))


def rmm_leadbin_skill(fcst, obs, start, *, bins=DEFAULT_LEAD_BINS, amp_thresh=1.0, min_samples=5) -> pd.DataFrame:
    """COR/RMSE pooled over the days of each lead bin (one run)."""
    f, o, lead = _paired(fcst, obs, start, amp_thresh)
    rows = []
    for name, a, b in bins:
        m = (lead >= a) & (lead <= b)
        args = (f.RMM1.values[m], f.RMM2.values[m], o.RMM1.values[m], o.RMM2.values[m])
        ok = m.sum() >= min_samples
        rows.append(dict(lead_bin=name, COR=bivariate_cor(*args) if ok else np.nan,
                         RMSE=bivariate_rmse(*args) if ok else np.nan, n_samples=int(m.sum())))
    return pd.DataFrame(rows)


def rmm_ensemble_scores(fcst: xr.Dataset, obs: xr.Dataset, init: str) -> xr.Dataset:
    """Per lead day of one hindcast ensemble (``fcst`` member x time).

    Ensemble-mean bivariate error, mean member error, ensemble spread
    (sqrt of the unbiased member variance of the RMM vector), amplitude and
    phase errors of the ensemble mean, and the observed amplitude.
    """
    f, o = xr.align(fcst, obs, join="inner", exclude=["member"])
    m1, m2 = f.RMM1.mean("member"), f.RMM2.mean("member")
    o1, o2 = o.RMM1, o.RMM2
    n = f.sizes["member"]
    member_err = np.sqrt((f.RMM1 - o1) ** 2 + (f.RMM2 - o2) ** 2)
    out = xr.Dataset({
        "error_ensmean": np.sqrt((m1 - o1) ** 2 + (m2 - o2) ** 2),
        "error_member_mean": member_err.mean("member"),
        "spread": np.sqrt((((f.RMM1 - m1) ** 2 + (f.RMM2 - m2) ** 2).sum("member")) / max(n - 1, 1)),
        "amp_error_ensmean": np.sqrt(m1 ** 2 + m2 ** 2) - o.RMM_amp,
        "phase_error_ensmean": xr.apply_ufunc(phase_error, m1, m2, o1, o2),
        "obs_amp": o.RMM_amp,
        "obs_phase": o.RMM_phase,
    })
    return out.assign_coords(lead=("time", lead_days(out.time.values, init))).swap_dims(time="lead").assign_attrs(
        init=init, n_members=n)


def rmm_window_skill(fcst: xr.Dataset, obs: xr.Dataset, init: str, windows=DEFAULT_LEAD_BINS) -> pd.DataFrame:
    """COR/RMSE pooled over lead windows: ensemble mean, and the mean over members."""
    f, o = xr.align(fcst, obs, join="inner", exclude=["member"])
    lead = lead_days(f.time.values, init)
    rows = []
    for name, a, b in windows:
        m = (lead >= a) & (lead <= b)
        if not m.any():
            continue
        o1, o2 = o.RMM1.values[m], o.RMM2.values[m]
        e1, e2 = f.RMM1.mean("member").values[m], f.RMM2.mean("member").values[m]
        mem = [(bivariate_cor(f.RMM1.values[k, m], f.RMM2.values[k, m], o1, o2),
                bivariate_rmse(f.RMM1.values[k, m], f.RMM2.values[k, m], o1, o2)) for k in range(f.sizes["member"])]
        rows.append(dict(window=name, n_days=int(m.sum()), COR_ensmean=bivariate_cor(e1, e2, o1, o2),
                         RMSE_ensmean=bivariate_rmse(e1, e2, o1, o2), COR_member_mean=np.mean([c for c, _ in mem]),
                         RMSE_member_mean=np.mean([r for _, r in mem]), obs_amp_mean=float(np.mean(np.hypot(o1, o2)))))
    return pd.DataFrame(rows)


def plot_running_skill(curves: dict[str, xr.Dataset], *, style_map=None, cor_ref=0.5) -> plt.Figure:
    style_map = style_map or STYLE_MAP
    fig, axes = plt.subplots(2, 1, figsize=(10.8, 7), sharex=True)
    for name, ds in curves.items():
        axes[0].plot(ds.lead_day, ds.COR, label=name, **style_map.get(name, {}))
        axes[1].plot(ds.lead_day, ds.RMSE, **style_map.get(name, {}))
    any_ds = next(iter(curves.values()))
    axes[0].axhline(cor_ref, color="k", linewidth=0.8, linestyle="--")
    axes[1].axhline(np.sqrt(2), color="k", linewidth=0.8, linestyle=":")
    axes[0].set_ylim(-1, 1)
    axes[0].set_ylabel("Bivariate COR")
    axes[1].set_ylabel("Bivariate RMSE")
    axes[1].set_xlabel("Lead day (start of forward window)")
    axes[0].set_title(f"(a) RMM skill, {any_ds.window_days}-day forward window, obs amplitude ≥ "
                      f"{any_ds.amp_thresh:g}", loc="left")
    axes[1].set_title(r"(b) dotted: $\sqrt{2}$ (climatological forecast)", loc="left")
    axes[0].legend(ncol=3, fontsize=8, frameon=False)
    for ax in axes:
        ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig


def plot_leadbin_skill(table: pd.DataFrame, *, style_map=None, cor_ref=0.5) -> plt.Figure:
    """``table`` with columns experiment, lead_bin, COR, RMSE, n_samples."""
    style_map = style_map or STYLE_MAP
    exps = list(dict.fromkeys(table.experiment))
    bins = list(dict.fromkeys(table.lead_bin))
    x = np.arange(len(bins))
    width = 0.8 / len(exps)
    fig, axes = plt.subplots(2, 1, figsize=(12.8, 7.2), sharex=True)
    for j, e in enumerate(exps):
        sub = table[table.experiment == e].set_index("lead_bin").loc[bins]
        for ax, col in zip(axes, ("COR", "RMSE")):
            ax.bar(x - 0.4 + width * (j + 0.5), sub[col], width, label=e, color=style_map.get(e, {}).get("color"),
                   alpha=0.85, hatch="//" if style_map.get(e, {}).get("linestyle") in ("--", ":") else None)
    axes[0].axhline(cor_ref, color="k", linewidth=0.8, linestyle="--")
    axes[0].set_ylim(-1, 1)
    axes[0].set_ylabel("Bivariate COR")
    axes[1].set_ylabel("Bivariate RMSE")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(bins)
    axes[0].set_title("(a) RMM skill by lead bin", loc="left")
    axes[1].set_title("(b)", loc="left")
    axes[0].legend(ncol=4, fontsize=8, frameon=False)
    for ax in axes:
        ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    return fig


def plot_ensemble_rmm(fcsts: dict[str, xr.Dataset], obs: xr.Dataset, init: str, *, colors=None,
                      max_lead: int | None = None, title: str = "") -> plt.Figure:
    """(a) phase diagram from ``init``: observed, members (thin) and ensemble mean;
    (b) ensemble-mean error and spread; (c) amplitude error; (d) phase error, vs lead day."""
    colors = colors or {}
    fig = plt.figure(figsize=(17, 9.5))
    gs = fig.add_gridspec(3, 2, width_ratios=(1.1, 1.4))
    ax0 = fig.add_subplot(gs[:, 0])
    axes = [fig.add_subplot(gs[i, 1]) for i in range(3)]
    t0 = np.datetime64(init)
    o = obs.sel(time=slice(t0, None))
    if max_lead is not None:
        o = o.isel(time=slice(0, max_lead + 1))
    ax0.plot(o.RMM1, o.RMM2, "k-", linewidth=2.2, label="OBS")
    ax0.plot(o.RMM1[0], o.RMM2[0], "ko", markersize=7)
    for k, (name, f) in enumerate(fcsts.items()):
        c = colors.get(name, f"C{k}")
        f = f.sel(time=slice(t0, o.time.values[-1]))
        for m in f.member.values:
            ax0.plot(f.RMM1.sel(member=m), f.RMM2.sel(member=m), color=c, linewidth=0.3, alpha=0.35)
        ax0.plot(f.RMM1.mean("member"), f.RMM2.mean("member"), color=c, linewidth=2, label=f"{name} mean")
        s = rmm_ensemble_scores(f, obs, init)
        if max_lead is not None:
            s = s.sel(lead=slice(0, max_lead))
        axes[0].plot(s.lead, s.error_ensmean, color=c, label=name)
        axes[0].plot(s.lead, s.spread, color=c, linestyle="--")
        axes[1].plot(s.lead, s.amp_error_ensmean, color=c)
        axes[2].plot(s.lead, s.phase_error_ensmean, color=c)
    th = np.linspace(0, 2 * np.pi, 200)
    ax0.plot(np.cos(th), np.sin(th), color="0.5", linewidth=0.8)
    lim = 3.2
    ax0.plot([-lim, lim], [-lim, lim], color="0.8", linewidth=0.5)
    ax0.plot([-lim, lim], [lim, -lim], color="0.8", linewidth=0.5)
    ax0.axhline(0, color="0.8", linewidth=0.5)
    ax0.axvline(0, color="0.8", linewidth=0.5)
    for p, (x, y) in enumerate([(-0.8, -0.35), (-0.35, -0.8), (0.35, -0.8), (0.8, -0.35),
                                (0.8, 0.35), (0.35, 0.8), (-0.35, 0.8), (-0.8, 0.35)], start=1):
        ax0.text(x * lim * 0.95, y * lim * 0.95, str(p), color="0.5", ha="center", va="center")
    ax0.set_xlim(-lim, lim)
    ax0.set_ylim(-lim, lim)
    ax0.set_aspect("equal")
    ax0.set_xlabel("RMM1")
    ax0.set_ylabel("RMM2")
    ax0.set_title(f"(a) RMM from {init} (dot: start)", loc="left")
    ax0.legend(fontsize=8, frameon=False, loc="lower left")
    axes[0].axhline(np.sqrt(2), color="k", linewidth=0.6, linestyle=":")
    axes[0].set_title(r"(b) Ensemble-mean error (solid) and spread (dashed); dotted: $\sqrt{2}$", loc="left")
    axes[1].set_title("(c) Amplitude error of the ensemble mean", loc="left")
    axes[2].set_title("(d) Phase error of the ensemble mean (deg; + = ahead)", loc="left")
    axes[1].axhline(0, color="k", linewidth=0.6)
    axes[2].axhline(0, color="k", linewidth=0.6)
    axes[2].set_xlabel("Lead day")
    axes[0].legend(fontsize=8, frameon=False, ncol=2)
    for ax in axes:
        ax.grid(True, alpha=0.3)
    fig.suptitle(title)
    fig.tight_layout()
    return fig
