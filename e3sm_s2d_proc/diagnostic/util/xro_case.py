"""XRO (extended nonlinear recharge oscillator) helpers for the RTENSO2026 7a-7d notebooks.

The XRO state vector of Zhao et al. (2024, Nature, doi:10.1038/s41586-024-07534-6) is
ENSO [Nino3.4 SST T, recharge h] plus eight other climate-mode SST indices. The boxes follow
``XRO.state.calc_XRO_indices``, except that the SASD north-east box (30S-20S, 40W-20E) wraps
through 0E here (XRO's 320-380E slice stops at 360E on a 0-360 grid). The recharge h is the
5S-5N, 120E-80W mean of the 0-700 m heat content (``rtenso_case`` field "OHC700"; EN4 for
observations) instead of the 20 C isotherm depth, which the hindcasts do not save. The variable
keeps XRO's name "WWV" in the datasets but is labelled "h (OHC700)" in figures.

The module holds the configuration shared by the notebooks (``OBS_SST``, ``EVENTS``, ``CASES``,
``XRO_SETUP``) so that every notebook fits, forecasts and verifies the same way, plus cached
builders:

* ``obs_state``        observed state vector (HadISST 1.1 SST, EN4 OHC700)
* ``hindcast_state``   drift-corrected state vector of one hindcast (member, lead)
* ``fit_xro``          detrend over a training period and fit XRO
* ``event_forecast``   XRO forecast aligned with a hindcast start (lead 0 = start month)
* ``rolling_hindcasts`` out-of-sample XRO hindcasts (refit each year on the preceding years)

Indices are area means of anomalies, so they work for observations (time, lat, lon) and
hindcasts (member, lead, lat, lon) alike.
"""

from __future__ import annotations

import warnings
from contextlib import contextmanager

import cftime
import numpy as np
import pandas as pd
import xarray as xr
from XRO import XRO

from esp_lab.utils.netcdf_utils import atomic_to_netcdf, load_netcdf

from . import rtenso_case as rc

CACHE_VERSION = 1

# -----------------------------------------------------------------------------
# Shared configuration of the 7a-7d notebooks
# -----------------------------------------------------------------------------
OBS_SST = "HadISST1.1"          # rc.OBS_ALT_PRODUCTS["TS"]; heat content is EN4.2.2 (rc.FIELDS["OHC700"])
OBS_YEARS = (1981, 2026)        # observed state vector: training, initial states and verification

# Control XRO of Zhao et al. (2024): annual + semi-annual cycle, quadratic IOD term, T^2 and T*h terms
# in the ENSO T equation. Trend orders per variable ("default" for the rest): linear for SST, quadratic
# for h, as in the XRO real-time notebook. "spread" is the percentile range shown for every ensemble
# (E3SM members and XRO stochastic copies alike). A few red-noise members run away through the
# nonlinear terms (about 0.1% of the validation members, mostly beyond lead 11); a member whose
# |Nino3.4| ever exceeds "diverge_threshold" (degC; observed 1981-2026 maximum 2.6) is set to NaN.
XRO_SETUP = {
    "ac_order": 2, "maskb": ("IOD",), "maskNT": ("T2", "TH"),
    "detrend_order": {"default": 1, "WWV": 2},
    "ncopy": 100, "seed": 2026,
    "spread": (5, 95),
    "diverge_threshold": 6.0,
}

# Out-of-sample validation (7a) and the cross-validated spread calibration derived from it, shared so
# that 7b/7c apply exactly the factors 7a evaluates.
VALIDATION = {"years": (1995, 2025), "train_start": 1981, "n_month": 24, "ncopy": 50, "seed": 1995}
CALIBRATION = {"pool": 1, "min_samples": 60, "vars": ("Nino34",)}   # h is not calibrated (7a section 8)

# Events: "init_months" are the hindcast starts, "train" the XRO training years (all before the event),
# "peak_season" the target season of the summaries and sensitivity experiments.
EVENTS = {
    2015: {"label": "2015/16 El Niño", "init_months": (2, 5, 8, 11), "train": (1981, 2014),
           "peak_season": {"label": "DJF 2015/16", "months": [(2015, 12), (2016, 1), (2016, 2)]}},
    2026: {"label": "2026/27 El Niño", "init_months": (2, 5, 8), "train": (1981, 2025),
           "peak_season": {"label": "DJF 2026/27", "months": [(2026, 12), (2027, 1), (2027, 2)]}},
}

# Hindcast cases (as in 1_enso_plume). Missing runs (e.g. Reanalysis-FOSIRL in 2026) are skipped.
MODEL_PREFIX = rc.MODEL_PREFIX
CASES = {
    "Reanalysis-CLIM": {"case_prefix": f"{MODEL_PREFIX}_RTENSO2026_clim", "source": "RTENSO2026_clim",
                        "init_months": (2, 5, 8, 11), "color": "#1f77b4", "short": "CLIM"},
    "Reanalysis-OBS": {"case_prefix": f"{MODEL_PREFIX}_RTENSO2026_obs", "source": "RTENSO2026_obs",
                       "init_months": (2, 5, 8, 11), "color": "#ff7f0e", "short": "OBS"},
    "Reanalysis-NSIDC": {"case_prefix": f"{MODEL_PREFIX}_RTENSO2026_obs_ice_nsidc", "source": "RTENSO2026_obs_ice_nsidc",
                         "init_months": (2, 5, 8, 11), "color": "#d62728", "short": "NSIDC"},
    "Reanalysis-FOSIRL": {"case_prefix": f"{MODEL_PREFIX}_BruteForce", "source": "Reanalysis",
                          "init_months": (2, 5, 8, 11), "color": "#2ca02c", "short": "FOSIRL"},
}
OBS_STYLE = {"color": "k", "label": "Observed"}
XRO_STYLE = {"color": "#7b2cbf", "ls": "--", "label": "XRO"}
MONTH_NAMES = {m: pd.Timestamp(2000, m, 1).strftime("%b") for m in range(1, 13)}

# [lon0, lon1, lat0, lat1] (0-360 longitudes; lon0 > lon1 wraps through 0E)
SST_BOXES = {
    "Nino34": (190.0, 240.0, -5.0, 5.0),
    "NPMM": (200.0, 240.0, 10.0, 25.0),
    "SPMM": (250.0, 270.0, -25.0, -15.0),
    "IOB": (40.0, 100.0, -20.0, 20.0),
    "IODw": (50.0, 70.0, -10.0, 10.0),
    "IODe": (90.0, 110.0, -10.0, 0.0),
    "SIODw": (65.0, 85.0, -25.0, -10.0),
    "SIODe": (90.0, 120.0, -30.0, -5.0),
    "TNA": (305.0, 345.0, 5.0, 25.0),
    "ATL3": (340.0, 360.0, -3.0, 3.0),
    "SASD_SW": (300.0, 360.0, -45.0, -35.0),
    "SASD_NE": (320.0, 20.0, -30.0, -20.0),
}
WWV_BOX = (120.0, 280.0, -5.0, 5.0)

# Order matters: XRO takes ENSO T and h as the first two states.
XRO_VARS = ["Nino34", "WWV", "NPMM", "SPMM", "IOB", "IOD", "SIOD", "TNA", "ATL3", "SASD"]
LABELS = {
    "Nino34": "Niño3.4 (ENSO T)", "WWV": "Recharge h (eq. Pacific OHC700)",
    "NPMM": "North Pacific Meridional Mode", "SPMM": "South Pacific Meridional Mode",
    "IOB": "Indian Ocean Basin", "IOD": "Indian Ocean Dipole", "SIOD": "Southern Indian Ocean Dipole",
    "TNA": "Tropical North Atlantic", "ATL3": "Atlantic Niño (ATL3)", "SASD": "South Atlantic Subtropical Dipole",
}
SHORT = {"Nino34": "Niño3.4", "WWV": "h (OHC700)"}
UNITS = {v: "°C" for v in XRO_VARS} | {"WWV": "GJ m$^{-2}$"}


@contextmanager
def quiet_xro():
    """Silence the xarray deprecation and divide warnings raised inside the XRO package only."""
    with warnings.catch_warnings():
        for category in (FutureWarning, DeprecationWarning, RuntimeWarning):
            warnings.filterwarnings("ignore", category=category)
        yield


def short(var):
    return SHORT.get(var, var)


def available_cases(year, init):
    return [case for case, cfg in CASES.items() if init in cfg["init_months"]]


def state_vector(sst, hc):
    """XRO state vector (Dataset of XRO_VARS) from SST (°C) and heat-content (GJ m-2) anomalies."""
    b = {name: rc.box_mean(sst, box) for name, box in SST_BOXES.items()}
    out = {
        "Nino34": b["Nino34"], "WWV": rc.box_mean(hc, WWV_BOX), "NPMM": b["NPMM"], "SPMM": b["SPMM"],
        "IOB": b["IOB"], "IOD": b["IODw"] - b["IODe"], "SIOD": b["SIODe"] - b["SIODw"],
        "TNA": b["TNA"], "ATL3": b["ATL3"], "SASD": b["SASD_SW"] - b["SASD_NE"],
    }
    ds = xr.Dataset(out)[XRO_VARS]
    for var in XRO_VARS:
        ds[var].attrs = {"long_name": LABELS[var], "units": UNITS[var]}
    return ds


# -----------------------------------------------------------------------------
# Time helpers
# -----------------------------------------------------------------------------
def decimal_year(times):
    """Mid-month decimal year of cftime/datetime64 monthly times."""
    idx = xr.DataArray(np.atleast_1d(times), dims="t")
    return idx.dt.year.values + (idx.dt.month.values - 0.5) / 12.0


def to_month_start(ds):
    """Same series with datetime64 month-start times (XRO reads ``time.dt.month``)."""
    t = ds["time"]
    times = pd.to_datetime({"year": t.dt.year.values, "month": t.dt.month.values, "day": 1})
    return ds.assign_coords(time=times.values)


def add_months(year, month, n):
    k = month - 1 + n
    return year + k // 12, k % 12 + 1


def month_index(ds, year, month):
    """Position of (year, month) along ``time``, or None."""
    hit = np.flatnonzero((ds["time"].dt.year.values == year) & (ds["time"].dt.month.values == month))
    return int(hit[0]) if hit.size else None


def obs_value(X, var, year, month):
    """Observed value of ``var`` in (year, month), NaN if not observed."""
    i = month_index(X, year, month)
    return np.nan if i is None else float(X[var].isel(time=i))


# -----------------------------------------------------------------------------
# Cached state vectors
# -----------------------------------------------------------------------------
def _cached(path, expected, mode, build):
    if mode != "rebuild" and rc._cache_ok(path, expected):
        return load_netcdf(path)
    if mode == "require":
        raise FileNotFoundError(f"missing or stale cache {path}")
    ds = build()
    if ds is None:
        return None
    ds.attrs.update(expected)
    path.parent.mkdir(parents=True, exist_ok=True)
    atomic_to_netcdf(ds, path)
    return ds


def obs_state(mode="auto"):
    """Observed state vector (time) of HadISST 1.1 SST and EN4 OHC700 anomalies from 1981-2010.

    Rebuild (``mode="rebuild"``) after new observations are added to the archive.
    """
    path = rc.diag_dir("observations", "xro",
                       f"{OBS_SST}_EN4.2.2_XRO_state_y{OBS_YEARS[0]}-{OBS_YEARS[1]}_"
                       f"clim{rc.CLIM_YEARS[0]}-{rc.CLIM_YEARS[1]}_{rc.GRID_TAG}.nc")
    expected = {"cache_version": CACHE_VERSION, "sst_product": OBS_SST, "heat_content_product": "EN4.2.2",
                "start_year": OBS_YEARS[0], "end_year": OBS_YEARS[1]}

    def build():
        sst = rc.obs_anomaly("TS", years=OBS_YEARS, mode=mode, product=OBS_SST)
        hc = rc.obs_anomaly("OHC700", years=OBS_YEARS, mode=mode)
        return state_vector(sst, hc).load()

    return _cached(path, expected, mode, build)


def last_observed(X, var):
    t = X[var].dropna("time")["time"].values[-1]
    return t.year, t.month


def hindcast_state(case, year, init, method="lead_drift", mode="auto"):
    """Drift-corrected state vector (member, lead) of one hindcast with a target-month ``time``
    coordinate along ``lead``, or None when the runs do not exist. SST is land-masked as the
    observations (EAM LANDFRAC < 0.5)."""
    cfg = CASES[case]
    path = rc.diag_dir(cfg["source"], "xro_state",
                       f"{cfg['source']}_init{init:02d}_y{year}_XRO_state_{method}_m{rc.NENS}_l{rc.NLEAD}_{rc.GRID_TAG}.nc")
    expected = {"cache_version": CACHE_VERSION, "case_prefix": cfg["case_prefix"], "init_year": year,
                "init_month": init, "method": method}

    def build():
        sst = rc.hindcast_anomaly(cfg["case_prefix"], init, "TS", year=year, method=method, mode=mode)
        if sst is None:
            return None
        hc = rc.hindcast_anomaly(cfg["case_prefix"], init, "OHC700", year=year, method=method, mode=mode)
        if hc is None:
            return None
        ocean = rc.land_fraction() < 0.5
        X = state_vector(sst.where(ocean), hc).load()
        X = X.drop_vars([c for c in X.coords if c not in ("member", "lead")])
        X.attrs["climatology_method"] = sst.attrs.get("climatology_method", "")
        return X

    if mode != "rebuild" and not rc._cache_ok(path, expected) and rc.member_file(cfg["case_prefix"], year, init, 0, "TS") is None:
        return None
    X = _cached(path, expected, mode, build)
    if X is None:
        return None
    X = X.drop_vars([c for c in X.coords if c not in ("member", "lead")])
    return X.assign_coords(time=("lead", rc._target_times(year, init, X.sizes["lead"])))


def hindcast_states(years=None, method="lead_drift", mode="auto", verbose=True):
    """{(year, case, init): Dataset} for every existing hindcast of ``EVENTS``."""
    out = {}
    for year, event in EVENTS.items():
        if years is not None and year not in years:
            continue
        for case in CASES:
            for init in event["init_months"]:
                if init not in CASES[case]["init_months"]:
                    continue
                X = hindcast_state(case, year, init, method=method, mode=mode)
                if X is None:
                    if verbose:
                        print(f"[{year} {case} init{init:02d}] not available")
                    continue
                out[(year, case, init)] = X
    return out


# -----------------------------------------------------------------------------
# Trend (fitted on the training period, extrapolated over the forecasts)
# -----------------------------------------------------------------------------
def fit_trend(X, period, orders):
    """Polynomial trend coefficients per variable over ``period`` (first_year, last_year)."""
    train = X.sel(time=slice(f"{period[0]}-01", f"{period[1]}-12"))
    t = decimal_year(train["time"].values)
    coeffs = {}
    for var in X.data_vars:
        y = train[var].values
        ok = np.isfinite(y)
        coeffs[var] = np.polyfit(t[ok], y[ok], orders.get(var, orders.get("default", 1)))
    return coeffs


def trend_values(coeffs, times, like=None):
    """Trend Dataset at ``times`` (dims follow ``like``'s time coordinate when given)."""
    t = decimal_year(times)
    dims = ("time",) if like is None else like["time"].dims
    data = {var: (dims, np.polyval(c, t).reshape([like.sizes[d] for d in dims] if like is not None else -1))
            for var, c in coeffs.items()}
    return xr.Dataset(data)


def _trend_on_init_lead(coeffs, fc):
    """Trend on an (init, lead) grid: target month = init month + lead."""
    init = pd.DatetimeIndex(fc["init"].values)
    t = (init.year.values + (init.month.values - 0.5) / 12.0)[:, None] + fc["lead"].values[None, :] / 12.0
    return xr.Dataset({var: (("init", "lead"), np.polyval(c, t)) for var, c in coeffs.items()},
                      coords={"init": fc["init"], "lead": fc["lead"]})


# -----------------------------------------------------------------------------
# XRO fit and forecasts
# -----------------------------------------------------------------------------
def fit(X_train, ac_order=2, maskb=("IOD",), maskNT=("T2", "TH")):
    """Fit XRO to a detrended (time) state vector starting in any month. Returns (model, fit_ds)."""
    model = XRO(ncycle=12, ac_order=ac_order)
    with quiet_xro():
        fit_ds = model.fit_matrix(to_month_start(X_train[XRO_VARS]), maskb=list(maskb), maskNT=list(maskNT))
    return model, fit_ds


def fit_xro(X, train, orders=None, setup=XRO_SETUP):
    """Detrend ``X`` with the ``train`` period trend and fit XRO on the training years.

    Returns {"model", "fit", "trend", "X_dt", "train", "orders"}.
    """
    orders = orders or setup["detrend_order"]
    coeffs = fit_trend(X, train, orders)
    X_dt = X - trend_values(coeffs, X["time"].values)
    sample = X_dt.sel(time=slice(f"{train[0]}-01", f"{train[1]}-12"))
    model, fit_ds = fit(sample, ac_order=setup["ac_order"], maskb=setup["maskb"], maskNT=setup["maskNT"])
    return {"model": model, "fit": fit_ds, "trend": coeffs, "X_dt": X_dt, "train": tuple(train), "orders": dict(orders)}


def forecast(model, fit_ds, X0, ic_year, ic_month, n_month, ncopy=1, noise_type="zero", seed=None):
    """XRO forecast from the state ``X0`` (Dataset of scalars) of month (ic_year, ic_month).

    Returns (lead[, member]) with lead 0 the initial month itself and a cftime ``time`` coordinate
    of target months (mid-month, as the hindcasts).
    """
    init = X0[XRO_VARS].expand_dims(time=[pd.Timestamp(ic_year, ic_month, 1)])
    with quiet_xro():
        fc = model.reforecast(fit_ds=fit_ds, init_ds=init, n_month=n_month, ncopy=ncopy,
                              noise_type=noise_type, seed=seed).isel(init=0, drop=True)
    times = [cftime.DatetimeNoLeap(*add_months(ic_year, ic_month, int(k)), 15) for k in fc["lead"].values]
    fc = fc.drop_vars([c for c in fc.coords if c not in ("lead", "member")])
    return fc.assign_coords(time=("lead", times))


def mask_diverged(ens, threshold=XRO_SETUP["diverge_threshold"]):
    """Set every variable of a stochastic member to NaN at all leads if its |Nino3.4| ever exceeds
    ``threshold``; the count is stored in attrs["diverged_members"]."""
    bad = (np.abs(ens["Nino34"]) > threshold).any("lead")
    out = ens.where(~bad)
    out.attrs = dict(ens.attrs, diverged_members=int(bad.sum()), diverge_threshold=threshold)
    return out


def initial_state(xfit, year, init):
    """Detrended observed state of the month before the hindcast start, or None if not fully observed."""
    ic_year, ic_month = add_months(year, init, -1)
    i = month_index(xfit["X_dt"], ic_year, ic_month)
    if i is None:
        return None
    X0 = xfit["X_dt"].isel(time=i, drop=True)
    return None if bool(X0.to_array().isnull().any()) else X0


def event_forecast(xfit, year, init, X0=None, stochastic=True, n_month=rc.NLEAD + 1, setup=XRO_SETUP):
    """XRO forecast aligned with a hindcast started on 1 ``init`` ``year``.

    XRO starts from the (detrended) observed state of the month before the start (or ``X0``),
    so XRO lead 1 is the start month; leads are shifted by one so that lead 0 is the start month,
    as in the hindcasts, and lead -1 is the initial month. Returns {"det", "det_dt"[, "stoch"]}
    (trend added back to "det" and "stoch"), or None without an observed initial state.
    """
    X0 = initial_state(xfit, year, init) if X0 is None else X0
    if X0 is None:
        return None
    ic_year, ic_month = add_months(year, init, -1)

    def run(**kwargs):
        fc = forecast(xfit["model"], xfit["fit"], X0, ic_year, ic_month, n_month, **kwargs)
        return fc.assign_coords(lead=fc["lead"].values - 1)

    det_dt = run()
    tr = trend_values(xfit["trend"], det_dt["time"].values, like=det_dt).assign_coords(det_dt.coords)
    out = {"det": det_dt + tr, "det_dt": det_dt}
    if stochastic:
        out["stoch"] = mask_diverged(run(ncopy=setup["ncopy"], noise_type="red", seed=setup["seed"]) + tr,
                                     setup["diverge_threshold"])
    return out


def removal_sensitivity(xfit, year, init):
    """Initial-state removal sensitivity: detrended deterministic forecast minus the forecast with one
    initial variable set to zero ("other modes": all eight non-ENSO modes together).
    XRO is nonlinear, so these sensitivities are not additive. Dataset with a ``removed`` dim."""
    X0 = initial_state(xfit, year, init)
    if X0 is None:
        return None
    full = event_forecast(xfit, year, init, X0=X0, stochastic=False)["det_dt"]
    cases = {var: X0.assign({var: 0.0 * X0[var]}) for var in XRO_VARS}
    cases["other modes"] = X0.assign({var: 0.0 * X0[var] for var in XRO_VARS[2:]})
    out = [full - event_forecast(xfit, year, init, X0=state, stochastic=False)["det_dt"] for state in cases.values()]
    return xr.concat(out, dim=pd.Index(list(cases), name="removed"))


def perturbed_forecast(xfit, year, init, delta):
    """Deterministic forecast (trend added) from the observed initial state plus ``delta``
    ({var: offset}); None without an observed initial state."""
    X0 = initial_state(xfit, year, init)
    if X0 is None:
        return None
    X0 = X0.assign({var: X0[var] + float(d) for var, d in delta.items()})
    return event_forecast(xfit, year, init, X0=X0, stochastic=False)["det"]


def simulate(xfit, nyear=100, ncopy=50, seed=2000):
    """Long stochastic XRO simulation (time, member) from the first training month (detrended)."""
    X0 = to_month_start(xfit["X_dt"].sel(time=slice(f"{xfit['train'][0]}-01", None)))[XRO_VARS].isel(time=0)
    with quiet_xro():
        return xfit["model"].simulate(fit_ds=xfit["fit"], X0_ds=X0, nyear=nyear, ncopy=ncopy, seed=seed,
                                      is_xi_stdac=False)


def ro_parameters(xfit):
    """Recharge-oscillator parameters of the ENSO (T, h) block, per calendar month and annual mean.

    XRO writes dT/dt = R T + F1 h, dh/dt = -F2 T - epsilon h (per year). The eigenvalues give the
    growth rate (R - epsilon)/2 and the angular frequency sqrt(F1 F2 - ((R + epsilon)/2)^2); the
    period does not depend on the units of h. Annual-mean growth and period use the annual-mean
    parameters.
    """
    with quiet_xro():
        p = xfit["model"].get_RO_parameters(xfit["fit"])
    p = p[["R", "F1", "F2", "epsilon"]]
    growth = (p["R"] - p["epsilon"]) / 2.0
    disc = p["F1"] * p["F2"] - ((p["R"] + p["epsilon"]) / 2.0) ** 2
    p["growth"] = growth
    p["period"] = 2 * np.pi / np.sqrt(disc.where(disc > 0))
    ann = p[["R", "F1", "F2", "epsilon"]].mean("cycle")
    ann_disc = float(ann["F1"] * ann["F2"] - ((ann["R"] + ann["epsilon"]) / 2.0) ** 2)
    annual = {k: float(ann[k]) for k in ann.data_vars}
    annual["growth"] = float((ann["R"] - ann["epsilon"]) / 2.0)
    annual["period"] = float(2 * np.pi / np.sqrt(ann_disc)) if ann_disc > 0 else np.nan
    return p, annual


# -----------------------------------------------------------------------------
# Out-of-sample validation
# -----------------------------------------------------------------------------
def rolling_hindcasts(X, years, train_start, n_month=24, orders=None, ncopy=0, seed=None,
                      setup=XRO_SETUP, mode="auto"):
    """Out-of-sample XRO hindcasts. For each year Y in ``years`` (first, last), XRO is detrended and
    fitted on ``train_start``..Y-1 only and run from every fully observed month of Y.

    Returns (det, stoch): Datasets (init, lead) and (init, lead, member) of full anomalies (the
    training trend added back), with XRO's lead convention (lead 0 = initial month); stoch is None
    when ``ncopy`` is 0.
    """
    orders = orders or setup["detrend_order"]
    tag = "_".join(f"{k}{v}" for k, v in sorted(orders.items()))
    base = rc.diag_dir("observations", "xro",
                       f"XRO_rolling_y{years[0]}-{years[1]}_train{train_start}_{tag}_n{n_month}")
    expected = {"cache_version": CACHE_VERSION, "first_year": years[0], "last_year": years[1],
                "train_start": train_start, "n_month": n_month, "orders": tag, "sst_product": OBS_SST,
                "ac_order": setup["ac_order"], "maskb": ",".join(setup["maskb"]), "maskNT": ",".join(setup["maskNT"])}
    stoch_expected = dict(expected, ncopy=ncopy, seed=str(seed), diverge_threshold=setup["diverge_threshold"])
    det_path, stoch_path = base.with_name(base.name + "_det.nc"), base.with_name(base.name + f"_stoch{ncopy}.nc")
    det_ok = mode != "rebuild" and rc._cache_ok(det_path, expected)
    stoch_ok = ncopy == 0 or (mode != "rebuild" and rc._cache_ok(stoch_path, stoch_expected))
    if det_ok and stoch_ok:
        return load_netcdf(det_path), (load_netcdf(stoch_path) if ncopy else None)
    if mode == "require":
        raise FileNotFoundError(f"missing or stale rolling hindcast cache {det_path}")

    det, stoch = [], []
    for Y in range(years[0], years[1] + 1):
        xfit = fit_xro(X, (train_start, Y - 1), orders=orders, setup=setup)
        ics = xfit["X_dt"].sel(time=slice(f"{Y}-01", f"{Y}-12"))
        ok = ics.to_array().notnull().all("variable").values
        if not ok.any():
            continue
        init_ds = to_month_start(ics.isel(time=np.flatnonzero(ok)))[XRO_VARS]
        with quiet_xro():
            fc = xfit["model"].reforecast(fit_ds=xfit["fit"], init_ds=init_ds, n_month=n_month, noise_type="zero")
        det.append(fc + _trend_on_init_lead(xfit["trend"], fc))
        if ncopy:
            with quiet_xro():
                fs = xfit["model"].reforecast(fit_ds=xfit["fit"], init_ds=init_ds, n_month=n_month, ncopy=ncopy,
                                              noise_type="red", seed=None if seed is None else seed + Y)
            stoch.append(fs + _trend_on_init_lead(xfit["trend"], fs))
    det = xr.concat(det, dim="init")
    det.attrs.update(expected)
    det_path.parent.mkdir(parents=True, exist_ok=True)
    atomic_to_netcdf(det, det_path)
    if not ncopy:
        return det, None
    stoch = mask_diverged(xr.concat(stoch, dim="init").transpose("init", "lead", "member"), setup["diverge_threshold"])
    stoch.attrs.update(stoch_expected)
    atomic_to_netcdf(stoch, stoch_path)
    return det, stoch


def target_obs(fc, X):
    """Observed values on a forecast's (init, lead) grid (target = init + lead months); NaN if unobserved."""
    init = pd.DatetimeIndex(fc["init"].values)
    leads = fc["lead"].values
    times = pd.DatetimeIndex([i + pd.DateOffset(months=int(L)) for i in init for L in leads])
    obs = to_month_start(X).reindex(time=times)
    return xr.Dataset({v: (("init", "lead"), obs[v].values.reshape(len(init), len(leads))) for v in fc.data_vars},
                      coords={"init": fc["init"], "lead": fc["lead"]})


def persistence(fc, obs):
    """Anomaly persistence forecast on the same grid: the observed initial-month value at every lead."""
    return obs.isel(lead=0, drop=True).broadcast_like(obs)


def skill(fc, obs, dim="init"):
    """ACC, RMSE and sample size along ``dim`` over pairs where both are finite."""
    ok = fc.notnull() & obs.notnull()
    f, o = fc.where(ok), obs.where(ok)
    return xr.Dataset({"acc": xr.corr(f, o, dim=dim), "rmse": np.sqrt(((f - o) ** 2).mean(dim)), "n": ok.sum(dim)})


# -----------------------------------------------------------------------------
# Probabilistic calibration of the stochastic ensemble (spread inflation)
# -----------------------------------------------------------------------------
def _month_number(times):
    t = pd.DatetimeIndex(np.atleast_1d(times))
    return t.year.values * 12 + t.month.values - 1


def required_inflation(ens, obs, interval=XRO_SETUP["spread"]):
    """Factor by which each forecast's central interval (``interval`` percentiles), scaled about the
    ensemble mean, must widen to just contain the observation (side of the observation only).
    ens: (init, lead, member); obs: (init, lead). NaN where unobserved."""
    mu = ens.mean("member")
    qlo = ens.quantile(interval[0] / 100, "member").drop_vars("quantile")
    qhi = ens.quantile(interval[1] / 100, "member").drop_vars("quantile")
    up, dn = (qhi - mu).clip(min=1e-8), (mu - qlo).clip(min=1e-8)
    return xr.where(obs >= mu, (obs - mu) / up, (mu - obs) / dn).where(obs.notnull())


def inflation_factors(r, issue, coverage=None, pool=1, min_samples=60, interval=XRO_SETUP["spread"]):
    """Lead-dependent inflation alpha(lead) for a forecast issued at ``issue`` (its initial month).

    Uses only historical forecasts whose target month (init + lead) precedes the issue month, i.e.
    verified when the forecast is issued, pooling leads within +-``pool``. alpha is the
    ``coverage`` quantile (default: the interval's nominal coverage) of the required inflation, so
    the inflated interval would have contained that fraction of the past observations. NaN with
    fewer than ``min_samples`` cases. Returns (alpha, n) DataArrays over lead.
    """
    coverage = (interval[1] - interval[0]) / 100 if coverage is None else coverage
    issue_n = _month_number(pd.Timestamp(issue))[0]
    target = _month_number(r["init"].values)[:, None] + r["lead"].values[None, :]
    vals = np.where(target < issue_n, r.values, np.nan)
    leads = r["lead"].values
    alpha, n = np.full(len(leads), np.nan), np.zeros(len(leads), int)
    for j, L in enumerate(leads):
        sample = vals[:, np.abs(leads - L) <= pool].ravel()
        sample = sample[np.isfinite(sample)]
        n[j] = sample.size
        if sample.size >= min_samples:
            alpha[j] = np.quantile(sample, coverage)
    coords = {"lead": r["lead"]}
    return xr.DataArray(alpha, coords=coords, dims="lead"), xr.DataArray(n, coords=coords, dims="lead")


def inflate(ens, alpha):
    """Calibrated ensemble: mean + alpha * (member - mean); the ensemble mean is unchanged."""
    mu = ens.mean("member")
    return mu + alpha * (ens - mu)


def crps_ensemble(ens, obs):
    """Empirical ensemble CRPS: E|X - o| - 0.5 E|X - X'| over the valid (non-NaN) members,
    NaN where unobserved."""
    x = np.sort(ens.transpose(..., "member").values, axis=-1)          # NaN members sort last
    m = np.isfinite(x).sum(-1)
    o = obs.broadcast_like(ens.isel(member=0, drop=True)).transpose(*ens.isel(member=0, drop=True).dims).values
    with warnings.catch_warnings():          # rows without an observation are all-NaN by design
        warnings.simplefilter("ignore", RuntimeWarning)
        term1 = np.nanmean(np.abs(x - o[..., None]), axis=-1)
    i = np.arange(1, x.shape[-1] + 1)
    w = np.where(i <= m[..., None], 2 * i - m[..., None] - 1, 0)
    term2 = np.nansum(x * w, axis=-1) * 2 / np.maximum(m, 1) ** 2
    crps = np.where(m > 0, term1 - 0.5 * term2, np.nan)
    dims = [d for d in ens.dims if d != "member"]
    return xr.DataArray(crps, coords={d: ens[d] for d in dims}, dims=dims).where(obs.notnull())


def event_inflation(X, year, init, var="Nino34", mode="auto"):
    """Cross-validated inflation factors for an event forecast (``event_forecast`` leads).

    Factors come from the 7a validation hindcasts verified before XRO's initial month (the month
    before the start); event lead L uses the factor of XRO lead L + 1, and lead -1 (the observed
    initial month) gets 1. Returns a DataArray over event leads (NaN where too few past cases).
    """
    v = VALIDATION
    _, stoch = rolling_hindcasts(X, v["years"], v["train_start"], v["n_month"], ncopy=v["ncopy"], seed=v["seed"], mode=mode)
    r = required_inflation(stoch[var], target_obs(stoch[[var]], X)[var])
    issue = pd.Timestamp(*add_months(year, init, -1), 1)
    alpha, _ = inflation_factors(r, issue, pool=CALIBRATION["pool"], min_samples=CALIBRATION["min_samples"])
    alpha = alpha.assign_coords(lead=alpha["lead"].values - 1)
    return xr.where(alpha["lead"] < 0, 1.0, alpha)


def calibrated_stoch(X, year, init, stoch, mode="auto"):
    """Stochastic event ensemble with CALIBRATION["vars"] inflated about the ensemble mean (other
    variables unchanged)."""
    out = stoch.copy()
    for var in CALIBRATION["vars"]:
        alpha = event_inflation(X, year, init, var, mode=mode).reindex(lead=stoch["lead"])
        out[var] = inflate(stoch[var], alpha)
    out.attrs["calibrated_vars"] = ",".join(CALIBRATION["vars"])
    return out


__all__ = [
    "VALIDATION", "CALIBRATION", "event_inflation", "calibrated_stoch",
    "OBS_SST", "OBS_YEARS", "XRO_SETUP", "EVENTS", "CASES", "OBS_STYLE", "XRO_STYLE", "MONTH_NAMES",
    "SST_BOXES", "WWV_BOX", "XRO_VARS", "LABELS", "SHORT", "UNITS", "quiet_xro", "short", "available_cases",
    "state_vector", "decimal_year", "to_month_start", "add_months", "month_index", "obs_value", "obs_state",
    "last_observed", "hindcast_state", "hindcast_states", "fit_trend", "trend_values", "fit", "fit_xro",
    "forecast", "initial_state", "event_forecast", "removal_sensitivity", "perturbed_forecast", "simulate",
    "ro_parameters", "rolling_hindcasts", "target_obs", "persistence", "skill",
    "required_inflation", "inflation_factors", "inflate", "crps_ensemble", "mask_diverged",
]
