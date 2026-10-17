"""QBO time-height sections and amplitudes for the AMIP hindcasts.

For each source (experiment or reanalysis) and variable (``ua`` zonal wind,
``ta`` temperature) the monthly field is reduced to the equatorial
(5S-5N) zonal mean on the model pressure levels and saved with its
amplitude profiles:

* ``qbo(time, plev)``: the time-height section (plev in hPa);
* ``amp_p2p``: half the peak-to-peak range over the record;
* ``amp_rms`` / ``amp_rms_detrend``: RMS of anomalies about the record mean
  (optionally after removing a linear trend).

For temperature the amplitudes are taken after removing the month-of-year
climatology (``remove_seasonal_cycle``), so the annual cycle does not count
as QBO signal; the saved section stays raw. FFT/PSD/wavelet amplitudes need
>= 8 years and are written as NaN placeholders for the 2-year record unless
``enable_spectral`` is set on a long record.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import xarray as xr
from scipy.signal import detrend as scipy_detrend

from configs.hind_experiment_config import (
    DEFAULT_INPUT_DATA_ROOT,
    END_DATE,
    MODEL_PLEV_PA,
    PERIOD,
    START_DATE,
    check_experiments,
)
from util.hind_model_data import open_model_monthly, select_plev
from util.hind_observation_data import (
    DEFAULT_OBSERVATION_ROOT,
    REANALYSES,
    open_reference_monthly,
)

QBO_VARIABLES = {
    "ua": dict(label="U", units=r"m s$^{-1}$", remove_seasonal_cycle=False, long_name="zonal wind"),
    "ta": dict(label="T", units="K", remove_seasonal_cycle=True, long_name="temperature"),
}


@dataclass(frozen=True)
class QboConfig:
    variable: str = "ua"
    lat_band: tuple[float, float] = (-5.0, 5.0)
    remove_seasonal_cycle: bool = False   # for the amplitudes only
    pow_spec_lev_hpa: float = 20.0
    min_months_for_spectral: int = 96

    @classmethod
    def for_variable(cls, variable: str, **kwargs) -> "QboConfig":
        """Config with the variable's ``remove_seasonal_cycle`` from ``QBO_VARIABLES``."""
        return cls(variable=variable, remove_seasonal_cycle=QBO_VARIABLES[variable]["remove_seasonal_cycle"], **kwargs)


class QboAnalyzer:
    """Equatorial time-height section and amplitude metrics of one monthly field."""

    def __init__(self, cfg: QboConfig):
        self.cfg = cfg

    def reduce(self, da: xr.DataArray) -> xr.DataArray:
        """Lazy zonal mean over the latitude band (time, plev, lat): the expensive part of ``section``."""
        da = da.sel(lat=slice(*self.cfg.lat_band))
        return da.mean("lon") if "lon" in da.dims else da

    def section(self, da: xr.DataArray) -> xr.DataArray:
        """5S-5N, cos-latitude weighted zonal mean (time, plev[hPa]); computed. Accepts ``reduce`` output."""
        lat0, lat1 = self.cfg.lat_band
        da = self.reduce(da)
        weights = np.cos(np.deg2rad(da.lat))
        qbo = da.weighted(weights).mean("lat").transpose("time", "plev").compute()
        qbo = qbo.assign_coords(plev=qbo.plev.values / 100.0).sortby("plev")
        qbo.plev.attrs["units"] = "hPa"
        qbo.name = "qbo"
        qbo.attrs = {
            "description": f"Equatorial ({abs(lat0):g}S-{lat1:g}N) zonal-mean {self.cfg.variable} time-height section",
            "units": da.attrs.get("units", ""),
        }
        return qbo

    @staticmethod
    def remove_monthly_clim(x: xr.DataArray) -> xr.DataArray:
        g = x.groupby("time.month")
        out = g - g.mean("time")
        return out.drop_vars("month", errors="ignore")

    @staticmethod
    def _amp_p2p(qbo: xr.DataArray) -> xr.DataArray:
        amp = (0.5 * (qbo.max("time") - qbo.min("time"))).rename("amp_p2p")
        amp.attrs["long_name"] = "Half peak-to-peak amplitude over analysis window"
        return amp

    @staticmethod
    def _amp_rms(qbo: xr.DataArray, detrend: bool) -> xr.DataArray:
        x = qbo
        if detrend:
            x = xr.DataArray(scipy_detrend(qbo.transpose("time", ...).values, axis=0),
                             coords=qbo.coords, dims=qbo.dims)
        amp = np.sqrt(((x - x.mean("time")) ** 2).mean("time"))
        amp = amp.rename("amp_rms_detrend" if detrend else "amp_rms")
        amp.attrs.update(long_name="RMS amplitude of anomalies over analysis window", detrended=str(detrend))
        return amp

    def _spectral(self, qbo: xr.DataArray) -> dict[str, xr.DataArray]:
        """20-40 month FFT amplitude, PSD and wavelet spectrum at ~20 hPa (long records only)."""
        nan_p = xr.full_like(qbo.isel(time=0, drop=True), np.nan).rename("amp_20_40_fft")
        out = {
            "amp_20_40_fft": nan_p,
            "psd_period": xr.DataArray([np.nan], dims="period"),
            "psd_amp": xr.DataArray([np.nan], dims="period"),
            "wave_period": xr.DataArray([np.nan], dims="period_wavelet"),
            "wavelet": xr.DataArray([np.nan], dims="period_wavelet"),
        }
        if qbo.sizes["time"] < self.cfg.min_months_for_spectral:
            return out
        q = qbo.transpose("time", "plev").values
        months = np.arange(q.shape[0]) % 12
        deseason = q - np.stack([np.nanmean(q[months == m], axis=0) for m in range(12)])[months]
        n = q.shape[0]
        freq = np.arange(n // 2) / n
        with np.errstate(divide="ignore"):
            period = 1.0 / freq
        spec = np.fft.fft(deseason, axis=0)[: n // 2]
        band = (period >= 20) & (period <= 40)
        psd = 2 * np.nansum((spec * np.conj(spec)).real[band], axis=0)
        out["amp_20_40_fft"] = nan_p.copy(data=np.sqrt(2 * psd) * (freq[1] - freq[0]))
        ilev = int(np.argmin(np.abs(qbo.plev.values - self.cfg.pow_spec_lev_hpa)))
        x = deseason[:, ilev]
        nfft = int(2 ** np.ceil(np.log2(n)))
        x0 = np.fft.fft(x, n=nfft) / n
        f = np.arange(0, nfft / 2 + 1) / nfft
        with np.errstate(divide="ignore"):
            per = 1 / f
        amp0 = 2 * np.abs(x0[: nfft // 2 + 1])
        period_new = np.concatenate((np.arange(2.0, 33.0), np.arange(34.0, 100.0, 2.0)))
        out["psd_period"] = xr.DataArray(period_new, dims="period")
        out["psd_amp"] = xr.DataArray(np.interp(period_new, per[::-1][:-1], amp0[::-1][:-1]), dims="period")
        try:
            import pywt
            data = scipy_detrend(qbo.transpose("time", "plev").values[:, ilev] - np.nanmean(q[:, ilev]))
            widths = np.arange(1, 56)
            cfs, _ = pywt.cwt(data, scales=widths, wavelet="cmor1.5-1.0")
            out["wave_period"] = xr.DataArray(widths, dims="period_wavelet")
            out["wavelet"] = xr.DataArray(np.sqrt(np.mean(np.abs(cfs) ** 2, axis=1)), dims="period_wavelet")
        except ImportError:
            pass
        return out

    def compute_metrics(self, name: str, da: xr.DataArray, enable_spectral: bool = False) -> xr.Dataset:
        qbo = self.section(da)
        qbo_amp = self.remove_monthly_clim(qbo) if self.cfg.remove_seasonal_cycle else qbo
        ds = xr.Dataset({
            "qbo": qbo,
            "amp_p2p": self._amp_p2p(qbo_amp),
            "amp_rms": self._amp_rms(qbo_amp, detrend=False),
            "amp_rms_detrend": self._amp_rms(qbo_amp, detrend=True),
        })
        if enable_spectral:
            ds = ds.assign(self._spectral(qbo))
        ds.attrs = {
            "name": name,
            "variable": self.cfg.variable,
            "lat_band": str(self.cfg.lat_band),
            "amplitude_seasonal_cycle_removed": int(self.cfg.remove_seasonal_cycle),
        }
        return ds


def qbo_section_path(data_dir, source: str, variable: str, period: str = PERIOD) -> Path:
    return Path(data_dir) / f"qbo_{source}_{variable}_{period}.nc"


def open_qbo_input(
    source: str,
    variable: str,
    *,
    data_root=DEFAULT_INPUT_DATA_ROOT,
    observation_root=DEFAULT_OBSERVATION_ROOT,
    start: str = START_DATE,
    end: str = END_DATE,
    plev_pa=MODEL_PLEV_PA,
) -> xr.DataArray:
    """Lazy monthly ``variable`` of an experiment or reanalysis (``ERA5``, ``ERAI``) on the model levels."""
    if source in REANALYSES:
        da = open_reference_monthly(source, variable, start=start, end=end, observation_root=observation_root)
    else:
        check_experiments([source])
        da = open_model_monthly(source, variable, data_root=data_root, start=start, end=end)
    return select_plev(da, plev_pa)
