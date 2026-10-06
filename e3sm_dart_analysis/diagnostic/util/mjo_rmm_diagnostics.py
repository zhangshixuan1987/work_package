"""RMM/EOF diagnostics used by the maintained MJO notebooks."""

import datetime
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from scipy.signal import butter, filtfilt
from sklearn.decomposition import PCA


class MJODiagnosticAnalyzer:
    def __init__(self, mjoindx_dir, olr_data, u850_data, v850_data, u200_data, v200_data,
                 time, frequency, output_dir, figure_dir=None, rmm_cache_dir=None,
                 region_bounds=(-15, 15), months=None, years=None,
                 reload_mjo=False, overwrite_eofs=False):
        self.olr = olr_data
        self.u850 = u850_data
        self.v850 = v850_data
        self.u200 = u200_data
        self.v200 = v200_data
        self.time = time
        self.output_dir = output_dir
        self.figure_dir = figure_dir or output_dir
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.figure_dir, exist_ok=True)
        self.region_bounds = region_bounds
        self.mjoindx_dir = mjoindx_dir
        self.overwrite_eofs = overwrite_eofs
        self.dt_hours = 6 if frequency == "6hourly" else 24
        self.obs_rmm = self.get_mjo(
            mjoindx_dir, cache_dir=rmm_cache_dir, months=months, years=years,
            reload_data=reload_mjo,
        )

    @staticmethod
    def get_mjo(mjoindx_dir, cache_dir=None, months=None, years=None, reload_data=False):
        cache_dir = cache_dir or mjoindx_dir
        os.makedirs(cache_dir, exist_ok=True)
        source_txt = os.path.join(mjoindx_dir, 'rmm.74toRealtime.txt')
        local_data_txt = os.path.join(cache_dir, 'rmm.74toRealtime.txt')
        local_data_nc = os.path.join(cache_dir, 'rmm.74toRealtime.nc')
        col_names = ['year', 'month', 'day', 'RMM1', 'RMM2', 'phase', 'amplitude', 'Missing Value']

        if os.path.exists(local_data_nc) and not reload_data:
            ds = xr.open_dataset(local_data_nc)
        else:
            source = source_txt if os.path.exists(source_txt) and not reload_data else local_data_txt
            if os.path.exists(source) and not reload_data:
                df = pd.read_csv(source, skiprows=2, names=col_names, sep=r'\s+')
            else:
                url = 'http://www.bom.gov.au/climate/mjo/graphics/rmm.74toRealtime.txt'
                df = pd.read_csv(url, skiprows=2, names=col_names, sep=r'\s+')
                df.to_csv(local_data_txt, sep=' ', index=False, header=False)

            df.index = [datetime.datetime(int(y), int(m), int(d), 12) for y, m, d in zip(df.year, df.month, df.day)]
            df = df[['RMM1', 'RMM2', 'phase', 'amplitude']]
            df[df >= 999] = np.nan
            ds = xr.Dataset.from_dataframe(df).rename({'index': 'time'})
            ds.to_netcdf(local_data_nc)

        if months:
            ds = ds.sel(time=ds.time.dt.month.isin(months))
        if years:
            ds = ds.sel(time=ds.time.dt.year.isin(years))
        return ds

    def _bandpass_filter(self, data, low=1 / 100, high=1 / 20):
        data = data.interpolate_na(dim="time", method="linear", fill_value="extrapolate")
        nyq = 0.5 / (self.dt_hours / 24)
        Wn = [low / nyq, high / nyq]
        b, a = butter(N=2, Wn=Wn, btype='band')
        padlen = 3 * max(len(b), len(a))
        if data.sizes["time"] <= padlen:
            print(f"[WARN] Skipping filter: time length ({data.sizes['time']}) < padlen ({padlen})")
            return data
        return xr.apply_ufunc(
            lambda x: filtfilt(b, a, x, axis=0),
            data,
            input_core_dims=[["time"]],
            output_core_dims=[["time"]],
            vectorize=True,
            dask="parallelized",
            output_dtypes=[data.dtype]
        )

    def _load_or_compute_eofs(self, data_stack, member_idx=0):
        varnames = ["olr", "u850", "v850", "u200", "v200"]
        eof1_paths = [os.path.join(self.output_dir, f"eof1_en{member_idx:02d}_{v}.nc") for v in varnames]
        eof2_paths = [os.path.join(self.output_dir, f"eof2_en{member_idx:02d}_{v}.nc") for v in varnames]

        if all(os.path.exists(p) for p in eof1_paths + eof2_paths) and not self.overwrite_eofs:
            eof1 = {v: xr.open_dataarray(f) for v, f in zip(varnames, eof1_paths)}
            eof2 = {v: xr.open_dataarray(f) for v, f in zip(varnames, eof2_paths)}
            return eof1, eof2

        ordered = data_stack.transpose("time", "var", "lat", "lon")
        X = ordered.stack(features=("var", "lat", "lon")).fillna(0).values
        pca = PCA(n_components=2)
        pca.fit(X)
        eofs = pca.components_.reshape((2, 5, data_stack.sizes['lat'], data_stack.sizes['lon']))

        eof1 = {}
        eof2 = {}
        for i, var in enumerate(varnames):
            da1 = xr.DataArray(eofs[0][i], coords=[self.olr.lat, self.olr.lon], dims=["lat", "lon"], name=f"eof1_{var}")
            da2 = xr.DataArray(eofs[1][i], coords=[self.olr.lat, self.olr.lon], dims=["lat", "lon"], name=f"eof2_{var}")
            da1.to_netcdf(eof1_paths[i])
            da2.to_netcdf(eof2_paths[i])
            eof1[var], eof2[var] = da1, da2
        return eof1, eof2

    def compute_rmm_projection(self, member_idx=0):
        lat_mask = (self.olr.lat >= self.region_bounds[0]) & (self.olr.lat <= self.region_bounds[1])
        olr = self._bandpass_filter(self.olr.isel(member=member_idx).sel(lat=lat_mask))
        u850 = self._bandpass_filter(self.u850.isel(member=member_idx).sel(lat=lat_mask))
        v850 = self._bandpass_filter(self.v850.isel(member=member_idx).sel(lat=lat_mask))
        u200 = self._bandpass_filter(self.u200.isel(member=member_idx).sel(lat=lat_mask))
        v200 = self._bandpass_filter(self.v200.isel(member=member_idx).sel(lat=lat_mask))

        data_stack = xr.concat([olr, u850, v850, u200, v200], dim='var')
        eof1, eof2 = self._load_or_compute_eofs(data_stack, member_idx)

        rmm1 = sum((v - v.mean("time")) * eof1[k] for v, k in zip([olr, u850, v850, u200, v200], eof1)) \
               .sum(dim=["lat", "lon"])
        rmm2 = sum((v - v.mean("time")) * eof2[k] for v, k in zip([olr, u850, v850, u200, v200], eof2)) \
               .sum(dim=["lat", "lon"])
        return rmm1, rmm2

    def plot_phase_space(self, rmm1, rmm2, label="forecast"):
        amp = np.sqrt(rmm1 ** 2 + rmm2 ** 2)
        plt.figure(figsize=(6, 6))
        plt.plot(rmm1, rmm2, label=label, color='blue')
        if self.obs_rmm is not None:
            plt.plot(self.obs_rmm["RMM1"], self.obs_rmm["RMM2"], label="OBS", color='black', linestyle='--')
        plt.xlabel("RMM1")
        plt.ylabel("RMM2")
        plt.title(f"MJO Phase Space: {label}")
        plt.grid(True)
        plt.axhline(0, color='gray', linewidth=0.5)
        plt.axvline(0, color='gray', linewidth=0.5)
        plt.legend()
        outpath = os.path.join(self.figure_dir, f"mjo_phase_{label}.png")
        plt.savefig(outpath)
        plt.close()
        return outpath

    def plot_amplitude_distribution(self, amplitudes):
        plt.figure(figsize=(8, 5))
        for i, amp in enumerate(amplitudes):
            plt.hist(amp, bins=30, alpha=0.3, label=f"member{i:02d}", density=True)
        if self.obs_rmm is not None:
            obs_amp = np.sqrt(self.obs_rmm["RMM1"]**2 + self.obs_rmm["RMM2"]**2)
            plt.hist(obs_amp, bins=30, color="black", alpha=0.5, label="OBS", density=True, histtype='step', linewidth=2)
        plt.xlabel("RMM Amplitude")
        plt.ylabel("Probability Density")
        plt.title("RMM Amplitude Distribution (All Members)")
        plt.legend()
        plt.grid(True)
        plt.savefig(os.path.join(self.figure_dir, "rmm_amplitude_distribution.png"))
        plt.close()

    def analyze_all_members(self):
        num_members = self.olr.sizes['member']
        result_paths = []
        all_rmm1 = []
        all_rmm2 = []
        all_amp = []
        acc_scores = []

        for m in range(num_members):
            rmm1, rmm2 = self.compute_rmm_projection(member_idx=m)
            rmm_ds = xr.Dataset({"RMM1": rmm1, "RMM2": rmm2})
            nc_path = os.path.join(self.output_dir, f"rmm_member{m:02d}.nc")
            rmm_ds.to_netcdf(nc_path)
            fig_path = self.plot_phase_space(rmm1, rmm2, label=f"member{m:02d}")
            result_paths.append((nc_path, fig_path))
            all_rmm1.append(rmm1)
            all_rmm2.append(rmm2)
            all_amp.append(np.sqrt(rmm1 ** 2 + rmm2 ** 2))

            if self.obs_rmm is not None:
                acc1 = xr.corr(rmm1, self.obs_rmm["RMM1"], dim="time")
                acc2 = xr.corr(rmm2, self.obs_rmm["RMM2"], dim="time")
                acc_scores.append({"member": m, "ACC_RMM1": float(acc1), "ACC_RMM2": float(acc2)})
                print(f"[INFO] Member {m:02d} ACC: RMM1={float(acc1):.3f}, RMM2={float(acc2):.3f}")

        self.plot_amplitude_distribution(all_amp)

        if acc_scores:
            acc_df = pd.DataFrame(acc_scores)
            acc_path = os.path.join(self.output_dir, "mjo_skill_scores.csv")
            acc_df.to_csv(acc_path, index=False)
            print(f"[INFO] Saved ACC skill scores to: {acc_path}")

        return result_paths
