import os
import glob
import numpy as np
import xarray as xr
import xcdat
import pandas as pd
import xskillscore as xs
from typing import Dict, Tuple
from properscoring import crps_ensemble
from datetime import datetime

class ModelDataReader:
    def __init__(self, base_path, exp_base, regnam, frequency, period=None, component="atm"):
        self.base_path = base_path
        self.exp_base = exp_base
        self.regnam = regnam
        self.frequency = frequency
        self.period = period
        self.component = component
        self.region = self._define_region(regnam)
        self.exp_dict = self._extract_exp_info(exp_base, base_path)

    @staticmethod
    def compute_total_soil_moisture(ds, h2osoi_var='H2OSOI', dzsoi_var='DZSOI') -> xr.DataArray:
        """
        Compute total column soil moisture by integrating volumetric H2OSOI over depth.

        Returns
        -------
        total_sm : xr.DataArray
            Total soil moisture in meters (can multiply by 1000 for mm).
        """
        h2osoi = ds[h2osoi_var]  # [time, levgrnd, ncol]
        dzsoi = ds[dzsoi_var]

        if dzsoi.dims != h2osoi.dims:
            dzsoi = dzsoi.broadcast_like(h2osoi)

        total_sm = (h2osoi * dzsoi).sum(dim='levgrnd')
        total_sm.name = 'TotalSoilMoisture'
        total_sm.attrs['units'] = 'm'
        total_sm.attrs['description'] = 'Total soil water content integrated over soil column'

        return total_sm

    @staticmethod
    def compute_swc_top_depth(ds: xr.Dataset, h2osoi_var='H2OSOI', dzsoi_var='DZSOI', depth_threshold=0.05) -> xr.DataArray:
        """
        Compute volumetric soil water content averaged over the top soil layers
        up to a specified depth (e.g., 5 cm = 0.05 m).

        Parameters
        ----------
        ds : xr.Dataset
            Dataset containing H2OSOI (volumetric soil moisture, mm3/mm3) and DZSOI (layer thickness, m).
        h2osoi_var : str
            Name of volumetric soil moisture variable (default 'H2OSOI').
        dzsoi_var : str
            Name of layer thickness variable (default 'DZSOI').
        depth_threshold : float
            Maximum depth (in meters) to average over (default 0.05 for 5 cm).

        Returns
        -------
        swc_top : xr.DataArray
            Depth-weighted average soil water content over the top `depth_threshold` meters [unit: m3/m3].
        """

        h2osoi = ds[h2osoi_var]  # shape: [time, levgrnd, ncol] or [levgrnd, ncol]
        dzsoi = ds[dzsoi_var]    # shape: [levgrnd] or broadcastable to h2osoi

        # Broadcast DZSOI to H2OSOI shape if needed
        if dzsoi.dims != h2osoi.dims:
            dzsoi = dzsoi.broadcast_like(h2osoi)

        # Compute cumulative depth
        depth = dzsoi.cumsum(dim='levgrnd')

        # Mask layers above threshold and compute weighted average
        mask = depth <= depth_threshold
        weighted_swc = (h2osoi * dzsoi).where(mask, 0.0)
        total_depth = dzsoi.where(mask, 0.0).sum(dim='levgrnd')

        swc_top = weighted_swc.sum(dim='levgrnd') / total_depth
        swc_top.name = f'swc_top{int(depth_threshold*100)}cm'
        swc_top.attrs['units'] = 'm3 m-3'
        swc_top.attrs['description'] = f'Depth-weighted volumetric soil moisture in top {int(depth_threshold*100)} cm'

        return swc_top

    @staticmethod
    def _define_region(regnam: str = 'global') -> Tuple[Tuple[float, float], Tuple[float, float]]:
        reg_dict = {
            'global':     [(-90, 90), (-180, 180)],
            'NHMidLat':   [(25, 50), (-60, 150)],
            'Tropics':    [(-10, 10), (-90, 60)],
            'Atlantic':   [(5, 55), (-95, -40)],
            'CONUS':      [(25, 50), (-125, -95)],
            'Antarctic':  [(-90, -50), (-180, 180)],
            'PolarN':     [(50, 90), (-180, 180)],
            'Greenland':  [(60, 85), (-75, -10)],
        }
        if regnam not in reg_dict:
            raise ValueError(f"Region '{regnam}' not defined. Available: {list(reg_dict.keys())}")
        return reg_dict[regnam]

    @staticmethod
    def _extract_exp_info(base_name, data_path) -> dict:
        exps = {
            'CTRLEN10':       {'name': 'ctrl_en10', 'nens': 10, 'period': '201112'},
            'DARTEN10':       {'name': 'dart_en10', 'nens': 10, 'period': '201112'},
            'DARTEN20':       {'name': 'dart_en20', 'nens': 20, 'period': '201112'},
            'DARTEN40':       {'name': 'dart_en40', 'nens': 40, 'period': '201112'},
            'CAPTEN10_15day': {'name': 'capt_en10', 'nens': 10, 'period': '201201-201202'},
            'DARTEN20_15day': {'name': 'dart_en20', 'nens': 20, 'period': '201201-201202'},
            'DARTEN40_15day': {'name': 'dart_en40', 'nens': 40, 'period': '201201-201202'},
        }
        return {
            exp_name: {
                'run': f"{exp_name}_{base_name}" if base_name else exp_name,
                'key': exp_data['name'],
                'nens': exp_data['nens'],
                'atm': 'archive/post/atm/180x360_aave',
                'lnd': 'archive/post/lnd/180x360_aave',
                'period': exp_data['period']
            }
            for exp_name, exp_data in sorted(exps.items())
        }

    def _ensure_datetime64(self, da: xr.DataArray) -> xr.DataArray:
        if not isinstance(da.indexes['time'], pd.DatetimeIndex):
            da['time'] = da.indexes['time'].to_datetimeindex(time_unit='ns')
        return da

    def _convert_model_units(self, var, data):
        if var == 'PRECT':
            return data * 86400.0 * 1000.0  # m/s to mm/day
        elif var in ['TREFHT', 'TS', 'T850']:
            return data - 273.15  # K to C
        elif var in ['PSL', 'PS']:
            return data / 100.0  # Pa to hPa
        return data

    def _extract_year(self, period_str):
        return ''.join(filter(str.isdigit, period_str))[:4]

    def _get_search_pattern(self, var, member, run, freq, period):
        if freq in ['6hourly', 'daily']:
            return f"{var}.{member}.{period}.nc"
        elif freq == 'monthly':
            return f"{run}.{member}.*.{period}.nc"
        else:
            raise ValueError(f"Unsupported frequency: {freq}")

    def read_variable(self, var, exp, frequency=None, regional_mean=False):
        run = self.exp_dict[exp]['run']
        component_path = self.exp_dict[exp][self.component]
        nens = self.exp_dict[exp]['nens']
        if not self.period:
            period_raw = self.exp_dict[exp]['period']
        else:
            period_raw = self.period
        freq = frequency or self.frequency
        period = self._extract_year(period_raw)

        if freq in ['6hourly']:
            prefix = 'ts/6hourly'
        elif freq in ['daily']:
            prefix = 'ts/daily'
        elif freq == 'monthly':
            prefix = 'clim'
            period = period_raw
        else:
            raise ValueError(f"Unsupported frequency: {freq}")

        (lat1, lat2), (lon1, lon2) = self.region
        ensemble_data = []

        for i in range(1, nens + 1):
            member = f"EN{i:02d}"
            search_dir = os.path.join(self.base_path, run, component_path, prefix)
            file_pattern = self._get_search_pattern(var, member, run, freq, period)
            search_path = os.path.join(search_dir, file_pattern)
            rpath = sorted(glob.glob(search_path))

            if not rpath:
                raise FileNotFoundError(f"No files for {member}: {search_path}")

            ds = xcdat.open_mfdataset(
                rpath,
                combine='nested',
                concat_dim='time',
                coords='minimal',
                compat='override',
                parallel=True
            )

            if "lon" in ds and ds.lon.min() >= 0:
                ds = ds.assign_coords(lon=((ds.lon + 180) % 360 - 180)).sortby("lon")

            if var not in ds:
                if var in ['PRECT']:
                    ds[var] = ds['PRECC'] + ds['PRECL']
                elif var in ['PRECST']:
                    ds[var] = ds['PRECSC'] + ds['PRECSL']
                else:
                    raise KeyError(f"Variable '{var}' missing in {rpath[0]}")

            data = self._convert_model_units(var, ds[var])
            data = data.sel(lat=slice(lat1, lat2), lon=slice(lon1, lon2))

            if regional_mean:
                weights = np.cos(np.deg2rad(data.lat))
                weights /= weights.mean()
                data = (data * weights).mean(dim=["lat", "lon"])

            data = self._ensure_datetime64(data)
            ensemble_data.append(data.expand_dims(member=[i]))

        combined = xr.concat(ensemble_data, dim="member")
        return combined.chunk({"member": -1, "time": -1}) if regional_mean else combined.chunk({"member": -1, "time": -1, "lat": -1, "lon": -1})

    def read_variable_across_experiments(self, var: str, model_list: list, **kwargs) -> Dict[str, xr.DataArray]:
        data_dict = {}
        for model_name in model_list:
            da = self.read_variable(var, model_name, **kwargs)
            data_dict[model_name] = da
        return data_dict

class EnsembleMetricEvaluator:
    def __init__(self, obs, model_dict, component, ref_dataset, var, vunit, output_path):
        self.obs_data = obs
        self.model_dict = model_dict
        self.component = component
        self.reference = ref_dataset
        self.var_name = var
        self.var_unit = vunit
        self.output_path = output_path

        self.results = self.derive_metrics_data()

    def _generate_coslat_weight(self, da):
        weights = np.cos(np.deg2rad(da['lat']))
        weights, _ = xr.broadcast(weights, da)
        return weights

    def _get_ensemble_dim(self, da: xr.DataArray) -> str:
        for dim in da.dims:
            if dim in ['ensemble', 'member']:
                return dim
        raise ValueError("No ensemble dimension found (expected 'ensemble' or 'member').")

    def compute_rank_histogram(self, members: xr.DataArray, obs: xr.DataArray, normalize=False):
        ens_sorted = np.sort(members.values, axis=0)
        obs_val = obs.values
        n_ens = ens_sorted.shape[0]

        ranks = np.sum(ens_sorted < obs_val, axis=0).flatten()
        valid_ranks = ranks[(~np.isnan(ranks)) & (ranks >= 0) & (ranks <= n_ens)]
        hist, _ = np.histogram(valid_ranks, bins=np.arange(n_ens + 2) - 0.5, density=normalize)
        return xr.DataArray(hist, dims=['rank'], coords={'rank': np.arange(n_ens + 1)})

    def compute_ensemble_coverage(self, ens_dim: str, members: xr.DataArray, obs: xr.DataArray):
        ens_min = members.min(dim=ens_dim)
        ens_max = members.max(dim=ens_dim)
        within = ((obs >= ens_min) & (obs <= ens_max)).astype(float)

        return xr.DataArray(
            within.mean().compute().item(),
            attrs={
                'units': '1',
                'long_name': 'Fraction of domain where obs within ensemble envelope'
            }
        )

    def compute_crps(self, members: xr.DataArray, obs: xr.DataArray) -> xr.DataArray:
        try:
            ens = members.transpose('lat', 'lon', members.dims[0]).values
            obs_arr = obs.values
            crps = crps_ensemble(obs_arr, ens)
            crps_da = xr.DataArray(crps, dims=['lat', 'lon'], coords={'lat': obs.lat, 'lon': obs.lon})
            return xr.where(np.isfinite(obs), crps_da, np.nan)
        except Exception as e:
            print(f"[WARN] CRPS computation failed: {e}")
            return xr.full_like(obs, np.nan)

    def _derive_hor_metrics_data(self, dso1, dsm1):
        mec_dict = {}
        dat_mean, dat_std = None, None

        dsm_mean = dsm1.mean(dim='time') if 'time' in dsm1.dims else dsm1
        weights = self._generate_coslat_weight(dsm_mean)

        mec_dict['mean'] = (dsm_mean * weights).sum(dim=["lat", "lon"]) / weights.sum(dim=["lat", "lon"])
        mec_dict['rmsm'] = np.sqrt(((dsm_mean ** 2) * weights).sum(dim=["lat", "lon"]) / weights.sum(dim=["lat", "lon"]))

        if isinstance(dso1, xr.DataArray):
            if 'time' in dsm1.dims:
                obs_broadcast = dso1.expand_dims({'time': dsm1.time}, axis=0)
            else:
                obs_broadcast = dso1

            dsm1 = dsm1.chunk(dict(lat=-1, lon=-1))
            obs_broadcast = obs_broadcast.chunk(dict(lat=-1, lon=-1))

            mec_dict['rmso'] = np.sqrt(((dso1 ** 2) * weights).sum(dim=["lat", "lon"]) / weights.sum(dim=["lat", "lon"]))
            mec_dict['rmse'] = xs.rmse(obs_broadcast, dsm1, dim=["lat", "lon"], weights=weights, skipna=True)
            mec_dict['pcor'] = xs.pearson_r(obs_broadcast, dsm1, dim=["lat", "lon"], weights=weights, skipna=True)

            diff = dsm1 - obs_broadcast
            dat_mean = diff.mean(dim='time') if 'time' in diff.dims else diff
            dat_std = diff.std(dim='time') if 'time' in diff.dims else xr.zeros_like(diff)

        return mec_dict, dat_mean, dat_std

    def derive_metrics_data(self):
        results = {}
        for model in self.model_dict:
            print(f"[INFO] Processing {model}")

            ens_dim = self._get_ensemble_dim(self.model_dict[model])
            mod_full = self.model_dict[model]

            if 'time' in mod_full.dims and mod_full.sizes['time'] == 1:
                mod_full = mod_full.squeeze('time')

            ens_mean = mod_full.mean(dim=ens_dim)
            ens_std = mod_full.std(dim=ens_dim)
            metrics, mean_bias, std_bias = self._derive_hor_metrics_data(self.obs_data, ens_mean)

            spread = ens_std.std(dim='time') if 'time' in ens_std.dims else ens_std
            spread_mean = spread.mean(dim=["lat", "lon"])
            rmse_mean = metrics['rmse'].mean()
            bias_mean = mean_bias.mean(dim=["lat", "lon"])
            mse = rmse_mean ** 2
            bias_sq = bias_mean ** 2
            spread_sq = spread_mean ** 2
            sem = spread_mean / np.sqrt(mod_full.sizes[ens_dim])

            metrics.update({
                'spread': spread_mean,
                'spread_skill_ratio': spread_mean / rmse_mean,
                'bias_squared': bias_sq,
                'spread_squared': spread_sq,
                'mse': mse,
                'sem': sem
            })

            rank_hist = self.compute_rank_histogram(mod_full, self.obs_data)
            coverage = self.compute_ensemble_coverage(ens_dim, mod_full, self.obs_data)
            crps_map = self.compute_crps(mod_full, self.obs_data)

            metrics.update({
                'rank_histogram': rank_hist,
                'coverage': coverage,
                'crps': crps_map,
                'crps_mean': crps_map.mean(dim=['lat', 'lon'])
            })

            results[model] = {
                'metrics': metrics,
                'mod': mod_full,
                'mean_map': ens_mean,
                'spread_map': ens_std,
                'obs': self.obs_data,
                'bias': mean_bias,
                'stddev': std_bias
            }

            filepath_nc = os.path.join(self.output_path, f"{model}_ensemble_mean_bias.nc")
            self.save_to_netcdf(filepath_nc, model, results[model])
            filepath_csv = os.path.join(self.output_path, f"{model}_ensemble_bias_summary.csv")
            self.save_summary_csv(filepath_csv, model, results[model])

        return results

    def save_to_netcdf(self, filepath, model, model_data):
        ds_out = {}
        metadata = {
            'mean': {'units': self.var_unit, 'long_name': f'{self.var_name} area-weighted mean (model)'},
            'rmsm': {'units': self.var_unit, 'long_name': f'{self.var_name} RMS (model)'},
            'rmso': {'units': self.var_unit, 'long_name': f'{self.var_name} RMS (observation)'},
            'rmse': {'units': self.var_unit, 'long_name': f'{self.var_name} RMSE (model vs {self.reference})'},
            'pcor': {'units': '1', 'long_name': f'{self.var_name} pattern correlation (model vs {self.reference})'},
            'bias': {'units': self.var_unit, 'long_name': f'{self.var_name} mean bias (model - {self.reference})'},
            'stddev': {'units': self.var_unit, 'long_name': f'{self.var_name} stddev of bias over time'},
            'spread_map': {'units': self.var_unit, 'long_name': f'{self.var_name} ensemble spread'},
            'mean_map': {'units': self.var_unit, 'long_name': f'{self.var_name} ensemble mean'},
            'spread': {'units': self.var_unit, 'long_name': f'{self.var_name} spatial mean ensemble spread'},
            'spread_skill_ratio': {'units': '1', 'long_name': 'Spread-to-RMSE ratio'},
            'bias_squared': {'units': f'{self.var_unit}^2', 'long_name': 'Squared ensemble bias'},
            'spread_squared': {'units': f'{self.var_unit}^2', 'long_name': 'Squared ensemble spread'},
            'mse': {'units': f'{self.var_unit}^2', 'long_name': 'Mean square error (bias² + spread²)'},
            'sem': {'units': self.var_unit, 'long_name': 'Standard error of the ensemble mean'},
            'rank_histogram': {'units': 'count', 'long_name': 'Rank histogram of observation within ensemble'},
            'coverage': {'units': '1', 'long_name': 'Fraction of domain where obs within ensemble envelope'},
            'crps': {'units': self.var_unit, 'long_name': f'CRPS: Continuous Ranked Probability Score for {self.var_name}'},
            'crps_mean': {'units': self.var_unit, 'long_name': f'Mean CRPS over spatial domain for {self.var_name}'}
        }

        for key, val in model_data['metrics'].items():
            da = val if isinstance(val, xr.DataArray) else xr.DataArray(val)

            # Clean dimension conflicts
            if 'time' in da.coords and 'time' not in da.dims:
                da = da.reset_coords('time', drop=True)

            if da.ndim == 0:
                da = da.expand_dims({'scalar_dim': [0]})

            if key in metadata:
                da.attrs.update(metadata[key])
            ds_out[key] = da

        for key in ['bias', 'stddev', 'mean_map', 'spread_map']:
            val = model_data.get(key)
            if isinstance(val, xr.DataArray):
                if 'time' in val.coords and 'time' not in val.dims:
                    val = val.reset_coords('time', drop=True)
                val.attrs.update(metadata.get(key, {}))
                ds_out[key] = val

        ds = xr.Dataset(ds_out)
        ds.attrs.update({
            'model': model,
            'component': self.component,
            'reference_dataset': self.reference,
            'variable': self.var_name,
            'units': self.var_unit,
            'created_by': 'EnsembleMetricEvaluator',
            'creation_time': datetime.now().isoformat()
        })

        ds.to_netcdf(filepath)
        print(f"[INFO] Saved metrics for {model} to {filepath}")

    def safe_mean(self, val):
        if isinstance(val, xr.DataArray):
            try:
                return float(val.mean().values)
            except Exception:
                return np.nan
        elif isinstance(val, (int, float, np.number)):
            return float(val)
        else:
            return np.nan

    def save_summary_csv(self, filepath, model, model_data):
        summary = []
        metrics = model_data['metrics']
        row = {
            'model': model,
            'RMSE': self.safe_mean(metrics.get('rmse')),
            'PCOR': self.safe_mean(metrics.get('pcor')),
            'Spread': self.safe_mean(metrics.get('spread')),
            'SSR': self.safe_mean(metrics.get('spread_skill_ratio')),
            'Bias^2': self.safe_mean(metrics.get('bias_squared')),
            'Spread^2': self.safe_mean(metrics.get('spread_squared')),
            'MSE': self.safe_mean(metrics.get('mse')),
            'SEM': self.safe_mean(metrics.get('sem')),
            'Coverage': self.safe_mean(metrics.get('coverage')),
            'CRPS': self.safe_mean(metrics.get('crps_mean'))
        }

        if all(pd.isna(v) for k, v in row.items() if k != 'model'):
            print(f"[WARN] No valid metrics found for {model}, CSV will be empty.")
            return

        summary.append(row)
        df = pd.DataFrame(summary)
        if os.path.exists(filepath):
            os.remove(filepath)
        df.to_csv(filepath, index=False)
        print(f"[INFO] Saved summary CSV to {filepath}")
