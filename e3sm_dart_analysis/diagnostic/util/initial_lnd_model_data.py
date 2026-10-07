import os
import glob
import numpy as np
import xarray as xr

from xcdat.dataset import open_dataset
from xcdat.bounds import create_bounds
from xcdat.dataset import open_mfdataset

import pandas as pd
import xskillscore as xs
from typing import Dict, Tuple
from datetime import datetime

from util.array_stats import nanquantile_numpy

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

    def compute_total_soil_moisture(self, ds, h2osoi_var='H2OSOI', dzsoi_var='DZSOI') -> xr.DataArray:
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

    def compute_swc_top_depth(self, ds: xr.Dataset, h2osoi_var='H2OSOI', dzsoi_var='DZSOI', depth_threshold=0.05) -> xr.DataArray:
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
            'CTRLEN10':       {'name': 'ctrl_en10', 'nens': 10, 'period': '201112-201112'},
            'DARTEN10':       {'name': 'dart_en10', 'nens': 10, 'period': '201112-201112'},
            'DARTEN20':       {'name': 'dart_en20', 'nens': 20, 'period': '201112-201112'},
            'DARTEN40':       {'name': 'dart_en40', 'nens': 40, 'period': '201112-201112'},
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

    def _ensure_datetime64(self, da: xr.DataArray, frequency: str) -> xr.DataArray:
        """Ensure 'time' coordinate is datetime64[ns] and normalize if daily or coarser."""
        # Convert CFTimeIndex or any non-datetime64 to datetime64[ns]
        if not isinstance(da.indexes['time'], pd.DatetimeIndex):
            da['time'] = da.indexes['time'].to_datetimeindex(time_unit='ns')

        # Normalize time for daily and coarser frequencies
        if frequency in ['daily', 'day', 'D', 'monthly', 'mon', 'M', 'yearly', 'annual', 'Y']:
            # Normalize to midnight to remove any hour/minute/second inconsistencies
            da['time'] = pd.to_datetime(da['time'].values).normalize()

        return da

    def _convert_model_units(self, var, data):
        if var == 'PRECT':
            return data * 86400.0 * 1000.0  # m/s to mm/day
        elif var in ['PSL', 'PS']:
            return data / 100.0  # Pa to hPa
        return data

    def _extract_year(self, period_str):
        return ''.join(filter(str.isdigit, period_str))[:4]

    def _get_search_pattern(self, var, member, run, freq):
        if freq in ['6hourly', 'daily']:
            return f"{var}.{member}.*.nc"
        elif freq == 'monthly':
            return f"{run}.{member}.*.*.nc"
        else:
            raise ValueError(f"Unsupported frequency: {freq}")

    def parse_time_range(self,time):
        # Check if the input is a time range string like '201101-201202'
        if isinstance(time, str) and '-' in time:
            start_date, end_date = time.split('-')
            start = pd.to_datetime(start_date)
            #end = pd.to_datetime(end_date) + pd.offsets.MonthEnd(0)
            end = pd.to_datetime(end_date)
            end = end + pd.Timedelta(days=1) - pd.Timedelta(seconds=1)  # Make the end time inclusive

            start_year, start_month = int(start_date[:4]), int(start_date[4:])
            end_year, end_month = int(end_date[:4]), int(end_date[4:])

            # Create a date range from the start month to the end month
            time_range = slice(start, end)
        elif isinstance(time, slice):
            time_range = time
            if time.start is None or time.stop is None:
                raise ValueError("Both start and end must be specified in the slice.")
            start = pd.to_datetime(time.start)
            end = pd.to_datetime(time.stop)
        else:
            raise TypeError("Only string ('YYYYMM-YYYYMM') or slice inputs are supported.")

        # Extract the years from the time range
        years = list(range(start.year, end.year + 1))

        return time_range, years

    def read_variable(self, var, exp, frequency=None, regional_mean=False):
        run = self.exp_dict[exp]['run']
        component_path = self.exp_dict[exp][self.component]
        nens = self.exp_dict[exp]['nens']

        # identify time range and list of years
        if not self.period:
            time_range, years = self.parse_time_range(self.exp_dict[exp]['period'])
        else:
            time_range, years = self.parse_time_range(self.period)

        freq = frequency or self.frequency

        if freq in ['6hourly']:
            prefix = 'ts/6hourly'
        elif freq in ['daily']:
            prefix = 'ts/daily'
        elif freq == 'monthly':
            prefix = 'clim'
        else:
            raise ValueError(f"Unsupported frequency: {freq}")

        (lat1, lat2), (lon1, lon2) = self.region
        ensemble_data = []

        for i in range(1, nens + 1):
            member = f"EN{i:02d}"
            search_dir = os.path.join(self.base_path, run, component_path, prefix)
            file_pattern = self._get_search_pattern(var, member, run, freq)
            search_path = os.path.join(search_dir, file_pattern)
            rpath = sorted(glob.glob(search_path))

            if not rpath:
                raise FileNotFoundError(f"No files for {member}: {search_path}")

            ds = open_mfdataset(
                rpath,
                combine='nested',
                concat_dim='time',
                coords='minimal',
                compat='override',
                parallel=True
            )

            # Normalize the time to ignore the hour, minute, and second part
            ds = self._ensure_datetime64(ds,freq)

            if "lon" in ds and ds.lon.min() >= 0:
                ds = ds.assign_coords(lon=((ds.lon + 180) % 360 - 180)).sortby("lon")

            if var in ['H2OSOI']:
                ds[var] = self.compute_total_soil_moisture(ds)
            elif var not in ds:
                if var in ['PRECT']:
                    ds[var] = ds['PRECC'] + ds['PRECL']
                elif var in ['PRECST']:
                    ds[var] = ds['PRECSC'] + ds['PRECSL']
                else:
                    raise KeyError(f"Variable '{var}' missing in {rpath[0]}")

            data = self._convert_model_units(var, ds[var])
            data = data.sel(lat=slice(lat1, lat2), lon=slice(lon1, lon2), time=time_range)

            if regional_mean:
                weights = np.cos(np.deg2rad(data.lat))
                weights /= weights.mean()
                data = (data * weights).mean(dim=["lat", "lon"])

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

    def _generate_coslat_weight(self, da: xr.DataArray) -> xr.DataArray:
        """
        Generate cosine-latitude weights broadcasted to the shape of the input DataArray.

        Parameters
        ----------
        da : xr.DataArray
            Input data array with 'lat' dimension or coordinate.

        Returns
        -------
        xr.DataArray
            Cosine-latitude weights broadcast to the shape of da.
        """
        if 'lat' not in da.coords:
            raise ValueError("Latitude coordinate 'lat' not found in DataArray.")

        weights = np.cos(np.deg2rad(da['lat']))

        # If 'lat' is not a dimension but just a coordinate, make it a dimension
        if 'lat' not in da.dims:
            weights = weights.expand_dims(dim='lat')

        # Broadcast to match da shape
        weights, _ = xr.broadcast(weights, da)

        return weights

    def _get_ensemble_dim(self, da: xr.DataArray) -> str:
        for dim in da.dims:
            if dim in ['ensemble', 'member']:
                return dim
        raise ValueError("No ensemble dimension found (expected 'ensemble' or 'member').")

    def bootstrap_ci_map(self, model, obs, metric_fn, ens_dim = "ensemble", n_boot=1000, alpha=0.05):
        """
        Compute metric with bootstrap significance test.

        Parameters:
            model (xr.DataArray): [ensemble, lat, lon]
            obs (xr.DataArray): [lat, lon]
            metric_fn (callable): Function(model_subset, obs) -> xr.DataArray[lat, lon]
            n_boot (int): Number of bootstrap resamples
            alpha (float): Significance level

        Returns:
            metric (xr.DataArray): Mean metric
            mask (xr.DataArray): Boolean mask where metric is significant
            ci (tuple): Lower and upper confidence interval (xr.DataArray)
        """
        bootstrap_metrics = []

        for i in range(n_boot):
            # Resample with replacement along the ensemble dimension
            resample = model.sel({ens_dim: np.random.choice(model[ens_dim].values, size=model.sizes[ens_dim], replace=True)})
            metric = metric_fn(resample, obs, ens_dim = ens_dim)
            bootstrap_metrics.append(metric)

        boot_maps = np.stack([m.values for m in bootstrap_metrics], axis=0)
        lower = nanquantile_numpy(boot_maps, alpha / 2, axis=0)
        upper = nanquantile_numpy(boot_maps, 1 - alpha / 2, axis=0)
        mean_metric = metric_fn(model, obs, ens_dim = ens_dim)
        mask = ((mean_metric < lower) | (mean_metric > upper)).astype(int)

        lower_da = xr.DataArray(lower, dims=['lat', 'lon'], coords={'lat': obs.lat, 'lon': obs.lon})
        upper_da = xr.DataArray(upper, dims=['lat', 'lon'], coords={'lat': obs.lat, 'lon': obs.lon})
        signif_mask = xr.DataArray(mask, dims=['lat', 'lon'], coords={'lat': obs.lat, 'lon': obs.lon})

        return lower_da, upper_da, signif_mask

    def bias_map_fn(self, model, obs, ens_dim="member"):
        # Drop ensemble dim from obs if accidentally present (e.g., from broadcasting)
        if ens_dim in obs.dims:
            obs = obs.mean(dim=ens_dim, skipna=True)  # Or .isel({ens_dim: 0}) if you expect one value
        model = model.transpose(ens_dim,'lat','lon')
        obs = obs.transpose('lat','lon')
        return model.mean(dim=ens_dim, skipna=True) - obs

    def rmse_map_fn(self, model, obs, ens_dim="member"):
        model = model.transpose(ens_dim, "lat", "lon")
        # Drop ensemble dim from obs if accidentally present (e.g., from broadcasting)
        if ens_dim in obs.dims:
            obs = obs.mean(dim=ens_dim, skipna=True)  # Or .isel({ens_dim: 0}) if you expect one value
        # Broadcast obs to model's shape
        obs = obs.expand_dims({ens_dim: model[ens_dim]})
        obs = obs.transpose(ens_dim, "lat", "lon")
        return xs.rmse(model, obs, dim=ens_dim, skipna=True)

    def spread_map_fn(self, model, obs=None, ens_dim='ensemble'):
        model = model.transpose(ens_dim,'lat','lon')
        return model.std(dim=ens_dim, skipna=True)

    def crps_map_fn(self, model, obs, ens_dim='ensemble'):
        # If obs accidentally has an ensemble dimension, remove it
        if ens_dim in obs.dims:
            obs = obs.mean(dim=ens_dim, skipna=True)

        # Ensure both are lat-lon aligned
        model = model.transpose(ens_dim, "lat", "lon")
        obs = obs.transpose("lat", "lon")
        # Compute CRPS (no dim specified; shape is [lat, lon])
        crps = xs.crps_ensemble(obs, model, member_dim=ens_dim, dim=[])
        return crps.where(obs.notnull())

    def compute_ci_and_significance(
           self, metric_name, metric_fn, members, obs, ens_dim='ensemble',
           alpha=0.05, n_boot=1000
        ):
        """
        Compute 2D CI and significance mask for a given metric map (e.g., bias, RMSE).
        Returns (lower_CI, upper_CI, significance_mask)
        """
        try:
            ci_lower, ci_upper, signif_mask = self.bootstrap_ci_map(
                members, obs, metric_fn, ens_dim, n_boot, alpha
            )
            signif_mask.name = f'{metric_name}_significance_mask'
            ci_lower.name = f'{metric_name}_ci_lower'
            ci_upper.name = f'{metric_name}_ci_upper'

            signif_mask.attrs.update({'units': '1', 'long_name': f'Significance mask for {metric_name.replace("_", " ")}'})
            ci_lower.attrs.update({'units': self.var_unit, 'long_name': f'Lower bound of 95% CI for {metric_name.replace("_", " ")}'})
            ci_upper.attrs.update({'units': self.var_unit, 'long_name': f'Upper bound of 95% CI for {metric_name.replace("_", " ")}'})

            return ci_lower, ci_upper, signif_mask
        except Exception as e:
            print(f"[WARN] Failed CI + significance for {metric_name}: {e}")
            return None, None, None

    def bootstrap_confidence_interval(self, metric_fn, members, obs, ens_dim='ensemble', n_boot=1000, alpha=0.05):
        """Compute bootstrap confidence interval for an ensemble-based metric."""
        rng = np.random.default_rng(seed=42)
        n_ens = members.sizes[ens_dim]
        boot_vals = []

        for _ in range(n_boot):
            sample_idx = rng.integers(0, n_ens, size=n_ens)
            sample = members.isel({ens_dim: sample_idx})
            try:
                val = metric_fn(sample.mean(dim=ens_dim), obs, ens_dim = ens_dim)
                if np.isfinite(val):
                    boot_vals.append(val)
            except Exception:
                continue

        if not boot_vals:
            return np.nan, np.nan, np.nan

        boot_vals = np.array(boot_vals)
        lower = np.percentile(boot_vals, 100 * alpha / 2)
        upper = np.percentile(boot_vals, 100 * (1 - alpha / 2))
        return lower, upper, np.mean(boot_vals)

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
            ens_dim = self._get_ensemble_dim(members)
            # Ensure dims are ordered correctly and data types are preserved
            members = members.transpose(ens_dim, "lat", "lon")
            obs = obs.transpose("lat", "lon")
            # Compute CRPS using xskillscore (returns xarray.DataArray)
            crps = xs.crps_ensemble(obs, members, member_dim=ens_dim, dim=[])
            return crps.where(obs.notnull())
        except Exception as e:
            print(f"[WARN] CRPS computation failed: {e}")
            return xr.full_like(obs, np.nan)

    def _derive_hor_metrics_data(self, dso1, dsm1, ens_dim):
        # Chunk data for performance (if using dask-backed xarray)
        dsm1 = dsm1.chunk({dim: -1 for dim in dsm1.dims if dim in ['time', ens_dim, 'lat', 'lon']})
        dso1 = dso1.chunk({dim: -1 for dim in dso1.dims if dim in ['time', 'lat', 'lon']})

        mec_dict = {}
        bias, stddev, rmse, spread, dsm_mean = None, None, None, None, None

        ens_mean = dsm1.mean(dim=ens_dim, skipna=True)
        ens_std = dsm1.std(dim=ens_dim, skipna=True)

        dsm_mean = ens_mean.mean(dim='time', skipna=True) if 'time' in ens_mean.dims else ens_mean
        weights = self._generate_coslat_weight(dsm_mean)

        # Model global mean and RMS
        mec_dict['mod_mean'] = (dsm_mean * weights).sum(dim=["lat", "lon"], skipna=True) / weights.sum(dim=["lat", "lon"], skipna=True)
        mec_dict['mod_rms'] = np.sqrt((((dsm_mean - mec_dict['mod_mean']) ** 2) * weights).sum(dim=["lat", "lon"], skipna=True) / weights.sum(dim=["lat", "lon"], skipna=True))

        if isinstance(dso1, xr.DataArray):
            dso_mean = dso1.mean(dim='time', skipna=True) if 'time' in dso1.dims else dso1

            mec_dict['obs_mean'] = (dso_mean * weights).sum(dim=["lat", "lon"], skipna=True) / weights.sum(dim=["lat", "lon"], skipna=True)
            mec_dict['obs_rms'] = np.sqrt((((dso_mean - mec_dict['obs_mean']) ** 2) * weights).sum(dim=["lat", "lon"], skipna=True) / weights.sum(dim=["lat", "lon"], skipna=True))

            # Global RMSE and pattern correlation
            mec_dict['rmse_glb'] = xs.rmse(dso_mean, dsm_mean, dim=["lat", "lon"], weights=weights, skipna=True)
            mec_dict['pcor'] = xs.pearson_r(
                dso_mean - mec_dict['obs_mean'],
                dsm_mean - mec_dict['mod_mean'],
                dim=["lat", "lon"],
                weights=weights,
                skipna=True
            )

            # Bias (model mean - obs) averaged over time and space
            diff = ens_mean - dso1
            bias = diff.mean(dim='time', skipna=True) if 'time' in diff.dims else diff
            mec_dict['bias'] = (bias * weights).sum(dim=["lat", "lon"], skipna=True) / weights.sum(dim=["lat", "lon"], skipna=True)

            # Std dev of difference (observation error) and spread of ensemble
            stddev = diff.std(dim='time', skipna=True) if 'time' in diff.dims else xr.zeros_like(diff)
            spread = ens_std.mean(dim='time', skipna=True) if 'time' in ens_std.dims else ens_std
            mec_dict['spread'] = (spread * weights).sum(dim=["lat", "lon"], skipna=True) / weights.sum(dim=["lat", "lon"], skipna=True)

            # RMSE per ensemble member and spatial average
            dso1_expanded = dso1.expand_dims({ens_dim: dsm1[ens_dim]})
            rerr = xs.rmse(dsm1, dso1_expanded, dim=ens_dim, skipna=True)
            rmse = rerr.mean(dim='time', skipna=True) if 'time' in rerr.dims else rerr
            mec_dict['rmse'] = (rmse * weights).sum(dim=["lat", "lon"], skipna=True) / weights.sum(dim=["lat", "lon"], skipna=True)

        return mec_dict, bias, stddev, rmse, spread, dsm_mean

    def derive_metrics_data(self):
        results = {}
        for model in self.model_dict:
            print(f"[INFO] Processing {model}")

            ens_dim = self._get_ensemble_dim(self.model_dict[model])
            mod_full = self.model_dict[model]
            obs_full = self.obs_data

            if 'time' in mod_full.dims and mod_full.sizes['time'] == 1:
                mod_full = mod_full.squeeze('time')

            if 'time' in obs_full.dims and obs_full.sizes['time'] == 1:
                obs_full = obs_full.squeeze('time')

            metrics, bias_map, bias_std_map, rmse_map, spread_map, mean_map = self._derive_hor_metrics_data(
                obs_full, mod_full, ens_dim
            )

            mse = metrics['rmse'] ** 2
            bias_sq = metrics['bias'] ** 2
            spread_sq = metrics['spread'] ** 2
            sem = metrics['spread'] / np.sqrt(mod_full.sizes[ens_dim])
            metrics.update({
                'ssr': metrics['spread'] / metrics['rmse'],
                'bias_squared': bias_sq,
                'spread_squared': spread_sq,
                'mse': mse,
                'sem': sem
            })

            rank_hist = self.compute_rank_histogram(mod_full, obs_full)
            coverage = self.compute_ensemble_coverage(ens_dim, mod_full, obs_full)
            crps_map = self.compute_crps(mod_full, obs_full)
            weights = self._generate_coslat_weight(crps_map)
            crps_mean = (crps_map * weights).sum(dim=["lat", "lon"], skipna=True) / weights.sum(dim=["lat", "lon"], skipna=True)
            metrics.update({
                'rank_histogram': rank_hist,
                'coverage': coverage,
                'crps': crps_mean
            })

            # Bootstrap significance testing
            bias_fn = lambda mod, obs: float((mod - obs).mean(dim=["lat", "lon"], skipna=True).values)
            rmse_fn = lambda mod, obs: float(xs.rmse(obs, mod, dim=["lat", "lon"], skipna=True).values)
            spread_fn = lambda mod, obs: float(mod.std(dim=ens_dim, skipna=True).mean(dim=["lat", "lon"], skipna=True).values)
            pcor_fn = lambda mod, obs: float(xs.pearson_r(obs, mod, dim=["lat", "lon"], skipna=True).values)
            ssr_fn = lambda mod, obs: spread_fn(mod, obs) / rmse_fn(mod, obs)

            try:
                bias_lb, bias_ub, _ = self.bootstrap_confidence_interval(bias_fn, mod_full, self.obs_data, ens_dim)
                rmse_lb, rmse_ub, _ = self.bootstrap_confidence_interval(rmse_fn, mod_full, self.obs_data, ens_dim)
                spread_lb, spread_ub, _ = self.bootstrap_confidence_interval(spread_fn, mod_full, self.obs_data, ens_dim)
                pcor_lb, pcor_ub, _ = self.bootstrap_confidence_interval(pcor_fn, mod_full, self.obs_data, ens_dim)
                ssr_lb, ssr_ub, _ = self.bootstrap_confidence_interval(ssr_fn, mod_full, self.obs_data, ens_dim)

                metrics.update({
                    'bias_ci_lower': bias_lb,
                    'bias_ci_upper': bias_ub,
                    'rmse_ci_lower': rmse_lb,
                    'rmse_ci_upper': rmse_ub,
                    'spread_ci_lower': spread_lb,
                    'spread_ci_upper': spread_ub,
                    'pcor_ci_lower': pcor_lb,
                    'pcor_ci_upper': pcor_ub,
                    'ssr_ci_lower': ssr_lb,
                    'ssr_ci_upper': ssr_ub
                })
            except Exception as e:
                print(f"[WARN] Bootstrap CI failed: {e}")

            # Bias Map
            bias_ci_lo, bias_ci_hi, bias_ci_mask = self.compute_ci_and_significance(
                "bias_map", self.bias_map_fn, mod_full, obs_full, ens_dim)
            # RMSE Map
            rmse_ci_lo, rmse_ci_hi, rmse_ci_mask = self.compute_ci_and_significance(
                "rmse_map", self.rmse_map_fn, mod_full, obs_full, ens_dim)
            # Spread Map
            spread_ci_lo, spread_ci_hi, spread_ci_mask = self.compute_ci_and_significance(
                "spread_map", self.spread_map_fn, mod_full, obs_full, ens_dim)
            # CRPS Map
            crps_ci_lo, crps_ci_hi, crps_ci_mask = self.compute_ci_and_significance(
                "crps_map", self.crps_map_fn, mod_full, obs_full, ens_dim)

            results[model] = {
                'metrics': metrics,
                'mod_map': mod_full.mean(dim='time', skipna=True) if 'time' in mod_full.dims else mod_full,
                'obs_map': obs_full.mean(dim='time', skipna=True) if 'time' in obs_full.dims else obs_full,
                'mean_map': mean_map,
                'spread_map': spread_map,
                'bias_map': bias_map,
                'bias_stddev_map': bias_std_map,
                'rmse_map': rmse_map,
                'crps_map': crps_map,
                'bias_map_ci_lower': bias_ci_lo,
                'bias_map_ci_upper': bias_ci_hi,
                'bias_map_significance_mask': bias_ci_mask,
                'rmse_map_ci_lower': rmse_ci_lo,
                'rmse_map_ci_upper': rmse_ci_hi,
                'rmse_map_significance_mask': rmse_ci_mask,
                'spread_map_ci_lower': spread_ci_lo,
                'spread_map_ci_upper': spread_ci_hi,
                'spread_map_significance_mask': spread_ci_mask,
                'crps_map_ci_lower': crps_ci_lo,
                'crps_map_ci_upper': crps_ci_hi,
                'crps_map_significance_mask': crps_ci_mask
            }

            filepath_nc = os.path.join(
                self.output_path,
                f"{model}_{self.var_name}_{self.ref_dataset}_ensemble_mean_bias.nc"
            )
            self.save_to_netcdf(filepath_nc, model, results[model])

            filepath_csv = os.path.join(
                self.output_path,
                f"{model}_{self.var_name}_{self.ref_dataset}_ensemble_bias_summary.csv"
            )
            self.save_summary_csv(filepath_csv, model, results[model])

        return results

    def save_to_netcdf(self, filepath, model, model_data):
        ds_out = {}
        metadata = {
            'mod_map': {'units': self.var_unit, 'long_name': f'{self.var_name} ensemble member (model)'},
            'obs_map': {'units': self.var_unit, 'long_name': f'{self.var_name} ensemble member ({self.reference})'},
            'pcor_map': {'units': '1', 'long_name': f'{self.var_name} pattern correlation (model vs {self.reference})'},
            'rmse_map': {'units': self.var_unit, 'long_name': f'{self.var_name} ensemble rmse (model vs {self.reference})'},
            'bias_map': {'units': self.var_unit, 'long_name': f'{self.var_name} ensemble mean bias (model - {self.reference})'},
            'spread_map': {'units': self.var_unit, 'long_name': f'{self.var_name} ensemble spread'},
            'mean_map': {'units': self.var_unit, 'long_name': f'{self.var_name} ensemble mean'},
            'bias_stddev_map': {'units': self.var_unit, 'long_name': f'{self.var_name} stddev of bias over time'},
            'bias_map_ci_lower': {'units': self.var_unit, 'long_name': 'Lower bound of CI for mean bias'},
            'bias_map_ci_upper': {'units': self.var_unit, 'long_name': 'Upper bound of CI for mean bias'},
            'bias_map_significance_mask': {'units': '1', 'long_name': 'Significance mask for mean bias'},
            'rmse_map_ci_lower': {'units': self.var_unit, 'long_name': 'Lower bound of CI for RMSE map'},
            'rmse_map_ci_upper': {'units': self.var_unit, 'long_name': 'Upper bound of CI for RMSE map'},
            'rmse_map_significance_mask': {'units': '1', 'long_name': 'Significance mask for RMSE'},
            'crps_map_ci_lower': {'units': self.var_unit, 'long_name': f'Lower bound of 95% CI for CRPS'},
            'crps_map_ci_upper': {'units': self.var_unit, 'long_name': f'Upper bound of 95% CI for CRPS'},
            'crps_map_significance_mask': {'units': '1', 'long_name': 'Significance mask for CRPS'},
            'spread_map_ci_lower': {'units': self.var_unit, 'long_name': f'Lower bound of 95% CI for CRPS'},
            'spread_map_ci_upper': {'units': self.var_unit, 'long_name': f'Upper bound of 95% CI for CRPS'},
            'spread_map_significance_mask': {'units': '1', 'long_name': 'Significance mask for ensemble spread'},
            'mean': {'units': self.var_unit, 'long_name': f'{self.var_name} area-weighted mean (model)'},
            'rmsm': {'units': self.var_unit, 'long_name': f'{self.var_name} global RMS (model)'},
            'rmso': {'units': self.var_unit, 'long_name': f'{self.var_name} global RMS (observation)'},
            'rmse': {'units': self.var_unit, 'long_name': f'{self.var_name} global mean ensemble RMSE (model vs {self.reference})'},
            'rmse_glb': {'units': self.var_unit, 'long_name': f'{self.var_name} ensemble mean global RMSE (model vs {self.reference})'},
            'spread': {'units': self.var_unit, 'long_name': f'{self.var_name} spatial mean ensemble spread'},
            'ssr': {'units': '1', 'long_name': 'Spread-to-RMSE ratio'},
            'bias': {'units': f'{self.var_unit}^2', 'long_name': 'Spatial mean ensemble mean bias'},
            'bias_squared': {'units': f'{self.var_unit}^2', 'long_name': 'Squared ensemble mean bias'},
            'spread_squared': {'units': f'{self.var_unit}^2', 'long_name': 'Squared ensemble spread'},
            'mse': {'units': f'{self.var_unit}^2', 'long_name': 'Mean square error (bias² + spread²)'},
            'sem': {'units': self.var_unit, 'long_name': 'Standard error of the ensemble mean'},
            'bias_ci_lower': {'units': self.var_unit, 'long_name': 'Lower bound of 95% CI for bias'},
            'bias_ci_upper': {'units': self.var_unit, 'long_name': 'Upper bound of 95% CI for bias'},
            'rmse_ci_lower': {'units': self.var_unit, 'long_name': 'Lower bound of 95% CI for RMSE'},
            'rmse_ci_upper': {'units': self.var_unit, 'long_name': 'Upper bound of 95% CI for RMSE'},
            'spread_ci_lower': {'units': self.var_unit, 'long_name': 'Lower bound of 95% CI for spread'},
            'spread_ci_upper': {'units': self.var_unit, 'long_name': 'Upper bound of 95% CI for spread'},
            'pcor_ci_lower': {'units': '1', 'long_name': 'Lower bound of 95% CI for pattern correlation'},
            'pcor_ci_upper': {'units': '1', 'long_name': 'Upper bound of 95% CI for pattern correlation'},
            'ssr_ci_lower': {'units': '1', 'long_name': 'Lower bound of 95% CI for spread-to-RMSE ratio'},
            'ssr_ci_upper': {'units': '1', 'long_name': 'Upper bound of 95% CI for spread-to-RMSE ratio'},
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

        for key, val in model_data.items():
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
                mean_val = val.mean(skipna=True).values
                return float(mean_val) if not np.isnan(mean_val) else np.nan
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
            'RMSE_GLB': self.safe_mean(metrics.get('rmse_glb')),
            'PCOR': self.safe_mean(metrics.get('pcor')),
            'Spread': self.safe_mean(metrics.get('spread')),
            'SSR': self.safe_mean(metrics.get('ssr')),
            'Bias^2': self.safe_mean(metrics.get('bias_squared')),
            'Spread^2': self.safe_mean(metrics.get('spread_squared')),
            'MSE': self.safe_mean(metrics.get('mse')),
            'SEM': self.safe_mean(metrics.get('sem')),
            'Coverage': self.safe_mean(metrics.get('coverage')),
            'CRPS': self.safe_mean(metrics.get('crps_mean')),
            'Bias_CI_Lower': metrics.get('bias_ci_lower'),
            'Bias_CI_Upper': metrics.get('bias_ci_upper'),
            'RMSE_CI_Lower': metrics.get('rmse_ci_lower'),
            'RMSE_CI_Upper': metrics.get('rmse_ci_upper'),
            'Spread_CI_Lower': metrics.get('spread_ci_lower'),
            'Spread_CI_Upper': metrics.get('spread_ci_upper'),
            'PCOR_CI_Lower': metrics.get('pcor_ci_lower'),
            'PCOR_CI_Upper': metrics.get('pcor_ci_upper'),
            'SSR_CI_Lower': metrics.get('ssr_ci_lower'),
            'SSR_CI_Upper': metrics.get('ssr_ci_upper')
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
