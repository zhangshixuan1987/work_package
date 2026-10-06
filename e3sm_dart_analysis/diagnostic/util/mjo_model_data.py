import os
import glob
import numpy as np
import xarray as xr


import pandas as pd
from typing import Dict, Tuple
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

    def parse_time_range(self, time):
        """Return an inclusive slice for a YYYYMM-YYYYMM period or input slice."""
        if isinstance(time, str) and '-' in time:
            start_date, end_date = time.split('-', 1)
            if len(start_date) != 6 or len(end_date) != 6:
                raise ValueError(f"Expected YYYYMM-YYYYMM, got: {time}")
            start = pd.Period(start_date, freq='M').start_time
            end = pd.Period(end_date, freq='M').end_time
            time_range = slice(start, end)
        elif isinstance(time, slice):
            if time.start is None or time.stop is None:
                raise ValueError("Both start and end must be specified in the slice.")
            start = pd.to_datetime(time.start)
            end = pd.to_datetime(time.stop)
            time_range = time
        else:
            raise TypeError("Only string ('YYYYMM-YYYYMM') or slice inputs are supported.")

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

            ds = xr.open_mfdataset(
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
