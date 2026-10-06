import os
import xarray as xr
import pandas as pd
import numpy as np

class ObsDataReader:
    """
    Class to read and preprocess observational data with internal dataset registry.
    Supports regional mean extraction.
    """

    def __init__(self):
        self.obs_dict = self._define_obs_group()

    def _define_obs_group(self):
        return {
            'ERA5': {
                'monthly': {
                    'path': '/compyfs/zhan391/acme_init/Observations/ERA5/monthly',
                    'template': 'ERA5_analysis_monthly_%(year)s.nc',
                    'period': '2001-2018'
                },
                '6hourly': {
                    'path': '/compyfs/zhan391/acme_init/Observations/ERA5/6hourly',
                    'template': 'ERA5.6hourly.en00.%(var)s.%(year)s01-%(year)s12.nc',
                    'period': '1979-2018'
                }
            },
            'GPCP': {
                'monthly': {
                    'path': '/compyfs/zhan391/acme_init/Observations/GPCP/monthly',
                    'template': 'PRECT.monthly.%(year)s.nc',
                    'period': '1979-2017'
                },
                'daily': {
                    'path': '/compyfs/zhan391/acme_init/Observations/GPCP/daily',
                    'template': 'PRECT.daily.%(year)s.nc',
                    'period': '2007-2011'
                }
            },
            'IMERG': {
                'daily': {
                    'path': '/compyfs/zhan391/acme_init/Observations/IMERG/daily',
                    'template': 'PRECT.daily.%(year)s.nc',
                    'period': '2007-2011'
                }
            },
            'CERES-OAFlux': {
                'monthly': {
                    'path': '/compyfs/zhan391/acme_init/Observations/CERES-OAFlux/monthly',
                    'template': 'CERES-OAFlux_%(year)s.nc',
                    'period': '2001-2018',
                }
            },
            'MODIS_LST': {
                'monthly': {
                    'path': '/compyfs/zhan391/acme_init/Observations/MODIS_LST/monthly',
                    'template': 'MOD11C3.monthly.%(year)s.nc',
                    'period': '2001-2024'
                }
            }
        }

    @staticmethod
    def define_region(regnam='global'):
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

    def get_group(self, name=None):
        return self.obs_dict.get(name) if name else self.obs_dict['ERA5']

    def list_datasets(self):
        return list(self.obs_dict.keys())

    def _ensure_datetime64(self, da: xr.DataArray) -> xr.DataArray:
        """Convert time to datetime64 if cftime is used."""
        if not isinstance(da.indexes['time'], pd.DatetimeIndex):
            da['time'] = da.indexes['time'].to_datetimeindex(time_unit='ns')
        return da

    def read(self, var, obsname, time, frequency=None, regnam='global', regional_mean=False):
        obs_group = self.get_group(obsname)
        if obs_group is None:
            raise ValueError(f"Observational dataset '{obsname}' not found.")
        if frequency not in obs_group:
            raise ValueError(f"Frequency '{frequency}' not found in '{obsname}'")

        group_info = obs_group[frequency]

        # Convert time to list of years
        if isinstance(time, slice):
            time_range = pd.date_range(start=time.start, end=time.stop, freq="D")
        elif isinstance(time, (list, np.ndarray, pd.DatetimeIndex)):
            time_range = pd.to_datetime(time)
        else:
            time_range = pd.to_datetime([time])
        years = sorted(set(str(y) for y in time_range.year))

        # Generate file paths
        paths = []
        for y in years:
            try:
                fname = group_info["template"] % {'year': y, 'var': var}
            except KeyError:
                fname = group_info["template"] % {'year': y}
            paths.append(os.path.join(group_info["path"], fname))

        # Open and merge datasets
        try:
            refds = xr.open_mfdataset(paths, decode_times=True, combine='by_coords', parallel=True)
        except Exception as e:
            raise RuntimeError(f"Failed to open datasets: {paths}") from e

        # Normalize longitude to [-180, 180]
        if refds.lon.min() >= 0:
            refds = refds.assign_coords(lon=((refds.lon + 180) % 360 - 180)).sortby("lon")

        if var not in refds:
            available = list(refds.data_vars)
            raise KeyError(f"Variable '{var}' not found. Available: {available}")

        data = self.convert_obs_units(var, refds[var])
        data = self._ensure_datetime64(data)

        # Apply region
        (lat1, lat2), (lon1, lon2) = self.define_region(regnam)
        data = data.sel(lat=slice(lat1, lat2), lon=slice(lon1, lon2))

        # Apply time selection
        if isinstance(time, slice):
            data = data.sel(time=time)
        elif isinstance(time, (list, np.ndarray, pd.DatetimeIndex)):
            data = data.sel(time=time)
        else:
            data = data.sel(time=time, method='nearest', tolerance=np.timedelta64(1, 'D'))

        # Compute regional mean if requested
        if regional_mean:
            weights = np.cos(np.deg2rad(data.lat))
            weights /= weights.mean()
            data = (data * weights).mean(dim=["lat", "lon"])

        return data.chunk({"time": -1}) if 'time' in data.dims else data

    def convert_obs_units(self, var, data_array):
        conversions = {
            'PRECT': lambda x: x * 86400.0 * 1000.0,  # mm/day
            'TREFHT': lambda x: x - 273.15,           # °C
            'TS': lambda x: x - 273.15,               # °C
            'T850': lambda x: x - 273.15,             # °C
            'PSL': lambda x: x / 100.0,               # hPa
            'PS': lambda x: x / 100.0,                # hPa
        }

        if 'time' not in data_array.dims:
            raise ValueError(f"Expected 'time' dimension in data for '{var}'")

        return conversions[var](data_array) if var in conversions else data_array
