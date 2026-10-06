import os
import xarray as xr
import pandas as pd
import numpy as np

from xcdat.dataset import open_dataset
from xcdat.bounds import create_bounds
from xcdat.dataset import open_mfdataset

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
                    'period': '2010-2015'
                },
                'daily': {
                    'path': '/compyfs/zhan391/acme_init/Observations/GPCP/daily',
                    'template': 'PRECT.daily.%(year)s.nc',
                    'period': '2010-2015'
                }
            },
            'IMERG': {
                'daily': {
                    'path': '/compyfs/zhan391/acme_init/Observations/IMERG/daily',
                    'template': 'PRECT.daily.%(year)s.nc',
                    'period': '2007-2011'
                },
                'monthly': {
                    'path': '/compyfs/zhan391/acme_init/Observations/IMERG/monthly',
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
            'CPC_SOM': {
                'monthly': {
                    'path': '/compyfs/zhan391/acme_init/Observations/CPC_SOM/monthly',
                    'template': 'SOILWATER_10CM.monthly.%(year)s.nc',
                    'period': '2001-2020'
                }
            },
            'ESA_CCI': {
                'daily': {
                    'path': '/compyfs/zhan391/acme_init/Observations/ESA_CCI/daily',
                    'template': 'H2OSOI.monthly.%(year)s.nc',
                    'period': '2001-2020'
                },
                'monthly': {
                    'path': '/compyfs/zhan391/acme_init/Observations/ESA_CCI/monthly',
                    'template': 'H2OSOI.monthly.%(year)s.nc',
                    'period': '2001-2020'
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

    def parse_time_range(self,time):
        # Check if the input is a time range string like '201101-201202'
        if isinstance(time, str) and '-' in time:
            start_date, end_date = time.split('-')
            start = pd.to_datetime(start_date)
            end = pd.to_datetime(end_date) + pd.offsets.MonthEnd(0)

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

    def read(self, var, obsname, time, frequency=None, regnam='global', regional_mean=False):
        obs_group = self.get_group(obsname)
        if obs_group is None:
            raise ValueError(f"Observational dataset '{obsname}' not found.")
        if frequency not in obs_group:
            raise ValueError(f"Frequency '{frequency}' not found in '{obsname}'")

        group_info = obs_group[frequency]

        # identify time range and list of years
        time_range, years = self.parse_time_range(time)

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
            refds = open_mfdataset(
                paths,
                combine='nested',
                concat_dim='time',
                coords='minimal',
                compat='override',
                parallel=True
            )
            #refds = open_mfdataset(paths, decode_times=True, combine='by_coords', parallel=True)
        except Exception as e:
            raise RuntimeError(f"Failed to open datasets: {paths}") from e

        # Normalize longitude to [-180, 180]
        if refds.lon.min() >= 0:
            refds = refds.assign_coords(lon=((refds.lon + 180) % 360 - 180)).sortby("lon")

        if var not in refds:
            available = list(refds.data_vars)
            raise KeyError(f"Variable '{var}' not found. Available: {available}")

        data = self.convert_obs_units(obsname, var, refds[var])

        data = self._ensure_datetime64(data,frequency)

        # Apply region +  time selection
        (lat1, lat2), (lon1, lon2) = self.define_region(regnam)
        data = data.sel(lat=slice(lat1, lat2), lon=slice(lon1, lon2), time=time_range)

        # Compute regional mean if requested
        if regional_mean:
            weights = np.cos(np.deg2rad(data.lat))
            weights /= weights.mean()
            data = (data * weights).mean(dim=["lat", "lon"])

        return data.chunk({"time": -1}) if 'time' in data.dims else data

    def convert_obs_units(self, obsname, var, data_array):
        conversions = {
            'PRECT': lambda x: x * 86400.0 * 1000.0,  # mm/day to m/s
            'PSL': lambda x: x / 100.0,               # Pa to hPa
            'PS': lambda x: x / 100.0,                # Pa to hPa
        }

        # Ensure 'time' dimension is present
        if 'time' not in data_array.dims:
            raise ValueError(f"Expected 'time' dimension in data for '{var}'")

        # If variable is 'PRECT' and 'gpcp' is in the name, return as is
        if var == 'PRECT' and 'gpcp' in obsname.lower():
            return data_array

        # Perform conversion if a valid conversion exists for the variable
        if var in conversions:
            return conversions[var](data_array)

        # If no conversion exists for the variable, return data as is
        return data_array
