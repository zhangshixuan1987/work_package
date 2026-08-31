import os
import numpy as np
import pandas as pd
import xarray as xr
from tc_basins import BASINS
from collections import defaultdict

def get_tc_basin(): 
    BASINS = {
        "North Atlantic":          {"lat1": 0,   "lat2": 45,  "lon1": -100, "lon2": -10},
        "Eastern Pacific":         {"lat1": 0,   "lat2": 35,  "lon1": -140, "lon2": -100},
        "Western Pacific":         {"lat1": 0,   "lat2": 45,  "lon1": 100,  "lon2": 180},
        "North Indian Ocean":      {"lat1": 0,   "lat2": 30,  "lon1": 40,   "lon2": 100},
        "Southwest Indian Ocean":  {"lat1": -40, "lat2": 0,   "lon1": 30,   "lon2": 100},
        "Southeast Indian Ocean":  {"lat1": -40, "lat2": 0,   "lon1": 100,  "lon2": 135},
        "South Pacific":           {"lat1": -40, "lat2": 0,   "lon1": 135,  "lon2": -100},
    
        # Added hemispheric regions
        "Northern Hemisphere":  {"lat1": 0,   "lat2": 60,  "lon1": -180, "lon2": 180},
        "Southern Hemisphere":  {"lat1": -60, "lat2": 0,   "lon1": -180, "lon2": 180}
    }
    
    return BASINS
    
def load_ibtracs_genesis_density(ibtracs_file, syear, eyear, res=1.0):
    """
    Loads TC genesis density from IBTrACS and returns a monthly gridded count: [time, lat, lon].

    Parameters:
    - ibtracs_file: str, path to IBTrACS NetCDF
    - syear, eyear: int, time window
    - res: float, lat-lon grid resolution (default 1.0)

    Returns:
    - genesis_da: xarray.DataArray with dims [time, lat, lon], monthly genesis counts
    """
    ds = xr.open_dataset(ibtracs_file)

    # Parse time and genesis info
    times = pd.to_datetime(ds["time"].values)
    lats = ds["lat"].values
    lons = ds["lon"].values
    storm_ids = ds["sid"].values

    # Convert lon to [-180, 180]
    lons = ((lons + 180) % 360) - 180

    # Keep only genesis points (first time per storm)
    df = pd.DataFrame({"time": times, "lat": lats, "lon": lons, "sid": storm_ids})
    df = df.dropna()
    df = df.sort_values(["sid", "time"])
    genesis_df = df.groupby("sid").first().reset_index()

    # Filter by year
    genesis_df = genesis_df[
        (genesis_df["time"].dt.year >= syear) & (genesis_df["time"].dt.year <= eyear)
    ]

    # Define lat/lon bins
    lat_bins = np.arange(-90, 90 + res, res)
    lon_bins = np.arange(-180, 180 + res, res)

    # Bin monthly counts
    time_range = pd.date_range(f"{syear}-01", f"{eyear}-12", freq="MS")
    data = np.zeros((len(time_range), len(lat_bins)-1, len(lon_bins)-1), dtype=int)

    for i, t in enumerate(time_range):
        month_df = genesis_df[
            (genesis_df["time"].dt.year == t.year) & (genesis_df["time"].dt.month == t.month)
        ]
        H, _, _ = np.histogram2d(
            month_df["lat"], month_df["lon"], bins=[lat_bins, lon_bins]
        )
        data[i] = H

    # Create coordinates
    lat_centers = 0.5 * (lat_bins[:-1] + lat_bins[1:])
    lon_centers = 0.5 * (lon_bins[:-1] + lon_bins[1:])
    genesis_da = xr.DataArray(
        data,
        coords={"time": time_range, "lat": lat_centers, "lon": lon_centers},
        dims=["time", "lat", "lon"],
        name="genesis_density"
    )

    return genesis_da
    
class IBTrACSAnalyzer:
    def __init__(self, ibtracs_path, syear, eyear):
        self.ibtracs_path = ibtracs_path
        self.syear = syear
        self.eyear = eyear
        self.ibtracs_df = None

    def assign_basin(self, lat, lon):
        for basin, b in BASINS.items():
            lon1, lon2 = b["lon1"], b["lon2"]
            if lon1 > lon2:
                in_lon = (lon >= lon1 or lon <= lon2)  # crosses dateline
            else:
                in_lon = (lon >= lon1 and lon <= lon2)
            in_lat = (lat >= b["lat1"] and lat <= b["lat2"])
            if in_lat and in_lon:
                return basin
        return None
        
    def __repr__(self):
        return f"IBTrACSAnalyzer(years={self.syear}-{self.eyear})"

    def load_and_filter_ibtracs(self):
        ds = xr.open_dataset(self.ibtracs_path)
        df = ds.to_dataframe().reset_index()

        # Basic QC and decoding
        df = df.dropna(subset=["lat", "lon", "iso_time", "wmo_wind", "wmo_pres"])
        df["datetime"] = pd.to_datetime(df["iso_time"].apply(lambda x: x.decode("utf-8") if isinstance(x, bytes) else x))
        df["storm_id"] = df["sid"]
        df["storm_name"] = df["name"].apply(lambda x: x.decode("utf-8") if isinstance(x, bytes) else x).str.strip()
        df["lon"] = df["lon"].apply(lambda x: x - 360 if x > 180 else x)
        df["wmo_wind"] = pd.to_numeric(df["wmo_wind"], errors="coerce")
        df["wmo_pres"] = pd.to_numeric(df["wmo_pres"], errors="coerce")

        # Temporal and quality filtering
        df = df[(df["datetime"].dt.year >= self.syear) & (df["datetime"].dt.year <= self.eyear)]
        df = df[df["datetime"].dt.hour.isin([0, 6, 12, 18])]
        df = df[(df["wmo_wind"] >= 34) & (df["lat"].between(-50.0, 50.0))]

        # Apply duration + time gap filter
        valid_ids = []
        for sid, sdf in df.groupby("storm_id"):
            sdf = sdf.sort_values("datetime")
            duration = (sdf["datetime"].max() - sdf["datetime"].min()).total_seconds() / 3600.0
            time_diffs = sdf["datetime"].diff().dt.total_seconds().dropna() / 3600.0
            if duration >= 24 and (time_diffs <= 6).all():
                valid_ids.append(sid)

        self.ibtracs_df = df[df["storm_id"].isin(valid_ids)]
        return self.ibtracs_df

    def save_yearly_density(self, out_dir="./", bins=(4.0, 4.0)):
        os.makedirs(out_dir, exist_ok=True)
        lat_bins = np.arange(-90,  91, bins[0])
        lon_bins = np.arange(-180, 181, bins[1])
        lat_centers = 0.5 * (lat_bins[:-1] + lat_bins[1:])
        lon_centers = 0.5 * (lon_bins[:-1] + lon_bins[1:])

        for year in range(self.syear, self.eyear + 1):
            df_year = self.ibtracs_df[self.ibtracs_df["datetime"].dt.year == year]
            if df_year.empty:
                print(f"[SKIP] No data for {year}")
                continue

            out_file = os.path.join(out_dir, f"tc_metrics_IBTrACS_{year}_bin{int(bins[0])}x{int(bins[1])}.nc")
            if os.path.exists(out_file):
                print(f"[SKIP] Already exists: {out_file}")
                continue

            track_H, _, _ = np.histogram2d(df_year["lat"], df_year["lon"], bins=[lat_bins, lon_bins])
            genesis_df = df_year.sort_values("datetime").groupby("storm_id").first().reset_index()
            genesis_H, _, _ = np.histogram2d(genesis_df["lat"], genesis_df["lon"], bins=[lat_bins, lon_bins])
            lysis_df = df_year.sort_values("datetime").groupby("storm_id").last().reset_index()
            lysis_H, _, _ = np.histogram2d(lysis_df["lat"], lysis_df["lon"], bins=[lat_bins, lon_bins])

            ds = xr.Dataset(
                data_vars={
                    "track_density": (("lat", "lon"), track_H),
                    "genesis_density": (("lat", "lon"), genesis_H),
                    "lysis_density": (("lat", "lon"), lysis_H),
                },
                coords={"lat": lat_centers, "lon": lon_centers},
                attrs={
                    "title": f"TC Metrics (IBTrACS) {year}",
                    "description": "TC track, genesis, and lysis density from filtered IBTrACS",
                    "year": year,
                    "bin_size_degrees": f"{bins[0]}x{bins[1]}"
                }
            )
            ds.to_netcdf(out_file)
            print(f"[SAVED] {out_file}")

    def save_multiyear_climatology(self, metric_dir, out_file, bins=(5.0, 5.0)):
        files = [
            os.path.join(metric_dir, f"tc_metrics_IBTrACS_{y}_bin{int(bins[0])}x{int(bins[1])}.nc")
            for y in range(self.syear, self.eyear + 1)
        ]
        valid_files = [f for f in files if os.path.exists(f)]
        if not valid_files:
            print("[ERROR] No valid yearly NetCDF files found.")
            return

        ds_all = xr.concat([xr.open_dataset(f) for f in valid_files], dim="year")
        ds_clim = ds_all.mean(dim="year")
        ds_clim.attrs.update({
            "title": f"IBTrACS TC Metrics Climatology ({self.syear}-{self.eyear})",
            "description": "Mean of yearly track, genesis, and lysis density",
        })
        ds_clim.to_netcdf(out_file)
        print(f"[SAVED] {out_file}")
        
    def get_genesis_climatology_da(self, metric_dir, bins=(1.0, 1.0), save_path=None):
        """
        Computes and returns multi-year mean gridded genesis density from filtered IBTrACS data.

        Parameters:
        - metric_dir: directory to store/load yearly density files
        - bins: tuple (lat_bin_deg, lon_bin_deg), default 1x1
        - save_path: optional NetCDF path to save the climatology

        Returns:
        - genesis_density: xarray.DataArray [lat, lon]
        """
        if self.ibtracs_df is None:
            self.load_and_filter_ibtracs()

        self.save_yearly_density(out_dir=metric_dir, bins=bins)

        clim_file = os.path.join(
            metric_dir,
            f"tc_metrics_IBTrACS_climatology_{self.syear}-{self.eyear}_bin{int(bins[0])}x{int(bins[1])}.nc"
        )
        self.save_multiyear_climatology(metric_dir, clim_file, bins=bins)

        ds = xr.open_dataset(clim_file)
        genesis_da = ds["genesis_density"]

        if save_path:
            genesis_da.to_netcdf(save_path)
            print(f"[SAVED] Genesis density climatology: {save_path}")

        return genesis_da
    
    def save_monthly_genesis_climatology(self, out_path, bins=(1.0, 1.0)):
        lat_bins = np.arange(-90, 91, bins[0])
        lon_bins = np.arange(-180, 181, bins[1])
        lat_centers = 0.5 * (lat_bins[:-1] + lat_bins[1:])
        lon_centers = 0.5 * (lon_bins[:-1] + lon_bins[1:])

        monthly_genesis = []
        times = []

        for month in range(1, 13):
            df_month = self.ibtracs_df[self.ibtracs_df["datetime"].dt.month == month]
            genesis_df = df_month.sort_values("datetime").groupby("storm_id").first().reset_index()
            H, _, _ = np.histogram2d(genesis_df["lat"], genesis_df["lon"], bins=[lat_bins, lon_bins])
            monthly_genesis.append(H[np.newaxis, ...])
            times.append(pd.Timestamp(f"2000-{month:02d}-01"))

        da = xr.DataArray(
            np.concatenate(monthly_genesis, axis=0),
            coords={"time": times, "lat": lat_centers, "lon": lon_centers},
            dims=["time", "lat", "lon"],
            name="genesis_density"
        )
        da.attrs["description"] = "Monthly TC genesis density from IBTrACS"
        da.to_netcdf(out_path)
        print(f"[SAVED] Monthly genesis climatology → {out_path}")

    def get_basin_monthly_genesis_and_density_timeseries(self, bins=(1.0, 1.0), normalize=False):
        """
        Compute monthly TC genesis and track density counts per basin using histogram + basin mask.

        Parameters:
            bins (tuple): lat/lon bin size in degrees
            normalize (bool): if True, normalize counts by number of grid cells in basin

        Returns:
            genesis_df: pd.DataFrame [time, basin]
            density_df: pd.DataFrame [time, basin]
        """
        if self.ibtracs_df is None:
            self.load_and_filter_ibtracs()

        # Define lat/lon bins and centers
        lat_bins = np.arange(-90, 91, bins[0])
        lon_bins = np.arange(-180, 181, bins[1])
        lat_centers = 0.5 * (lat_bins[:-1] + lat_bins[1:])
        lon_centers = 0.5 * (lon_bins[:-1] + lon_bins[1:])
        lon2d, lat2d = np.meshgrid(lon_centers, lat_centers)

        # Precompute basin masks and their normalization factors
        basin_masks = {}
        basin_norm = {}  # number of grid cells per basin
        for basin, b in BASINS.items():
            lon1, lon2 = b["lon1"], b["lon2"]
            if lon1 > 180 or lon2 > 180:
                lon1 = ((lon1 + 180) % 360) - 180
                lon2 = ((lon2 + 180) % 360) - 180

            if lon1 > lon2:
                lon_mask = (lon2d >= lon1) | (lon2d <= lon2)
            else:
                lon_mask = (lon2d >= lon1) & (lon2d <= lon2)

            lat_mask = (lat2d >= b["lat1"]) & (lat2d <= b["lat2"])
            mask = lon_mask & lat_mask
            basin_masks[basin] = mask
            basin_norm[basin] = np.sum(mask)

        # Time series containers
        all_times = []
        genesis_series = {basin: [] for basin in BASINS}
        track_series = {basin: [] for basin in BASINS}

        for year in range(self.syear, self.eyear + 1):
            for month in range(1, 13):
                time = pd.Timestamp(f"{year}-{month:02d}-01")
                df_month = self.ibtracs_df[self.ibtracs_df["datetime"].dt.to_period("M") == time.to_period("M")]

                if df_month.empty:
                    for basin in BASINS:
                        genesis_series[basin].append(0)
                        track_series[basin].append(0)
                    all_times.append(time)
                    continue

                # Histograms
                track_H, _, _ = np.histogram2d(df_month["lat"], df_month["lon"], bins=[lat_bins, lon_bins])
                genesis_df = df_month.sort_values("datetime").groupby("storm_id").first().reset_index()
                genesis_H, _, _ = np.histogram2d(genesis_df["lat"], genesis_df["lon"], bins=[lat_bins, lon_bins])

                for basin, mask in basin_masks.items():
                    gval = np.sum(genesis_H[mask])
                    tval = np.sum(track_H[mask])
                    if normalize and basin_norm[basin] > 0:
                        gval /= basin_norm[basin]
                        tval /= basin_norm[basin]
                    genesis_series[basin].append(gval)
                    track_series[basin].append(tval)

                all_times.append(time)

        # Convert to DataFrames
        genesis_df = pd.DataFrame(genesis_series, index=all_times)
        track_df = pd.DataFrame(track_series, index=all_times)

        return genesis_df, track_df


    def save_basin_timeseries_to_netcdf(self, genesis_df, density_df, out_path):
        """
        Save monthly TC genesis and track density counts by basin to NetCDF.

        Parameters:
            genesis_df (pd.DataFrame): Monthly genesis counts by basin
            density_df (pd.DataFrame): Monthly track density counts by basin
            out_path (str): Output NetCDF path
        """
        # Convert index to datetime
        time = pd.to_datetime(genesis_df.index)

        # Create xarray Dataset
        ds = xr.Dataset(
            {
                "genesis_count": (["time", "basin"], genesis_df.values),
                "track_density_count": (["time", "basin"], density_df.values),
            },
            coords={
                "time": time,
                "basin": genesis_df.columns.tolist(),
            },
            attrs={
                "description": "Monthly TC genesis and track density counts by basin",
                "source": "Filtered IBTrACS",
            }
        )
        ds.to_netcdf(out_path)
        print(f"[SAVED] Basin TC monthly timeseries → {out_path}")
