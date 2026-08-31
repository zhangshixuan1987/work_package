from typing import Optional
import os
import numpy as np
import xarray as xr


class AROccurrenceFrequencyCalculator:
    def __init__(self, data_dir, file_template="CLIM_ARS_{}.nc",
                 varname="AR_binary_tag", syear=2008, eyear=2017):
        self.data_dir = data_dir
        self.file_template = file_template
        self.varname = varname
        self.syear = syear
        self.eyear = eyear
        self.monthly_counts = None
        self.annual_days = None
        self._lat = None
        self._lon = None
        self.composited_data = {}
        
    def _normalize_coords(self, da):
        rename_dict = {}
        if "latitude" in da.coords:
            rename_dict["latitude"] = "lat"
        if "longitude" in da.coords:
            rename_dict["longitude"] = "lon"
        da = da.rename(rename_dict)

        if "lon" in da.coords:
            da = da.assign_coords(lon=(((da.lon + 180) % 360) - 180))
            da = da.sortby("lon")

        return da

    def process_single_year(self, year, output_path):
        fname = os.path.join(self.data_dir, self.file_template.format(year))
        if not os.path.exists(fname):
            print(f"[WARN] File not found: {fname}")
            return

        ds = xr.open_dataset(fname)
        if self.varname not in ds:
            print(f"[WARN] Variable '{self.varname}' not found in {fname}")
            return

        ds = self._normalize_coords(ds)
        mask = ds[self.varname]  # [time, lat, lon]
        if self._lat is None or self._lon is None:
            self._lat = ds["lat"]
            self._lon = ds["lon"]

        mask["time"] = xr.decode_cf(ds).time
        mask["month"] = mask["time"].dt.month

        monthly_sum = np.zeros((12,) + mask.shape[1:], dtype=np.float32)
        monthly_npts = np.zeros(12, dtype=int)
        monthly_days = np.zeros(12, dtype=int)

        grouped = mask.groupby("month")
        for m, da_m in grouped:
            monthly_sum[m - 1] += da_m.sum(dim="time").values
            monthly_npts[m - 1] += da_m.sizes["time"]
            monthly_days[m - 1] = len(np.unique(da_m["time"].dt.day.values))

        # === Monthly frequency (fraction of 6-hourly timesteps with AR)
        monthly_mean = monthly_sum / monthly_npts[:, np.newaxis, np.newaxis]
        monthly_counts = xr.DataArray(
            monthly_mean,
            dims=("month", "lat", "lon"),
            coords={"month": np.arange(1, 13), "lat": self._lat, "lon": self._lon},
            name="monthly_ar_frequency"
        )
        
        monthly_counts.attrs.update({
            "description": "Monthly mean AR occurrence frequency (fraction of 6-hourly timesteps)",
            "units": "fraction (0–1)",
            "source": "6-hourly AR_binary_tag mask"
        })

        # === Monthly AR days
        ar_days_per_month = (monthly_mean * monthly_days[:, np.newaxis, np.newaxis])
        monthly_days_da = xr.DataArray(
            ar_days_per_month,
            dims=("month", "lat", "lon"),
            coords={"month": np.arange(1, 13), "lat": self._lat, "lon": self._lon},
            name="monthly_ar_days"
        )
        monthly_days_da.attrs.update({
            "description": "Monthly mean number of days with AR occurrence",
            "units": "days/month",
            "source": "Derived from 6-hourly AR_binary_tag mask"
        })

        # === Annual AR days
        annual_days = (mask.sum(dim="time") / mask.sizes["time"]) * 365 * 4
        annual_days = annual_days.assign_coords(lat=self._lat, lon=self._lon)
        annual_days.name = "annual_ar_days"
        annual_days.attrs.update({
            "description": "Annual mean number of days with AR occurrence",
            "units": "days/year",
            "source": "Derived from 6-hourly AR_binary_tag mask"
        })

        # === Save to NetCDF
        ds_out = xr.Dataset({
            "monthly_ar_frequency": monthly_counts,
            "monthly_ar_days": monthly_days_da,
            "annual_ar_days": annual_days
        })
        ds_out.to_netcdf(output_path)
        print(f"[SAVED] {output_path}")


    def compute_composited_means(self, year, varname, varsub=None, output_path=None):
        if "ERA5" in self.data_dir:
            var_alias = {"TUQ": "VIWVE", "TVQ": "VIWVN"}
        else:
            var_alias = {}

        # Determine variable names
        input_varname = var_alias.get(varname, varname)
        actual_varname = f"{input_varname}_{varsub}" if varsub else input_varname

        fname = os.path.join(self.data_dir, self.file_template.format(year).replace("nff", input_varname))

        if not os.path.exists(fname):
            print(f"[WARN] File not found: {fname}")
            return

        ds = xr.open_dataset(fname)
        ds = self._normalize_coords(ds)

        # === Handle IVT derivation if needed ===
        if varname == "IVT" and actual_varname not in ds:
            # Derive IVT from TUQ and TVQ (or their aliases)
            tuq_input = var_alias.get("TUQ", "TUQ") 
            tvq_input = var_alias.get("TVQ", "TVQ") 
            tuq_name = var_alias.get("TUQ", "TUQ") + (f"_{varsub}" if varsub else "")
            tvq_name = var_alias.get("TVQ", "TVQ") + (f"_{varsub}" if varsub else "")
            f1 = os.path.join(self.data_dir, self.file_template.format(year).replace("nff", tuq_input))
            ds1 = xr.open_dataset(f1)
            ds1 = self._normalize_coords(ds1)
            f2 = os.path.join(self.data_dir, self.file_template.format(year).replace("nff", tvq_input))
            ds2 = xr.open_dataset(f2)
            ds2 = self._normalize_coords(ds2)
            if tuq_name not in ds1 or tvq_name not in ds2:
                print(f"[WARN] TUQ ({tuq_name}) or TVQ ({tvq_name}) not found in {fname}")
                return
            da_tuq = ds1[tuq_name]
            da_tvq = ds2[tvq_name]
            da = np.sqrt(da_tuq**2 + da_tvq**2)
            da.name = "IVT"
        else:
            if actual_varname not in ds:
                print(f"[WARN] Variable '{actual_varname}' not found in {fname}")
                return
            da = ds[actual_varname]

        # Time decoding and coordinate check
        da["time"] = xr.decode_cf(ds).time
        da["month"] = da["time"].dt.month
        if self._lat is None or self._lon is None:
            self._lat = da["lat"]
            self._lon = da["lon"]

        # Compute means
        monthly_mean = da.groupby("month").mean("time")
        annual_mean = da.mean("time")

        ds_out = xr.Dataset({
            f"ar_composite_monthly_{varname}": monthly_mean,
            f"ar_composite_annual_{varname}": annual_mean
        })
        ds_out.to_netcdf(output_path)
        print(f"[SAVED] {output_path}")

