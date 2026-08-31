import os
import xarray as xr

class GPIReader:
    def __init__(self, data_dir, experiment, syear, eyear, gpi_vars=["GPI"], file_template=None):
        """
        Parameters:
        - data_dir: str, base directory where files are stored
        - experiment: str, name of the experiment (used in subfolder and filename)
        - syear, eyear: int, start and end year
        - gpi_vars: list of str, variable names to read (e.g., ["GPI", "GPI2010"])
        - file_template: optional override of filename template
        """
        self.data_dir = data_dir
        self.experiment = experiment
        self.syear = syear
        self.eyear = eyear
        self.gpi_vars = gpi_vars
        self.file_template = file_template  # allow optional override

        self.dataset = None        # Full xarray.Dataset
        self.gpi_data = {}         # Dict of {varname: DataArray}
        
    def normalize_longitude(self,da):
        if da.lon.max() > 180:
            da = da.assign_coords(lon=(((da.lon + 180) % 360) - 180))
            da = da.sortby('lon')
        return da

    def _get_file_path(self, year):
        """Constructs the file path for a given year based on experiment and filename convention."""
        if self.file_template:
            fname = self.file_template.format(exp=self.experiment, year=year)
        elif self.experiment == "ERA5":
            fname = f"{self.experiment}_monthly_{year}_h0_ens00.nc"
        else:
            fname = f"{self.experiment}_monthly_{year}_h0.nc"
        return os.path.join(self.data_dir, self.experiment, fname)

    def load_data(self):
        """Loads all requested GPI variants across years and stores them in self.gpi_data."""
        for varname in self.gpi_vars:
            datasets = []
            for year in range(self.syear, self.eyear + 1):
                fpath = self._get_file_path(year)
                if os.path.exists(fpath):
                    ds = xr.open_dataset(fpath)
                    ds = self.normalize_longitude(ds)
                    
                    if varname in ds:
                        datasets.append(ds[varname])
                    else:
                        print(f"⚠️ Variable '{varname}' not found in file: {fpath}")
                else:
                    print(f"⚠️ File not found: {fpath}")
            if datasets:
                self.gpi_data[varname] = xr.concat(datasets, dim="time")
            else:
                print(f"⚠️ No data loaded for variable: {varname}")
        
        # Store the dataset (first varname only)
        if self.gpi_data:
            first_var = self.gpi_vars[0]
            self.dataset = self.gpi_data[first_var].to_dataset(name=first_var)
        return self.gpi_data

    def compute_monthly_climatology(self, varname):
        """Compute 12-month climatology for a given variable name."""
        if varname not in self.gpi_data:
            raise ValueError(f"Variable '{varname}' not loaded.")
        return self.gpi_data[varname].groupby("time.month").mean("time")

    def get_data(self, varname):
        """Return full time series for a variable."""
        return self.gpi_data.get(varname, None)

    def get_dataset(self):
        """Return the base xarray dataset (uses first varname)."""
        return self.dataset

