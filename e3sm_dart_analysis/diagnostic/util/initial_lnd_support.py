"""Support-file generation for initial-land diagnostics."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


def _sample_model_file(data_dict, sample_key: str, target_date: str, ensemble: str) -> Path:
    info = data_dict[sample_key]
    year = pd.to_datetime(target_date).year
    filename = info["template"] % {"year": year, "ensemble": ensemble}
    return Path(info["path"]) / filename


def _layer_thickness_from_midpoints(levgrnd: xr.DataArray) -> xr.DataArray:
    z = np.asarray(levgrnd.values, dtype=float)
    if z.ndim != 1 or z.size < 2:
        raise ValueError("levgrnd must be a 1-D coordinate with at least two levels.")

    interfaces = np.empty(z.size + 1, dtype=float)
    interfaces[0] = 0.0
    interfaces[1:-1] = 0.5 * (z[:-1] + z[1:])
    interfaces[-1] = z[-1] + 0.5 * (z[-1] - z[-2])
    dz = np.diff(interfaces)

    out = xr.DataArray(
        dz.astype("float32"),
        dims=(levgrnd.dims[0],),
        coords={levgrnd.dims[0]: levgrnd.values},
        name="DZSOI",
        attrs={
            "long_name": "soil layer thickness derived from levgrnd midpoints",
            "units": "m",
            "source": "Generated from model levgrnd coordinate for initial-land diagnostics.",
        },
    )
    return out


def ensure_initial_land_support_files(
    data_dict,
    *,
    landmask_file,
    soilayer_file,
    sample_key="CTRLEN10",
    target_date="2012-01-01",
    ensemble="EN01",
    overwrite=False,
):
    """Create landmask and DZSOI support files if they do not already exist."""

    landmask_path = Path(landmask_file)
    soilayer_path = Path(soilayer_file)
    landmask_path.parent.mkdir(parents=True, exist_ok=True)
    soilayer_path.parent.mkdir(parents=True, exist_ok=True)

    need_landmask = overwrite or not landmask_path.exists()
    need_soilayer = overwrite or not soilayer_path.exists()
    if not need_landmask and not need_soilayer:
        return str(landmask_path), str(soilayer_path)

    sample_path = _sample_model_file(data_dict, sample_key, target_date, ensemble)
    if not sample_path.exists():
        raise FileNotFoundError(
            f"Cannot generate initial-land support files because sample model file is missing: {sample_path}"
        )

    with xr.open_dataset(sample_path) as ds:
        if need_landmask:
            if "landfrac" not in ds:
                raise KeyError(f"Cannot generate {landmask_path}; 'landfrac' not found in {sample_path}")
            landmask = (ds["landfrac"] > 0.1).astype("float32").rename("landmask")
            landmask.attrs.update({
                "long_name": "land mask derived from landfrac > 0.1",
                "units": "1",
                "source": str(sample_path),
            })
            landmask.to_dataset().to_netcdf(landmask_path)

        if need_soilayer:
            if "DZSOI" in ds:
                dz = ds["DZSOI"]
                if "time" in dz.dims:
                    dz = dz.isel(time=0)
                dz = dz.rename("DZSOI")
            elif "levgrnd" in ds:
                dz = _layer_thickness_from_midpoints(ds["levgrnd"])
            else:
                raise KeyError(f"Cannot generate {soilayer_path}; neither 'DZSOI' nor 'levgrnd' found in {sample_path}")
            dz.to_dataset().to_netcdf(soilayer_path)

    return str(landmask_path), str(soilayer_path)
