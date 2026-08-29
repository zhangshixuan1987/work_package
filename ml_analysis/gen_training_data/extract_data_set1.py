import os
import sys
import numpy as np
import xarray as xr
from dask.distributed import Client

# OPTIONAL: Start Dask distributed client for parallelism
client = Client()
print("Dask client started:", client)

# --- Input Arguments ---
if len(sys.argv) > 2:
    syear = int(sys.argv[1])
    eyear = int(sys.argv[2])
else:
    print("Must specify year range!")
    sys.exit(1)

topdir = "/global/cfs/cdirs/e3sm/www/zhan391/darpa_temporary_data_share"
outdir = "/global/cfs/cdirs/e3sm/www/zhan391/SEA_CROGS"
exp = "E3SMv2_NDGUV_tau6"

groups = ["before_nudging", "nudging_tendency"]

for grp in groups:
    if grp == "before_nudging":
        smdh = "01010000"
        emdh = "12312100"
        var_list = [
            "WS_bf_ndg", "WD_bf_ndg", "PRESSURE_bf_ndg",
            "U_bf_ndg", "V_bf_ndg", "T_bf_ndg",
            "Q_bf_ndg", "PS_bf_ndg"
        ]
        var_outs = ["WS", "WD", "PRESSURE", "U", "V", "T", "Q", "PS"]
    else:
        smdh = "01010300"
        emdh = "01010000"
        var_list = ["Nudge_U", "Nudge_V", "Nudge_T", "Nudge_Q"]
        var_outs = ["UTEND", "VTEND", "TEND", "QTEND"]

    datadir = os.path.join(topdir, "SE_PG2", grp)
    outpath = os.path.join(outdir, "New_Training", exp)
    os.makedirs(outpath, exist_ok=True)

    for year in range(syear, eyear + 1):
        for var, vou in zip(var_list, var_outs):
            print(f"\nProcessing variable: {var}")

            if "Nudge" in var:
                fname = f"{exp}_3hourly_{year:04d}{smdh}-{year + 1:04d}{emdh}.nc"
            else:
                fname = f"{exp}_3hourly_{year:04d}{smdh}-{year:04d}{emdh}.nc"

            fpath = os.path.join(datadir, fname)
            if not os.path.exists(fpath):
                print(f"File not found: {fpath}")
                continue

            # Open with Dask chunking for memory efficiency
            try:
                ds = xr.open_dataset(fpath, chunks={"time": 100})
            except Exception as e:
                print(f"Failed to open {fname}: {e}")
                continue

            lat = ds['lat'].values
            lon = ds['lon'].values
            lev = ds['lev'].values
            time = ds['time'].values

            print(f"latitude range: {lat.min()} to {lat.max()}")
            print(f"longitude range: {lon.min()} to {lon.max()}")
            print(f"level range: {lev.min()} to {lev.max()}")
            print(f"time range: {time.min()} to {time.max()}")

            try:
                # --- Derived or direct variable extraction ---
                if var == "PRESSURE_bf_ndg":
                    hyam = ds['hyam'].values
                    hybm = ds['hybm'].values
                    P0 = ds['P0'].values
                    PS = ds['PS_bf_ndg']
                    pressure = []
                    for k in range(len(lev)):
                        pressure.append(hyam[k] * P0 + hybm[k] * PS)
                    data = xr.concat(pressure, dim="lev")
                    data_np = data.transpose("time", "lev", "lat").compute().astype(np.float32)

                elif var == "PS_bf_ndg":
                    data_np = ds[var].compute().astype(np.float32)

                elif var in ["WS_bf_ndg", "WD_bf_ndg"]:
                    U = ds["U_bf_ndg"]
                    V = ds["V_bf_ndg"]
                    WS = (U**2 + V**2) ** 0.5
                    WD = np.arctan2(V, U)
                    data = WS if var == "WS_bf_ndg" else WD
                    data_np = data.compute().astype(np.float32)

                else:
                    data_np = ds[var].compute().astype(np.float32)

            except Exception as e:
                print(f"Error reading variable {var}: {e}")
                ds.close()
                continue

            # --- Save output ---
            outname = f"{vou}_3hourly_{year:04d}.npy"
            fout = os.path.join(outpath, outname)

            if os.path.exists(fout):
                os.remove(fout)

            np.save(fout, data_np)
            print(f"Saved: {fout}")

            ds.close()
            del data_np
    del smdh, emdh, var_list, var_outs

