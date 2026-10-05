"""
Latitude-band e3sm_diags lat_lon metrics (Test/Ref mean, STD, RMSE, PCC) for one
model climatology window, written as e3sm_diags-style ``{season}_metrics_table.csv``
files so ``SpinupMetricAnalyzer`` reads them exactly like the global tables.

Uses the xarray/xCDAT e3sm_diags (3.x) lat_lon driver pieces, so it runs under the
e3sm-unified environment (not the notebook kernel, which has no e3sm_diags):

    PY=/lcrc/soft/climate/e3sm-unified/e3smu_1_12_0/chrysalis/conda/envs/e3sm_unified_1.12.0_login/bin/python
    $PY -m util.regional_latlon_metrics \\
        --case 20240707.v3.LR.piControl.Rerun.GC2BC.chrysalis \\
        --climo-dir <MODEL_ROOT>/<case>/post/atm/180x360_aave/clim/50yr --period 0251-0300 \\
        --out-root <OUTPUT_ROOT>/data/model_error_compare/regional_metrics \\
        --catalog config/pltvar_region_sources.json

Each row reuses the default e3sm_diags lat_lon model_vs_obs parameters (obs file,
derived variable, pressure level, regrid tool/method, land/ocean mask). After the
driver's own subset + regrid to the lower resolution (``subset_and_align_datasets``),
cells outside the latitude band are set to NaN; the area-weighted metrics skip NaN.
Bands use grid-cell centre latitude and combine both hemispheres:

    global    no band mask (sanity check against the e3sm_diags tables)
    tropics   |lat| <  30
    midlat    30 <= |lat| < 60
    highlat   60 <= |lat| <= 90

Output: <out-root>/<band>/<case>/e3sm_diags/atm_monthly_180x360_aave/
        model_vs_obs_<period>/viewer/table-data/<season>_metrics_table.csv
"""

import argparse
import csv
import fnmatch
import glob
import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

BANDS = {
    "global":  None,
    "tropics": (0.0, 30.0),
    "midlat":  (30.0, 60.0),
    "highlat": (60.0, 90.0),
}
SEASONS = ["DJF", "MAM", "JJA", "SON"]
OBS_PATH = "/lcrc/soft/climate/e3sm_diags_data/obs_for_e3sm_diags/climatology"
TABLE_SUBDIR = "e3sm_diags/atm_monthly_180x360_aave/model_vs_obs_{period}/viewer/table-data"
COLS = ["Variables", "Unit", "Test_mean", "Ref._mean", "Mean_Bias",
        "Test_STD", "Ref._STD", "RMSE", "Correlation"]
MISSING = 999.999
# e3sm-unified environment with the xarray/xCDAT e3sm_diags (3.x)
E3SM_PYTHON = ("/lcrc/soft/climate/e3sm-unified/e3smu_1_12_0/chrysalis/conda/envs/"
               "e3sm_unified_1.12.0_login/bin/python")


def table_path(out_root, band, case, period, season):
    return os.path.join(out_root, band, case, TABLE_SUBDIR.format(period=period),
                        f"{season}_metrics_table.csv")


def _allowed_rows(catalog_path, exclude):
    """{(var_key, region, ref_name)} from the plot catalog, minus excluded variables."""
    with open(catalog_path) as f:
        cat = json.load(f)
    allowed = set()
    for var, val in cat.items():
        if any(fnmatch.fnmatch(var, pat) for pat in exclude):
            continue
        regmap = val["regions"] if isinstance(val, dict) and "regions" in val else val
        for region, srcs in regmap.items():
            for src in (srcs if isinstance(srcs, list) else [srcs]):
                allowed.add((var, region, src))
    return allowed


def _band(ds, var_key, band):
    """Copy of ds with var_key set to NaN outside the |lat| band (both hemispheres)."""
    import numpy as np
    import xcdat as xc

    if ds is None or BANDS[band] is None:
        return ds
    lo, hi = BANDS[band]
    alat = np.abs(xc.get_dim_coords(ds[var_key], axis="Y"))
    keep = (alat >= lo) & ((alat < hi) if hi < 90.0 else (alat <= hi))
    out = ds.copy()
    out[var_key] = ds[var_key].where(keep)
    return out


def _row(label, unit, var_key, ds_test, ds_ref, band):
    """One metrics-table row, computed as e3sm_diags does for the regridded fields."""
    from e3sm_diags.metrics.metrics import correlation, rmse, spatial_avg, std

    t = _band(ds_test, var_key, band)
    tm, ts = float(spatial_avg(t, var_key)), float(std(t, var_key))
    if ds_ref is None:
        return [label, unit, round(tm, 3), MISSING, MISSING, round(ts, 3), MISSING, MISSING, MISSING]
    r = _band(ds_ref, var_key, band)
    rm = float(spatial_avg(r, var_key))
    return [label, unit, round(tm, 3), round(rm, 3), round(tm - rm, 3), round(ts, 3),
            round(float(std(r, var_key)), 3),
            round(float(rmse(t, r, var_key)), 3),
            round(float(correlation(t, r, var_key)), 3)]


def _run_parameters(test_dir, case):
    """Default lat_lon model_vs_obs parameters, merged exactly as the e3sm_diags runner does."""
    from e3sm_diags.parameter.core_parameter import CoreParameter
    from e3sm_diags.run import runner

    param = CoreParameter()
    param.test_data_path = test_dir
    param.test_name = case
    param.reference_data_path = OBS_PATH
    param.results_dir = tempfile.mkdtemp(prefix="regional_latlon_")
    param.run_type = "model_vs_obs"
    param.multiprocessing = False
    runner.sets_to_run = ["lat_lon"]
    return runner.get_run_parameters([param])


def compute(case, climo_dir, period, out_root, catalog, exclude=(), bands=None,
            seasons=SEASONS, overwrite=False):
    from e3sm_diags.driver.lat_lon_driver import _get_ref_dataset
    from e3sm_diags.driver.utils.dataset_xr import Dataset
    from e3sm_diags.driver.utils.regrid import (
        get_z_axis, has_z_axis, regrid_z_axis_to_plevs, subset_and_align_datasets)

    bands = list(bands or BANDS)
    todo = [s for s in seasons
            if overwrite or not all(os.path.isfile(table_path(out_root, b, case, period, s)) for b in bands)]
    if not todo:
        print(f"[skip] {case} {period}: tables exist")
        return

    # e3sm_diags finds climo files by "<case>_<season>*" in test_data_path: expose
    # only this window's files through a symlink directory.
    # ncclimo names carry the season's own months, e.g. _JJA_000106_005008_climo.nc
    y0, y1 = (int(t) for t in period.split("-"))
    test_dir = tempfile.mkdtemp(prefix=f"climo_{period}_")
    for s in todo:
        pattern = os.path.join(climo_dir, f"{case}_{s}_{y0:04d}[0-9][0-9]_{y1:04d}[0-9][0-9]_climo.nc")
        found = sorted(glob.glob(pattern))
        if len(found) != 1:
            raise FileNotFoundError(f"expected one file for {pattern}, found {found}")
        os.symlink(found[0], os.path.join(test_dir, os.path.basename(found[0])))

    allowed = _allowed_rows(catalog, exclude)
    params = _run_parameters(test_dir, case)

    for season in todo:
        rows = {b: [] for b in bands}
        for p in params:
            p.seasons = [season]
            plevs = list(getattr(p, "plevs", []) or [])
            jobs = []
            for var_key in p.variables:
                keys = [f"{var_key}-{int(pl)}mb" for pl in plevs] if plevs else [var_key]
                wanted = [(k, rg) for k in keys for rg in p.regions if (k, rg, p.ref_name) in allowed]
                if wanted:
                    jobs.append((var_key, wanted))
            if not jobs:
                continue

            test_ds = Dataset(p, data_type="test")
            ref_ds = Dataset(p, data_type="ref")
            for var_key, wanted in jobs:
                p.var_id = var_key
                ds_test = test_ds.get_climo_dataset(var_key, season)
                ds_mask = test_ds._get_land_sea_mask(season)
                ds_ref = _get_ref_dataset(ref_ds, var_key, season)

                # {row key: (test, ref)} -- one entry per pressure level for 3-D variables
                if has_z_axis(ds_test[var_key]) and ds_ref is not None and has_z_axis(ds_ref[var_key]):
                    t_rg = regrid_z_axis_to_plevs(ds_test, var_key, plevs)
                    r_rg = regrid_z_axis_to_plevs(ds_ref, var_key, plevs)
                    z = get_z_axis(t_rg[var_key]).name
                    fields = {f"{var_key}-{int(pl)}mb": (t_rg.sel({z: pl}), r_rg.sel({z: pl})) for pl in plevs}
                else:
                    fields = {var_key: (ds_test, ds_ref)}

                for key, region in wanted:
                    if key not in fields:
                        continue
                    t, r = fields[key]
                    if r is None:   # model-only (no obs): test is not regridded
                        t_reg, r_reg = t, None
                    else:
                        _, t_reg, _, r_reg, _ = subset_and_align_datasets(p, t, r, ds_mask, var_key, region)
                    label = f"{key} {region} {p.ref_name}"
                    unit = t[var_key].attrs.get("units", "")
                    for b in bands:
                        rows[b].append(_row(label, unit, var_key, t_reg, r_reg, b))
                    print(f"[{case} {period} {season}] {label}", flush=True)

        for b in bands:
            # e3sm_diags keys table rows by label, so a later cfg entry with the same
            # label (e.g. TREFMNAV/TREFMXAV MERRA2) replaces the earlier one
            table = {row[0]: row for row in rows[b]}
            path = table_path(out_root, b, case, period, season)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "w") as f:
                w = csv.writer(f, delimiter=",", lineterminator="\n", quoting=csv.QUOTE_NONE)
                w.writerow(COLS)
                w.writerows(table.values())
        print(f"[done] {case} {period} {season}: {len(table)} rows", flush=True)


def run_windows(windows, out_root, catalog, exclude=(), bands=None, seasons=SEASONS,
                n_workers=8, python=E3SM_PYTHON, overwrite=False):
    """
    Compute the band tables for many (case, climo_dir, period) windows in parallel,
    one subprocess per window and season (the notebook kernel has no e3sm_diags).
    Windows/seasons whose tables already exist are skipped. Returns failed tasks.
    """
    bands = list(bands or BANDS)
    tasks = [(c, d, p, s) for c, d, p in windows for s in seasons
             if overwrite or not all(os.path.isfile(table_path(out_root, b, c, p, s)) for b in bands)]
    if not tasks:
        print(f"[regional] all {len(windows) * len(seasons)} window/season tables exist under {out_root}")
        return []

    log_dir = os.path.join(out_root, "logs")
    os.makedirs(log_dir, exist_ok=True)
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # workflow root (has util/)
    env = dict(os.environ, HDF5_USE_FILE_LOCKING="FALSE", OMP_NUM_THREADS="1")

    def _one(task):
        case, climo_dir, period, season = task
        cmd = [python, "-m", "util.regional_latlon_metrics", "--case", case,
               "--climo-dir", climo_dir, "--period", period, "--out-root", out_root,
               "--catalog", catalog, "--seasons", season, "--bands", *bands]
        if exclude:
            cmd += ["--exclude", *exclude]
        if overwrite:
            cmd += ["--overwrite"]
        log = os.path.join(log_dir, f"{case}_{period}_{season}.log")
        with open(log, "w") as f:
            rc = subprocess.run(cmd, cwd=root, env=env, stdout=f, stderr=subprocess.STDOUT).returncode
        print(f"[regional] {'ok  ' if rc == 0 else 'FAIL'} {case} {period} {season}"
              + ("" if rc == 0 else f"  (see {log})"), flush=True)
        return None if rc == 0 else task

    print(f"[regional] computing {len(tasks)} window/season tables with {n_workers} workers ...", flush=True)
    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        return [t for t in pool.map(_one, tasks) if t is not None]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", required=True)
    ap.add_argument("--climo-dir", required=True)
    ap.add_argument("--period", required=True, help="e.g. 0251-0300")
    ap.add_argument("--out-root", required=True)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--exclude", nargs="*", default=[])
    ap.add_argument("--bands", nargs="*", default=list(BANDS))
    ap.add_argument("--seasons", nargs="*", default=SEASONS)
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args(argv)
    compute(a.case, a.climo_dir, a.period, a.out_root, a.catalog, a.exclude,
            a.bands, a.seasons, a.overwrite)


if __name__ == "__main__":
    sys.exit(main())
