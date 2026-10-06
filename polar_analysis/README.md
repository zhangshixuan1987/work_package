# E3SM v3 polar analysis

Data-processing and plotting workflows for E3SM v3 polar diagnostics:
regional time series, mean-state biases, surface-temperature budgets, and Arctic
amplification (AA). The package serves two projects:

| Project | Output tree (current machine) | Origin |
| --- | --- | --- |
| v3 large-ensemble paper | `/lcrc/group/e3sm/public_html/diagnostic_output/ac.szhang/v3LE_paper` | `v3_le_paper/1_*`–`5_*`, `script/`, `data/` |
| v3 polar regional analysis | `/lcrc/group/e3sm/public_html/diagnostic_output/ac.szhang/v3_polar_analysis` | `v3_polar_analysis/<region>/`, `data/` |

Both trees use `data/` (processed inputs), `diag_data/` (diagnostic products),
`figures/`, and `analysis_v0/` (the original folders, archived unchanged).
Nothing in the package hard-codes these locations: every notebook takes its
directories as parameters, and every shell script takes `DATA_DIR`.

## Layout

```text
polar_analysis/
├── config/regions.json  per-region settings (mask, reference, colour scales, ...)
├── jupyter/
│   ├── process/         process_*.ipynb: read inputs, write data/ or diag_data/
│   └── plot/            plot_*.ipynb: read processed data, write figures/
├── scripts/             shared modules (exp_info, regions, data_utils)
└── shell_script/        shell/Slurm preprocessing (climatologies, regrid maps, MPAS)
```

## Running a notebook

1. Use a kernel with xarray, xcdat, cartopy, global_land_mask, and xskillscore
   (on Chrysalis: the `zppy-pcmdi-diags` conda environment).
2. Run the **setup cell** (identical in every notebook). It finds the package
   from the working directory and puts `scripts/` on `sys.path`.
3. Edit the **parameters cell** (tagged `parameters`) for your machine and
   experiment: input directories (`CLIMO_DIR`, `REGMN_DIR`, `INPUT_DIR`, ...),
   outputs (`DIAG_DIR`, `FIG_DIR`), and `REGION` where relevant. Its defaults are
   this machine's paths.
4. Run the remaining cells.

To run a notebook non-interactively with other values, use papermill, which
injects the overrides after the parameters cell:

```bash
papermill jupyter/plot/plot_bias_zonal.ipynb out.ipynb \
  -p REGION Greenland \
  -p CLIMO_DIR /lcrc/group/e3sm/public_html/diagnostic_output/ac.szhang/v3_polar_analysis/data/climo \
  -p DIAG_DIR  /lcrc/group/e3sm/public_html/diagnostic_output/ac.szhang/v3_polar_analysis/diag_data \
  -p FIG_DIR   /lcrc/group/e3sm/public_html/diagnostic_output/ac.szhang/v3_polar_analysis/figures/mean_bias
```

## Regions

The climatology-based notebooks (`plot_bias_*`, `plot_budget_ts_{sea,zonal}`,
`plot_budget_regmean`) run for any region in `config/regions.json` by setting
`REGION`. Region bounds come from `scripts/exp_info.py` (`REGION_CATALOG`); the
JSON adds the per-region choices that used to live in separate copies:

| Region | Surface mask | Reference | Other |
| --- | --- | --- | --- |
| Arctic, Greenland | land | MODIS | — |
| Antarctic | land | MODIS | drops HadISST; ±10 colour scales; budget band −90…−65° |
| Atlantic | ocean | HadISST | TS reference HadISST |

`scripts/regions.py` provides `region_settings(REGION)` and
`region_catalogs(REGION)` (run and variable catalogs with the region's
exclusions and overrides). They reproduce the per-region `exp_info.py` copies
from `v3_polar_analysis` exactly. Add a region by adding its bounds to
`REGION_CATALOG` and an entry to `config/regions.json`.

## Notebooks

Names follow `<stage>_<diagnostic>_<detail>_<variant>.ipynb`.

### Data processing (`jupyter/process/`)

| Step | Notebooks | Purpose |
| --- | --- | --- |
| 1 | `process_3d_{e3sm,era5,hadisst,noaa20c,noaav5,picontrol}_data`, `process_3d_picontrol_trend`, `process_seaice_ts` | Regrid/subset 3D fields, the piControl trend, sea-ice time series → `data/global/`. |
| 2 | `process_ts_{e3sm,era5,noaa20c}` | Regional-mean time series → `data/regmn/`. |
| 3 | `process_3d_e3sm_driftcorr`, `process_ts_e3sm_driftcorr` | piControl drift correction → `data/*_driftcorr/`. |
| — | `process_ts_{atm,sea_ice}_imsk`, `process_ts_nmsk` | Regional anomaly statistics (ice mask / no mask) → `diag_data/`. |
| — | `process_aa_{ts,ac,2d}_{65n,65n_driftcorr,50n}` | Arctic amplification time series, annual cycle, maps → `diag_data/`. |

Climatologies (`data/climo/`) and regrid maps (`data/map/`) are built by
`shell_script/run_process_*.bash` and `run_gen_map.sh`.

### Plots (`jupyter/plot/`)

| Notebooks | Figures |
| --- | --- |
| `plot_ts_sea_ice_imsk` | Seasonal sea-ice time-series panels. |
| `plot_bias_map`, `plot_bias_zonal` | Regional bias maps and zonal-mean biases (any `REGION`). |
| `plot_budget_ts_sea`, `plot_budget_ts_zonal` | TS budget with land/sea mask, zonal-mean TS budget (any `REGION`). |
| `plot_budget_regmean` | Regional-mean monthly TS budget and per-model/season CSVs (any `REGION`). |
| `plot_budget_seb`, `plot_budget_ts_mean` | Surface energy budget closure; ensemble time-mean TS budget (v3 LE). |
| `plot_aa_{stat_ts,stat_ac,stat_pdf,2d}_{65n,65n_driftcorr,50n}` | Arctic amplification statistics and maps. |

AA variants: `_65n` (poleward of 65°N), `_65n_driftcorr` (65°N, drift-corrected
E3SM), `_50n` (poleward of 50°N). Run `process_aa_*` before `plot_aa_*`.

## Shell scripts

Every script starts with a `DATA_DIR` parameter that defaults to this machine's
location and can be overridden from the environment:

```bash
DATA_DIR=/other/machine/data bash shell_script/run_process_era5.bash
sbatch --export=ALL,DATA_DIR=/other/machine/data shell_script/sbatch_mpas_analysis_ice_2024.bash
```

- `run_process_*.bash`, `derive.bash`: climatologies with derived surface-energy
  fields into `${DATA_DIR}/climo/` (default: the v3_polar_analysis data tree,
  which v3LE_paper/data/climo links to). The `_v3le` variants of `HadSST` and
  `e3sm_hist` are the v3 LE versions (CMIP6 SST/ice forcing 1869–2022; other
  member list). `run_process_e3sm_hist_v3le.bash` still has a debugging `exit`.
- `run_gen_map.sh`: ERA5/MERRA2 → 1° regrid maps in `${DATA_DIR}/map/`.
- `sbatch_mpas_analysis_*.bash`, `sbatch_process_ice_*.bash`: MPAS sea-ice/ocean
  time series and regridded sea-ice fields (default: the v3LE_paper data tree).

## Checks

After editing notebooks or scripts, run from `work_package/`:

```bash
python tools/check_package.py polar_analysis
```

It fails on syntax errors, notebooks without exactly one `parameters` cell right
after the shared setup cell, absolute `/lcrc` paths outside the parameters cell,
invalid `config/regions.json` entries, and shell scripts without a `DATA_DIR`
parameter. It warns about unused parameters and input paths that do not exist
on the current machine.

## Shared modules

- `scripts/exp_info.py`: run, region, variable, and budget catalogs.
- `scripts/regions.py`: per-region settings from `config/regions.json`.
- `scripts/data_utils.py`: run/variable metadata, region definitions, longitude
  normalization, and land/ice masks (formerly `function.py`).

## Known input issues on this machine

- `process_ts_nmsk`, `plot_budget_seb`, and `plot_budget_ts_mean` read annual
  `v3.LR.historical.enNN.area_mean_ts.185001-202412.nc` files from
  `DATA_DIR/Arctic/`. The originals (in the removed `large_ensemble/` tree) are
  gone and no notebook here produces them. Checks against the values printed in
  the notebooks' saved outputs show they were annual, cos-latitude-weighted
  65–90°N means: atmosphere fields from `data/global/*.monthly_base.*` reproduce
  the original SST exactly, but the sea-ice fields (`SICONC`, `SIMASS`,
  `SITHICK`, and the `_imsk` masks built from them) came from regridded
  MPAS-SeaIce output (`shell_script/sbatch_process_ice_*.bash`), which must be
  regenerated before the files can be rebuilt.
- `diag_data/` has no `v3.LR.historical.Arctic.DJF.ensemble_stats.*.nc`, so
  `plot_ts_sea_ice_imsk` stops at DJF until `process_ts_sea_ice_imsk` regenerates it.

## Not carried over

Only in the `analysis_v0/` archives: `trash/` folders, `*_bak`/`*-Copy1`
notebooks, the superseded `arctic_65N_v0/` AA configuration, the older
`energy_balance` notebooks and `thermo.py`/`constants.py`/`energy_budget.py`
helpers, the per-region notebook copies from `v3_polar_analysis`, its `scripts/`
v2 energy files, and previously generated figures.
