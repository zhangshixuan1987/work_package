# E3SM v3 large-ensemble Arctic analysis

This package contains the data-processing and plotting workflows for the E3SM v3
large-ensemble (v3 LE) Arctic paper: regional time series, mean-state biases,
surface-temperature budgets, and Arctic amplification (AA). The notebooks under
`jupyter/` hold the experiment configuration and analysis. The shared run,
variable, and region catalogs live under `scripts/`.

The workflows were collected from the `1_*`–`5_*` directories of
`/lcrc/group/e3sm/ac.szhang/acme_scratch/e3sm_project/v3_le_paper`.

## Layout

```text
arctic_analysis/
├── jupyter/   notebooks (process_* writes data, plot_* makes figures)
└── scripts/   shared modules imported by the notebooks
```

All figures are written under the fixed paper location
`/lcrc/group/e3sm/public_html/diagnostic_output/ac.szhang/v3LE_paper`
(`V3LE_FIG_ROOT` in each notebook's first code cell), one subdirectory per
workflow: `time_series_imsk/`, `time_series_nmsk/`, `mean_bias/`, `budget/`, and
`arctic_amplification/<suffix>/`.

## Notebook organization

Notebook names follow:

```text
<stage>_<diagnostic>_<detail>_<variant>.ipynb
```

`process_*` notebooks read raw model or observational output and write
intermediate NetCDF/JSON products. `plot_*` notebooks read those products and
make figures. Run each group in the order listed below.

### 1. Data processing

| Step | Notebooks | Purpose |
| --- | --- | --- |
| 1 | `process_3d_{e3sm,era5,hadisst,noaa20c,noaav5,picontrol}_data.ipynb`, `process_3d_picontrol_trend.ipynb`, `process_seaice_ts.ipynb` | Regrid/subset 3D fields, the piControl trend, and sea-ice time series. |
| 2 | `process_ts_{e3sm,era5,noaa20c}.ipynb` | Build regional-mean time series. |
| 3 | `process_3d_e3sm_driftcorr.ipynb`, `process_ts_e3sm_driftcorr.ipynb` | Apply the piControl drift correction to E3SM fields and time series. |

Outputs go to `v3_le_paper/data` (set by `top_dir`/`out_dir` in each notebook).

### 2. Regional time series

| Notebook | Purpose |
| --- | --- |
| `process_ts_atm_imsk.ipynb` | Atmospheric anomaly statistics with the ice mask (`imsk`). |
| `process_ts_sea_ice_imsk.ipynb` | Sea-ice anomaly statistics with the ice mask. |
| `plot_ts_sea_ice_imsk.ipynb` | Seasonal sea-ice time-series panels. Figures go to `time_series_imsk/`. |
| `process_ts_nmsk.ipynb` | Annual anomaly mean/spread panels without the mask (`nmsk`). Figures go to `time_series_nmsk/<region>/`. |

### 3. Mean-state bias

- `plot_bias_map.ipynb`: regional bias maps against observations.
- `plot_bias_zonal.ipynb`: zonal-mean biases.

Figures go to `mean_bias/<region>/`.

### 4. Surface-temperature budget

- `plot_budget_seb.ipynb`: surface energy budget terms and closure check.
- `plot_budget_ts_mean.ipynb`: time-mean TS budget decomposition.
- `plot_budget_ts_sea.ipynb`: TS budget with the land/sea mask.
- `plot_budget_ts_zonal.ipynb`: zonal-mean TS budget.

These notebooks use `RUN_CATALOG`, `REGION_CATALOG`, `VARIABLE_CATALOG`, and
`Budget_CATALOG` from `scripts/exp_info.py`. Figures go to `budget/<region>/`.

### 5. Arctic amplification

Each AA configuration has the same seven notebooks, distinguished by suffix:

| Suffix | Configuration (original folder) |
| --- | --- |
| `_65n` | Arctic poleward of 65°N (`arctic_65N`) |
| `_65n_driftcorr` | 65°N with drift-corrected E3SM data (`arctic_65N_drifcorr`) |
| `_50n` | Arctic poleward of 50°N (`arctic_50N`) |

Run order within a configuration:

1. `process_aa_ts_<suffix>.ipynb`: AA time series.
2. `process_aa_ac_<suffix>.ipynb`: AA annual cycle.
3. `process_aa_2d_<suffix>.ipynb`: AA spatial fields.
4. `plot_aa_stat_ts_<suffix>.ipynb`, `plot_aa_stat_ac_<suffix>.ipynb`,
   `plot_aa_stat_pdf_<suffix>.ipynb`, `plot_aa_2d_<suffix>.ipynb`: figures.

Figures go to `arctic_amplification/<suffix>/`.

## Running a notebook

1. Start Jupyter from anywhere inside `arctic_analysis/` with an environment
   that provides xarray, xcdat, matplotlib, and cartopy (for example, E3SM
   Unified).
2. Run the first code cell. It locates `PROJECT_ROOT` (the directory that
   contains `scripts/`), adds `scripts/` to `sys.path`, and, where needed,
   defines `V3LE_FIG_ROOT` and the workflow's `FIG_DIR_ROOT`.
3. Review the configuration cell (`TOP_DIR`, `DATA_DIR`, `OUT_DIR`, run lists,
   periods) before running the rest. Input and intermediate data paths are
   absolute LCRC paths, mostly under `v3_le_paper/` or `large_ensemble/`.

## Shared modules

- `scripts/data_utils.py`: run metadata (`RUN_DICT`), variable metadata and
  contour levels, region definitions, longitude normalization, and land/ice
  masks. It was `function.py` in the original folders.
- `scripts/exp_info.py`: dataclass-based run, region, variable, and budget
  catalogs for the bias and budget notebooks. It merges the two earlier
  copies: the `3_mean_bias` version plus the `Arctic` region from
  `4_budget_analysis`.

## Not carried over

These items remain in the original `v3_le_paper/` directories:

- `trash/` folders, `*_bak` and `*-Copy1` notebooks, and the superseded
  `arctic_65N_v0/` AA configuration.
- The `function.py` copy in `3_mean_bias/` and `4_budget_analysis/`. No notebook
  imports it, and it depends on the `thermo`/`constants` modules in
  `4_budget_analysis/trash/`.
- Generated PDFs, `__pycache__/`, and `.ipynb_checkpoints/`.
