# E3SM v3 large-ensemble polar analysis

This package contains the data-processing and plotting workflows for the polar
diagnostics of the E3SM v3 large-ensemble (v3 LE) paper (Arctic-focused so far):
regional time series, mean-state biases, surface-temperature budgets, and Arctic
amplification (AA). The notebooks under
`jupyter/` hold the experiment configuration and analysis. The shared run,
variable, and region catalogs live under `scripts/`.

The workflows were collected from the `1_*`–`5_*` directories of
`/lcrc/group/e3sm/ac.szhang/acme_scratch/e3sm_project/v3_le_paper`. Those
originals, together with `6_pcmdi_diag/` and `script/`, are archived unchanged in
`/lcrc/group/e3sm/public_html/diagnostic_output/ac.szhang/v3LE_paper/analysis_v0`.

## Layout

```text
polar_analysis/
├── config/paths.json  the only place the v3LE_paper root is written
├── jupyter/           notebooks (process_* writes data, plot_* makes figures)
├── scripts/           shared modules; paths.py exposes the v3LE locations
└── shell_script/      shell/Slurm preprocessing; paths.sh exposes the same locations
```

All data and diagnostic output goes under the fixed paper location
`/lcrc/group/e3sm/public_html/diagnostic_output/ac.szhang/v3LE_paper`:

```text
v3LE_paper/
├── data/          V3LE_DATA_DIR: processed inputs (former v3_le_paper/data)
│   ├── global/            monthly_base 3D fields per member/dataset (~410 GB)
│   ├── global_driftcorr/  drift-corrected monthly_base fields (~860 GB)
│   ├── regmn/             region_means time series (global, Arctic, ...)
│   ├── regmn_driftcorr/   region_means of drift-corrected fields (+ backup/)
│   ├── sea_ice/           per-member mpas_ts sea-ice/atm time series (imsk)
│   ├── mpas_ts/           MPAS-Analysis mpassi/mpaso time series (shell_script/sbatch_mpas_*)
│   ├── pcmdi_diags/       symlinks to e3sm-pcmdi-le/climo member directories
│   ├── climo/             created by shell_script/run_process_*.bash (bias/budget inputs)
│   └── logs/              missing-file lists from the MPAS batch scripts
├── diag_data/     V3LE_DIAG_DIR: diagnostic NetCDF/JSON products (AA_*, AA2_*,
│                  *.ensemble_stats.*.nc, mean_bias/budget <region>/ outputs, ...)
├── figures/       V3LE_FIG_ROOT: time_series_imsk/, time_series_nmsk/, mean_bias/,
│                  budget/, arctic_amplification/<suffix>/, pcmdi/
└── analysis_v0/   archive of the original v3_le_paper workflow folders
```

The pipeline is `data/` → `process_*` → `diag_data/` → `plot_*` → `figures/`.
`diag_data/` replaces the former `v3_le_paper/figure_data` (same flat file names).

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

Outputs go to `V3LE_DATA_DIR` (`global/`, `regmn/`, and the `*_driftcorr` variants).

### 2. Regional time series

| Notebook | Purpose |
| --- | --- |
| `process_ts_atm_imsk.ipynb` | Atmospheric anomaly statistics with the ice mask (`imsk`). |
| `process_ts_sea_ice_imsk.ipynb` | Sea-ice anomaly statistics with the ice mask. |
| `plot_ts_sea_ice_imsk.ipynb` | Seasonal sea-ice time-series panels. Figures go to `figures/time_series_imsk/`. |
| `process_ts_nmsk.ipynb` | Annual anomaly mean/spread panels without the mask (`nmsk`). Figures go to `figures/time_series_nmsk/<region>/`. |

### 3. Mean-state bias

- `plot_bias_map.ipynb`: regional bias maps against observations.
- `plot_bias_zonal.ipynb`: zonal-mean biases.

Figures go to `figures/mean_bias/<region>/`.

### 4. Surface-temperature budget

- `plot_budget_seb.ipynb`: surface energy budget terms and closure check.
- `plot_budget_ts_mean.ipynb`: time-mean TS budget decomposition.
- `plot_budget_ts_sea.ipynb`: TS budget with the land/sea mask.
- `plot_budget_ts_zonal.ipynb`: zonal-mean TS budget.

These notebooks use `RUN_CATALOG`, `REGION_CATALOG`, `VARIABLE_CATALOG`, and
`Budget_CATALOG` from `scripts/exp_info.py`. Figures go to `figures/budget/<region>/`.

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

Figures go to `figures/arctic_amplification/<suffix>/`.

## Running a notebook

1. Start Jupyter from anywhere inside `polar_analysis/` with an environment
   that provides xarray, xcdat, matplotlib, and cartopy (for example, E3SM
   Unified).
2. Run the first code cell. Every notebook has the same setup cell: it finds the
   package from the working directory, adds `scripts/` to `sys.path`, and imports
   `V3LE_DATA_DIR`, `V3LE_DIAG_DIR`, `V3LE_FIG_ROOT`, `fig_dir`, and `require`
   from `scripts/paths.py`; plotting notebooks then set
   `FIG_DIR_ROOT = fig_dir("<workflow>")`.
3. Review the configuration cell (`DATA_DIR`, `OUT_DIR`, run lists, periods)
   before running the rest. Raw model and observation inputs are
   absolute LCRC paths outside this package (for example `CVDP_RGD`).

## Paths and checks

- To move the v3LE output tree, edit `v3le_root` in `config/paths.json`. For a
  single session, set `V3LE_ROOT=/other/v3LE_paper` before starting Jupyter or
  submitting a script; notebooks and shell scripts both honor it.
- `paths.require(path, ...)` raises one error listing every missing input.
- Shell scripts `source shell_script/paths.sh`. Under `sbatch`, they locate it
  through `scontrol show job` because Slurm runs a spooled copy of the script.
- After editing notebooks or scripts, run from `work_package/`:

  ```bash
  python tools/check_package.py polar_analysis
  ```

  It fails on syntax errors, hard-coded v3LE roots, path names used before the
  `from paths import`, differing setup cells, and shell scripts that do not
  source `paths.sh`. It warns about input paths that do not exist; these
  warnings match the "Known input-path issues" below.

## Shared modules

- `scripts/paths.py`: v3LE locations from `config/paths.json` (see above).
- `scripts/data_utils.py`: run metadata (`RUN_DICT`), variable metadata and
  contour levels, region definitions, longitude normalization, and land/ice
  masks. It was `function.py` in the original folders.
- `scripts/exp_info.py`: dataclass-based run, region, variable, and budget
  catalogs for the bias and budget notebooks. It merges the two earlier
  copies: the `3_mean_bias` version plus the `Arctic` region from
  `4_budget_analysis`.

## Not carried over

These items are only in the `analysis_v0/` archive:

- `trash/` folders, `*_bak` and `*-Copy1` notebooks, and the superseded
  `arctic_65N_v0/` AA configuration.
- The `function.py` copy in `3_mean_bias/` and `4_budget_analysis/`. No notebook
  imports it, and it depends on the `thermo`/`constants` modules in
  `4_budget_analysis/trash/`.
- Generated PDFs, `__pycache__/`, and `.ipynb_checkpoints/`.
- `script/`: older versions of the time-series and budget workflows
  (`plot_time_series.ipynb`, `energy_balance.ipynb`, `function.py`,
  `energy_budget.py`) and the unused `thermo.py`/`constants.py` thermodynamics
  helpers.

## Shell scripts

`shell_script/` holds the shell preprocessing that feeds `V3LE_DATA_DIR`:

- `run_process_{e3sm_hist,e3sm_amip,era5,noaa2c,HadSST,era5.e3smdiag}.bash`:
  build 30-year/2001–2018 climatologies with derived surface-energy-budget
  fields into `data/climo/` (the `plot_bias_*`/`plot_budget_ts_{sea,zonal}`
  inputs). Scripts that regrid expect a `map/map_721x1440_to_180x360_conserve.nc`
  file relative to the working directory. `run_process_e3sm_hist.bash` still has
  a debugging `exit` inside its loop.
- `sbatch_mpas_analysis_{ice_2024,ice_2050,ocn_2050}.bash`: extract MPAS-Analysis
  sea-ice/ocean time series into `data/mpas_ts/`; missing inputs are appended to
  `data/logs/missing.txt`.
- `sbatch_process_ice_{2024,2050}.bash`: regrid MPAS-SeaIce history to 1° under
  `data/v3.LR.historical/`.

## Known input-path issues

These notebooks still read inputs from project directories that no longer exist,
and their data was not found in `V3LE_DATA_DIR`. Update their
`TOP_DIR`/`top_dir`/`top_path` before running them (the checker lists them):

- `large_ensemble/`: the `*_50n` AA notebooks (left unmapped because `data/regmn`
  stores one "Arctic" region, used by the 65°N notebooks), `process_ts_nmsk`
  (`data/Arctic/*.area_mean_ts.*`), `process_seaice_ts` (`data/seaice`), and
  `plot_budget_seb`/`plot_budget_ts_mean` (`data/<region>`).

`data/climo/` does not exist yet; run the `shell_script/run_process_*.bash` scripts
before the bias and budget notebooks. `diag_data/` has no
`v3.LR.historical.Arctic.DJF.ensemble_stats.*.nc`, so `plot_ts_sea_ice_imsk`
stops at DJF until `process_ts_sea_ice_imsk` regenerates it.
