# Diagnostics

This directory contains diagnostic notebooks, helper scripts, legacy NCL workflows,
and generated figures used for E3SM-DART analysis.

## Repository Layout Convention

Going forward on this branch, keep workflow code organized by role:

- `jupyter/` contains driver notebooks only. Notebooks should configure a run,
  call reusable modules, inspect intermediate results, and make figures; avoid
  adding importable helper modules or large reusable classes here.
- `util/` contains reusable helper/library modules, such as
  `util/dask_helpers.py`, `util/dart_observation_diagnostics.py`, and data-reader or analysis
  classes shared by multiple notebooks or scripts.
- `script/` contains runnable workflow scripts and CLI-style batch drivers.
  Use this for command-line entry points, batch processing wrappers, and
  repeatable production runs that should not require opening a notebook.
- `configs/` contains experiment dictionaries and configuration-building
  helpers shared by notebooks, utilities, and scripts.

If a notebook cell grows into reusable logic, move that logic into `util/`. If a
workflow needs to run unattended or from a scheduler, put that entry point in
`script/` and keep the notebook as an interactive driver or demonstration.

## Diagnostic Notebooks

Every maintained notebook begins with a tagged parameter cell defining `INPUT_DATA_ROOT`, `WORK_DIR`, and `DIAGNOSTIC_OUTPUT_ROOT`. Data products are written under `DIAGNOSTIC_OUTPUT_ROOT/data/<workflow>/`; figures are written under `DIAGNOSTIC_OUTPUT_ROOT/figure/<workflow>/`. The repository defaults are `/compyfs/zhan391/v3_dart_cda_scratch`, `/compyfs/zhan391/work_package/e3sm_dart_analysis`, and `/compyfs/www/zhan391/e3sm_dart/diag_dart_2026`, respectively.

Active land-atmosphere workflow drivers live in `jupyter/analysis_lac/`:

1. `01_prepare_metrics.ipynb` - configure and compute lag-correlation and spatial-feedback metrics.
2. `02_compute_correlation_coupling.ipynb` - compute the LHFLX-TREFHT squared-correlation coupling diagnostic.
3. `03_compute_dirmeyer_tci.ipynb` - compute model and observational Dirmeyer TCI for top-5-cm, full-column, or 10-cm soil moisture.
4. `04_plot_tci_differences.ipynb` - plot the original model TCI comparison workflow.
5. `05_plot_coupling_diagnostics.ipynb` - plot model TCI/TCC differences, observational TCI maps, and regional time series.

Analysis-initialization diagnostics live in `jupyter/analysis_atm_init/`:

1. `01_compute_and_plot_horizontal_maps.ipynb` handles surface and selected-level bias, RMSE, and spread maps.
2. `02_compute_and_plot_vertical_cross_sections.ipynb` handles latitude–level or longitude–level bias, RMSE, and spread sections for 3-D variables.
3. `03_compute_and_plot_analysis_increments.ipynb` handles cached DA-only posterior-minus-prior maps and time-mean vertical increment/posterior-spread sections without a reference dataset.

Each notebook reads raw fields only for missing or forced products. Existing request-keyed NetCDF products are reused unless `force_compute = True`.

Bias verification is split by period and component:

1. `jupyter/analysis_atm_bias/01_compute_and_plot_bias_metrics.ipynb` - atmospheric DA-analysis-period verification.
2. `jupyter/analysis_lnd_bias/01_compute_and_plot_bias_metrics.ipynb` - land DA-analysis-period verification.
3. `jupyter/fcst_atm_bias/01_compute_and_plot_bias_metrics.ipynb` - atmospheric forecast-period verification.
4. `jupyter/fcst_lnd_bias/01_compute_and_plot_bias_metrics.ipynb` - land forecast-period verification.

All four use `util/bias_metrics.py`, keep analysis and forecast caches separate, and reuse existing NetCDF products unless `force_compute = True`.
The broader legacy variable metadata is preserved in `configs/bias_variable_catalog.py`.

Active data-assimilation diagnostic workflow drivers live in `jupyter/analysis_atm_da/`:

1. `01_obs_distribution.ipynb` - inspect and plot DART obs_seq observation distributions.
2. `02_obs_diag_check.ipynb` - compare overlapping obs_diag products numerically and visually.
3. `03_obs_diag_compare.ipynb` - compare obs_diag bias, RMSE, total spread, spread/RMSE, and rejection-rate time series; switch modes with `PLOT_MODE`.
4. `04_obs_profile_diagnostics.ipynb` - plot bias, RMSE, total spread, spread/RMSE, and rejection-rate profiles; switch modes with `PROFILE_MODE`.
5. `05_obs_multilevel_diagnostics.ipynb` - compare diagnostics across multiple pressure layers and regions.

Initial-land diagnostics follow the same notebook-driver layout:

1. `jupyter/analysis_lnd_init/01_soil_moisture_distribution.ipynb` drives the soil-moisture map workflow.
2. `jupyter/analysis_lnd_init/02_surface_variable_distribution.ipynb` processes and plots PRECT/TREFHT bias and spread.
3. `jupyter/analysis_lnd_init/03_restart_processing.ipynb` drives column-to-gridcell restart aggregation and optional conservative regridding.
4. `jupyter/analysis_lnd_init/04_process_initial_condition_errors.ipynb` computes ensemble soil-moisture bias, spread, and significance products.
5. `jupyter/analysis_lnd_init/05_plot_initial_condition_error_maps.ipynb` plots those processed initial-condition products.
6. `jupyter/analysis_lnd_init/06_plot_soil_point_timeseries.ipynb` extracts cached or direct point time series.
7. `configs/initial_land_experiment_config.py`, `configs/initial_land_experiment_registry.py`, and `configs/initial_land_observation_registry.py` store paths and experiment/observation definitions.
8. `util/initial_lnd_model_data.py`, `util/initial_lnd_observation_data.py`, `util/initial_lnd_restart_aggregation.py`, and `util/initial_lnd_restart_regridding.py` hold reusable processing.
9. `script/run_process_land_rest.bash` drives the complete batch workflow; `script/initial_land/aggregate_restart.py` and `script/initial_land/regrid_restart.py` expose its individual stages.
10. Regrid reference grids and map files are expected under `/compyfs/zhan391/v3_dart_cda_scratch/reference/regrid_maps/` by default.
11. Initial-land support files such as `dzsoi_elm.nc` and `landmask_1x1.nc` are generated under `/compyfs/zhan391/v3_dart_cda_scratch/reference/lnd_sea_mask/` when missing.
12. Generated initial-land NetCDF/log outputs live under `/compyfs/www/zhan391/e3sm_dart/diag_dart_2026/data/analysis_lnd_init/`; generated figures live under `/compyfs/www/zhan391/e3sm_dart/diag_dart_2026/figure/analysis_lnd_init/`.

The former `7_mjo_analysis/` workflow is consolidated under `jupyter/analysis_mjo/`:
`01_hovmoller_phase_composites.ipynb` handles configurable PRECT/FLUT
Hovmöller and phase composites, while `02_rmm_diagnostics.ipynb` handles
ensemble RMM/EOF diagnostics. Reusable implementations live in `util/`, and
all generated MJO data and figures use the shared diagnostic output roots. See
`jupyter/README.md` for the complete notebook index and workflow order.


S2S forecast verification is consolidated into two parameterized workflows:

1. `jupyter/fcst_s2s_pcc/` computes regional PCC/ACC and RMSE, horizontal ACC maps, and plots cached regional skill.
2. `jupyter/fcst_s2s_tcc/` computes lead-window TCC/RMSE/bias maps, reduces them regionally, and plots cached map or regional products.
3. Shared experiment, observation, metric, and plotting implementations live in `configs/s2s_config.py` and `util/s2s_*.py`.
4. `jupyter/fcst_atm_timeseries/` provides cached regional PRECT/FLUT ensemble time series, and `jupyter/fcst_atm_initial_shock/` diagnoses experiment-minus-baseline decay timescales.

The TCC workflow also includes the Dec2011 `da` analysis selection inherited from `6_s2s_skills/`; weekly, biweekly, 15-day, and monthly caches remain separate. The former seasonal and ensemble-slice notebook copies in `analysis_s2s_pcc/` and `analysis_s2s_tcc/` are superseded by notebook parameters. Cached NetCDF products are reused unless `force_compute = True`.

The legacy `cross_correlation/` workflow is consolidated into `jupyter/analysis_lac/06_compute_and_plot_ensemble_cross_correlation.ipynb`, and the generic `ensemble_plot/` notebook is now `jupyter/analysis_atm_da/06_compute_and_plot_ensemble_ranges.ipynb`, which compares cached CTRL and DART hindcast spatial-range ensembles. The duplicate `fcst_initial_shock/` notebook is superseded by the cached `jupyter/fcst_atm_initial_shock/` workflow.

Atmospheric restart stability checking is maintained in `jupyter/analysis_atm_init/04_check_restart_stability.ipynb` with reusable processing in `util/atm_restart_stability.py`. It supersedes the broken one-off `sanity_check/` notebook and caches request-specific range, non-finite, spatial, and transition diagnostics.

Model evaluation is split by period and component under `jupyter/analysis_*_eval/` and `jupyter/fcst_*_eval/`. Atmospheric workflows include the separate surface-energy-budget diagnostic. Extended reusable metrics live in `util/model_evaluation.py`; cached energy-budget processing lives in `util/atm_surface_energy_budget.py`. These workflows supersede the duplicated `3_model_evaluation/` and `model_evaluation/` trees and reuse cached products unless `force_compute = True`.

## Configs

`configs/` contains experiment dictionaries and configuration-building helpers
shared by the diagnostic notebooks. Keep reusable importable modules here; keep
standalone commands in `script/`.

## NCL Scripts

Actively maintained NCL workflow bundles live in `script/ncl/`. Historical
paper-figure sources are stored outside the repository at
`/compyfs/www/zhan391/e3sm_dart/diag_dart_2026/paper/ncl_paper_figures/`. The bundles retain
their internal directory structure because several NCL programs use relative
`load "./..."` statements.

- `script/ncl/obs_distribution/` contains scripts from the old `1_obs_distribution/` workflow.
- `script/ncl/obs_diagnostics/time_series/` contains helper scripts from the old `2_obs_diagnostics/time_series/` workflow.
- `script/ncl/obs_diagnostics/profile/` contains helper scripts from the old `2_obs_diagnostics/profile/` workflow.
- `script/ncl/rmse_*/` contains the retained RMSE, bias, and increment workflows.

Run these programs with `script/run_process_ncl.bash`. The launcher selects the correct
working directory and moves newly generated figures from both maintained and
historical scripts to the fixed
`DIAGNOSTIC_OUTPUT_ROOT/figure/ncl/` tree. Generated NetCDF files and logs go to
`DIAGNOSTIC_OUTPUT_ROOT/data/ncl/`. The default output root is
`/compyfs/www/zhan391/e3sm_dart/diag_dart_2026`.

For example:

```bash
diagnostic/script/run_process_ncl.bash \
  diagnostic/script/ncl/rmse_profile/plot_vertical_preshgt.ncl
```

The RMSE/bias profile backup directories document four regional scripts:

- Southern Hemisphere: `plot_rmse_bias_profil_eamv0_vs_eam_v1_T_SH.ncl`
- Northern Hemisphere: `plot_rmse_bias_profil_eamv0_vs_eam_v1_T_NH.ncl`
- Tropics: `plot_rmse_bias_profil_eamv0_vs_eam_v1_T_TP.ncl`
- North America: `plot_rmse_bias_profil_eamv0_vs_eam_v1_T_NA.ncl`

## Helper Scripts

`script/` contains helper scripts for running or batching diagnostic workflows.
