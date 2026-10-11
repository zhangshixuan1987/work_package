# Diagnostic notebooks

The notebooks are grouped by diagnostic area and numbered in their recommended
workflow order. They are interactive drivers: shared implementation belongs in
`../util/`, shared configuration in `../configs/`, and unattended entry points
in `../script/`.

## Top-level path parameters

Every notebook starts with a code cell tagged `parameters`. Input data depend on
the machine and analysis case, so these paths are notebook parameters rather
than values in `../configs/`. Edit them per run, or pass them with papermill
(`papermill nb.ipynb out.ipynb -p INPUT_DATA_ROOT /other/scratch`):

```python
INPUT_DATA_ROOT = Path("/compyfs/zhan391/v3_dart_cda_scratch")
WORK_DIR = Path("/compyfs/zhan391/work_package/e3sm_dart_analysis")
DIAGNOSTIC_OUTPUT_ROOT = Path("/compyfs/www/zhan391/e3sm_dart/diag_dart_2026")
# only where used:
OBSERVATION_ROOT = Path("/compyfs/zhan391/acme_init/Observations")   # observations
LAND_SUPPORT_DIR = Path("/compyfs/zhan391/acme_init/lnd_sea_mask")   # land mask, soil-layer thickness
REGRID_MAP_DIR = Path("/compyfs/zhan391/acme_init/map_file")         # regrid grids and maps
```

The parameter cell holds plain values only. `sys.path` setup and
`workflow_output_dirs(...)` live in the following cell, so overridden
parameters are honored. Shared modules accept paths as arguments
(e.g. `get_experiment_dict(..., data_root=INPUT_DATA_ROOT)`) instead of
reading fixed defaults.

Outputs are derived consistently as
`DIAGNOSTIC_OUTPUT_ROOT/data/<workflow>/...` for NetCDF and intermediate data,
and `DIAGNOSTIC_OUTPUT_ROOT/figure/<workflow>/...` for figures. Model output
is resolved beneath `INPUT_DATA_ROOT`; observations beneath `OBSERVATION_ROOT`
(the shared archive, one copy per product; ESA CCI soil moisture is in
`ESA_CCI/H20SOI/`); land masks and soil-layer thicknesses in `LAND_SUPPORT_DIR`;
regrid grids and maps in `REGRID_MAP_DIR`. Module defaults for these locations
are in `../configs/input_paths.py` (overridable with `E3SM_DART_OBS_ROOT`,
`E3SM_DART_LAND_REFERENCE_DIR`, `E3SM_DART_REGRID_MAP_DIR`). Repository modules
are loaded from `WORK_DIR/diagnostic`.

## `analysis_atm_da/`

1. `01_obs_distribution.ipynb` - map DART observation distributions.
2. `02_obs_diag_check.ipynb` - compare two observation-diagnostic products.
3. `03_obs_diag_compare.ipynb` - compare bias, RMSE, total spread, spread/RMSE, and rejection-rate time series across experiments.
4. `04_obs_multilevel_diagnostics.ipynb` - analyze multiple pressure layers and regions.
5. `05_prior_innovation_profiles.ipynb` - separate systematic and random components of prior observation-space mismatch.
6. `06_state_update_profiles.ipynb` - compute long-window pressure-profile increments, normalized increments, and posterior/prior spread ratios from DART ensemble summaries.
7. `07_state_update_spatial_patterns.ipynb` - compute representative increment maps and spatial increment/spread-ratio sections as a one-day case-study companion to notebook 06.
8. `08_obs_profile_diagnostics.ipynb` - compare bias, RMSE, total spread, spread/RMSE, and rejection-rate vertical profiles.
9. `09_compute_and_plot_ens_ranges.ipynb` - cache and compare CTRL/DART ensemble-member spatial minima, means, maxima, envelopes, and ensemble means.

## `analysis_atm_init/`

1. `01_compute_and_plot_horizontal_maps.ipynb` - compute and plot bias, RMSE, and spread maps for surface variables or selected 3-D levels.
2. `02_compute_and_plot_vertical_cross_sections.ipynb` - compute and plot latitude–level or longitude–level bias, RMSE, and spread sections for 3-D variables.
3. `04_check_restart_stability.ipynb` - validate native EAM restarts across members and timestamps using cached range, non-finite, spatial-flag, and inter-snapshot change diagnostics.

Each notebook is independently runnable and lists its available variables and metrics. Request-keyed NetCDF files are stored under `DIAGNOSTIC_OUTPUT_ROOT/data/analysis_atm_init/atmosphere_metrics/`; existing files are overwritten only when `force_compute = True`.

## `analysis_atm_bias/`

1. `01_compute_and_plot_bias_metrics.ipynb` - compute cached atmospheric bias, RMSE, spread, and CRPS maps for DA analysis periods.

## `fcst_atm_bias/`

1. `01_compute_and_plot_bias_metrics.ipynb` - compute the same atmospheric metrics for forecast periods only.

The land counterparts are `analysis_lnd_init/01_compute_and_plot_bias_metrics.ipynb` (DA analysis periods) and `fcst_lnd_eval/01_compute_and_plot_bias_metrics.ipynb` (forecast periods). The four bias workflows are intentionally separated by period and component. Each notebook reads, processes, caches, and plots independently; cached products are reused unless `force_compute = True`.


## `fcst_atm_timeseries/`

1. `01_compute_and_plot_regional_timeseries.ipynb` - compute, cache, and plot regional daily PRECT and FLUT forecast ensemble time series against observations.

## `fcst_atm_initial_shock/`

1. `01_compute_and_plot_decay_timescales.ipynb` - compute, cache, and plot experiment-minus-baseline initial-shock anomalies and exponential decay timescales. This supersedes both legacy initial-shock notebook copies.

## `fcst_s2s_pcc/`

1. `01_compute_regional_skill.ipynb` - compute and cache lead-dependent regional PCC/ACC and RMSE metrics for selected groups, runs, observations, and variables.
2. `02_compute_acc_maps.ipynb` - compute and cache horizontal ACC maps for configurable lead windows.
3. `03_plot_regional_skill.ipynb` - plot cached regional PCC/ACC and RMSE summaries without rereading raw forecasts.

## `fcst_s2s_tcc/`

1. `01_compute_window_metrics.ipynb` - compute and cache per-member and ensemble-mean TCC, RMSE, and bias maps by lead window for forecast groups or the Dec2011 DA analysis group.
2. `02_compute_regional_metrics.ipynb` - reduce cached window maps to weighted regional ensemble statistics.
3. `03_plot_window_maps.ipynb` - plot cached TCC, RMSE, or bias maps across experiments and lead windows.
4. `04_plot_regional_skill.ipynb` - plot cached regional TCC/RMSE ensemble summaries.

Both S2S workflows expose January/June forecast groups, observation-variable selections, regions, and metrics as notebook parameters. The TCC workflow also accepts the Dec2011 analysis group with the `da` run. Weekly, biweekly, 15-day, and monthly TCC caches are kept in separate subdirectories. Existing NetCDF products are reused unless `force_compute = True`.

## Model evaluation workflows

Analysis-period verification:

1. `analysis_atm_eval/01_compute_and_plot_ensemble_metrics.ipynb` - atmospheric maps, scalar verification metrics, rank histograms, and optional bootstrap confidence intervals.
2. `analysis_atm_eval/02_compute_and_plot_surface_energy_budget.ipynb` - atmospheric surface-energy-budget terms and regional ensemble time series.
3. `analysis_lnd_eval/01_land_forcing_scorecard.ipynb` - land-only DA scorecard on the 64 scorecard cycles: % RMSE change of DART10/20/40 vs CTRL10 for TS, TREFHT, surface fluxes (LHFLX, SHFLX, FSNS), and atmospheric forcing (PRECT, FSDS, FLDS) (`util/da_scorecard.py`).
4. `analysis_lnd_eval/02_soil_moisture_daily_skill.ipynb` - daily top-5 cm soil-moisture anomaly correlation and ubRMSE of DART40-S0/S1 against ESA CCI over the scorecard windows (`util/lnd_sm_daily_skill.py`).
5. `analysis_lnd_eval/03_plot_atm_da_land_impact.ipynb` - spread ratio, mean shift, and error change of DART vs CTRL10 top-5 cm soil moisture at the end of the S0 window; reads the snapshot from `analysis_lnd_init/02` (`util/initial_lnd_da_impact.py`).
6. `analysis_lnd_eval/04_surface_variable_distribution.ipynb` - PRECT/TREFHT ensemble bias, spread, and significance vs GPCP/ERA5 on selected January forecast dates (`util/initial_lnd_surface_variables.py`).

Forecast-period verification:

1. `fcst_atm_eval/01_compute_and_plot_ensemble_metrics.ipynb` - forecast atmospheric ensemble verification.
2. `fcst_atm_eval/02_compute_and_plot_surface_energy_budget.ipynb` - forecast surface-energy-budget terms and regional ensemble time series.
3. `fcst_lnd_eval/01_compute_and_plot_ensemble_metrics.ipynb` - forecast land ensemble verification.
4. `fcst_lnd_eval/02_plot_soil_point_timeseries.ipynb` - cached member spaghetti or mean ± spread of forecast soil moisture against ESA CCI at selected points (`util/initial_lnd_point_timeseries.py`).
5. `fcst_lnd_eval/03_process_initial_condition_errors.ipynb` - top-10 cm soil-moisture bias, spread, and significance vs ESA CCI on the forecast start dates (January/June 2012) (`util/initial_lnd_file_collection.py`, `util/initial_lnd_soil_error.py`).
6. `fcst_lnd_eval/04_plot_initial_condition_error_maps.ipynb` - 2x2 bias/spread maps of two experiments from notebook 03.

Each notebook reads, processes, caches, and plots independently. The notes list every supported variable and metric; existing request-keyed NetCDF products are reused unless `force_compute = True`.

## `analysis_lnd_init/`

Land state during the DA (analysis) periods. Implementation lives in `../util/`; notebooks are thin drivers. Verification over the DA windows is in `analysis_lnd_eval/`, forecast-period land diagnostics in `fcst_lnd_eval/` (see Model evaluation workflows).

1. `01_compute_and_plot_bias_metrics.ipynb` - cached land bias, RMSE, spread, and CRPS maps for the DA analysis periods (`util/bias_metrics.py`).
2. `02_soil_moisture_distribution.ipynb` - top-5 cm soil-moisture bias and spread vs ESA CCI at the end of the DA window; writes the processed snapshot used by `analysis_lnd_eval/03` (`util/initial_lnd_soil_moisture.py`).
3. `03_restart_processing.ipynb` - ELM restart column-to-gridcell aggregation and optional conservative regridding (`util/initial_lnd_restart_aggregation.py`, `util/initial_lnd_restart_regridding.py`).

## `analysis_lac/`

1. `01_prepare_metrics.ipynb` - lag-correlation and spatial-feedback metrics.
2. `02_compute_correlation_coupling.ipynb` - LHFLX-TREFHT squared-correlation coupling.
3. `03_compute_dirmeyer_tci.ipynb` - model and observational Dirmeyer TCI.
4. `04_plot_tci_differences.ipynb` - original model TCI comparison plots.
5. `05_plot_coupling_diagnostics.ipynb` - TCI/TCC differences and observational products.
6. `06_compute_and_plot_ensemble_cross_correlation.ipynb` - compute, cache, and plot ensemble atmosphere-land restart cross correlations on collocated grids or from one land point to all atmosphere columns.

## `analysis_mjo/`

1. `01_hovmoller_phase_composites.ipynb` - configurable PRECT/FLUT Hovmöller diagrams, RMM phase space, and phase composites.
2. `02_rmm_diagnostics.ipynb` - ensemble EOF/RMM projections, skill scores, phase-space plots, and amplitude distributions.

Run notebooks from the repository root when practical. Each maintained notebook
also discovers the repository through `E3SM_DART_ANALYSIS_ROOT` or its parent
directories, so opening it from its topic directory is supported.
