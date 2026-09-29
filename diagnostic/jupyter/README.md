# Diagnostic notebooks

The notebooks are grouped by diagnostic area and numbered in their recommended
workflow order. They are interactive drivers: shared implementation belongs in
`../util/`, shared configuration in `../configs/`, and unattended entry points
in `../script/`.

## Top-level path parameters

Every notebook starts with a code cell tagged `parameters`. Edit only these three
paths when relocating a workflow:

```python
INPUT_DATA_ROOT = Path("/compyfs/zhan391/v3_dart_cda_scratch")
WORK_DIR = Path("/compyfs/zhan391/e3sm_dart_analysis")
DIAGNOSTIC_OUTPUT_ROOT = Path("/compyfs/www/zhan391/e3sm_dart/diag_out")
```

Outputs are derived consistently as
`DIAGNOSTIC_OUTPUT_ROOT/data/<workflow>/...` for NetCDF and intermediate data,
and `DIAGNOSTIC_OUTPUT_ROOT/figure/<workflow>/...` for figures. Input products
are resolved beneath `INPUT_DATA_ROOT`; repository modules are loaded from
`WORK_DIR/diagnostic`. The regional forecast time-series notebook additionally
exposes `OBSERVATION_ROOT` because its IMERG and NOAA products live in the shared
observation archive.

## `analysis_atm_da/`

1. `01_obs_distribution.ipynb` - map DART observation distributions.
2. `02_obs_diag_check.ipynb` - compare two observation-diagnostic products.
3. `03_obs_diag_compare.ipynb` - compare bias, RMSE, total spread, spread/RMSE, and rejection-rate time series across experiments.
4. `04_obs_profile_diagnostics.ipynb` - compare bias, RMSE, total spread, spread/RMSE, and rejection-rate vertical profiles.
5. `05_obs_multilevel_diagnostics.ipynb` - analyze multiple pressure layers and regions.
6. `06_compute_and_plot_ensemble_ranges.ipynb` - cache and plot ensemble-member spatial minimum, mean, and maximum time series from NetCDF, Zarr, or Intake-ESM inputs.

## `analysis_atm_init/`

1. `01_compute_and_plot_horizontal_maps.ipynb` - compute and plot bias, RMSE, and spread maps for surface variables or selected 3-D levels.
2. `02_compute_and_plot_vertical_cross_sections.ipynb` - compute and plot latitude–level or longitude–level bias, RMSE, and spread sections for 3-D variables.
3. `03_compute_and_plot_analysis_increments.ipynb` - compute cached posterior-minus-prior maps and time-mean vertical increment/posterior-spread sections without a reference dataset.
4. `04_check_restart_stability.ipynb` - validate native EAM restarts across members and timestamps using cached range, non-finite, spatial-flag, and inter-snapshot change diagnostics.

Each notebook is independently runnable and lists its available variables and metrics. Request-keyed NetCDF files are stored under `DIAGNOSTIC_OUTPUT_ROOT/data/analysis_atm_init/atmosphere_metrics/`; existing files are overwritten only when `force_compute = True`.

## `analysis_atm_bias/`

1. `01_compute_and_plot_bias_metrics.ipynb` - compute cached atmospheric bias, RMSE, spread, and CRPS maps for DA analysis periods.

## `analysis_lnd_bias/`

1. `01_compute_and_plot_bias_metrics.ipynb` - compute cached land bias, RMSE, spread, and CRPS maps for DA analysis periods.

## `fcst_atm_bias/`

1. `01_compute_and_plot_bias_metrics.ipynb` - compute the same atmospheric metrics for forecast periods only.

## `fcst_lnd_bias/`

1. `01_compute_and_plot_bias_metrics.ipynb` - compute the same land metrics for forecast periods only.

The four bias workflows are intentionally separated by period and component. Each notebook reads, processes, caches, and plots independently; cached products are reused unless `force_compute = True`.


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
3. `analysis_lnd_eval/01_compute_and_plot_ensemble_metrics.ipynb` - land maps and scalar ensemble verification metrics.

Forecast-period verification:

1. `fcst_atm_eval/01_compute_and_plot_ensemble_metrics.ipynb` - forecast atmospheric ensemble verification.
2. `fcst_atm_eval/02_compute_and_plot_surface_energy_budget.ipynb` - forecast surface-energy-budget terms and regional ensemble time series.
3. `fcst_lnd_eval/01_compute_and_plot_ensemble_metrics.ipynb` - forecast land ensemble verification.

Each notebook reads, processes, caches, and plots independently. The notes list every supported variable and metric; existing request-keyed NetCDF products are reused unless `force_compute = True`.

## `analysis_lnd_init/`

1. `01_soil_moisture_distribution.ipynb` - soil-moisture bias and spread maps.
2. `02_surface_variable_distribution.ipynb` - PRECT/TREFHT bias and spread.
3. `03_restart_processing.ipynb` - cached restart aggregation and optional conservative regridding by calling reusable `util/` implementations directly.
4. `04_process_initial_condition_errors.ipynb` - process ensemble soil-moisture bias, spread, and significance.
5. `05_plot_initial_condition_error_maps.ipynb` - plot processed initial-condition error maps.
6. `06_plot_soil_point_timeseries.ipynb` - extract and plot soil-moisture time series at selected points.

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
