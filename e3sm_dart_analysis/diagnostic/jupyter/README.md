# Diagnostic notebooks

The notebooks are grouped by diagnostic area and numbered in their recommended
workflow order. They are interactive drivers: shared implementation belongs in
`../util/`, shared configuration in `../configs/`, and unattended entry points
in `../script/`.

## Top-level path parameters

Each notebook starts with a path-parameter code cell (tagged `parameters` for papermill). Input data depend on
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

## Layout

Directories are prefixed by period and component: `analysis_*` covers the DA
(analysis) windows, S0 = Dec 2011 and S1 = May 2012; `fcst_*` covers the
hindcasts started from them. `atm`, `lnd`, and `lac` (land-atmosphere coupling)
name the component; `*_init` describes the initial/analysis state, `*_eval`
verifies it against observations. Notebooks are numbered in workflow order
within a directory; notebooks sharing a number are independent workflows.

`hind_*` covers the 2012-2013 AMIP hindcasts: eight continuous E3SMv3 LR AMIP
runs started on 2012-01-01 from different atmosphere/land initial states
(CTRL from the piControl; DART, ERA5, and mixed atmosphere/land states;
`configs/hind_experiment_config.py`). They run on Perlmutter, so their
parameter cells default to

```python
INPUT_DATA_ROOT = Path("/global/cfs/cdirs/m4849/zhan391/eam_hindcast")
WORK_DIR = Path("/global/cfs/cdirs/m4849/zhan391/work_package/e3sm_dart_analysis")
DIAGNOSTIC_OUTPUT_ROOT = Path("/global/cfs/cdirs/e3sm/www/zhan391/e3sm_dart/diag_dart_2026")
OBSERVATION_ROOT = Path("/global/cfs/cdirs/e3sm/zhan391/data")  # reference-data archive
```

The archive holds every reference these notebooks read: ERA5 6-hourly 1-degree
fields (`ERA5/1x1/6hourly`, 1979-2019, averaged to days), ERA5 single-level
monthly fields (`ERA5/1x1/monthly`), the E3SM Diags time series
(`obs_for_e3sm_diags/time-series`: NOAA OLR and IMERG daily; ERA5, ERA-Interim,
NOAA-20C monthly pressure-level fields), and the BoM RMM index (`MJO_INDEX`).

with outputs under `DIAGNOSTIC_OUTPUT_ROOT/{data,figure}/<workflow>/<topic>/`.
Model and reference readers (`util/hind_model_data.py`,
`util/hind_observation_data.py`) open files lazily with Dask, label daily means
by the day they cover and monthly means by year-month, and subset before
computing (the ERA5 pressure-level references are 0.25-degree, ~75 GB files).

## `analysis_atm_da/`

DART observation-space and state-space diagnostics of the assimilation.

1. `01_obs_distribution.ipynb` - map DART observation distributions.
2. `02_obs_diag_check.ipynb` - compare two observation-diagnostic products.
3. `03_obs_diag_compare.ipynb` - compare bias, RMSE, total spread, spread/RMSE, and rejection-rate time series across experiments.
4. `04_obs_multilevel_diagnostics.ipynb` - analyze multiple pressure layers and regions.
5. `05_prior_innovation_profiles.ipynb` - separate systematic and random components of prior observation-space mismatch.
6. `06_state_update_profiles.ipynb` - compute long-window pressure-profile increments, normalized increments, and posterior/prior spread ratios from DART ensemble summaries.
7. `07_state_update_spatial_patterns.ipynb` - compute representative increment maps and spatial increment/spread-ratio sections as a one-day case-study companion to notebook 06.
8. `08_obs_profile_diagnostics.ipynb` - compare bias, RMSE, total spread, spread/RMSE, and rejection-rate vertical profiles.

## `analysis_atm_init/`

1. `01_compute_and_plot_horizontal_maps.ipynb` - compute and plot bias, RMSE, and spread maps for surface variables or selected 3-D levels.
2. `02_compute_and_plot_vertical_cross_sections.ipynb` - compute and plot latitude–level or longitude–level bias, RMSE, and spread sections for 3-D variables.
3. `04_check_restart_stability.ipynb` - validate native EAM restarts across members and timestamps using cached range, non-finite, spatial-flag, and inter-snapshot change diagnostics (`util/atm_restart_stability.py`).

Each notebook is independently runnable and lists its available variables and metrics. Request-keyed NetCDF files are stored under `DIAGNOSTIC_OUTPUT_ROOT/data/analysis_atm_init/atmosphere_metrics/`; existing files are overwritten only when `force_compute = True`.

## `analysis_atm_bias/`

1. `01_compute_and_plot_bias_metrics.ipynb` - compute cached atmospheric bias, RMSE, spread, and CRPS maps for DA analysis periods (`util/bias_metrics.py`).

The bias workflows are separated by period and component: `analysis_atm_bias/01` and `analysis_lnd_init/01` (DA analysis periods), `fcst_atm_eval/04_compute_and_plot_bias_metrics` and `fcst_lnd_eval/01_compute_and_plot_bias_metrics` (forecast periods). Each reads, processes, caches, and plots independently; cached products are reused unless `force_compute = True`.

## `analysis_atm_eval/`

Atmospheric verification over the DA windows (workflow `analysis_atm_model_evaluation`).

1. `01_compute_and_plot_scalar_metrics.ipynb` - global (cos-latitude weighted) ensemble verification metrics and rank histograms per case (`util/model_evaluation.py`).
2. `02_compute_and_plot_map_metrics.ipynb` - maps of the same metrics (bias, RMSE, spread, CRPS, ...) per case, variable, and metric, from the same verification step as notebook 01.
3. `03_compute_and_plot_surface_energy_budget.ipynb` - member-level surface-energy-budget caches and regional ensemble summaries (`util/atm_surface_energy_budget.py`).
4. `04_compute_da_scorecard_metrics.ipynb` - score variables against references at every verification time of the DA window and write the scorecard tables: S0 % RMSE change of DART10/20/40 vs CTRL10, DART40 S0 vs S1 normalised RMSE, and the update impact (`util/da_scorecard.py`).
5. `05_plot_da_scorecards.ipynb` - scorecard heat maps from the notebook 04 tables; no model data are read.

## `analysis_lnd_init/`

Land state during the DA (analysis) periods. Implementation lives in `../util/`; notebooks are thin drivers.

1. `01_compute_and_plot_bias_metrics.ipynb` - cached land bias, RMSE, spread, and CRPS maps for the DA analysis periods (`util/bias_metrics.py`).
2. `02_soil_moisture_distribution.ipynb` - top-5 cm soil-moisture bias and spread vs ESA CCI at the end of the DA window; writes the processed snapshot used by `analysis_lnd_eval/03` (`util/initial_lnd_soil_moisture.py`).
3. `03_restart_processing.ipynb` - ELM restart column-to-gridcell aggregation and optional conservative regridding (`util/initial_lnd_restart_aggregation.py`, `util/initial_lnd_restart_regridding.py`).

## `analysis_lnd_eval/`

Land verification over the DA windows.

1. `01_land_forcing_scorecard.ipynb` - land-only DA scorecard on the 64 scorecard cycles: % RMSE change of DART10/20/40 vs CTRL10 for TS, TREFHT, net surface shortwave (FSNS), and atmospheric forcing (PRECT, FSDS, FLDS) (`util/da_scorecard.py`).
2. `02_soil_moisture_daily_skill.ipynb` - daily top-5 cm soil-moisture anomaly correlation and ubRMSE of DART40-S0/S1 and the CTRF10-S0/S1 free-run (no-DA) control against ESA CCI over the scorecard windows (`util/lnd_sm_daily_skill.py`).
3. `03_plot_atm_da_land_impact.ipynb` - spread ratio, mean shift, and error change of DART vs CTRL10 top-5 cm soil moisture at the end of the S0 window; reads the snapshot from `analysis_lnd_init/02` (`util/initial_lnd_da_impact.py`).
4. `04_surface_variable_distribution.ipynb` - PRECT/TREFHT ensemble bias, spread, and significance vs GPCP/ERA5 on selected January forecast dates (`util/initial_lnd_surface_variables.py`).

## `fcst_atm_eval/`

Atmospheric verification of the hindcasts. Every notebook loops over `cases` (S0: Jan-Feb 2012, S1: Jun-Jul 2012 by default) and can run on its own; 02 and 03 share caches.

1. `01_compute_and_plot_ens_ts.ipynb` - two modes: `daily_performance` (6-hourly member fields aggregated to 24-h lead bins, then daily ensemble-mean RMSE, finite-ensemble-corrected spread/RMSE, bias, and fair CRPS of PSL, TREFHT, PRECT, FLUT against ERA5, GPM IMERG, and NOAA OLR, with a 10-member-subset fairness test; `util/forecast_daily_skill.py`) and `shock` (6-hourly initial-shock tendency of PSL, TREFHT, TS; `util/forecast_ensemble_ts.py`) (`fcst_atm_ensemble_ts`).
2. `02_compute_and_plot_scalar_metrics.ipynb` - global ensemble verification metrics, summary table, and rank histograms (`fcst_atm_model_evaluation`; `util/model_evaluation.py`).
3. `03_compute_and_plot_map_metrics.ipynb` - bias, RMSE, spread, CRPS, and coverage maps from the same verification step (`fcst_atm_model_evaluation`; `util/model_evaluation.py`).
4. `04_compute_and_plot_bias_metrics.ipynb` - bias, RMSE, spread, and CRPS maps with bootstrap significance hatching (`fcst_atm_bias`; `util/bias_metrics.py`).
5. `05_compute_and_plot_surface_energy_budget.ipynb` - member-level surface-energy-budget terms and regional ensemble time series (`fcst_atm_model_evaluation/surface_energy_budget`; `util/atm_surface_energy_budget.py`).
6. `06_compute_and_plot_decay_timescales.ipynb` - shock-minus-baseline anomalies and exponential decay timescales (`fcst_atm_initial_shock`; `util/forecast_initial_shock.py`).
7. `07_compute_and_plot_window_scorecard.ipynb` - heat maps of anomaly PCC, MSSS vs the climatology forecast (optionally with the mean bias removed), and fair-CRPS skill vs the climatological ensemble per biweekly lead window (four triangles per cell) for atmosphere (PSL; U, V, T, Q, Z at 200/500/850 hPa; TREFHT, PRECT, FLUT) and land variables, CTRL/CAPT/DART40 in S0 and S1, against ERA5, GPM IMERG, NOAA OLR, and ESA CCI with a shared 2001-2017 climatology (`fcst_window_scorecard`; `util/forecast_window_scorecard.py`).

## `fcst_lnd_eval/`

Land verification of the hindcasts.

1. `01_compute_and_plot_bias_metrics.ipynb` - forecast-period land bias, RMSE, spread, and CRPS maps (`fcst_lnd_bias`; `util/bias_metrics.py`).
2. `01_compute_and_plot_ensemble_metrics.ipynb` - extended forecast land ensemble verification (`fcst_lnd_model_evaluation`; `util/model_evaluation.py`).
3. `02_plot_soil_point_timeseries.ipynb` - cached member spaghetti or mean ± spread of forecast soil moisture against ESA CCI at selected points (`util/initial_lnd_point_timeseries.py`).
4. `03_process_initial_condition_errors.ipynb` - top-10 cm soil-moisture bias, spread, and significance vs ESA CCI on the forecast start dates (January/June 2012) (`util/initial_lnd_file_collection.py`, `util/initial_lnd_soil_error.py`).
5. `04_plot_initial_condition_error_maps.ipynb` - 2x2 bias/spread maps of two experiments from notebook 03.

## `fcst_lac_eval/`

Land-atmosphere coupling in the hindcasts.

1. `01_compute_and_plot_lag_coupling.ipynb` - lagged correlation of land and atmosphere ensemble departures (SM5-LHFLX, SM5-PRECT, TS-LHFLX): significance maps and regional lag curves for S0/S1 hindcasts (`fcst_lac_feedback`; `util/lnd_atm_feedback.py`).
2. `02_compute_correlation_coupling.ipynb` - LHFLX-TREFHT squared-correlation coupling.
3. `03_compute_dirmeyer_tci.ipynb` - model and observational Dirmeyer TCI.
4. `04_plot_tci_differences.ipynb` - original model TCI comparison plots.
5. `05_plot_coupling_diagnostics.ipynb` - TCI/TCC differences and observational products.
6. `06_compute_and_plot_ensemble_cross_correlation.ipynb` - compute, cache, and plot ensemble atmosphere-land restart cross correlations on collocated grids or from one land point to all atmosphere columns.

## `fcst_atm_mjo/`

1. `01_hovmoller_phase_composites.ipynb` - configurable PRECT/FLUT Hovmöller diagrams, RMM phase space, and phase composites.
2. `02_rmm_diagnostics.ipynb` - RMM of every S0/S1 hindcast member and of the observations, all projected onto the observed WH04 EOFs of `hind_atm_mjo/03` (aligned with the BoM index; previous-120-day mean from the observations before each start date); ensemble-mean error, spread, amplitude and phase errors per lead day, and bivariate COR/RMSE per lead window (`fcst_atm_mjo`; `util/mjo_rmm_projection.py`).

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

## `hind_atm_mjo/`

MJO in the AMIP hindcasts (same diagnostics as `fcst_atm_mjo/`, for continuous runs; `util/mjo_hovmoller_coupling.py`, `util/mjo_rmm_projection.py`).

1. `01_compute_and_plot_hovmoller.ipynb` - 10S-10N Hovmöller diagrams of OLR and precipitation vs NOAA OLR and IMERG, with pattern correlation and RMSE.
2. `02_compute_and_plot_convection_wind_coupling.ipynb` - 20-100 day band-passed OLR-U850 lag correlation over 60-150E and its minimum and lag.
3. `03_compute_and_plot_rmm_skill.ipynb` - observed WH04 EOFs of OLR/U850/U200 (1986-2017; previous-120-day mean removed), order and signs aligned with the BoM index, and their checks (exact self-projection, orthogonality, agreement with BoM); the EOF file is shared with `fcst_atm_mjo/02`. Then the projection of observations and experiments on those EOFs (each run's previous-120-day mean starts from the observations before 2012-01-01) and bivariate COR/RMSE (Rashid et al. 2011) in running windows and lead bins. Reads ERA5 U850/U200 from `ERA5/1x1/6hourly`.

## `hind_atm_eval/`

Stratosphere, extratropical modes, and tropical winds in the AMIP hindcasts.

1. `01_compute_and_plot_qbo_wind.ipynb` - 5S-5N zonal-mean U time-height sections and amplitudes for every experiment, ERA5, and ERA-Interim (`util/qbo_diagnostics.py`); sections, differences, amplitudes, RMSE/PCC skill bars over both reanalyses with block-bootstrap CI, persistence (`util/qbo_plotting.py`).
2. `02_compute_and_plot_qbo_temperature.ipynb` - the same for T: anomalies, amplitudes, vertical coherence, and T′-wind-shear (thermal-wind) coupling (uses the U sections of notebook 01).
3. `03_compute_and_plot_ssw_diagnostics.ipynb` - monthly U(55-65N), polar-cap Z and T, biases and window RMSEs vs ERA5 and ERA-Interim (`util/ssw_diagnostics.py`); 2012-13 SSW bias and running RMSE, RMSE attribution by initial-state component, evolution, RMSE profiles (`util/ssw_plotting.py`).
4. `04_compute_and_plot_nao_index.ipynb` - NAO EOF base from NOAA-20C SLP (NDJFM, 1985-2014) and its check, then the daily NAO index of every experiment and ERA5 projected on it; 14-day lead-window means and correlations (`util/nao_index.py`).
5. `05_compute_and_plot_tropical_wind_skill.ipynb` - daily 20S-20N U850/U200 spatial ACC and RMSE vs ERA5 by lead day, window means with block-bootstrap CI (`util/tropical_wind_skill.py`).

Every `hind_atm_*` notebook covers one metric in the same order: read input, derive, save (then reload every source from disk), scores or checks, plot. Saved diagnostics are reused: input fields are opened lazily, and only for sources whose output file is missing, unless `FORCE_COMPUTE = True` (a rebuilt NAO or RMM EOF file recomputes every index projected on it). The long training records behind the NAO and RMM EOF checks are read only when the EOFs are built or `RECHECK = True`. In the derive step the expensive reductions of all sources are computed in one Dask graph (`util.dask_helpers.compute_all`); `DASK_WORKERS` (default 8) caps the threads, which bounds memory on many-core login nodes; `USE_DISTRIBUTED = True` runs the graph on a local `dask.distributed` cluster of `DASK_WORKERS` processes instead, which read NetCDF files in parallel.

These notebooks replace the former `v3_dart_paper/` notebooks (`mjo_analysis`, `qbo_analysis`, `ssw_analysis`, `nao_analysis`, `tropicwind_analysis`), whose duplicated classes and `ExperimentUtil.py` copies are now the modules above.

## Running

Notebooks load repository modules from `WORK_DIR/diagnostic`, so they can be
opened from any directory; set `WORK_DIR` in the parameter cell to use another
checkout.
