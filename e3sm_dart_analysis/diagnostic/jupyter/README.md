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

The archive holds every reference these notebooks read: ERA5 1-degree daily means
(`ERA5/1x1/daily`, 1979-2026, from `ERA5/1x1/6hourly`), ERA5 single-level monthly
fields (`ERA5/1x1/monthly`), NOAA OLR CDR daily (`NOAA_CDR/1x1/daily`), GPCP V3.3
daily (`GPCPv3.3/1x1/daily`) and IMERG daily (`NASA_GPM/1x1/daily`), the BoM RMM index
(`MJO_INDEX`), and, from the E3SM Diags time series (`obs_for_e3sm_diags/time-series`),
only the ERA5, ERA-Interim and NOAA-20C monthly pressure-level fields.

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

Atmospheric verification of the hindcasts. Every notebook loops over `cases` (S0: Jan-Feb 2012, S1: Jun-Jul 2012 by default) and can run on its own; the two 03 notebooks share caches.

1. `01_compute_and_plot_scorecard.ipynb` - heat maps of anomaly PCC, MSSS vs the climatology forecast (optionally with the mean bias removed), and fair-CRPS skill vs the climatological ensemble per biweekly lead window (four triangles per cell) for surface (PSL, TREFHT, TS, 10-m wind, PRECT, FLUT, surface fluxes), upper-air (U, V, T, Q, Z at 200/500/700/850/925 hPa), and land variables, CTRL/CAPT/DART40 in S0 and S1 in one figure (`full_list`: every variable; `sub_list`: key variables), against ERA5, GPM IMERG, NOAA OLR, and ESA CCI with a shared 2001-2017 climatology (`fcst_window_scorecard`; `util/forecast_window_scorecard.py`).
2. `02_compute_and_plot_ens_ts.ipynb` - daily performance (6-hourly member fields aggregated to 24-h lead bins, then daily ensemble-mean RMSE, finite-ensemble-corrected spread/RMSE, bias, and fair CRPS of PSL, TREFHT, PRECT, FLUT against ERA5, GPCP, and NOAA OLR, with a 10-member-subset fairness test; `util/forecast_daily_skill.py`) and initial shock (6-hourly tendency of the same variables; `util/forecast_ensemble_ts.py`) (`fcst_atm_ensemble_ts`).
3. `03_compute_and_plot_scalar_metrics.ipynb` - global ensemble verification metrics, summary table, and rank histograms (`fcst_atm_model_evaluation`; `util/model_evaluation.py`).
4. `03_compute_and_plot_map_metrics.ipynb` - bias, RMSE, spread, CRPS, and coverage maps from the same verification step (`fcst_atm_model_evaluation`; `util/model_evaluation.py`).
5. `04_compute_and_plot_bias_metrics.ipynb` - bias, RMSE, spread, and CRPS maps with bootstrap significance hatching (`fcst_atm_bias`; `util/bias_metrics.py`).
6. `05_compute_and_plot_surface_energy_budget.ipynb` - member-level surface-energy-budget terms and regional ensemble time series (`fcst_atm_model_evaluation/surface_energy_budget`; `util/atm_surface_energy_budget.py`).
7. `06_compute_and_plot_decay_timescales.ipynb` - shock-minus-baseline anomalies and exponential decay timescales (`fcst_atm_initial_shock`; `util/forecast_initial_shock.py`).

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

## `hind_lead_eval/`

Data and methods, with the reasons for each choice (sampling, climatologies, bias treatments, MJO winds, limitations): [`hind_lead_eval/METHODS.md`](hind_lead_eval/METHODS.md).

Initialization-strategy hindcasts on Perlmutter (`INPUT_DATA_ROOT/F20TR_ne30pg2_r05_IcoswISC30E3r5_<CONFIG>_<YYYYMMDD>00/<ENxx>`): five configurations, CLIM_BASE (free-running climate state), HYB_BASE (ERA5-based atmosphere + DATM-forced land), DA_ATM (DART atmosphere), DA_LND (DART land), and DA_COUPLED (DART atmosphere and land). Each runs 20 members for 2 months from 4 starts per season, 16-31 Dec 2011 and 16-31 May 2012, 5 days apart (`configs/fcst_init_experiment_config.py`). Readers: `util/fcst_init_data.py`, which handles lead-day alignment, sampling (00Z states for pressure levels) and units as in the window scorecard, and reads ERA5, GPCP V3.3 daily (precipitation), NOAA OLR and ESA CCI references and their climatology, and the monthly climatology of the free-running E3SMv3 AMIP runs (`v3.LR.amip_0101`, `0151`) used for the forecast anomalies. The ERA5 1-degree archive is built independently of this repository by `/global/cfs/cdirs/e3sm/zhan391/data/ERA5/scripts/` (README there) (6-hourly files from the hourly 0.25-degree ERA5, and `1x1/daily`, their daily means: exact for fluxes, trapezoid rule over 00Z..00Z for states). Daily-mean fields are verified against `1x1/daily`, the 00Z pressure-level fields against `1x1/6hourly`. Scores: `util/fcst_init_skill.py`. Figures: `util/fcst_init_plotting.py`. The model time series must first be regenerated with `post_process/submit_regen_f20tr_ts.sbatch`; the zppy `ncclimo --split` files stop after 4 days. Each notebook only uses configuration-seasons whose members all report `regen_ts_<stream>.status` OK (`REQUIRE_ALL_MEMBERS`). Every notebook's parameter cell holds `CONFIGS`, a dictionary from configuration key (`CLIM_BASE`, `HYB_BASE`, ...) to raw case name, the case of a start being `INPUT_DATA_ROOT/<case name>_<YYYYMMDD>00`; `set_case_names(CONFIGS)` (`configs/fcst_init_experiment_config.py`) makes every reader use these names, and saved products record the case they were built from, so renaming a case recomputes them.

1. `01_compute_and_plot_lead_skill.ipynb` - daily ensemble-mean RMSE, bias, spread, spread/RMSE, and fair CRPS by lead day (1-56), pooled over the starts of a season, for each region.
2. `02_compute_and_plot_s2s_weekly_skill.ipynb` - weekly S2S anomaly skill following the ESP-Lab `s2s_skill` protocol (weeks 1-8, ACC with a t-test, ΔACC vs `BASELINE`). Both climatologies are built the same way for every variable: calendar-month means over `CLIM_YEARS` (2001-2014), interpolated linearly between mid-month values to the verification days. Observed anomalies use the monthly means of the daily reference (ERA5 `1x1/daily`, GPCP V3.3 daily, NOAA OLR), so the climatology comes from the same product as the verifying days; forecast anomalies are computed two ways (`BIASES`): against the observed climatology (`"none"`: initialization drift and model mean bias count as error, the full-field forecast) and against the AMIP monthly climatology (`"amip"`: model mean bias removed, the S2S convention; 4 starts per season are too few for a hindcast climatology). Neither removes the drift; `s2s_bias_impact_*` shows their difference (the effect of the model mean bias on skill, negative while the forecast is still near the observed state, positive once it has drifted to the model climate). Both are all-hours means, so the 00Z anomalies of the pressure-level fields carry each system's 00Z departure from the daily mean (ERA5 RMS about 6 m Z500, 0.5 K T850). The leave-one-start-out drift removal (`"loo"`) is kept for comparison only, since its starts verify on overlapping days. Outputs: regional lead-time skill curves in the layout of ESP-Lab `2a_regional_acc_skill_ts` (area-mean gridpoint ACC and regional nRMSE by lead week; rows Dec, May and all starts, columns ACC and nRMSE), regional-index skill and series in the layout of `3a_sst_skill_ts` (ACC with filled markers where significant and nRMSE of the regional-mean anomaly; the index by start date at fixed lead windows, for `INDEX_REGIONS`), regional pattern ACC and MSSS, plus gridpoint ACC over the starts of the available seasons (8 when both are ready). ΔACC uses only the seasons that both the configuration and `BASELINE` have. Every variable and region is scored and saved; only the selection is plotted: a named setup of `PLOT_SETUPS` (`PLOT_SETUP`: bias treatments, seasons, figure types), `PLOT_VARIABLES` and `PLOT_REGIONS`.
   `02b_compute_and_plot_s2s_leadtime_skill_maps.ipynb` - lead-time skill maps of two-week means (weeks 1-2, 3-4, 5-6, 7-8), the S2S analogue of the ESP-Lab `1a`/`1b` lead-time ACC and RMSE maps, from the weekly anomalies saved by 02: gridpoint ACC over the starts (stippled above the 95% critical value) and nRMSE = RMSE / climatology-forecast RMSE (below 1 beats climatology), Dec and May starts stacked in one figure (optionally both pooled), configurations in columns; differences from `BASELINE` (CLIM_BASE) and a regional summary (`s2s_window_summary.csv`).
3. `03_compute_and_plot_land_skill_and_coupling.ipynb` - top-5-cm soil moisture vs ESA CCI and TSA vs ERA5 over land (lead-day and weekly, anomalies as in 02), and ensemble-based soil-moisture memory and SM-TSA/SM-LHFLX coupling per configuration.
4. `04_compute_and_plot_mjo_rmm_skill.ipynb` - RMM projection of every member on the observed EOFs of `hind_atm_mjo/03` (previous-120-day mean taken from the observations), then bivariate COR/RMSE, spread, and amplitude error by lead pooled over starts, plus lead-window skill.
5. `05a`-`05c` - atmospheric and land DART-initialization benefits for 2-m temperature over 30-60N land. A 2 x 2 of the atmosphere and land initial states: HYB_BASE (E00, CAPT-like: ERA5 + DATM-forced land) and DA_COUPLED (E11, coupled DART) are the consistent initializations, DA_ATM (E10) and DA_LND (E01) the mixed ones. Benefits are pattern-ACC differences centred on DA_COUPLED: atmosphere E11 - E01, land E11 - E10, and the land-atmosphere consistency (interaction) E11 - E10 - E01 + E00 (members exchangeable, not paired across experiments; `util/init_benefit.py`), with a paired bootstrap over starts and member indices. DA_LND starts every member from the unperturbed ERA5 atmosphere (no initial atmospheric spread). `05a_compute_and_plot_init_benefit_curves.ipynb` the three-panel figure on daily running means (schematic hypothesis, ERA5 verification, held-out-member verification); `05b_compute_and_plot_init_benefit_weekly.ipynb` weekly benefits with intervals and absolute skill; `05c_compute_and_plot_init_benefit_maps.ipynb` gridpoint-ACC benefit maps, weeks 3-4 (one figure per bias treatment: DJF, JJA and both seasons by atmosphere and land benefit).
6. `06a`-`06h` - case study of the summer 2012 central U.S. heat wave (May 2012 starts): `06a` observed event, thresholds and simulation coverage (design, event definition and reference states); `06b` experiment QC; `06c` start-resolved regional trajectories; `06d` fixed-event-window metrics and maps; `06e` S2S lead-window scores; `06f` event probability and Brier score; `06g` unpaired ensemble contrasts; `06h` sensitivity of the event definition (`util/heatwave_case.py`).
7. `07a`-`07d` - case study of the January 2012 Pacific Northwest storm (Dec 2011 starts): `07a` regional series (design and event definition); `07b` member storm metrics; `07c` event maps of precipitation and temperature anomaly; `07d` storm environment (IVT, Z500, PSL) (`util/winter_storm_case.py`).

The case notebooks 05-07 are split into small notebooks that each read the data they need, compute, save and plot, independently of each other. Their case definitions are dataclass defaults in `configs/hind_lead_eval_cases.py` (one per series, changed per notebook with `CASE_OVERRIDES`), and their shared products come from `util/init_benefit_workflow.py`, `util/heatwave_workflow.py` and `util/winter_storm_workflow.py`. Every product is saved with the settings it was computed from and read back when they are unchanged (`util/hind_case_io.py`); `FORCE_COMPUTE` recomputes a notebook's results, and `FORCE_INPUTS` also re-reads the model output and observations.

## `hind_atm_mjo/`

MJO in the AMIP hindcasts (same diagnostics as `fcst_atm_mjo/`, for continuous runs; `util/mjo_hovmoller_coupling.py`, `util/mjo_rmm_projection.py`).

1. `01_compute_and_plot_hovmoller.ipynb` - 10S-10N Hovmöller diagrams of OLR and precipitation vs NOAA OLR (CDR) and GPCP V3.3, with pattern correlation and RMSE.
2. `02_compute_and_plot_convection_wind_coupling.ipynb` - 20-100 day band-passed OLR-U850 lag correlation over 60-150E and its minimum and lag.
3. `03_compute_and_plot_rmm_skill.ipynb` - observed WH04 EOFs of OLR/U850/U200 (1986-2017; previous-120-day mean removed), order and signs aligned with the BoM index, and their checks (exact self-projection, orthogonality, agreement with BoM); the EOF file is shared with `fcst_atm_mjo/02`. Then the projection of observations and experiments on those EOFs (each run's previous-120-day mean starts from the observations before 2012-01-01) and bivariate COR/RMSE (Rashid et al. 2011) in running windows and lead bins. Reads ERA5 U850/U200 from `ERA5/1x1/6hourly`.

## `hind_atm_eval/`

Stratosphere, extratropical modes, and tropical winds in the AMIP hindcasts.

1. `01_compute_and_plot_qbo.ipynb` - 5S-5N zonal-mean U and T time-height sections and amplitudes for every experiment, ERA5, and ERA-Interim (`util/qbo_diagnostics.py`). Wind: sections, differences, amplitudes, RMSE/PCC skill bars over both reanalyses with block-bootstrap CI, persistence, memory and onset timing (`util/qbo_plotting.py`). Temperature: anomalies, amplitudes, vertical coherence, and T′-wind-shear (thermal-wind) coupling; the paper figure `qbo_thermal_wind_consistency.png` splits each run's temperature error into the part balanced by its shear error and the unbalanced rest, by lead month and atmospheric initial state.
2. `02_compute_and_plot_ssw_diagnostics.ipynb` - monthly U(55-65N), polar-cap Z and T, biases and window RMSEs vs ERA5 and ERA-Interim (`util/ssw_diagnostics.py`); 2012-13 SSW bias and running RMSE, RMSE attribution by initial-state component, evolution, RMSE profiles (`util/ssw_plotting.py`); uses the QBO U sections of notebook 01.
3. `03_compute_and_plot_nao_index.ipynb` - NAO EOF base from NOAA-20C SLP (NDJFM, 1985-2014) and its check, then the daily NAO index of every experiment and ERA5 projected on it; 14-day lead-window means and correlations (`util/nao_index.py`).
4. `04_compute_and_plot_tropical_wind_skill.ipynb` - daily 20S-20N U850/U200 spatial ACC and RMSE vs ERA5 by lead day, window means with block-bootstrap CI (`util/tropical_wind_skill.py`).

Experiment set: every `hind_atm_*` parameter cell maps run keys to figure labels in `EXPERIMENTS_FULL`, all eight runs of the atmosphere x land factorial, and in `EXPERIMENTS_FCST`, the three runs whose initial states match the `fcst_atm_eval` forecasts and which carry their names (CTRL -> CTRL10-S0, ERA5OP -> CAPT10-S0, DARTOP -> DART40-S0); the mixed runs are labelled `<atm>A-<lnd>L-S0` by the experiments their states come from (e.g. ERA5DL -> CAPTA-DARTL-S0); `EXPERIMENTS = EXPERIMENTS_FCST` by default. Files and tables keep the run keys; `util.hind_factorial.relabel_figure` applies the labels when a figure is saved. Per-source files (sections, indices, projections) are shared between the sets; the tables and figures of a subset go to a `fcst_aligned/` subdirectory of the data and figure directories. The factorial effects and scorecards, and the noise floor of the SSW RMSE attribution, need the full set and are skipped otherwise (the SSW attribution then shows only the bars whose two runs are present).

Every `hind_atm_*` notebook covers one metric in the same order: read input, derive, save (then reload every source from disk), scores or checks, plot. Saved diagnostics are reused: input fields are opened lazily, and only for sources whose output file is missing, unless `FORCE_COMPUTE = True` (a rebuilt NAO or RMM EOF file recomputes every index projected on it). The long training records behind the NAO and RMM EOF checks are read only when the EOFs are built or `RECHECK = True`. In the derive step the expensive reductions of all sources are computed in one Dask graph (`util.dask_helpers.compute_all`); `DASK_WORKERS` (default 8) caps the threads, which bounds memory on many-core login nodes; `USE_DISTRIBUTED = True` runs the graph on a local `dask.distributed` cluster of `DASK_WORKERS` processes instead, which read NetCDF files in parallel.

These notebooks replace the former `v3_dart_paper/` notebooks (`mjo_analysis`, `qbo_analysis`, `ssw_analysis`, `nao_analysis`, `tropicwind_analysis`), whose duplicated classes and `ExperimentUtil.py` copies are now the modules above.

## Running

Notebooks load repository modules from `WORK_DIR/diagnostic`, so they can be
opened from any directory; set `WORK_DIR` in the parameter cell to use another
checkout.
