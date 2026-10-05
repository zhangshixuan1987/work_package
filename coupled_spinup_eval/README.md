# Coupled spin-up evaluation

Notebook-driven analyses of the E3SM v3 coupled spin-up: the fully coupled spin-up (Full-CPL / FC) with its
piControl, and the alternating fully coupled / forced ocean–sea ice spin-ups (FC-FOSI-S1, S2).

* `jupyter/` — one notebook per analysis. **All setup lives in the notebook**: paths, the experiment table,
  labels/colors, years, and each step's figure parameters and plot style. Change them there (or copy the
  notebook for a new scenario).
* `util/` — only reusable code: one shared module per kind of analysis. Notebooks import the classes and drive
  them; nothing in `util/` holds experiment- or scenario-specific settings.
* `config/pltvar_region_sources.json` — e3sm_diags variable/region catalog (notebooks 8–9).
* `scripts/submit_slurm.sh` — execute notebooks headlessly on a compute node.

## Notebooks

| notebook | shows | classes used (`util/`) |
|---|---|---|
| `1_enso_compare_analysis` | Niño3.4 SST, index and power spectrum, FC vs FC-FOSI | `NinoSSTPlotter`, `NinoIndexPlotter`, `OCNEnsoSpectrum` (regional_timeseries) |
| `2_ocn_sfc_compare_analysis` | global surface ocean time series | `FigureDataCollector`, `OCNDiagnosticsPlotter` (alt_global_timeseries) |
| `3_tke_amoc_compare_analysis` | surface TKE and AMOC time series | same as 2 |
| `4_ohc_compare_analysis` | global OHC by layer and depth–time Hovmöller | `OHCFigureDataCollector`, `OHCDiagnosticsPlotter` (alt_global_timeseries) |
| `5_ohu_compare_analysis` | ocean heat uptake vs surface fluxes | `OHUFigureDataCollector`, `OHUDiagnosticsPlotter` (alt_global_timeseries) |
| `6_computation_cost_analysis` | cost of alternating FC / FOSI segments | `SimulationRun`, `AlternatingSimulationEvaluator` (computation_cost) |
| `7_ocn_map_compare_analysis` | SST, SSS and AMOC climatologies (+ combined figure), global MOC | `Surface2DClimoPlotter`, `AMOCClimoPlotter` (climo_maps) |
| `8_model_error_drift_analysis` | drift of e3sm_diags PCC/RMSE, Full-CPL | `SpinupMetricAnalyzer` (model_error) |
| `9_model_error_compare_analysis` | e3sm_diags PCC/RMSE, FC vs FC-FOSI (global and regional) | `SpinupMetricAnalyzer` (model_error) |
| `10_land_ctrl_analysis` | land BGC spin-up stability, Full-CPL + piControl | `LandCtrlStability` (land_bgc) |
| `11_land_compare_analysis` | land BGC, FC vs FC-FOSI | `LandCompareStability` (land_bgc) |
| `12_atm_flux_ctrl_analysis` | atmospheric energy and water balance trends | `AtmosphereBalanceAnalysis` (ctrl_processing), `AtmosphereBalancePlotter` (ctrl_plotting) |
| `13_ocn_conserve_ctrl_analysis` | ocean volume, salt and freshwater conservation | `OceanBudgetBuilder`, `CoupledBudgetPlotter` — needs the output of 12 |
| `14_ocn_ts_ctrl_analysis` | SST, SSS, SSH, AMOC and sea-ice volume drift | `SSTBuilder`, `SSSBuilder`, `SSHBuilder`, `AMOCDiagnosticsBuilder`, `SeaIceBuilder`, `OceanDriftPlotter` |
| `15_ohc_ts_ctrl_analysis` | OHC drift by layer | `OHCBuilder`, `OHCDriftPlotter` |
| `16_recoupling_atm_analysis` | recoupling shock, annual global means 201–350 | `load_or_read`, `plot_annual_timeseries` (area_mean_drift) |
| `17_ohc_corr_compare_analysis` | OHC Hovmöller and lagged flux–OHC correlation (0001–0220, 0221–0350) | `OHCHovmollerPlotter`, `LaggedOHCCorrelation`, `write_annual_flux` (regional_timeseries) — needs the output of 5 |
| `18_ohc_corr_full_period_analysis` | the same lag correlation over 0001–0350 (backup supplement to 17) | `LaggedOHCCorrelation` (regional_timeseries) — needs the output of 17 |
| `19_ohc_2d_compare_analysis` | OHC in 2-D, FC-FOSI − FC: 1 horizontal ΔOHC maps per layer (MPAS-Analysis `deltaOHC`, 0301–0350, piControl significance) + combined figure; 2 zonally integrated OHC, latitude × depth (decadal `ocn_2d` climatologies, years ≤ 250, no significance) | `Surface2DClimoPlotter` (climo_maps), `ZonalOHCSection` (ohc_zonal) |
| `20_cmip6_rmse_compare_analysis` | e3sm_diags seasonal RMSE of 9 atmospheric fields against the CMIP6 historical spread (box plots per window; ANN RMSE through the spin-up) | `read_cmip6_csv`, `read_e3sm_diags_rmse`, `plot_rmse_boxes` (cmip6_rmse) |
| `21_atm_drift_obs_analysis` | drift of the global-mean surface climate (TS, TREFHT, PSL, PRECT, SHFLX, LHFLX, FLUT) against NOAA-20C 1871–1890: sliding 20-yr bias and RMSE | `load_or_read` (area_mean_drift), `obs_annual_global_mean`, `sliding_drift`, `plot_drift` (obs_drift) |
| `22_fw_flux_compare_analysis` | ocean surface freshwater flux by source (rain, evaporation, snow, runoff, sea ice, frazil, net), FC vs FC-FOSI; FOSI-segment shading and period-mean differences | `read_fw_annual`, `plot_fw_components`, `plot_fw_differences` (fw_flux) |
| `23_natl_density_flux_analysis` | North Atlantic deep-water formation: surface density flux (thermal / haline), March mixed-layer depth and convective area in the Labrador, Irminger, GIN Seas and subpolar North Atlantic, per decade, with AMOC at 26°N | `regional_decades`, `plot_region_panels` (natl_density_flux) |

Every notebook has the same structure:

1. **Setup** — `MODEL_ROOT` (case directories), `OUTPUT_ROOT`/`OUTPUT_URL` (public_html), `FIG_DIR` (figures),
   `DATA_DIR` (cached/derived data), `EXPERIMENTS` (case names and periods), `LABEL`,
   and the analysis settings the notebook uses: `EXPS`, `COLORS`, `YEARS`, `CLIMO_LEN`,
   `SWITCH_YEAR` (alternate spin-up comparisons) or `EXPERIMENT`, `CTRL_YEARS`, `SPINUP_MARKERS`
   (single-experiment spin-up + piControl analyses).
2. **Shared inputs** — `run_dict = make_run_dict(EXPERIMENTS, MODEL_ROOT)`, output folders.
3. Per analysis step: **Settings** (its parameters and plot style) and **Compute & plot** (instantiate the
   class, extract/cache, plot, `publish(FIG_DIR, DATA_DIR)`).

Run from `jupyter/` with the E3SM Unified kernel. Batch:
`sbatch scripts/submit_slurm.sh jupyter/<notebook>.ipynb [...]` (no argument = all, in order);
executed copies go to `jupyter/executed/`.

## Outputs and inputs

* Everything is written to public_html (made world-readable),
  `OUTPUT_ROOT = /lcrc/group/e3sm/public_html/diagnostic_output/ac.szhang/v3spinup_paper`: figures to
  `figures/<topic>/` and cached/derived data to `data/<topic>/` — the same `<topic>` name, so each figure folder
  has a matching data folder (notebooks 10 and 11 share `data/land/`). `data/README.md` indexes every folder.
* Inputs: the raw model output under `MODEL_ROOT`, and in `OUTPUT_ROOT/data/`: `obs/` (HadSST Niño3.4,
  notebook 1), `grid/depth.nc` (MPAS-Ocean level depths), `land_mask/landmask.nc`, `diag_data/` (frozen copy of
  the e3sm_diags ANN and seasonal metrics tables, notebooks 8–9 and 20), `cmip6/` (CMIP6 benchmark RMSE tables,
  notebook 20); NOAA-20C observational time series from the e3sm_diags data on LCRC (notebook 21).
* Re-running a notebook reuses its caches in `data/<topic>/`; delete a file to recompute it. Several ctrl
  processing settings default to `force_reprocess = True`, i.e. they rebuild from raw model output; set them to
  `False` in the settings cell to reuse caches.

## util/

| module | contents |
|---|---|
| `experiments.py` | `make_run_dict`, `build_exp_subdirs_from_run_dict` |
| `figures.py` | `publish` (world-readable outputs for public_html), `compare_pdfs` (pixel difference of two PDFs) |
| `ctrl_processing.py` | `MPASDiagnosticsBuilder` base → `SSTBuilder`, `SSSBuilder`, `SSHBuilder`, `SeaIceBuilder`, `OHCBuilder`, `OceanBudgetBuilder`; `AMOCDiagnosticsBuilder`; `AtmosphereBalanceAnalysis` |
| `ctrl_plotting.py` | `MPASTimeDiagnosticsPlotter` base → `OceanDriftPlotter`, `OHCDriftPlotter`; `AtmosphereBalancePlotter`; `CoupledBudgetPlotter` |
| `alt_global_timeseries.py` | `FigureExtractSpec`, `FigureDataCollector`, `OCNDiagnosticsPlotter`, `PlotRequest` (+ `OHU…`, `OHC…` variants) |
| `regional_timeseries.py` | `OCNTimeSeriesPlotter` base → `NinoSSTPlotter`, `NinoIndexPlotter`, `OHCHovmollerPlotter`; `OCNEnsoSpectrum`; `LaggedOHCCorrelation` |
| `climo_maps.py` | `Surface2DClimoPlotter`, `AMOCClimoPlotter` |
| `ohc_zonal.py` | `ZonalOHCSection` (zonally integrated heat content per metre depth per degree latitude from the decadal `ocn_2d` climatologies), `decadal_climo_files` |
| `cmip6_rmse.py` | `FIELDS`, `read_cmip6_csv` (CMIP6 benchmark table), `read_e3sm_diags_rmse` (e3sm_diags metrics tables), `plot_rmse_boxes` |
| `fw_flux.py` | `COMPONENTS`, `ocean_area`, `read_fw_annual` (MPAS global freshwater fluxes, mm yr⁻¹; net = sum of sources), `plot_fw_components`, `plot_fw_differences` |
| `natl_density_flux.py` | `REGIONS`, `regional_decades` (decadal `ocn_2d` climatologies → regional density flux, SST/SSS, March MLD), `plot_region_panels` |
| `obs_drift.py` | `obs_annual_global_mean`, `model_annual_global_mean`, `sliding_drift`, `plot_drift` |
| `computation_cost.py` | `SimulationRun`, `AlternatingSimulationEvaluator` |
| `model_error.py` | `SpinupMetricAnalyzer` |
| `area_mean_drift.py` | `read_area_means`/`load_or_read` (post/time_series/atm area means), `AreaMeanDrift` (bias vs piControl, σ-scaled drift, difference from a baseline run, annual-cycle corr/RMSE), `plot_window_heatmaps`, `plot_annual_timeseries` |
| `land_bgc.py` | `SpinupStabilityBGC` base → `LandCtrlStability`, `LandCompareStability` |
| `ocn_timeseries.py`, `ohc_timeseries.py`, `flux_timeseries.py` | MPAS-O time-series extractors (flux v1 + `…V2` used by the ocean-heat-uptake analysis) |

Each class family is one **base class** holding the common code plus a thin **subclass per analysis** that
overrides only what differs. Smaller differences are options:

* `SpinupMetricAnalyzer(bias_ratio_metric=...)` — `"Mean_Bias"` (notebook 9, default) or `"Bias"` (notebook 8).
* `Surface2DClimoPlotter.plot_map_compare(ref_cmap_n, ref_darken, diff_cmap_n)` — SST defaults 200/1.0/200; SSS uses 256/0.8/256.
* `AMOCClimoPlotter.plot_amoc_2d_compare` (contours) and `plot_amoc_2d_compare_mesh` (pcolormesh).

## Notes

* Notebook dependencies: 13 needs `data/atm_flux/` from 12; 17 needs `data/ohu/` from 5; 18 needs
  `data/ohc_corr/` from 17. A batch run with no arguments runs the notebooks in numeric order.
* The MPAS `netMassFlux` and `netEnergyFlux` series are 0 in the first month after every restart (accumulator
  reset; every 5th year in the FC-FOSI runs). The net fluxes are therefore built from their components:
  `Fmass` = sum of the freshwater sources (notebook 2 panel (c), notebook 22) and `Qnet` = heat fluxes +
  ρ₀c_p × temperature fluxes, which reproduces `netEnergyFlux` exactly in all other months (notebooks 5, 17, 18).
  The raw MPAS fields remain available as `Mass_net` and `Qnet_mpas`.
* Q̄source in notebooks 5, 17 and 18 is Qsw + Qlw_dn. Lag-correlation files carry the flux name:
  `lagcorr_<flux>_ohc_bydepth_<years>.nc`.
* `sh: getfattr: command not found` messages come from the netCDF-C library and are harmless.
