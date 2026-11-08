# RTENSO2026 case-study diagnostics

Diagnostics for the RTENSO2026 hindcasts (`../hindcast_run/RTENSO2026/`,
post-processed by `../postprocess/RTENSO2026/`), treated as case studies of the
2015/16 El Nino (2015 Feb/May/Aug/Nov starts) and the 2026/27 El Nino (2026
Feb/May/Aug starts; no Reanalysis-FOSIRL runs). Each notebook's `EVENTS` lists the events. The notebooks follow the layout and conventions of the ESP-Lab
suite (`~/code/ESP-Lab/jupyter`): a run-control block at the top, cache-first
products under `S2D_DIAG_ROOT`, and figures named with `figure_filename`. ESP-Lab
measures skill over many start years; here one event is followed through 4 starts.

| Experiment | Starts (2015) | Members | Differs by |
| --- | --- | --- | --- |
| `Reanalysis-CLIM` (`RTENSO2026_clim`) | Feb, May, Aug, Nov | 10 | land/sea-ice/river ICs from the BruteForce year closest to the same-date climatology |
| `Reanalysis-OBS` (`RTENSO2026_obs`) | Feb, May, Aug, Nov | 10 | land/sea-ice/river ICs from the BruteForce year closest to the ERA5 anomaly |
| `Reanalysis-NSIDC` (`RTENSO2026_obs_ice_nsidc`) | placeholder | 10 | `Reanalysis-OBS` with NSIDCv6 sea ice; enabled by setting its `init_months` once post-processed |
| `Reanalysis-FOSIRL` (BruteForce) | Feb, May, Aug, Nov | 10 | the actual 2015 restarts |

All share the ERA5 atmosphere and CFSR ocean ICs. Cases whose `init_months` is empty
are placeholders and are dropped from the notebooks.

## 1 Feb 2026 start: `Reanalysis-CLIM` and `Reanalysis-OBS` are the same runs

For this one start the `clim` and `obs` IC selections both chose the BruteForce **2010**
land, river and sea-ice state, so the two hindcasts started from identical initial
conditions and are bit-identical member for member (Nino3.4 difference 0.0 everywhere).
Every other start uses different years for the two methods; see "Selected years" in
`../ensemble_init/RTENSO2026/README.md`. The runs were kept as they are.

What this means in the diagnostics:

* **Feb 2026 has no CLIM vs OBS contrast.** In plots the `Reanalysis-CLIM` line or map is
  hidden under `Reanalysis-OBS`; any CLIM-OBS difference for that start is zero by
  construction, not a result.
* **Pooled ensembles count those members once.** The notebook 4 inter-member regression
  and the notebook 6 outlook pool `Reanalysis-CLIM`, `-OBS` and `-NSIDC`, and skip a case
  whose members repeat an already pooled case exactly (each run prints a note). The Feb
  2026 pool is therefore 20 members (CLIM + NSIDC) instead of 30.
* `Reanalysis-NSIDC` is unaffected: its NSIDC sea ice differs from the 2010 state.

## Drift correction

Model anomalies remove the Reanalysis (BruteForce) 1981-2010 lead-time climatology,
the same IC system:

- May and Nov starts use the Reanalysis climatology for that start month and lead.
- Feb and Aug starts have no Reanalysis climatology (only 1980 and 1985 exist). The
  default `CLIM_METHOD = "lead_drift"` uses the observed 1981-2010 climatology of each
  target month plus the Reanalysis drift at the same lead, averaged over the May and
  Nov starts. `"target_month"` instead borrows the May/Nov climatology for the same
  target month at the nearest lead; it is biased at short leads because the
  equatorial cold tongue drifts about 1 C in the first months (Nino3.4 error at the
  Aug start, lead 0: +0.82 C with `target_month`, -0.10 C with `lead_drift`).

Observed anomalies are relative to the 1981-2010 monthly climatology: HadISST2 SST
(land-masked, verified 60S-60N), ERA5 `TREFHT`, `PSL`, `U850`, and GPCP v2.3 `PRECT`,
all conservatively regridded to the model 1x1 grid.

These e3sm_diags time series end in 2022 (HadISST2), 2019 (ERA5) and 2017 (GPCP). For
later years, such as the 2026 starts, `obs_anomaly(field, years=(2026, 2027))` switches to
the 1x1 archive products under `/global/cfs/cdirs/e3sm/zhan391/data` (`OBS_ARCHIVE`), which
then also supply the 1981-2010 climatology: NOAA OISST v2.1 for `TS` (to Oct 2026; Jan-Aug
climatology 1982-2010), GPCP v2.3 for `PRECT` (to Jun 2026) and NOAA/NSIDC CDR v6 for `SICE`
(to Sep 2026), and NOAA CDR for `FLUT` (OLR, to Sep 2026). There is no archive product yet for
`PSL`, `U850`, `Z500`, `T925` or the 200-hPa winds after 2019; `TREFHT` over land comes from
ERA5-Land (`product="ERA5_Land"`, to Sep 2026).

The observation caches record the latest modification time of the archive files, so
`derivation_mode = "auto"` rebuilds them after a product is extended (no manual cleanup needed).

Alternative products are picked by name with `obs_anomaly(field, product=...)` from
`OBS_ALT_PRODUCTS`: for `TS`, `HadISST1.1` (`data/HadSST`, 1870 on), `OISSTv2.1`, `ERSSTv5`
(1854 on) and `ESA_CCI_SST` (ESA CCI SST v3.0, 1980 on), each with its own 1981-2010
climatology. `1_enso_plume.ipynb` plots HadISST 1.1, ERSSTv5, ESA CCI SST and NOAA OISST v2.1 and keeps
HadISST 1.1 as the reference for errors. ERA5-Land has no ocean
points, so it gives no Nino index.

## Files

| File | Purpose |
| --- | --- |
| `util/rtenso_case.py` | Shared helpers, imported by every notebook as `from util import rtenso_case as rc`: case registry, hindcast loading, lead-time climatologies, observation products and anomalies, cache builds (`build_caches`), pattern ACC/RMSE, regions, velocity potential, output paths. |
| `0_build_caches.ipynb` | Builds the Reanalysis lead-time climatologies and observation caches ahead of the other notebooks (`FIELDS`, `derivation_mode`), through `rc.build_caches`. |
| `1_enso_plume.ipynb` | Nino3.4/3/4/1+2/ONI forecast plumes per start against HadISST 1.1, ERSSTv5, ESA CCI SST and NOAA OISST v2.1 with the observed peak period shaded, peak-season summary, ensemble-mean error and spread (after 3a), and a map of the Nino regions. |
| `2_anomaly_maps_acc.ipynb` | Peak-season (DJF) anomaly maps from every start, and spatial pattern ACC and RMSE by lead season and region, per member and ensemble mean (after 1a/1b); only fully observed seasons are scored. Set `FIELD` to `TS`, `PRECT` (both events), `TREFHT`, `PSL` or `U850` (2015 only). |
| `3_ic_sensitivity.ipynb` | `Reanalysis-CLIM` vs `Reanalysis-OBS` vs `Reanalysis-FOSIRL` differences: land-mean differences with member-based intervals, fraction of land with significant differences by lead, stippled difference maps, and hemispheric sea-ice area against HadISST2. |
| `4_teleconnections.ipynb` | Peak-season teleconnections of `PRECT`/`TREFHT`/`PSL` per degC of Nino3.4 against the observed event and the 1982-2010 observed regression, an inter-member regression, and CPC-style equatorial Hovmollers of SST, U850, OLR and 200-hPa velocity potential (`rc.velocity_potential`) (after 3b). |
| `5_npac_namerica.ipynb` | CPC-style North Pacific / North America maps (120E-60W, 20-70N) for Nov-Feb: 500-hPa height and anomaly, 200-hPa wind, 925-hPa temperature anomaly, observed (ERA5, through 2019) and each case's ensemble mean. |
| `6_us_temp_precip.ipynb` | CPC-style U.S. maps: 1- and 3-month percent of average precipitation and temperature departures, and a tercile-probability seasonal outlook (NDJ) from the pooled members with the observed category. |
| `util/xro_case.py` | XRO state vector (Nino3.4, recharge *h*, NPMM, SPMM, IOB, IOD, SIOD, TNA, ATL3, SASD boxes) and the configuration shared by the 7a-7d notebooks (`EVENTS`, `CASES`, `XRO_SETUP`); cached observed and hindcast state vectors, training-period trend, XRO fit, event forecasts aligned with the hindcast starts, initial-state sensitivity experiments and out-of-sample rolling hindcasts. |
| `util/pubfig.py` | Publication figure style for 7a-7d: figures drawn at journal column width (3.5/7.2 in) with 7-9 pt text, saved as vector PDF (embedded fonts) and 600-dpi PNG. |
| `7a_xro_validation.ipynb` | XRO (Zhao et al. 2024) validated against observations only: training-period trends (linear vs quadratic *h*), fitted recharge-oscillator growth rate and period, stochastic simulations against observed seasonal synchronization, *h*-*T* lag correlation and persistence, out-of-sample skill (refit every year 1995-2025) by lead and initial month against persistence, DJF amplitude (does XRO underpredict strong El Ninos?), and calibration of the stochastic ensemble. |
| `7b_xro_events.ipynb` | XRO forecasts of the 2015/16 and 2026/27 events from the months before each start: Nino3.4 and *h* plumes, peak-season summary, peak-season Nino3.4 sensitivity to each initial-state component (removal experiments, with the non-additive remainder shown), and to the initial *h* alone. |
| `7c_e3sm_xro_state.ipynb` | E3SM hindcasts in the XRO state space against observations and XRO: state-vector plumes (common 5-95% ranges), error decomposition into the start-month error B(0) and the error accumulated afterwards G(L) = B(L) - B(0), summary heat maps, growth of Nino3.4 and *h* from the start month, peak-season summary. |
| `7d_e3sm_recharge.ipynb` | *T*-*h* phase planes, and XRO started from E3SM's start-month offsets (*h* only, and the full state) to test whether the E3SM-XRO peak difference comes from the initial state or from E3SM's evolution. |
| `8_enso2026_animation.ipynb` | Animated month-by-month evolution of the 2026/27 El Nino: global maps (Robinson; `DOMAIN` for a regional box) of SST (OISST v2.1), 0-700 m heat content (EN4.2.2), precipitation (GPCP v2.3) and land 2-m air temperature (ERA5-Land) anomalies, observed and the pooled (`clim`/`obs`/`NSIDC`) ensemble mean of each 2026 start, above a Nino3.4 plume drawn up to the frame month, with the deterministic XRO forecast from each start (7b control setup) as a dashed reference. Saves a GIF and a video per field (H.264 MP4 if ffmpeg has libx264, e.g. via `imageio-ffmpeg`; else VP9 WebM) and a static storyboard of selected months; the pooled ensemble means are cached under `data/`. |
| `util/strat_case.py` | Stratosphere readers: daily zonal means and daily 10-hPa maps from the native `eam.h1` files, monthly zonal means of the post-processed U/T/Z3, ERA5 daily zonal means and 00Z maps, lead-time climatologies of zonal means (`direct`, `lead_drift`), the daily-minus-monthly Z offset, SSW central dates (Charlton and Polvani 2007), final warming with coverage status, event details, lag composites. |
| `util/strat_diag.py` | Shared inputs of the 9a-9i notebooks: `StratConfig` (data roots and run controls), `load_era5` (cached ERA5 indices and climatologies), `load_references` (Reanalysis zonal means of every start month), `load_model` (pooled 2026 starts with their daily and monthly indices), `ssw_catalogue`, `djf_vortex`, DJF helpers. |
| `util/vortex_plot.py` | North-polar map helpers of 9g-9i: projection and extent, fixed-level anomaly shading and total contours, U60 badges, video writer. |
| `9a_strat_build_caches.ipynb` | Builds and checks every stratosphere cache (the only 9x notebook that reads raw data; the others use `derivation_mode = "require"`): ERA5 daily zonal means and indices, model daily/monthly zonal means and daily 10-hPa maps, Reanalysis climatologies of every start month (`direct`: Feb/Aug from the six 1985-2010 starts); cache inventory; the daily-minus-monthly height offset check (38 m at 10 hPa, constant over leads and starts). |
| `9b_strat_model_climatology.ipynb` | E3SM's stratosphere against ERA5: seasonal cycle of U(60N, 10 hPa) and polar-cap temperature, DJF zonal-mean U bias, DJF vortex strength and interannual variability. |
| `9c_vortex_forecast.ipynb` | The vortex forecast: Feb 2026 start against the observed end of winter 2025/26, daily U(60N, 10 hPa) plumes for 2026/27 and 2027/28, the vortex index plume (U60, polar-cap T and height, polar-cap SLP), DJF tercile probabilities with their sensitivity to the climatology method. |
| `9d_ssw_occurrence.ipynb` | SSW frequency by start x case x winter (raw wind; bias-adjusted as sensitivity) against ERA5 1979/80 on, central-date months, final warming with censoring; event catalogue CSV. |
| `9e_enso_qbo_background.ipynb` | QBO forecasts (time-height and 50 hPa), the ENSO/QBO background state by start and case, within-group member covariation with the DJF vortex. |
| `9f_ssw_dynamics.ipynb` | Season-matched SSW lag composites (U60, 100-hPa heat flux, polar-cap height) with cluster-bootstrap intervals, event-property distributions, vertical evolution of polar-cap anomalies and SLP. |
| `9g_ssw_mechanism.ipynb` | What distinguishes an SSW member from comparable members that do not reverse? A featured SSW member (`EVENT`) against a matched non-SSW member and a near-miss: aligned U60 / polar-cap T / heat-flux timelines, ensemble ranks, the ensemble relationship of preconditioning and pre-onset wave activity with the outcome, pre-onset 10-hPa maps, side-by-side animation. |
| `9h_vortex_ssw_animation.ipynb` | Day-by-day 10-hPa polar maps of one winter: each start's ensemble mean and a featured SSW member, the all-start mean (or ERA5 for an observed winter), with a U60 timeline; GIF, video, storyboard. |
| `9i_vortex_evolution_animation.ipynb` | Seasonal evolution Oct 2026-Apr 2027 without event-based selection: one fixed member and the ensemble mean of each start, raw 10-hPa height contours with optional anomaly shading, reference climatologies, U60 timeline; GIF, video, preview, storyboard. |

## Running

Use the `e3sm_analysis` kernel (ESP-Lab installed in editable mode) and start Jupyter in
this directory so the `util` package (`util/rtenso_case.py`) imports. Build the caches once with
`0_build_caches.ipynb` (about 4 minutes per field and start month; the U850 reads are the slowest),
interactively or headless:

```bash
conda activate e3sm_analysis
jupyter nbconvert --to notebook --execute --inplace 0_build_caches.ipynb
```

Each notebook's run-control block holds `derivation_mode` ("auto" reuses valid caches,
"rebuild" recomputes, "require" fails on a missing cache) and `CLIM_METHOD`.

Inputs: the RTENSO2026 cases read their time series from
`/global/cfs/cdirs/e3sm/S2S2D/RTENSO2026/post_process` (`rc.RTENSO2026_POST_ROOT`);
the Reanalysis (BruteForce) cases read from `RAW_MODEL_ROOT`
(`/global/cfs/cdirs/e3sm/S2S2D/post_process`).

Outputs:

- Notebook products: one folder per notebook under the fixed
  `/global/cfs/cdirs/e3sm/www/zhan391/e3sm-rtenso2026/` (`rc.OUTPUT_ROOT`), kept apart from
  the ESP-Lab gallery. Each notebook sets `NOTEBOOK` and writes through `rc.save_figure` and
  `rc.output_dir`:

  ```text
  e3sm-rtenso2026/
    1_enso_plume/{figures,data}/          data: Nino index netCDF per case/start, peak summary CSV
    2_anomaly_maps_acc/{figures,data}/    data/<FIELD>/: pattern-skill netCDF, ACC CSV
    3_ic_sensitivity/{figures,data}/
    4_teleconnections/{figures,data}/     data: teleconnection ACC CSV
    5_npac_namerica/figures/
    6_us_temp_precip/figures/
    7a_xro_validation/{figures,data}/     data: RO parameters, skill, DJF amplitude and calibration CSVs
    7b_xro_events/{figures,data}/         data: XRO event forecasts (netCDF), peak and sensitivity CSVs
    7c_e3sm_xro_state/{figures,data}/     data: start-month and evolution error CSVs, peak summary
    7d_e3sm_recharge/{figures,data}/      data: matched-initial-state CSV
    8_enso2026_animation/{figures,data}/  figures: GIF/video animations and storyboards; data: pooled ensemble-mean caches
    9a_strat_build_caches/ ... 9i_vortex_evolution_animation/   one folder per 9x notebook: figures, CSV tables
  ```
- Reusable caches stay outside it: `$S2D_DIAG_ROOT/Reanalysis/leadtime_clim/`,
  `$S2D_DIAG_ROOT/<source>/strat/{daily_zonal,monthly_zonal,daily_maps}/` and `observations/strat/daily_maps/` (notebooks 9a-9i) and
  `$S2D_DIAG_ROOT/observations/case_study/` (`S2D_DIAG_ROOT` defaults to
  `/global/cfs/cdirs/e3sm/S2S2D/s2d_diag`).

## XRO recharge proxy (`OHC700`)

The hindcasts save no subsurface temperature, so the 7a-7d notebooks use the 5S-5N, 120E-80W
mean 0-700 m heat content anomaly as the recharge *h* (instead of the 20 C isotherm depth).
MPAS-O writes `oceanHeatContentSfcTo700m` per native cell (J, integrated with T in K);
`rc.open_member` divides it by the native cell area remapped with the file's own `map_file`
weights (`rc.mpas_cell_area`), giving GJ m-2. Observations are EN4.2.2 OHC700
(`data/ENv4.2.2/1x1`, 1979 to Jul 2026). Anomalies remove the 1981-2010 climatology as for
the other fields; the constant offset between K- and C-based integrals cancels at fixed depth
(if MPAS integrates from the moving sea surface, sea-level anomalies add about 0.1 GJ m-2 per
0.1 m in phase with recharge; not yet checked). The XRO state vectors and rolling hindcasts are
cached under `$S2D_DIAG_ROOT/observations/xro/` and `$S2D_DIAG_ROOT/<source>/xro_state/`; rerun
with `derivation_mode = "rebuild"` after new observations are added. The XRO
package (`pip install XRO`, examples in `~/code/XRO`) must be installed in `e3sm_analysis`.

Land and sea-ice diagnostics in `3_ic_sensitivity.ipynb` need the RTENSO2026 `lnd_mon1`
and `ice_mon1` post-processing; cases not yet processed are skipped.

## Stratosphere series (9a-9i)

Run `9a_strat_build_caches` first (about 1 h the first time: the native daily model reads); the other 9x
notebooks then only read caches (`derivation_mode = "require"`) and can run in any order. Each starts with
its question and findings. The polar-vortex notebooks need `ERA5/1x1/daily_zonal` in the reference
archive (built by `ERA5/scripts/daily_zonal`). On a login node (30 GB per user) run one 9x notebook at a time;
`.split_work/run_strat_series.sbatch` runs the whole series on a compute node.

The notebooks are generated by `.make_strat_series.py` from the cell sources in `.split_work/` (the former
notebooks 9, 10 and 11, split by question) and the shared loaders in `util/strat_diag.py`; edit those and
regenerate (`python .make_strat_series.py`) rather than the notebooks, or edit a notebook directly and
carry the change back.

Model anomalies use each start's own Reanalysis climatology (`CLIM_METHOD = "direct"`): May/Nov
1981-2010, Feb/Aug the six 1985-2010 starts. The earlier `lead_drift` climatology produced a spurious
weak-vortex signal for the Feb/Aug starts (9c keeps it as a sensitivity test).
