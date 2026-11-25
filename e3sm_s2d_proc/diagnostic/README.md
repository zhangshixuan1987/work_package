# RTENSO2026 case-study diagnostics

Diagnostics for the RTENSO2026 hindcasts (`../hindcast_run/RTENSO2026/`,
post-processed by `../postprocess/RTENSO2026/`), treated as case studies of the
2015/16 El Nino (2015 Feb/May/Aug/Nov starts) and the 2026/27 El Nino (2026
Feb/May/Aug/Nov starts; no Reanalysis-FOSIRL runs). Each notebook's `EVENTS` lists the events. The notebooks follow the layout and conventions of the ESP-Lab
suite (`~/code/ESP-Lab/jupyter`): a run-control block at the top, cache-first
products under `S2D_DIAG_ROOT`, and figures named with `figure_filename`. ESP-Lab
measures skill over many start years; here one event is followed through 4 starts.

| Experiment | 2015 starts | 2026 starts | Members | Differs by |
| --- | --- | --- | --- | --- |
| `Reanalysis-CLIM` (`RTENSO2026_clim`) | Feb, May, Aug, Nov | Feb, May, Aug, Nov | 10 | land/sea-ice/river ICs from the BruteForce year closest to the same-date climatology |
| `Reanalysis-OBS` (`RTENSO2026_obs`) | Feb, May, Aug, Nov | Feb, May, Aug, Nov | 10 | land/sea-ice/river ICs from the BruteForce year closest to the ERA5 anomaly |
| `Reanalysis-NSIDC` (`RTENSO2026_obs_ice_nsidc`) | Feb, May, Aug, Nov | Feb, May, Aug, Nov | 10 | `Reanalysis-OBS` with NSIDCv6 sea ice |
| `Reanalysis-FOSIRL` (BruteForce) | Feb, May, Aug, Nov | none | 10 | the actual 2015 restarts (the Reanalysis set ends in 2015) |

All share the ERA5 atmosphere and CFSR ocean ICs. Every set is complete: 24 months for all members,
atmosphere, ocean, land and sea ice, except the 1 Nov 2026 start of the three RTENSO2026 cases, which was
extended to 36 months (Nov 2026 - Oct 2029; `rc.EXTENDED_STARTS`, series in `ts/monthly/3yr`). The
notebooks read each start at its own length (`rc.start_nlead`).

The events, cases and their starts are defined once in `util/rtenso_case.py`: `rc.EVENTS` (label, starts,
DJF peak season), `rc.CASES` (case prefix, cache `source`, `years`, `init_months`, colour, short label),
`rc.POOLED_CASES` (CLIM, OBS and NSIDC, pooled in `enso/03_teleconnections`, `enso/05_us_temp_precip`, `enso/10_enso2026_animation` and the stratosphere notebooks) and `rc.has_start(case, year,
init)`. The notebooks, `util/xro_case.py` (which adds the XRO training years) and `util/strat_diag.py` read
them from there; change a case or start there, not in a notebook.

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
* **Pooled ensembles count those members once.** The `enso/01_enso_plume` RONI outlook, the `enso/05_us_temp_precip`
  outlook, `enso/10_enso2026_animation` and the stratosphere notebooks pool `Reanalysis-CLIM`, `-OBS` and `-NSIDC`, and skip a case
  whose members repeat an already pooled case exactly (each run prints a note). The Feb
  2026 pool is therefore 20 members (CLIM + NSIDC) instead of 30.
* `Reanalysis-NSIDC` is unaffected: its NSIDC sea ice differs from the 2010 state.

## Drift correction

Every notebook uses the same lead-dependent mean-bias correction (ICPO 2011),
`CLIM_METHOD = "start_drift"`, built from the Reanalysis (BruteForce) hindcasts, the same IC
system. For start month *m* and lead *L*:

    drift_m(L) = mean over the Reanalysis starts of month m in 1981-2010 of
                 [ensemble mean(L) - observed value at the same verification month]
    clim_m(L)  = observed 1981-2010 climatology of the target month + drift_m(L)
    anomaly    = forecast - clim_m(L)

- May and Nov: 30 start years (1981-2010). The observed terms then average back to the
  observed climatology, so the climatology equals the Reanalysis ensemble-mean climatology
  (to within a few hundredths of a degree for Nino3.4; leads that verify after 2010 differ
  slightly).
- Feb and Aug: the Reanalysis has these starts every fifth year, 1980-2015; the six inside
  1981-2010 are 1985, 1990, 1995, 2000, 2005 and 2010 (10 members each). Measuring the drift
  against observations of the same years removes the ENSO state of those years (several were
  La Nina years, so a plain six-year mean of the hindcasts is about 0.5 C too cold in Nino3.4).
- Fields without an observational product (land `SOILWATER_10CM`, `H2OSNO`, `TWS`, `TSA`), and
  points without observations (e.g. `TS` over land), use the Reanalysis ensemble mean of the same
  starts (`"direct"`).
- The observed reference is each field's default product (`rc.FIELDS[field]["obs"]`: HadISST 1.1 for
  `TS`, ERA5, GPCP v2.3, EN4.2.2); the stratosphere uses the ERA5 daily zonal means.
- Uncertainty: `rc.start_drift(field, init)` returns the standard error of the drift over the
  start years (members of one start share their initial state, so the sample size is the number
  of years). For Nino3.4 (against HadISST 1.1) it is 0.13-0.38 C at leads 0-11 of the Feb/Aug
  starts (0.27-0.38 C at the DJF-peak leads) and 0.06-0.18 C for May/Nov. Subsampling the May/Nov
  starts to the same six years changes their Nino3.4 drift by up to 0.36 C (May) and 0.52 C (Nov)
  at leads 0-11, and 0.86 C at Nov leads 20-23, consistent with that standard error.

Alternatives kept for sensitivity tests (`matched_climatology`, `strat_case.lead_climatology`):
`"direct"` (the Reanalysis ensemble mean of the same starts), `"lead_drift"` (observed target-month
climatology plus the mean May/Nov drift; misses the Aug-start drift of about -1.4 C in Nino3.4 at
leads 1-2 and gives a spurious weak polar vortex for Feb/Aug), and `"target_month"` (the May/Nov
climatology of the same target month at the nearest lead).

**Leads past 24 months.** The Reanalysis starts run 24 months, so the third year of the 36-month
1 Nov 2026 start (leads 24-35) has no drift estimate of its own: it reuses the climatology of the lead 12
months earlier, which has the same target month (`rc.extend_leads`; the XRO spread calibration and the
stratosphere tercile bounds follow the same rule). This assumes the drift no longer changes after the
second year; the step it leaves at lead 24 (e.g. about +0.4 C in the Nov 2026 Nino3.4 ensemble mean) is
part of that assumption, and the figures mark where it starts.

Observed anomalies are relative to the 1981-2010 monthly climatology: HadISST 1.1 SST
(`data/HadSST`, 1870 to Sep 2026, land-masked where the model land fraction is 0.5 or more; the
SST product everywhere: drift reference, observed climatology, Nino indices, maps), ERA5 `TREFHT`,
`PSL`, `U850`, and GPCP v2.3 `PRECT`, all conservatively regridded to the model 1x1 grid.

The e3sm_diags time series end in 2019 (ERA5) and 2017 (GPCP). For later years, such as the
2026 starts, `obs_anomaly(field, years=(2026, 2027))` switches to the 1x1 archive products under
`/global/cfs/cdirs/e3sm/zhan391/data` (`OBS_ARCHIVE`), which then also supply the 1981-2010
climatology: GPCP v2.3 for `PRECT` (to Oct 2026) and NOAA/NSIDC CDR v6 for `SICE`
(to Oct 2026; also `product="NSIDCv6"` for any year, used by `initialization/01_ic_sensitivity` for both events), and NOAA CDR for `FLUT` (OLR, to Oct 2026). There is no archive product yet for
`PSL`, `U850`, `Z500`, `T925` or the 200-hPa winds after 2019; `TREFHT` over land comes from
ERA5-Land (`product="ERA5_Land"`, to Sep 2026).

The observation caches record the latest modification time of the archive files, so
every notebook rebuilds them after a product is extended (no manual cleanup needed).

Alternative products are picked by name with `obs_anomaly(field, product=...)` from
`OBS_ALT_PRODUCTS`: for `TS`, `HadISST1.1` (`data/HadSST`, 1870 on), `OISSTv2.1`, `ERSSTv5`
(1854 on) and `ESA_CCI_SST` (ESA CCI SST v3.0, 1980 on), each with its own 1981-2010
climatology. `enso/01_enso_plume.ipynb` plots HadISST 1.1, ERSSTv5, ESA CCI SST and NOAA OISST v2.1 and keeps
HadISST 1.1 as the reference for errors. ERA5-Land has no ocean
points, so it gives no Nino index.

## Files

The notebooks are in `jupyter/`, grouped by topic (notebook paths below are relative to `jupyter/`); within a
folder the number is the reading order. Every notebook is self-contained: it reads the model output it
needs, builds and saves the caches that are missing or stale, reuses valid ones, and runs its diagnostics
(see "Caches" under Running); there is no build step to run first.

Each notebook finds `util/` by itself (its first code cell adds `diagnostic/` to `sys.path`), so it runs from its own folder.
Batch scripts are in `scripts/`; run output (logs, executed copies) goes to `/pscratch/sd/z/zhan391/run_logs/`.


### Shared code (`util/`)

| File | Purpose |
| --- | --- |
| `util/rtenso_case.py` | Shared helpers, imported by every notebook as `from util import rtenso_case as rc`: case registry, hindcast loading, lead-time climatologies, observation products and anomalies, an optional pre-build of the reference caches (`build_caches`), pattern ACC/RMSE, regions, velocity potential, output paths. |
| `util/xro_case.py` | XRO state vector (Nino3.4, recharge *h*, NPMM, SPMM, IOB, IOD, SIOD, TNA, ATL3, SASD boxes) and the configuration shared by the XRO notebooks (`enso/06`-`09`) (`EVENTS`, `CASES`, `XRO_SETUP`); cached observed and hindcast state vectors, training-period trend, XRO fit, event forecasts aligned with the hindcast starts, initial-state sensitivity experiments and out-of-sample rolling hindcasts. |
| `util/pubfig.py` | Publication figure style for the XRO notebooks: figures drawn at journal column width (3.5/7.2 in) with 7-9 pt text, saved as vector PDF (embedded fonts) and 600-dpi PNG. |
| `util/strat_case.py` | Stratosphere readers: daily zonal means and daily 10-hPa maps from the native `eam.h1` files, monthly zonal means of the post-processed U/T/Z3, ERA5 daily zonal means and 00Z maps, lead-time climatologies of zonal means (default `start_drift`; `direct`, `lead_drift`), the daily-minus-monthly Z offset, SSW central dates (Charlton and Polvani 2007), final warming with coverage status, event details, lag composites. |
| `util/strat_diag.py` | Shared inputs of the stratosphere notebooks: `StratConfig` (data roots and run controls), `load_era5` (cached ERA5 indices and climatologies), `load_references` (Reanalysis zonal means of every start month), `load_model` (pooled 2026 starts with their daily and monthly indices), `ssw_catalogue`, `djf_vortex`, DJF helpers. |
| `util/vortex_plot.py` | North-polar map helpers of `stratosphere/06`-`08`: projection and extent, fixed-level anomaly shading and total contours, U60 badges, video writer. |

### `enso/`: the El Niño forecasts, their verification and ENSO dynamics

| File | Purpose |
| --- | --- |
| `enso/01_enso_plume.ipynb` | Nino3.4/3/4/1+2/ONI forecast plumes per start against HadISST 1.1, ERSSTv5, ESA CCI SST and NOAA OISST v2.1 with the observed peak period shaded, IRI-style seasonal plumes, a CPC-style RONI outlook (Relative Oceanic Nino Index: 3-month Nino3.4 minus 20S-20N tropical-mean anomaly, rescaled to the ONI's standard deviation with HadISST 1.1 1981-2010; pooled-ensemble 5-95/15-85/25-75% ranges), peak-season summary, ensemble-mean error and spread (after ESP-Lab 3a), and a map of the Nino regions. |
| `enso/02_anomaly_maps_acc.ipynb` | Peak-season (DJF) anomaly maps from every start, and spatial pattern ACC and RMSE by lead season and region, per member and ensemble mean (after ESP-Lab 1a/1b); only fully observed seasons are scored. Set `FIELD` to `TS`, `PRECT` (both events), `TREFHT`, `PSL` or `U850` (2015 only). |
| `enso/03_teleconnections.ipynb` | Peak-season teleconnections of `PRECT`/`TREFHT`/`PSL` per degC of Nino3.4 against the observed event and the 1982-2010 observed regression, and CPC-style equatorial Hovmollers of SST, U850, OLR and 200-hPa velocity potential (`rc.velocity_potential`) (after ESP-Lab 3b). |
| `enso/04_npac_namerica.ipynb` | CPC-style North Pacific / North America maps (120E-60W, 20-70N) for Nov-Feb: 500-hPa height and anomaly, 200-hPa wind, 925-hPa temperature anomaly, observed (ERA5, through 2019) and each case's ensemble mean. |
| `enso/05_us_temp_precip.ipynb` | CPC-style U.S. maps: 1- and 3-month percent of average precipitation and temperature departures, and a tercile-probability seasonal outlook (NDJ) from the pooled members with the observed category. |
| `enso/06_xro_validation.ipynb` | XRO (Zhao et al. 2024) validated against observations only: training-period trends (linear vs quadratic *h*), fitted recharge-oscillator growth rate and period, stochastic simulations against observed seasonal synchronization, *h*-*T* lag correlation and persistence, out-of-sample skill (refit every year 1995-2025) by lead and initial month against persistence, DJF amplitude (does XRO underpredict strong El Ninos?), and calibration of the stochastic ensemble. |
| `enso/07_xro_events.ipynb` | XRO forecasts of the 2015/16 and 2026/27 events from the months before each start: Nino3.4 and *h* plumes, peak-season summary, peak-season Nino3.4 sensitivity to each initial-state component (removal experiments, with the non-additive remainder shown), and to the initial *h* alone. |
| `enso/08_e3sm_xro_state.ipynb` | E3SM hindcasts in the XRO state space against observations and XRO: state-vector plumes (common 5-95% ranges), error decomposition into the start-month error B(0) and the error accumulated afterwards G(L) = B(L) - B(0), summary heat maps, growth of Nino3.4 and *h* from the start month, peak-season summary. |
| `enso/09_e3sm_recharge.ipynb` | *T*-*h* phase planes, and XRO started from E3SM's start-month offsets (*h* only, and the full state) to test whether the E3SM-XRO peak difference comes from the initial state or from E3SM's evolution. |
| `enso/10_enso2026_animation.ipynb` | Animated month-by-month evolution of the 2026/27 El Nino: global maps (Robinson; `DOMAIN` for a regional box) of SST (HadISST 1.1), 0-700 m heat content (EN4.2.2), precipitation (GPCP v2.3) and land 2-m air temperature (ERA5-Land) anomalies, observed and the pooled (`clim`/`obs`/`NSIDC`) ensemble mean of each 2026 start, above a Nino3.4 plume drawn up to the frame month, with the deterministic XRO forecast from each start (`07_xro_events` control setup) as a dashed reference. Saves a GIF and a video per field (H.264 MP4 if ffmpeg has libx264, e.g. via `imageio-ffmpeg`; else VP9 WebM) and a static storyboard of selected months; the pooled ensemble means are cached under `data/`. |

### `stratosphere/`: polar vortex and SSWs in the 2026 starts

| File | Purpose |
| --- | --- |
| `stratosphere/01_model_climatology.ipynb` | E3SM's stratosphere against ERA5: seasonal cycle of U(60N, 10 hPa) and polar-cap temperature, DJF zonal-mean U bias, DJF vortex strength and interannual variability. |
| `stratosphere/02_vortex_forecast.ipynb` | The vortex forecast: Feb 2026 start against the observed end of winter 2025/26, daily U(60N, 10 hPa) plumes for 2026/27 and 2027/28, the vortex index plume (U60, polar-cap T and height, polar-cap SLP), DJF tercile probabilities with their sensitivity to the climatology method. |
| `stratosphere/03_ssw_occurrence.ipynb` | SSW frequency by start x case x winter (raw wind; bias-adjusted as sensitivity) against ERA5 1979/80 on, central-date months, final warming with censoring; event catalogue CSV. |
| `stratosphere/04_enso_qbo_background.ipynb` | QBO forecasts (time-height and 50 hPa), the ENSO/QBO background state by start and case, within-group member covariation with the DJF vortex. |
| `stratosphere/05_ssw_dynamics.ipynb` | Season-matched SSW lag composites (U60, 100-hPa heat flux, polar-cap height) with cluster-bootstrap intervals, event-property distributions, vertical evolution of polar-cap anomalies and SLP. |
| `stratosphere/06_ssw_mechanism.ipynb` | What distinguishes an SSW member from comparable members that do not reverse? A featured SSW member (`EVENT`) against a matched non-SSW member and a near-miss: aligned U60 / polar-cap T / heat-flux timelines, ensemble ranks, the ensemble relationship of preconditioning and pre-onset wave activity with the outcome, pre-onset 10-hPa maps, side-by-side animation. |
| `stratosphere/07_vortex_ssw_animation.ipynb` | Day-by-day 10-hPa polar maps of one winter: each start's ensemble mean and a featured SSW member, the all-start mean (or ERA5 for an observed winter), with a U60 timeline; GIF, video, storyboard. |
| `stratosphere/08_vortex_evolution_animation.ipynb` | Seasonal evolution Oct 2026-Apr 2027 without event-based selection: one fixed member and the ensemble mean of each start, raw 10-hPa height contours with optional anomaly shading, reference climatologies, U60 timeline; GIF, video, preview, storyboard. |

### `nov2026_outlook/`: the 1 Nov 2026 case (36-month start) and its La Niña outlook

| File | Purpose |
| --- | --- |
| `nov2026_outlook/01_year3_drift_sensitivity.ipynb` | Year 3 of the 36-month start: the drift by lead and its year-1-to-2 change, where the lead-24 step comes from, year-3 drift extrapolations (held, per-month, adjusted for the initial-adjustment leads; r = 0 to 1), standard vs relative Niño3.4, 1981-2010 vs 1991-2018/2020 bases, and P(La Niña) for DJF 2027/28, 2028/29 and both from the 30 pooled members. |
| `nov2026_outlook/02_outlook_checks.ipynb` | Checks behind the outlook: F-test (and Levene) of the DJF means across the clim/obs/NSIDC sets (pooling), Niño3.4 DJF anomaly correlation of the Reanalysis 1 Nov starts in years 1 and 2 against persistence, and the peak and decay timing of the forecast against observed strong El Niños and the other 2026 starts. |
| `nov2026_outlook/03_drift_trend_and_test.ipynb` | Trend-aware drift correction (T1: linear trend in the drift; T2: Kharin et al. 2012, trend in the model climatology) against fixed 1981-2010, 1991-2018 and recent-window drifts; a test in time (drift from the 1981-2008 starts, verified on 2009-2018: bias, RMSE, La Niña Brier score); P(La Niña) of the 2026 forecast for each correction; and the 2026 forecast's sensitivity to the fixed drift window (1981-2010 default against sliding 15-year windows, 1981-1995 to 2004-2018); and the trend correction's sensitivity to the training starts (T1 held at 2018, fitted on 1981-2018 against windows starting later and ending in 2018), with its test in time (windows ending in 2008, verified on 2009-2018). Member-level index cache in the same folder. |

The notebooks share an index-space cache (Niño3.4 and the tropical ocean mean of the Reanalysis Nov starts, the 30 members and HadISST) under `$S2D_DIAG_ROOT/RTENSO2026_pool/nino_drift/`.

### `initialization/`: sensitivity to the land, sea-ice and river ICs

| File | Purpose |
| --- | --- |
| `initialization/01_ic_sensitivity.ipynb` | `Reanalysis-CLIM` vs `-OBS`, `-OBS` vs `-NSIDC`, and each vs `Reanalysis-FOSIRL` (2015 only), for both events: land-mean differences with member-based intervals, fraction of land with significant differences by lead, stippled difference maps, and hemispheric sea-ice area against NOAA/NSIDC CDR v6. |

## Running

Use the `e3sm_analysis` kernel (ESP-Lab installed in editable mode); start Jupyter in this directory
or any subfolder (each notebook adds `diagnostic/` to `sys.path` so `util` imports). Run a notebook
interactively, or headless with its executed copy written outside this tree:

```bash
conda activate e3sm_analysis
jupyter nbconvert --to notebook --execute --output-dir=/pscratch/sd/z/zhan391/run_logs/rtenso2026_notebook_runs \
    jupyter/enso/01_enso_plume.ipynb
```

**Caches.** Each notebook's run-control block holds `FORCE_COMPUTE` and `CLIM_METHOD`:

- `FORCE_COMPUTE = False` (default): read valid caches, build missing or stale ones. A cache is stale when
  its key (inputs, start length, observation file times, ...) no longer matches, so new data or a longer
  run is picked up without forcing.
- `FORCE_COMPUTE = True`: recompute this notebook's own products (the hindcast-derived caches it uses).
- Shared reference caches (Reanalysis 1981-2010 lead-time climatologies, observations, ERA5, the XRO
  observed state and rolling hindcasts) are reused even with `FORCE_COMPUTE`; set
  `rc.REBUILD_REFERENCE = True` in the notebook to recompute them too (`rc.reference_mode`).
  `rc.build_caches(mode="rebuild")` pre-builds all of them at once (about 4 minutes per field and start
  month).
- Caches are shared between notebooks (`$S2D_DIAG_ROOT`); whichever notebook needs one first builds it.

Inputs: the RTENSO2026 cases read their time series from
`/global/cfs/cdirs/e3sm/S2S2D/RTENSO2026/post_process` (`rc.RTENSO2026_POST_ROOT`);
the Reanalysis (BruteForce) cases read from `RAW_MODEL_ROOT`
(`/global/cfs/cdirs/e3sm/S2S2D/post_process`).

Outputs:

- Notebook products: one folder per notebook (mirroring the folders here) under the fixed
  `/global/cfs/cdirs/e3sm/www/zhan391/e3sm-rtenso2026/` (`rc.OUTPUT_ROOT`), kept apart from
  the ESP-Lab gallery. Each notebook sets `NOTEBOOK` and writes through `rc.save_figure` and
  `rc.output_dir`:

  ```text
  e3sm-rtenso2026/
    enso/01_enso_plume/{figures,data}/            data: Nino index netCDF per case/start, peak summary CSV
    enso/02_anomaly_maps_acc/{figures,data}/      data/<FIELD>/: pattern-skill netCDF, ACC CSV
    enso/03_teleconnections/{figures,data}/       data: teleconnection ACC CSV
    enso/04_npac_namerica/figures/
    enso/05_us_temp_precip/figures/
    enso/06_xro_validation/{figures,data}/        data: RO parameters, skill, DJF amplitude and calibration CSVs
    enso/07_xro_events/{figures,data}/            data: XRO event forecasts (netCDF), peak and sensitivity CSVs
    enso/08_e3sm_xro_state/{figures,data}/        data: start-month and evolution error CSVs, peak summary
    enso/09_e3sm_recharge/{figures,data}/         data: matched-initial-state CSV
    enso/10_enso2026_animation/{figures,data}/    figures: GIF/video animations and storyboards; data: pooled ensemble-mean caches
    stratosphere/01_model_climatology/ ... stratosphere/08_vortex_evolution_animation/   figures, CSV tables
    initialization/01_ic_sensitivity/{figures,data}/
  ```
- Reusable caches stay outside it: `$S2D_DIAG_ROOT/Reanalysis/leadtime_clim/`,
  `$S2D_DIAG_ROOT/<source>/strat/{daily_zonal,monthly_zonal,daily_maps}/` and `observations/strat/daily_maps/` (the stratosphere notebooks) and
  `$S2D_DIAG_ROOT/observations/case_study/` (`S2D_DIAG_ROOT` defaults to
  `/global/cfs/cdirs/e3sm/S2S2D/s2d_diag`).

## XRO recharge proxy (`OHC700`)

The hindcasts save no subsurface temperature, so the XRO notebooks (`enso/06`-`09`) use the 5S-5N, 120E-80W
mean 0-700 m heat content anomaly as the recharge *h* (instead of the 20 C isotherm depth).
MPAS-O writes `oceanHeatContentSfcTo700m` per native cell (J, integrated with T in K);
`rc.open_member` divides it by the native cell area remapped with the file's own `map_file`
weights (`rc.mpas_cell_area`), giving GJ m-2. Observations are EN4.2.2 OHC700
(`data/ENv4.2.2/1x1`, 1979 to Sep 2026). Anomalies remove the 1981-2010 climatology as for
the other fields; the constant offset between K- and C-based integrals cancels at fixed depth
(if MPAS integrates from the moving sea surface, sea-level anomalies add about 0.1 GJ m-2 per
0.1 m in phase with recharge; not yet checked). The XRO state vectors and rolling hindcasts are
cached under `$S2D_DIAG_ROOT/observations/xro/` and `$S2D_DIAG_ROOT/<source>/xro_state/`; rerun
with `rc.REBUILD_REFERENCE = True` to force them (they also rebuild by themselves when the observation files change). The XRO
package (`pip install XRO`, examples in `~/code/XRO`) must be installed in `e3sm_analysis`.

## Stratosphere series (`stratosphere/`)

Each notebook builds the stratosphere caches it needs; the first one to run (or one with `FORCE_COMPUTE`)
reads the native daily model output, about 1 h, so run that one on a compute node. After that they only
read caches and can run in any order. `stratosphere/02_vortex_forecast` also checks the one correction
applied to the model data (the daily-minus-monthly 10-hPa height offset, about 38 m). Each starts with
its question and findings. The polar-vortex notebooks need `ERA5/1x1/daily_zonal` in the reference
archive (built by `ERA5/scripts/daily_zonal`). On a login node (30 GB per user) run one stratosphere notebook at a time;
`scripts/run_strat_series.sbatch` runs the whole series on a compute node (executed copies and logs go to
`/pscratch/sd/z/zhan391/run_logs/rtenso2026_notebook_runs/`, not into this tree).

The notebooks share their loaders in
`util/strat_diag.py`; edit the notebooks directly.

Model anomalies use `CLIM_METHOD = "start_drift"` as in the other notebooks (see "Drift correction"):
the ERA5 target-month climatology plus the drift of the Reanalysis starts of the same month (May/Nov
1981-2010, Feb/Aug the six 1985-2010 starts). The earlier `lead_drift` climatology produced a spurious
weak-vortex signal for the Feb/Aug starts; `stratosphere/02_vortex_forecast` compares the methods.
