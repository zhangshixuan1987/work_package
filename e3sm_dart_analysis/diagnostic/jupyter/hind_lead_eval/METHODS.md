# `hind_lead_eval`: data and methods

How the notebooks turn the F20TR initialization hindcasts and the observations into skill scores,
and why each choice was made. Code: `util/fcst_init_data.py` (fields and climatologies),
`util/fcst_init_skill.py` (scores), `util/fcst_init_plotting.py` (figures),
`configs/fcst_init_experiment_config.py` (experiments), `util/init_benefit.py` (05a-05c), `util/heatwave_case.py`
(06a-06h), `util/winter_storm_case.py` (07a-07d); the case settings of 05-07 are in `configs/hind_lead_eval_cases.py`
and their shared, cached products in `util/*_workflow.py`. Last revised 2026-11-08.

## 1. Hindcasts

Five initialization configurations, 20 members each, run for about 62 days from 4 start dates per season
(16, 21, 26, 31 December 2011 and 16, 21, 26, 31 May 2012):

| Configuration | Atmosphere | Land |
|---|---|---|
| CLIM_BASE | free-running climate state | free-running |
| HYB_BASE | ERA5-based | DATM-forced |
| DA_ATM | DART analysis | as HYB_BASE |
| DA_LND | as HYB_BASE | DART analysis |
| DA_COUPLED | DART analysis | DART analysis |

Both start sets hold the same experiments: 5 configurations x 4 starts x 20 members. The members of each
ensemble are exchangeable and their numbers are labels only. The BASE configurations are EN01-EN20, with EN01
the unperturbed control and EN02-EN20 the control plus the centred perturbations of 19 of a season's 20
selected DART members (`ensemble_perturbations/02_ic_free_clim_dartpert.ipynb` and `03_ic_era5_dartfill_dartpert.ipynb`). The DA
configurations keep the DART numbers of the same 20 members, a different subset of the 40 in each season
(Dec 2011, `DARTEN40`: EN01 EN03 EN04 EN07 EN08 EN09 EN10 EN14 EN17 EN18 EN20 EN21 EN25 EN26 EN28 EN31 EN35
EN36 EN37 EN40; May 2012, `DARTEN40S1`: EN02 EN03 EN05 EN07 EN11 EN12 EN13 EN14 EN15 EN16 EN20 EN22 EN24
EN26 EN27 EN28 EN31 EN32 EN34 EN38), so each configuration's list comes from `members(config, season)`.
Members are never paired across configurations by number.

Lead day 1 is the start date (00-24Z of the initial day); scores use lead days 1-56 (8 weeks). A
configuration-season is used only when all 20 members report a finished post-processing
(`regen_ts_<stream>.status` OK; zppy's own files stop after 4 days, so `post_process/submit_process_e3sm_post.sbatch`
makes full-length ones). Every configuration must be ready in both seasons (`balanced_runs`); otherwise the
notebooks stop rather than pool an unbalanced set.

### Ocean and sea-ice boundary

All runs are AMIP-type (interactive atmosphere and land) hybrid starts with the same prescribed SST: DOCN
`SSTDATA` mode, the NOAA AVHRR 1x1 file, identical stream settings and records; no run reads a DOCN restart.
Open-ocean TS on the first days is identical across configurations in 99 % of the fully ocean cells within
45S-45N. The reference case differs: the DATM case for HYB_BASE, DA_ATM and DA_LND, the DART case for
DA_COUPLED, so DA_COUPLED's sea-ice state also comes from DART; the remaining TS differences are in
sea-ice and marginal-ice cells (prognostic sea-ice surface temperature).

### What the model output contains

The output is set in each case's `atm_in` (identical for all configurations):

| Stream | Frequency | Contents relevant here |
|---|---|---|
| h1 | daily, one record | **instantaneous** 00Z state (`:I`, `cell_methods = "time: point"`, stamped 00Z at the end of the day) of every pressure-level field (U, V, T, Z, Q, OMEGA at 19 levels) and of TREFHT, PSL, TS, TMQ; daily **means** (`:A`) only of FLUT, PRECT, TCO, SCO |
| h2 | 6-hourly | 6-hour **means** of TREFHT, TS, TMQ, surface fluxes, clouds, radiation; no pressure-level fields |
| h0 | monthly | full 3-D U, V, T, Z3 as monthly means |

So a true daily mean exists for the surface and column fields (from h2) and for FLUT and PRECT (h1), but
**not for any pressure-level field**: those exist only as one 00Z state per day.

## 2. Observations

All observations come from the reference-data archive. Its location is a notebook input
(`OBSERVATION_ROOT` in the setup cell, with `LAND_SUPPORT_DIR` and `AMIP_ROOT` below it); the code holds only
paths relative to it and has no default location.

| Fields | Product | Archive path |
|---|---|---|
| TREFHT, TS, PSL, TMQ, surface fluxes (daily means) | ERA5, daily means of the 6-hourly 1-degree archive: exact for fluxes, trapezoid rule over 00Z..00Z of the next day for states | `ERA5/1x1/daily` |
| Z500, T850, U850, U200 (00Z) | ERA5 00Z analyses | `ERA5/1x1/6hourly` |
| PRECT | GPCP V3.3 daily (2 x 2 area-weighted mean of the 0.5-degree files) | `GPCPv3.3/1x1/daily` |
| FLUT | NOAA OLR Climate Data Record v02r00 | `NOAA_CDR/1x1/daily` |
| Top-5-cm soil moisture | ESA CCI | `ESA_CCI_H2OSOI/1x1/daily` |
| RMM check | BoM RMM index | `MJO_INDEX/rmm.74toRealtime.txt` |

Not used: the 1-degree E3SM Diags file `obs_for_e3sm_diags/time-series/GPCPDAY_3.3/precip_199801_202312.nc`,
whose latitudes are reversed (r = 0.02 with IMERG as stored, 0.86 flipped), and the older GPCP V3.1 files
`GPCPv3.3/1x1/daily/GPCP_daily_<YYYY>_1x1.nc`. GES DISC describes GPCP V3.3 as not yet fully validated and
asks users to contact George Huffman before publishing with it.

## 3. Daily values on the lead days

Every field is put on (start date, lead day, lat, lon), 1-degree grid:

* **Daily-mean fields** (TREFHT, TS, PSL, TMQ, fluxes): model = mean of the four 6-hour means of the day (h2);
  ERA5 = `1x1/daily`. FLUT and PRECT: model daily means (h1); GPCP and NOAA OLR daily.
* **Pressure-level fields** (Z500, T850, U850, U200): the 00Z state ending each day, model and ERA5 at the
  same instant. Instant against instant, so the verification itself is consistent.
* **Land**: ELM daily TSA and top-5-cm soil moisture (layer-thickness weighted) against ERA5 TREFHT and
  ESA CCI.

## 4. Climatologies (notebooks 02 and 03)

Both climatologies are built the same way, for every variable:

1. calendar-month means over `CLIM_YEARS` = 2001-2014;
2. linear interpolation between mid-month values (day 15) to each verification day (midday for daily means,
   the 00Z instant for 00Z fields). Interpolation avoids steps at month boundaries; the same function is
   applied to both sides.

* **Observed**: monthly means of the same daily observations as above (ERA5 `1x1/daily`, GPCP, NOAA OLR,
  ESA CCI), so the climatology comes from the same product as the verifying days.
* **Model**: monthly means of the free-running E3SMv3 AMIP runs `v3.LR.amip_0101` and `v3.LR.amip_0151`
  (same grid, observed SST and sea ice; `AMIP_ROOT`). Four starts per season are too few for the hindcasts'
  own lead-dependent climatology, the S2S standard.

Both are all-hours means, also for the 00Z fields (the AMIP runs have no 00Z output), so each 00Z anomaly
carries its own system's 00Z departure from the daily mean (ERA5: RMS about 6 m Z500, 0.5 K T850,
0.5 m s-1 U850/U200). This shared pattern slightly raises pattern ACC; it does not affect gridpoint ACC
within a season.

## 5. Weekly anomalies and the two bias treatments

Week w is the mean of lead days 7(w-1)+1 to 7w (weeks 1-8; a week needs all 7 days). The forecast uses the
ensemble mean. The observed anomaly is always observations minus the observed climatology. The forecast
anomaly is computed both ways (`BIASES`):

| `BIAS` | Forecast anomaly | Counts as error |
|---|---|---|
| `none` | ensemble mean minus the **observed** climatology | initialization drift and the model mean bias (the full-field forecast) |
| `amip` | ensemble mean minus the **AMIP** climatology | the drift, but not the model mean bias (the S2S convention of separate climatologies) |

Neither removes the drift. Under `amip` a forecast still close to the observed state (weeks 1-2) is penalized
by the opposite of the model bias. The difference `amip - none` (figures `s2s_bias_impact_*`,
`land_bias_impact_*`) shows the effect of the model mean bias on skill: positive once the forecast has drifted
to the model climate, negative in the first weeks; where it crosses zero indicates each configuration's drift
time scale. `loo` (mean error of the other starts) is kept only for comparison: starts 5 days apart verify on
overlapping days, so it leaks their observed anomalies into the forecast.

## 6. Scores

| Notebook | Scores |
|---|---|
| 01 lead-day performance | per lead day and region, pooled over the 4 starts of a season: ensemble-mean RMSE and bias, spread (unbiased member variance x (N+1)/N), spread/RMSE, fair CRPS; window means (days 1-3, 4-10, weeks 3-4, 5-8). Full fields, no anomalies. |
| 02 weekly S2S skill | per bias treatment: pattern ACC (centred, cos-lat weighted) and MSSS by region and week; gridpoint ACC over the starts (pooled over the available seasons) with a one-sided t-test; ΔACC against `HYB_BASE`; scorecard; bias impact. |
| 03 land | lead-day scores over land; weekly skill as in 02 for both bias treatments; soil-moisture memory (week-1 with week-w soil moisture) and SM-TSA, SM-LHFLX coupling, correlations over the 80 member-starts of a season (departures from the ensemble mean, so independent of the climatology). |
| 04 MJO | RMM index of each member (section 7); bivariate COR and RMSE of the ensemble mean by lead (Rashid et al. 2011), spread, amplitude error, lead-window skill. |
| 05a-05c initialization benefits | 2 x 2 roles (atmosphere, land DART-initialized: E00 HYB_BASE, E10 DA_ATM, E01 DA_LND, E11 DA_COUPLED): pattern ACC of TREFHT over 30-60N land and the contrasts E11-E01 (atmosphere), E11-E10 (land) and E11-E10-E01+E00 (land-atmosphere consistency), weekly and on 7-day running means, both bias treatments; paired bootstrap over starts and members; panel against a held-out E11 member. Incremental benefits of DART over the reference initialization, not total component contributions. |
| 06a-06h heat wave 2012 | case study of the central U.S. heat wave (35-45N, 105-85W land; May starts): event from ERA5 daily-mean TREFHT above its calendar-day 90th percentile (2001-2014, +-7 days) for 3+ days (27 Jun-7 Jul 2012, peak 4 Jul, +5.1 K); coverage (leads 28-53); regional series of all members; member metrics (event-mean and peak anomaly, peak-date error, exceedance days, longest run, accumulated exceedance, area fraction, event precipitation), both bias treatments; maps. Organized as Q1 (HYB_BASE, DA_COUPLED vs CLIM_BASE), Q2 (DA_COUPLED vs HYB_BASE), Q3 (DA_ATM, DA_LND). |
| 07a-07d winter storm 2012 | case study of the 16-20 Jan 2012 Pacific Northwest storm (Dec starts; WEST_WA_OR and INLAND_PNW land): regional series of the 80 pooled members, member metrics (event and maximum daily precipitation, peak date, event-mean and minimum temperature, precipitation on freezing and near-freezing days, pattern scores), event maps, storm environment (IVT, Z500, PSL); both bias treatments. |

## 7. MJO: winds for the RMM index

The Wheeler-Hendon RMM index combines daily-mean OLR, U850 and U200 (15S-15N), projected on the observed
EOFs built by `hind_atm_mjo/03_compute_and_plot_rmm_skill.ipynb` (1986-2017), with the previous-120-day mean
taken from the observations before each start. The winds are required (they carry most of the RMM variance)
and need daily means.

* **Observed index**: true ERA5 daily means (`ERA5/1x1/daily`) and NOAA OLR CDR daily means.
* **Model index**: OLR daily means (h1); U850/U200 exist only as 00Z states (section 1), so the daily mean is
  estimated as the mean of the two 00Z states bounding the day (`open_member_daily(pair_mean=True)`; lead day 1
  uses the initial state). This estimate is centred on midday, like a true daily mean.

How much this matters, tested on the observed index for 2012-2013 (730 days), against the index built from
true daily-mean winds:

| Daily U850/U200 in the index | Bivariate COR | RMSE |
|---|---|---|
| mean of the two 00Z states (used) | 0.9995 | 0.030 |
| single 00Z state ending the day | 0.994 | 0.105 |
| single 00Z state starting the day | 0.995 | 0.098 |

At individual points the two-00Z estimate differs from the true daily mean by RMS 0.9 m s-1 (850 hPa) and
1.8 m s-1 (200 hPa) in 15S-15N, against daily variability of 4.0 and 9.9 m s-1, but with a near-zero mean;
the projection on the large-scale EOFs averages it out. An index RMSE of 0.03 is negligible next to forecast
RMM errors of 0.5-1.5 and the usual skill thresholds (COR 0.5, RMSE about 1.4), so the evaluation is sound.
The observed index keeps the standard true daily means. For future hindcasts, `U850:A` and `U200:A` (and V,
Z500, T850) in the daily stream would give true model daily means at little cost.

## 8. Limitations

* **Incomplete ERA5 initial atmosphere (HYB_BASE, DA_LND): their results are not valid as they stand.** The
  ERA5 file `DART_INIT/ATM/v3.LR.ne30pg2.ERA5.eam.i.<date>.nc`, and the `HYB_DARTPERT20` files built from it,
  hold 13 variables (U, V, T, Q, PS, CLDLIQ, CLDICE and coordinates) against 129 in the DART and FREE_CLIM
  files: no ozone, no linoz gases (H2OLNZ, CH4LNZ, N2OLNZ), no chemistry, no aerosols, no sea-ice surface
  fields. EAM starts these near zero, and the radiation reads its gases from those tracers
  (`rad_climate = 'A:H2OLNZ:H2O', 'A:O3:O3', 'A:N2OLNZ:N2O', 'A:CH4LNZ:CH4'`), so the runs start
  radiatively dry and without ozone (day 1 against DA_ATM: -190 W m-2 downward longwave over land, -250 DU
  ozone). Global clear-sky downward longwave at the surface (FLDSC) is then
  94 W m-2 on day 1 against 318 in DA_ATM, recovers only over about 8 weeks, and the lower stratosphere
  (50 hPa) cools by 11 K. Over land TREFHT is 6-7 K too cold for 2-4 weeks in both seasons (water vapour
  is normal; the ocean is held by the prescribed SST). Contrasts with HYB_BASE or DA_LND, including
  notebook 05's contrasts involving E00 or E01 (atmosphere, consistency, total), mostly measure this defect, not the initialization.
  Found 2026-11-06. The completed files (ERA5 fields plus the DART-mean of every other field, H2OLNZ from
  ERA5 Q) are built by `ensemble_perturbations/03_ic_era5_dartfill_dartpert.ipynb` (step 3a); these two
  configurations need new runs from them.
* **Drift**: not removed by either bias treatment (section 5); week 1-2 comparisons between configurations
  include it.
* **Verifying analysis**: ERA5 also provides the initial atmosphere of HYB_BASE and DA_LND, which favours them
  in week 1.
* **Significance**: the t-test treats the 8 starts as independent, but starts 5 days apart verify on
  overlapping periods, so significance is overstated. Lead with gridpoint ACC within a season and ΔACC.
* **00Z fields**: verified at 00Z only (no model daily means), with all-hours climatologies (section 4).
* **DA_COUPLED**: started from the DART case, so its sea-ice state also differs (polar only; section 1).
* **Precipitation**: GPCP V3.3 is built from IMERG V07B within 55N-55S, so it is not independent of IMERG
  there.

## 9. Outputs and caches

Under `DIAGNOSTIC_OUTPUT_ROOT/{data,figure}/hind_lead_eval/`: `lead_skill/` (01), `s2s_weekly/` (02),
`land/` (03), `mjo/` (04), `init_benefit/` (05), `heatwave_2012/` (06), `winter_storm_2012/` (07). In 05-07 every product
(model and observed fields, metrics, scores) is saved with the settings it was computed from and read back when they are
unchanged; `FORCE_COMPUTE` recomputes a notebook's results and `FORCE_INPUTS` also re-reads the model output and
observations. In 01-04 climatologies and anomalies are cached and reused unless `FORCE_COMPUTE`; a cached
observed climatology is reused only if it was built from monthly means of the daily series for the current
`CLIM_YEARS`, and an anomaly file only if its `bias` and `clim_years` attributes match. Results built from the
observations used before 2026-11-05 were moved to `data/_stale_obs_20261105/`.
