# `hind_lead_eval`: data and methods

How the four notebooks turn the F20TR initialization hindcasts and the observations into skill scores,
and why each choice was made. Code: `util/fcst_init_data.py` (fields and climatologies),
`util/fcst_init_skill.py` (scores), `util/fcst_init_plotting.py` (figures),
`configs/fcst_init_experiment_config.py` (experiments). Last revised 2026-11-05.

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

Lead day 1 is the start date (00-24Z of the initial day); scores use lead days 1-56 (8 weeks). A
configuration-season is used only when all 20 members report a finished post-processing
(`regen_ts_<stream>.status` OK; the zppy files stop after 4 days, so `post_process/regen_f20tr_ts.sh`
regenerates them).

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

* **Drift**: not removed by either bias treatment (section 5); week 1-2 comparisons between configurations
  include it.
* **Verifying analysis**: ERA5 also provides the initial atmosphere of HYB_BASE and DA_LND, which favours them
  in week 1.
* **Significance**: the t-test treats the 8 starts as independent, but starts 5 days apart verify on
  overlapping periods, so significance is overstated. Lead with gridpoint ACC within a season and ΔACC.
* **00Z fields**: verified at 00Z only (no model daily means), with all-hours climatologies (section 4).
* **Precipitation**: GPCP V3.3 is built from IMERG V07B within 55N-55S, so it is not independent of IMERG
  there.

## 9. Outputs and caches

Under `DIAGNOSTIC_OUTPUT_ROOT/{data,figure}/hind_lead_eval/`: `lead_skill/` (01), `s2s_weekly/` (02),
`land/` (03), `mjo/` (04). Climatologies and anomalies are cached and reused unless `FORCE_COMPUTE`; a cached
observed climatology is reused only if it was built from monthly means of the daily series for the current
`CLIM_YEARS`, and an anomaly file only if its `bias` and `clim_years` attributes match. Results built from the
observations used before 2026-11-05 were moved to `data/_stale_obs_20261105/`.
