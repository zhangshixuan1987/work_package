# ERA5 initial-condition tests

Evaluates the test hindcasts that address the issues found in the ERA5-initialized F20TR hindcasts.
Every test is a single member started on 16 Dec 2011 00Z and run for 2 months. Each is compared with the
production members (`INPUT_DATA_ROOT`, eam_hindcast) that share its other initial states.

## Issues

1. **Incomplete ERA5 atmosphere.** `DART_INIT/ATM/v3.LR.ne30pg2.ERA5.eam.i.<date>.nc` holds only T, U, V, Q, PS,
   CLDLIQ and CLDICE, and EAM starts every missing field from zero. The radiation reads its gases from tracers
   (`rad_climate = 'A:H2OLNZ:H2O', 'A:O3:O3', 'A:N2OLNZ:N2O', 'A:CH4LNZ:CH4'`), so HYB_BASE and DA_LND start
   radiatively dry and without ozone, aerosols or hydrometeor numbers, about 6 K too cold over land.
2. **Source of the filled fields.** The completed ERA5 states (`ensemble_perturbations/util/complete_era5_ic.py`)
   take the missing fields from the DART ensemble mean or from the free-running CLIM_BASE (FREE_CLIM) state.
3. **Component contrasts confounded by the defect.** The production DA_LND (bare ERA5 atmosphere) and DA_ATM
   (DATM-forced land) mix the IC defect, or a different land spin-up, into the atmosphere and land contrasts.
   The CLIMO tests replace one component of DA_COUPLED EN40 with the complete CLIM_BASE state (Berner-style).

## Runs (`ic_util/ic_test_runs.py`)

| Key | Case (`_2011121600`) / member | Atmosphere IC | Land IC |
|---|---|---|---|
| ERA5_BARE | HYB_BASE / EN01 (production) | bare ERA5 | DATM EN01 |
| ERA5_DARTFILL | ERA5COMPLETE_TEST / EN01 | ERA5 + DART-mean fields | DATM EN01 |
| ERA5_CLIMFILL | ERA5CLIMFILL_TEST / EN01 | ERA5 + CLIM_BASE fields | DATM EN01 |
| DA_COUPLED | DA_COUPLED / EN40 (production) | DART EN40 | DART EN40 |
| DA_ATM | DA_ATM / EN40 (production) | DART EN40 | DATM EN40 |
| DA_LND | DA_LND / EN40 (production) | bare ERA5 | DART EN40 |
| CLIMO_ATM | CLIMO_ATM_TEST / EN40 | CLIM_BASE + (DART EN40 - DART mean) | DART EN40 |
| CLIMO_LND | CLIMO_LND_TEST / EN40 | DART EN40 | CLIM_BASE |

Comparisons (A - B): `tracer_fix` ERA5_DARTFILL - ERA5_BARE, `fill_source` ERA5_CLIMFILL - ERA5_DARTFILL,
`atm_clim` CLIMO_ATM - DA_COUPLED, `atm_bare_era5` DA_LND - DA_COUPLED, `lnd_clim` CLIMO_LND - DA_COUPLED,
`lnd_datm` DA_ATM - DA_COUPLED. A single-member difference is judged against the noise scale
sqrt(2) x the member spread of the 20 DA_COUPLED members of the same start.

## Notebooks

1. `01_initial_state_completeness.ipynb` - which fields each atmosphere initial file holds, their global
   means and zonal profiles (total ozone, radiative water vapour, N2O/CH4, aerosol burdens, droplet number,
   50-hPa temperature), and a check that the completed files keep the ERA5 fields unchanged.
2. `02_drift_and_bias_timeseries.ipynb` - daily regional bias against ERA5 (lead days 1-56) of every run
   with the DA_COUPLED ensemble range, and the window means of each comparison with their ratio to the noise scale.
3. `03_lead_window_difference_maps.ipynb` - maps of each comparison for days 1-3, 4-10, 11-28 and 29-56,
   stippled where |A - B| > 2 x noise.

## Results (16 Dec 2011 start)

1. **Completing the ERA5 state removes the cold start** (`tracer_fix`). The ERA5 fields are unchanged in the
   completed files. Over land the completed run is 6.1 K warmer (2-m temperature, days 1-10), with surface
   downward longwave 99 W m-2 higher in days 1-3. Total ozone is 240 DU higher globally. The bare run's land
   cold bias lasts until about day 30, and its 50-hPa temperature falls to 200 K. The completed run stays
   inside the DA_COUPLED ensemble range.
2. **The fill source matters little, except for ozone** (`fill_source`). CLIM_BASE fill instead of the DART
   mean changes TREFHT and T850 by less than the noise. Total ozone is 8-15 DU lower and T050 is about
   1.5 K warmer later. The CLIM_BASE fill also carries more CH4 (column mean 1722 vs 1048 ppbv).
3. **The production DA_LND contrast is the IC defect** (`atm_bare_era5`): -6.3 K over land in days 1-10,
   the same signal as item 1. With complete states, CLIMO_ATM - DA_COUPLED is +0.18 K globally in days 1-3
   (9 x noise) and within noise later. CLIMO_LND and DA_ATM minus DA_COUPLED are within noise in every window.
   One member and one start are too few to detect the land-initialization benefit.

## Before running

The notebooks read the 1-degree time series that `post_process/regen_f20tr_ts.sh` writes under
`<member>/post/`. For the tests (on scratch), run:

```bash
cd post_process
ROOT=/pscratch/sd/z/zhan391/e3sm_dart/simulation STREAMS="atm_6hourly atm_daily lnd_daily" \
CASES="ERA5COMPLETE_TEST_2011121600 ERA5CLIMFILL_TEST_2011121600 CLIMO_ATM_TEST_2011121600 CLIMO_LND_TEST_2011121600" \
sbatch --nodes=1 --qos=debug --time=00:30:00 submit_regen_f20tr_ts.sbatch
```

Paths are set in each notebook's user-setup cell: `INPUT_DATA_ROOT` (production), `TEST_DATA_ROOT` (tests),
`IC_ROOT` (initial-file library), `OBSERVATION_ROOT`. Outputs go to
`DIAGNOSTIC_OUTPUT_ROOT/{data,figure}/era5_ic_test/<topic>/`. Products are cached with the settings they were
built from and recomputed when those change or `FORCE_COMPUTE = True`.

## Code

* `ic_util/ic_test_runs.py`: the runs, their initial files and the comparisons. It keeps names and relative
  layouts only; the roots are notebook inputs.
* `ic_util/ic_test_eval.py`: initial-file summaries, the daily readers (model and ERA5), regions, regional
  means, lead windows and the noise scale. The model time handling and the product caching are reused,
  unchanged, from `diagnostic/util/model_daily.py`, `time_labels.py` and `hind_case_io.py`.

Add a test by adding a `Run` (and, where useful, a comparison and a group) in `ic_util/ic_test_runs.py`.
