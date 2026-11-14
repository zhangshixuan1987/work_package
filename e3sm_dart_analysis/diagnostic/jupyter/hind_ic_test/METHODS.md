# ERA5 initial-condition tests: data and methods

Evaluates the test hindcasts that address the issues found in the ERA5-initialized F20TR hindcasts.
Every test is a single member started on 16 Dec 2011 00Z and run for 2 months. Each is compared with the
production members (`INPUT_DATA_ROOT`, eam_hindcast) that share its other initial states.

## Focus: ERA5 + CLIM_BASE aerosol and ozone initialization vs ERA5 initialization (HYB_BASE)

The notebooks verify, by default, one pair of runs from the same member EN01: HYB_BASE EN01 (`ERA5_BARE`, the
ERA5 file as used in production) and `ERA5_CLIMFILL` EN01. Both start from the same ERA5 state with the same
DART perturbation EN01, over the same DATM-forced land EN01. ERA5_CLIMFILL keeps the ERA5 fields unchanged and
takes every field the ERA5 file lacks (aerosols, ozone and the other radiative gases, hydrometeor numbers)
from the CLIM_BASE (FREE_CLIM) state of the same date. Their difference `clim_fill` = ERA5_CLIMFILL - ERA5_BARE
is therefore the effect of that initialization. It is judged against sqrt(2) x the member spread of the 20
HYB_BASE members of the same start (`ENSEMBLE_RUN = "ERA5_BARE"`), and both runs are verified against ERA5.
Aerosol and ozone effects are traced through the ozone columns (TOZ, SCO, TCO), stratospheric temperature
(T100, T050, T010), the clear-sky shortwave (FSDSC, FSNTOAC) and longwave (FLDSC, FLUTC) and the surface
temperature. The other runs and comparisons below stay in `configs/hind_ic_test_runs.py`; add them to
`RUN_KEYS`, `PLOT_GROUPS` and `COMPARISON_KEYS` to include them.

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

## Runs (`configs/hind_ic_test_runs.py`)

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

Comparisons (A - B): `clim_fill` ERA5_CLIMFILL - ERA5_BARE (the focus), `tracer_fix` ERA5_DARTFILL - ERA5_BARE, `fill_source` ERA5_CLIMFILL - ERA5_DARTFILL,
`atm_clim` CLIMO_ATM - DA_COUPLED, `atm_bare_era5` DA_LND - DA_COUPLED, `lnd_clim` CLIMO_LND - DA_COUPLED,
`lnd_datm` DA_ATM - DA_COUPLED. A single-member difference is judged against the noise scale
sqrt(2) x the member spread of the 20 members of `ENSEMBLE_RUN` of the same start (HYB_BASE for the focus;
DA_COUPLED for the component swaps).

## Notebooks

1. `01_compute_and_plot_initial_state_completeness.ipynb` - which fields each atmosphere initial file holds, their global
   means and zonal profiles (total ozone, radiative water vapour, N2O/CH4, aerosol burdens, droplet number,
   50-hPa temperature), and a check that the completed files keep the ERA5 fields unchanged.
2. `02_compute_and_plot_drift_and_bias.ipynb` - daily regional bias against ERA5 (lead days 1-56) of each run
   with the range of the `ENSEMBLE_RUN` members, and the window means of each comparison with their ratio to the
   noise scale. Section 7 scores a test member against the whole ensemble of another configuration
   (`ENSEMBLE_COMPARISONS`; by default ERA5_CLIMFILL against the 20 DA_COUPLED members): ensemble-mean
   difference, z = difference / (spread x sqrt(1 + 1/N)), and the run's rank among the members.
3. `03_compute_and_plot_lead_window_maps.ipynb` - maps of each comparison for days 1-3, 4-10, 11-28 and 29-56,
   stippled where |A - B| > 2 x noise.

## Results (16 Dec 2011 start)

**ERA5 + CLIM_BASE fill vs HYB_BASE, member EN01** (`clim_fill`; noise from the 20 HYB_BASE members):

1. The filled run starts with the ozone, radiative water vapour and aerosols the HYB_BASE member lacks:
   total ozone +231 DU globally in days 1-3 (stratospheric +201, tropospheric +28), still +106 DU in days 29-56
   as HYB_BASE recovers only slowly. Clear-sky surface downward longwave is +206 W m-2 and clear-sky surface
   downward shortwave -58 W m-2 in days 1-3 (absorption by water vapour, ozone and aerosols).
2. This removes the cold start. Land 2-m temperature bias against ERA5 goes from -6.0 K (HYB_BASE EN01, and the
   HYB_BASE ensemble mean) to +0.0 K in days 1-3 and +0.2 K in days 4-10; the difference is 10-85 x noise
   through day 28, and the filled run is warmer than all 20 HYB_BASE members. T850 bias goes from -2.4 to -5.7 K
   to -0.7 to -1.6 K.
3. After about day 28 the near-surface fields converge (within 1-2 x noise), but HYB_BASE overshoots warm
   (land TREFHT bias +1.6 K, Z500 +37 m vs +1.0 K and +7 m). The stratosphere stays apart: 50 hPa is 9 K
   warmer in the filled run in days 29-56, 10 hPa 2.6 K colder after day 10.

**ERA5 + CLIM_BASE fill vs the DA_COUPLED ensemble** (`climfill_vs_dacoupled`; atmosphere and land differ):

4. Near the surface the filled run behaves like a DA_COUPLED member after day 3 (TREFHT, T850, FLDS, PRECT
   within the member spread). In days 1-3 it is closer to ERA5 (land TREFHT bias +0.02 vs +0.28 K), as it
   starts from ERA5.
5. Ozone and the stratosphere differ throughout: total ozone bias -14 DU globally vs -5 DU for DA_COUPLED
   (NH land -23 vs +5 DU; better only in the tropics), the CLIM_BASE ozone being lower than DART's. 50 hPa is
   1-2.5 K warmer after day 10 and Z500 about 17 m higher in every window (rank 19-20 of 20). These come from
   the fill source, not from ERA5 vs DART, and will carry into the rerun HYB_BASE/DA_LND vs DA_COUPLED contrasts.

**Other tests** (not run by default; noise from DA_COUPLED):

6. The DART-mean fill gives the same near-surface result as the CLIM_BASE fill (difference within noise);
   ozone is 8-15 DU higher and 50 hPa about 1.5 K colder. Its CH4 column mean is 1048 ppbv vs 1722 for CLIM_BASE.
7. The production DA_LND - DA_COUPLED contrast is the same defect (-6.3 K over land in days 1-10), so
   HYB_BASE and DA_LND must both be rerun from completed states. With complete states, CLIMO_ATM - DA_COUPLED
   is +0.18 K globally in days 1-3 (9 x noise) and within noise later; CLIMO_LND and DA_ATM minus DA_COUPLED
   are within noise in every window. One member and one start are too few to detect the land-initialization benefit.

## Before running

The notebooks read the 1-degree time series that `post_process/submit_process_e3sm_post.sbatch` writes under
`<member>/post/`. For the tests (on scratch), run:

```bash
cd ../../../post_process   # e3sm_dart_analysis/post_process
ROOT=/pscratch/sd/z/zhan391/e3sm_dart/simulation \
KINDS="ts_atm_6hourly_180x360_aave ts_atm_daily_180x360_aave ts_lnd_daily_180x360_aave" \
CASES="ERA5COMPLETE_TEST_2011121600 ERA5CLIMFILL_TEST_2011121600 CLIMO_ATM_TEST_2011121600 CLIMO_LND_TEST_2011121600" \
sbatch --nodes=1 --qos=debug --time=00:30:00 submit_process_e3sm_post.sbatch
```

Paths are set in each notebook's user-setup cell: `INPUT_DATA_ROOT` (production), `TEST_DATA_ROOT` (tests),
`IC_ROOT` (initial-file library), `OBSERVATION_ROOT`. Outputs go to
`DIAGNOSTIC_OUTPUT_ROOT/{data,figure}/hind_ic_test/<topic>/`. Products are cached with the settings they were
built from and recomputed when those change or `FORCE_COMPUTE = True`.

## Code

* `configs/hind_ic_test_runs.py`: the runs, their initial files and the comparisons. It keeps names and relative
  layouts only; the roots are notebook inputs.
* `util/hind_ic_test_eval.py`: initial-file summaries, the daily readers (model and ERA5), regions, regional
  means, lead windows and the noise scale. The model time handling and the product caching are reused,
  unchanged, from `diagnostic/util/model_daily.py`, `time_labels.py` and `hind_case_io.py`.

Add a test by adding a `Run` (and, where useful, a comparison and a group) in `configs/hind_ic_test_runs.py`.
