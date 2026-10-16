# v3_dart_cda run scripts

Snapshot of the E3SMv3–DART workflow scripts from
`/qfs/people/zhan391/e3sm_dart_work/run_script/v3_dart_cda`, organized by
experiment group. Each case folder holds its numbered step scripts,
`create_and_setup_case.sh` and `workflow_lib/`; runtime state (`runtmp/`), job
logs and diagnostic outputs are not included. Paths inside the configs still
point to the original run directories.

| Folder | Case | Period | Description |
|---|---|---|---|
| `1_ne30pg2_ctrl_atm_en10`  | CTRLEN10   | 2011-12 | Atmosphere control, DART cycling with all obs evaluate-only |
| `1_ne30pg2_ctrl1_atm_en10` | CTRFEN10   | 2011-12 | Atmosphere control, continuous free run (no DART cycling) |
| `1_ne30pg2_dart_atm_en10`  | DARTEN10   | 2011-12 | EAM DA, 10 members |
| `1_ne30pg2_dart_atm_en20`  | DARTEN20   | 2011-12 | EAM DA, 20 members |
| `1_ne30pg2_dart_atm_en40`  | DARTEN40   | 2011-12 | EAM DA, 40 members |
| `2_ne30pg2_ctrl_atm_en10`  | CTRLEN10S1 | 2012-05 | Atmosphere control, DART cycling evaluate-only |
| `2_ne30pg2_ctrl1_atm_en10` | CTRFEN10S1 | 2012-05 | Atmosphere control, continuous free run |
| `2_ne30pg2_dart_atm_en40`  | DARTEN40S1 | 2012-05 | EAM DA, 40 members |
| `3_ne30pg2_ctrl_cpl_en10`  | CTRLEN10 (WCYCL20TR) | 2011-12 | Coupled control (not run) |
| `3_ne30pg2_dart_cpl_en40`  | DARTEN40 (coupled)   | 2011-12 | Coupled EAM DA, 40 members |
| `4_ne30pg2_dart_scp_en10`  | DARTSCPEN10          | 2011-12-01 | Strongly coupled DA test case |
| `e3sm_dart_compare`        | —          | 2011-12 | Observation-space comparison across experiments |

Earlier hindcast experiments (CAPT, DART, control) are in `../v3_dart_hindcast/`.

## Workflow

DA and evaluate-only control cases run Steps 1-8 in order (`1_run_e3sm_ensemble_setup.sh`
through `8_run_post_init.sh`); see each case's `README.md`.

The `ctrl1` free-run cases run Steps 1-3 and then, on a login node,
`./4_run_e3sm_continuous.sh` (submit / `status` / `finalize`): each member runs
as its own CIME job chain, so no DART cycling or `eam.i` output is produced.

## History output and post-processing (Step 7)

Every DA cycle is a separate 6-h model run started from the analysis, so each
history file holds an nstep-0 record (the analysis for instantaneous fields; a
single-timestep value for mean/max/min fields) and, on tapes written every 6 h,
a record spanning the cycle (6-h forecast means/extremes; background state).

- `7_run_post_hist.sh` is the same in every case. Cycled cases run
  `workflow_lib/post/cycled_split.py`; free runs run `free_run_hfreq.py`. Both
  write, under `ENxx/archive/post/<atm|lnd>/180x360_aave/`, the layout used by
  the hindcasts: `ts/6hourly`, `ts/6hourly_analysis` (cycled only),
  `ts/daily_mean`, `ts/daily_max`, `ts/daily_min`, `ts/daily_00z` and
  `monthly/`. Records are chosen by time bounds and `cell_methods`; daily values
  cover calendar days 00Z-24Z.
- Cases set up with `nhtfrq=-24` (or `mfilt=1`) kept no 6-h statistics on those
  tapes; they are listed in `ANALYSIS_ONLY_TAPES` and give analysis states only.
  In `3_ne30pg2_dart_cpl_en40` the EAM h2 6-h means were overwritten
  (`mfilt=1`), so it has no 6-h flux or precipitation means.
- New cases: the run template writes every EAM/ELM/MOSART tape at `-6` with
  `mfilt=4`, and h1 holds per-member U/V/T/Q/Z on 16 pressure levels (analysis
  and background). Step 4 refuses to cycle unless every tape has
  `nhtfrq=-<cycle hours>` and `mfilt>=2`
  (`workflow_lib/cycle/check_history_namelists.sh`). ELM-DART finds its prior
  record in either history layout, so its streams can keep `mfilt>=2` too.
