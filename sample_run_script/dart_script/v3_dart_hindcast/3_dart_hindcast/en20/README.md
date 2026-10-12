# E3SM ensemble hindcast workflow

This directory runs one ensemble hindcast: every member is initialized once
from a DART analysis (`my_refcase`, at `my_casedate`/`my_casetod`) and then
runs a single forecast of `my_stop_n` `my_stop_option`. There is no cycling.

The numbered scripts are the entry points. Internal workers, templates and
libraries live under `workflow_lib/` and are not run directly. Runtime state
(logs, locks, completion records) lives under `runtmp/`.

## Quick start

1. Edit `create_and_setup_case.sh` (ensemble size, start date, initial
   conditions, run type, nodes, queue). For `my_ic_source=perturbed`, run
   Steps 0a and 0b first (see Initial conditions). Namelists and history output are in
   `workflow_lib/run_template/compile_and_setup_e3sm.${my_compset}.sh`.
2. Check the inputs. This changes nothing and also creates `runtmp/logs`,
   which the `#SBATCH --output` lines need:

   ```bash
   PREFLIGHT_ONLY=TRUE ./1_run_e3sm_ensemble_setup.sh
   ```

3. Run the stages in order, waiting for each to finish:

   ```bash
   sbatch 1_run_e3sm_ensemble_setup.sh     # build once, create + configure members
   sbatch 2_run_e3sm_hindcast.sh           # submit each member's hindcast
   ./3_check_hindcast.sh                   # repeat until all members are COMPLETE
   ./4_post_process_all.sh                 # submits one post job per member/product
   sbatch 5_post_6hr_to_daily.sh           # after the post jobs finished
   sbatch 6_run_compress_data.sh           # optional
   ```

Every script takes an optional member list as its first argument: `all`
(default), `3`, `1,4,7-9`, `EN05`, and for Step 2 also `failed` and
`unfinished`.

## Workflow overview

```text
1. Ensemble setup  (build e3sm.exe once, create/configure every member)
        |
2. Submit the hindcast of each member  (one model run per member)
        |
3. Check / validate the finished runs  -> hindcast_complete.ENxx
        |
        +-- 4. Post-processing  (one batch job per member and product)
        |        +-- 5. 6-hourly -> daily averages
        +-- 6. Compression (optional)
```

## Directory structure

```text
.
├── 1_run_e3sm_ensemble_setup.sh
├── 2_run_e3sm_hindcast.sh
├── 3_check_hindcast.sh
├── 4_post_process_all.sh
├── 5_post_6hr_to_daily.sh
├── 6_run_compress_data.sh
├── create_and_setup_case.sh
├── workflow_lib/
│   ├── common.sh        # shared functions (config, members, records, case setup)
│   ├── run_template/    # compile_and_setup_e3sm.<COMPSET>.sh (namelists, history)
│   ├── post/            # post-processing job scripts + post_lib.sh
│   ├── compress/        # Step 6 worker
│   ├── env/             # machine environment (nccopy etc.)
│   ├── namelists/eam/   # DART perturbation namelist (Step 0b)
│   └── maintenance/     # runtmp cleanup utility
└── runtmp/
    ├── logs/            # Slurm and per-member logs
    ├── locks/           # flock lock files
    ├── run_scripts/     # the materialized compile_and_setup_e3sm.sh that was run
    └── status/          # completion and in-progress records
```

## Initial conditions

`my_ic_source` in `create_and_setup_case.sh` selects where each member's
initial conditions come from:

- `dart_analysis` (DART hindcasts): the analysis restarts of the DA case
  `my_refcase`, under `${my_refdir}/ENxx/archive/rest/<date>`.
- `perturbed` (control and CAPT hindcasts): one unperturbed state is staged
  for every member and then perturbed with DART, before Step 1:

  ```bash
  sbatch 0a_stage_initial_conditions.sh     # copy EAM/ELM/MOSART/MPAS/cpl ICs per member
  sbatch 0b_perturb_initial_conditions.sh   # DART perturbation of the EAM state
  ```

  The EAM state is `my_refeam_in` with `PS,U,V,T,Q,CLDLIQ,CLDICE` from
  `my_refeam_ic`; the other components come from `my_refelm_in`,
  `my_refrof_in`, `my_refocn_in`, `my_refice_in`, `my_refcpl_in`. Both steps
  write into the hindcast case itself (`my_refcase` = `my_casename`,
  `my_refdir` = `my_modeldir`), and Step 1 requires the Step 0b record
  `runtmp/status/perturb_complete.<date>`. The perturbation namelist is
  `workflow_lib/namelists/eam/perturb.nml`.

## Configuration

`create_and_setup_case.sh` is the only file to edit for a new experiment.

- `my_workflow_root` is the directory holding the file, so a copied workflow
  needs no path change.
- `my_casename` is derived from the DA case (`DARTEN10_F20TR_...` ->
  `DARTEN10_15day_F20TR_...`). Set `my_casetag` to run a second hindcast from
  the same analysis.
- `my_expected_daily_streams` / `my_expected_monthly_streams` list the history
  files Step 3 requires before a member is COMPLETE. Keep them in line with the
  history settings in the run template.

## Stage guide

### Step 1 — Ensemble setup

```bash
sbatch 1_run_e3sm_ensemble_setup.sh [members]
```

Copies the run template for `my_compset` to `runtmp/run_scripts/` and runs it
for the base case (compiles `e3sm.exe`), then for each member, reusing that
executable. Each member is then pointed at its DA restarts and verified.
The restarts are read from `REF/ENxx/archive/rest/<date>`, or from the older
shared layout. Settings: `REBUILD`, `RECREATE_EXISTING`,
`ADOPT_EXISTING_CASE`, `PREFLIGHT_ONLY`.

Writes `build_complete`, `setup_complete.ENxx` and, once every member is
ready, `setup_complete`. Members already set up are skipped.

### Step 2 — Submit the hindcast

```bash
sbatch 2_run_e3sm_hindcast.sh [members]
sbatch 2_run_e3sm_hindcast.sh failed
```

Requires `setup_complete.ENxx`. For each member it sets the hybrid start
(`RUN_REFCASE`/`RUN_REFDIR`, `ncdata`, `finidat`, `finidat_rtm`, MPAS restart
streams, `rpointer.*`). It runs `preview_namelists`, checks that the generated
`atm_in`, `lnd_in`, `mosart_in`, streams and pointers reference the intended
files, and then runs `case.submit`. A member that fails verification is not
submitted.

Writes `hindcast_submitted.ENxx` and removes that member's downstream records.
Queued, running, finished and complete members are skipped
(`RESUBMIT_COMPLETED=TRUE` to rerun). `DRY_RUN=TRUE` shows what would be
submitted.

### Step 3 — Check and validate

```bash
./3_check_hindcast.sh [members]
```

Prints each member's state:

| State | Meaning |
|---|---|
| NO_CASE | not set up |
| READY | set up, never submitted |
| PENDING / RUNNING | queued / running |
| FINISHED | run ended, not yet validated |
| COMPLETE | outputs validated (`hindcast_complete.ENxx`) |
| INVALID | run ended but outputs incomplete (reason shown) |
| FAILED | run failed or was cancelled |
| FAILED_ARCHIVE | model finished, short-term archiving failed |

Validation requires the final restart (`rest/<end>`) plus the expected daily
and monthly history files. For `FAILED_ARCHIVE`, run `./case.st_archive` in
the member's `case_scripts` and check again. Resubmit failures with
`sbatch 2_run_e3sm_hindcast.sh failed`.

### Step 4 — Post-processing

```bash
./4_post_process_all.sh [members] [products...]
sbatch 5_post_6hr_to_daily.sh [members]
```

This submits one batch job per member and
per post script in `workflow_lib/post/`. The variable lists are edited in
those scripts.

| Product | Script | Input |
|---|---|---|
| eam_6hourly | `post_process_eam_6hourly.sh` | `eam.h2` |
| eam_daily | `post_process_eam_daily.sh` | `eam.h1` |
| eam_clim | `post_process_eam_monthly.sh` | `eam.h0` |
| elm_daily | `post_process_elm_daily.sh` | `elm.h1` |
| elm_clim | `post_process_elm_monthly.sh` | `elm.h0` |

Before submitting anything, the script checks each member's
`hindcast_complete` record and every history file the selected products need
for `POST_START`..`POST_END` (defaults from the config). Each job writes
`postprocess_complete.<product>.ENxx.<start>_<end>` after checking its
outputs. Verified products and jobs already in the queue are not submitted
again (`OVERWRITE_EXISTING=TRUE` to redo).

`5_post_6hr_to_daily.sh` averages the regridded 6-hourly series to
daily once `eam_6hourly` is complete. Daily files already written from
`eam.h1` are kept.

Outputs: `ENxx/archive/post/{atm,lnd}/180x360_aave/{ts/<freq>,clim}/`.

### Step 6 — Compression (optional)

```bash
sbatch 6_run_compress_data.sh [members]
```

Requires `hindcast_complete.ENxx`. Compresses to netCDF-4 classic with
deflate level 5. Each file is written to a temporary copy, checked with
`ncdump`, then renamed into place. Already-compressed files are skipped, and
leftovers of an interrupted run are cleaned up. Settings: `COMPRESS_HISTORY`,
`COMPRESS_RESTARTS`, `COMPRESS_MAX_PARALLEL`. Uses
`workflow_lib/env/env_${my_machine}_specific.sh`.

## Completion records, locks and reruns

- `runtmp/status/*_complete*`: validated work. Each record carries the
  workflow root, case, ensemble size and start time; a record that does not
  match the current configuration is ignored.
- `runtmp/status/*_in_progress*`: an attempt started and has not finished.
- `runtmp/locks/*.lock`: `flock` locks per stage and per member, which stop
  two scripts working on the same thing.

Resubmitting a stage only does what is left. Do not create or delete records
by hand to skip a stage; fix the cause and rerun that stage.

## New DA experiment

```bash
rsync -a --exclude runtmp en10/ en_new/
cd en_new && vi create_and_setup_case.sh    # my_refcase, my_ensnum, my_casedate, ...
PREFLIGHT_ONLY=TRUE ./1_run_e3sm_ensemble_setup.sh
```

The case directory records which workflow directory created it
(`.hindcast_owner`). A copy that still points at another workflow's case is
refused instead of modifying it.

## Logs and troubleshooting

```bash
squeue -u "$USER"
./3_check_hindcast.sh
ls -t runtmp/logs | head        # setup_ENxx.log, run_ENxx.log, post_<product>_ENxx.*, ...
ls runtmp/status
```

Each driver ends with a one-line `INFO`/`ERROR` summary. A failing member is
listed together with the command to retry it.

## Runtime cleanup

```bash
workflow_lib/maintenance/cleanup_runtmp.sh --dry-run
workflow_lib/maintenance/cleanup_runtmp.sh --apply
```

Removes logs older than 30 days. Records, run scripts and locks are kept. It
refuses to run while any workflow lock is held.

## Validation after editing scripts

```bash
find . -name '*.sh' -not -path './runtmp/*' -print0 | xargs -0 -n1 bash -n
```
