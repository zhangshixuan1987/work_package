# E3SM–DART Atmospheric Ensemble Workflow

This directory contains the restart-safe, Slurm-driven workflow for the 10-member AMIP control ensemble (`CTRLEN10`). The numbered scripts are the user-facing entry points. This source tree was refactored after 200 legacy six-hour cycles, with archived output through 2012-01-20-00000; existing model output is not moved automatically. Internal workers and templates live under `workflow_lib/` and should not normally be executed directly.

## Quick start

1. Edit `create_and_setup_case.sh` and verify the machine, project, ensemble, model, reference restart, DART, and date settings.
   Create the Slurm log directory before the first submission if it does not already exist: `mkdir -p runtmp/logs`.
2. Confirm the required E3SM/DART code, input data, observations, map files, and reference restarts exist.
3. Submit Steps 1–4 in order, waiting for each stage to complete successfully:

   ```bash
   sbatch 1_run_e3sm_ensemble_setup.sh
   sbatch 2_run_dart_e3sm_icbc.sh
   sbatch 3_run_dart_eam_pertrub.sh
   sbatch 4_run_dart_e3sm_cycleda.sh
   ```

4. Run the optional compression, diagnostics, history processing, and restart processing stages as needed:

   ```bash
   sbatch 5_run_dart_compress.sh START_TIMESTAMP END_TIMESTAMP
   sbatch 6_run_dart_diag.sh
   sbatch 7_run_post_hist.sh
   sbatch 8_run_post_init.sh
   ```

Use timestamps in `YYYY-MM-DD-SSSSS` form, where `SSSSS` is seconds since midnight. For example, midnight is `00000` and 18 UTC is `64800`.

Do not bypass a failed stage by manually creating completion records. Downstream stages use those records as part of their safety checks.

## Workflow overview

```text
1. Ensemble setup
        ↓
2. Initial and boundary conditions
        ↓
3. Initial ensemble perturbation
        ↓
4. Atmospheric forecast–assimilation cycling
        ├── 5. Compression (optional, independent range driver)
        ├── 6. DART diagnostics
        ├── 7. EAM/ELM history post-processing
        └── 8. EAM/ELM restart-variable extraction
```

Steps 1–4 form the core cycling workflow. Steps 5–8 are post-cycle utilities and may be run when their required upstream cycle records and data are available.

## Directory structure

```text
.
├── 1_run_e3sm_ensemble_setup.sh
├── 2_run_dart_e3sm_icbc.sh
├── 3_run_dart_eam_pertrub.sh
├── 4_run_dart_e3sm_cycleda.sh
├── 5_run_dart_compress.sh
├── 6_run_dart_diag.sh
├── 7_run_post_hist.sh
├── 8_run_post_init.sh
├── create_and_setup_case.sh
├── workflow_lib/
│   ├── compress/       # Step 5 worker
│   ├── cycle/          # Step 4 cycle, assimilation, and handoff logic
│   ├── diagnostics/    # Step 6 workers
│   ├── namelists/      # Authoritative DART runtime namelist templates
│   ├── maintenance/    # Safe runtime cleanup utility
│   ├── post/           # Step 7 workers
│   └── run_template/   # Maintained E3SM run-script templates
├── runtmp/
│   ├── handoff/        # Transaction recovery snapshots
│   ├── locks/          # flock lock files
│   ├── logs/           # Slurm and worker logs
│   ├── run_scripts/    # Generated executable run scripts
│   └── status/         # Completion and in-progress records
└── deprecated/         # Historical scripts; not part of the active workflow
```

Only the numbered scripts should normally be submitted with `sbatch`. Files in `workflow_lib/` are workers, sourced libraries, templates, or maintenance utilities.

## Shared configuration

`create_and_setup_case.sh` is the authoritative configuration. Review it before starting a new experiment.

Important groups include:

- Runtime environments: `my_conda_setup_file` and `my_analysis_conda_env` for Steps 2–3, plus `my_dart_env_file` for DART-dependent stages.
- Slurm resources: `my_task_per_node`, `my_job_nnodes`, `my_project`, `my_jobqueue`, and `my_walltime`.
- Ensemble configuration: `my_ensnum`, `my_nodes_per_member`, setup concurrency, and forecast retry settings.
- Model configuration: `my_e3sm_code`, `my_runtype`, `my_compset`, `my_resolution`, `my_runpath`, and `my_casename`.
- Initial state: `my_casedate`, `my_casetod`, `my_refcase`, `my_refdate`, `my_reftod`, `my_refdir`, and component restart paths.
- Timeline and DART configuration: `my_e3sm_cycle_hours`, the shared E3SM start/end time, component-specific `my_eam_dart_cycle_hours` and `my_elm_dart_cycle_hours`, DART code and run directories, observation paths, and diagnostic ranges.
- Cycling behavior: `my_cycles_per_job`, minimum cycle runtime, shutdown margin, and handoff concurrency.

`my_conda_setup_file` and `my_analysis_conda_env` explicitly select the analysis environment that supplies NCO and related tools for initial-condition generation and perturbation.

`my_dart_env_file` names the workflow-owned, machine-specific environment used by cycling, compression, and diagnostics. It resolves from `my_machine` to `workflow_lib/env/env_${my_machine}_specific.sh` (for example, `env_compy_specific.sh`) rather than to a generated file beneath a DART `work/` directory. The numbered drivers validate it before starting substantive work.

Runtime DART namelists are also workflow-owned. The explicit
`my_eam_filter_nml`, `my_eam_perturb_nml`, `my_eam_diag_nml`, and
`my_elm_filter_nml` settings select the templates under
`workflow_lib/namelists/`. DART model `work/input.nml` files remain reserved for
build-time utilities such as `preprocess`.

Exceptional assimilation cycles can override localization cutoff and inflation damping in the `my_eam_dart_cycle_overrides` table. Unlisted cycles continue to use the ensemble-size defaults.

The configuration defines these canonical workflow paths from the explicit `my_workflow_root` setting:

```bash
my_workflow_root
my_workflow_lib
my_runtime_dir
my_log_dir
my_status_dir
my_lock_dir
my_handoff_dir
my_run_script_dir
```

This allows the scripts to resolve their workers consistently from Slurm, an interactive shell, or another directory.

## Stage guide

### Step 1 — Ensemble setup

```bash
sbatch 1_run_e3sm_ensemble_setup.sh
```

Step 1 selects the model template matching `my_compset` from `workflow_lib/run_template/`, materializes an executable script under `runtmp/run_scripts/`, builds/configures the base case, and creates the ensemble cases.

It writes:

```text
runtmp/status/setup_complete
```

only after all requested ensemble case directories pass validation. A rerun removes the old completion record before modifying setup state.

### Step 2 — Initial and boundary conditions

```bash
sbatch 2_run_dart_e3sm_icbc.sh
```

Step 2 requires a matching Step 1 completion record. It prepares restart inputs for all ensemble members and validates the required components before declaring completion:

- EAM
- ELM
- MOSART
- MPAS-I
- Coupler
- MPAS-O for Full-CPL runs

AMIP does not require MPAS-O, but this workflow still requires MPAS-I and coupler files.

The completion record is named:

```text
runtmp/status/icbc_complete.<valid-time>
```

### Step 3 — Initial perturbation

```bash
sbatch 3_run_dart_eam_pertrub.sh
```

Step 3 requires the matching Step 2 record. It calculates DART tasks from the actual Slurm node allocation and `my_task_per_node`, perturbs the initial ensemble, validates the result, and writes:

```text
runtmp/status/perturb_complete.<valid-time>
```

An archive-level `.dart_perturb_in_progress` marker protects partially perturbed ensembles.

### Step 4 — Atmospheric cycling DA

```bash
sbatch 4_run_dart_e3sm_cycleda.sh
```

Step 4 requires the initial perturbation record and verifies that the Slurm allocation matches `my_job_nnodes`. Before running a cycle, it checks whether the configured DA end time has already been reached.

Each cycle performs ensemble forecasts, validates restart products and observations, runs DART assimilation, transactionally commits the next-cycle state, updates the cycle counter, and writes:

```text
runtmp/status/cycle_complete.<valid-time>
```

The driver can run multiple cycles per allocation using `my_cycles_per_job`. It stops conservatively when insufficient wall time remains and can recover or submit a continuation job without duplicating active continuations.

For an explicit migration recovery where all case XML clocks are uniformly
older than a known completed-cycle archive, set
`my_allow_forward_xml_recovery="TRUE"`. Step 4 then validates every member
restart and archived EAM DART evidence before reconstructing auditable recovery
metadata. The override is reset to `FALSE` after a cycle commits successfully.

Handoff snapshots use visible names such as:

```text
runtmp/handoff/cycle_109.job_774000/
```

The newest successful snapshot is retained. Older snapshots are removed only when their cycle and job IDs match durable completion records. Failed or unverified snapshots remain for recovery.

### Step 5 — Compression

Step 5 is optional and is not automatically submitted by Step 4.

```bash
sbatch 5_run_dart_compress.sh 2011-12-20-00000 2011-12-26-00000
```

The range must be aligned with `COMPRESS_INTERVAL_HOURS`, which defaults to 6. The worker validates upstream completion state, prevents conflicts with active cycling, checks available space, validates NetCDF files, and replaces compressed files atomically.

Compression behavior can be controlled through environment variables such as:

```bash
export COMPRESS_HISTORY=TRUE
export COMPRESS_RESTARTS=FALSE
export COMPRESS_MAX_PARALLEL=4
export COMPRESS_RESTART_MAX_PARALLEL=1
export COMPRESS_SPACE_MARGIN_MB=1024
```

### Step 6 — DART diagnostics

Edit the user settings at the top of `6_run_dart_diag.sh`:

```bash
DIAG_MODE="obs"   # seq, obs, common, closest, or all
DIAG_START=""     # empty uses my_eam_dart_diag_start
DIAG_END=""       # empty uses my_eam_dart_diag_end
```

The modes are independent:

- `seq` preprocesses the cycle `obs_seq_final` files.
- `obs` runs the standard observation-space diagnostics.
- `common` runs diagnostics for observations common across cycles and writes to `dart_diagnostics/cam6_common/`.
- `closest` runs `closest_member_tool` and writes to `dart_diagnostics/closest_member/`.
- `all` runs `seq`, `obs`, and `common`; it also runs `closest` when `my_eam_dart_diag_run_closest_member` is `.true.`.

Enable closest-member processing in `create_and_setup_case.sh`:

```bash
export my_eam_dart_diag_run_closest_member=".true."
```

Closest mode validates every ensemble-member restart over the requested range. It creates fresh per-quantity summaries in temporary files and atomically replaces existing `*.closest_member.txt` summaries after successful completion, so rerunning a range does not append duplicate records.

Then submit:

```bash
sbatch 6_run_dart_diag.sh
```

Step 6 validates the entire requested cycle range before starting. Observation modes require completed observation sequences; closest mode requires all ensemble restart files. Mode-specific workers live under `workflow_lib/diagnostics/` and cannot be executed directly. Each mode has a separate completion record containing its mode, range, case, and ensemble size. A matching completion record skips only that mode, allowing failed or newly enabled modes to be run independently.

### Step 7 — History post-processing

Edit the settings at the top of `7_run_post_hist.sh`:

```bash
POST_START="2011-12-01-00000"
POST_END="2011-12-31-00000"      # last day processed (inclusive)
MAX_CONCURRENT_WORKERS=4
OVERWRITE_EXISTING="FALSE"
RUN_EAM="TRUE"
RUN_ELM="TRUE"
ANALYSIS_ONLY_TAPES="eam.h1 elm.h1"   # tapes this archive holds without 6-h records
```

Submit with:

```bash
sbatch 7_run_post_hist.sh
```

Each DA cycle is a separate 6-h run started from the analysis, so every
history file holds a nstep-0 record (zero-width bounds) and, on tapes written
every 6 h, a record spanning the cycle. Point fields in the nstep-0 record are
the analysis state; mean/max/min fields there are a single-timestep value and
are never used. A tape without 6-h records stops Step 7 unless it is listed in
`ANALYSIS_ONLY_TAPES` (runs set up with `nhtfrq=-24`, or with `mfilt=1`, which
lets the next cycle overwrite them). Step 4 refuses to start cycling unless every
tape has `nhtfrq=-<cycle hours>` and `mfilt>=2`
(`workflow_lib/cycle/check_history_namelists.sh`). Tapes with `nhtfrq=-24` hold the nstep-0 record only, and the
monthly tape covers only the last cycle, so it is not processed.
`workflow_lib/post/cycled_split.py` (products `eam_cycled`, `elm_cycled`)
selects records by their time bounds and `cell_methods` and writes, under
`ENxx/archive/post/<atm|lnd>/180x360_aave/`:

| Output | Contents |
|---|---|
| `ts/6hourly` | forecast: 6-h means/max/min labelled at the window middle with bounds; point fields (background) at the window end |
| `ts/6hourly_analysis` | analysis state at every cycle start (point fields, plus state variables labelled `time: mean`, e.g. `U850`, `TSOI`) |
| `ts/daily_mean` | calendar day 00Z-24Z: mean of the four 6-h forecast means; fields without one use the mean of the 00/06/12/18Z analysis states (stated in `cell_methods`) |
| `ts/daily_max`, `ts/daily_min` | from 6-h maxima/minima; only for runs with `nhtfrq=-6` on that tape |
| `ts/daily_00z` | analysis state at 00Z |
| `monthly/VAR.ENxx.YYYY-MM.nc` | mean of `ts/daily_mean` over the days of the month in range |

The first run on a member moves the outputs of the earlier version
(`ts/6hourly`, `ts/6hourly_da`, `ts/daily`, `monthly`, `clim`) to
`180x360_aave/legacy/` and leaves a `LAYOUT` marker. Each member/product has
its own lock and completion record; the worker lists every file it wrote in
`runtmp/status/postprocess_manifest.*`, and each must be a readable NetCDF
file before completion is recorded.

Free-run cases (`my_free_layout` set; `eam_hfreq`, `elm_hfreq`,
`free_run_hfreq.py`) write the same layout without `ts/6hourly_analysis` (a
free run has no analysis; `ts/daily_00z` is the model state at 00Z). Their
earlier `ts/daily`, `monthly` and `clim` are moved to `legacy/` on the first run.

### Step 8 — Restart-variable extraction

Edit the settings at the top of `8_run_post_init.sh`:

```bash
INTERP_START="2011-12-15-00000"
INTERP_END="2011-12-28-00000"
MAX_CONCURRENT_WORKERS=4
OVERWRITE_EXISTING="FALSE"
PROCESS_EAM="TRUE"
PROCESS_ELM="TRUE"
REQUIRE_UPSTREAM_MARKERS="TRUE"
```

Then submit:

```bash
sbatch 8_run_post_init.sh
```

Step 8 operates on midnight restart files. It preflights every requested date and member before writing output, extracts the configured EAM and ELM variables, validates results, and moves them into place atomically.

Keep `REQUIRE_UPSTREAM_MARKERS=TRUE` for normal operation. Set it to `FALSE` only when deliberately processing a legacy restart archive whose files have been independently validated but whose old completion markers were not retained.

## Dependencies and submission practice

The scripts enforce data dependencies through completion records, so submitting a downstream stage prematurely fails safely. The clearest operating practice is still to submit stages after verifying the previous stage completed:

```bash
squeue -u "$USER"
tail -f runtmp/logs/<relevant-log>
```

A Slurm `afterok` dependency can be used for Steps 1–3, but it does not replace completion-record validation:

```bash
job1=$(sbatch --parsable 1_run_e3sm_ensemble_setup.sh)
job2=$(sbatch --parsable --dependency="afterok:${job1}" 2_run_dart_e3sm_icbc.sh)
job3=$(sbatch --parsable --dependency="afterok:${job2}" 3_run_dart_eam_pertrub.sh)
sbatch --dependency="afterok:${job3}" 4_run_dart_e3sm_cycleda.sh
```

Step 4 manages its own continuation behavior after it starts cycling.

## Completion records, locks, and reruns

The workflow uses three kinds of runtime control data:

- `runtmp/status/*_in_progress*`: a stage or member/product attempt has started.
- `runtmp/status/*_complete*`: validated work completed successfully.
- `runtmp/locks/*.lock`: `flock` coordination prevents concurrent conflicting operations.

Most stages remove or invalidate their old completion record when a rerun begins. If the rerun fails, downstream stages cannot mistake the old result for a new success.

Do not delete status records merely to force a downstream stage to run. Correct the failed upstream condition and rerun the owning stage.

## Logs and troubleshooting

Slurm and worker logs are written under:

```text
runtmp/logs/
```

Useful checks include:

```bash
squeue -u "$USER"
find runtmp/status -maxdepth 1 -type f -print | sort
find runtmp/handoff -mindepth 1 -maxdepth 1 -type d -print
find runtmp/logs -type f -printf '%TY-%Tm-%Td %TT %p\n' | sort | tail
```

When a stage fails:

1. Read its Slurm log and any member-specific failure/retry logs.
2. Check the corresponding in-progress record.
3. Confirm required upstream completion records still match the case, valid time, and ensemble size.
4. Correct missing or invalid inputs.
5. Rerun the same numbered stage. Do not run its internal worker directly.

## Safe runtime cleanup

Preview cleanup first:

```bash
workflow_lib/maintenance/cleanup_runtmp.sh --dry-run
```

Apply exactly the reported cleanup:

```bash
workflow_lib/maintenance/cleanup_runtmp.sh --apply
```

Default cleanup behavior:

- Removes regular log files older than 30 days.
- Preserves all completion records.
- Preserves all in-progress records.
- Preserves generated run scripts.
- Preserves lock files.
- Preserves failed handoffs.
- Refuses to run when another runtime cleanup or a held workflow lock is detected.
- Validates that cleanup targets remain beneath `runtmp/`.

Optional stale-progress and failed-handoff cleanup is disabled by default. Review the settings at the top of `workflow_lib/maintenance/cleanup_runtmp.sh` before enabling either option.

Never remove `runtmp/` wholesale while jobs are active or while the workflow may need to resume.

## Validation after editing scripts

After changing shell scripts, run:

```bash
find . -type f -name '*.sh' ! -path './deprecated/*' ! -path './runtmp/logs/*' -print0 \
  | xargs -0 bash -n
```

This catches shell syntax errors but does not replace runtime preflight or a controlled test on the target machine.

## Deprecated scripts

Files under `deprecated/` are retained only for historical reference. They are not part of the active workflow and may lack current safety checks. Do not submit them as replacements for numbered stages.

### Component assimilation scheduling in Step 4

The E3SM forecast advances on the shared `my_e3sm_cycle_hours` timeline. EAM
and ELM assimilation have independent enable switches, cadences, and end times.
At each forecast valid time, a component is classified as `on`, `not_due`,
`ended`, or `off`. Cycle handoff and the shared cycle counter advance only after
every component due at that time completes successfully and validates all
ensemble members.

When EAM and ELM are both due, their execution mode is derived from the strongly
coupled settings:

- Direct mode is used unless both `strongly_coupled_on=on` and
  `lnd_da_use_sequential_prior_post=.true.`. EAM and ELM run concurrently and
  split the 80-node allocation equally.
- Sequential mode is used when both settings above are enabled. EAM runs first
  on all 80 nodes and produces the sequential prior; after EAM succeeds and
  the dependent inputs validate, ELM runs on all 80 nodes.

If only one component is due, it receives all 80 nodes. If neither component
is due, Step 4 performs a forecast-only cycle. A failed component assimilation
leaves an in-progress marker, and the retry path rebuilds every forecast member
before assimilation is attempted again.

In direct mode, land observations are read from `my_elm_dart_obsdir` using
`YYYYMM_6H/obs_seq.YYYY-MM-DD-SSSSS`. In sequential mode, ELM consumes the EAM
sequential-prior output. This experiment uses separate working copies of the
gridded ELM `h1` stream for history and vector-history input.

Component assimilation is controlled in `create_and_setup_case.sh`:

```bash
export my_eam_dart_da="on"   # on or off
export my_elm_dart_da="on"   # on or off
```

The enable switches do not determine cadence: an enabled component runs only
when its component-specific interval is due and its configured end time has not
been exceeded. In the current configuration EAM is enabled every 6 hours and
ELM is disabled.
