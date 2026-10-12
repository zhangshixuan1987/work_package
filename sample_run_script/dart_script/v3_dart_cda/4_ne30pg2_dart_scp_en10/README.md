# E3SM–DART Strongly Coupled AMIP Ensemble Workflow

This directory contains a restart-safe, Slurm-driven E3SM–DART cycling data-assimilation workflow for the 10-member strongly coupled AMIP case. The numbered scripts are the user-facing entry points. Internal workers and templates live under `workflow_lib/` and should not normally be executed directly.

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
4. Coupled forecast–assimilation cycling
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

- Slurm resources: `my_task_per_node`, `my_job_nnodes`, `my_project`, `my_jobqueue`, and `my_walltime`.
- Ensemble configuration: `my_ensnum`, `my_nodes_per_member`, setup concurrency, and forecast retry settings.
- Model configuration: `my_e3sm_code`, `my_runtype`, `my_compset`, `my_resolution`, `my_runpath`, and `my_casename`.
- Initial state: `my_casedate`, `my_casetod`, `my_refcase`, `my_refdate`, `my_reftod`, `my_refdir`, and component restart paths.
- DART configuration: EAM and ELM enable switches, component code and observation paths, `my_e3sm_cycle_hours`, `my_eam_dart_cycle_hours`, `my_elm_dart_cycle_hours`, component DA end times, diagnostic range, and the shared E3SM timeline.
- Cycling behavior: `my_cycles_per_job`, minimum cycle runtime, shutdown margin, and handoff concurrency.

Exceptional assimilation cycles can override localization cutoff and inflation damping in the `my_eam_dart_cycle_overrides` table. Unlisted cycles continue to use the ensemble-size defaults.

Set `strongly_coupled_on="on"` to apply both the ATM and LND DA profiles. With `lnd_da_use_sequential_prior_post=".true."`, EAM runs first and writes the current sequential-prior observation sequence; ELM validates and consumes that file afterward. A static ELM preflight and transaction marker are completed before EAM starts, so any later EAM or ELM failure forces a complete forecast rebuild. With sequential-prior use set to `.false.`, ELM retains the direct-observation path and may run in parallel with EAM; this direct mode requires `lnd_da_strongly_coupled=".false."`, `lnd_da_state_model="Land"`, and `lnd_da_obs_model="Land"`. Set `strongly_coupled_on="off"` to retain the source DART namelist values and the direct ELM observation stream.

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

### Step 4 — Coupled cycling DA

```bash
sbatch 4_run_dart_e3sm_cycleda.sh
```

Step 4 requires the initial perturbation record and verifies that the Slurm allocation matches `my_job_nnodes`. Before running a cycle, it checks whether the configured DA end time has already been reached.

Each cycle advances the coupled E3SM ensemble by `my_e3sm_cycle_hours`, validates the forecast products, and independently determines whether EAM-DART and ELM-DART are due. A component runs when its configured interval divides the elapsed time from the DA start. After every due analysis validates, the workflow transactionally commits the next-cycle state, updates the global cycle counter, and writes:

```text
runtmp/status/cycle_complete.<valid-time>
```

The shared timeline is tracked by `my_e3sm_completed_cycles`, `my_e3sm_cycle_hours`, `my_e3sm_start_date`, and `my_e3sm_end_date`/`my_e3sm_end_tod`. The counter advances only after a complete forecast/handoff transaction. Component assimilation is controlled by `my_eam_dart_da` and `my_elm_dart_da`, with frequencies set by `my_eam_dart_cycle_hours` and `my_elm_dart_cycle_hours`. Each component interval must be an integer multiple of `my_e3sm_cycle_hours`. With the default 6-hour EAM and 24-hour ELM settings, EAM runs at every forecast endpoint and ELM runs at midnight relative to a midnight DA start. When both are due in direct mode, each filter receives 10 of the configured 20 nodes and runs concurrently. In sequential-prior mode, EAM and then ELM each receive all 20 nodes because their executions cannot overlap. When only one component is due, it receives all 20 nodes. A component that is enabled but not scheduled records `da_mode=not_due`; after its `my_eam_dart_end_*` or `my_elm_dart_end_*` timestamp it records `da_mode=ended`. Component end times are inclusive and must not exceed the shared E3SM end time. Handoff proceeds only after every due component validates successfully.

The shared model archive remains `my_dart_archive_dir`. EAM and ELM DART setup, inflation, and run products are isolated under `my_eam_dart_run_dir` and `my_elm_dart_run_dir`, respectively.

ELM prior and posterior inflation restarts are carried forward independently from the most recent earlier ELM directory containing every required file. When the configured land inflation flavors are both zero, `fill_inflation_restart` is neither required nor staged.

The driver can run multiple cycles per allocation using `my_cycles_per_job`. It stops conservatively when insufficient wall time remains and can recover or submit a continuation job without duplicating active continuations.

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

Step 6 validates upstream completion records across the entire requested cycle range before starting. Closest mode additionally requires every ensemble restart file. Its EAM assimilation count is derived from the requested end time and `my_eam_dart_cycle_hours`, rather than from the shared forecast-cycle counter. Mode-specific workers live under `workflow_lib/diagnostics/` and cannot be executed directly. Each mode has a separate completion record containing its mode, range, case, and ensemble size. A matching completion record skips only that mode, allowing failed or newly enabled modes to be run independently.

### Step 7 — History post-processing

Edit the settings at the top of `7_run_post_hist.sh`:

```bash
POST_MODE="all"                # base, monthly, or all
POST_START="2011-12-01-00000"
POST_END="2011-12-26-00000"
MAX_CONCURRENT_WORKERS=4
OVERWRITE_EXISTING="FALSE"
```

Individual EAM/ELM products can be enabled or disabled using the `RUN_*` switches. Submit with:

```bash
sbatch 7_run_post_hist.sh
```

In `all` mode, base products complete before monthly aggregation begins. In `monthly` mode, existing base completion records and outputs must validate. Every member/product uses its own lock and completion record. Outputs must be nonempty, readable NetCDF files before completion is recorded.

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

Keep `REQUIRE_UPSTREAM_MARKERS=TRUE` for normal operation. Set it to `off` only when deliberately processing a legacy restart archive whose files have been independently validated but whose old completion markers were not retained.

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
