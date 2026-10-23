# E3SM S2D Bundled Ensemble Hindcast Workflow

## Overview

This workflow sets up and runs E3SM S2D ensemble hindcast simulations on a
Slurm-based HPC system. The top-level driver is:

```bash
run_e3sm_ens_bundle.sh
```

The workflow is designed for experiments with multiple initialization dates and
multiple ensemble members. It uses a bundled execution strategy: one large Slurm
allocation runs many already-configured E3SM cases concurrently with
`case.submit --no-batch`, reducing queue overhead.

The workflow has three main modes:

1. Compile and initial setup
2. Continuation setup
3. Bundled ensemble execution

Configuration is centralized in:

```bash
create_and_setup_case.sh
```

## Sample run sets

Each subdirectory is a complete, self-contained run workflow for one
initial-condition set from `../ensemble_init/`, copied from a production
hindcast batch. All use 10 members, 2-year forecasts, `custom-7_1x2_nyears`
(7 nodes per member), and a 280-node bundle that runs 40 members at once.

| Directory | `my_leadref` | Sample batch | Lead dates | `my_refcase1` (atm) / `my_refcase2` (other) | Flags as saved |
| --- | --- | --- | --- | --- | --- |
| `4DEnVarOcn.v0/` | `4DEnVarOcn` | 1980-1981 | May, Nov | `assm_anom_ocn_4DEnVar_loc100_v3` / same | initial compile/setup |
| `BruteForce/` | `BruteForce` | 1990-1995 | Feb, Aug | `v3.LR.amip_0101.HICCUP.atm_era5` / `20251024_s2d_spinup` | bundled run, restarts one month after the lead date |
| `JRA55_FOSIRL/` | `JRA55-FOSIRL` | 2017-2018 | May, Nov | `v3.LR.amip_0101.HICCUP.atm_era5` / `20251024_s2d_spinup` | bundled run |
| `RTENSO2026/` | `RTENSO2026` | 2015 | Feb, May, Aug, Nov | `v3.LR.amip_0101.HICCUP.atm_era5` / `20251024_s2d_spinup` | initial compile/setup |

`4DEnVarOcn.v0/` and `BruteForce/` use the current driver; their
`run_e3sm_ens_bundle.sh` differ only in `#SBATCH` account and wall time, and
share the same case template. `JRA55_FOSIRL/` uses an earlier driver and case
template with the compile/setup, continuation, and run modes but no
short-spinup phase (`do_short_spinup` and `do_skip_completed_spinup` are not
used).

`RTENSO2026/` is `BruteForce/` reading the RTENSO2026 initial conditions
(`/global/cfs/cdirs/e3sm/S2S2D/RTENSO2026/initial_conditions/`, see
`../ensemble_init/RTENSO2026/`), which keep one directory per component:
`YYYYMMDD/atm/`, `ocn/`, and `lnd/`, `rof/`, `ice/`, `cpl/<method>/`. Each batch
uses one method, `my_ic_method` (`clim`, `obs`, or `hybrid_ic`; in the case
name), and optionally another sea ice, `my_ice_method` (for example `nsidc`).
The 2026 dates are added to `my_leadymd` once their `atm/` exists.

To run a new batch, copy the matching set to a working directory, edit
`my_leadymd` (and `my_restart_ymd` for continuations), the `#SBATCH` account and
time, and the mode flags below.

`spinup/` holds the standalone one-month spin-up bundle, and `stream_file/`
holds reference MPAS stream files; both are shared by the sets.

## Quick Start

Always run commands from the run-set directory (for example `BruteForce/`)
because `my_workdir` is set from `$PWD`. Choose exactly one workflow stage at a time:

| Stage | `do_compile_setup_only` | `do_continue_setup` | `do_short_spinup` | `do_ensemble_run` |
| --- | --- | --- | --- | --- |
| Initial setup | `true` | `false` | `false` or `true` | `false` |
| Continuation setup | `false` | `true` | `false` | `false` |
| Standard or continuation run | `false` | `false` | `false` | `true` |
| Initial short-spinup chain | `false` | `false` | `true` | `true` |

For setup stages, run:

```bash
./run_e3sm_ens_bundle.sh
```

After setup succeeds, change the flags to the appropriate run stage and submit:

```bash
sbatch run_e3sm_ens_bundle.sh
```

Do not enable `do_continue_setup` and `do_ensemble_run` together. Continuation
setup also requires `do_short_spinup=false`.

## Main Files

* `run_e3sm_ens_bundle.sh` - top-level driver and Slurm batch script.
* `create_and_setup_case.sh` - user-facing configuration file.
* `compile_and_setup_e3sm_<resolution>.<compset>.sh` - E3SM case creation,
  setup, build, and runtime configuration template.
* `e3sm_boundle_extra.sh` - placeholder runtime hook for advanced monitoring or
  control logic during bundled runs.
* `e3sm_boundle_cycling.sh` - placeholder post-run hook for future cycling or
  post-processing steps.
* `ssp245_user_nl_eam.txt` and `ssp245_user_nl_elm.txt` - forcing overrides
  appended for simulations initialized after 2014.
* `stream_file/` - reference MPAS-Ocean and MPAS-Seaice stream files from
  component experts. These are not directly required by the workflow.

Note: these hook scripts are currently placeholders. The filenames currently use
`boundle`, and the driver expects those exact names.

## Configuration

Edit `create_and_setup_case.sh` before running. The most important settings are:

### Ensemble and Dates

* `my_ensnum` - number of ensemble members.
* `my_leadymd` - initialization dates, for example `1980-05-01`.
* `my_leadtod` - initialization time of day in seconds.
* `my_leadref` - initialization source, for example `BruteForce`.

### E3SM Case

* `my_e3sm_code` - E3SM source tree.
* `my_compset` - E3SM compset.
* `my_resolution` - E3SM resolution.
* `my_runtype` - usually `hybrid` for this S2D workflow.

### Paths and Restart Sources

* `my_runpath` - root directory where lead-date/member cases are created.
* `my_refdir` - root directory containing initial restart files.
* `my_refcase1` - atmosphere initial condition case name.
* `my_refcase2` - non-atmosphere restart case name.

### Run Length and Archiving

* `my_runopt`, `my_runn` - forecast segment length.
* `my_restopt`, `my_restn` - restart output frequency.
* `my_short_archive` - controls short-term archiving.
* `my_walltime` - case walltime for the forecast segment.

### HPC Layout

* `my_machine` - E3SM machine name.
* `my_project` - allocation/project.
* `my_jobqueue` - Slurm queue/QOS setting used by CIME.
* `my_job_nnodes` - total nodes in the bundled Slurm allocation.
* `my_min_nodes_per_sim` - nodes required by one ensemble member.
* `my_layout` - PE layout string used by the template script.

The bundled run concurrency is computed as:

```bash
NMAXPS = my_job_nnodes / my_min_nodes_per_sim
```

For example, `my_job_nnodes=280` and `my_min_nodes_per_sim=7` allows up to 40
concurrent member runs inside one Slurm allocation.

The total number of cases is:

```text
number of lead dates x my_ensnum
```

If the total case count is larger than `NMAXPS`, the driver runs cases in
waves. Ensure the `#SBATCH --nodes` value in `run_e3sm_ens_bundle.sh` agrees
with `my_job_nnodes`.

## Workflow Control Flags

The user normally edits only these flags in `create_and_setup_case.sh`:

* `do_compile_setup_only`
* `do_continue_setup`
* `do_short_spinup`
* `do_ensemble_run`
* `my_continue_recompile`

The driver derives the internal flags `do_e3sm_compile` and
`do_ensemble_setup` from these choices.

`my_continue_recompile` is used only during continuation setup. Set it to
`true` when `my_layout` changed; it is ignored during bundled execution.

### Mode 1: Compile and Initial Setup

Use this mode to create all cases, build `e3sm.exe` once, copy/reuse that
executable for the other members, copy initial restart files, patch namelists,
and run `case.setup`.

```bash
export do_compile_setup_only=true
export do_continue_setup=false
export do_short_spinup=false
export do_ensemble_run=false
```

Set `do_short_spinup=true` here only when the first run should use the configured
short-spinup phase.

Run from the workflow directory:

```bash
./run_e3sm_ens_bundle.sh
```

This mode performs:

* run preflight checks for required configuration, templates, and input files
* create and build `EN00` for the first lead date
* verify `EN00/build/e3sm.exe`
* create the remaining ensemble member cases
* copy the shared executable into the other build directories
* set up all configured lead dates and members
* copy required restart and rpointer files into each member archive
* update `user_nl_eam`, `user_nl_elm`, `user_nl_mosart`,
  `user_nl_mpaso`, and `user_nl_mpasi`
* patch MPAS-Ocean and MPAS-Seaice stream files

Setup is parallelized with an internal limit of 10 background setup jobs.

### Short Spinup

Short spinup is an initial-run option, not a continuation option. When
`do_short_spinup=true`, initial setup uses `my_spinup_*` time-step, stop, and
restart settings. During bundled execution, each member runs the short-spinup
phase and then the remaining main phase. If
`do_skip_completed_spinup=true`, members with a complete spinup restart set skip
the first phase.

Set `do_short_spinup=false` for continuation setup and continuation execution.

### Mode 2: Continuation Setup

Use this mode when a previous segment has completed and you want to prepare the
same cases for a continuation run.

```bash
export do_compile_setup_only=false
export do_continue_setup=true
export do_short_spinup=false
export do_ensemble_run=false
export my_continue_recompile=false
```

Set `my_restart_ymd` so that it corresponds element-by-element to
`my_leadymd`, and set `my_restart_tod`.

`do_short_spinup` must be `false` during continuation setup. The driver exits
with an error if continuation setup and short spinup are enabled together.

Run:

```bash
./run_e3sm_ens_bundle.sh
```

This mode performs:

* run preflight checks for restart-date configuration and required files
* locate each member restart directory under its archive
* copy restart files into the run directory
* copy recent history files needed by continuation runs when available
* set `CONTINUE_RUN=TRUE`
* align `RUN_REFDATE` and `RUN_REFTOD` with each case's original lead date
* update continuation run length, restart frequency, walltime, and queue
* rerun `case.setup`

If `my_continue_recompile=true`, the driver also enables the compile path before
continuation setup. Use this when changing `my_layout`. The driver rebuilds
`e3sm.exe` once using the requested layout and reuses it for the other cases.
Each lead date starts from a fresh template so later dates do not inherit the
first lead date's restart settings.

### Mode 3: Bundled Ensemble Run

Use this mode only after compile/setup or continuation setup has completed.

```bash
export do_compile_setup_only=false
export do_continue_setup=false
export do_short_spinup=false
export do_ensemble_run=true
```

Submit the bundled Slurm job:

```bash
sbatch run_e3sm_ens_bundle.sh
```

This mode performs:

* run preflight checks for expected case directories and hook scripts
* loop over all configured lead dates and ensemble members
* enter each member's `case_scripts` directory
* launch the case with `./case.submit --no-batch`
* run up to `NMAXPS` members concurrently
* warn and continue if an individual member run fails
* print a final summary of failed member jobs and log paths, if any
* source placeholder hook `e3sm_boundle_extra.sh` during the run-monitoring
  section
* source placeholder hook `e3sm_boundle_cycling.sh` after the bundled run
  section

## Directory Layout

For each lead date, the driver creates a lead-specific run path:

```bash
${my_runpath}/${my_leadcase}_${YYYYMMDDHH}
```

Each ensemble member is placed under:

```bash
${my_runpath}/${my_leadcase}_${YYYYMMDDHH}/EN##
```

Important member subdirectories are:

* `case_scripts/` - CIME case scripts.
* `build/` - executable and build products.
* `run/` - active run directory.
* `archive/` - short-term archive and restart output.

## Initial Restart File Expectations

For initial setup, restart files are read from:

```bash
${my_refdir}/${START_DATE}-${START_TOD}
```

The workflow expects component restart filenames based on these rules:

* Atmosphere uses `my_refcase1` and member-specific `EN##.eam.i` files.
* Ocean and sea ice use `my_refcase2` with MPAS-style `_` date separators.
* Land, river, and coupler use `my_refcase2` with standard `-` date separators.
* Non-atmosphere components also require `rpointer.<component>` files.

The copied restart files are staged into:

```bash
${ARCHIVE_DIR}/rest/${START_DATE}-${START_TOD}
```

The atmosphere initial condition is written to `ncdata`, land to `finidat`, and
river to `finidat_rtm` through the corresponding `user_nl_*` files.

## MPAS Stream and Namelist Customization

During initial setup, the driver patches:

* `streams.ocean`
* `streams.seaice`

The patched versions are also copied into:

```bash
SourceMods/src.mpaso
SourceMods/src.mpassi
```

The driver also updates MPAS time steps in:

* `user_nl_mpaso`
* `user_nl_mpasi`

When `do_short_spinup=false`, the workflow currently resets these to the
default values:

```bash
my_mpaso_dt="00:30:00"
my_mpassi_dt=1800
```

## Simulations Beyond 2014

For initialization years after 2014, the workflow appends:

```bash
ssp245_user_nl_eam.txt
ssp245_user_nl_elm.txt
```

This switches prescribed forcing beyond the historical forcing period. Users
should treat diagnostics around the 2014-2015 boundary carefully because the
forcing source changes from historical to SSP245 scenario data.

## Typical Usage

### First Segment

1. Edit `create_and_setup_case.sh`.
2. Set:

```bash
export do_compile_setup_only=true
export do_continue_setup=false
export do_short_spinup=false
export do_ensemble_run=false
```

3. Run setup:

```bash
./run_e3sm_ens_bundle.sh
```

4. After setup completes, set:

```bash
export do_compile_setup_only=false
export do_continue_setup=false
export do_short_spinup=false
export do_ensemble_run=true
```

5. Submit the bundled run:

```bash
sbatch run_e3sm_ens_bundle.sh
```

### Continuation Segment

1. Set `my_restart_ymd` and `my_restart_tod`.
2. Set:

```bash
export do_compile_setup_only=false
export do_continue_setup=true
export do_short_spinup=false
export do_ensemble_run=false
```

3. Run continuation setup:

```bash
./run_e3sm_ens_bundle.sh
```

4. Then set:

```bash
export do_compile_setup_only=false
export do_continue_setup=false
export do_short_spinup=false
export do_ensemble_run=true
```

5. Submit the bundled run:

```bash
sbatch run_e3sm_ens_bundle.sh
```

## Before Running

Check the following before each stage:

* `my_leadymd` contains the intended initialization dates.
* `my_restart_ymd` has the same number of entries for continuation setup.
* `my_layout` and `my_min_nodes_per_sim` describe the same per-case layout.
* `my_job_nnodes` and the Slurm `#SBATCH --nodes` request agree.
* Setup completed successfully for every lead date and ensemble member before
  enabling `do_ensemble_run=true`.
* During a run, member logs are written under each case's `case_scripts/`
  directory as `e3sm_EN##.log.o<jobid>`.

## Notes

* The setup and bundled execution steps are intentionally separate. Submitting
  the top-level script with `sbatch` runs only the mode selected by the flags.
* `scripts/` is generated by the workflow and stores temporary helper scripts
  used to create/setup member cases.
* The compile-once strategy assumes the same executable is valid for all members
  and lead dates in the configured experiment.
* `my_restart_ymd` must contain exactly one restart date for every entry in
  `my_leadymd`.
* A continuation run is a two-step workflow: first run continuation setup with
  `do_continue_setup=true`, then run the cases with `do_ensemble_run=true`.
* The bundled run assumes all target cases already exist and have completed
  setup successfully.
