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
`./4_submit_member_jobs.sh` (submit / `status` / `finalize`): each member runs
as its own CIME job chain, so no DART cycling or `eam.i` output is produced.
