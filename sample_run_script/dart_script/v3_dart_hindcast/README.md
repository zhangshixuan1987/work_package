# v3_dart_hindcast run scripts

Snapshot of the E3SMv3 hindcast scripts from
`/qfs/people/zhan391/e3sm_dart_work/run_script/v3_dart_hindcast`. Only scripts
(`*.sh`), namelists (`*.nml`) and READMEs are included; runtime state, job logs
and `trash/` folders are not. Paths inside the configs still point to the
original run directories.

The hindcast experiments in `1_ctrl_hindcast`, `2_capt_hindcast` and
`3_dart_hindcast` run the same workflow: the step scripts, `README.md` and
`workflow_lib/` are identical, and only `create_and_setup_case.sh` differs. The workflow is documented in
`3_dart_hindcast/en10/README.md`.

| Folder | Experiments |
|---|---|
| `1_ctrl_baseline` | `en01`, `en01cpl` |
| `1_ctrl_hindcast` | `en10/{20111216,20111221,20111226,20120101}`, `en10cpl`, `en10s1`, `en10s1cpl` |
| `2_capt_hindcast` | `en10`, `en10cpl`, `en10s1`, `en10s1cpl` |
| `3_dart_hindcast` | `en10`, `en20`, `en40`, `en40s1`, `en40cpl`, `en40s1cpl` |
