# Control hindcasts

Ensemble hindcasts from a perturbed control state (`my_ic_source=perturbed`).
Every experiment runs the shared hindcast workflow documented in
`../3_dart_hindcast/en10/README.md`, with two extra steps first:

```bash
sbatch 0a_stage_initial_conditions.sh     # stage the control state for every member
sbatch 0b_perturb_initial_conditions.sh   # DART perturbation of the EAM state
sbatch 1_run_e3sm_ensemble_setup.sh       # then Steps 1-6 as for the DART hindcasts
```

| Experiment | Case | Start | Run type | EAM initial state |
|---|---|---|---|---|
| `en10/20111216` | CTRLEN10_30day   | 2011-12-16 | AMIP     | CTRL_F20TR `init/` amip_0101 (not on disk) |
| `en10/20111221` | CTRLEN10_25day   | 2011-12-21 | AMIP     | CTRL_F20TR `init/` amip_0101 (not on disk) |
| `en10/20111226` | CTRLEN10_20day   | 2011-12-26 | AMIP     | CTRL_F20TR `init/` amip_0101 (not on disk) |
| `en10/20120101` | CTRLEN10_15day   | 2012-01-01 | AMIP     | CTRL_F20TR `init/` amip_0101 |
| `en10cpl`       | CTRLEN10_15day   | 2012-01-01 | Full-CPL | CTRL_F20TR `init/` amip_0101 |
| `en10s1`        | CTRLEN10s1_15day | 2012-06-01 | AMIP     | CTRL_F20TR EN01 restart (no longer in archive/rest) |
| `en10s1cpl`     | CTRLEN10s1_15day | 2012-06-01 | Full-CPL | CTRL_WCYCL20TR EN01 restart (no longer in archive/rest) |

All run 2 months with monthly restarts; land, river, sea-ice, ocean and
coupler states come from the `20241109.v3-LR.DATM` spin-up at the start date.
The three dated December starts were never run.
