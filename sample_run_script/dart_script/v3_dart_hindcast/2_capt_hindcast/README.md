# CAPT hindcasts

Ensemble hindcasts initialized from ERA5 (`my_ic_source=perturbed`): the EAM
state is the 2012-01-01 `v3.LR.amip_0101` file with `PS,U,V,T,Q,CLDLIQ,CLDICE`
from `v3.LR.ne30pg2.ERA5.eam.i` at the start date, then perturbed with DART.
Every experiment runs the shared hindcast workflow documented in
`../3_dart_hindcast/en10/README.md`, with two extra steps first:

```bash
sbatch 0a_stage_initial_conditions.sh
sbatch 0b_perturb_initial_conditions.sh
sbatch 1_run_e3sm_ensemble_setup.sh       # then Steps 1-6
```

| Experiment | Case | Start | Run type |
|---|---|---|---|
| `en10`      | CAPTEN10_15day   | 2012-01-01 | AMIP     |
| `en10cpl`   | CAPTEN10_15day   | 2012-01-01 | Full-CPL |
| `en10s1`    | CAPTEN10S1_15day | 2012-06-01 | AMIP     |
| `en10s1cpl` | CAPTEN10S1_15day | 2012-06-01 | Full-CPL |

All run 2 months with monthly restarts.
