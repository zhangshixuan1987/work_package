# DART hindcasts

One directory per experiment, all running the same workflow (step scripts,
`workflow_lib/` and README are identical; only `create_and_setup_case.sh`
differs). `en10/README.md` documents the workflow; `en10` is the template for
new experiments.

| Experiment | DA case (initial conditions) | Start | Run type | Nodes/member |
|---|---|---|---|---|
| `en10`      | DARTEN10_F20TR   | 2012-01-01 | AMIP     | 10 |
| `en20`      | DARTEN20_F20TR   | 2012-01-01 | AMIP     | 10 |
| `en40`      | DARTEN40_F20TR   | 2012-01-01 | AMIP     | 10 |
| `en40s1`    | DARTEN40S1_F20TR | 2012-06-01 | AMIP     | 10 |
| `en40cpl`   | DARTEN40_F20TR   | 2012-01-01 | Full-CPL | 40 |
| `en40s1cpl` | DARTEN40S1_F20TR | 2012-06-01 | Full-CPL | 40 |

All run 2 months with monthly restarts.
