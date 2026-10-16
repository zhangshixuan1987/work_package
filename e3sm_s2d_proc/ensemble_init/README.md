# S2D ensemble initial conditions

Each subdirectory is a self-contained working directory for one
initial-condition set:

| Directory | `workflow_mode` | Non-atmosphere state |
| --- | --- | --- |
| `jra55_fosirl/` | `jra55-fosirl` | Ocean, ice, land, river, and coupler state from the matching JRA55-forced (FOSI) restart. |
| `bruteforce/luke_cfsr/` | `bruteforce-cfsr` | Same restart, with the MPAS-Ocean state replaced by Luke Van Roekel's CFSRv2 analysis. |
| `4DEnVarOcn/` | `4denvar-ocn` | 4DEnVar atmosphere, ocean, land, river, sea-ice, coupler, and atmosphere-restart files (May/November 1980-2015, selected from the matching decade directory). |

`gen_ensemble_s2d_init.sh` and `create_and_setup_case.sh` are identical across
the three sets; `submit_all_year.bash` (and, for 4DEnVarOcn,
`process_diag.eami.sh`) carry the set-specific settings.

`bruteforce/ben_kirtman/` is the earlier brute-force generator, kept for
reference. `bfb_test/` checks that the generator reproduces established
initial conditions bit-for-bit with `cprnc`.

## Output

Generated initial conditions are under
`/global/cfs/cdirs/e3sm/zhan391/v3.LR.S2D.ENSINT/<set>/<YYYY-MM-DD>-00000/`:

| `<set>` | Dates available |
| --- | --- |
| `BruteForce` | Feb, May, Aug, Nov 01, 1980-2018 (156 dates) |
| `JRA55-FOSIRL` | Feb, May, Aug, Nov 01, 1980-02 to 2019-05 (157 dates) |
| `4DEnVarOcn` | May, Nov 01, 1980-2015 (72 dates) |

Each date directory holds the hindcast inputs for one initialization, for
example in `BruteForce/` (20 files):

- `<refcase1>.EN00-EN09.eam.i.<date>.nc`: atmosphere members (`EN00`
  unperturbed)
- `<refcase2>.{cpl.r,elm.r,mosart.r}.<date>.nc` and
  `<refcase2>.{mpaso,mpassi}.rst.<date with _>.nc`: coupler, land, river,
  ocean, and sea-ice restarts
- `rpointer.{drv,ice,lnd,ocn,rof}`

`hindcast_run/` reads these through `my_refdir`, `my_refcase1`, and
`my_refcase2`.

## In preparation

`rtenso2026/` is in preparation: initial conditions for the 2026 hindcasts
(Feb, May, Aug, and Nov 01; the same 2015 dates test the sensitivity), with
ERA5 atmosphere and CFSR ocean from collaborators (as in the brute-force set),
sea-ice, land, and river restarts taken here from a BruteForce year closest
to the same-date climatology (`clim`) or to the ERA5 anomaly (`obs`), and DART member-subset selection; see its
README.

## Atmospheric ensemble perturbations with DART

The three established sets build the ensemble members by calling DART: `gen_ensemble_s2d_init.sh`
copies the base EAM initial condition (`eam.i`) to every member, stages the
members in a DART run directory, and runs DART's `perturb_single_instance`
from `${my_dart_code}/models/eam-se/work`, configured by
`perturb_single_instance_nml` in `eam_pert_input.nml`, to perturb the
temperature field. No observations or inflation are involved. The DART code, environment, and task count are set in the
`my_dart_*` block of `create_and_setup_case.sh`. Member 0 (`EN00`) is the
unperturbed control.

## Running

`submit_all_year.bash` is the only driver. Edit its `USER SETTINGS` block to
choose the workflow, years, months, concurrency, and overwrite policy. Then run:

```bash
./submit_all_year.bash
```

Dates with unavailable inputs can be omitted in the same settings block. For
example, this skips only February 2011:

```bash
skip_years=(2011)
skip_months=(2)
```

If only `skip_years` is populated, every selected month in those years is
skipped. If only `skip_months` is populated, those months are skipped in every
selected year. Leave both arrays empty to skip nothing.

It checks every configured date and submits itself as a Slurm array. To check
without submitting, set `check_only=true` in the same settings block and still
run `./submit_all_year.bash`.

Existing nonempty output directories are protected by default. To deliberately
replace them, set `overwrite_existing=true` in the `USER SETTINGS` block.
The generator still validates all source files before removing any old output.

For a single array mapping/configuration dry run without generating files:

```bash
IC_MODE=bruteforce-cfsr SLURM_ARRAY_TASK_ID=77 DRY_RUN=true \
  bash submit_all_year.bash
```
