# S2D ensemble initial conditions

Each subdirectory is a self-contained working directory for one
initial-condition set:

| Directory | `workflow_mode` | Non-atmosphere state |
| --- | --- | --- |
| `JRA55_FOSIRL/` | `jra55-fosirl` | Ocean, ice, land, river, and coupler state from the matching JRA55-forced (FOSI) restart. |
| `BruteForce/luke_cfsr/` | `bruteforce-cfsr` | Same restart, with the MPAS-Ocean state replaced by Luke Van Roekel's CFSRv2 analysis (`20251024_s2d_spinup.mpaso.rst.<model date>.cfsrv2.nc`). |
| `Reanalysis/` | `reanalysis-cfsr` | Same restart, with the MPAS-Ocean state replaced by Luke Van Roekel's calendar-date CFSR ocean ICs (`mpas_<YYYYMMDD>.cfsr.ic.nc`). |
| `4DEnVarOcn.v1/` | `4denvar-ocn` | 4DEnVar atmosphere, ocean, land, river, sea-ice, coupler, and atmosphere-restart files (`assm_full_ocn_4DEnVar_loc_m100_v3`, May 1980-1995, selected from the matching decade directory). |

`gen_ensemble_s2d_init.sh` and `create_and_setup_case.sh` are identical in
`JRA55_FOSIRL/`, `BruteForce/luke_cfsr/`, and `Reanalysis/`; `4DEnVarOcn.v1/`
has its own copies set up for the `assm_full` 4DEnVar sources.
`submit_all_year.bash` (and, for 4DEnVarOcn, `process_diag.eami.sh`) carry the
set-specific settings.

`Reanalysis/` is **on hold** until its inputs are complete: the spin-up
restarts in `my_fosdir1` (`/global/cfs/cdirs/e3sm/lvroekel/archive/rest`) are
no longer at that path, 2026 has no spin-up restarts or ERA5 atmosphere ICs,
and the sea-ice files are not yet on the `IcoswISC30E3r5` mesh. The scripts
pass syntax and preflight path checks but have not generated any output.

`Reanalysis/` reads `v3.LR.S2D.ENSINT/Reanalysis_ICs/mpas_<YYYYMMDD>.cfsr.ic.nc`
(2015 Feb/May/Aug/Nov and 2026 Feb/May/Aug so far) through `my_ocndir4` and
`my_ocncase4` (`my_ocnopt=LukeVanRoekel3`) and writes to
`v3.LR.S2D.ENSINT/Reanalysis/`. As in `BruteForce/luke_cfsr/`, the ocean file is
appended onto the matching spin-up `mpaso.rst` (`ncks -A`), and the land, river,
sea-ice, and coupler restarts come from the spin-up model year
(`year - 1980 + 353`).

Reanalysis-based sea ice is under development. Setting `ice="NSIDCv6"` in the
`Reanalysis/submit_all_year.bash` settings (`my_iceopt`) replaces the spin-up
`mpassi.rst` with
`v3.LR.S2D.ENSINT/Reanalysis_ICs/mpassi.EC30to60E2r2.iceArea_from_NSIDCv6_spunup.<YYYYMMDD>.nc`
(`my_icedir`, `my_icecase`). It is off by default: the current files are on
the `EC30to60E2r2` mesh (236,853 cells), not `IcoswISC30E3r5` (465,044 cells),
and the generator stops when the sea-ice and ocean meshes differ. Update
`my_icecase` once the `IcoswISC30E3r5` files are available.

`BruteForce/ben_kirtman/` is the earlier brute-force generator, kept for
reference.

`4DEnVarOcn.v0/` is the earlier 4DEnVar setup, kept unchanged as a record. It
read the `assm_anom_ocn_4DEnVar_loc100_v3` assimilation (May/November
1980-2015, with `WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_4DEnVar` atmosphere ICs
before 1990) and produced 72 dates. Those sources and outputs now live under
`v3.LR.S2D.ENSINT/4DEnVar_ICs_SHIPF.OLD/` and `v3.LR.S2D.ENSINT/4DEnVarOcn.OLD/`;
the scripts still point at the current `4DEnVar_ICs_SHIPF/` and
`4DEnVarOcn/` paths, so re-running them requires repointing `my_daroot` and
the output directory first.

| | `4DEnVarOcn.v0/` | `4DEnVarOcn.v1/` (current) |
| --- | --- | --- |
| Assimilation case | `assm_anom_ocn_4DEnVar_loc100_v3` | `assm_full_ocn_4DEnVar_loc_m100_v3` |
| Dates | May, Nov 01, 1980-2015 (72) | May 01, 1980-1995 (16) |
| Sources | `4DEnVar_ICs_SHIPF.OLD/` | `4DEnVar_ICs_SHIPF/` |
| Output | `4DEnVarOcn.OLD/` | `4DEnVarOcn/` |

`BFB_TEST/` checks that the generator reproduces established initial
conditions bit-for-bit with `cprnc`.

## Output

Generated initial conditions are under
`/global/cfs/cdirs/e3sm/zhan391/v3.LR.S2D.ENSINT/<set>/<YYYY-MM-DD>-00000/`:

| `<set>` | Dates available |
| --- | --- |
| `BruteForce` | Feb, May, Aug, Nov 01, 1980-2018 (156 dates) |
| `JRA55-FOSIRL` | Feb, May, Aug, Nov 01, 1980-02 to 2019-05 (157 dates) |
| `4DEnVarOcn` | May 01, 1980-1995 (16 dates) |
| `Reanalysis` | none yet (on hold; 2015 and 2026 ocean ICs available) |

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

`RTENSO2026/` is in preparation: initial conditions for the 2026 hindcasts
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
