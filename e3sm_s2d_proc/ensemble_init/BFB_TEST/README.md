# Bit-for-bit (BFB) regression test

Checks that `gen_ensemble_s2d_init.sh` still reproduces the archived initial
conditions after a change to the generator or `create_and_setup_case.sh`.
One array task per set.

| File | Purpose |
| --- | --- |
| `bfb_test.sh` | Slurm array driver: regenerates one date per set into an isolated scratch directory, then runs `bfb_compare.sh`. |
| `bfb_compare.sh` | Compares a candidate directory with its reference: `cmp` first, then `cprnc` for NetCDF files that differ only in metadata. Can be run on its own. |

## Running

Submit from this directory:

```bash
sbatch bfb_test.sh               # all three sets
sbatch --array=1 bfb_test.sh     # one set
BFB_DATE=2010-05-01 sbatch --array=0-1 bfb_test.sh
```

| Array task | Set directory | Reference set | Files | Default date |
| --- | --- | --- | --- | --- |
| 0 | `../BruteForce/luke_cfsr/` | `BruteForce` | 20 | 2015-11-01 |
| 1 | `../JRA55_FOSIRL/` | `JRA55-FOSIRL` | 20 | 2015-11-01 |
| 2 | `../4DEnVarOcn.v1/` | `4DEnVarOcn` | 23 | 1995-05-01 |

`BFB_DATE` overrides the date for every task; 4DEnVarOcn exists only for May
1980-1995.

Output goes to
`/pscratch/sd/z/zhan391/e3sm_project/bfb_s2d_init/job_<jobid>_<task>/<set>/<date>-00000/`;
the reference under `/global/cfs/cdirs/e3sm/zhan391/v3.LR.S2D.ENSINT/` is only
read. The test refuses to reuse an existing scratch directory.

To recheck an earlier run without regenerating:

```bash
./bfb_compare.sh \
  /global/cfs/cdirs/e3sm/zhan391/v3.LR.S2D.ENSINT/BruteForce/2015-11-01-00000 \
  /pscratch/sd/z/zhan391/e3sm_project/bfb_s2d_init/job_58599001_0/BruteForce/2015-11-01-00000
```

## Result

`BFB_REPORT.txt` in the candidate directory lists each file as
`STRICT_MATCH`, `CPRNC_MATCH`, `MISMATCH`, `MISSING`, or `EXTRA`. The test
passes when there are no mismatches, missing, or extra files; `PASS (strict)`
means every file is byte-identical, `PASS (cprnc)` means some differ only in
metadata. A failing file's cprnc output is kept as `<file>.cprnc.txt`.
