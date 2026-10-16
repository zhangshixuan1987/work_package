# Ensemble perturbation analysis

`verify_ens_perturb.ipynb` checks a set of perturbed EAM initial conditions
against the unperturbed base state. It is independent of how the perturbations
were made, so it applies to the `ensemble_init/` sets and to either method in
`ensemble_generation/`.

## What it checks

1. **Control member:** `EN00` against the base file. It compares SHA-256 hashes
   and reports the maximum absolute difference, which should be zero.
2. **Perturbation statistics:** pooled and per-member moments of
   `member - base`, using a fixed-seed random sample of grid points
   (`sample_per_member`, `seed`).
3. **Figures:** a histogram with a normal Q-Q plot of the perturbations, and
   the vertical profile of the perturbation standard deviation.

The notebook currently analyzes `T` only (`variables=("T",)`). Other
variables need matching input files.

## Inputs

The notebook expects single-variable files for the base state and each member
in one directory:

```text
${data_path}/${init_date}/diag/T.base.eam.i.<date>_SE.nc
${data_path}/${init_date}/diag/T.EN00.eam.i.<date>_SE.nc
${data_path}/${init_date}/diag/T.EN##.eam.i.<date>_SE.nc
```

For the S2D initial conditions, `ensemble_init/*/process_diag.eami.sh` writes
these files, on both the native SE grid (`_SE.nc`) and a 360x180 regridded grid
(`_360x180.nc`). For other sources, extract the variable into files that follow
the same naming.

Set these in the driver cell (`if __name__ == "__main__":`):

| Setting | Meaning |
| --- | --- |
| `data_path` | Root of the ensemble IC output, for example `.../v3.LR.S2D.ENSINT` |
| `init_date` | Initialization date directory, for example `1980-05-01-00000` |
| `grid_tag` | `_SE.nc` (native) or `_360x180.nc` (regridded) |
| `diag_path` | Where `ens_ic_diag_out/<init_date>/` is written |
| `fig_path` | Optional extra copy of the figures (they are always written to the output directory) |
