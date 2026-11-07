# Ensemble initial conditions

Builds the 20-member initial ensembles of the F20TR hindcasts (start dates 16, 21, 26 and 31 Dec 2011 and
16, 21, 26 and 31 May 2012). Every ensemble is a reference state plus the same DART-derived atmospheric
perturbations; the three sets differ only in the reference state.

## Workflow

```text
step 1  01_select_dart_20_member_subset.ipynb
        20 of the 40 DART members per season, matching the 40-member mean, spread and covariance
          |
          +--> step 2  02_ic_free_clim_dartpert.ipynb
          |            CLIM_BASE (FREE_CLIM) state                 + DART perturbations  -> FREECLIM_DARTPERT20
          |
          +--> step 3  03_ic_era5_dartfill_dartpert.ipynb
          |      3a    ERA5 atmosphere, other fields = DART mean   -> DART_INIT/ATM_COMPLETE/*.ERA5c.*
          |      3b                                                + DART perturbations  -> HYB_DARTPERT20C
          |
          +--> step 4  04_ic_era5_climfill_dartpert.ipynb
                 4a    ERA5 atmosphere, other fields = CLIM_BASE   -> DART_INIT/ATM_COMPLETE_CLIMFILL/*.ERA5climfill.*
                 4b                                                + DART perturbations  -> HYB_CLIMFILL_DARTPERT20
```

Steps 2-4 are independent of each other; each needs only step 1 (its member lists are recorded in
`util/dart_perturbations.py`, `DART_SEASONS`, and checked against the step-1 CSV when present).

| Step | Reference atmosphere | Aerosols, ozone, chemistry, ... | Land / sea ice | Output (`e3sm_dart/`) | Used by |
|---|---|---|---|---|---|
| 2 | CLIM_BASE | CLIM_BASE | CLIM_BASE | `free_clim_dartpert20_atm_multidate` | CLIM_BASE |
| 3 | ERA5 | DART mean | HYB (DATM) | `hyb_dartpert20c_atm_multidate` | HYB_BASE (DA_LND uses the 3a file) |
| 4 | ERA5 | CLIM_BASE | HYB (DATM) | `hyb_climfill_dartpert20_atm_multidate` | test of the fill fields |

Steps 3 and 4 complete the ERA5 states first: the files `DART_INIT/ATM/v3.LR.ne30pg2.ERA5.eam.i.<date>.nc`
hold only T, U, V, Q, PS, CLDLIQ and CLDICE, and EAM starts every missing field from zero, including the
radiation's water vapour, ozone, N2O and CH4 tracers. The first HYB build from those bare files
(`hyb_dartpert20_atm_multidate`, `HYB_DARTPERT20`) started about 6 K too cold over land and is kept for
reference only.

## Perturbations

For the season's 20 selected DART members S, output EN01 is the reference itself and

    EN(i+1) = reference + alpha * (DART member s_i - mean over S),   i = 1..19

so the 20th selected member enters the mean only. Steps 3 and 4 perturb T, U, V, Q, PS and H2OLNZ; step 2's
existing library perturbed T, U, V, Q and PS. Each date directory holds the 20 `eam.i` files, links to the
reference's land, river, sea-ice and coupler restarts, the member mapping, diagnostics, validation and
metadata.

## Code

* `util/dart_subset_selection.py`: step 1.
* `util/complete_era5_ic.py`: completion of the ERA5 states (steps 3a and 4a).
* `util/dart_perturbations.py`: perturbation, writing and validation (steps 2, 3b and 4b).

Paths are set in each notebook's "User setup" cell. Existing outputs are validated and reused unless
`OVERWRITE = True`.
