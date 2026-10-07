# Maintained NCL workflows

This tree contains the actively maintained NCL workflow bundles. Their original
internal directory structure is retained because many scripts use relative
`load` paths.

Run an NCL program through the repository launcher instead of invoking `ncl`
directly:

```bash
diagnostic/script/run_process_ncl.bash \
  diagnostic/script/ncl/rmse_profile/plot_vertical_preshgt.ncl
```

The launcher runs from the program's own directory and then moves new figures
to `${DIAGNOSTIC_OUTPUT_ROOT}/figure/ncl/active/<bundle>/`. It moves generated
NetCDF products and logs to `${DIAGNOSTIC_OUTPUT_ROOT}/data/ncl/active/<bundle>/`.
`DIAGNOSTIC_OUTPUT_ROOT` defaults to
`/compyfs/www/zhan391/e3sm_dart/diag_dart_2026`.
