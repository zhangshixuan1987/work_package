# Initial-land command-line drivers

This directory contains optional thin entry points for unattended or scheduled
restart processing. Reusable implementations live in `diagnostic/util/`;
the Jupyter workflow calls those utilities directly.

## Complete workflow

Use the Bash driver to aggregate and regrid in one reproducible command:

```bash
diagnostic/script/run_process_land_rest.bash \
  --restart /path/to/ELM_restart.nc \
  --source-domain /path/to/land_domain.nc \
  --variables H2OSOI,TSA,SOILM
```

Products default to
`/compyfs/www/zhan391/e3sm_dart/diag_dart_2026/data/analysis_lnd_init/restart/`.
Existing aggregation, SCRIP, weight-map, and regridded files are reused. Pass
`--force` to regenerate them, `--aggregate-only` to skip regridding, or
`--dry-run` to inspect the commands without writing files.

For Slurm, submit the same driver with site-specific resources, for example:

```bash
sbatch --job-name=initial_land_restart \
  --output=initial_land_restart_%j.log \
  --wrap="diagnostic/script/run_process_land_rest.bash \
    --restart /path/to/ELM_restart.nc \
    --source-domain /path/to/land_domain.nc"
```

The job environment must provide NCO/`ncremap` and
`ESMF_RegridWeightGen`. The default Python executable is
`/qfs/people/zhan391/.conda/envs/e3sm_analysis/bin/python`; override it with
`--python` or `PYTHON_BIN`.

## Individual stages

Aggregate restart columns to gridcells:

```bash
python diagnostic/script/initial_land/aggregate_restart.py \
  --restart ELM_restart.nc \
  --src-domain land_domain.nc \
  --out /path/to/agg_gridcell.nc \
  --vars H2OSOI TSA SOILM
```

Conservatively regrid an existing aggregated product:

```bash
python diagnostic/script/initial_land/regrid_restart.py \
  --input /path/to/agg_gridcell.nc \
  --destination-scrip /path/to/cmip6_180x360_scrip.nc \
  --source-scrip /path/to/src_unstructured.scrip.nc \
  --map-file /path/to/map_unstructured_to_180x360.nc \
  --output /path/to/elm_180x360_aave.nc
```
