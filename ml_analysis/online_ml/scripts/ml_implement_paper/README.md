# Online ML bias-correction paper workflow

This package provides a Jupyter-driven, reproducible comparison of `CTRL`,
`UNET`, and `UNETXTR` against a common reference dataset. The notebook is the
scientific driver; reusable I/O, preprocessing, metrics, and plotting logic is
kept in `src/`.

## Setup

1. Edit `config/paths.yaml` with external file patterns, analysis dates, and
   model-specific variable/dimension names.
2. Ensure the environment provides `xarray`, `dask`, `netCDF4` or `h5netcdf`,
   `numpy`, `pandas`, `matplotlib`, `PyYAML`, and optionally `cartopy`.
3. Run the figure-specific notebooks independently:

   - `notebooks/fig01_global_rmse.ipynb`
   - `notebooks/fig02_vertical_rmse.ipynb`
   - `notebooks/fig03_spatial_error_reduction.ipynb`
   - `notebooks/fig04_stability_ml_forcing.ipynb`
   - `notebooks/fig05_monthly_ml_forcing.ipynb`
   - `notebooks/fig06_monthly_vertical_bias.ipynb`

Each notebook loads only the variables needed for its diagnostic, runs
Dask-backed quality control, writes its processed NetCDF product, and reloads
that product in a separate plotting section. Figure 3 also writes the combined
summary table after loading the Figure 1 diagnostic product.

The workflow assumes experiment and reference fields are already on a common
horizontal grid. It matches common timestamps and interpolates three-dimensional
fields to configured pressure levels. E3SM grid-cell area is used when present;
cosine-latitude weights are the fallback.

## External outputs

Generated artifacts are written outside Git under:

```text
/global/cfs/cdirs/m1867/zhan391/paper_work/online_ml/ml_implement_paper/
├── processed/
├── figures/
└── tables/
```

The notebooks write per-case NetCDF diagnostic caches, publication
figures, and `tables/summary_metrics.csv`. No analysis output should be written
inside this source directory.

## Pressure interpolation

If the data have a one-dimensional pressure coordinate on `lev`, the workflow
uses it directly. For native E3SM hybrid levels, it calculates pressure from
`hyam`, `hybm`, `P0`, and surface pressure. Update these names in the YAML file
when the input convention differs.

## Monthly vertical bias cross sections (Figure 6)

Figure 6 compares model U, V, T, Q against the processed ERA5 monthly reanalysis
specified by `ERA5_PATH`. Only selected months and fields (plus surface pressure)
are loaded. ERA5 and model grids must match exactly. ERA5 standard levels are
selected directly; model fields are interpolated as needed to the common
1000–1 hPa levels. Both are masked below their surfaces before experiment-minus-ERA5
biases are averaged over longitude. Equivalent wind/humidity unit labels are
normalized. Model inputs prefer monthly time series and fall back to exact
single-month products under `post/atm/180x360_aave/clim/1yr` when needed.

Set `PLOT_MONTH` to YYYYMM for one month, or `None` for the period mean weighted
by ERA5 calendar days (February 2012 has 29 days). Output caches are separated
from the old nudged-reference products in
`OUTPUT_ROOT/processed/monthly_vertical_bias/ERA5/`; figures use
`fig06_monthly_vertical_bias_ERA5_*.png` and explicitly label ERA5 as reference. Each
panel also shows RMSE and spatial PCC of the displayed experiment and ERA5 mean
sections, weighted by cosine latitude and pressure-layer thickness over paired
finite cells. The caches retain paired model/reference zonal means for these
statistics; older bias-only caches are recomputed automatically. Q RMSE uses g/kg.
These annotations describe mean cross sections, not daily RMSE or temporal PCC.

## Figure 3 variable checkpoints

Figure 3 saves each completed case/variable pair immediately under
`processed/spatial_error_reduction/by_variable/<case>/<variable>_<period>.nc`.
With `FORCE_COMPUTE=False`, adding an entry to `REQUIRED_VARIABLES` computes only
missing pairs. Rerun the notebook from its settings cell after changing selections;
rerunning the processing cell alone also resumes interrupted work. Plotting reloads
the selected checkpoints directly, without rebuilding a combined case file.

`REUSE_LEGACY_CACHES=True` copies compatible variables from the original combined
case files on the next run, leaving those files intact. Old caches lack complete
provenance; migration checks their recorded case/CTRL/REF identifiers and period
and requires the full-period dates and original RMSE threshold. Disable legacy
reuse if the original settings differed. New caches record per-variable settings;
changing dates, source paths, or the RMSE threshold invalidates affected pairs.
If source file contents change at the same paths, use `FORCE_COMPUTE=True`.

Inputs for a new pair must have identical timestamps and horizontal grids to REF,
so adding a case does not silently change the averaging period of cached results.
Finite-value QC runs during field loading to avoid a separate full-data read.
