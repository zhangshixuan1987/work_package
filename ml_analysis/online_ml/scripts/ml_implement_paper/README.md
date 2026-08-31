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
3. Open `notebooks/ml_online_bias_correction_analysis.ipynb` and run it from
   top to bottom.

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

The expected products are four NetCDF diagnostic files, four publication
figures, and `tables/summary_metrics.csv`. No analysis output should be written
inside this source directory.

## Pressure interpolation

If the data have a one-dimensional pressure coordinate on `lev`, the workflow
uses it directly. For native E3SM hybrid levels, it calculates pressure from
`hyam`, `hybm`, `P0`, and surface pressure. Update these names in the YAML file
when the input convention differs.

