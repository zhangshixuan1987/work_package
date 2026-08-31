# External input data

Input datasets are not stored in this repository. Supply their directory to
each analysis script with `--data-dir`.

The current scripts expect these files directly inside the directory passed to
`--data-dir`:

- abrupt4CO2 analysis: `Anomalies_of_abrupt4CO2_in_FOM.pkl` and
  `cam_grid_info.nc`
- GFexp analysis: `qdp.PatchEXP.108.nc` and `SOM_Patches_TS.nc`

The exact external directory can differ between users and machines. Do not
commit symlinks or copies of large datasets here.
