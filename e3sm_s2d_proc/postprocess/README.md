# S2D hindcast post-processing

Archive, transfer, and post-process the bundled hindcast output from
`../hindcast_run/`. Copied from the production workflow in
`/global/cfs/cdirs/e3sm/zhan391/s2d_diag` (`post_*` and `scripts/`).

```
postprocess/
  4DEnVarOcn.v0/    per-set drivers (from post_4DEnVarOcn)
  BruteForce/       per-set drivers (from post_bruteforce)
  JRA55_FOSIRL/     per-set drivers (from post_fosirl)
  RTENSO2026/       per-set drivers (copied from BruteForce/)
  scripts/          shared ncclimo / ncremap / e3sm_to_cmip templates
  zstash_S2D.sh     zstash create / update / check for one case
  zstash_extract.bash  restore cases from HPSS
```

The per-set directories hold the same drivers; they differ in `exp`
(`..._4DEnVarOcn`, `..._BruteForce`, `..._JRA55_FOSIRL`,
`..._RTENSO2026_<ic_method>`), source and output
roots, lead-date lists, and a few variable lists. Edit the settings block at the
top of each driver (`exp`, `leads`, `run_dir`, `out_dir`, `work_dir`,
`max_jobs`) before running. All drivers use E3SM Unified on `pm-cpu`.

`RTENSO2026/` covers the 2015 hindcasts from `../hindcast_run/RTENSO2026/`
(leads 2015020100, 2015050100, 2015080100, 2015110100; 10 members; 2 years
including the 1-month short spin-up). One set of drivers serves every IC
method: `IC_METHOD` (default `clim`) selects
`WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_RTENSO2026_${IC_METHOD}`, e.g.
`IC_METHOD=obs sbatch run_post_atm_mon1.bash`.
Before stage 1 moves the history files to CFS, set `RUN_DIR` to read them
from `$PSCRATCH`, e.g.
`RUN_DIR=$PSCRATCH/e3sm_project/E3SMv3_S2D IC_METHOD=obs sbatch run_post_atm_mon1.bash`.

Case directories follow the hindcast layout:
`${run_dir}/${exp}_${YYYYMMDDHH}/EN##/archive/<comp>/hist`.

## Stage 1: archive and transfer

Run in order for each batch of lead dates:

| Script | Purpose |
| --- | --- |
| `0_run_archive.bash` | Run `case.st_archive` for every member of the selected cases (`CASE_DATE` empty = all, `"1980*"` = one year, or one date); logs to `archive_logs/`. |
| `1_run_zstash_multi.bash` | Archive each case to HPSS with `../zstash_S2D.sh` (copied to `work_dir`). |
| `2_run_mov_archive.bash` | Move `archive/<comp>/hist` from `$PSCRATCH` to the CFS simulation directory with `rsync --remove-source-files`. |
| `3_run_check_clean.bash` | Compare source and destination history files with `cmp` and remove the scratch copies that match; failures go to `post_process_logs/`. |

`JRA55_FOSIRL/` has no `0_run_archive.bash`. `BruteForce/run_rsync.bash` (also in
`RTENSO2026/`) is an alternative copy from `$PSCRATCH` to the
post-processing root (`--dry-run` supported).

`zstash_S2D.sh` usage:

```bash
./zstash_S2D.sh [dry-run] {create|update|check} CASE_NAME
```

## Stage 2: time series and climatologies

Each `run_post_<comp>_<freq>.bash` loops over lead dates and members, renders
`scripts/${group}_${comp}_${freq}_${grid}.bash` into a per-member script, and
runs up to `max_jobs` at once. Output goes to `out_dir`
(`/global/cfs/cdirs/e3sm/S2S2D/post_process`).

| Component | Drivers | Products |
| --- | --- | --- |
| Atmosphere (EAM) | `atm_hour1`, `atm_hour3`, `atm_hour6`, `atm_day1`, `atm_mon1`, `atm_glb`, `atm_clim_mon`, `atm_clim_diurnal` | Time series regridded to 180x360 (`aave`), global means, monthly and diurnal (3-hourly) climatologies |
| Land (ELM) | `lnd_2d_day1`, `lnd_3d_day1`, `lnd_mon1`, `lnd_native_day1`, `lnd_native_mon1`, `lnd_clim_mon`, `lnd_clim_native` | Regridded and native-grid time series and climatologies |
| River (MOSART) | `rof_day1`, `rof_mon1`, `rof_native_day1`, `rof_native_mon1` | Regridded and native-grid time series |
| Ocean (MPAS-O) | `ocn_day1`, `ocn_mon1`, `ocn_rgn` | Regridded time series and regional means |
| Sea ice (MPAS-SI) | `ice_day1`, `ice_mon1`, `ice_rgn` | Regridded time series and regional means |

`length` in each driver is the number of forecast months to process, counted
from the lead month (1 = lead month only). File
names differ slightly between sets: 4DEnVarOcn uses `run_post_atm_month1.bash`
(others `run_post_atm_mon1.bash`) and keeps an older `run_post_lnd_day1.bash`;
JRA55-FOSIRL adds `run_post_atm_climo.bash`.
