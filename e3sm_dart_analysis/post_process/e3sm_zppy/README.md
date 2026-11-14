# e3sm_zppy templates

zppy-style post-processing scripts for the `F20TR_ne30pg2_r05_IcoswISC30E3r5_<CONFIG>_<YYYYMMDDHH>/<ENxx>`
hindcasts, one per product (`<kind>.bash`). They are the zppy scripts that made the `post/` tree of the
existing configurations (taken from `CLIM_BASE_2011121600/EN01/post/scripts`), with the run-specific
header values replaced by `@ROOT@ @CASE@ @ENS@ @YR_SRT@ @YR_END@ @MTH_SRT@ @MTH_END@` and MPI
parallelism switched to serial (`--parallel=serial --jobs=1`, `ncremap -p serial`).

They are run by `../run_process_e3sm_post.sh` (one product of one member) and
`../submit_process_e3sm_post.sbatch` (many cases and members).

- Monthly products (`ts_*_monthly_*`): unchanged zppy commands. Filled in, they reproduce the scripts of
  all 4 existing configurations x 8 starts exactly (apart from the two changes above; checked 2026-11-14).
- Daily products (`ts_{atm_daily,atm_6hourly,lnd_daily,rof_daily}_180x360_aave`, `ts_{lnd,rof}_daily_native`):
  zppy's `ncclimo --split` reads only the first few history files of these 2-3 month runs (atm stops
  after 4 days, lnd/rof at the end of the start month), so these templates instead concatenate every
  history file of the stream (`ncrcat`), remap once (180x360 only) and write one file per variable,
  `<VAR>_<YYYYMM>_<YYYYMM>.nc` named by the data's time range. Header, variables, maps and output
  directories are zppy's. The 180x360 ones also write `post/scripts/regen_ts_<stream>.status`
  (`OK files=... records=... range=...`), which the hindcast notebooks check. This replaces the former
  `post_process/regen_f20tr_ts.sh`: identical files for CLIM_BASE_2011123100/EN01 (atm, rof bit-identical;
  2 lnd fields differ by at most 1 float32 ulp from threaded remapping).
- Climatologies (`clim_*`): kept but not run; `ncclimo` needs whole years of monthly files.
