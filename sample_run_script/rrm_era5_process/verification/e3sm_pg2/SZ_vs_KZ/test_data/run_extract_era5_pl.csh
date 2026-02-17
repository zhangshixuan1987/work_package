#!/bin/csh
set date = "2009-04-01_00:00:00"
set stim = '2009-04-01 00:00:00'
set etim = '2009-04-01 00:00:00'
set era5_dir = "/global/cfs/projectdirs/m3522/cmip6/ERA5"
set map_file = "/global/cscratch1/sd/zhan391/SciDAC_m3089/acme_init/era5_reanalysis/era5_to_e3sm_rrm_map/map_721x1440_ne30pg2.nc"
set out_file = "ERA5_analysis_pl_721x1440_${date}.nc" 

rm -rvf $out_file

set var = "U"
set file = "${era5_dir}/e5.oper.an.pl/200904/e5.oper.an.pl.128_131_u.ll025uv.2009040100_2009040123.nc"
ncks -v  $var -d time,"${stim}","${etim}" $file  $out_file

set var = "V"
set file = "${era5_dir}/e5.oper.an.pl/200904/e5.oper.an.pl.128_132_v.ll025uv.2009040100_2009040123.nc"
ncks -A -v  $var -d time,"${stim}","${etim}" $file  $out_file

set var = "T"
set file = "${era5_dir}/e5.oper.an.pl/200904/e5.oper.an.pl.128_130_t.ll025sc.2009040100_2009040123.nc"
ncks -A -v  $var -d time,"${stim}","${etim}" $file  $out_file

set var = "Q"
set file = "${era5_dir}/e5.oper.an.pl/200904/e5.oper.an.pl.128_133_q.ll025sc.2009040100_2009040123.nc"
ncks -A -v  $var -d time,"${stim}","${etim}" $file  $out_file

set var = "SP"
set file = "${era5_dir}/e5.oper.an.sfc/200904/e5.oper.an.sfc.128_134_sp.ll025sc.2009040100_2009043023.nc"
ncks -A -v  $var -d time,"${stim}","${etim}" $file  $out_file
ncrename  -v $var,PS $out_file
#ncrename  -d latitude,lat   $out_file
#ncrename  -d longitude,lon  $out_file
#ncrename  -d level,lev      $out_file
#ncrename  -v latitude,lat   $out_file
#ncrename  -v longitude,lon  $out_file
#ncrename  -v level,lev      $out_file
