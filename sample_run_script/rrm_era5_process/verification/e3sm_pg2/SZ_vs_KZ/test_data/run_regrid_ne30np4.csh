#!/bin/csh

module load nco
set script_path = /global/u1/z/zender/bin_cori
set path        = ( ${script_path}  ${path} )

set work_dir  = "/global/cscratch1/sd/zhan391/SciDAC_m3089/acme_init/era5_reanalysis"
set time_tag  = "2009-04-01-00000"
set stim      = '2009-04-01 00:00:00'
set etim      = '2009-04-01 00:00:00'

set era5_grid        = "721x1440"
set e3sm_grid        = "ne30np4"
set era5_to_e3sm_dir = "${work_dir}/era5_to_e3sm"
set map_file         = "${work_dir}/era5_to_e3sm_map/era5_pl/map_${e3sm_grid}_to_era5_conserve.20191206.nc"

set ndg_file =  "./era5_ne30L72_2009-04-01.nc" #"${era5_to_e3sm_dir}/ERA5_ne30np4_L72_${e3sm_grid}.${time_tag}.nc"
set rgd_file = "ZK_${e3sm_grid}_to_ERA5_${era5_grid}_${time_tag}.nc"

ncks -d time,"${stim}","${etim}" $ndg_file  tmp.nc

$script_path/ncremap -m ${map_file} tmp.nc ${rgd_file}

rm -rvf tmp.nc

