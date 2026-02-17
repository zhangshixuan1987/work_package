#!/bin/csh

module load nco
set script_path = /global/u1/z/zender/bin_cori
set path        = ( ${script_path}  ${path} )

set work_dir         = "/global/cscratch1/sd/zhan391/SciDAC_m3089/acme_init/era5_reanalysis"
set time_tag         = "2009-04-01-00000"

set era5_grid        = "721x1440"
set e3sm_grid        = "ne30pg2"
set era5_to_e3sm_dir = "${work_dir}/era5_pl_to_e3sm"
set map_file         = "${work_dir}/era5_to_e3sm_map/era5_pl/map_ne30pg2_to_era5_conserve.20210331.nc"

set ndg_file = "${era5_to_e3sm_dir}/ERA5_${e3sm_grid}_L72.${time_tag}.nc"
set rgd_file = "SZ_${e3sm_grid}_to_ERA5_${era5_grid}_${time_tag}.nc"

$script_path/ncremap -m ${map_file} ${ndg_file} ${rgd_file}

