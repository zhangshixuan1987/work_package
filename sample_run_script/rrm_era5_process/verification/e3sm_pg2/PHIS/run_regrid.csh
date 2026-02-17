#!/bin/csh

module load nco

set era5_grid        = "640x1280"
set e3sm_grid        = "northamericax4v1pg2"
set work_dir         = "/global/cscratch1/sd/zhan391/SciDAC_m3089/acme_init/era5_reanalysis"
set era5_data_dir    = "${work_dir}/era5_an_ml" 
set era5_to_e3sm_dir = "${work_dir}/era5_to_e3sm_rrm"
set map_file         = "${work_dir}/era5_to_e3sm_rrm_map/map_${e3sm_grid}_${era5_grid}.nc"

set time_tag         = "2009-04-01-00000"

set ndg_file = "${era5_to_e3sm_dir}/ERA5_ne30pg2_L72_${e3sm_grid}.${time_tag}.nc"
set rgd_file = "ERA5_${e3sm_grid}_${time_tag}_${era5_grid}.nc"

ncremap -m ${map_file} ${ndg_file} ${rgd_file}

