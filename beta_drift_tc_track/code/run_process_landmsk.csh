#!/bin/csh
# Submit this script as : sbatch ./[script-name]
#SBATCH -A esmd
#SBATCH -p short
#SBATCH -t 2:00:00
#SBATCH -N 4
#SBATCH  --job-name=ncclimo_ctrl

module load nco 
set script_name = ncremap
set remap       = ${script_name}
set MAP_FILE    = "/global/cscratch1/sd/zhan391/DARPA_project/post_process/map_file/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc"

set expnam      = "e3sm_v2_CLIM_F20TR_ne30pg2_EC30to60E2r2"

set run_dir     = "/global/cscratch1/sd/zhan391/DARPA_project/nudged_e3smv2_simulation/${expnam}/run"
set out_dir     = "/global/cscratch1/sd/zhan391/DARPA_project/Track_model/data"

set syear       = 2007
set eyear       = 2017

#process land mask file 
set files = "${run_dir}/${expnam}.eam.h0.2010-01.nc"

if ( ! -d ${out_dir} ) then 
  mkdir -p ${out_dir}
endif 

set varlist = ( "LANDFRAC"  "OCNFRAC" ) 
set len = $#varlist

foreach i (`seq 1 1 $len`)
 set outfile = ${out_dir}/$varlist[$i]_E3SM_180x360.nc
 rm -rvf msk_out.nc
 ncwa -a time -v $varlist[$i] $files msk_out.nc
 $remap -m $MAP_FILE -i msk_out.nc -o $outfile
 rm -rvf msk_out.nc
end
