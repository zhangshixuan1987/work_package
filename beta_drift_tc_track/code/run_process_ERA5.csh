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
set outnam      = "CLIM"

set run_dir     = "/global/cscratch1/sd/zhan391/DARPA_project/nudged_e3smv2_simulation/${expnam}/run"
set out_dir     = "/global/cscratch1/sd/zhan391/DARPA_project/Track_model/data"

set syear       = 2007
set eyear       = 2017

set year = $syear
while ( $year <= $eyear ) 

 #process 3-hour data 
 setenv hfile "h2"
 set files =
 @ syearm1 = $year - 1
 @ eyearp1 = $year + 1 
 set iyear = $syearm1
 while ( $iyear <= $eyearp1 ) 
  echo $iyear 
  set stim  = `printf "%04d" $iyear`
  set ftmp  = `ls ${run_dir}/${expnam}*eam.${hfile}*${stim}*.nc`
  set files = (${files} ${ftmp})
  @ iyear++ 
 end

 if ( ! -d ${out_dir}/${outnam} ) then 
   mkdir -p ${out_dir}/${outnam}
 endif 

 set varlist = ("U850"  "V850" "U200" "V200")
 set len = $#varlist

 foreach i (`seq 1 1 $len`)
  set outfile = ${out_dir}/${outnam}/$varlist[$i]_${year}.nc
  ncrcat -d time,"${year}-01-01 00:00:0.0","${year}-12-31 21:00:0.0" -v $varlist[$i] $files tmp_out_${outnam}.nc
  $remap -m $MAP_FILE -i tmp_out_${outnam}.nc -o $outfile
  rm -rvf tmp_out_${outnam}.nc
 end
 
 @ year++ 
end 
