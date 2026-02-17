#!/bin/csh
# Submit this script as : sbatch ./[script-name]
#SBATCH -N 1
#SBATCH -C amd
#SBATCH -t 6:00:00
#SBATCH -q bigmem
#SBATCH -A m3089

module load nco

setenv OMP_NUM_THREADS 1

set workdir = `pwd`

cd $workdir

set Cases  = (e3sm_v2_NDGERA5_PL_UV3_tau06_F20TR_ne30pg2_EC30to60E2r2)
set Outnam = (ERA5)
set ncase  = $#Cases

set i = 1
while ( $i <= $ncase )
 
 set CASE_NAME = $Cases[$i]
 set OUT_NAME  = $Outnam[$i]
 echo $CASE_NAME

 #location of model history file
 set OUT_DIR_SE    = /global/cfs/cdirs/e3sm/www/zhan391/darpa_temporary_data_share/SE_PG2/reference_data
 set OUT_DIR_FV    = /global/cfs/cdirs/e3sm/www/zhan391/darpa_temporary_data_share/FV_180x360/reference_data

 #Mapping file
 set MAP_FILE      = /global/cscratch1/sd/zhan391/DARPA_project/post_process/map_file/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc

 if ( ! -d $OUT_DIR_SE ) then
   mkdir -p $OUT_DIR_SE
 endif

 if ( ! -d $OUT_DIR_FV ) then
   mkdir -p $OUT_DIR_FV
 endif

 ### state variable for input of ML #######
 set year = 2012
 while ($year <= 2012)

  set timstr  = "${year}01010000-${year}12311500"
  set SE_file = $OUT_DIR_SE/${OUT_NAME}_3hourly_${timstr}.nc
  set FV_file = $OUT_DIR_FV/ML_REF_${OUT_NAME}_3hourly_${timstr}.nc
  rm -rvf ${FV_file}
  ncremap -m $MAP_FILE -i $SE_file -o $FV_file

  @ year++
 end

 @ i++
end
