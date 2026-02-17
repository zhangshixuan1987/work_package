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

set Outnam = (E3SMv2_NDGUV_tau6)
set ncase  = $#Outnam

set i = 1
while ( $i <= $ncase )
 
 set OUT_NAME  = $Outnam[$i]

 #location of model history file
 set RUN_FILE_DIR  = /global/cfs/cdirs/e3sm/www/zhan391/darpa_temporary_data_share/SE_PG2/before_nudging
 set OUT_DIR_SE    = /global/cfs/cdirs/e3sm/www/zhan391/darpa_temporary_data_share/SE_PG2/clim


 ### state variable for input of ML #######
 set year = 2007
 while ($year <= 2017) 

  set timstr  = "${year}01010000-${year}12312100"
  set SE_file = $RUN_FILE_DIR/${OUT_NAME}_3hourly_${timstr}.nc 
  echo $SE_file 

  rm -rvf ${OUT_DIR_SE}/time_${year}.nc
  ncks -v time $SE_file ${OUT_DIR_SE}/time_${year}.nc 

  @ year++ 
 end 

 @ i++
end 
