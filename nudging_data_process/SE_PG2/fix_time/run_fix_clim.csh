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

set Outnam = (E3SMv2_CLIM) 
set ncase  = $#Outnam

set i = 1
while ( $i <= $ncase )
 
 set OUT_NAME  = $Outnam[$i]

 #location of model history file
 set RUN_FILE_DIR  = /global/cfs/cdirs/e3sm/www/zhan391/darpa_temporary_data_share/SE_PG2/clim_old
 set OUT_DIR_SE    = /global/cfs/cdirs/e3sm/www/zhan391/darpa_temporary_data_share/SE_PG2/clim
 set OUT_DIR_FV    = /global/cfs/cdirs/e3sm/www/zhan391/darpa_temporary_data_share/FV_180x360/clim

 if ( ! -d $OUT_DIR_SE ) then
   mkdir -p $OUT_DIR_SE
 endif

 if ( ! -d $OUT_DIR_FV ) then
   mkdir -p $OUT_DIR_FV
 endif

 rm -rvf ${OUT_DIR_SE}/tmp_clim_${Outnam}_*.nc

 ### state variable for input of ML #######
 set year = 2007
 while ($year <= 2017) 

  set timstr      = "${year}01010000-${year}12312100"
  set SE_file_old = $RUN_FILE_DIR/${OUT_NAME}_3hourly_${timstr}.nc 
  set SE_file     = $OUT_DIR_SE/${OUT_NAME}_3hourly_${timstr}.nc
  echo $SE_file_old
  echo $SE_file
  rm -rvf $SE_file
  ncks -x -v PSL,TS,TREFHT,ZBOT,UBOT,VBOT,QBOT ${SE_file_old} ${SE_file}
  ncks -C -A -v time ./time_${year}.nc ${SE_file} 

 #set FV_file = $OUT_DIR_FV/${OUT_NAME}_3hourly_${timstr}.nc
 #rm -rvf ${FV_file}
 #ncremap -m $MAP_FILE -i $SE_file -o $FV_file

  @ year++ 
 end 

 @ i++
end 
