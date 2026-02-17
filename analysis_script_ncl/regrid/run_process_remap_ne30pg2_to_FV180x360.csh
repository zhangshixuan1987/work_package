#!/bin/csh
# Submit this script as : sbatch ./[script-name]
#SBATCH -N 1
#SBATCH -C amd
#SBATCH -t 6:00:00
#SBATCH -q bigmem
#SBATCH -A m3089

module load nco
module load cdo 

setenv OMP_NUM_THREADS 1

set workdir = `pwd`

cd $workdir

set Expnam = (E3SMv2_CLIM)
set ncase  = $#Expnam

#mapping file for horizonal regriding 
set MAP_FILE      = "./map_ne30pg2_to_cmip6_180x360_aave.20200201.nc" 

#location of model history file
set RUN_FILE_DIR  = "/global/cfs/cdirs/e3sm/www/zhan391/darpa_temporary_data_share/SE_PG2/clim"
set OUT_DIR       = "./"

# only extract a specific period 
set filname       = "E3SMv2_CLIM_3hourly_201301010000-201312311500.nc"
set start_date    = "2013-01-01 00:00:0.0"
set end_date      = "2013-01-01 18:00:0.0"
set tim_str       = "2013010100-2013010118"

set i = 1
while ( $i <= $ncase )
 
  set CASE    = $Expnam[$i]
  set outfile = "$OUT_DIR/$CASE.$tim_str.nc"

  set infiles = $RUN_FILE_DIR/$filname

  echo $infiles
  rm -rvf SE_NE30PG2.nc
  ncks -d time,"$start_date","${end_date}" $infiles SE_NE30PG2.nc

  rm -rvf ${outfile}
  ncremap -m $MAP_FILE -i SE_NE30PG2.nc -o $outfile
  rm -rvf SE_NE30PG2

 @ i++
end 
