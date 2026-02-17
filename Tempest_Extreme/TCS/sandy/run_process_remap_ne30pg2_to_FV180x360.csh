#!/bin/csh
# Submit this script as : sbatch ./[script-name]
#SBATCH -N 1
#SBATCH -C amd
#SBATCH -t 6:00:00
#SBATCH -q bigmem
#SBATCH -A m3089

module load nco
module load cdo 
module swap nco/5.0.1-gcc8 

setenv OMP_NUM_THREADS 1

set workdir = `pwd`

cd $workdir

set Expnam = ("ERA5" "E3SMv2_CLIM" "E3SMv2_NDGUV_tau6" "E3SMv2_NDGUVT_tau6" "E3SMv2_NDGUVTQ_tau6")
set Config = ("reference_data" "clim" "before_nudging" "before_nudging" "before_nudging")
set ncase  = $#Expnam

#mapping file for horizonal regriding 
set MAP_FILE      = "../../regrid/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc"

#location of model history file
set RUN_FILE_DIR  = "/global/cfs/cdirs/e3sm/www/zhan391/darpa_temporary_data_share/SE_PG2"
set OUT_DIR       = "./"

# extract data for sandy 
set year          = 2012
set mon           = 10
set dd1           = 20
set dd2           = 31
set start_date    = "${year}-${mon}-${dd1} 00:00:0.0"
set end_date      = "${year}-${mon}-${dd2} 18:00:0.0"
set tim_str       = "${year}${mon}${dd1}00-${year}${mon}${dd2}18"

set i = 1
while ( $i <= 1 ) #$ncase )
 
  set CASE    = $Expnam[$i]
  set Key     = $Config[$i]

  set filname = "${CASE}_3hourly_${year}01010000-${year}12311500.nc"
  set infiles = $RUN_FILE_DIR/${Key}/$filname
  echo $infiles

  set outfile = "$OUT_DIR/$CASE.$tim_str.nc"

  rm -rvf SE_NE30PG2.nc
  ncks -d time,"$start_date","${end_date}" $infiles SE_NE30PG2.nc

  rm -rvf ${outfile}
  ncremap -m $MAP_FILE -i SE_NE30PG2.nc -o $outfile
  rm -rvf SE_NE30PG2
   
 @ i++
end 
