#!/bin/csh
# Submit this script as : sbatch ./[script-name]
#SBATCH -N 1
#SBATCH -C amd
#SBATCH -t 6:00:00
#SBATCH -q bigmem
#SBATCH -A m3089

#module load nco
#module load cdo 
#module swap nco/5.0.1-gcc8 

setenv OMP_NUM_THREADS 1

set workdir = `pwd`

cd $workdir
set datadir = "/pscratch/sd/z/zhan391/darpa_scratch"
set Runnam  = ( "FREE_SANDY2012_F20TR_ne30pg2_EC30to60E2r2_pm-cpu" ) 
set Expnam  = ( "CTRL" ) 
set ncase   = $#Expnam

#mapping file for horizonal regriding 
set MAP_FILE   = "map_ne30pg2_to_cmip6_180x360_aave.20200201.nc"
set OUT_DIR    = "./"

# extract data  sandy 
set start_date = "2012-10-21 00:00:0.0"
set end_date   = "2012-11-03 00:00:0.0"
set tim_str    = "2012102100-2012110300"
set var_List   = ( "PRECT"     "PRECL"    "PRECC"   "CLDLOW"  "CLDMED"  "CLDHGH"  "CLDTOT"  \
                   "TGCLDLWP"  "TGCLDIWP" "CAPE"    "CIN"     "TAUX"    "TAUY"    \
                   "SHFLX"     "LHFLX"    "FLUT"    "FLUTC"   "SWCF"    "LWCF"    ) 
set nvars  = $#var_List
set eam_hf = "eam.h1"

set i = $ncase 
while ( $i <= $ncase )
 
  set EXPS    = $Runnam[$i]
  set CASE    = $Expnam[$i]

  #location of model history file
  set RUN_FILE_DIR = "${datadir}/${EXPS}/run" 
  set outdir       = "$OUT_DIR/$CASE"

  if ( ! -d $outdir ) then
    mkdir -p $outdir
  endif

  set j = 1
  while ( $j <= $nvars)
    set var = $var_List[$j]
    echo $var
    set routfile = "$outdir/${CASE}.${var}.${tim_str}.nc"
    set infiles  = "$RUN_FILE_DIR/${EXPS}*${eam_hf}*.nc"
    echo $infiles 

    if ( ! -d $outdir/SE ) then
     mkdir -p $outdir/SE
    endif

    set outfile = "$outdir/SE/${EXPS}.${var}.${tim_str}.nc"
    rm -rvf $outfile
    ncrcat -d time,"$start_date","${end_date}" -v $var $infiles $outfile
    rm -rvf $routfile
    ncremap -m $MAP_FILE -i $outfile  -o $routfile
  @ j++
 end   
 @ i++
end 
