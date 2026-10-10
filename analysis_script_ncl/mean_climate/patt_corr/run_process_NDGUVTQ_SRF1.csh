#!/bin/csh
#SBATCH -A e3sm
#SBATCH -q slurm
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH  --job-name=nudge

set BASE_DIR = "/pscratch/sd/z/zhan391/DARPA_project/post_process/data/model_output"
set OBS_DIR  = "/global/cfs/cdirs/e3sm/zhan391/data"

set models   = ( "NDGUVTQ_SRF1" )
set lenm     = $#models

#set var      = ("U" "V" "T" "Q" "OMEGA" "OLR" "OLR" "PRECT" "PRECT" "PRECT")
#set lenv     = $#var

#set obdat    = ("ERA5" "ERA5" "ERA5" "ERA5" "ERA5" "NOAA_AVHRR" "NOAA_HIRS" "ERA5" "GPCP" "TRMM")
#set leno     = $#obdat

set var      = ("OLR" "OLR" "PRECT" "PRECT" "PRECT")
set lenv     = $#var

set obdat    = ("NOAA_AVHRR" "NOAA_HIRS" "ERA5" "GPCP" "TRMM")
set leno     = $#obdat

set plevs    = (850 700 600 500 400 300 200 100)
set lenp     = $#plevs

set outdir = "../fig_data" 
if( ! -d $outdir ) then
  mkdir -p  $outdir
endif

set yst  = 2007
set yed  = 2017
set nyr  = 0
set yyyy = 0

@ nyr =  $yed - $yst  + 1

foreach i (`seq 1 1 $lenm`)

 foreach j (`seq 1 1 $nyr `) 

  @ yyyy = $yst + $j - 1

  setenv year       $yyyy
  setenv MODEL_NAME $models[$i]
  setenv INPUT      $BASE_DIR
  setenv OUTPUT     $outdir
  setenv OBSDIR     $OBS_DIR
  setenv model      $models[$i]

  foreach j (`seq 1 1 $lenv`)

    setenv OBS_NAME $obdat[$j]
    setenv VAR      $var[$j]
    
    if ( $VAR == "OLR" ) then
      setenv PLEV 1000
      ncl make_sac.ncl & 
      ncl make_sac_3regions.ncl & 
      ncl make_tac.ncl
    else if ( $VAR == "PRECT" ) then
      setenv PLEV 1000
      ncl make_sac.ncl & 
      ncl make_sac_3regions.ncl & 
      ncl make_tac.ncl
    else
      foreach k (`seq 1 1 $lenp`)
        setenv PLEV $plevs[$k]
        ncl make_sac.ncl & 
        ncl make_sac_3regions.ncl &
        ncl make_tac.ncl 
      end
    endif

  end 

 end 

 echo "Calculate ACC done..."

end 

