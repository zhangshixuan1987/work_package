#!/bin/csh
#SBATCH -N 1
#SBATCH -C amd
#SBATCH -t 30:00
#SBATCH -q bigmem
#SBATCH -A m3089

#conda activate e3sm_analysis 
#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.csh

set command0 = "StitchNodes"
set command1 = "HistogramNodes"
set varnam   = "ETCS"
set syear    = 2007
set eyear    = 2017

set dirlist  = ("ppe_clim" "ppe_ndguv"     "ppe_ndguvt"     "ppe_ndguvtq")
set explist  = ("CLIM"     "NDGUV3_tau06"  "NDGUVT3_tau06"  "NDGUVTQ3_tau06")
set nexp     = $#explist

set parlist  = ( "set1" )
set npar     = $#parlist

set grdfil   = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/grid_info/outCS_ne30pg2_connect.txt"

set trajRange     = 9.0    # unit: degree Maximum travel distance for cyclone in 6 h [4,12]degree 
set trajMinLength = "24h"  # unit: hour, Minimum cyclone lifetime [12,72]h
set trajMaxGap    = "1"    # unit: hour, Maximum allowable gap in trajectory [0,30]h  

set i = 1
while ($i <= $nexp ) 

 set j = 1
 while ($j <= $npar ) 

  set expnam = $explist[$i]
  set dirnam = $dirlist[$i]
  set parset = $parlist[$j]
  echo $expnam  $parset 
 
  #post-process the track data 
  set tdn_file = "../${dirnam}/${expnam}_${parset}_${varnam}_${syear}-${eyear}.txt"
  set trk_file = "${expnam}_${varnam}_${parset}_Track_${syear}-${eyear}.txt"
  set hst_file = "${expnam}_${varnam}_${parset}_hist_gram_${syear}-${eyear}.nc"
  rm -f $trk_file; touch $trk_file

  if ( $expnam == "ERA5" ) then 

    $command0  --in  $tdn_file \
               --out $trk_file \
               --in_fmt "lon,lat,slp,wind,zs" \
               --range ${trajRange} \
               --mintime ${trajMinLength} \
               --maxgap ${trajMaxGap}  \
               --min_endpoint_dist 12.0

    $command1  --in $trk_file \
               --iloncol 3 \
               --ilatcol 4 \
               --out $hst_file

  else

    $command0  --in  $tdn_file \
               --out $trk_file \
               --in_connect ${grdfil} \
               --in_fmt "lon,lat,slp,wind,zs" \
               --range ${trajRange} \
               --mintime ${trajMinLength} \
               --maxgap ${trajMaxGap}  \
               --min_endpoint_dist 12.0

    $command1  --in $trk_file \
               --iloncol 2 \
               --ilatcol 3 \
               --out $hst_file

  endif 

  sleep 1

  @ j++ 
 end 

 @ i++
end 

exit

