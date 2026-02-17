#!/bin/csh
#SBATCH -N 1
#SBATCH -C amd
#SBATCH -t 30:00
#SBATCH -q bigmem
#SBATCH -A m3089

# Here the StitchNodes thresholds require that 
# (1) storms persist for at least 60 h with 
# (2) a maximum gap (time between sequential detections satisfying 
#                    the DetectNodes criteria) of at most 1-h. 
# (3) ETCs move at least 12-degree GCD from the start to the end of
#     the trajectory (min_endpoint_dist), in order to eliminate stationary 
#     features (e.g., the Icelandic Low) and spurious shallow lows 
#     generated over regions of high topography.
# (4) Optional: at least one point must pass through a geographic region
#     Global: 90S-90N, 0-360   longitude
#     CONUS:  24-52N, 234-294E longitude

conda activate e3sm_analysis 
#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-knl.csh
#conda activate e3sm_unified_1.6.0

set command0 = "StitchNodes"
set command1 = "HistogramNodes"
set varnam   = "ETCS"
set syear    = 2007
set eyear    = 2017

set explist  = ("ERA5" "CLIM" "NDGUV3_tau06" "NDGUVT3_tau06" "NDGUVTQ3_tau06")
set nexp     = $#explist

set connect  = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/grid_info/outCS_ne30pg2_connect.txt"

foreach expnam ($explist) 

  #process the data 
  set tdn_file = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/${varnam}/${expnam}_${varnam}_${syear}-${eyear}.txt"
  set trk_file = "${expnam}_${varnam}_Track_${syear}-${eyear}.txt"
  set hst_file = "${expnam}_${varnam}_hist_gram_${syear}-${eyear}.nc"
  rm -f $trk_file; touch $trk_file

  if ( $expnam == "ERA5" ) then 
    $command0  --in  $tdn_file \
               --out $trk_file \
               --in_fmt "lon,lat,slp,wind,zs" \
               --range 6.0 \
               --mintime "24h"\
               --maxgap "6"  \
               --min_endpoint_dist 12.0

    $command1  --in $trk_file \
               --iloncol 3 \
               --ilatcol 4 \
               --out $hst_file

  else

    $command0  --in  "$tdn_file" \
              --out "$trk_file" \
              --in_connect "${connect}" \
              --in_fmt "lon,lat,slp,wind,zs" \
              --range 6.0 \
              --mintime "24h" \
              --maxgap "6" \
              --min_endpoint_dist 12.0 

    $command1  --in $trk_file \
               --iloncol 2 \
               --ilatcol 3 \
               --out $hst_file

  endif 

end 

exit

