#!/bin/csh
#SBATCH -N 1
#SBATCH -C amd
#SBATCH -t 30:00
#SBATCH -q bigmem
#SBATCH -A m3089

# Here the StitchNodes thresholds require that 
# (1) storms persist for at least 10 h (mintime, i.e. the minimum 
#     persistence time of each trajectory, calculated as the time 
#     between initiation and termination) 
# (2) a maximum distance of 8-degree (in CGD) for a feature that can 
#     move between subsequent detections
# (3) a maximum gap (time between sequential detections satisfying 
#                    the DetectNodes criteria) of at most 3-h. 
# (4) field-dependent thresholds: 
#     a. "wind,>=,10.0,10;" : the wind magnitude (derived from the 
#        "wind" column in the nodefile. this ensures that these 
#        features are sufficiently intense to be classified as tropical storms
#     b, "lat,<=,50.0,10;lat,>=,-50.0,10;": the latitude of the feature must be 
#         between 50S and 50N for at least 10 time slices, so as to eliminate 
#         any extratropical features that could not have existed as 
#         tropical storms
#     c. "zs,<=,150.0,10": the feature must exist at an elevation below 150 m 
#         for at least 10 time slices; this removes false alarms that can 
#         often appear in regions of rough topography that are associated 
#         with the sea level pressure correction.

conda activate e3sm_analysis 
#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-knl.csh
#conda activate e3sm_unified_1.6.0

set command0 = "StitchNodes"
set command1 = "HistogramNodes"
set varnam   = "TCS"
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
               --range 8.0 \
               --mintime "10"\
               --maxgap "3"  \
               --threshold "wind,>=,10.0,10;lat,<=,50.0,10;lat,>=,-50.0,10;zs,<=,150.0,10" 

    $command1  --in $trk_file \
               --iloncol 3 \
               --ilatcol 4 \
               --out $hst_file

  else

    $command0  --in  "$tdn_file" \
              --out "$trk_file" \
              --in_connect "${connect}" \
              --in_fmt "lon,lat,slp,wind,zs" \
              --range 8.0 \
              --mintime "10" \
              --maxgap "3" \
              --threshold "wind,>=,10.0,10;lat,<=,50.0,10;lat,>=,-50.0,10;zs,<=,150.0,10" \
              --allow_repeated_times 

    $command1  --in $trk_file \
               --iloncol 2 \
               --ilatcol 3 \
               --out $hst_file

  endif 

end 

exit

