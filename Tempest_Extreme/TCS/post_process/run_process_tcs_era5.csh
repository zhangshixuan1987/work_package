#!/bin/csh
#SBATCH -N 1
#SBATCH -C amd
#SBATCH -t 30:00
#SBATCH -q bigmem
#SBATCH -A m3089

set command0 = "StitchNodes"
set command1 = "HistogramNodes"
set varnam   = "TCS"
set syear    = 2007
set eyear    = 2017
set parset   = "set4"

set dirlist  = ("ppe_era5")
set explist  = ("ERA5")
set nexp     = $#explist

set connect  = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/grid_info/outCS_ne30pg2_connect.txt"

set wind_thresh  = 10.0 
set lat_thresh   = 50.0 
set zs_thresh    = 15.0 
set range_thresh = 8.0 
set time_thresh  = "10"
set time_gap     = "3"
set time_int     = 10

set i = 1
while ($i <= $nexp ) 

  set expnam = $explist[$i]
  set dirnam = $dirlist[$i]

  #process the data 
  set tdn_file = "../${dirnam}/${expnam}_${parset}_${varnam}_${syear}-${eyear}.txt"
  set trk_file = "${expnam}_${varnam}_${parset}_Track_${syear}-${eyear}.txt"
  set hst_file = "${expnam}_${varnam}_${parset}_hist_gram_${syear}-${eyear}.nc"
  rm -f $trk_file; touch $trk_file

  if ( $expnam == "ERA5" ) then 
    $command0  --in  $tdn_file \
               --out $trk_file \
               --in_fmt "lon,lat,slp,wind,zs" \
               --range ${range_thresh} \
               --mintime ${time_thresh} \
               --maxgap ${time_gap}  \
               --threshold "wind,>=,${wind_thresh},${time_int};lat,<=,${lat_thresh},${time_int};lat,>=,-${lat_thresh},${time_int};zs,<=,${zs_thresh},${time_int}" 

    $command1  --in $trk_file \
               --iloncol 3 \
               --ilatcol 4 \
               --out $hst_file

  else

    $command0  --in  "$tdn_file" \
               --out "$trk_file" \
               --in_connect "${connect}" \
               --in_fmt "lon,lat,slp,wind,zs" \
               --range ${range_thresh} \
               --mintime ${time_thresh} \
               --maxgap ${time_gap}  \
               --threshold "wind,>=,${wind_thresh},${time_int};lat,<=,${lat_thresh},${time_int};lat,>=,-${lat_thresh},${time_int};zs,<=,${zs_thresh},${time_int}"

    $command1  --in $trk_file \
               --iloncol 2 \
               --ilatcol 3 \
               --out $hst_file
  endif 

  sleep 1

 @ i++
end 

exit

