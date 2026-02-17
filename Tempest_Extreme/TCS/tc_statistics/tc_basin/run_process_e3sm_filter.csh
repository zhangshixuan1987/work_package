#!/bin/csh

#conda activate e3sm_analysis 
#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.csh

set command0 = "DetectNodes"
set command1 = "StitchNodes"
set command2 = "HistogramNodes"

set runnam   = "e3sm_v2_CLIM_F20TR_ne30pg2_EC30to60E2r2"
set expnam   = "CLIM"
set freq     = "6hourly"

set data_dir = "/global/cscratch1/sd/zhan391/DARPA_project/post_process/data/model_output/${expnam}/${freq}"
set zsfil    = "/global/cscratch1/sd/zhan391/DARPA_project/post_process/data/model_output/landmsk_coord/landmsk_coord.nc"

set fvList   = ("PSL" "UBOT" "VBOT" "Z200" "Z500" "T200" "T500")
set nfvs     = $#fvList 
#echo $nfvs

set varnam   = "ETCS"
set syear    = 2007
set eyear    = 2017
set parset   = "set1"

#process the data 
set tdn_file = ${expnam}_${parset}_${varnam}_${syear}-${eyear}.txt
set trk_file = "${expnam}_${parset}_${varnam}_NA_Track_${syear}-${eyear}.txt"
rm -f $trk_file; touch $trk_file

$command1  --in  $tdn_file \
           --out $trk_file \
           --in_fmt "lon,lat,slp,wind,zs" \
           --range ${trajRange} \
           --mintime ${trajMinLength} \
           --maxgap ${trajMaxGap}  \
           --min_endpoint_dist 12.0 

exit
