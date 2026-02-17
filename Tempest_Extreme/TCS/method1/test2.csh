#!/bin/csh

#conda activate e3sm_analysis 
source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.csh

set command1 = "DetectNodes"
set command2 = "StitchNodes"
set data_dir = "/global/cscratch1/sd/zhan391/DARPA_project/nudged_e3smv2_simulation"
set runnam   = "e3sm_v2_CLIM_F20TR_ne30pg2_EC30to60E2r2"
set expnam   = "CLIM"

set varnam   = "TC"
set syear    = 2007
set eyear    = 2017

#process the data 
set outfile  = ${expnam}_TCS_${syear}-${eyear}.txt
set trk_file = ${expnam}_TCS_Track_${syear}-${eyear}.txt
rm -f $trk_file; touch $trk_file
$command2 --in_list $outfile \
          --out $trk_file \
          --in_fmt "lon,lat,slp,wind,zs"\
          --range 8.0 --mintime "54h"\
          --maxgap "24h"\
          --threshold "wind,>=,10.0,10;lat,<=,50.0,10;lat,>=,-50.0,10;zs,<=,150.0,10"

