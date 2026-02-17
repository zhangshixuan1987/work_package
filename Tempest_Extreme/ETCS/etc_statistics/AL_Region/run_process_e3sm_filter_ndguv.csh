#!/bin/csh

conda activate e3sm_analysis 

set runnam   = "e3sm_v2_NDGUV_F20TR_ne30pg2_EC30to60E2r2"
set expnam   = "NDGUV"
set trknam   = "ppe_ndguv"
set feature  = "ETCS"
set parset   = "set1"

set freq     = "6hourly"
set syear    = 2007
set eyear    = 2017

set top_dir  = "/global/cscratch1/sd/zhan391/DARPA_project"
set data_dir = "${top_dir}/post_process/data/model_output/${expnam}/${freq}"
set zsfil    = "${top_dir}/post_process/data/model_output/landmsk_coord/landmsk_coord.nc"
set out_dir  = "./filter_data"
set trk_dir  = "${top_dir}/TempestExtremes/${feature}/script_FV/${trknam}"

set varList  = ("PRECT" "PSL" "Q850" "U850" "V850" "RH850" "Z200" "Z500" "TS" "TREFHT" "T500" "T850" "OMEGA500")
set nvars    = $#varList

set basin    = "AL"
set lats     = 24
set latn     = 52
set lonw     = 234
set lone     = 294

if ( ! -d $out_dir ) then 
  mkdir -p $out_dir
endif 

#process the etc track data
set command1 = "StitchNodes"
set tdn_file = ${trk_dir}/${expnam}_${parset}_${feature}_${syear}-${eyear}.txt
set trk_file = "${expnam}_${basin}_track_nf_${syear}-${eyear}.txt"
rm -f $trk_file; touch $trk_file
$command1  --in  $tdn_file \
           --out $trk_file \
           --in_fmt "lon,lat,slp,wind,zs" \
           --range   6.0   \
           --mintime "60h" \
           --maxgap  "24h" \
           --min_endpoint_dist 12.0 \
           --threshold "lat,>,${lats},1;lon,>,${lonw},1;lat,<,${latn},1;lon,<,${lone},1"

set iv = 1
while ( $iv <= $nvars ) 

  set varnam = $varList[$iv]

  #generate composite file list 
  rm -f infile_${parset}_${expnam}_list; touch infile_${parset}_${expnam}_list
  rm -f filter_${parset}_${expnam}_list; touch filter_${parset}_${expnam}_list

  set iyear  = $syear
  while ( $iyear <= $eyear )
     set inpfil = `ls ${data_dir}/*${varnam}*${iyear}*.nc`
     set outfil = $out_dir/${expnam}_${feature}_${varnam}_filter_${iyear}.nc
     echo $inpfil >> infile_${parset}_${expnam}_list
     echo $outfil >> filter_${parset}_${expnam}_list
   @ iyear++
  end

   set command2 = "NodeFileEditor"
   set trk_filter = "${expnam}_${basin}_track_wf_${syear}-${eyear}.txt"
   rm -f $trk_filter; touch $trk_filter
   $command2  --in_nodefile  $trk_file \
           --in_data_list infile_${parset}_${expnam}_list \
           --in_fmt "lon,lat,slp,wind,zs" \
           --out_fmt "lon,lat,slp,wind,zs" \
           --out_nodefile $trk_filter \
           --colfilter "slp,<=,99000."

   #filter data within 25deg of storm center and composite
   set command3 = "NodeFileFilter"
   $command3  --in_nodefile  $trk_file \
              --in_fmt "lon,lat,slp,wind,zs" \
              --in_data_list  infile_${parset}_${expnam}_list \
              --out_data_list filter_${parset}_${expnam}_list \
              --var "$varnam" \
              --bydist 35.0 \
              --maskvar "mask"

   set compfil = "$out_dir/${expnam}_all_${feature}_${varnam}_${syear}-${eyear}_composite.nc"
   set command4 = "NodeFileCompose"
   $command4  --in_nodefile  $trk_file \
              --in_fmt "lon,lat,slp,wind,zs" \
              --in_data_list  infile_${parset}_${expnam}_list \
              --out_data      "$compfil" \
              --var "$varnam" \
              --max_time_delta "2h" \
              --op "mean" \
              --dx 1.0 \
             --resx 80

   set compfil = "$out_dir/${expnam}_strong_${feature}_${varnam}_${syear}-${eyear}_composite.nc"
   $command4  --in_nodefile  $trk_filter \
              --in_fmt "lon,lat,slp,wind,zs" \
              --in_data_list  infile_${parset}_${expnam}_list \
              --out_data      "$compfil" \
              --var "$varnam" \
              --max_time_delta "2h" \
              --op "mean" \
              --dx 1.0 \
             --resx 80
 @ iv ++
end  
