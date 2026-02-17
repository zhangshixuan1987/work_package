#!/bin/csh

conda activate e3sm_analysis 

set runnam   = "e3sm_v2_NDGUVTQ_F20TR_ne30pg2_EC30to60E2r2"
set expnam   = "NDGUVTQ"
set feature  = "TCS"
set parset   = "set1"

set freq     = "6hourly"
set syear    = 2007
set eyear    = 2017

set top_dir  = "/global/cscratch1/sd/zhan391/DARPA_project"
set data_dir = "${top_dir}/post_process/data/model_output/${expnam}/${freq}"
set zsfil    = "${top_dir}/post_process/data/model_output/landmsk_coord/landmsk_coord.nc"
set out_dir  = "./filter_data"

set varList  = ("PRECT" "PSL" "Q850" "U850" "V850" "RH850" "Z200" "Z500" "TS" "TREFHT" "T500" "T850" "OMEGA500")
set nvars    = $#varList

if ( ! -d $out_dir ) then 
  mkdir -p $out_dir
endif 

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

   #filter data within r8 radius of storm center and composite
   set trk_file = "${out_dir}/${expnam}_${feature}_radprofs_${syear}-${eyear}.txt"
   set command3 = "NodeFileFilter"
   $command3  --in_nodefile  $trk_file \
              --in_fmt "lon,lat,rsize,rprof" \
              --in_data_list  infile_${parset}_${expnam}_list \
              --out_data_list filter_${parset}_${expnam}_list \
              --var "$varnam" \
              --bydist "rsize" \
#              --maskvar "mask"

 @ iv ++
end  
