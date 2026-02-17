#!/bin/csh

conda activate e3sm_analysis 

set runnam   = "ERA5"
set expnam   = "ERA5"
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
set vnmList  = ("PRECT" "PSL" "Q850(850hPa)" "U850(850hPa)" "V850(850hPa)" "RH850(850hPa)" "Z200(200hPa)" "Z500(500hPa)" "TS" "TREFHT" "T500(500hPa)" "T850(850hPa)" "OMEGA500(500hPa)")
set nvars    = $#varList

if ( ! -d $out_dir ) then 
  mkdir -p $out_dir
endif 

set iv = 3
while ( $iv <= $nvars ) 

  set varnam = $varList[$iv]
  set varin  = $vnmList[$iv]

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
              --var "$varin" \
              --bydist "rsize" \
#              --maskvar "mask"

 @ iv ++
end  
