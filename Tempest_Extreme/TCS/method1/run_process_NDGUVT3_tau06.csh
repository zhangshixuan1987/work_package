#!/bin/csh

#conda activate e3sm_analysis 
#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.csh

set command = "DetectNodes"
set zsfile   = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/data/e3sm_phis.nc"
set data_dir = "/global/cscratch1/sd/zhan391/DARPA_project/nudged_e3smv2_simulation"

set runnam   = "e3sm_v2_NDGERA5_PL_UVT3_tau06_F20TR_ne30pg2_EC30to60E2r2"
set expnam   = "NDGUVT3_tau06"

set varnam   = "TCS"
set syear    = 2007
set eyear    = 2017

#generate file list 
rm -f file_${expnam}_list; touch file_${expnam}_list
set iyear = $syear
while ( $iyear <= $eyear )
   echo ${data_dir}/${runnam}/run/*eam.h2*${iyear}*.nc >> file_${expnam}_list
 @ iyear++
end

#process the data 
set outfile  = "${expnam}_${varnam}_${syear}-${eyear}.txt"
rm -f $outfile; touch $outfile

set files = `cat file_${expnam}_list`
set nfile = $#files 
echo $nfile 

set i = 1
while ( $i <= $nfile ) 

  set file   = $files[$i]
  set finput = "input_${expnam}_${varnam}_list.txt"
  rm -f $finput; touch $finput

  echo  "$file;$zsfile" >> ${finput}

  rm out_${expnam}_${varnam}.dat

  $command --in_data_list $finput \
            --verbosity 0 \
            --in_connect "outCS_ne30pg2_connect.txt" \
            --closedcontourcmd "PSL,200.0,5.5,0;_DIFF(Z200,Z500),-11.0,6.5,1.0" \
            --mergedist 6.0 \
            --searchbymin "PSL" \
            --outputcmd "PSL,min,0;_VECMAG(UBOT,VBOT),max,2;_DIV(PHIS,9.81),min,0" \
            --timefilter "6hr" \
            --latname lat  \
            --lonname lon \
            --out out_${expnam}_${varnam}.dat

  cat out_${expnam}_${varnam}.dat >> $outfile
  rm out_${expnam}_${varnam}.dat

 @ i++
end 

exit
