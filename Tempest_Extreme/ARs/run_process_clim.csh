#!/bin/csh

#conda activate e3sm_analysis 
#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.csh

set command1 = "DetectBlobs"
set command2 = "NodeFileFilter"
set command3 = "VariableProcessor"

set data_dir = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes"
set zsfil    = "${data_dir}/e3sm_zs.nc"
set ftctrk   = "${data_dir}/TCS2/ERA5_TCS_Track_2007-2017.txt"

set runnam   = "E3SMv2_CLIM"
set expnam   = "CLIM"
set out_dir  = "./$expnam"


set varnam   = "ARS"
set syear    = 2007
set eyear    = 2017

if( ! -d $out_dir ) then
  mkdir -p $out_dir
endif

#process the data 
set iyear = $syear
while ( $iyear <= $eyear )

 set files = `echo ${data_dir}/data/${runnam}/${runnam}*${iyear}*.nc`
 set nfile = $#files 
 echo $nfile 

 set outfile  = $out_dir/${expnam}_${varnam}_${iyear}.nc
 rm -f $outfile

 set i = 1
 while ( $i <= $nfile )

  set file   = $files[$i]
  set finput = "input_${expnam}_${varnam}_list.txt"
  rm -f $finput; touch $finput

  set foutput = "output_${expnam}_${varnam}_list.txt"
  rm -f $foutput; touch $foutput

  set tmpstr = $file
  set tmpstr = "$tmpstr"
  echo "$tmpstr" >> ${finput}

  set foutput = "output_${expnam}_${varnam}_list.txt"
  rm -f $foutput; touch $foutput
  echo "$outfile" >> ${foutput}

  $command1 --in_data_list $finput \
            --verbosity 0 \
            --in_connect "outCS_ne30pg2_connect.txt" \
            --thresholdcmd "_LAPLACIAN{8,10}(_VECMAG(TUQ,TVQ)),<=,-20000,0" \
            --minabslat 15 \
            --geofiltercmd "area,>=,4e5km2" \
            --timefilter "6hr" \
            --latname lat  \
            --lonname lon \
            #--out_list $foutput \
            --out $outfile \
            --tagvar "AR_binary_tag"

  @ i++
 end
 @ iyear++
end

