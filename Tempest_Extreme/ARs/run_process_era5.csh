#!/bin/csh

#conda activate e3sm_analysis 
#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.csh

set command1 = "DetectBlobs"
set command2 = "NodeFileFilter"
set command3 = "VariableProcessor"

set expnam   = "ERA5"
set data_dir = "/global/cfs/projectdirs/m3522/cmip6/${expnam}"
set zsfil    = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/data/era5_zs.nc"
set ftctrk   = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/TCS2/ERA5_TCS_Track_2007-2017.txt"
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

set im = 1
while ( $im <= 12 ) 
 
 set mstr  = `printf "%02d" $im` 

 set files = `echo ${data_dir}/e5.oper.an.vinteg/${iyear}${mstr}/e5.oper.an.vinteg.*_viwve.*.nc`
 set nfile = $#files 

 set vfils = `echo ${data_dir}/e5.oper.an.vinteg/${iyear}${mstr}/e5.oper.an.vinteg.*_viwvn.*.nc`
 set nvfil = $#vfils

 echo $nfile $nvfil

 set outfile  = $out_dir/${expnam}_${varnam}_${iyear}-${mstr}.nc
 rm -f $outfile

 set i = 1
 while ( $i <= $nfile ) 

  set file   = $files[$i]
  set finput = "input_${expnam}_${varnam}_list.txt"
  rm -f $finput; touch $finput

  set foutput = "output_${expnam}_${varnam}_list.txt"
  rm -f $foutput; touch $foutput

  set tmpstr = $file
  set tmpstr = "$tmpstr;$vfils[$i]"
  echo "$tmpstr" >> ${finput}

  set foutput = "output_${expnam}_${varnam}_list.txt"
  rm -f $foutput; touch $foutput
  echo "$outfile" >> ${foutput}

  $command1 --in_data_list $finput \
            --verbosity 0 \
            --thresholdcmd "_LAPLACIAN{8,10}(_VECMAG(VIWVE,VIWVN)),<=,-20000,0" \
            --minabslat 15 \
            --geofiltercmd "area,>=,4e5km2" \
            --timefilter "6hr" \
            --latname latitude  \
            --lonname longitude \
            #--out_list $foutput \
            --out $outfile \
            --tagvar "AR_binary_tag"

  @ i++
 end 
 @ im ++ 
end 
 @ iyear++
end 

#set finput = "input_${expnam}_${varnam}_list.txt"
#rm -f $finput; touch $finput
#ls $out_dir/*.nc >> $finput

#set ar_nff_file = ${expnam}_${varnam}_NFF_${syear}-${eyear}.txt
#rm -f $ar_nff_file; touch $ar_nff_file
#$command2 --in_nodefile ERA5_TC_tracks.txt \
#          --in_fmt "lon,lat,slp,wind,zs" \
#          --in_data_list $finput \
#          --out_data_list $ar_nff_file \
#          --var "binary_tag" \
#          --bydist 5.0 \
#          --invert \
#          --var "TC_binary_tag"

exit
