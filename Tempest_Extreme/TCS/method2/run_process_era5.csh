#!/bin/csh

#conda activate e3sm_analysis 
#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.csh

set command1 = "DetectNodes"
set command2 = "StitchNodes"

set expnam   = "ERA5"
set data_dir = "/global/cfs/projectdirs/m3522/cmip6/${expnam}"
set zsfil    = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/data/era5_zs.nc"

set varnam   = "TCS"
set syear    = 2007
set eyear    = 2017

#process the data 
set outfile  = ${expnam}_${varnam}_${syear}-${eyear}.txt
rm -f $outfile; touch $outfile

set iyear = $syear
while ( $iyear <= $eyear )

set im = 1
while ( $im <= 12 ) 
 
 set mstr  = `printf "%02d" $im` 

 set files = `echo ${data_dir}/e5.oper.an.sfc/${iyear}${mstr}/e5.oper.an.sfc.*_msl.*.nc`
 set nfile = $#files 

 set ufils  = `echo ${data_dir}/e5.oper.an.sfc/${iyear}${mstr}/e5.oper.an.sfc.*_10u.*.nc`
 set vfils  = `echo ${data_dir}/e5.oper.an.sfc/${iyear}${mstr}/e5.oper.an.sfc.*_10v.*.nc`
 set nuvfil = $#ufils

 set tfils = `echo ${data_dir}/e5.oper.an.pl/${iyear}${mstr}/e5.oper.an.pl.*_t.*.nc`
 set ntfil = $#tfils

 echo $nfile $nuvfil $ntfil

 set i = 1
 while ( $i <= $ntfil ) 
  set tfile  = $tfils[$i]

  set finput = "input_${expnam}_${varnam}_list.txt"
  rm -f $finput; touch $finput

  set tmpstr = "$tfile;$files;$ufils;$vfils;$zsfil"
  echo "$tmpstr" >> ${finput}

  rm out_${expnam}_${varnam}.dat

  $command1 --in_data_list $finput \
            --verbosity 0 \
            --closedcontourcmd "MSL,300.0,4.0,0;_AVG(T(200hPa),T(500hPa)),-0.6,4,0.30" \
            --mergedist 6.0 \
            --searchbymin "MSL" \
            --outputcmd "MSL,min,0;_VECMAG(VAR_10U,VAR_10V),max,2;_DIV(zs,9.81),min,0" \
            --timefilter "6hr" \
            --latname latitude  \
            --lonname longitude \
            --out out_${expnam}_${varnam}.dat

  cat out_${expnam}_${varnam}.dat >> $outfile
  rm out_${expnam}_${varnam}.dat

  @ i++
 end 
 @ im ++ 
end 
 @ iyear++
end 

#set trk_file = ${expnam}_TCS_Track_${syear}-${eyear}.txt
#rm -f $trk_file; touch $trk_file
#
#$command2 --in  $outfile \
#          --out $trk_file \
#          --in_fmt "lon,lat,slp,wind,zs" \
#          --range 8.0 --mintime "54h"\
#          --maxgap "24h" \
#          --threshold "wind,>=,10.0,10;lat,<=,50.0,10;lat,>=,-50.0,10;zs,<=,150.0,10"

exit
