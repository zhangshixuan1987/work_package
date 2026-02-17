#!/bin/csh

#conda activate e3sm_analysis 
#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.csh

set command = "DetectNodes"
set expnam   = "ERA5"
set data_dir = "/global/cfs/projectdirs/m3522/cmip6/${expnam}"
set zsfil    = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/data/era5_zs.nc"

set varnam   = "ETCS"
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

 set zfils = `echo ${data_dir}/e5.oper.an.pl/${iyear}${mstr}/e5.oper.an.pl.*_z.*.nc`
 set nzfil = $#zfils

 echo $nfile $nuvfil $nzfil

 set i = 1
 while ( $i <= $nzfil )
  set zfile  = $zfils[$i]
  set finput = "input_${expnam}_${varnam}_list.txt"
  rm -f $finput; touch $finput

  set tmpstr = "$zfile;$files;$ufils;$vfils;$zsfil"
  echo "$tmpstr" >> ${finput}

  rm out_${expnam}_${varnam}.dat

  $command --in_data_list $finput \
            --verbosity 0 \
            --closedcontourcmd "MSL,200.0,5.5,0" \
            --noclosedcontourcmd "_AVG(Z(200hPa),Z(500hPa)),-58.8,6.5,1.0" \
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

exit

