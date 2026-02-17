#!/bin/csh

module load cdo

set expnam   = "ERA5"
set data_dir = "/global/cfs/projectdirs/m3522/cmip6/${expnam}"
set out_dir  = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/data/${expnam}"
set syear    = 2007
set eyear    = 2017

set var_List  = ( "U10" "V10" "PSL" "PHIS" "Z" "U" "V" "T" "Q" )
set var_ERA5  = ( "10u" "10v" "msl"    "z" "z" "u" "v" "t" "q" )
set var_key   = ( "an.sfc" "an.sfc" "an.sfc" "invariant" "an.pl" "an.pl" "an.pl" "an.pl" )
set nvars     = $#var_List

if ( ! -d  $out_dir )then 
  mkdir -p $out_dir
endif 

set iyear = $syear
while ($iyear <= $eyear) 

  set iv = 1
  while ( $iv <= $nvars ) 

   set vout = $var_List[$iv]
   set vint = $var_ERA5[$iv]
   set vkey = $var_key[$iv]

   echo $vout $vint $vkey

   rm -rvf $out_dir/${vout}_${iyear}.nc
   if ($vout == "PHIS") then
     set file_list = `echo ${data_dir}/e5.oper.${vkey}/197901/e5.oper.${vkey}.*_z.*.nc` 
     cdo select,name=Z    $file_list $out_dir/${vout}_${iyear}.nc 
     ncrename -v Z,PHIS $out_dir/${vout}_${iyear}.nc
   else
     set file_list = `echo ${data_dir}/e5.oper.${vkey}/${iyear}*/e5.oper.${vkey}.*_${vint}.*.nc`   
     echo $file_list
     if($vout == "PSL") then 
       ncrcat -d time,,,6  -v MSL $file_list $out_dir/${vout}_${iyear}.nc
     else if($vout == "U10") then
       ncrcat -d time,,,6 -v VAR_10U $file_list $out_dir/${vout}_${iyear}.nc
       ncrename -v VAR_10U,U10 $out_dir/${vout}_${iyear}.nc
     else if($vout == "V10") then
       ncrcat -d time,,,6  -v VAR_10V $file_list $out_dir/${vout}_${iyear}.nc
       ncrename -v VAR_10V,V10 $out_dir/${vout}_${iyear}.nc
     else
       ncrcat -d time,,,6 -v $vout $file_list $out_dir/${vout}_${iyear}.nc
     endif 
   endif 

   @ iv++  
  end 

  @ iyear++
end 
