#!/bin/csh

conda activate e3sm_analysis 

set runnam   = "e3sm_v2_NDGUVTQ_F20TR_ne30pg2_EC30to60E2r2"
set expnam   = "NDGUVTQ"
set trknam   = "ppe_ndguvtq"
set feature  = "TCS"
set parset   = "set1"

set freq     = "6hourly"
set syear    = 2007
set eyear    = 2017

set top_dir  = "/global/cscratch1/sd/zhan391/DARPA_project"
set data_dir = "${top_dir}/post_process/data/model_output/${expnam}/${freq}"
set zsfil    = "${top_dir}/post_process/data/model_output/landmsk_coord/landmsk_coord.nc"
set out_dir  = "./filter_data"
set trk_dir  = "${top_dir}/TempestExtremes/${feature}/script_FV/${trknam}"

set fvList   = ("PSL" "UBOT" "VBOT" "Z200" "Z500" "T200" "T500")
set nfvs     = $#fvList

set basin    = "AL"
set lats     = 24
set latn     = 52
set lonw     = 234
set lone     = 294

if ( ! -d $out_dir ) then 
  mkdir -p $out_dir
endif 

# compute the outer radius of each tracked TC
# Following Schenkel et al. (2017) and Stansfield et al. (2020),
# the largest radius outside of the eyewall where the azimuthally averaged 
# wind speed exceeds 8 m s^-1 (r8) is used to measure the size of a TC. 
#generate TC file list 
rm -f tcfile_${parset}_${expnam}_list; touch tcfile_${parset}_${expnam}_list
set iyear = $syear
while ( $iyear <= $eyear )
   set ifv  = 1
   set fstr =
   while ( $ifv <= $nfvs )
     set tmpfil = `ls ${data_dir}/*$fvList[$ifv]_${iyear}*.nc`
     if( $ifv == 1 ) then
       set fstr = $tmpfil
     else
       set fstr = "$fstr;$tmpfil"
     endif
    @ ifv++
   end
   #set fstr = "$fstr;$zsfil"
   echo $fstr >> tcfile_${parset}_${expnam}_list
 @ iyear++
end

set command    = "NodeFileEditor"
set trk_file   = ${trk_dir}/${expnam}_${parset}_${feature}_Track_${syear}-${eyear}.txt
set out_file   = "${out_dir}/${expnam}_${feature}_radprofs_${syear}-${eyear}.txt"
rm -f $out_file; touch $out_file 
set var_psl  = "PSL"
set var_U10  = "UBOT"
set var_V10  = "VBOT"
set var_ZS   = "PHIS"
set var_lat  = "lat"
set var_lon  = "lon"
$command  --in_nodefile  $trk_file \
          --in_data_list tcfile_${parset}_${expnam}_list \
          --in_fmt  "lon,lat,slp,wind,zs" \
          --out_nodefile $out_file \
          --out_fmt "lon,lat,rsize,rprof" \
          --timefilter "6hr" \
          --calculate "rprof=radial_wind_profile($var_U10,$var_V10,39,1.0);rsize=lastwhere(rprof,>,8)"
exit
