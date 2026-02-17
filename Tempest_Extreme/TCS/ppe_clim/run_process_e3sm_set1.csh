#!/bin/csh

#conda activate e3sm_analysis 
#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.csh

set command0 = "DetectNodes"
set command1 = "StitchNodes"
set command2 = "HistogramNodes"

set runnam   = "e3sm_v2_CLIM_F20TR_ne30pg2_EC30to60E2r2"
set expnam   = "CLIM"

set data_dir = "/global/cscratch1/sd/zhan391/DARPA_project/nudged_e3smv2_simulation"
set zsfil    = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/data/e3sm_phis.nc"
set grdfil   = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/grid_info/outCS_ne30pg2_connect.txt"

set varnam   = "TCS"
set syear    = 2007
set eyear    = 2017
set parset   = "set1"

#define the parameter for detection of TempestExtreme 
#reference: Colin M. Zarzycki and Paul A. Ullrich: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016GL071606
#Input variables: 
#VAR_10U, VAR_10U: zonal and meridional wind speeds at 10 m  (reanalysis) or at model bottom level 
#VAR_MSL: mean sea level pressure 
#ZS: surface height 
#Z: geopotential height  
set var_psl  = "PSL"
set var_U10  = "UBOT"
set var_V10  = "VBOT"
set var_ZS   = "PHIS"       
set var_lat  = "lat"
set var_lon  = "lon"

set zs_factor     = 9.81   # unit: m/s2, convert the surface geopotential to m, set to 1.0 if no conversion is needed 
set pslFOmag      = 200.0  # unit: Pa, Strength of local sea level pressure (PSL) minimum [0.5,6]hPa
set pslFOdist     = 5.5    # unit: degree, Allowable distance for PSL closed contour (great circle distance) [1,8]degree
set wcFOmag_T     = -0.6   # unit: K, Strength of warm core anomaly (T) [-2.0,-0.1]K
set wcFOmag_DZ    = -6.0   # unit: m, Strength of warm core anomaly (DZ)[-50,-1]m  
set wcFOdist      = 6.5    # unit: degree, Allowable distance for T/DZ closed contour [1,8]degree
set wcMaxOffset   = 1.0    # unit: degree, Maximum separation between PSL minimum, T/DZ maximum [0,3]degree 
set mergeDist     = 6.0    # unit: degree, Minimum allowable distance between two candidates [0,10]degree
set trajRange     = 8.0    # unit: degree Maximum travel distance for cyclone in 6 h [4,12]degree 
set trajMinLength = "10"   # unit: hour, Minimum cyclone lifetime [12,72]h
set trajMaxGap    = "3"    # unit: hour, Maximum allowable gap in trajectory [0,30]h  
set maxTopo       = 150.0  # unit: m, Maximum topography directly under PSL minimum [10,6000]m
set maxLat        = 50.0   # unit: degree, Maximum latitude of PSL minimum [40,90] degree
set minWind       = 10.0   # unit: m/s , Minimum lowest model level wind speed [10,17.5] m/s
set scidist       = 10     # unit: trajMinLength
set iloncol       = 2      # colomn number for longitude 
set ilatcol       = 3      # column number for latitude 
set tfilter       = "6hr"  # sampling frequency  

if ( ${parset} == "set1" ) then 
  set var_WC1  = "Z200"
  set var_WC2  = "Z500"
  set wcFOmag  = ${wcFOmag_DZ} 
  set varcmd   = "_DIFF"
endif 
 
if ( ${parset} == "set2" ) then
  set var_WC1  = "T200"
  set var_WC2  = "T500"
  set wcFOmag  = ${wcFOmag_T}
  set varcmd   = "_AVG"
endif

if ( ${parset} == "set3" ) then
  set var_WC1  = "Z300"
  set var_WC2  = "Z500"
  set wcFOmag  = ${wcFOmag_DZ}
  set varcmd   = "_DIFF"
endif

if ( ${parset} == "set4" ) then
  set var_WC1  = "T300"
  set var_WC2  = "T500"
  set wcFOmag  = ${wcFOmag_T}
  set varcmd   = "_AVG"
endif

#process the data 
set outfile  = ${expnam}_${parset}_${varnam}_${syear}-${eyear}.txt
rm -f $outfile; touch $outfile

#generate file list 
rm -f file_${parset}_${expnam}_list; touch file_${parset}_${expnam}_list
set iyear = $syear
while ( $iyear <= $eyear )
   echo ${data_dir}/${runnam}/run/*eam.h2*${iyear}*.nc >> file_${parset}_${expnam}_list
 @ iyear++
end

set zfils = `cat file_${parset}_${expnam}_list`
set nzfil = $#zfils
echo $nzfil

set i = 1
while ( $i <= $nzfil ) 
  set zfile  = $zfils[$i]
  set finput = "input_${parset}_${expnam}_${varnam}_list.txt"
  rm -f $finput; touch $finput

  set tmpstr = "$zfile;$zsfil"
  echo "$tmpstr" >> ${finput}

  set tmpout = out_${parset}_${expnam}_${varnam}.dat
  rm -rvf $tmpout 

  $command0 --in_data_list $finput \
            --verbosity 0 \
            --in_connect ${grdfil} \
            --closedcontourcmd "${var_psl},${pslFOmag},${pslFOdist},0;${varcmd}(${var_WC1},${var_WC2}),${wcFOmag},${wcFOdist},${wcMaxOffset}" \
            --mergedist ${mergeDist} \
            --searchbymin "${var_psl}" \
            --outputcmd "${var_psl},min,0;_VECMAG(${var_U10},${var_V10}),max,2;_DIV(${var_ZS},${zs_factor}),min,0" \
            --timefilter ${tfilter}  \
            --latname ${var_lat}  \
            --lonname ${var_lon} \
            --out $tmpout 

  cat $tmpout >> $outfile
  rm -rvf $tmpout

 @ i++
end 

#post-process the track data 

set tdn_file = "${outfile}"
set trk_file = "${expnam}_${parset}_${varnam}_Track_${syear}-${eyear}.txt"
set hst_file = "${expnam}_${parset}_${varnam}_hist_gram_${syear}-${eyear}.nc"
rm -f $trk_file; touch $trk_file

$command1  --in  $tdn_file \
           --out $trk_file \
           --in_connect ${grdfil} \
           --in_fmt "lon,lat,slp,wind,zs" \
           --range ${trajRange} \
           --mintime ${trajMinLength} \
           --maxgap ${trajMaxGap}  \
           --threshold "wind,>=,${minWind},${scidist};lat,<=,${maxLat},${scidist};lat,>=,-${maxLat},${scidist};zs,<=,${maxTopo},${scidist}"

$command2  --in $trk_file \
           --iloncol ${iloncol} \
           --ilatcol ${ilatcol} \
           --out $hst_file

exit

