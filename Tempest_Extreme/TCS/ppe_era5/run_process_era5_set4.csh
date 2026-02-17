#!/bin/csh

#conda activate e3sm_analysis 
#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.csh

set command0 = "DetectNodes"
set command1 = "StitchNodes"
set command2 = "HistogramNodes"
set expnam   = "ERA5"

set data_dir = "/global/cfs/projectdirs/m3522/cmip6/${expnam}"
set zsfil    = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/data/era5_zs.nc"

set varnam   = "TCS"
set syear    = 2007
set eyear    = 2017
set parset   = "set4"

#define the parameter for detection of TempestExtreme 
#reference: Colin M. Zarzycki and Paul A. Ullrich: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016GL071606
#Input variables: 
#VAR_10U, VAR_10U: zonal and meridional wind speeds at 10 m  (reanalysis) or at model bottom level 
#VAR_MSL: mean sea level pressure 
#ZS: surface height 
#Z: geopotential height  
set var_psl  = "MSL"
set var_U10  = "VAR_10U"
set var_V10  = "VAR_10V"
set var_ZS   = "zs"       
set var_WC1  = "Z(200hPa)"
set var_WC2  = "Z(500hPa)"

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
set iloncol       = 3      # colomn number for longitude 
set ilatcol       = 4      # column number for latitude 
set tfilter       = "6hr"  # sampling frequency  

if ( ${parset} == "set1" ) then 
  set var_WC1  = "Z(200hPa)"
  set var_WC2  = "Z(500hPa)"
  set wcFOmag  = ${wcFOmag_DZ} 
  set varcmd   = "_DIFF"
  ##Note: the setup below depends on the data 
  if ($expnam == "ERA5") then 
    set wcFOmag = `echo "${wcFOmag_DZ} * 9.81" | bc -l` 
    echo $wcFOmag
    set fvar0 = "msl"
    set fvar1 = "10u"
    set fvar2 = "10v"
    set fvar3 = "z"
  endif 
endif 
 
if ( ${parset} == "set2" ) then
  set var_WC1  = "T(200hPa)"
  set var_WC2  = "T(500hPa)"
  set wcFOmag  = ${wcFOmag_T}
  set varcmd   = "_AVG"
  ##Note: the setup below depends on the data 
  if ($expnam == "ERA5") then
    set fvar0 = "msl"
    set fvar1 = "10u"
    set fvar2 = "10v"
    set fvar3 = "t"
  endif
endif

if ( ${parset} == "set3" ) then
  set var_WC1  = "Z(300hPa)"
  set var_WC2  = "Z(500hPa)"
  set wcFOmag  = ${wcFOmag_DZ}
  set varcmd   = "_DIFF"
  ##Note: the setup below depends on the data 
  if ($expnam == "ERA5") then
    set wcFOmag = `echo "${wcFOmag_DZ} * 9.81" | bc -l`
    echo $wcFOmag
    set fvar0 = "msl"
    set fvar1 = "10u"
    set fvar2 = "10v"
    set fvar3 = "z"
  endif
endif

if ( ${parset} == "set4" ) then
  set var_WC1  = "T(300hPa)"
  set var_WC2  = "T(500hPa)"
  set wcFOmag  = ${wcFOmag_T}
  set varcmd   = "_AVG"
  ##Note: the setup below depends on the data 
  if ($expnam == "ERA5") then
    set fvar0 = "msl"
    set fvar1 = "10u"
    set fvar2 = "10v"
    set fvar3 = "t"
  endif
endif

#process the data 
set outfile  = ${expnam}_${parset}_${varnam}_${syear}-${eyear}.txt
rm -f $outfile; touch $outfile

set iyear = $syear
while ( $iyear <= $eyear )

set im = 1
while ( $im <= 12 ) 
 
 set mstr  = `printf "%02d" $im` 

 set files  = `echo ${data_dir}/e5.oper.an.sfc/${iyear}${mstr}/e5.oper.an.sfc.*_${fvar0}.*.nc`
 set nfile  = $#files 

 set ufils  = `echo ${data_dir}/e5.oper.an.sfc/${iyear}${mstr}/e5.oper.an.sfc.*_${fvar1}.*.nc`
 set vfils  = `echo ${data_dir}/e5.oper.an.sfc/${iyear}${mstr}/e5.oper.an.sfc.*_${fvar2}.*.nc`
 set nuvfil = $#ufils

 set zfils = `echo ${data_dir}/e5.oper.an.pl/${iyear}${mstr}/e5.oper.an.pl.*_${fvar3}.*.nc`
 set nzfil = $#zfils

 echo $nfile $nuvfil $nzfil

 set i = 1
 while ( $i <= $nzfil ) 
  set zfile  = $zfils[$i]
  set finput = "input_${parset}_${expnam}_${varnam}_list.txt"
  rm -f $finput; touch $finput

  set tmpstr = "$zfile;$files;$ufils;$vfils;$zsfil"
  echo "$tmpstr" >> ${finput}

  set tmpout = out_${parset}_${expnam}_${varnam}.dat
  rm -rvf $tmpout 

  $command0 --in_data_list $finput \
            --verbosity 0 \
            --closedcontourcmd "${var_psl},${pslFOmag},${pslFOdist},0;${varcmd}(${var_WC1},${var_WC2}),${wcFOmag},${wcFOdist},${wcMaxOffset}" \
            --mergedist ${mergeDist} \
            --searchbymin "${var_psl}" \
            --outputcmd "${var_psl},min,0;_VECMAG(${var_U10},${var_V10}),max,2;_DIV(${var_ZS},${zs_factor}),min,0" \
            --timefilter ${tfilter}  \
            --latname latitude  \
            --lonname longitude \
            --out $tmpout 

  cat $tmpout >> $outfile
  rm -rvf $tmpout

  @ i++
 end 
 @ im ++ 
end 
 @ iyear++
end 


#post-process the track data 

set tdn_file = "${outfile}"
set trk_file = "${expnam}_${parset}_${varnam}_Track_${syear}-${eyear}.txt"
set hst_file = "${expnam}_${parset}_${varnam}_hist_gram_${syear}-${eyear}.nc"
rm -f $trk_file; touch $trk_file

$command1  --in  $tdn_file \
           --out $trk_file \
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

