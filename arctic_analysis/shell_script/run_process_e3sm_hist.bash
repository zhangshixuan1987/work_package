#!/bin/bash

map="map/map_721x1440_to_180x360_conserve.nc"
outdir="/lcrc/group/e3sm/public_html/diagnostic_output/ac.szhang/v3LE_paper/data/climo"
if [ ! -d ${outdir} ];then 
  mkdir -p ${outdir}
fi 

syear=1985   #2001
eyear=2014   #2018
inty="30yr"


exps=( "v3.LR.historical_0051" "v3.LR.historical_0101" 
       "v3.LR.historical_0101_bcdt15m" "v3.LR.historical_0151" 
       "v3.LR.historical_0201" "v3.LR.historical_0251" )

data_dir="/lcrc/group/e3sm/ac.wlin/E3SMv3"

sub_dir="post/atm/180x360_aave/clim"

for igp in `seq 1 2`;do
  for exp in ${exps[@]};do
    indir="${data_dir}/${exp}/${sub_dir}/${inty}"
    echo $indir
    exit

    period="${syear}-${eyear}"
    files=
    unset files
    for imon in `seq 1 12`;do
      month=`printf "%02d" $imon`
      files=(${files[@]} `echo ${indir}/*_${month}_*${syear}*_*${eyear}*.nc`)
    done
    echo ${files[@]}
    outfile=${outdir}/${exp}.climo.${period}.nc
    rm -rvf ${outfile}
    ncrcat -d time,0, ${files[@]} ${outfile}

    ncap2 -s "FLDSC=FLDS+FLNS-FLNSC"     ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
    ncap2 -s "FLUS=FLDS+FLNS"            ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
    ncap2 -s "FLUSC=FLDSC+FLNSC"         ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
    ncap2 -s "FSUS=FSDS-FSNS"            ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
    ncap2 -s "FSUSC=FSDSC-FSNSC"         ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
    ncap2 -s "FSUTOA=SOLIN-FSNTOA"       ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
    ncap2 -s "FSUTOAC=SOLIN-FSNTOAC"     ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
    ncap2 -s "RESTOM=FSNT-FLNT"          ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
    ncap2 -s "RESTOMC=FSNTC-FLNTC"       ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
    ncap2 -s "RESTOA=FSNTOA-FLUT"        ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
    ncap2 -s "RESTOAC=FSNTOAC-FLUTC"     ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
    ncap2 -s "NETCF=SWCF+LWCF"           ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
    ncap2 -s "SWCF_SRF=FSNS-FSNSC"       ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
    ncap2 -s "LWCF_SRF=-(FLNS-FLNSC)"    ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
    ncap2 -s "NETWCF_SRF=SWCF_SRF+LWCF_SRF"   ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
    ncap2 -s "NETFLX4_SRF=FSNS-FLNS-LHFLX-SHFLX"              ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
    ncap2 -s "NETFLX6_SRF=FSDS-FSUS+(FLDS-FLUS)-LHFLX-SHFLX"  ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
    ncap2 -s "LHQFLX=(2.501e6 + 3.337e5)*QFLX-3.337e5*1.0e3*(PRECC+PRECL-PRECSC-PRECSL)"  ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
    ncap2 -s "NETFLX8_SRF=FSNS-FLNS-LHQFLX-SHFLX"  ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
    ncap2 -s "ALBEDO_SRF=FSUS/FSDS"        ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
    ncap2 -s "ALBEDOC_SRF=FSUSC/FSDSC"     ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
    ncap2 -s "ALBEDO=FSUTOA/SOLIN"         ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile} 
    ncap2 -s "ALBEDOC=FSUTOAC/SOLIN"       ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
  done
done  
