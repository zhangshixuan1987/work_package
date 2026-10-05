#!/bin/bash

outdir="/lcrc/group/e3sm/public_html/diagnostic_output/ac.szhang/v3LE_paper/data/climo"
if [ ! -d ${outdir} ];then 
  mkdir -p ${outdir}
fi 

data_dir="/lcrc/group/e3sm/ac.szhang/acme_scratch/data/era5"
exp="ERA5"

for igp in `seq 2 2`;do 
  if [[ $igp == 1 ]];then 
    syear=1985
    eyear=2014
  else
    syear=2001
    eyear=2018
  fi 
  period="${syear}-${eyear}"
  for ens in analysis ens01 ens02 ens03 ens04 ens05 ens06 ens07 ens08 ens09 ens10;do 
    file=`echo ${data_dir}/*_${ens}_*.nc`
    echo $file
    for year in `seq $syear $eyear`;do 
       tmpfil="${outdir}/${exp}_${year}.nc"
       if [ ! -f "${tmpfil}" ];then 
         ncks -d time,"${year}-01-01 0:00:0.0","${year}-12-31 00:00:0.0" ${file} ${tmpfil}  
       fi   
       #ncdump -t -v time ${tmpfil}
    done

    echo ${outdir}/${exp}_*.nc
    outfile="${outdir}/${exp}.${ens}.climo.${period}.nc"
    rm -rvf ${outfile}
    ncea ${outdir}/${exp}_*.nc ${outfile} 
    rm -rvf ${outdir}/${exp}_*.nc
    ncap2 -s "FLNT=FLNTOA"               ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
    ncap2 -s "FLNTC=FLNTOAC"             ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
    ncap2 -s "FSNT=FSNTOA"               ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
    ncap2 -s "FSNTC=FSNTOAC"             ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
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
    ncap2 -s "LHQFLX=(2.501e6+3.337e5)*MEP-3.337e5*1.0e3*(PRECT-PRECST)"  ${outfile} ${outfile}.tmp
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
