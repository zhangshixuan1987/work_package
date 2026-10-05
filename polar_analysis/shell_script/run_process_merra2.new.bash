#!/bin/bash

# Parameters: override from the environment for another machine/experiment,
# e.g.  DATA_DIR=/other/data bash run_process_merra2.new.bash   (or: sbatch --export=ALL,DATA_DIR=... run_process_merra2.new.bash)
DATA_DIR="${DATA_DIR:-/lcrc/group/e3sm/public_html/diagnostic_output/ac.szhang/v3_polar_analysis/data}"

outdir="${DATA_DIR}/climo"
if [ ! -d ${outdir} ];then 
  mkdir -p ${outdir}
fi 

data_dir="/lcrc/group/acme/ac.szhang/acme_scratch/data/merra2/rgd_data"
syear=1985
eyear=2014
syear=2001
eyear=2018
exp="MERRA2"
ens="analysis"
period="${syear}-${eyear}"
files=
for year in `seq ${syear} ${eyear}`;do 
  ystr=`printf "%04d" $year`
  files=( ${files[@]} `echo ${data_dir}/${exp}.1x1.${ystr}.nc` )
done 
echo ${files[@]}
outfile="${outdir}/${exp}.${ens}.climo.${period}.nc"
rm -rvf $outfile
nces ${files[@]} ${outfile} 

#net downward radiative fluxes
ncap2 -s "FLNS=FLNS*-1.0"           ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "FLNSC=FLNSC*-1.0"          ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "FLUSC=FLNSC+FLDSC"         ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}

ncap2 -s "FLNT=FLUT"                 ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "FLNTC=FLUTC"               ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "FLNTOA=FLUT"                 ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "FLNTOAC=FLUTC"               ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "FSNTOA=FSNT"               ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "FSNTOAC=FSNTC"             ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}

ncap2 -s "FSUTOA=SOLIN-FSNTOA"       ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "FSUTOAC=SOLIN-FSNTOAC"     ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}

ncap2 -s "FSUS=FSDS-FSNS"            ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "FSUSC=FSDSC-FSNSC"         ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}

ncap2 -s "RESTOM=FSNT-FLNT"          ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "RESTOMC=FSNTC-FLNTC"       ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "RESTOA=FSNTOA-FLUT"        ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "RESTOAC=FSNTOAC-FLUTC"     ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "SWCF=FSNTOA-FSNTOAC"       ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "LWCF=FLNTOA-FLNTOAC"       ${outfile} ${outfile}.tmp
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
ncap2 -s "LHQFLX=(2.501e6+3.337e5)*QFLX-3.337e5*1.0e3*(PRECT-PRECST)"  ${outfile} ${outfile}.tmp
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

