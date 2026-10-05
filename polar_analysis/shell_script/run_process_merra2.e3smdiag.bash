#!/bin/bash

# Parameters: override from the environment for another machine/experiment,
# e.g.  DATA_DIR=/other/data bash run_process_merra2.e3smdiag.bash   (or: sbatch --export=ALL,DATA_DIR=... run_process_merra2.e3smdiag.bash)
DATA_DIR="${DATA_DIR:-/lcrc/group/e3sm/public_html/diagnostic_output/ac.szhang/v3_polar_analysis/data}"

outdir="${DATA_DIR}/climo"
if [ ! -d ${outdir} ];then 
  mkdir -p ${outdir}
fi 

data_dir="/lcrc/soft/climate/e3sm_diags_data/obs_for_e3sm_diags/climatology"
map="${DATA_DIR}/map/map_361x576_to_180x360_conserve.nc"

syear=1980
eyear=2016
exp="MERRA2"
period="${syear}-${eyear}"
files=
for mon in `seq 1 12`;do 
  mstr=`printf "%02d" $mon`
  files=( ${files[@]} `echo ${data_dir}/${exp}/*_${mstr}_*.nc` )
done 
echo ${files[@]}

ens="analysis"
outfile="${outdir}/${exp}.${ens}.climo.${period}.nc"
rm -vf ${outfile}.tmp ${outfile}
ncrcat -d time,0, ${files[@]} ${outfile}.tmp
ncremap -m ${map} -i ${outfile}.tmp -o ${outfile}
rm -rvf ${outfile}.tmp

ncrename -v "uas,U10"         ${outfile}
ncrename -v "vas,V10"         ${outfile}
ncrename -v "huss,QREFHT"     ${outfile}
ncrename -v "tas,TREFHT"      ${outfile}
ncrename -v "pr,PRECT"        ${outfile}
ncrename -v "prw,TMQ"         ${outfile}
ncrename -v "psl,PSL"         ${outfile}
ncrename -v "tasmin,TREFMNAV" ${outfile}
ncrename -v "tasmax,TREFMXAV" ${outfile}
ncrename -v "tauu,TAUX"       ${outfile}
ncrename -v "tauv,TAUY"       ${outfile}
ncrename -v "ta,T"            ${outfile}
ncrename -v "ua,U"            ${outfile}
ncrename -v "va,V"            ${outfile}
ncrename -v "hus,Q"           ${outfile}
ncrename -v "wap,OMEGA"       ${outfile}
ncrename -v "zg,Z3"           ${outfile}
ncrename -v "hur,RELHUM"      ${outfile}
ncrename -v "hfss,SHFLX"      ${outfile}
ncrename -v "hfls,LHFLX"      ${outfile}
ncrename -v "rlut,FLUT"       ${outfile}
ncrename -v "rlutcs,FLUTC"    ${outfile}
ncrename -v "rsds,FSDS"       ${outfile}
ncrename -v "rsdscs,FSDSC"    ${outfile}
ncrename -v "rlds,FLDS"       ${outfile}
ncrename -v "rsdt,SOLIN"      ${outfile}
ncrename -v "rsut,FSUTOA"     ${outfile}
ncrename -v "rsutcs,FSUTOAC"  ${outfile}

ncap2 -s "FLNS=rlus-FLDS"        ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile} 
ncap2 -s "FSNS=FSDS-rsus"        ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "FSNTOA=SOLIN-FSUTOA"   ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "FSNTOAC=SOLIN-FSUTOAC" ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "RESTOA=FSNTOA-FLUT"    ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "RESTOM=FSNTOA-FLUT"    ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "FLNT=FSNTOA-RESTOM"    ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}

#ncap2 -s "FSNSC=FSDSC-rsuscs"     ${outfile} ${outfile}.tmp
#mv ${outfile}.tmp ${outfile}
#ncap2 -s "FLNSC=FLDS+FLNS-rldscs" ${outfile} ${outfile}.tmp
#mv ${outfile}.tmp ${outfile}
#ncap2 -s "RESTOM=SOLIN-FSUT-FLUT"  ${outfile} ${outfile}.tmp
#mv ${outfile}.tmp ${outfile}

#FLUS = FLNS + FLDS
#FLUSC = FLNSC + FLDSC
#rtmt = FSNT - FLNT
#rlut = FSNTOA - (FSNT - FLNT)
