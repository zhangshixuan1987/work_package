#!/bin/bash

# Load shared output locations (V3LE_DATA_DIR, ...) from paths.sh next to this
# script. Under sbatch the job runs a spooled copy, so ask Slurm for the original.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  _self="$(scontrol show job "${SLURM_JOB_ID}" | awk -F= '/ Command=/{print $2; exit}')"
else
  _self="${BASH_SOURCE[0]}"
fi
source "$(cd "$(dirname "${_self}")" && pwd)/paths.sh"

data_dir="/lcrc/group/e3sm2/ac.wlin/E3SMv3/AMIP/sstice-ext"
syear=1869
eyear=2022
exp="HadISST"
period="${syear}-${eyear}"
file="${data_dir}/sst_ice_CMIP6_DECK_E3SM_1x1_c20221024.nc"

cd ${data_dir}

outdir="${V3LE_DATA_DIR}/climo"
if [ ! -d ${outdir} ];then
  mkdir -p ${outdir}
fi

for year in `seq $syear $eyear`;do 
  time1=`printf "%04d" $year`"-01-01 00:00:0.0"
  time2=`printf "%04d" $year`"-12-31 23:59:59.0"
  rm -rvf ${outdir}/${exp}_${year}.nc
  ncrcat -d time,"${time1}","${time2}" ${file} ${outdir}/${exp}_${year}.nc 
done 

outfile="${outdir}/${exp}.observation.climo.${period}.nc"
rm -rvf ${outfile}
echo ${outdir}/${exp}_*.nc
if [ ! -f "${outfile}" ];then 
  rm -rvf ${outfile}
  ncea ${outdir}/${exp}_*.nc ${outfile} 
  rm -rvf ${outdir}/${exp}_*.nc
fi 
#ncap2 -s "lon=lon+180.0" ${outfile} ${outfile}.tmp 
#mv ${outfile}.tmp ${outfile}
