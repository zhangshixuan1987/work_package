#!/bin/bash

# Load shared output locations (V3LE_DATA_DIR, ...) from paths.sh next to this
# script. Under sbatch the job runs a spooled copy, so ask Slurm for the original.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  _self="$(scontrol show job "${SLURM_JOB_ID}" | awk -F= '/ Command=/{print $2; exit}')"
else
  _self="${BASH_SOURCE[0]}"
fi
source "$(cd "$(dirname "${_self}")" && pwd)/paths.sh"

map="map/map_721x1440_to_180x360_conserve.nc"
outdir="${V3LE_DATA_DIR}/climo"
if [ ! -d ${outdir} ];then 
  mkdir -p ${outdir}
fi 

#for exp in "ERA5","v3.LR.amip_0101","v3.LR.historical_0101";do 
exp="ERA5"
data_dir="/lcrc/soft/climate/e3sm_diags_data/obs_for_e3sm_diags/climatology/ERA5"
syear=1979
eyear=2019
period="${syear}-${eyear}"
for imon in `seq 1 12`;do 
  month=`printf "%02d" $imon`
  file=`echo ${data_dir}/*_${month}_*${syear}*_*${eyear}*.nc`
  ncremap -m ${map}  -i ${file} -o ${outdir}/${exp}_${month}.nc 
done 
echo ${outdir}/${exp}_*.nc
rm -rvf ${outdir}/${exp}.climo.${period}.nc
ncrcat -d time,0, ${outdir}/${exp}_*.nc ${outdir}/${exp}.e3sm.diag.${period}.nc
rm -rvf ${outdir}/${exp}_*.nc
