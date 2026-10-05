#!/bin/bash

map="map/map_721x1440_to_180x360_conserve.nc"
outdir="/lcrc/group/e3sm/public_html/diagnostic_output/ac.szhang/v3LE_paper/data/climo"
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
