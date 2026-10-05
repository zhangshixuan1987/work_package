#!/bin/bash

# Parameters: override from the environment for another machine/experiment,
# e.g.  DATA_DIR=/other/data bash run_process_modis.bash   (or: sbatch --export=ALL,DATA_DIR=... run_process_modis.bash)
DATA_DIR="${DATA_DIR:-/lcrc/group/e3sm/public_html/diagnostic_output/ac.szhang/v3_polar_analysis/data}"

outdir="${DATA_DIR}/climo"
if [ ! -d ${outdir} ];then 
  mkdir -p ${outdir}
fi 
data_dir="/lcrc/group/e3sm/ac.szhang/acme_scratch/data/modis/rgd_data"
#syear=2001
#eyear=2020
syear=2001
eyear=2018
exp="MODIS"
period="${syear}-${eyear}"
file=`echo ${data_dir}/*_1x1_*.nc`
for year in `seq $syear $eyear`;do 
   tmpfil="${outdir}/${exp}_${year}.nc"
   if [ ! -f "${tmpfil}" ];then 
     ncrcat -v "LST_INTERP" -d time,"${year}-01-01 0:00:0.0","${year}-12-31 00:00:0.0" ${file[@]} ${tmpfil}  
     ncrename -v "LST_INTERP,TS" ${tmpfil} ${tmpfil}.tmp 
     mv ${tmpfil}.tmp ${tmpfil}
     ncatted -O -a history,global,o,c,"" ${tmpfil} ${tmpfil}.tmp 
     mv ${tmpfil}.tmp ${tmpfil}
   fi
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
