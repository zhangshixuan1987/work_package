#!/bin/bash

# Parameters: override from the environment for another machine/experiment,
# e.g.  DATA_DIR=/other/data bash run_process_HadSST.bash   (or: sbatch --export=ALL,DATA_DIR=... run_process_HadSST.bash)
DATA_DIR="${DATA_DIR:-/lcrc/group/e3sm/public_html/diagnostic_output/ac.szhang/v3_polar_analysis/data}"

data_dir="/lcrc/group/e3sm/ac.szhang/acme_scratch/data/HadISST"
#syear=2001
#eyear=2020
syear=2001
eyear=2018
exp="HadISST"
period="${syear}-${eyear}"
file="${data_dir}/rgd_data/ts_186901_202212.nc"

cd ${data_dir}

outdir="${DATA_DIR}/climo"
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
