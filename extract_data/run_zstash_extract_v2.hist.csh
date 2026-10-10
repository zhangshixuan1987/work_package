#!/bin/csh
# Submit this script as : sbatch ./[script-name]
#SBATCH -A e3sm
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH  --job-name=zstash_v21
#SBATCH  --output=job%j 

source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.csh

set workdir = `pwd`
set years = 1850
set yeare = 2014

set hpssdir = "/home/g/golaz/E3SMv2_1"
set model   = "v2_1.LR"
set type    = "historical"
set rundir  = "archive/archive/ice/hist"

cd $workdir

foreach ens ( "0101" "0151" "0201" "0251" "0301")
  set case  = "${model}.${type}" 
  set ddir  = "${hpssdir}/${case}_${ens}"
  set iy    = $years
  while ( $iy <= $yeare )
    set yst = `printf "%04d" $iy`
    zstash extract --hpss=$ddir -v "*mpassi.hist.am.timeSeriesStatsMonthly.${yst}*.nc" &
   @ iy++ 
  end 
  wait
end 

