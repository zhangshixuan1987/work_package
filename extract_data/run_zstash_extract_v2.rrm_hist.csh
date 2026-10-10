#!/bin/csh
# Submit this script as : sbatch ./[script-name]
#SBATCH -A m3525
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 12:00:00
#SBATCH -N 2
#SBATCH  --job-name=ncclimo_ctrl
#SBATCH  --output=job%j 

source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.csh
#List of simulations 
#v2.LR.amip_0101_bonus/           v2.NARRM.amip_0201/              v2.NARRM.historical_0151/        
#v2.LR.historical_0101_bonus/     v2.NARRM.amip_0301/              v2.NARRM.historical_0201/        
#v2.NARRM.amip_0101/              v2.NARRM.historical_0101/        v2.NARRM.historical_0251/        
#v2.NARRM.amip_0101_bonus/        v2.NARRM.historical_0101_bonus/  v2.NARRM.piControl/

set workdir = `pwd`
set years = 1979
set yeare = 2014

set hpssdir = "/home/t/tang30/E3SMv2"
set model   = "v2.NARRM"
set type    = "historical"
set rundir  = "archive/atm/hist" 

cd $workdir

foreach ens ( "0101" "0151" "0201" "0251" "0301")
  set case  = "${model}.${type}" 
  set ddir  = "${hpssdir}/${case}_${ens}"
  set iy    = $years
  while ( $iy <= $yeare )
    set yst = `printf "%04d" $iy`
    zstash extract --hpss=$ddir -v "*eam.h2.${yst}-*-*-*.nc" &
    zstash extract --hpss=$ddir -v "*eam.h3.${yst}-*-*-*.nc" & 
    zstash extract --hpss=$ddir -v "*eam.h4.${yst}-*-*-*.nc" & 
    wait
   @ iy++ 
  end 
end 

