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

set workdir = `pwd`
set syear = 1850
set eyear = 2014

set hpssdir = "/home/projects/e3sm/www/WaterCycle/E3SMv2/LR"
set case    = "v2.LR.piControl"
#set rundir  = "archive/archive/ice/hist"

cd $workdir
set ddir  = "${hpssdir}/${case}"
echo $ddir
zstash ls --hpss=$ddir

