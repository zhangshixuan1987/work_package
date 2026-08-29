#!/bin/csh
# Submit this script as : sbatch ./[script-name]
#SBATCH -A m3525
#SBATCH -q regular
#SBATCH -t 48:00:00
#SBATCH -N 2
#SBATCH  --job-name=ncclimo_ctrl
#SBATCH  --output=job%j 
#SBATCH  --exclusive 
#SBATCH  --constraint=knl,quad,cache

#ssh dtn01.nersc.gov
#screen
#bash
#source ~/.bashrc.ext
#conda activate e3sm_unified_latest 
#source /global/common/software/e3sm/anaconda_envs/base/envs/e3sm_unified_latest.csh 
#exit

set workdir = "./zstash"
set hpssdir = "/home/z/zhan391/DARPA_DIAG_Data/NDGUVTQ_SRF1"
cd $workdir
zstash extract --hpss=$hpssdir 
wait

