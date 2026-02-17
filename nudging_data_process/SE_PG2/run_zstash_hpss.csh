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
#source ~/.bashrc.ext
#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.sh

set workdir = `pwd`
set hpssdir = "/home/z/zhan391/DARPA_ML_Train_Data/SE_PG2"
cd $workdir

foreach dir ( visualization fix_time clim_sfc clim \
              before_nudging after_nudging nudging_tendency \
              clim_30yr_phase2 SSP245_30yr_phase2 SSP585_30yr_phase2 \
              reference_data reference_data_nudge ) # reference_data_tempest visualization )

cd $workdir/$dir 

#zstash create --hpss=$hpssdir/$dir --maxsize 128 . >& zstash_create.log 
zstash update --hpss=$hpssdir/$dir  >& zstash_create.log

wait

end

