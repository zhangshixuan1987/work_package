#!/bin/bash -l 
ulimit -S -c 0
#SBATCH -A phy220062      # Allocation name 
#SBATCH --nodes=1         # Total # of nodes (must be 1 for serial job)
#SBATCH --ntasks=128      # Total # of MPI tasks (should be 1 for serial job) 
#SBATCH --time=8:00:00    # Total run time limit (hh:mm:ss)
#SBATCH -J data_process   # Job name
#SBATCH -o myjob.o%j      # Name of stdout output file
#SBATCH -e myjob.e%j      # Name of stderr error file
#SBATCH -p shared        # Queue (partition) name

source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh

#ssh dtn01.nersc.gov
#screen
#bash
#source ~/.bashrc.ext
#conda activate e3sm_unified_latest 
#source /global/common/software/e3sm/anaconda_envs/base/envs/e3sm_unified_latest.csh 
#exit
#screen
data_path="/compyfs/zhan391/initial_test_data"
hpss_path="globus://9cd89cfd-6d04-11e5-ba46-22000b92c6ec/~/EAM_DART/e3sm_experiment/${casename}"
archive_dir="F20TR_ne11_oQU240_DART80_compy F20TR_ne11_oQU240_DART80_DT7200_compy F20TR_ne16_oQU240_DART80_compy F20TR_ne30pg2_EC30to60E2r2_DART80_compy" #"map_file SST_forcing topo"
archive_dir="F20TR_ne30pg2_EC30to60E2r2_DART80_compy F20TR_ne11_oQU240_DART80_compy F20TR_ne16_oQU240_DART80_compy"

for dir in $archive_dir; do 
  cd ${data_path}/${dir}
  if [ ! -d ./zstash ]; then
    mkdir zstash
  fi
  if [ ! -f zstash/index.db ];then 
    zstash create --hpss=none --maxsize 128 . >& zstash_create.log & 
  else
    zstash update --hpss=none >& zstash_create.log & 
  fi 
done 
wait 
