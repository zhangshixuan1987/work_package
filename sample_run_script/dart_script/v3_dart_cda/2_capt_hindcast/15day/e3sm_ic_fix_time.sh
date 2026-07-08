#!/bin/bash -l
#################################################################################
# This script attempts to prepare initial condition files for startup of DART DA 
#################################################################################
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --job-name=e3sm_dart_ensda_init 
#SBATCH  --nodes=1
#SBATCH  --output=e3sm_dart_ensda_init.%j 
#SBATCH  --exclusive 
#SBATCH  --account=esmd
#SBATCH  --time=02:00:00
#SBATCH  --qos=short

environment_ext=/share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh
source ${environment_ext}

source create_and_setup_case.sh

my_wkdir=${PWD}

# Loop over members
if [ -f "${my_refeam_ic}" ];then
  echo "reference file: ${my_refeam_ic}"
  ncatted -O -a fv_nphys,global,c,l,2 ${my_refeam_ic} ${my_refeam_ic}
  ncatted -O -a np,global,c,l,4       ${my_refeam_ic} ${my_refeam_ic}
  ncatted -O -a ne,global,c,l,30      ${my_refeam_ic} ${my_refeam_ic}
fi

wait
 
exit
