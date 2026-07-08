#!/bin/csh
# Submit this script as : sbatch ./[script-name]

#SBATCH -A m3089
#SBATCH -N 1
##SBATCH -q regular
#SBATCH -q scavenger
#SBATCH -C knl,quad,cache
#SBATCH -o DT1800.cam.h1
#SBATCH -t 04:00:00

#Set the minimal environment
#module purge
#unsetenv LD_LIBRARY_PATH
#module load intel/16.0.2
#module load mvapich2/2.2b

setenv OMP_NUM_THREADS 1
module load anaconda3
source /share/apps/anaconda3/2019.03/etc/profile.d/conda.csh
conda activate eralib

set dtimeList = (1800)
set irS = 1
set irE = 1
set casestr = "newLscale.alpha5_59.F2010.ne30pg2_r05_oECv3.compy"
set rgdstr  = ${casestr} #"v2_Baseline"

#Script name and path
set script_name1 = ncremap
set script_name2 = ncrcat

set CLIMO_SCRIPT1 = $script_name1
set CLIMO_SCRIPT2 = $script_name2

#Mapping fil
set MAP_FILE      = /qfs/people/zender/data/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc
#location of model history file
set RUN_FILE_DIR  =  /compyfs/zhan391/clubb_test_run_nudg/${casestr}/run/

#Output climo files in FV grid
set REGRIDDED_DIR = /compyfs/zhan391/clubb_test_run_nudg/post/${rgdstr}

if ( ! -d $REGRIDDED_DIR ) then
  mkdir -p $REGRIDDED_DIR
endif

cd ${RUN_FILE_DIR}/

foreach file (*eam.h0* *eam.h1*)
  set rgdfil = `basename $file`
  rm -rvf ${REGRIDDED_DIR}/${rgdfil}
  $CLIMO_SCRIPT1 -i ${file} -m ${MAP_FILE} -o ${REGRIDDED_DIR}/${rgdfil}
end  #----------------------------------------------

