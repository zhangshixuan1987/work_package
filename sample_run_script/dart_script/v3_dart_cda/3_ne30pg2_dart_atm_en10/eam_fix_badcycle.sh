#!/bin/bash -el 
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --account=esmd
#SBATCH  --time=2:00:00
#SBATCH  --partition=short
#SBATCH  --job-name=e3sm_dart_ensda_cyc 
#SBATCH  --nodes=20
#SBATCH  --output=e3sm_dart_ensda_cyc.%j 
#SBATCH  --exclusive
#SBATCH  --no-kill
#SBATCH  --requeue

#For cshell:
#limit stacksize unlimited
#limit datasize unlimited

#For bash 
#ulimit -s unlimited
#ulimit -d unlimited

#export SLURM_NNODES=20
#export SLURM_NTASKS=800

echo == Start of eam_fix_badcycle.sh ==
date
echo ============================================

#source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh
#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.sh
source /qfs/people/zhan391/e3sm_dart_work/code/DART/models/eam-se/work/env_mach_specific.sh

my_wkdir=${PWD}

cd ${my_wkdir}
source ./create_and_setup_case.sh

export ENVIRONMENT_RUNNING_E3SM_UNIFIED_USE_ANOTHER_TERMINAL='e3sm'

DATA_ASSIMILATION_ATM=TRUE
DATA_ASSIMILATION_CYCLES=${my_dart_cycle}
DATA_ASSIMILATION_WINDOW=${my_dart_window}

CASE_ROOT=${my_modeldir}/EN01/case_scripts
RUN_ROOT=${my_modeldir}/EN01/run
ARCHIVE_DIR="${my_modeldir}/archive"

#my_fixdate="2011-12-27"
#my_fixtod="43200"
my_fixdate="2011-12-28"
my_fixtod="43200"
if [ -d "${CASE_ROOT}" ]; then
  cd ${CASE_ROOT}
  CUR_YMD=${my_fixdate}
  CUR_TOD=${my_fixtod}
else
  echo "ERROR: Case directory does not exist: ${my_modelcase}"
  exit
fi

#determine the time for previous DA cycle 
CUR_DATE=( `echo ${CUR_YMD}-${CUR_TOD} | sed -e "s#-# #g"` )
CUR_YEAR=`echo "${CUR_DATE[0]}" | bc`
CUR_MONTH=`echo "${CUR_DATE[1]}" | bc`
CUR_DAY=`echo "${CUR_DATE[2]}" | bc`
CUR_HOUR=`echo "${CUR_DATE[3]}" / 3600 | bc`
CUR_SECONDS=`echo "${CUR_DATE[3]}" | bc`
echo "valid time for eam forecast cycle is $CUR_YEAR $CUR_MONTH $CUR_DAY $CUR_SECONDS (seconds)"

# First Step: Loop over members and run 6-hour ensembel forecast
for i in `seq 1 ${my_ensnum}`;do
  echo === Starting member ${i} ===
  ENSTR=EN`printf "%02d" ${i}`
  CASE_NAME=${my_casename}.${ENSTR}
  CASE_DIR=`echo ${CASE_ROOT} | sed "s/EN01/${ENSTR}/g"`
  RUN_DIR=`echo ${RUN_ROOT} | sed "s/EN01/${ENSTR}/g"`
  REF_DIR="${ARCHIVE_DIR}/rest/${CUR_YMD}-${CUR_TOD}"
  cd ${CASE_DIR}
  ./xmlchange run_exe="--kill-on-bad-exit=1 --job-name=${CASE_NAME} \${EXEROOT}/e3sm.exe "
  ./xmlchange RUN_TYPE="hybrid"
  ./xmlchange CONTINUE_RUN=FALSE
  ./xmlchange BFBFLAG=FALSE
  ./xmlchange RUN_STARTDATE="${CUR_YMD}"
  ./xmlchange START_TOD="${CUR_TOD}"
  ./xmlchange REST_OPTION="nhours"
  ./xmlchange REST_N="${DATA_ASSIMILATION_WINDOW}"
  ./xmlchange STOP_OPTION="nhours"
  ./xmlchange STOP_N="${DATA_ASSIMILATION_WINDOW}"
  ./xmlchange GET_REFCASE=FALSE
  ./xmlchange RUN_REFCASE="${CASE_NAME}"
  ./xmlchange RUN_REFDATE="${CUR_YMD}"
  ./xmlchange RUN_REFTOD="${CUR_TOD}"
  ./xmlchange RUN_REFDIR="${REF_DIR}"
  ./xmlchange DOUT_S=True  #FALSE
  ./xmlchange DOUT_S_ROOT="${ARCHIVE_DIR}"

  echo ============================

done

# That's all folks!
sleep 10

echo ===== End of eam_fix_badcycle.sh =====
date
echo =====================================
