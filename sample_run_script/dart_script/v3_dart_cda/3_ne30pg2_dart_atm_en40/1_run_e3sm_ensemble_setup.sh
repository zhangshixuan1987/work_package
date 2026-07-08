#!/bin/bash -l
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --job-name=e3sm_dart_ensda_setup 
#SBATCH  --nodes=1
#SBATCH  --output=e3sm_dart_ensda_setup.%j 
#SBATCH  --exclusive 
#SBATCH  --account=esmd
#SBATCH  --time=02:00:00
#SBATCH  --qos=slurm

#source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh
#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.sh

source create_and_setup_case.sh
export do_fetch_code=false
export do_create_newcase=true
export do_case_setup=true
export do_case_build=true
export do_case_submit=false

readonly my_wkdir=${PWD}
readonly CASE_NAME=${my_casename}
readonly CODE_ROOT="${my_e3sm_code}"
readonly CASE_ROOT="${my_runpath}/${my_casename}"
readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/case_scripts
readonly CASE_BUILD_DIR=${CASE_ROOT}/build
readonly CASE_ARCHIVE_DIR=${CASE_ROOT}/archive
readonly CASE_RUN_DIR=${CASE_ROOT}/run
readonly RUN_REFDIR=${CASE_ARCHIVE_DIR}/rest/${my_refdate}-${my_reftod}

if [ ! -d ${my_wkdir}/scripts ] ; then
  mkdir -p ${my_wkdir}/scripts
fi
cd ${my_wkdir}/scripts

run_script="compile_and_setup_e3sm.sh"
cp -rp ${my_wkdir}/compile_and_setup_e3sm_${my_resolution}.${my_compset}.sh  ${run_script}

sed -i "s#MACHINE=.*#MACHINE=\"${my_machine}\"#"                         ${run_script}
sed -i "s#PROJECT=.*#PROJECT=\"${my_project}\"#"                         ${run_script}
sed -i "s#WALLTIME=.*#WALLTIME=\"${my_walltime}\"#"                      ${run_script}
sed -i "s#JOB_SLURM=.*#JOB_SLURM=\"${my_jobqueue}\"#"                    ${run_script}
sed -i "s#JOB_NTASKS=.*#JOB_NTASKS=\"${my_job_ntasks}\"#"                ${run_script}
sed -i "s#TASK_PER_NODE=.*#TASK_PER_NODE=\"${my_task_per_node}\"#"       ${run_script}
sed -i "s#COMPSET=.*#COMPSET=\"${my_compset}\"#"                         ${run_script}
sed -i "s#RESOLUTION=.*#RESOLUTION=\"${my_resolution}\"#"                ${run_script}
sed -i "s#run=.*#run=\"${my_layout}\"#"                                  ${run_script}

sed -i "s#CASE_ARCHIVE_DIR=.*#CASE_ARCHIVE_DIR=\"${CASE_ARCHIVE_DIR}\"#" ${run_script}
sed -i "s#CODE_ROOT=.*#CODE_ROOT=\"${CODE_ROOT}\"#"                      ${run_script}
sed -i "s#CASE_ROOT=.*#CASE_ROOT=\"${CASE_ROOT}\"#"                      ${run_script}
sed -i "s#CASE_NAME=.*#CASE_NAME=\"${CASE_NAME}\"#"                      ${run_script}
sed -i "s#CASE_BUILD_DIR=.*#CASE_BUILD_DIR=\"${CASE_BUILD_DIR}\"#"       ${run_script}
sed -i "s#CASE_RUN_DIR=.*#CASE_RUN_DIR=\"${CASE_RUN_DIR}\"#"             ${run_script}
sed -i "s#CASE_SCRIPTS_DIR=.*#CASE_SCRIPTS_DIR=\"${CASE_SCRIPTS_DIR}\"#" ${run_script}

sed -i "s#START_DATE=.*#START_DATE=\"${my_casedate}\"#"                  ${run_script}
sed -i "s#START_TOD=.*#START_TOD=\"${my_casetod}\"#"                     ${run_script}
sed -i "s#GET_REFCASE=.*#GET_REFCASE=TRUE#"                              ${run_script}
sed -i "s#RUN_REFDATE=.*#RUN_REFDATE=\"${my_casedate}\"#"                ${run_script}
sed -i "s#RUN_REFTOD=.*#RUN_REFTOD=\"${my_casetod}\"#"                   ${run_script}
sed -i "s#RUN_REFCASE=.*#RUN_REFCASE=\"${CASE_NAME}\"#"                  ${run_script}
sed -i "s#RUN_REFDIR=.*#RUN_REFDIR=\"${RUN_REFDIR}\"#"                   ${run_script}

if [ "${do_create_newcase,,}" != "false"  ]; then
  if [ -d ${CASE_SCRIPTS_DIR} ]; then
    rm -rvf ${CASE_SCRIPTS_DIR}
  fi
fi

if [ "${do_case_build,,}" != "false"  ]; then
  if [ -d ${CASE_BUILD_DIR} ]; then
    rm -rvf ${CASE_BUILD_DIR}
  fi
  ./${run_script}
else
  cd ${CASE_SCRIPTS_DIR} 
  my_modelexe="${CASE_BUILD_DIR}/e3sm.exe"
  if [ ! -f ${my_modelexe} ];then
    echo $'\n----- e3sm.exe does not exit, please compile model first-----\n'
    return
  fi
  echo 'WARNING: Setting BUILD_COMPLETE = TRUE.  This is a little risky, but trusting the user.'
  ./xmlchange BUILD_COMPLETE=TRUE
fi

# Loop over members and run ensembles
export do_case_build=false
for i in `seq 1 $my_ensnum`;do 
  my_enscase=EN`printf "%02d" ${i}`
  echo "ens: $i  case: $my_enscase"
  SUB_CASE_NAME=${my_casename}.${my_enscase}
  SUB_CASE_DIR=${CASE_ROOT}/${my_enscase}/case_scripts
  SUB_BUILD_DIR=${CASE_ROOT}/${my_enscase}/build
  SUB_RUN_DIR=${CASE_ROOT}/${my_enscase}/run
  SUB_REFDIR=${RUN_REFDIR}
  if [ -d ${SUB_CASE_DIR} ]; then
     rm -rvf  ${SUB_CASE_DIR}
  fi
  sed -i "s#RUN_REFCASE=.*#RUN_REFCASE=\"${SUB_CASE_NAME}\"#"            ${run_script}
  sed -i "s#RUN_REFDIR=.*#RUN_REFDIR=\"${SUB_REFDIR}\"#"                 ${run_script}
  sed -i "s#CASE_NAME=.*#CASE_NAME=\"${SUB_CASE_NAME}\"#"                ${run_script}
  sed -i "s#CASE_SCRIPTS_DIR=.*#CASE_SCRIPTS_DIR=\"${SUB_CASE_DIR}\"#"   ${run_script}
  sed -i "s#CASE_BUILD_DIR=.*#CASE_BUILD_DIR=\"${SUB_BUILD_DIR}\"#"      ${run_script}
  sed -i "s#CASE_RUN_DIR=.*#CASE_RUN_DIR=\"${SUB_RUN_DIR}\"#"            ${run_script}
  sed -i "s#old_modelexe#\"${my_modelexe}\"#"                            ${run_script}
  echo $run_script 
  ./${run_script}
done

wait 

exit
