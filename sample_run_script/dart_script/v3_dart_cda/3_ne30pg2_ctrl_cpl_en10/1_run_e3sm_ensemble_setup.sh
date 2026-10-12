#!/bin/bash -el
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --job-name=e3sm_dart_ensda_setup 
#SBATCH  --nodes=1
#SBATCH  --output=e3sm_dart_ensda_setup.%j 
#SBATCH  --exclusive 
#SBATCH  --account=esmd
#SBATCH  --time=02:00:00
#SBATCH  --qos=short

set -Eeuo pipefail
#source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh
#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.sh

fail() {
  echo "ERROR: $*"
  exit 1
}

check_command() {
  local cmd="$1"
  command -v "${cmd}" >/dev/null 2>&1 || fail "required command not found: ${cmd}"
}

validate_positive_int() {
  local name="$1"
  local value="$2"
  if ! [[ "${value}" =~ ^[0-9]+$ ]] || (( value < 1 )); then
    fail "${name} must be a positive integer, got: ${value}"
  fi
}

SCRIPT_NAME=$(basename "$0")
RUN_START_EPOCH=$(date +%s)
on_exit_summary() {
  local exit_code=$?
  local run_end_epoch
  local elapsed
  run_end_epoch=$(date +%s)
  elapsed=$((run_end_epoch - RUN_START_EPOCH))
  if (( exit_code == 0 )); then
    echo "INFO: ${SCRIPT_NAME} completed successfully in ${elapsed}s at $(date '+%Y-%m-%d %H:%M:%S')"
  else
    echo "ERROR: ${SCRIPT_NAME} failed with exit ${exit_code} after ${elapsed}s at $(date '+%Y-%m-%d %H:%M:%S')"
  fi
}
trap on_exit_summary EXIT
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  WORK_DIR=$(readlink -f "${SLURM_SUBMIT_DIR}") || fail "cannot resolve SLURM_SUBMIT_DIR"
else
  SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || fail "cannot resolve script path"
  WORK_DIR=$(dirname "${SCRIPT_PATH}")
fi
cd "${WORK_DIR}"
source "${WORK_DIR}/create_and_setup_case.sh"
mkdir -p "${my_log_dir}" "${my_status_dir}" "${my_lock_dir}" "${my_run_script_dir}"

for cmd in sed cp dirname flock mkdir readlink rm; do
  check_command "${cmd}"
done

validate_positive_int "my_ensnum" "${my_ensnum}"
exec 8>"${my_lock_dir}/step1_setup.lock"
flock -n 8 || fail "another step-1 setup is already running"

mkdir -p "${my_status_dir}"
setup_status="${my_status_dir}/setup_complete"
setup_in_progress="${my_status_dir}/setup_in_progress"
rm -f -- "${setup_status}"
printf 'valid_time=%s-%s\ncase=%s\nensemble_size=%s\nslurm_job_id=%s\nstarted_at=%s\n' "${my_refdate}" "${my_reftod}" "${my_casename}" "${my_ensnum}" "${SLURM_JOB_ID:-none}" "$(date '+%Y-%m-%d %H:%M:%S')" > "${setup_in_progress}"

export do_fetch_code=false
export do_create_newcase=true
export do_case_setup=true
export do_case_build=true
export do_case_submit=false

readonly CASE_NAME=${my_casename}
readonly my_wkdir="${WORK_DIR}"
readonly CODE_ROOT="${my_e3sm_code}"
readonly CASE_ROOT="${my_runpath}/${my_casename}"
readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/case_scripts
readonly CASE_BUILD_DIR=${CASE_ROOT}/build
readonly CASE_ARCHIVE_DIR=${CASE_ROOT}/archive
readonly CASE_RUN_DIR=${CASE_ROOT}/run
readonly RUN_REFDIR=${CASE_ARCHIVE_DIR}/rest/${my_refdate}-${my_reftod}

safe_remove_case_dir() {
  local target="$1"
  case "${target}" in
    "${CASE_ROOT}/case_scripts"|"${CASE_ROOT}/build"|"${CASE_ROOT}"/EN[0-9][0-9]/case_scripts) ;;
    *) fail "refusing to remove unexpected setup path: ${target}" ;;
  esac
  rm -rvf -- "${target}"
}

mkdir -p "${my_run_script_dir}"
cd "${my_run_script_dir}"

run_script="compile_and_setup_e3sm.sh"
template_dir="${my_workflow_lib}/run_template"
template_script="${template_dir}/compile_and_setup_e3sm.${my_compset}.sh"
[[ -f "${template_script}" ]] || fail "missing template script: ${template_script}"
cp -p -- "${template_script}" "${run_script}"
chmod +x "${run_script}"

sed -i "s#MACHINE=.*#MACHINE=\"${my_machine}\"#"                         ${run_script}
sed -i "s#PROJECT=.*#PROJECT=\"${my_project}\"#"                         ${run_script}
sed -i "s#WALLTIME=.*#WALLTIME=\"${my_walltime}\"#"                      ${run_script}
sed -i "s#JOB_SLURM=.*#JOB_SLURM=\"${my_jobqueue}\"#"                    ${run_script}
sed -i "s#JOB_NTASKS=.*#JOB_NTASKS=\"${my_task_per_node}\"#"                ${run_script}
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
    safe_remove_case_dir "${CASE_SCRIPTS_DIR}"
  fi
fi

if [ "${do_case_build,,}" != "false" ]; then
  if [ -d "${CASE_BUILD_DIR}" ]; then
    safe_remove_case_dir "${CASE_BUILD_DIR}"
  fi
  ./${run_script}
else
  cd ${CASE_SCRIPTS_DIR} 
  my_modelexe="${CASE_BUILD_DIR}/e3sm.exe"
  if [ ! -f ${my_modelexe} ];then
    echo $'\n----- e3sm.exe does not exit, please compile model first-----\n'
    exit 1
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
  SUB_ARCHIVE_DIR=${CASE_ROOT}/${my_enscase}/archive
  SUB_REFDIR=${RUN_REFDIR}
  if [ -d "${SUB_CASE_DIR}" ]; then
     safe_remove_case_dir "${SUB_CASE_DIR}"
  fi
  sed -i "s#RUN_REFCASE=.*#RUN_REFCASE=\"${SUB_CASE_NAME}\"#"            ${run_script}
  sed -i "s#RUN_REFDIR=.*#RUN_REFDIR=\"${SUB_REFDIR}\"#"                 ${run_script}
  sed -i "s#CASE_NAME=.*#CASE_NAME=\"${SUB_CASE_NAME}\"#"                ${run_script}
  sed -i "s#CASE_SCRIPTS_DIR=.*#CASE_SCRIPTS_DIR=\"${SUB_CASE_DIR}\"#"   ${run_script}
  sed -i "s#CASE_BUILD_DIR=.*#CASE_BUILD_DIR=\"${SUB_BUILD_DIR}\"#"      ${run_script}
  sed -i "s#CASE_RUN_DIR=.*#CASE_RUN_DIR=\"${SUB_RUN_DIR}\"#"            ${run_script}
  sed -i "s#CASE_ARCHIVE_DIR=.*#CASE_ARCHIVE_DIR=\"${SUB_ARCHIVE_DIR}\"#" ${run_script}
  sed -i "s#old_modelexe#\"${my_modelexe}\"#"                            ${run_script}
  echo $run_script 
  ./${run_script}
done

[[ -x "${my_modelexe}" ]] || fail "model executable missing after setup: ${my_modelexe}"
for i in `seq 1 ${my_ensnum}`; do
  member=$(printf "EN%02d" "${i}")
  member_case="${CASE_ROOT}/${member}/case_scripts"
  [[ -x "${member_case}/xmlquery" && -s "${member_case}/env_run.xml" ]] || fail "incomplete case setup for ${member}: ${member_case}"
done
setup_tmp="${setup_status}.tmp.${SLURM_JOB_ID:-$$}"
printf 'valid_time=%s-%s\ncase=%s\nensemble_size=%s\nslurm_job_id=%s\ncompleted_at=%s\n' "${my_refdate}" "${my_reftod}" "${CASE_NAME}" "${my_ensnum}" "${SLURM_JOB_ID:-none}" "$(date '+%Y-%m-%d %H:%M:%S')" > "${setup_tmp}"
mv -f "${setup_tmp}" "${setup_status}"
rm -f -- "${setup_in_progress}"
echo "Wrote step-1 completion record: ${setup_status}"
wait 

exit
