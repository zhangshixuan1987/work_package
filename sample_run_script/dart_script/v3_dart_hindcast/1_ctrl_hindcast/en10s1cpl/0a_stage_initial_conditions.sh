#!/bin/bash -el
###############################################################################
# Step 0a (hindcasts with my_ic_source=perturbed) -- stage the unperturbed initial conditions of every member.
# Derived from the v3_dart_cda Step 2 script; writes runtmp/status/icbc_complete.<date>.
###############################################################################
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --job-name=e3sm_hindcast_icbc
#SBATCH  --nodes=1
#SBATCH  --output=runtmp/logs/e3sm_hindcast_icbc.%j 
#SBATCH  --exclusive 
#SBATCH  --account=esmd
#SBATCH  --time=02:00:00
#SBATCH  --partition=short

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
  my_wkdir=$(readlink -f "${SLURM_SUBMIT_DIR}") || fail "cannot resolve SLURM_SUBMIT_DIR"
else
  script_path=$(readlink -f "${BASH_SOURCE[0]}") || fail "cannot resolve script path"
  my_wkdir=$(dirname "${script_path}")
fi
cd "${my_wkdir}"
source "${my_wkdir}/create_and_setup_case.sh"

[[ -n "${my_conda_setup_file:-}" ]] || fail "my_conda_setup_file is unset"
[[ -r "${my_conda_setup_file}" ]] || fail "configured Conda setup is not readable: ${my_conda_setup_file}"
[[ -n "${my_analysis_conda_env:-}" ]] || fail "my_analysis_conda_env is unset"
echo "Activating configured analysis environment: ${my_analysis_conda_env}"
source "${my_conda_setup_file}"
conda activate "${my_analysis_conda_env}" || fail "could not activate Conda environment: ${my_analysis_conda_env}"
mkdir -p "${my_log_dir}" "${my_status_dir}" "${my_lock_dir}"

for cmd in awk bc ncks ncap2 ncrename ncdump dirname flock mkdir cp mv readlink rm touch; do
  check_command "${cmd}"
done

validate_positive_int "my_ensnum" "${my_ensnum}"
exec 8>"${my_lock_dir}/step0a_icbc.lock"
flock -n 8 || fail "another Step 0a IC/BC staging is already running"

# Hindcast case ownership (same rule as workflow_lib/common.sh): one experiment
# directory per case directory. Step 0a claims a new case; both steps refuse a
# case that belongs to another experiment directory.
hindcast_owner_file="${my_modeldir}/.hindcast_owner"
hindcast_owner=$(cat "${hindcast_owner_file}" 2>/dev/null || true)
if [[ -n "${hindcast_owner}" && "${hindcast_owner}" != "${my_workflow_root}" ]]; then
  fail "case ${my_modeldir} belongs to experiment directory ${hindcast_owner}"
fi
if [[ -z "${hindcast_owner}" ]]; then
  [[ ! -d "${my_modeldir}" ]] || fail "case ${my_modeldir} exists but was not created by this workflow; move it aside or claim it by writing ${my_workflow_root} to ${hindcast_owner_file}"
  mkdir -p "${my_modeldir}"
  echo "${my_workflow_root}" > "${hindcast_owner_file}"
fi

mkdir -p "${my_status_dir}"
icbc_status="${my_status_dir}/icbc_complete.${my_refdate}-${my_reftod}"
icbc_in_progress="${my_status_dir}/icbc_in_progress.${my_refdate}-${my_reftod}"
rm -f -- "${icbc_status}"
printf 'valid_time=%s-%s\ncase=%s\nensemble_size=%s\nslurm_job_id=%s\nstarted_at=%s\n' "${my_refdate}" "${my_reftod}" "${my_casename}" "${my_ensnum}" "${SLURM_JOB_ID:-none}" "$(date '+%Y-%m-%d %H:%M:%S')" > "${icbc_in_progress}"

case "${my_runtype}" in
  Full-CPL|AMIP)
    ;;
  *)
    fail "my_runtype must be either Full-CPL or AMIP, got: ${my_runtype:-unset}"
    ;;
esac

for req in my_refeam_in my_refeam_ic my_refelm_in my_refrof_in my_refocn_in my_refice_in my_refcpl_in; do
  [[ -n "${!req:-}" ]] || fail "required variable is unset or empty: ${req}"
done


# Run options
REF_DATE=${my_refdate}
REF_TOD=${my_reftod}
REF_CASE=${my_refcase}
REF_DIR=${my_refdir}
REF_HOUR=`echo "${REF_TOD}" / 3600 | bc`
REF_HOUR=`printf "%02d" $REF_HOUR`
echo "valid time is $REF_DATE $REF_TOD (seconds) $REF_HOUR (hours)"

# Each member owns its raw CIME initial-condition archive.

# ==============================================================================
# machine-specific dereferencing
# suppress "rm" warnings if wildcard does not match anything
VERBOSE='-v'
MOVE='/usr/bin/mv'
COPY='/usr/bin/cp --preserve=timestamps'
LINK='/usr/bin/ln -fs'
LINKV=TRUE
LIST='/usr/bin/ls'
REMOVE='/usr/bin/rm -fr'
# ==============================================================================

prepare_mpas_restart() {
  local source_file="$1"
  local destination_file="$2"
  local tmp_file="${destination_file}.tmp.${SLURM_JOB_ID:-$$}"
  rm -f -- "${tmp_file}"
  ncks -O --hdr_pad=10000 "${source_file}" "${tmp_file}" || return 1
  ncrename -v xtime,xtime.orig "${tmp_file}" || return 1
  ncdump -h "${tmp_file}" >/dev/null 2>&1 || return 1
  touch --reference="${source_file}" "${tmp_file}" || return 1
  mv -f "${tmp_file}" "${destination_file}"
}
# Loop over members
for i in `seq 1 ${my_ensnum}`;do
  echo === Member ${i} ===
  ENSTR=EN`printf "%02d" ${i}`
  DART_CASE=${my_casename}.${ENSTR}
  MEMBER_ARCHIVE_DIR="${my_modeldir}/${ENSTR}/archive"
  CASE_ARCHIVE_DIR="${MEMBER_ARCHIVE_DIR}/rest/${REF_DATE}-${REF_TOD}"
  echo "Run Case: ${DART_CASE}"
  echo "Run Directory: ${CASE_ARCHIVE_DIR}"
  if [ ! -d "${CASE_ARCHIVE_DIR}" ];then
    mkdir -p "${CASE_ARCHIVE_DIR}"
  fi
  for scomp in "atm" "lnd" "rof" "ocn" "ice" "drv"; do 
     echo === E3SM component ${scomp} ===
     cd ${CASE_ARCHIVE_DIR}
     if [[ ${scomp} == "atm" ]]; then 
        smod="eam"
        REF_DATE_EXT=${REF_DATE}-${REF_TOD}
        ATM_INITIAL_FILENAME=${DART_CASE}.${smod}.i.${REF_DATE_EXT}.nc
        echo "process ${smod} ic: ${ATM_INITIAL_FILENAME}"
        if [ -f "${my_refeam_in}" ];then
          ${COPY} ${my_refeam_in} ${ATM_INITIAL_FILENAME} || exit 01
          echo "reference file: ${my_refeam_ic}"
          ncks -A -v PS,U,V,T,Q,CLDLIQ,CLDICE ${my_refeam_ic} ${ATM_INITIAL_FILENAME} || exit 02
          ATM_DATE=( `echo $REF_DATE_EXT | sed -e "s#-# #g"` )
          mdate=`echo ${REF_DATE} | sed "s/-//g"`
          mtods=${REF_TOD}
          ncap2 -O -s "date=${mdate};datesec=${mtods}" ${ATM_INITIAL_FILENAME} ${ATM_INITIAL_FILENAME} || exit 03
        else
          echo "ERROR: initial condition file not found: ${my_refeam_in}"
          exit 1
        fi
        ATM_REST_FILENAME="${DART_CASE}.${smod}.r.${REF_DATE_EXT}.nc"
        echo "${ATM_REST_FILENAME}"   >  rpointer.atm
        if [ ! -f "${ATM_REST_FILENAME}" ];then 
          touch ${ATM_REST_FILENAME}
        fi 
     elif [[ ${scomp} == "lnd" ]]; then
        smod="elm"
        REF_DATE_EXT=${REF_DATE}-${REF_TOD}
        LND_INITIAL_FILENAME=${DART_CASE}.${smod}.r.${REF_DATE_EXT}.nc
        echo "process ${smod} ic: ${LND_INITIAL_FILENAME}"
        if [ -f "${my_refelm_in}" ];then
          ${COPY} ${my_refelm_in} ${LND_INITIAL_FILENAME} || exit 04
        else
          echo "ERROR: initial condition file not found: ${my_refelm_in}"
          exit 1
        fi 
        echo "./${LND_INITIAL_FILENAME}"   >  rpointer.lnd
     elif [[ ${scomp} == "rof" ]]; then
        smod="mosart"
        REF_DATE_EXT=${REF_DATE}-${REF_TOD}
        ROF_INITIAL_FILENAME=${DART_CASE}.${smod}.r.${REF_DATE_EXT}.nc
        echo "process ${smod} ic: ${ROF_INITIAL_FILENAME}"
        if [ -f "${my_refrof_in}" ];then
          ${COPY} ${my_refrof_in} ${ROF_INITIAL_FILENAME} || exit 05
        else
          echo "ERROR: initial condition file not found: ${my_refrof_in}"
          exit 1
        fi
        echo "./${ROF_INITIAL_FILENAME}"   >  rpointer.rof
     elif [[ ${scomp} == "ocn"  &&  ${my_runtype} == "Full-CPL" ]]; then
        smod="mpaso"
        REF_DATE_EXT=${REF_DATE}_${REF_TOD}
        OCN_INITIAL_FILENAME=${DART_CASE}.${smod}.rst.${REF_DATE_EXT}.nc
        echo "process ${smod} ic: ${OCN_INITIAL_FILENAME}"
        if [ -f "${my_refocn_in}" ];then
          prepare_mpas_restart "${my_refocn_in}" "${OCN_INITIAL_FILENAME}" || exit 07
        else
          echo "ERROR: initial condition file not found: ${my_refocn_in}"
          exit 1
        fi
        echo "${REF_DATE}_`printf "%02d" ${REF_HOUR}`:00:00"  > rpointer.ocn
     elif [[ ${scomp} == "ocn"  && ${my_runtype} == "AMIP" ]]; then
        smod="docn"
        REF_DATE_EXT=${REF_DATE}_${REF_TOD}
        OCN_INITIAL_FILENAME1="${DART_CASE}.${smod}.r.${REF_DATE_EXT}.nc"
        OCN_INITIAL_FILENAME2="${DART_CASE}.${smod}.rs1.${REF_DATE_EXT}.bin"
        echo "process ${smod} ic: ${OCN_INITIAL_FILENAME1}"
        echo "process ${smod} ic: ${OCN_INITIAL_FILENAME2}"
        if [ ! -f "${OCN_INITIAL_FILENAME1}" ];then
          touch ${OCN_INITIAL_FILENAME1}
        fi
        if [ ! -f "${OCN_INITIAL_FILENAME2}" ];then
          touch ${OCN_INITIAL_FILENAME2}
        fi
        echo "${OCN_INITIAL_FILENAME1}"  >  rpointer.ocn
        echo "${OCN_INITIAL_FILENAME2}"  >> rpointer.ocn
     elif [[ ${scomp} == "ice" ]]; then
        smod="mpassi"
        REF_DATE_EXT=${REF_DATE}_${REF_TOD}
        ICE_INITIAL_FILENAME=${DART_CASE}.${smod}.rst.${REF_DATE_EXT}.nc
        echo "process ${smod} ic: ${ICE_INITIAL_FILENAME}"
        if [ -f "${my_refice_in}" ];then
          prepare_mpas_restart "${my_refice_in}" "${ICE_INITIAL_FILENAME}" || exit 10
        else
          echo "ERROR: initial condition file not found: ${my_refice_in}"
          exit 1
        fi
        echo "${REF_DATE}_`printf "%02d" ${REF_HOUR}`:00:00"  > rpointer.ice
     else 
        smod="cpl"
        REF_DATE_EXT=${REF_DATE}-${REF_TOD}
        CPL_INITIAL_FILENAME=${DART_CASE}.${smod}.r.${REF_DATE_EXT}.nc
        echo "process ${smod} ic: ${CPL_INITIAL_FILENAME}"
        if [ -f "${my_refcpl_in}" ];then
          ${COPY} ${my_refcpl_in} ${CPL_INITIAL_FILENAME} || exit 12
        else
          echo "ERROR: initial condition file not found: ${my_refcpl_in}"
          exit 1
        fi
        echo "${CPL_INITIAL_FILENAME}"  >  rpointer.drv
     fi
  done
done

# Validate every required component before declaring step 2 complete.
for i in `seq 1 ${my_ensnum}`; do
  ENSTR=EN`printf "%02d" ${i}`
  MEMBER_ARCHIVE_DIR="${my_modeldir}/${ENSTR}/archive"
  member_prefix="${MEMBER_ARCHIVE_DIR}/rest/${REF_DATE}-${REF_TOD}/${my_casename}.${ENSTR}"
  required_files=("${member_prefix}.eam.i.${REF_DATE}-${REF_TOD}.nc" "${member_prefix}.elm.r.${REF_DATE}-${REF_TOD}.nc" "${member_prefix}.mosart.r.${REF_DATE}-${REF_TOD}.nc" "${member_prefix}.mpassi.rst.${REF_DATE}_${REF_TOD}.nc" "${member_prefix}.cpl.r.${REF_DATE}-${REF_TOD}.nc")
  [[ "${my_runtype}" == "Full-CPL" ]] && required_files+=("${member_prefix}.mpaso.rst.${REF_DATE}_${REF_TOD}.nc")
  for component_file in "${required_files[@]}"; do
    if [ ! -s "${component_file}" ] || ! ncdump -h "${component_file}" >/dev/null 2>&1; then
      fail "invalid generated restart for ${ENSTR}: ${component_file}"
    fi
  done
done

# A successful regeneration makes any failed member-local Step 0b marker stale.
for i in $(seq 1 "${my_ensnum}"); do
  ENSTR=$(printf 'EN%02d' "${i}")
  MEMBER_ARCHIVE_DIR="${my_modeldir}/${ENSTR}/archive"
  perturb_marker="${MEMBER_ARCHIVE_DIR}/rest/${REF_DATE}-${REF_TOD}/.dart_perturb_in_progress"
  if [[ -e "${perturb_marker}" ]]; then
    rm -f -- "${perturb_marker}" || fail "could not clear stale perturbation marker for ${ENSTR}"
    echo "Cleared stale perturbation marker for ${ENSTR}"
  fi
done
icbc_tmp="${icbc_status}.tmp.${SLURM_JOB_ID:-$$}"
printf 'valid_time=%s-%s\ncase=%s\nensemble_size=%s\narchive_layout=%s\ndart_root=%s\nslurm_job_id=%s\ncompleted_at=%s\n' "${REF_DATE}" "${REF_TOD}" "${my_casename}" "${my_ensnum}" "per_member" "${my_dart_root}" "${SLURM_JOB_ID:-none}" "$(date '+%Y-%m-%d %H:%M:%S')" > "${icbc_tmp}"
mv -f "${icbc_tmp}" "${icbc_status}"
rm -f -- "${icbc_in_progress}"
echo "Wrote Step 0a completion record: ${icbc_status}"
echo "Next: sbatch 0b_perturb_initial_conditions.sh"
wait
 
exit
