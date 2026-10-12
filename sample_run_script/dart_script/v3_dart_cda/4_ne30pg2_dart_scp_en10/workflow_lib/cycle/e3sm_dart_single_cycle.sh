#!/bin/bash -el 
if [[ "${CYCLE_DRIVER_ACTIVE:-FALSE}" != "TRUE" ]]; then
  echo "ERROR: internal cycle worker; run 4_run_dart_e3sm_cycleda.sh" >&2
  exit 1
fi

#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------

#For cshell:
#limit stacksize unlimited
#limit datasize unlimited

#For bash 
#ulimit -s unlimited
#ulimit -d unlimited

#export SLURM_NNODES=20
#export SLURM_NTASKS=800

echo == Start of e3sm_dart_single_cycle.sh ==
date
echo ============================================

#source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh
#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.sh
my_wkdir=${PWD}

fail() {
  echo "ERROR: $*"
  exit 1
}

check_command() {
  local cmd="$1"
  command -v "${cmd}" >/dev/null 2>&1 || fail "required command not found: ${cmd}"
}

check_required_commands() {
  local cmd
  for cmd in awk bc cp date ex find grep mv ncdump patch rm sort srun stat tail; do
    check_command "${cmd}"
  done
}

compute_cycle_state_time() {
  local base_ymd="${my_e3sm_start_date:-${my_casedate:-}}"
  local base_tod="${my_e3sm_start_tod:-${my_casetod:-}}"
  local base_seconds base_hours base_minutes base_secs
  local offset_hours base_epoch target_epoch resolved_time

  [[ -n "${base_ymd}" ]] || fail "missing DA start date: set my_e3sm_start_date or my_casedate"
  [[ -n "${base_tod}" ]] || fail "missing DA start time: set my_e3sm_start_tod or my_casetod"
  [[ "${base_ymd}" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]] || fail "invalid DA start date: ${base_ymd}"
  [[ "${base_tod}" =~ ^[0-9]{5}$ ]] || fail "invalid DA start time: ${base_tod}"

  base_seconds=$((10#${base_tod}))
  (( base_seconds >= 0 && base_seconds < 86400 )) || fail "DA start time is out of range: ${base_tod}"
  base_hours=$((base_seconds / 3600))
  base_minutes=$(((base_seconds % 3600) / 60))
  base_secs=$((base_seconds % 60))
  offset_hours=$((DATA_ASSIMILATION_CYCLES * DATA_ASSIMILATION_WINDOW))
  base_epoch=$(date -d "${base_ymd} $(printf '%02d:%02d:%02d' "${base_hours}" "${base_minutes}" "${base_secs}")" +%s) || {
    fail "could not convert DA start time ${base_ymd}-${base_tod} to epoch seconds"
  }
  target_epoch=$((base_epoch + offset_hours * 3600))

  resolved_time=$(date -d "@${target_epoch}" +"%Y-%m-%d %H %M %S") || {
    fail "could not resolve DA state time from ${base_ymd}-${base_tod} and cycle ${DATA_ASSIMILATION_CYCLES}"
  }

  read -r CUR_YMD CUR_HOUR CUR_MINUTE CUR_SECOND <<< "${resolved_time}"
  CUR_TOD=$(printf '%05d' $((10#${CUR_HOUR} * 3600 + 10#${CUR_MINUTE} * 60 + 10#${CUR_SECOND})))
}

check_required_vars() {
  local var
  local missing=0
  for var in my_job_nnodes my_ensnum my_modeldir my_casename my_e3sm_cycle_hours my_eam_dart_cycle_hours my_elm_dart_cycle_hours my_eam_dart_run_dir my_elm_dart_run_dir my_eam_dart_end_date my_eam_dart_end_tod my_elm_dart_end_date my_elm_dart_end_tod my_e3sm_completed_cycles; do
    if [[ -z "${!var:-}" ]]; then
      echo "ERROR: required variable is unset or empty: ${var}"
      missing=1
    fi
  done

  if (( missing != 0 )); then
    exit 1
  fi
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

on_signal() {
  echo "ERROR: ${SCRIPT_NAME} received a termination signal"
  if declare -p setup_pids >/dev/null 2>&1; then
    for child_pid in "${setup_pids[@]}"; do kill -TERM "${child_pid}" 2>/dev/null || true; done
  fi
  if declare -p handoff_pids >/dev/null 2>&1; then
    for child_pid in "${handoff_pids[@]}"; do kill -TERM "${child_pid}" 2>/dev/null || true; done
  fi
  if declare -p member_pids >/dev/null 2>&1; then
    for child_pid in "${member_pids[@]}"; do kill -TERM "${child_pid}" 2>/dev/null || true; done
  fi
  wait 2>/dev/null || true
  if declare -F handoff_abort_on_signal >/dev/null 2>&1; then
    handoff_abort_on_signal || true
  fi
  if declare -F restore_retry_bfbflag >/dev/null 2>&1; then
    restore_retry_bfbflag || true
  fi
  exit 143
}
trap on_signal INT TERM
check_required_commands

cd ${my_wkdir}
source ./create_and_setup_case.sh
[[ -n "${my_dart_env_file:-}" ]] || fail "my_dart_env_file is unset"
[[ -r "${my_dart_env_file}" ]] || fail "configured DART environment is not readable: ${my_dart_env_file}"
echo "Using configured DART machine environment: ${my_dart_env_file}"
source "${my_dart_env_file}"
EAM_DART_DA="${my_eam_dart_da,,}"
ELM_DART_DA="${my_elm_dart_da,,}"
[[ "${EAM_DART_DA}" == "on" || "${EAM_DART_DA}" == "off" ]] || fail "my_eam_dart_da must be on or off"
[[ "${ELM_DART_DA}" == "on" || "${ELM_DART_DA}" == "off" ]] || fail "my_elm_dart_da must be on or off"
echo "Component DA modes: EAM=${EAM_DART_DA}, ELM=${ELM_DART_DA}"


LOG_DIR="${my_log_dir}"
mkdir -p "${LOG_DIR}"
check_required_vars

case "${my_runtype}" in
  Full-CPL|AMIP)
    ;;
  *)
    fail "my_runtype must be either Full-CPL or AMIP, got: ${my_runtype:-unset}"
    ;;
esac

RETRY_FORECAST_TIMEOUT_SEC="${RETRY_FORECAST_TIMEOUT_SEC:-${my_retry_forecast_timeout_sec:-0}}"
FORECAST_READY_FOR_DA="${my_forecast_ready_for_da:-FALSE}"
FORECAST_READY_FOR_DA="$(echo "${FORECAST_READY_FOR_DA}" | tr '[:lower:]' '[:upper:]')"
NODES_PER_MEMBER="${my_nodes_per_member:-2}"
WAIT_POLL_INTERVAL_SEC="${WAIT_POLL_INTERVAL_SEC:-${my_wait_poll_interval_sec:-10}}"
SKIP_COMPLETED_MEMBERS="${my_skip_completed_members:-TRUE}"
SKIP_COMPLETED_MEMBERS="$(echo "${SKIP_COMPLETED_MEMBERS}" | tr '[:lower:]' '[:upper:]')"
MAX_PARALLEL_SETUP="${MAX_PARALLEL_SETUP:-${my_max_parallel_setup:-8}}"

validate_positive_int "my_ensnum" "${my_ensnum}"
validate_positive_int "my_job_nnodes" "${my_job_nnodes}"
validate_positive_int "NODES_PER_MEMBER" "${NODES_PER_MEMBER}"
validate_positive_int "WAIT_POLL_INTERVAL_SEC" "${WAIT_POLL_INTERVAL_SEC}"
validate_positive_int "MAX_PARALLEL_SETUP" "${MAX_PARALLEL_SETUP}"
if ! [[ "${RETRY_FORECAST_TIMEOUT_SEC}" =~ ^[0-9]+$ ]]; then
  fail "RETRY_FORECAST_TIMEOUT_SEC must be a non-negative integer, got: ${RETRY_FORECAST_TIMEOUT_SEC}"
fi
if [[ "${SKIP_COMPLETED_MEMBERS}" != "TRUE" && "${SKIP_COMPLETED_MEMBERS}" != "FALSE" ]]; then
  fail "my_skip_completed_members must be TRUE or FALSE, got: ${SKIP_COMPLETED_MEMBERS}"
fi

DATA_ASSIMILATION_ATM=TRUE
DATA_ASSIMILATION_CYCLES=${my_e3sm_completed_cycles}
DATA_ASSIMILATION_WINDOW=${my_e3sm_cycle_hours}
validate_positive_int "my_e3sm_cycle_hours" "${my_e3sm_cycle_hours}"
validate_positive_int "my_eam_dart_cycle_hours" "${my_eam_dart_cycle_hours}"
validate_positive_int "my_elm_dart_cycle_hours" "${my_elm_dart_cycle_hours}"
(( my_eam_dart_cycle_hours % my_e3sm_cycle_hours == 0 )) || fail "my_eam_dart_cycle_hours must be an integer multiple of my_e3sm_cycle_hours"
(( my_elm_dart_cycle_hours % my_e3sm_cycle_hours == 0 )) || fail "my_elm_dart_cycle_hours must be an integer multiple of my_e3sm_cycle_hours"
TARGET_ELAPSED_HOURS=$(((DATA_ASSIMILATION_CYCLES + 1) * my_e3sm_cycle_hours))
EAM_DART_RUN="off"
ELM_DART_RUN="off"
if [[ "${EAM_DART_DA}" == "on" ]]; then
  EAM_DART_RUN="not_due"
  (( TARGET_ELAPSED_HOURS % my_eam_dart_cycle_hours != 0 )) || EAM_DART_RUN="on"
fi
if [[ "${ELM_DART_DA}" == "on" ]]; then
  ELM_DART_RUN="not_due"
  (( TARGET_ELAPSED_HOURS % my_elm_dart_cycle_hours != 0 )) || ELM_DART_RUN="on"
fi
echo "Component DA schedule at +${TARGET_ELAPSED_HOURS}h: EAM=${EAM_DART_RUN}, ELM=${ELM_DART_RUN}"

CASE_ROOT=${my_modeldir}/EN01/case_scripts
RUN_ROOT=${my_modeldir}/EN01/run

if [ -d "${CASE_ROOT}" ]; then
  cd ${CASE_ROOT}
  compute_cycle_state_time
  if [[ -x ./xmlquery ]]; then
    XML_CUR_YMD=`./xmlquery RUN_STARTDATE --value`
    XML_CUR_TOD=`./xmlquery START_TOD --value`
    XML_CUR_TOD=`printf %05d ${XML_CUR_TOD}`
    if [[ "${XML_CUR_YMD}" != "${CUR_YMD}" || "${XML_CUR_TOD}" != "${CUR_TOD}" ]]; then
      echo "WARNING: case XML time ${XML_CUR_YMD}-${XML_CUR_TOD} does not match cycle-derived DA state ${CUR_YMD}-${CUR_TOD}; using cycle-derived time"
    fi
  fi
else
  fail "Case directory does not exist: ${my_modelcase}"
fi

#determine the time for previous DA cycle 
CUR_DATE=( `echo ${CUR_YMD}-${CUR_TOD} | sed -e "s#-# #g"` )
CUR_YEAR=`echo "${CUR_DATE[0]}" | bc`
CUR_MONTH=`echo "${CUR_DATE[1]}" | bc`
CUR_DAY=`echo "${CUR_DATE[2]}" | bc`
CUR_HOUR=`echo "${CUR_DATE[3]}" / 3600 | bc`
CUR_SECONDS=`echo "${CUR_DATE[3]}" | bc`
echo "valid time for eam forecast cycle is $CUR_YEAR $CUR_MONTH $CUR_DAY $CUR_SECONDS (seconds)"

# Precompute target DA valid time for this cycle and the archive directory to check.
DA_TARGET_DATE=`date -d "$CUR_HOUR:00 ${CUR_YEAR}-${CUR_MONTH}-${CUR_DAY} + ${DATA_ASSIMILATION_WINDOW} hours" +"%Y-%m-%d %H"`
DA_TARGET_DATE=( `echo $DA_TARGET_DATE | sed -e "s#-# #g"` )
DA_TARGET_YEAR=`echo "${DA_TARGET_DATE[0]}" | bc`
DA_TARGET_MONTH=`echo "${DA_TARGET_DATE[1]}" | bc`
DA_TARGET_DAY=`echo "${DA_TARGET_DATE[2]}" | bc`
DA_TARGET_HOUR=`echo "${DA_TARGET_DATE[3]}" | bc`
DA_TARGET_SECONDS=`echo "${DA_TARGET_DATE[3]}" \* 3600 | bc`
DA_TARGET_YMD=`printf "%04d" ${DA_TARGET_YEAR}`-`printf "%02d" ${DA_TARGET_MONTH}`-`printf "%02d" ${DA_TARGET_DAY}`
DA_TARGET_TOD=`printf "%05d" ${DA_TARGET_SECONDS}`
TARGET_STAMP="${DA_TARGET_YMD//-/}${DA_TARGET_TOD}"
EAM_END_STAMP="${my_eam_dart_end_date//-/}${my_eam_dart_end_tod}"
ELM_END_STAMP="${my_elm_dart_end_date//-/}${my_elm_dart_end_tod}"
if [[ "${EAM_DART_DA}" == "on" ]] && (( 10#${TARGET_STAMP} > 10#${EAM_END_STAMP} )); then
  EAM_DART_RUN="ended"
fi
if [[ "${ELM_DART_DA}" == "on" ]] && (( 10#${TARGET_STAMP} > 10#${ELM_END_STAMP} )); then
  ELM_DART_RUN="ended"
fi
echo "Component DA window at ${DA_TARGET_YMD}-${DA_TARGET_TOD}: EAM=${EAM_DART_RUN}, ELM=${ELM_DART_RUN}"
DA_TRANSACTION_DIR="${my_dart_root}/transactions/${DA_TARGET_YMD}-${DA_TARGET_TOD}"
STALE_DART_MARKER="${DA_TRANSACTION_DIR}/.dart_filter_in_progress"
STALE_ELM_DART_MARKER="${DA_TRANSACTION_DIR}/.dart_elm_filter_in_progress"
mapfile -t STALE_DART_MARKERS < <(find "${DA_TRANSACTION_DIR}" -maxdepth 1 -type f -name '.dart_*_in_progress' -print 2>/dev/null || true)
if (( ${#STALE_DART_MARKERS[@]} > 0 )); then
  echo "WARNING: previous EAM or ELM DART attempt did not finish cleanly"
  echo "WARNING: forcing a rebuild of every forecast member before retrying DART"
  SKIP_COMPLETED_MEMBERS="FALSE"
fi
echo "target DA valid time is ${DA_TARGET_YMD}-${DA_TARGET_TOD}"
echo "skip completed members: ${SKIP_COMPLETED_MEMBERS}"

#function to modify namlist 
user_eam_nl() {
  local file="user_nl_eam"
  local ncdata_path="$1"
  local hist_freq="$2"
  local inithist_freq="${hist_freq}-HOURLY"

  if [[ ! -f "$file" ]]; then
    echo "[ERROR] Namelist file '$file' not found!"
    return 1
  fi

  ex "$file" <<EOF
g/^ *ncdata *=/s@=.*@= "${ncdata_path}"@
g/^ *inithist *=/s@=.*@= '${inithist_freq}'@
g/^ *inithist_all *=/s@=.*@= .true.@
wq
EOF
}

user_elm_nl() {
  local file="user_nl_elm"
  local finidat="$1"
  local yr_check="$2"
  local dynpft_check="$3"
  local fsurdat_check="$4"
  local pct_check="$5"

  ex "$file" <<EOF
g/^ *finidat *=/s@=.*@= "${finidat}"@
g/^ *check_finidat_year_consistency *=/s@=.*@= ${yr_check}@
g/^ *check_dynpft_consistency *=/s@=.*@= ${dynpft_check}@
g/^ *check_finidat_fsurdat_consistency *=/s@=.*@= ${fsurdat_check}@
g/^ *check_finidat_pct_consistency *=/s@=.*@= ${pct_check}@
wq
EOF
}

user_mosart_nl() {
  local file="user_nl_mosart"
  local finidat_rtm="$1"

  ex "$file" <<EOF
g/^ *finidat_rtm *=/s@=.*@= "${finidat_rtm}"@
wq
EOF
}

user_mpassi_nl() {
  local ymd="$1"
  local hour="$2"
  local ref_time="$3"
  local file="user_nl_mpassi"

  ex "$file" <<EOF
g/^ *config_start_time *=/s@=.*@= '${ymd}_${hour}'@
g/^ *config_calendar_type *=/s@=.*@= 'gregorian'@
g/^ *config_initial_condition_type *=/s@=.*@= 'restart'@
g/^ *config_do_restart *=/s@=.*@= .true.@
g/^ *config_restart_timestamp_name *=/s@=.*@= 'rpointer.ice'@
wq
EOF
}

user_mpaso_nl() {
  local ymd="$1"
  local hour="$2"
  local ref_time="$3"
  local file="user_nl_mpaso"

  ex "$file" <<EOF
g/^ *config_start_time *=/s@=.*@= '${ymd}_${hour}'@
g/^ *config_calendar_type *=/s@=.*@= 'gregorian'@
g/^ *config_do_restart *=/s@=.*@= .true.@
g/^ *config_output_reference_time *=/s@=.*@= '${ref_time}'@
g/^ *config_restart_timestamp_name *=/s@=.*@= 'rpointer.ocn'@
wq
EOF
}

# =====================================
# Customize MPAS stream files if needed
# =====================================
patch_mpaso_streams() {
echo
echo 'Modifying MPAS (OCEAN) streams files'
pushd ${1}

rline=`sed -n -e 12p streams.ocean`
rline=`echo ${rline}| sed "s/filename_template=//g"`

patch streams.ocean << EOF
--- streams.ocean
+++ streams.ocean
@@ -12,1 +12,1 @@
-                  filename_template=${rline}
+                  filename_template="${2}"
EOF

# copy to SourceMods
cp streams.ocean  ${3}/SourceMods/src.mpaso/

popd

}

patch_mpassi_streams() {
echo
echo 'Modifying MPAS streams files'
pushd ${1}

rlin1=`sed -n -e 11p streams.seaice`
rlin1=`echo ${rlin1}| sed "s/filename_template=//g"`

rlin2=`sed -n -e 38p streams.seaice`
rlin2=`echo ${rlin2}| sed "s/filename_template=//g"`

patch streams.seaice << EOF
--- streams.seaice    
+++ streams.seaice    
@@ -11,1 +11,1 @@
-                  filename_template=${rlin1}
+                  filename_template="${2}"
@@ -38,1 +38,1 @@
-                  filename_template=${rlin2}
+                  filename_template="${2}"
EOF

# copy to SourceMods
cp streams.seaice ${3}/SourceMods/src.mpassi/

popd

}


stage_file_atomic() {
  local source_file="$1"
  local destination_file="$2"
  local tmp_file="${destination_file}.tmp.${SLURM_JOB_ID:-$$}"

  [[ -s "${source_file}" ]] || fail "restart source is missing or empty: ${source_file}"
  if [[ "${source_file}" == *.nc ]] && command -v ncdump >/dev/null 2>&1; then
    ncdump -h "${source_file}" >/dev/null 2>&1 || fail "restart source is invalid NetCDF: ${source_file}"
  fi
  cp -p "${source_file}" "${tmp_file}"
  mv -f "${tmp_file}" "${destination_file}"
}

is_valid_da_eam_file() {
  local file="$1"
  local min_size_bytes=4096

  # Must exist and be non-empty.
  [[ -f "${file}" ]] || return 1
  [[ -s "${file}" ]] || return 1

  # Guard against tiny/truncated placeholders.
  if (( $(stat -c%s "${file}") < min_size_bytes )); then
    return 1
  fi

  # If ncks is available, enforce core EAM variable presence.
  if command -v ncks >/dev/null 2>&1; then
    ncks -m -v PS,U,V,T,Q "${file}" >/dev/null 2>&1 || return 1
  fi

  return 0
}

wait_for_member_batch() {
  local label="$1"
  local timeout_sec="${2:-0}"
  local failed=0
  local deadline=0
  local idx
  local pid
  local member
  local status

  if (( ${#member_pids[@]} == 0 )); then
    return 0
  fi

  echo "Waiting for ${label} batch: ${member_names[*]}"
  if (( timeout_sec > 0 )); then
    deadline=$(( $(date +%s) + timeout_sec ))
  fi

  for idx in "${!member_pids[@]}"; do
    pid="${member_pids[$idx]}"
    member="${member_names[$idx]}"
    if (( timeout_sec > 0 )); then
      if wait_for_pid_until_deadline "${pid}" "${deadline}"; then
        echo "${label} launcher exited with status 0 for ${member} (pid ${pid}); output validation follows"
      else
        status=$?
        if (( status == 124 )); then
          echo "ERROR: ${label} timed out for ${member} (pid ${pid})"
          kill "${pid}" 2>/dev/null || true
          wait "${pid}" 2>/dev/null || true
        fi
        echo "ERROR: ${label} failed for ${member} (pid ${pid}, status ${status})"
        failed=1
      fi
    elif wait "${pid}"; then
      echo "${label} launcher exited with status 0 for ${member} (pid ${pid}); output validation follows"
    else
      status=$?
      echo "ERROR: ${label} failed for ${member} (pid ${pid}, status ${status})"
      failed=1
    fi
  done

  member_pids=()
  member_names=()

  if (( failed != 0 )); then
    return 1
  fi
}

restore_retry_bfbflag() {
  local member case_dir
  local failed=0

  if ! declare -p retry_names >/dev/null 2>&1; then
    return 0
  fi
  for member in "${retry_names[@]}"; do
    case_dir="${CASE_ROOT/EN01/${member}}"
    echo "Restoring BFBFLAG=TRUE after retry for ${member}"
    if ! (cd "${case_dir}" && ./xmlchange BFBFLAG=TRUE); then
      echo "ERROR: could not restore BFBFLAG=TRUE for retry member ${member}"
      failed=1
    fi
  done
  retry_names=()
  (( failed == 0 ))
}

wait_for_pid_until_deadline() {
  local pid="$1"
  local deadline="$2"
  local interval="${WAIT_POLL_INTERVAL_SEC}"
  local now

  while kill -0 "${pid}" 2>/dev/null; do
    now=$(date +%s)
    if (( now >= deadline )); then
      return 124
    fi
    sleep "${interval}"
  done

  wait "${pid}"
}

repair_case_locked_files() {
  local case_dir="$1"
  local enstr="$2"
  local locked_dir="${case_dir}/LockedFiles"
  local locked_file active_file tmp_file

  [[ -d "${case_dir}" ]] || fail "missing case directory: ${case_dir}"
  [[ -d "${locked_dir}" ]] || return 0
  while IFS= read -r locked_file; do
    active_file="${case_dir}/$(basename "${locked_file}")"
    [[ -s "${active_file}" ]] || fail "cannot repair ${locked_file}: active file is missing or empty: ${active_file}"
    if [[ "${active_file}" == *.xml ]] && command -v xmllint >/dev/null 2>&1; then
      xmllint --noout "${active_file}" || fail "cannot repair ${locked_file}: active XML is invalid"
    fi
    tmp_file="${locked_file}.tmp.${SLURM_JOB_ID:-$$}"
    echo "Repairing zero-byte locked file for ${enstr}: $(basename "${locked_file}")"
    cp -p "${active_file}" "${tmp_file}"
    mv -f "${tmp_file}" "${locked_file}"
  done < <(find "${locked_dir}" -maxdepth 1 -type f -size 0 -print)
}

repair_invalid_locked_files() {
  local i enstr case_dir
  echo "Checking CIME LockedFiles for all ensemble members"
  for i in $(seq 1 "${my_ensnum}"); do
    enstr=$(printf "EN%02d" "${i}")
    case_dir="${CASE_ROOT/EN01/${enstr}}"
    repair_case_locked_files "${case_dir}" "${enstr}"
  done
}

preflight_cycle_inputs() {
  local dart_workdir="${my_eam_dart_code}/models/${my_eam_dart_model}/work"
  local yyyymm obs_file i enstr case_name ref_dir run_dir input_file
  local required_files=()

  echo "Preflighting all ensemble restart inputs"
  if [[ "${EAM_DART_RUN}" == "on" ]]; then
    echo "Preflighting EAM DART dependencies and observations"
    [[ -x "${dart_workdir}/filter" ]] || fail "missing DART filter executable: ${dart_workdir}/filter"
    [[ -s "${my_eam_filter_nml}" ]] || fail "missing workflow EAM filter namelist: ${my_eam_filter_nml}"
    [[ -s "${my_eam_topography_file}" ]] || fail "missing EAM topography file: ${my_eam_topography_file}"
    [[ -s "${my_eam_se_mapping_file}" || -s "${my_eam_cs_grid_file}" ]] || fail "both EAM mapping/grid files are missing"
    yyyymm=$(printf "%04d%02d" "${DA_TARGET_YEAR}" "${DA_TARGET_MONTH}")
    obs_file="${my_eam_dart_obsdir}/${yyyymm}_6H_CESM/obs_seq.${DA_TARGET_YMD}-${DA_TARGET_TOD}"
    [[ -s "${obs_file}" ]] || fail "missing target EAM observation sequence: ${obs_file}"
  else
    obs_file="EAM DA ${EAM_DART_RUN}"
  fi
  for i in $(seq 1 "${my_ensnum}"); do
    enstr=$(printf "EN%02d" "${i}")
    case_name="${my_casename}.${enstr}"
    ref_dir="${my_modeldir}/${enstr}/archive/rest/${CUR_YMD}-${CUR_TOD}"
    run_dir="${RUN_ROOT/EN01/${enstr}}"
    required_files=(
      "${ref_dir}/${case_name}.eam.i.${CUR_YMD}-${CUR_TOD}.nc"
      "${ref_dir}/${case_name}.elm.r.${CUR_YMD}-${CUR_TOD}.nc"
      "${ref_dir}/${case_name}.mosart.r.${CUR_YMD}-${CUR_TOD}.nc"
      "${ref_dir}/${case_name}.cpl.r.${CUR_YMD}-${CUR_TOD}.nc"
      "${ref_dir}/${case_name}.mpassi.rst.${CUR_YMD}_${CUR_TOD}.nc"
    )
    if [[ "${my_runtype}" != "AMIP" ]]; then
      required_files+=("${ref_dir}/${case_name}.mpaso.rst.${CUR_YMD}_${CUR_TOD}.nc")
    fi

    for input_file in "${required_files[@]}"; do
      [[ -s "${input_file}" ]] || fail "missing or empty cycle input for ${enstr}: ${input_file}"
      if [[ "${input_file}" == *.nc ]]; then
        ncdump -h "${input_file}" >/dev/null 2>&1 || fail "invalid NetCDF cycle input for ${enstr}: ${input_file}"
      fi
    done
  done
  echo "Preflight passed for all ${my_ensnum} members and ${obs_file}"
}

repair_invalid_locked_files
preflight_cycle_inputs

# First step: prepare members in bounded parallel batches, then launch forecasts.
if (( my_job_nnodes % NODES_PER_MEMBER != 0 )); then
  fail "my_job_nnodes (${my_job_nnodes}) must be divisible by ${NODES_PER_MEMBER}"
fi
BATCH_SIZE=$((my_job_nnodes / NODES_PER_MEMBER))
(( BATCH_SIZE >= 1 )) || fail "invalid forecast batch size: ${BATCH_SIZE}"

prepare_member() (
  local i="$1"
  local ENSTR CASE_NAME CASE_DIR RUN_DIR MEMBER_ARCHIVE_DIR REF_DIR
  local atm_in lnd_in rof_in ocn_in ice_in drv_in
  local atm_rst lnd_rst rof_rst drv_rst ice_rst ocn1_rst ocn2_rst ocn_rst
  local file fname OCN_DATE OCN_HOUR OCN_TIME

  ENSTR=$(printf "EN%02d" "${i}")
  CASE_NAME="${my_casename}.${ENSTR}"
  CASE_DIR="${CASE_ROOT/EN01/${ENSTR}}"
  RUN_DIR="${RUN_ROOT/EN01/${ENSTR}}"
  MEMBER_ARCHIVE_DIR="${my_modeldir}/${ENSTR}/archive"
  REF_DIR="${MEMBER_ARCHIVE_DIR}/rest/${CUR_YMD}-${CUR_TOD}"
  echo "Preparing ${ENSTR}"

  cd "${CASE_DIR}"
  ./xmlchange run_exe="--kill-on-bad-exit=1 --job-name=${CASE_NAME} \${EXEROOT}/e3sm.exe "
  ./xmlchange RUN_TYPE="hybrid"
  ./xmlchange CONTINUE_RUN=FALSE
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
  ./xmlchange DOUT_S=True
  ./xmlchange DOUT_S_ROOT="${MEMBER_ARCHIVE_DIR}"

  atm_in="${REF_DIR}/${CASE_NAME}.eam.i.${CUR_YMD}-${CUR_TOD}.nc"
  lnd_in="${REF_DIR}/${CASE_NAME}.elm.r.${CUR_YMD}-${CUR_TOD}.nc"
  rof_in="${REF_DIR}/${CASE_NAME}.mosart.r.${CUR_YMD}-${CUR_TOD}.nc"
  ocn_in="${REF_DIR}/${CASE_NAME}.mpaso.rst.${CUR_YMD}_${CUR_TOD}.nc"
  ice_in="${REF_DIR}/${CASE_NAME}.mpassi.rst.${CUR_YMD}_${CUR_TOD}.nc"
  drv_in="${REF_DIR}/${CASE_NAME}.cpl.r.${CUR_YMD}-${CUR_TOD}.nc"

  cd "${RUN_DIR}"
  rm -rvf rpointer*
  atm_rst="${CASE_NAME}.eam.r.${CUR_YMD}-${CUR_TOD}.nc"
  lnd_rst="./${CASE_NAME}.elm.r.${CUR_YMD}-${CUR_TOD}.nc"
  rof_rst="./${CASE_NAME}.mosart.r.${CUR_YMD}-${CUR_TOD}.nc"
  drv_rst="${CASE_NAME}.cpl.r.${CUR_YMD}-${CUR_TOD}.nc"
  ice_rst="${CASE_NAME}.mpassi.rst.${CUR_YMD}_${CUR_TOD}.nc"
  echo "${atm_rst}" > rpointer.atm
  echo "${lnd_rst}" > rpointer.lnd
  echo "${rof_rst}" > rpointer.rof
  echo "${drv_rst}" > rpointer.drv
  echo "${CUR_YMD}_$(printf "%02d" "${CUR_HOUR}"):00:00" > rpointer.ice

  if [[ "${my_runtype}" == "AMIP" ]]; then
    ocn1_rst="${CASE_NAME}.docn.r.${CUR_YMD}_${CUR_TOD}.nc"
    ocn2_rst="${CASE_NAME}.docn.rs1.${CUR_YMD}_${CUR_TOD}.bin"
    echo "${ocn1_rst}" > rpointer.ocn
    echo "${ocn2_rst}" >> rpointer.ocn
    if [[ -s "${REF_DIR}/${ocn1_rst}" ]]; then
      stage_file_atomic "${REF_DIR}/${ocn1_rst}" "${ocn1_rst}"
    fi
    if [[ -s "${REF_DIR}/${ocn2_rst}" ]]; then
      stage_file_atomic "${REF_DIR}/${ocn2_rst}" "${ocn2_rst}"
    fi
  else
    ocn_rst="${CASE_NAME}.mpaso.rst.${CUR_YMD}_${CUR_TOD}.nc"
    echo "${CUR_YMD}_$(printf "%02d" "${CUR_HOUR}"):00:00" > rpointer.ocn
    stage_file_atomic "${REF_DIR}/${ocn_rst}" "${ocn_rst}"
  fi
  # EAM starts from ncdata in hybrid mode, so its eam.i file is the required
  # atmospheric input; the eam.r placeholder is not staged.
  for file in ${lnd_rst} ${rof_rst} ${drv_rst} ${ice_rst}; do
    fname=$(basename "${file}")
    stage_file_atomic "${REF_DIR}/${fname}" "${file}"
  done

  OCN_DATE=${CUR_YMD}
  OCN_HOUR=$(echo "${CUR_TOD} / 3600" | bc)
  OCN_TIME="${OCN_DATE}_$(printf "%02d" "${OCN_HOUR}"):00:00"
  cd "${CASE_DIR}"
  user_eam_nl "${atm_in}" "${DATA_ASSIMILATION_WINDOW}"
  if (( DATA_ASSIMILATION_CYCLES == 0 )); then
    user_elm_nl "${lnd_in}" .false. .false. .false. .false.
  else
    user_elm_nl "${lnd_in}" .false. .false. .true. .true.
  fi
  user_mosart_nl "${rof_in}"
  user_mpassi_nl "${OCN_DATE}" "${OCN_HOUR}" "${OCN_TIME}"
  if [[ "${my_runtype}" == "Full-CPL" ]]; then
    user_mpaso_nl "${OCN_DATE}" "${OCN_HOUR}" "${OCN_TIME}"
  fi

  ./case.setup
  repair_case_locked_files "${CASE_DIR}" "${ENSTR}"
  patch_mpassi_streams "${RUN_DIR}" "${ice_in}" "${CASE_DIR}"
  if [[ "${my_runtype}" == "Full-CPL" ]]; then
    patch_mpaso_streams "${RUN_DIR}" "${ocn_in}" "${CASE_DIR}"
  fi
  echo "Preparation completed for ${ENSTR}"
)

capture_member_failure_diagnostics() {
  local enstr="$1"
  local case_dir="$2"
  local run_dir="$3"
  local diagnostic_log latest_e3sm_log member_archive

  member_archive="${my_modeldir}/${enstr}/archive"
  diagnostic_log="${LOG_DIR}/failure.${enstr}.${SLURM_JOB_ID:-$$}.cycle${DATA_ASSIMILATION_CYCLES}.log"
  latest_e3sm_log=$(find "${run_dir}" -maxdepth 1 -type f -name 'e3sm.log.*' -printf '%T@ %p\n' 2>/dev/null |
    sort -nr | awk 'NR == 1 {$1=""; sub(/^ /, ""); print; exit}')

  {
    echo "Failure diagnostics for ${enstr}"
    echo "captured_at=$(date '+%Y-%m-%d %H:%M:%S')"
    echo "slurm_job_id=${SLURM_JOB_ID:-none}"
    echo "case_dir=${case_dir}"
    echo "run_dir=${run_dir}"
    echo "expected_output=${member_archive}/rest/${DA_TARGET_YMD}-${DA_TARGET_TOD}/${my_casename}.${enstr}.eam.i.${DA_TARGET_YMD}-${DA_TARGET_TOD}.nc"
    echo "===== CaseStatus tail ====="
    if [[ -r "${case_dir}/CaseStatus" ]]; then
      tail -n 80 "${case_dir}/CaseStatus"
    else
      echo "CaseStatus unavailable"
    fi
    echo "===== latest E3SM log ====="
    if [[ -n "${latest_e3sm_log}" && -r "${latest_e3sm_log}" ]]; then
      echo "log=${latest_e3sm_log}"
      echo "===== matched failure lines ====="
      grep -nEi 'error|fatal|abort|input/output|segmentation|signal|killed|bad termination|mpi_abort' \
        "${latest_e3sm_log}" | tail -n 120 || true
      echo "===== E3SM log tail ====="
      tail -n 160 "${latest_e3sm_log}"
    else
      echo "No readable e3sm.log.* file found"
    fi
  } > "${diagnostic_log}" 2>&1

  echo "Captured ${enstr} failure diagnostics: ${diagnostic_log}"
}

wait_for_setup_batch() {
  local failed=0 idx pid member status
  for idx in "${!setup_pids[@]}"; do
    pid="${setup_pids[$idx]}"
    member="${setup_names[$idx]}"
    if wait "${pid}"; then
      echo "setup completed successfully for ${member} (pid ${pid})"
    else
      status=$?
      echo "ERROR: setup failed for ${member} (pid ${pid}, status ${status})"
      failed=1
    fi
  done
  setup_pids=()
  setup_names=()
  (( failed == 0 ))
}

FORECAST_CYCLE_START_EPOCH=$(date +%s)
setup_pids=()
setup_names=()
members_to_run=()
for i in $(seq 1 "${my_ensnum}"); do
  ENSTR=$(printf "EN%02d" "${i}")
  CASE_NAME="${my_casename}.${ENSTR}"
  MEMBER_ARCHIVE_DIR="${my_modeldir}/${ENSTR}/archive"
  FORECAST_TARGET_FILE="${MEMBER_ARCHIVE_DIR}/rest/${DA_TARGET_YMD}-${DA_TARGET_TOD}/${CASE_NAME}.eam.i.${DA_TARGET_YMD}-${DA_TARGET_TOD}.nc"
  RUN_DIR="${RUN_ROOT/EN01/${ENSTR}}"
  if [[ "${SKIP_COMPLETED_MEMBERS}" == "TRUE" ]] && is_valid_da_eam_file "${FORECAST_TARGET_FILE}"; then
    echo "Skipping ${ENSTR}: target DA file already exists: ${FORECAST_TARGET_FILE}"
    continue
  fi
  if [[ "${SKIP_COMPLETED_MEMBERS}" == "FALSE" ]]; then
    # A failed CIME run can still execute st_archive and return success. Remove
    # only this member's current-cycle atmospheric target so an old run file
    # cannot be copied back into the archive as apparent new forecast data.
    RUN_TARGET_FILE="${RUN_DIR}/${CASE_NAME}.eam.i.${DA_TARGET_YMD}-${DA_TARGET_TOD}.nc"
    if [[ -e "${FORECAST_TARGET_FILE}" || -L "${FORECAST_TARGET_FILE}" ]]; then
      echo "Cleaning previous archived forecast target for ${ENSTR}: ${FORECAST_TARGET_FILE}"
      rm -f -- "${FORECAST_TARGET_FILE}"
    fi
    if [[ -e "${RUN_TARGET_FILE}" || -L "${RUN_TARGET_FILE}" ]]; then
      echo "Cleaning previous run forecast target for ${ENSTR}: ${RUN_TARGET_FILE}"
      rm -f -- "${RUN_TARGET_FILE}"
    fi
  fi
  members_to_run+=("${ENSTR}")
  setup_log="${LOG_DIR}/setup.${ENSTR}.${SLURM_JOB_ID:-$$}.cycle${DATA_ASSIMILATION_CYCLES}.log"
  prepare_member "${i}" > "${setup_log}" 2>&1 &
  setup_pids+=("$!")
  setup_names+=("${ENSTR}")
  if (( ${#setup_pids[@]} >= MAX_PARALLEL_SETUP )); then
    wait_for_setup_batch || fail "member setup batch failed; aborting before forecasts"
  fi
done
wait_for_setup_batch || fail "member setup failed; aborting before forecasts"
echo "All ${#members_to_run[@]} requested members completed setup"

TOTAL_WAVES=$(((${#members_to_run[@]} + BATCH_SIZE - 1) / BATCH_SIZE))
CURRENT_WAVE=1
STATUS_LOG="${LOG_DIR}/forecast_wave_status.${SLURM_JOB_ID:-$$}.cycle${DATA_ASSIMILATION_CYCLES}.log"
echo "# wave,status,members,elapsed_sec,timestamp" > "${STATUS_LOG}"
echo "Forecast wave status log: ${STATUS_LOG}"
echo "Forecast plan: ${#members_to_run[@]} members, ${my_job_nnodes} nodes, ${NODES_PER_MEMBER} nodes/member"
echo "Running up to ${BATCH_SIZE} members per wave for ${TOTAL_WAVES} wave(s)"

member_pids=()
member_names=()
for ENSTR in "${members_to_run[@]}"; do
  CASE_DIR="${CASE_ROOT/EN01/${ENSTR}}"
  if (( ${#member_pids[@]} == 0 )); then
    WAVE_START_EPOCH=$(date +%s)
  fi
  cd "${CASE_DIR}"
  ./case.submit --no-batch > "${LOG_DIR}/e3sm.${ENSTR}.o${SLURM_JOB_ID}.cycle${DATA_ASSIMILATION_CYCLES}.log" 2>&1 &
  member_pids+=("$!")
  member_names+=("${ENSTR}")

  if (( ${#member_pids[@]} >= BATCH_SIZE )); then
    WAVE_MEMBERS=${#member_pids[@]}
    echo "Waiting for forecast wave ${CURRENT_WAVE}/${TOTAL_WAVES}"
    if ! wait_for_member_batch "forecast"; then
      WAVE_END_EPOCH=$(date +%s)
      WAVE_ELAPSED=$((WAVE_END_EPOCH - WAVE_START_EPOCH))
      echo "${CURRENT_WAVE},failed,${WAVE_MEMBERS},${WAVE_ELAPSED},$(date '+%Y-%m-%d %H:%M:%S')" >> "${STATUS_LOG}"
      fail "forecast batch failed; aborting before DA step"
    fi
    WAVE_END_EPOCH=$(date +%s)
    WAVE_ELAPSED=$((WAVE_END_EPOCH - WAVE_START_EPOCH))
    echo "${CURRENT_WAVE},success,${WAVE_MEMBERS},${WAVE_ELAPSED},$(date '+%Y-%m-%d %H:%M:%S')" >> "${STATUS_LOG}"
    echo "Forecast wave ${CURRENT_WAVE}/${TOTAL_WAVES} completed"
    CURRENT_WAVE=$((CURRENT_WAVE + 1))
  fi
done

if (( ${#member_pids[@]} > 0 )); then
  WAVE_MEMBERS=${#member_pids[@]}
  echo "Waiting for forecast wave ${CURRENT_WAVE}/${TOTAL_WAVES}"
  if ! wait_for_member_batch "forecast"; then
    WAVE_END_EPOCH=$(date +%s)
    WAVE_ELAPSED=$((WAVE_END_EPOCH - WAVE_START_EPOCH))
    echo "${CURRENT_WAVE},failed,${WAVE_MEMBERS},${WAVE_ELAPSED},$(date '+%Y-%m-%d %H:%M:%S')" >> "${STATUS_LOG}"
    fail "forecast failed; aborting before DA step"
  fi
  WAVE_END_EPOCH=$(date +%s)
  WAVE_ELAPSED=$((WAVE_END_EPOCH - WAVE_START_EPOCH))
  echo "${CURRENT_WAVE},success,${WAVE_MEMBERS},${WAVE_ELAPSED},$(date '+%Y-%m-%d %H:%M:%S')" >> "${STATUS_LOG}"
  echo "Forecast wave ${CURRENT_WAVE}/${TOTAL_WAVES} completed"
fi
echo "All forecast waves completed"

# Second Step: run eam-dart data assimilation 
#-----------------------------------------------------------------------------------#
# Alway assume we perform forecast first, then DART data assimilation. At begining 
# of each forecast, we link the DART modified IC/BC files to the run directory to 
# enable the forecast for cycling data assimilation 
#-----------------------------------------------------------------------------------# 
echo "valid time of current DA cycle is $DA_TARGET_YEAR $DA_TARGET_MONTH $DA_TARGET_DAY $DA_TARGET_SECONDS (seconds)"

# check if needed data is ready for data assimilation 
TMP_DATE=${DA_TARGET_YMD}
TMP_TOD=${DA_TARGET_TOD}
member_pids=()
member_names=()
retry_names=()

for i in `seq 1 $my_ensnum`;do
  ENSTR=EN`printf "%02d" ${i}`
  CASE_NAME=${my_casename}.${ENSTR}
  CASE_DIR=`echo ${CASE_ROOT} | sed "s/EN01/${ENSTR}/g"`
  RUN_DIR=`echo ${RUN_ROOT} | sed "s/EN01/${ENSTR}/g"`
  DA_REF_DIR="${my_modeldir}/${ENSTR}/archive/rest/${DA_TARGET_YMD}-${DA_TARGET_TOD}"
  atm_in="${DA_REF_DIR}/${CASE_NAME}.eam.i.${TMP_DATE}-${TMP_TOD}.nc"
  atm_in1="${RUN_DIR}/${CASE_NAME}.eam.i.${TMP_DATE}-${TMP_TOD}.nc"

  if ! is_valid_da_eam_file "${atm_in}"; then
     if is_valid_da_eam_file "${atm_in1}"; then 
       cd ${CASE_DIR}
       srun -N 1 -n 1 ./case.st_archive
	   else
	     cd ${CASE_DIR}
	     capture_member_failure_diagnostics "${ENSTR}" "${CASE_DIR}" "${RUN_DIR}"
	     echo "Restaging cycle inputs and setting BFBFLAG=FALSE for retry forecast member ${ENSTR}"
	     retry_names+=("${ENSTR}")
	     member_number=$((10#${ENSTR#EN}))
	     (
	       prepare_member "${member_number}"
	       cd "${CASE_DIR}"
	       ./xmlchange BFBFLAG=FALSE
	       ./case.submit --no-batch
	     ) > "${LOG_DIR}/e3sm.${ENSTR}.o${SLURM_JOB_ID}.cycle${DATA_ASSIMILATION_CYCLES}.retry.log" 2>&1 &
       PID=$!
       member_pids+=("${PID}")
       member_names+=("${ENSTR}")
       echo "run model again to obtain ${atm_in1}"
     fi 
  fi
done

retry_forecast_failed=0
if ! wait_for_member_batch "retry forecast" "${RETRY_FORECAST_TIMEOUT_SEC}"; then
  retry_forecast_failed=1
fi
restore_retry_bfbflag || fail "failed to restore BFBFLAG=TRUE after retry forecasts"
if (( retry_forecast_failed != 0 )); then
  fail "retry forecast failed; aborting before DA step"
fi

missing_forecast_members=()
for i in `seq 1 $my_ensnum`;do
  ENSTR=EN`printf "%02d" ${i}`
  CASE_NAME=${my_casename}.${ENSTR}
  CASE_DIR=`echo ${CASE_ROOT} | sed "s/EN01/${ENSTR}/g"`
  RUN_DIR=`echo ${RUN_ROOT} | sed "s/EN01/${ENSTR}/g"`
  DA_REF_DIR="${my_modeldir}/${ENSTR}/archive/rest/${DA_TARGET_YMD}-${DA_TARGET_TOD}"
  atm_in="${DA_REF_DIR}/${CASE_NAME}.eam.i.${TMP_DATE}-${TMP_TOD}.nc"
  atm_in1="${RUN_DIR}/${CASE_NAME}.eam.i.${TMP_DATE}-${TMP_TOD}.nc"

  if ! is_valid_da_eam_file "${atm_in}" && is_valid_da_eam_file "${atm_in1}"; then
    cd ${CASE_DIR}
    srun -N 1 -n 1 ./case.st_archive
  fi

  if ! is_valid_da_eam_file "${atm_in}"; then
    missing_forecast_members+=("${ENSTR}")
    retry_log="${LOG_DIR}/e3sm.${ENSTR}.o${SLURM_JOB_ID}.cycle${DATA_ASSIMILATION_CYCLES}.retry.log"
    echo "ERROR: forecast output missing after archive/retry for ${CASE_NAME}"
    if [[ -s "${retry_log}" ]]; then
      echo "ERROR: ${ENSTR} retry log: ${retry_log}"
      grep -E 'Exception from case_run:|ERROR: RUN FAIL:|See log file for details:' "${retry_log}" | tail -n 3 || true
    fi
    continue
  fi
  if [[ "${SKIP_COMPLETED_MEMBERS}" == "FALSE" ]] && (( $(stat -c %Y "${atm_in}") < FORECAST_CYCLE_START_EPOCH )); then
    fail "forecast output is stale and was not replaced in this cycle: ${atm_in}"
  fi
done

if (( ${#missing_forecast_members[@]} > 0 )); then
  fail "forecast output validation failed for ${#missing_forecast_members[@]} member(s): ${missing_forecast_members[*]}"
fi

FORECAST_READY_FOR_DA="TRUE"

# Stale component markers are cleared only after every forecast member has
# been rebuilt and validated. Each assimilation creates its own fresh marker.
for stale_marker in "${STALE_DART_MARKERS[@]}"; do
  if [[ -e "${stale_marker}" ]]; then
    rm -f "${stale_marker}" || fail "could not clear stale DART marker after forecast recovery: ${stale_marker}"
    echo "Cleared stale marker after validating the regenerated forecast ensemble: ${stale_marker}"
  fi
done

if [[ "${FORECAST_READY_FOR_DA}" != "TRUE" ]]; then
  fail "forecast readiness check failed; aborting before DA step"
fi

# Run only component analyses due at this valid time. When both are due they use the
# equal halves of the allocation; a single enabled component receives all nodes.
DART_YEAR=${DA_TARGET_YEAR}
DART_MONTH=${DA_TARGET_MONTH}
DART_DAY=${DA_TARGET_DAY}
DART_HOUR=${DA_TARGET_HOUR}
DART_SECONDS=${DA_TARGET_SECONDS}
DA_VALID_TIME="${DA_TARGET_YMD}-${DA_TARGET_TOD}"
ELM_USES_SEQUENTIAL_PRIOR="FALSE"
if [[ "${strongly_coupled_on,,}" == "on" && "${lnd_da_use_sequential_prior_post,,}" == ".true." ]]; then
  ELM_USES_SEQUENTIAL_PRIOR="TRUE"
fi
if [[ "${EAM_DART_RUN}" == "on" && "${ELM_DART_RUN}" == "on" ]]; then
  if [[ "${ELM_USES_SEQUENTIAL_PRIOR}" == "TRUE" ]]; then
    my_eam_dart_nnodes=${my_job_nnodes}
    my_elm_dart_nnodes=${my_job_nnodes}
  else
    my_eam_dart_nnodes=$((my_job_nnodes / 2))
    my_elm_dart_nnodes=$((my_job_nnodes / 2))
  fi
elif [[ "${EAM_DART_RUN}" == "on" ]]; then
  my_eam_dart_nnodes=${my_job_nnodes}
  my_elm_dart_nnodes=0
elif [[ "${ELM_DART_RUN}" == "on" ]]; then
  my_eam_dart_nnodes=0
  my_elm_dart_nnodes=${my_job_nnodes}
else
  my_eam_dart_nnodes=0
  my_elm_dart_nnodes=0
fi
echo "DA valid time: ${DA_VALID_TIME}; EAM=${EAM_DART_RUN} (${my_eam_dart_nnodes} nodes); ELM=${ELM_DART_RUN} (${my_elm_dart_nnodes} nodes)"

BIDIRECTIONAL_RUN="FALSE"
if [[ "${dart_coupling_mode:-legacy}" == "bidirectional" ]]; then
  [[ "${EAM_DART_RUN}" == "on" && "${ELM_DART_RUN}" == "on" ]] || fail "bidirectional mode requires both component passes at ${DA_VALID_TIME}"
  BIDIRECTIONAL_RUN="TRUE"
  rm -f -- "${my_status_dir}/eam_assim_complete.${DA_VALID_TIME}" "${my_status_dir}/elm_assim_complete.${DA_VALID_TIME}" "${my_status_dir}/coupled_assim_complete.${DA_VALID_TIME}"
  ( set -Eeuo pipefail; cd "${my_wkdir}"; . "${my_workflow_lib}/cycle/bidirectional_dart_assimilation.sh" ) \
    || fail "bidirectional assimilation failed; cycle handoff is blocked"
fi

if [[ "${BIDIRECTIONAL_RUN}" != "TRUE" && "${ELM_DART_RUN}" == "on" && "${ELM_USES_SEQUENTIAL_PRIOR}" != "TRUE" ]]; then
  ( set -Eeuo pipefail; cd "${my_wkdir}"; ELM_DART_PREFLIGHT_ONLY=TRUE . "${my_workflow_lib}/cycle/elm_dart_assimilation.sh" ) || fail "ELM DART preflight failed; no component assimilation was started"
fi

if [[ "${BIDIRECTIONAL_RUN}" != "TRUE" && "${ELM_DART_RUN}" == "on" && "${ELM_USES_SEQUENTIAL_PRIOR}" == "TRUE" ]]; then
  ( set -Eeuo pipefail; cd "${my_wkdir}"; ELM_DART_PREFLIGHT_ONLY=TRUE ELM_DART_PREFLIGHT_STATIC_ONLY=TRUE . "${my_workflow_lib}/cycle/elm_dart_assimilation.sh" ) || fail "static ELM DART preflight failed; no component assimilation was started"
fi

[[ "${BIDIRECTIONAL_RUN}" == "TRUE" ]] || rm -f -- "${my_status_dir}/eam_assim_complete.${DA_VALID_TIME}" "${my_status_dir}/elm_assim_complete.${DA_VALID_TIME}" "${my_status_dir}/coupled_assim_complete.${DA_VALID_TIME}"
eam_da_log="${LOG_DIR}/assim.eam.${SLURM_JOB_ID:-$$}.cycle${DATA_ASSIMILATION_CYCLES}.log"
elm_da_log="${LOG_DIR}/assim.elm.${SLURM_JOB_ID:-$$}.cycle${DATA_ASSIMILATION_CYCLES}.log"
eam_da_pid=""
elm_da_pid=""
eam_da_status=0
elm_da_status=0
if [[ "${BIDIRECTIONAL_RUN}" != "TRUE" && "${ELM_DART_RUN}" == "on" && "${ELM_USES_SEQUENTIAL_PRIOR}" == "TRUE" ]]; then
  mkdir -p "${DA_TRANSACTION_DIR}" || fail "could not create DA transaction directory: ${DA_TRANSACTION_DIR}"
  [[ ! -e "${STALE_ELM_DART_MARKER}" ]] || fail "stale ELM transaction marker requires a full forecast rebuild: ${STALE_ELM_DART_MARKER}"
  printf 'cycle=%s\nvalid_time=%s\nslurm_job_id=%s\nstarted_at=%s\nphase=waiting_for_eam\n' \
    "${DATA_ASSIMILATION_CYCLES}" "${DA_VALID_TIME}" "${SLURM_JOB_ID:-none}" "$(date '+%F %T')" > "${STALE_ELM_DART_MARKER}" \
    || fail "could not create sequential ELM transaction marker"
fi
if [[ "${BIDIRECTIONAL_RUN}" != "TRUE" && "${EAM_DART_RUN}" == "on" ]]; then
  ( set -Eeuo pipefail; cd "${my_wkdir}"; . "${my_workflow_lib}/cycle/eam_dart_assimilation.sh" ) > "${eam_da_log}" 2>&1 &
  eam_da_pid=$!
else
  echo "EAM DART assimilation is ${EAM_DART_RUN}"
fi
if [[ "${BIDIRECTIONAL_RUN}" != "TRUE" && "${ELM_DART_RUN}" == "on" && "${ELM_USES_SEQUENTIAL_PRIOR}" == "TRUE" ]]; then
  [[ "${EAM_DART_RUN}" == "on" ]] || fail "ELM sequential-prior mode requires EAM DA at the same valid time"
  wait "${eam_da_pid}" || eam_da_status=$?
  eam_da_pid=""
  if (( eam_da_status != 0 )); then
    echo "ERROR: EAM assimilation failed before the sequential ELM pass (status ${eam_da_status})"
    tail -n 40 "${eam_da_log}" >&2 || true
    fail "ELM sequential-prior pass blocked because EAM did not complete"
  fi
  ( set -Eeuo pipefail; cd "${my_wkdir}"; ELM_DART_PREFLIGHT_ONLY=TRUE . "${my_workflow_lib}/cycle/elm_dart_assimilation.sh" ) || fail "ELM sequential-prior preflight failed after EAM completed"
fi
if [[ "${BIDIRECTIONAL_RUN}" != "TRUE" && "${ELM_DART_RUN}" == "on" ]]; then
  ( set -Eeuo pipefail; cd "${my_wkdir}"; ELM_DART_TRANSACTION_STARTED="${ELM_USES_SEQUENTIAL_PRIOR}" . "${my_workflow_lib}/cycle/elm_dart_assimilation.sh" ) > "${elm_da_log}" 2>&1 &
  elm_da_pid=$!
else
  echo "ELM DART assimilation is ${ELM_DART_RUN}"
fi
[[ -z "${eam_da_pid}" ]] || wait "${eam_da_pid}" || eam_da_status=$?
[[ -z "${elm_da_pid}" ]] || wait "${elm_da_pid}" || elm_da_status=$?
if (( eam_da_status != 0 || elm_da_status != 0 )); then
  echo "ERROR: component assimilation failed: EAM=${eam_da_status}, ELM=${elm_da_status}"
  [[ "${EAM_DART_RUN}" != "on" ]] || { echo "ERROR: EAM log: ${eam_da_log}"; tail -n 40 "${eam_da_log}" >&2 || true; }
  [[ "${ELM_DART_RUN}" != "on" ]] || { echo "ERROR: ELM log: ${elm_da_log}"; tail -n 40 "${elm_da_log}" >&2 || true; }
  fail "cycle handoff blocked because an enabled component analysis failed"
fi

write_da_record() {
  local component="$1" mode="$2"
  local record="${my_status_dir}/${component}_assim_complete.${DA_VALID_TIME}"
  printf 'valid_time=%s\ncase=%s\nensemble_size=%s\ncycle=%s\nda_mode=%s\nslurm_job_id=%s\ncompleted_at=%s\n' "${DA_VALID_TIME}" "${my_casename}" "${my_ensnum}" "${DATA_ASSIMILATION_CYCLES}" "${mode}" "${SLURM_JOB_ID:-none}" "$(date '+%F %T')" > "${record}.tmp.${SLURM_JOB_ID:-$$}"
  mv -f "${record}.tmp.${SLURM_JOB_ID:-$$}" "${record}"
}
[[ "${EAM_DART_RUN}" == "on" ]] || write_da_record eam "${EAM_DART_RUN}"
[[ "${EAM_DART_RUN}" != "on" ]] || write_da_record eam on
[[ "${ELM_DART_RUN}" == "on" ]] || write_da_record elm "${ELM_DART_RUN}"
if [[ "${ELM_DART_RUN}" == "on" ]]; then
  if [[ "${BIDIRECTIONAL_RUN}" == "TRUE" ]]; then
    [[ -s "${my_status_dir}/elm_assim_complete.${DA_VALID_TIME}.elm_atmosphere_cross" ]] || fail "bidirectional ELM cross-pass completion record is missing"
  else
    [[ -s "${my_status_dir}/elm_assim_complete.${DA_VALID_TIME}" ]] || fail "ELM completion record is missing"
  fi
fi
coupled_record="${my_status_dir}/coupled_assim_complete.${DA_VALID_TIME}"
printf 'valid_time=%s\ncase=%s\nensemble_size=%s\ncycle=%s\neam_da=%s\nelm_da=%s\nslurm_job_id=%s\ncompleted_at=%s\n' "${DA_VALID_TIME}" "${my_casename}" "${my_ensnum}" "${DATA_ASSIMILATION_CYCLES}" "${EAM_DART_RUN}" "${ELM_DART_RUN}" "${SLURM_JOB_ID:-none}" "$(date '+%F %T')" > "${coupled_record}.tmp.${SLURM_JOB_ID:-$$}"
mv -f "${coupled_record}.tmp.${SLURM_JOB_ID:-$$}" "${coupled_record}"
echo "Enabled component analyses completed; coupled handoff is permitted"

DATA_ASSIMILATION_CYCLES=$((DATA_ASSIMILATION_CYCLES+1))

# Post steps
cd ${my_wkdir}
. "${my_workflow_lib}/cycle/e3sm_bundle_cycling.sh"

# That's all folks!
sleep 10

echo ===== End of e3sm_dart_single_cycle.sh =====
date
echo =====================================
