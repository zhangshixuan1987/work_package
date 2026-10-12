#!/bin/bash
#------------------------------------------------------------------------------
# Multi-cycle Slurm driver. With my_cycles_per_job=1 this preserves the
# original one-cycle-per-job behavior.
#------------------------------------------------------------------------------
#SBATCH --account=esmd
#SBATCH --time=6:00:00
#SBATCH --partition=slurm
#SBATCH --job-name=e3sm_dart_ensda_cyc
#SBATCH --nodes=80
#SBATCH --output=runtmp/logs/e3sm_dart_ensda_cyc.%j
#SBATCH --exclusive
#SBATCH --no-kill
#SBATCH --requeue

set -Eeuo pipefail

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

validate_my_eam_cycle_overrides() {
  local key stamp parameter value ymd tod
  for key in "${!my_eam_cycle_overrides[@]}"; do
    if [[ ! "${key}" =~ ^([0-9]{4}-[0-9]{2}-[0-9]{2}-[0-9]{5}):(localization_cutoff|inflation_damping|no_obs_assim_above_level)$ ]]; then
      echo "invalid my_eam_cycle_overrides key: ${key}" >&2
      return 1
    fi
    stamp=${BASH_REMATCH[1]}
    parameter=${BASH_REMATCH[2]}
    ymd=${stamp:0:10}
    tod=${stamp:11:5}
    date -d "${ymd}" +%F >/dev/null 2>&1 || { echo "invalid override date: ${key}" >&2; return 1; }
    (( 10#${tod} < 86400 )) || { echo "override time is outside 00000-86399: ${key}" >&2; return 1; }
    value=${my_eam_cycle_overrides[${key}]}
    case "${parameter}" in
      localization_cutoff)
        [[ "${value}" =~ ^[0-9]+([.][0-9]+)?$ ]] && awk -v v="${value}" 'BEGIN {exit !(v > 0)}' || { echo "invalid localization cutoff for ${stamp}: ${value}" >&2; return 1; }
        ;;
      inflation_damping)
        [[ "${value}" =~ ^[0-9]+([.][0-9]+)?$ ]] && awk -v v="${value}" 'BEGIN {exit !(v >= 0 && v <= 1)}' || { echo "invalid inflation damping for ${stamp}: ${value}" >&2; return 1; }
        ;;
      no_obs_assim_above_level)
        [[ "${value}" =~ ^[1-9][0-9]*$ ]] && (( value <= 72 )) || { echo "invalid model-top cutoff level for ${stamp}: ${value}" >&2; return 1; }
        ;;
    esac
  done
}

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  [[ -n "${SLURM_SUBMIT_DIR:-}" ]] || fail "SLURM_SUBMIT_DIR is unavailable"
  my_wkdir=$(readlink -f "${SLURM_SUBMIT_DIR}") || fail "could not resolve SLURM_SUBMIT_DIR"
else
  SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || fail "could not resolve script path"
  my_wkdir=$(dirname "${SCRIPT_PATH}")
fi
CONFIG_FILE="${my_wkdir}/create_and_setup_case.sh"
SCRIPT_PATH="${my_wkdir}/4_run_dart_e3sm_cycleda.sh"
[[ -r "${SCRIPT_PATH}" ]] || fail "missing persistent cycle driver: ${SCRIPT_PATH}"
[[ -r "${CONFIG_FILE}" ]] || fail "missing cycle configuration: ${CONFIG_FILE}"
source "${CONFIG_FILE}"
[[ -n "${my_dart_env_file:-}" ]] || fail "my_dart_env_file is unset"
[[ -r "${my_dart_env_file}" ]] || fail "configured DART environment is not readable: ${my_dart_env_file}"
mkdir -p "${my_log_dir}" "${my_status_dir}" "${my_lock_dir}" "${my_handoff_dir}"

for cmd in awk basename cksum date dirname flock mkdir mv ncdump readlink sbatch squeue; do
  command -v "${cmd}" >/dev/null 2>&1 || fail "required command not found: ${cmd}"
done
validate_my_eam_cycle_overrides || fail "invalid my_eam_cycle_overrides configuration"
exec 9>"${my_lock_dir}/e3sm_dart_cycle.lock"
flock -n 9 || fail "another cycle driver is already running for ${my_wkdir}"
marker_field() { awk -F= -v key="$2" '$1 == key {sub(/^[^=]*=/, ""); print; found=1; exit} END {if (!found) exit 1}' "$1"; }
validate_marker() {
  local marker="$1" expected_time="$2" expected_case="$3" expected_size="$4"
  local actual_time actual_case actual_size actual_layout actual_dart_root
  [[ -s "${marker}" ]] || fail "missing upstream completion record: ${marker}"
  actual_time=$(marker_field "${marker}" valid_time) || fail "completion record lacks valid_time: ${marker}"
  actual_case=$(marker_field "${marker}" case) || fail "completion record lacks case: ${marker}"
  actual_size=$(marker_field "${marker}" ensemble_size) || fail "completion record lacks ensemble_size: ${marker}"
  actual_layout=$(marker_field "${marker}" archive_layout) || fail "completion record lacks archive_layout: ${marker}; rerun the upstream step"
  actual_dart_root=$(marker_field "${marker}" dart_root) || fail "completion record lacks dart_root: ${marker}; rerun the upstream step"
  [[ "${actual_time}" == "${expected_time}" && "${actual_case}" == "${expected_case}" && "${actual_size}" == "${expected_size}" && "${actual_layout}" == "per_member" && "${actual_dart_root}" == "${my_dart_root}" ]] || fail "upstream completion record does not match configured time, case, ensemble size, archive layout, or DART root: ${marker}"
}
validate_marker "${my_status_dir}/perturb_complete.${my_refdate}-${my_reftod}" "${my_refdate}-${my_reftod}" "${my_casename}" "${my_ensnum}"
validate_positive_int() {
  local name="$1"
  local value="$2"
  [[ "${value}" =~ ^[1-9][0-9]*$ ]] || fail "${name} must be a positive integer, got: ${value}"
}

validate_on_off() {
  local name="$1" value="${2,,}"
  [[ "${value}" == "on" || "${value}" == "off" ]] || fail "${name} must be on or off, got: ${2}"
}

validate_true_false() {
  local name="$1" value="${2^^}"
  [[ "${value}" == "TRUE" || "${value}" == "FALSE" ]] || fail "${name} must be TRUE or FALSE, got: ${2}"
}

remaining_allocation_seconds() {
  local value days=0 hours minutes seconds

  [[ -n "${SLURM_JOB_ID:-}" ]] || return 1
  value=$(squeue -h -j "${SLURM_JOB_ID}" -o '%L') || return 1
  value=${value//[[:space:]]/}
  [[ -n "${value}" && "${value}" != "UNLIMITED" && "${value}" != "NOT_SET" ]] || return 1

  if [[ "${value}" == *-* ]]; then
    days=${value%%-*}
    value=${value#*-}
  fi
  IFS=: read -r hours minutes seconds <<< "${value}"
  [[ "${days}" =~ ^[0-9]+$ && "${hours}" =~ ^[0-9]+$ && "${minutes}" =~ ^[0-9]+$ && "${seconds}" =~ ^[0-9]+$ ]] || return 1
  echo $((10#${days} * 86400 + 10#${hours} * 3600 + 10#${minutes} * 60 + 10#${seconds}))
}

expected_state_time() {
  local start_seconds start_hour start_minute start_second start_epoch state_epoch
  start_seconds=$((10#${my_e3sm_start_tod}))
  start_hour=$((start_seconds / 3600))
  start_minute=$(((start_seconds % 3600) / 60))
  start_second=$((start_seconds % 60))
  start_epoch=$(date -d "${my_e3sm_start_date} $(printf '%02d:%02d:%02d' "${start_hour}" "${start_minute}" "${start_second}")" +%s) || return 1
  state_epoch=$((start_epoch + my_e3sm_completed_cycles * my_e3sm_cycle_hours * 3600))
  date -d "@${state_epoch}" '+%Y-%m-%d %H %M %S' | awk '{printf "%s-%05d\n", $1, $2 * 3600 + $3 * 60 + $4}'
}

validate_runtime_state() {
  local expected marker marker_time marker_cycle case_dir member member_name member_archive
  local current_date current_tod actual common_actual="" expected_stamp actual_stamp
  local later_marker later_cycle rollback=FALSE forward_recovery=FALSE marker_tmp case_name file
  local dart_evidence dart_mean
  local -a required_files

  expected=$(expected_state_time) || fail "could not derive current state from the configured timeline"

  # First establish that every member has the configured archive and that all
  # live case clocks agree.
  for ((member=1; member<=my_ensnum; member++)); do
    member_name=$(printf 'EN%02d' "${member}")
    member_archive="${my_modeldir}/${member_name}/archive"
    [[ -d "${member_archive}/rest/${expected}" ]] || fail "missing restart archive for ${member_name} at configured state: ${member_archive}/rest/${expected}"
    case_dir="${my_modeldir}/${member_name}/case_scripts"
    [[ -x "${case_dir}/xmlquery" ]] || fail "missing xmlquery: ${case_dir}/xmlquery"
    current_date=$(cd "${case_dir}" && ./xmlquery RUN_STARTDATE --value) || fail "could not query RUN_STARTDATE from ${case_dir}"
    current_tod=$(cd "${case_dir}" && ./xmlquery START_TOD --value) || fail "could not query START_TOD from ${case_dir}"
    [[ "${current_tod}" =~ ^[0-9]+$ ]] || fail "invalid case XML time in ${case_dir}: ${current_date}-${current_tod}"
    actual=$(printf '%s-%05d' "${current_date}" "$((10#${current_tod}))")
    [[ -z "${common_actual}" || "${actual}" == "${common_actual}" ]] || fail "mixed case XML states: ${member_name} is ${actual}, earlier members are ${common_actual}"
    common_actual="${actual}"
  done

  expected_stamp=${expected//-/}
  actual_stamp=${common_actual//-/}
  if (( 10#${actual_stamp} < 10#${expected_stamp} )); then
    [[ "${my_allow_forward_xml_recovery^^}" == "TRUE" ]] || fail "case XML state ${common_actual} is older than configured state ${expected}; set my_allow_forward_xml_recovery=TRUE only for a deliberate, validated recovery"
    (( my_e3sm_completed_cycles > 0 )) || fail "forward XML recovery cannot target the pre-cycle state"
    forward_recovery=TRUE
  elif (( 10#${actual_stamp} > 10#${expected_stamp} )); then
    rollback=TRUE
  fi

  if (( my_e3sm_completed_cycles == 0 )); then
    if [[ "${rollback}" == "TRUE" ]]; then
      later_marker="${my_status_dir}/cycle_complete.${common_actual}"
      validate_marker "${later_marker}" "${common_actual}" "${my_casename}" "${my_ensnum}"
      later_cycle=$(marker_field "${later_marker}" cycle) || fail "later completion record lacks cycle: ${later_marker}"
      [[ "${later_cycle}" =~ ^[0-9]+$ ]] && (( 10#${later_cycle} > 0 )) \
        || fail "later marker does not prove progression beyond the configured perturbation state: ${later_marker}"
      echo "Validated configured rollback to the pre-cycle perturbation state ${expected}"
    fi
    [[ "${forward_recovery}" == "FALSE" ]] || fail "forward XML recovery cannot target the pre-cycle perturbation state"
    marker="${my_status_dir}/perturb_complete.${expected}"
    validate_marker "${marker}" "${expected}" "${my_casename}" "${my_ensnum}"
  else
    marker="${my_status_dir}/cycle_complete.${expected}"
    if [[ "${rollback}" == "TRUE" ]]; then
      later_marker="${my_status_dir}/cycle_complete.${common_actual}"
      validate_marker "${later_marker}" "${common_actual}" "${my_casename}" "${my_ensnum}"
      later_cycle=$(marker_field "${later_marker}" cycle) || fail "later completion record lacks cycle: ${later_marker}"
      [[ "${later_cycle}" =~ ^[0-9]+$ ]] && (( 10#${later_cycle} > my_e3sm_completed_cycles )) || fail "later marker does not prove progression beyond configured cycle ${my_e3sm_completed_cycles}: ${later_marker}"
    fi

    if [[ "${rollback}" == "TRUE" || "${forward_recovery}" == "TRUE" ]]; then
      # Validate every component input before creating recovery metadata.
      for ((member=1; member<=my_ensnum; member++)); do
        member_name=$(printf 'EN%02d' "${member}")
        member_archive="${my_modeldir}/${member_name}/archive"
        case_name="${my_casename}.${member_name}"
        required_files=(
          "${member_archive}/rest/${expected}/${case_name}.eam.i.${expected}.nc"
          "${member_archive}/rest/${expected}/${case_name}.elm.r.${expected}.nc"
          "${member_archive}/rest/${expected}/${case_name}.mosart.r.${expected}.nc"
          "${member_archive}/rest/${expected}/${case_name}.cpl.r.${expected}.nc"
          "${member_archive}/rest/${expected}/${case_name}.mpassi.rst.${expected:0:10}_${expected:11:5}.nc"
        )
        [[ "${my_runtype}" == "AMIP" ]] || required_files+=("${member_archive}/rest/${expected}/${case_name}.mpaso.rst.${expected:0:10}_${expected:11:5}.nc")
        for file in "${required_files[@]}"; do
          [[ -s "${file}" ]] || fail "missing recovery restart file: ${file}"
          ncdump -h "${file}" >/dev/null 2>&1 || fail "invalid recovery NetCDF file: ${file}"
        done
      done

      if [[ "${forward_recovery}" == "TRUE" && "${my_eam_dart_da,,}" == "on" ]]; then
        dart_evidence="${my_eam_dart_run_dir}/${expected}/${my_casename}.dart.e.eam_obs_seq_final.${expected}"
        dart_mean="${my_eam_dart_run_dir}/${expected}/${my_casename}.dart.i.eam_output_mean.${expected}.nc"
        [[ -s "${dart_evidence}" ]] || fail "missing archived EAM DART evidence for forward recovery: ${dart_evidence}"
        [[ -s "${dart_mean}" ]] || fail "missing archived EAM DART analysis mean for forward recovery: ${dart_mean}"
        ncdump -h "${dart_mean}" >/dev/null 2>&1 || fail "invalid EAM DART analysis mean for forward recovery: ${dart_mean}"
      fi

      if [[ ! -s "${marker}" ]]; then
        marker_tmp="${marker}.tmp.${SLURM_JOB_ID:-$$}"
        if [[ "${forward_recovery}" == "TRUE" ]]; then
          printf 'cycle=%s\nvalid_time=%s\ncase=%s\nensemble_size=%s\narchive_layout=%s\ndart_root=%s\nslurm_job_id=reconstructed\ncompleted_at=unknown\nreconstructed_at=%s\nreconstruction_reason=explicit_forward_xml_recovery\nprevious_xml_time=%s\n' \
            "${my_e3sm_completed_cycles}" "${expected}" "${my_casename}" "${my_ensnum}" "per_member" "${my_dart_root}" "$(date '+%Y-%m-%d %H:%M:%S')" "${common_actual}" > "${marker_tmp}" || fail "could not write reconstructed marker"
        else
          printf 'cycle=%s\nvalid_time=%s\ncase=%s\nensemble_size=%s\narchive_layout=%s\ndart_root=%s\nslurm_job_id=reconstructed\ncompleted_at=unknown\nreconstructed_at=%s\nreconstruction_reason=verified_counter_rollback_from_%s\n' \
            "${my_e3sm_completed_cycles}" "${expected}" "${my_casename}" "${my_ensnum}" "per_member" "${my_dart_root}" "$(date '+%Y-%m-%d %H:%M:%S')" "${common_actual}" > "${marker_tmp}" || fail "could not write reconstructed marker"
        fi
        mv -f "${marker_tmp}" "${marker}" || fail "could not commit reconstructed marker"
        echo "Reconstructed verified recovery marker: ${marker}"
      fi
      if [[ "${forward_recovery}" == "TRUE" ]]; then
        echo "Validated explicit forward XML recovery from ${common_actual} to configured cycle ${my_e3sm_completed_cycles} at ${expected}"
      else
        echo "Validated intentional rollback from XML state ${common_actual} to configured cycle ${my_e3sm_completed_cycles} at ${expected}"
      fi
    fi

    validate_marker "${marker}" "${expected}" "${my_casename}" "${my_ensnum}"
    marker_time=$(marker_field "${marker}" valid_time) || fail "completion record lacks valid_time: ${marker}"
    marker_cycle=$(marker_field "${marker}" cycle) || fail "completion record lacks cycle: ${marker}"
    [[ "${marker_time}" == "${expected}" && "${marker_cycle}" == "${my_e3sm_completed_cycles}" ]] || fail "cycle completion record does not match configured state: ${marker}"
  fi

  if [[ "${rollback}" == "TRUE" ]]; then
    echo "Rollback preflight passed; the cycle worker will reset all case XML clocks to ${expected}"
  elif [[ "${forward_recovery}" == "TRUE" ]]; then
    echo "Forward-recovery preflight passed; the cycle worker will advance all case XML clocks to ${expected}"
  else
    echo "Validated workflow state ${expected} across counter, completion record, ${my_ensnum} member archives, and case XML files"
  fi
}

reached_configured_end() {
  local current_time current_stamp end_stamp
  current_time=$(expected_state_time) || fail "could not derive current state from the configured timeline"
  current_stamp="${current_time//-/}"
  end_stamp="${my_e3sm_end_date//-/}${my_e3sm_end_tod}"
  (( 10#${current_stamp} >= 10#${end_stamp} ))
}

validate_cycle_end() {
  local prefix date_var tod_var date_value tod_value component_stamp
  local e3sm_stamp="${my_e3sm_end_date//-/}${my_e3sm_end_tod}"
  for prefix in my_e3sm my_eam_dart my_elm_dart; do
    date_var="${prefix}_end_date"
    tod_var="${prefix}_end_tod"
    date_value="${!date_var:-}"
    tod_value="${!tod_var:-}"
    [[ "${date_value}" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]] || fail "${date_var} must use YYYY-MM-DD format, got: ${date_value:-unset}"
    [[ "${tod_value}" =~ ^[0-9]{5}$ ]] || fail "${tod_var} must be five-digit seconds since midnight, got: ${tod_value:-unset}"
    (( 10#${tod_value} < 86400 )) || fail "${tod_var} is out of range: ${tod_value}"
    date -d "${date_value}" +%F >/dev/null 2>&1 || fail "invalid date in ${date_var}: ${date_value}"
    if [[ "${prefix}" != "my_e3sm" ]]; then
      component_stamp="${date_value//-/}${tod_value}"
      (( 10#${component_stamp} <= 10#${e3sm_stamp} )) || fail "${prefix} end time must not exceed the E3SM end time"
    fi
  done
}

write_continuation_record() {
  local record_file="$1" cycle="$2" job_id="$3" job_name="$4"
  local tmp_file="${record_file}.tmp.${SLURM_JOB_ID:-$$}"
  printf 'cycle=%s\njob_id=%s\njob_name=%s\nsubmitted_by=%s\nsubmitted_at=%s\n' \
    "${cycle}" "${job_id}" "${job_name}" "${SLURM_JOB_ID:-none}" "$(date '+%Y-%m-%d %H:%M:%S')" > "${tmp_file}"
  mv -f "${tmp_file}" "${record_file}"
}

ensure_continuation_job() {
  local cycle="${my_e3sm_completed_cycles}" status_dir record_file intent_file job_name
  local existing_job_id next_job_id scheduler_user workflow_id
  local -a dependency=()
  status_dir="${my_status_dir}"
  record_file="${status_dir}/continuation.cycle${cycle}"
  intent_file="${record_file}.intent"
  workflow_id=$(printf '%s' "${my_wkdir}" | cksum | awk '{print $1}')
  job_name="e3sm-${workflow_id}-c${cycle}"
  scheduler_user="${SLURM_JOB_USER:-${USER:-}}"
  [[ -n "${scheduler_user}" ]] || fail "cannot determine Slurm user for continuation lookup"
  mkdir -p "${status_dir}"

  printf 'cycle=%s\njob_name=%s\nrequested_by=%s\n' \
    "${cycle}" "${job_name}" "${SLURM_JOB_ID:-none}" > "${intent_file}.tmp.${SLURM_JOB_ID:-$$}"
  mv -f "${intent_file}.tmp.${SLURM_JOB_ID:-$$}" "${intent_file}"

  # Recover a successful sbatch if interruption occurred before the local
  # record was written. The workspace flock serializes this check/submit path.
  existing_job_id=$(squeue -h -u "${scheduler_user}" -n "${job_name}" -o '%A' |
    awk -v current="${SLURM_JOB_ID:-none}" 'NF && $1 != current {print; exit}')
  if [[ -n "${existing_job_id}" ]]; then
    [[ "${existing_job_id}" =~ ^[0-9]+$ ]] || fail "invalid existing continuation job ID: ${existing_job_id}"
    write_continuation_record "${record_file}" "${cycle}" "${existing_job_id}" "${job_name}"
    echo "Recovered existing continuation for cycle ${cycle}: job ${existing_job_id}"
    return 0
  fi

  [[ -n "${SLURM_JOB_ID:-}" ]] && dependency=(--dependency="afterok:${SLURM_JOB_ID}")
  next_job_id=$(sbatch --parsable --job-name="${job_name}" "${dependency[@]}" "${SCRIPT_PATH}")
  next_job_id=${next_job_id%%;*}
  [[ "${next_job_id}" =~ ^[0-9]+$ ]] || fail "could not parse continuation job ID: ${next_job_id}"
  write_continuation_record "${record_file}" "${cycle}" "${next_job_id}" "${job_name}"
  echo "Submitted continuation for cycle ${cycle} as job ${next_job_id}"
}

CYCLES_PER_JOB="${CYCLES_PER_JOB:-${my_cycles_per_job:-1}}"
MIN_CYCLE_TIME_SEC="${MIN_CYCLE_TIME_SEC:-${my_min_cycle_time_sec:-12600}}"
SHUTDOWN_MARGIN_SEC="${SHUTDOWN_MARGIN_SEC:-${my_cycle_shutdown_margin_sec:-900}}"
validate_positive_int "CYCLES_PER_JOB" "${CYCLES_PER_JOB}"
CYCLE_WORKER="${my_workflow_lib}/cycle/e3sm_dart_single_cycle.sh"
validate_on_off "my_eam_dart_da" "${my_eam_dart_da}"
validate_on_off "my_elm_dart_da" "${my_elm_dart_da}"
validate_true_false "my_allow_forward_xml_recovery" "${my_allow_forward_xml_recovery:-FALSE}"
validate_positive_int "my_e3sm_cycle_hours" "${my_e3sm_cycle_hours}"
[[ "${my_e3sm_completed_cycles:-}" =~ ^[0-9]+$ ]] || fail "my_e3sm_completed_cycles must be a non-negative integer"
validate_positive_int "my_eam_dart_cycle_hours" "${my_eam_dart_cycle_hours}"
validate_positive_int "my_elm_dart_cycle_hours" "${my_elm_dart_cycle_hours}"
(( my_eam_dart_cycle_hours % my_e3sm_cycle_hours == 0 )) || fail "my_eam_dart_cycle_hours must be an integer multiple of my_e3sm_cycle_hours"
(( my_elm_dart_cycle_hours % my_e3sm_cycle_hours == 0 )) || fail "my_elm_dart_cycle_hours must be an integer multiple of my_e3sm_cycle_hours"
validate_cycle_end
ELM_EXECUTION_MODE="direct"
if [[ "${strongly_coupled_on,,}" == "on" && "${lnd_da_use_sequential_prior_post,,}" == ".true." ]]; then
  ELM_EXECUTION_MODE="sequential"
fi
if [[ "${my_elm_dart_da,,}" == "on" && "${ELM_EXECUTION_MODE}" == "sequential" ]]; then
  [[ "${my_eam_dart_da,,}" == "on" ]] || fail "sequential-prior ELM requires my_eam_dart_da=on"
  [[ "${atm_da_output_sequential_prior_post,,}" == ".true." ]] || fail "sequential-prior ELM requires atm_da_output_sequential_prior_post=.true."
  [[ "${lnd_da_output_sequential_prior_post,,}" == ".false." ]] || fail "ELM cannot use and output sequential priors in the same filter pass"
  (( my_elm_dart_cycle_hours % my_eam_dart_cycle_hours == 0 )) || fail "ELM sequential-prior cadence must be an integer multiple of the EAM DA cadence"
  eam_end_stamp="${my_eam_dart_end_date//-/}${my_eam_dart_end_tod}"
  elm_end_stamp="${my_elm_dart_end_date//-/}${my_elm_dart_end_tod}"
  (( 10#${elm_end_stamp} <= 10#${eam_end_stamp} )) || fail "ELM sequential-prior end time must not exceed the EAM DA end time"
elif [[ "${my_elm_dart_da,,}" == "on" && "${strongly_coupled_on,,}" == "on" ]]; then
  [[ "${lnd_da_strongly_coupled,,}" == ".false." && "${lnd_da_state_model}" == "Land" && "${lnd_da_obs_model}" == "Land" ]] \
    || fail "direct-observation ELM requires lnd_da_strongly_coupled=.false., state_model=Land, and obs_model=Land"
fi
validate_positive_int "MIN_CYCLE_TIME_SEC" "${MIN_CYCLE_TIME_SEC}"
validate_positive_int "SHUTDOWN_MARGIN_SEC" "${SHUTDOWN_MARGIN_SEC}"
validate_positive_int "my_job_nnodes" "${my_job_nnodes}"
if [[ "${my_eam_dart_da,,}" == "on" && "${my_elm_dart_da,,}" == "on" && "${ELM_EXECUTION_MODE}" == "direct" ]] && (( my_job_nnodes % 2 != 0 )); then
  fail "my_job_nnodes must be even when concurrent EAM and ELM DART are both enabled"
fi
if [[ -n "${SLURM_JOB_NUM_NODES:-}" && "${SLURM_JOB_NUM_NODES}" != "${my_job_nnodes}" ]]; then
  fail "Slurm allocation has ${SLURM_JOB_NUM_NODES} nodes; configuration requires ${my_job_nnodes}"
fi
[[ -x "${CYCLE_WORKER}" ]] || fail "cycle worker is not executable: ${CYCLE_WORKER}"

SCRIPT_NAME=$(basename "${SCRIPT_PATH}")
RUN_START_EPOCH=$(date +%s)
on_exit_summary() {
  local exit_code=$?
  local elapsed=$(( $(date +%s) - RUN_START_EPOCH ))
  if (( exit_code == 0 )); then
    echo "INFO: ${SCRIPT_NAME} completed successfully in ${elapsed}s"
  else
    echo "ERROR: ${SCRIPT_NAME} failed with exit ${exit_code} after ${elapsed}s"
  fi
}
trap on_exit_summary EXIT

on_signal() {
  echo "ERROR: ${SCRIPT_NAME} received a termination signal"
  jobs -pr | while IFS= read -r child_pid; do kill -TERM "${child_pid}" 2>/dev/null || true; done
  wait 2>/dev/null || true
  exit 143
}
trap on_signal INT TERM
cd "${my_wkdir}"
echo "== Start of ${SCRIPT_NAME}: up to ${CYCLES_PER_JOB} cycle(s) in this allocation =="
date

validate_runtime_state
if reached_configured_end; then
  echo "Configured DA end time already reached; no cycle will be run"
  exit 0
fi

continuation_needed=FALSE
for ((cycle_index=1; cycle_index<=CYCLES_PER_JOB; cycle_index++)); do
  required_sec=$((MIN_CYCLE_TIME_SEC + SHUTDOWN_MARGIN_SEC))
  if remaining_sec=$(remaining_allocation_seconds); then
    echo "Remaining allocation time: ${remaining_sec}s; required for a cycle: ${required_sec}s"
    if (( remaining_sec < required_sec )); then
      echo "Insufficient time for a cycle; stopping cleanly"
      continuation_needed=TRUE
      break
    fi
  elif [[ -n "${SLURM_JOB_ID:-}" ]]; then
    echo "WARNING: unable to determine remaining Slurm time; stopping conservatively"
    continuation_needed=TRUE
    break
  else
    echo "INFO: no Slurm job detected; allowing a directly invoked cycle"
  fi

  echo "== Running cycle ${cycle_index}/${CYCLES_PER_JOB} in Slurm job ${SLURM_JOB_ID:-unknown} =="
  CYCLE_DRIVER_ACTIVE=TRUE DEFER_CYCLE_RESUBMIT=TRUE "${CYCLE_WORKER}"

  # Refresh the cycle counter and end-time settings atomically updated by the worker.
  source "${CONFIG_FILE}"
  validate_cycle_end
  validate_runtime_state
  if reached_configured_end; then
    echo "Configured DA end time reached; no continuation will be submitted"
    exit 0
  fi
  continuation_needed=TRUE
done

if [[ "${continuation_needed}" == "TRUE" ]]; then
  ensure_continuation_job
fi

echo "== End of ${SCRIPT_NAME} =="
date
