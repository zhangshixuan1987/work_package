#!/bin/bash
#SBATCH --account=esmd
#SBATCH --time=02:00:00
#SBATCH --partition=short
#SBATCH --job-name=e3sm_compress_diag
#SBATCH --nodes=1
#SBATCH --output=runtmp/logs/e3sm_compress_diag.%j
#SBATCH --exclusive
#SBATCH --no-kill
#SBATCH --requeue

set -Eeuo pipefail

if [[ "${STEP5_DRIVER_ACTIVE:-FALSE}" != "TRUE" ]]; then
  echo "ERROR: internal compression worker; run 5_run_dart_compress.sh" >&2
  exit 1
fi

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

normalize_bool() {
  local name="$1"
  local value="${2^^}"
  [[ "${value}" == "TRUE" || "${value}" == "FALSE" ]] || fail "${name} must be TRUE or FALSE, got: $2"
  printf '%s' "${value}"
}

if [[ -n "${WORKFLOW_ROOT:-}" ]]; then
  WORK_DIR=$(readlink -f "${WORKFLOW_ROOT}") || fail "cannot resolve WORKFLOW_ROOT"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  WORK_DIR=$(readlink -f "${SLURM_SUBMIT_DIR}") || fail "cannot resolve SLURM_SUBMIT_DIR"
else
  SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || fail "cannot resolve script path"
  WORK_DIR=$(dirname "${SCRIPT_PATH}")
fi
cd "${WORK_DIR}"

[[ -r "${WORK_DIR}/create_and_setup_case.sh" ]] || fail "missing workflow configuration"
source "${WORK_DIR}/create_and_setup_case.sh"
[[ -n "${my_dart_env_file:-}" ]] || fail "my_dart_env_file is unset"
[[ -r "${my_dart_env_file}" ]] || fail "configured DART environment is not readable: ${my_dart_env_file}"
echo "Using configured DART machine environment: ${my_dart_env_file}"
source "${my_dart_env_file}"
mkdir -p "${my_log_dir}" "${my_status_dir}" "${my_lock_dir}"

for cmd in awk chmod date df dirname flock mkdir mv nccopy ncdump readlink rm sort stat touch; do
  command -v "${cmd}" >/dev/null 2>&1 || fail "required command not found: ${cmd}"
done

[[ "${my_ensnum:-}" =~ ^[1-9][0-9]*$ ]] || fail "my_ensnum must be a positive integer"
[[ "${my_e3sm_completed_cycles:-}" =~ ^[0-9]+$ ]] || fail "my_e3sm_completed_cycles must be a non-negative integer"
[[ "${my_e3sm_cycle_hours:-}" =~ ^[1-9][0-9]*$ ]] || fail "my_e3sm_cycle_hours must be a positive integer"
[[ "${my_e3sm_start_date:-}" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]] || fail "my_e3sm_start_date must use YYYY-MM-DD format"
[[ "${my_e3sm_start_tod:-}" =~ ^[0-9]{5}$ ]] || fail "my_e3sm_start_tod must be five-digit seconds since midnight"
(( 10#${my_e3sm_start_tod} < 86400 )) || fail "my_e3sm_start_tod is out of range"
date -d "${my_e3sm_start_date}" +%F >/dev/null 2>&1 || fail "my_e3sm_start_date is not a valid date"

[[ $# -eq 2 ]] || fail "usage: $0 YYYY-MM-DD SSSSS"
ZIP_DATE="$1"
ZIP_TOD="$2"
[[ "${ZIP_DATE}" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]] || fail "invalid date: ${ZIP_DATE}"
[[ "${ZIP_TOD}" =~ ^[0-9]{5}$ ]] || fail "invalid seconds-of-day: ${ZIP_TOD}"
(( 10#${ZIP_TOD} < 86400 )) || fail "seconds-of-day out of range: ${ZIP_TOD}"
date -d "${ZIP_DATE}" +%F >/dev/null 2>&1 || fail "invalid date: ${ZIP_DATE}"

COMPRESS_HISTORY=$(normalize_bool COMPRESS_HISTORY "${COMPRESS_HISTORY:-TRUE}")
COMPRESS_RESTARTS=$(normalize_bool COMPRESS_RESTARTS "${COMPRESS_RESTARTS:-FALSE}")
HISTORY_MAX_PARALLEL="${COMPRESS_MAX_PARALLEL:-4}"
RESTART_MAX_PARALLEL="${COMPRESS_RESTART_MAX_PARALLEL:-1}"
SPACE_MARGIN_MB="${COMPRESS_SPACE_MARGIN_MB:-1024}"
[[ "${HISTORY_MAX_PARALLEL}" =~ ^[1-9][0-9]*$ ]] || fail "COMPRESS_MAX_PARALLEL must be positive"
[[ "${RESTART_MAX_PARALLEL}" =~ ^[1-9][0-9]*$ ]] || fail "COMPRESS_RESTART_MAX_PARALLEL must be positive"
[[ "${SPACE_MARGIN_MB}" =~ ^[0-9]+$ ]] || fail "COMPRESS_SPACE_MARGIN_MB must be non-negative"
[[ "${COMPRESS_HISTORY}" == "TRUE" || "${COMPRESS_RESTARTS}" == "TRUE" ]] || fail "both compression modes are disabled"

[[ "${COMPRESS_ENSTR:-}" =~ ^EN[0-9][0-9]$ ]] || fail "COMPRESS_ENSTR must identify one ensemble member"
ARCHIVE_DIR="${my_modeldir}/${COMPRESS_ENSTR}/archive"
STATUS_DIR="${my_status_dir}"
TIMESTAMP_DASH="${ZIP_DATE}-${ZIP_TOD}"
TIMESTAMP_UNDER="${ZIP_DATE}_${ZIP_TOD}"
JOB_TAG="${SLURM_JOB_ID:-$$}.${COMPRESS_ENSTR}"
LOCK_FILE="${my_lock_dir}/compress.${COMPRESS_ENSTR}.${TIMESTAMP_DASH}.lock"
COMPLETION_FILE="${STATUS_DIR}/compression_complete.${COMPRESS_ENSTR}.${TIMESTAMP_DASH}"
IN_PROGRESS_FILE="${STATUS_DIR}/compression_in_progress.${COMPRESS_ENSTR}.${TIMESTAMP_DASH}"
TEMP_FILES=()
CHILD_PIDS=()
mkdir -p "${STATUS_DIR}"
exec 8>"${LOCK_FILE}"
flock -n 8 || fail "another compression job is active for ${TIMESTAMP_DASH}"
if [[ -s "${COMPLETION_FILE}" ]]; then
  completed_history=$(awk -F= '$1 == "history" {print $2}' "${COMPLETION_FILE}")
  completed_restarts=$(awk -F= '$1 == "restarts" {print $2}' "${COMPLETION_FILE}")
  completed_case=$(awk -F= '$1 == "case" {print $2}' "${COMPLETION_FILE}")
  completed_size=$(awk -F= '$1 == "ensemble_size" {print $2}' "${COMPLETION_FILE}")
  if [[ "${completed_case}" == "${my_casename}" && "${completed_size}" == "${my_ensnum}" && ( "${COMPRESS_HISTORY}" == "FALSE" || "${completed_history}" == "TRUE" ) && ( "${COMPRESS_RESTARTS}" == "FALSE" || "${completed_restarts}" == "TRUE" ) ]]; then
    echo "Compression already completed for requested modes: ${COMPLETION_FILE}"
    exit 0
  fi
fi

cleanup_temps() {
  local file
  for file in "${TEMP_FILES[@]}"; do
    [[ -e "${file}" ]] && rm -f -- "${file}"
  done
}

on_signal() {
  local pid
  echo "ERROR: compression job received a termination signal" >&2
  for pid in "${CHILD_PIDS[@]}"; do
    kill -TERM "${pid}" 2>/dev/null || true
  done
  wait 2>/dev/null || true
  cleanup_temps
  exit 143
}
trap on_signal INT TERM
trap cleanup_temps EXIT

compute_current_cycle_stamp() {
  local base_seconds base_hour base_minute base_second base_epoch current_epoch current_date current_hour current_minute current_second current_tod
  base_seconds=$((10#${my_e3sm_start_tod}))
  base_hour=$((base_seconds / 3600))
  base_minute=$(((base_seconds % 3600) / 60))
  base_second=$((base_seconds % 60))
  base_epoch=$(date -d "${my_e3sm_start_date} $(printf '%02d:%02d:%02d' "${base_hour}" "${base_minute}" "${base_second}")" +%s) || return 1
  current_epoch=$((base_epoch + my_e3sm_completed_cycles * my_e3sm_cycle_hours * 3600))
  read -r current_date current_hour current_minute current_second < <(date -d "@${current_epoch}" '+%Y%m%d %H %M %S')
  current_tod=$((10#${current_hour} * 3600 + 10#${current_minute} * 60 + 10#${current_second}))
  printf '%s%05d\n' "${current_date}" "${current_tod}"
}

REQUESTED_STAMP="${ZIP_DATE//-/}${ZIP_TOD}"
CURRENT_CYCLE_STAMP=$(compute_current_cycle_stamp) || fail "cannot compute current workflow cycle"
if (( 10#${REQUESTED_STAMP} > 10#${CURRENT_CYCLE_STAMP} )); then
  fail "refusing to compress future cycle ${TIMESTAMP_DASH}; current workflow state is ${CURRENT_CYCLE_STAMP}"
fi

REST_DIR="${ARCHIVE_DIR}/rest/${TIMESTAMP_DASH}"
DART_MARKER="${my_dart_root}/transactions/${TIMESTAMP_DASH}/.dart_filter_in_progress"
INITIAL_STAMP="${my_e3sm_start_date}-${my_e3sm_start_tod}"
if [[ "${TIMESTAMP_DASH}" == "${INITIAL_STAMP}" ]]; then
  UPSTREAM_MARKER="${STATUS_DIR}/perturb_complete.${TIMESTAMP_DASH}"
else
  UPSTREAM_MARKER="${STATUS_DIR}/cycle_complete.${TIMESTAMP_DASH}"
fi
[[ -s "${UPSTREAM_MARKER}" ]] || fail "missing upstream completion record: ${UPSTREAM_MARKER}"
UPSTREAM_TIME=$(awk -F= '$1 == "valid_time" {sub(/^[^=]*=/, ""); print; exit}' "${UPSTREAM_MARKER}")
[[ "${UPSTREAM_TIME}" == "${TIMESTAMP_DASH}" ]] || fail "upstream completion valid_time mismatch: expected ${TIMESTAMP_DASH}, got ${UPSTREAM_TIME:-missing}"
if [[ "${COMPRESS_RESTARTS}" == "TRUE" && -e "${DART_MARKER}" ]]; then
  fail "refusing restart compression while DART is active or incomplete: ${DART_MARKER}"
fi
if [[ "${COMPRESS_RESTARTS}" == "TRUE" ]]; then
  exec 9>"${my_lock_dir}/e3sm_dart_cycle.lock"
  flock -n 9 || fail "refusing restart compression while the cycle driver is active"
fi

printf 'valid_time=%s\ncase=%s\nmember=%s\nensemble_size=%s\njob_id=%s\nstarted_at=%s\n' "${TIMESTAMP_DASH}" "${my_casename}" "${COMPRESS_ENSTR}" "${my_ensnum}" "${SLURM_JOB_ID:-none}" "$(date '+%Y-%m-%d %H:%M:%S')" > "${IN_PROGRESS_FILE}"

check_temp_space() {
  local source_file="$1"
  local target_dir source_bytes required_blocks available_blocks
  target_dir=$(dirname "${source_file}")
  source_bytes=$(stat -c %s "${source_file}") || return 1
  required_blocks=$(((source_bytes + 1023) / 1024 + SPACE_MARGIN_MB * 1024))
  available_blocks=$(df -Pk "${target_dir}" | awk 'NR == 2 {print $4}') || return 1
  [[ "${available_blocks}" =~ ^[0-9]+$ ]] || return 1
  if (( available_blocks < required_blocks )); then
    echo "ERROR: insufficient temporary space for ${source_file}: need ${required_blocks} KiB, have ${available_blocks} KiB" >&2
    return 1
  fi
}

compress_one() {
  local source_file="$1"
  local tmp_file="${source_file}.nccopy.${JOB_TAG}"
  [[ -s "${source_file}" ]] || { echo "ERROR: empty input: ${source_file}" >&2; return 1; }
  ncdump -h "${source_file}" >/dev/null 2>&1 || { echo "ERROR: invalid NetCDF input: ${source_file}" >&2; return 1; }
  check_temp_space "${source_file}" || return 1
  rm -f -- "${tmp_file}"
  nccopy -7 -d 5 "${source_file}" "${tmp_file}" || return 1
  [[ -s "${tmp_file}" ]] || return 1
  ncdump -h "${tmp_file}" >/dev/null 2>&1 || return 1
  chmod --reference="${source_file}" "${tmp_file}" || return 1
  touch --reference="${source_file}" "${tmp_file}" || return 1
  mv -f "${tmp_file}" "${source_file}" || return 1
  echo "Compressed ${source_file}"
}

check_set_space() {
  local max_parallel="$1"
  shift
  local files=("$@") source_file target_dir available_blocks required_blocks=0 count=0 source_bytes
  (( ${#files[@]} > 0 )) || return 0
  target_dir=$(dirname "${files[0]}")
  while read -r source_bytes; do
    required_blocks=$((required_blocks + (source_bytes + 1023) / 1024 + SPACE_MARGIN_MB * 1024))
    count=$((count + 1))
    (( count >= max_parallel )) && break
  done < <(for source_file in "${files[@]}"; do stat -c %s "${source_file}"; done | sort -nr)
  available_blocks=$(df -Pk "${target_dir}" | awk 'NR == 2 {print $4}')
  [[ "${available_blocks}" =~ ^[0-9]+$ ]] || return 1
  (( available_blocks >= required_blocks )) || { echo "ERROR: insufficient aggregate temporary space in ${target_dir}: need ${required_blocks} KiB, have ${available_blocks} KiB" >&2; return 1; }
}

compress_file_set() {
  local label="$1"
  local max_parallel="$2"
  shift 2
  local files=("$@")
  local paths=() pids=()
  local file pid index failed=0
  (( ${#files[@]} > 0 )) || { echo "No matching NetCDF files for ${label}"; return 0; }
  check_set_space "${max_parallel}" "${files[@]}" || return 1
  echo "Compressing ${#files[@]} file(s) for ${label}, up to ${max_parallel} concurrently"
  for file in "${files[@]}"; do
    TEMP_FILES+=("${file}.nccopy.${JOB_TAG}")
    compress_one "${file}" &
    pid=$!
    pids+=("${pid}")
    paths+=("${file}")
    CHILD_PIDS+=("${pid}")
    if (( ${#pids[@]} >= max_parallel )); then
      for index in "${!pids[@]}"; do
        if ! wait "${pids[$index]}"; then
          echo "ERROR: compression failed for ${paths[$index]}" >&2
          failed=1
        fi
      done
      pids=()
      paths=()
      CHILD_PIDS=()
    fi
  done
  for index in "${!pids[@]}"; do
    if ! wait "${pids[$index]}"; then
      echo "ERROR: compression failed for ${paths[$index]}" >&2
      failed=1
    fi
  done
  CHILD_PIDS=()
  (( failed == 0 ))
}

sum_file_bytes() {
  local total=0 file
  for file in "$@"; do
    total=$((total + $(stat -c %s "${file}")))
  done
  printf '%s' "${total}"
}

shopt -s nullglob
ALL_FILES=()
FAILED=0
if [[ "${my_runtype}" == "AMIP" ]]; then
  COMPONENTS=(atm cpl lnd rof)
else
  COMPONENTS=(atm cpl ice lnd ocn rof)
fi

HISTORY_FILE_COUNT=0
RESTART_FILE_COUNT=0

if [[ "${COMPRESS_HISTORY}" == "TRUE" ]]; then
  for component in "${COMPONENTS[@]}"; do
    history_dir="${ARCHIVE_DIR}/${component}/hist"
    [[ -d "${history_dir}" ]] || continue
    files=("${history_dir}"/*"${TIMESTAMP_DASH}"*.nc "${history_dir}"/*"${TIMESTAMP_UNDER}"*.nc)
    HISTORY_FILE_COUNT=$((HISTORY_FILE_COUNT + ${#files[@]}))
    ALL_FILES+=("${files[@]}")
    compress_file_set "${component}/hist at ${TIMESTAMP_DASH}" "${HISTORY_MAX_PARALLEL}" "${files[@]}" || FAILED=1
  done
  (( HISTORY_FILE_COUNT > 0 )) || fail "no history files found for requested timestamp ${TIMESTAMP_DASH}"
fi

if [[ "${COMPRESS_RESTARTS}" == "TRUE" ]]; then
  [[ -d "${REST_DIR}" ]] || fail "restart directory not found: ${REST_DIR}"
  restart_files=("${REST_DIR}"/*.nc)
  RESTART_FILE_COUNT=${#restart_files[@]}
  (( RESTART_FILE_COUNT > 0 )) || fail "no restart files found for requested timestamp ${TIMESTAMP_DASH}"
  ALL_FILES+=("${restart_files[@]}")
  compress_file_set "rest/${TIMESTAMP_DASH}" "${RESTART_MAX_PARALLEL}" "${restart_files[@]}" || FAILED=1
else
  echo "Restart compression disabled; set COMPRESS_RESTARTS=TRUE for an inactive cycle"
fi

(( FAILED == 0 )) || fail "one or more files failed compression"
FILE_COUNT=${#ALL_FILES[@]}
TOTAL_BYTES=$(sum_file_bytes "${ALL_FILES[@]}")
completion_tmp="${COMPLETION_FILE}.tmp.${JOB_TAG}"
printf 'valid_time=%s\ncase=%s\nmember=%s\nensemble_size=%s\njob_id=%s\ncompleted_at=%s\nfile_count=%s\nhistory_file_count=%s\nrestart_file_count=%s\ncompressed_bytes=%s\nhistory=%s\nrestarts=%s\n' "${TIMESTAMP_DASH}" "${my_casename}" "${COMPRESS_ENSTR}" "${my_ensnum}" "${SLURM_JOB_ID:-none}" "$(date '+%Y-%m-%d %H:%M:%S')" "${FILE_COUNT}" "${HISTORY_FILE_COUNT}" "${RESTART_FILE_COUNT}" "${TOTAL_BYTES}" "${COMPRESS_HISTORY}" "${COMPRESS_RESTARTS}" > "${completion_tmp}"
mv -f "${completion_tmp}" "${COMPLETION_FILE}"
rm -f -- "${IN_PROGRESS_FILE}"
echo "Compression completed for ${FILE_COUNT} file(s); record: ${COMPLETION_FILE}"
