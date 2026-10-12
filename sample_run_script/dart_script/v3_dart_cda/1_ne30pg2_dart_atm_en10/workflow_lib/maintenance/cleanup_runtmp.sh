#!/bin/bash

set -Eeuo pipefail

# Safe defaults. Run with --apply only after reviewing the dry-run output.
DRY_RUN="TRUE"
LOG_RETENTION_DAYS=30
CLEAN_STALE_PROGRESS="FALSE"
STALE_PROGRESS_DAYS=14
CLEAN_FAILED_HANDOFFS="FALSE"
FAILED_HANDOFF_RETENTION_DAYS=30

fail() { echo "ERROR: $*" >&2; exit 1; }
normalize_bool() {
  local name="$1" value="${2^^}"
  [[ "${value}" == "TRUE" || "${value}" == "FALSE" ]] || fail "${name} must be TRUE or FALSE"
  printf '%s' "${value}"
}

case "${1:-}" in
  ""|--dry-run) DRY_RUN="TRUE" ;;
  --apply) DRY_RUN="FALSE" ;;
  *) fail "usage: $0 [--dry-run|--apply]" ;;
esac

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || fail "cannot resolve script path"
WORKFLOW_ROOT=$(readlink -f "$(dirname "${SCRIPT_PATH}")/../..") || fail "cannot resolve workflow root"
CONFIG_FILE="${WORKFLOW_ROOT}/create_and_setup_case.sh"
[[ -r "${CONFIG_FILE}" ]] || fail "missing workflow configuration: ${CONFIG_FILE}"
source "${CONFIG_FILE}"
mkdir -p "${my_lock_dir}"

for cmd in awk du find flock readlink rm; do
  command -v "${cmd}" >/dev/null 2>&1 || fail "required command not found: ${cmd}"
done
for setting in CLEAN_STALE_PROGRESS CLEAN_FAILED_HANDOFFS; do
  printf -v "${setting}" '%s' "$(normalize_bool "${setting}" "${!setting}")"
done
for setting in LOG_RETENTION_DAYS STALE_PROGRESS_DAYS FAILED_HANDOFF_RETENTION_DAYS; do
  [[ "${!setting}" =~ ^[0-9]+$ ]] || fail "${setting} must be a non-negative integer"
done

for path in "${my_log_dir}" "${my_status_dir}" "${my_lock_dir}" "${my_handoff_dir}"; do
  case "${path}" in
    "${my_runtime_dir}"/*) ;;
    *) fail "runtime path escapes my_runtime_dir: ${path}" ;;
  esac
done

exec 8>"${my_lock_dir}/runtime_cleanup.lock"
flock -n 8 || fail "another runtime cleanup is active"

# Refuse cleanup when any workflow lock is currently held. Lock files themselves
# are intentionally retained: they consume negligible space and unlinking a lock
# can create a race with a new process.
for lock_file in "${my_lock_dir}"/*.lock; do
  [[ -e "${lock_file}" && "${lock_file}" != "${my_lock_dir}/runtime_cleanup.lock" ]] || continue
  exec {lock_fd}<>"${lock_file}"
  flock -n "${lock_fd}" || fail "workflow appears active; held lock: ${lock_file}"
  flock -u "${lock_fd}"
  eval "exec ${lock_fd}>&-"
done

remove_file() {
  local file="$1"
  [[ -f "${file}" && ! -L "${file}" ]] || fail "refusing unexpected cleanup target: ${file}"
  if [[ "${DRY_RUN}" == "TRUE" ]]; then
    echo "WOULD REMOVE file: ${file}"
  else
    rm -f -- "${file}"
    echo "REMOVED file: ${file}"
  fi
}

remove_handoff_dir() {
  local dir="$1" base
  base=${dir##*/}
  [[ "${dir}" == "${my_handoff_dir}/"* && "${base}" =~ ^cycle_[0-9]+\.job_[0-9]+$ && -d "${dir}" && ! -L "${dir}" ]] ||
    fail "refusing unexpected handoff target: ${dir}"
  if [[ "${DRY_RUN}" == "TRUE" ]]; then
    echo "WOULD REMOVE failed handoff: ${dir}"
  else
    rm -rf -- "${dir}"
    echo "REMOVED failed handoff: ${dir}"
  fi
}

before_size=$(du -sh "${my_runtime_dir}" | awk '{print $1}')
echo "Runtime cleanup mode: $([[ "${DRY_RUN}" == "TRUE" ]] && echo DRY-RUN || echo APPLY)"
echo "Runtime size before cleanup: ${before_size}"

echo "Scanning logs older than ${LOG_RETENTION_DAYS} day(s)"
while IFS= read -r -d '' file; do
  remove_file "${file}"
done < <(find "${my_log_dir}" -type f ! -name 'README*' -mtime "+${LOG_RETENTION_DAYS}" -print0)

if [[ "${CLEAN_STALE_PROGRESS}" == "TRUE" ]]; then
  echo "Scanning stale in-progress records older than ${STALE_PROGRESS_DAYS} day(s)"
  while IFS= read -r -d '' file; do
    remove_file "${file}"
  done < <(find "${my_status_dir}" -maxdepth 1 -type f -name '*in_progress*' -mtime "+${STALE_PROGRESS_DAYS}" -print0)
else
  echo "Preserving all in-progress records (CLEAN_STALE_PROGRESS=FALSE)"
fi

if [[ "${CLEAN_FAILED_HANDOFFS}" == "TRUE" ]]; then
  echo "Scanning failed handoffs older than ${FAILED_HANDOFF_RETENTION_DAYS} day(s)"
  while IFS= read -r -d '' dir; do
    remove_handoff_dir "${dir}"
  done < <(find "${my_handoff_dir}" -mindepth 1 -maxdepth 1 -type d -name 'cycle_*.job_*' -mtime "+${FAILED_HANDOFF_RETENTION_DAYS}" -print0)
else
  echo "Preserving failed handoffs (CLEAN_FAILED_HANDOFFS=FALSE)"
fi

after_size=$(du -sh "${my_runtime_dir}" | awk '{print $1}')
echo "Runtime size after cleanup: ${after_size}"
echo "Completion records, generated run scripts, lock files, and recent logs were preserved."
