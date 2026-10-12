#!/bin/bash
#------------------------------------------------------------------------------
# Safe cleanup of runtmp/. Preview first, then apply exactly what was reported:
#   workflow_lib/maintenance/cleanup_runtmp.sh --dry-run
#   workflow_lib/maintenance/cleanup_runtmp.sh --apply
#
# Default: removes log files older than LOG_RETENTION_DAYS. Completion and
# in-progress records, lock files and run scripts are preserved. Refuses to
# run while any workflow lock is held (a driver is active).
#------------------------------------------------------------------------------
set -Eeo pipefail

# Safe defaults.
LOG_RETENTION_DAYS=30
CLEAN_STALE_PROGRESS="FALSE"     # also remove *_in_progress* records older than STALE_PROGRESS_DAYS
STALE_PROGRESS_DAYS=14

fail() { echo "ERROR: $*" >&2; exit 1; }
case "${1:-}" in
  ""|--dry-run) DRY_RUN="TRUE" ;;
  --apply)      DRY_RUN="FALSE" ;;
  *) fail "usage: $0 [--dry-run|--apply]" ;;
esac

SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}")
WORKFLOW_ROOT=$(readlink -f "$(dirname "${SCRIPT_PATH}")/../..")
[[ -r ${WORKFLOW_ROOT}/create_and_setup_case.sh ]] || fail "missing workflow configuration in ${WORKFLOW_ROOT}"
source "${WORKFLOW_ROOT}/create_and_setup_case.sh"
[[ ${CLEAN_STALE_PROGRESS} == TRUE || ${CLEAN_STALE_PROGRESS} == FALSE ]] || fail "CLEAN_STALE_PROGRESS must be TRUE or FALSE"
for setting in LOG_RETENTION_DAYS STALE_PROGRESS_DAYS; do
  [[ ${!setting} =~ ^[0-9]+$ ]] || fail "${setting} must be a non-negative integer"
done
for path in "${my_log_dir}" "${my_status_dir}" "${my_lock_dir}"; do
  case "${path}" in "${my_runtime_dir}"/*) ;; *) fail "runtime path escapes my_runtime_dir: ${path}" ;; esac
done
[[ -d ${my_runtime_dir} ]] || { echo "nothing to clean: ${my_runtime_dir} does not exist"; exit 0; }
mkdir -p "${my_lock_dir}"

exec 8>"${my_lock_dir}/runtime_cleanup.lock"
flock -n 8 || fail "another runtime cleanup is active"

# Refuse while any workflow lock is held.
for lock_file in "${my_lock_dir}"/*.lock; do
  [[ -e ${lock_file} && ${lock_file} != "${my_lock_dir}/runtime_cleanup.lock" ]] || continue
  exec 9<>"${lock_file}"
  flock -n 9 || fail "workflow appears active; held lock: ${lock_file}"
  flock -u 9
  exec 9>&-
done

remove_file() {
  local file="$1"
  [[ -f ${file} && ! -L ${file} ]] || fail "refusing unexpected cleanup target: ${file}"
  if [[ ${DRY_RUN} == TRUE ]]; then echo "WOULD REMOVE file: ${file}"; else rm -f -- "${file}"; echo "REMOVED file: ${file}"; fi
}

echo "Runtime cleanup mode: $([[ ${DRY_RUN} == TRUE ]] && echo DRY-RUN || echo APPLY)"
echo "Runtime size before cleanup: $(du -sh "${my_runtime_dir}" | awk '{print $1}')"
if [[ -d ${my_log_dir} ]]; then
  while IFS= read -r -d '' file; do remove_file "${file}"; done \
    < <(find "${my_log_dir}" -type f -mtime "+${LOG_RETENTION_DAYS}" -print0)
fi
if [[ ${CLEAN_STALE_PROGRESS} == TRUE && -d ${my_status_dir} ]]; then
  while IFS= read -r -d '' file; do remove_file "${file}"; done \
    < <(find "${my_status_dir}" -maxdepth 1 -type f -name '*in_progress*' -mtime "+${STALE_PROGRESS_DAYS}" -print0)
else
  echo "Preserving all in-progress records (CLEAN_STALE_PROGRESS=FALSE)"
fi
echo "Runtime size after cleanup: $(du -sh "${my_runtime_dir}" | awk '{print $1}')"
echo "Completion records, run scripts, lock files and recent logs were preserved."
