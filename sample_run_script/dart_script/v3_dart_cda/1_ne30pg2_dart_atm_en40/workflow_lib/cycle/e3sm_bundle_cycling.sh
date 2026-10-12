#!/bin/bash -el

# Sourced by e3sm_dart_single_cycle.sh after DART completes.
echo "== Start of e3sm_bundle_cycling.sh =="
date
echo "============================================"

cycling_fail() {
  echo "ERROR: $*"
  return 1
}

valid_eam_initial_file() {
  local file="$1"
  [[ -s "${file}" ]] || return 1
  ncdump -h "${file}" >/dev/null 2>&1 || return 1
  if command -v ncks >/dev/null 2>&1; then
    ncks -m -v PS,U,V,T,Q "${file}" >/dev/null 2>&1 || return 1
  fi
}

update_cycle_counter() {
  local config_file="${my_wkdir}/create_and_setup_case.sh"
  local tmp_file="${config_file}.tmp.${SLURM_JOB_ID:-$$}"
  local match_count

  match_count=$(grep -c '^export my_e3sm_completed_cycles=' "${config_file}") || return 1
  if [[ "${match_count}" -ne 1 ]]; then
    echo "ERROR: expected exactly one my_e3sm_completed_cycles setting in ${config_file}, found ${match_count}"
    return 1
  fi
  awk -v cycle="${DATA_ASSIMILATION_CYCLES}" '
    /^export my_e3sm_completed_cycles=/ { print "export my_e3sm_completed_cycles=" cycle; next }
    /^export my_allow_forward_xml_recovery=/ { print "export my_allow_forward_xml_recovery=\"FALSE\""; next }
    { print }
  ' "${config_file}" > "${tmp_file}" || { rm -f "${tmp_file}"; return 1; }
  chmod --reference="${config_file}" "${tmp_file}" || { rm -f "${tmp_file}"; return 1; }
  mv -f "${tmp_file}" "${config_file}" || return 1
}

for cmd in awk basename chmod cp find grep mkdir mv ncdump rm; do
  command -v "${cmd}" >/dev/null 2>&1 || cycling_fail "required command not found: ${cmd}"
done

NEXT_DATE=$(printf "%04d-%02d-%02d" "${DART_YEAR}" "${DART_MONTH}" "${DART_DAY}")
NEXT_TOD=$(printf "%05d" "${DART_SECONDS}")
NEXT_STAMP="${NEXT_DATE//-/}${NEXT_TOD}"
END_STAMP="${my_e3sm_end_date//-/}${my_e3sm_end_tod}"
MAX_PARALLEL_HANDOFF="${MAX_PARALLEL_HANDOFF:-${my_max_parallel_handoff:-8}}"
[[ "${MAX_PARALLEL_HANDOFF}" =~ ^[1-9][0-9]*$ ]] || cycling_fail "MAX_PARALLEL_HANDOFF must be a positive integer"

DEFER_CYCLE_RESUBMIT="${DEFER_CYCLE_RESUBMIT:-TRUE}"
DEFER_CYCLE_RESUBMIT=$(printf "%s" "${DEFER_CYCLE_RESUBMIT}" | tr '[:lower:]' '[:upper:]')
[[ "${DEFER_CYCLE_RESUBMIT}" == "TRUE" || "${DEFER_CYCLE_RESUBMIT}" == "FALSE" ]] || cycling_fail "DEFER_CYCLE_RESUBMIT must be TRUE or FALSE"

handoff_id="cycle${DATA_ASSIMILATION_CYCLES}.job${SLURM_JOB_ID:-$$}"
handoff_dir="${my_handoff_dir}/cycle_${DATA_ASSIMILATION_CYCLES}.job_${SLURM_JOB_ID:-$$}"
mkdir -p "${handoff_dir}"
config_snapshot="${handoff_dir}/create_and_setup_case.sh.before"
cp -p "${my_wkdir}/create_and_setup_case.sh" "${config_snapshot}" || cycling_fail "could not snapshot cycle configuration"
prune_older_completed_handoffs() {
  local current_marker="$1" current_cycle current_time dir base old_cycle old_job marker marker_cycle marker_job
  [[ -d "${handoff_dir}" && -s "${current_marker}" ]] || return 1
  current_cycle=$(awk -F= '$1 == "cycle" {print $2; exit}' "${current_marker}") || return 1
  current_time=$(awk -F= '$1 == "valid_time" {print $2; exit}' "${current_marker}") || return 1
  [[ "${current_cycle}" == "${DATA_ASSIMILATION_CYCLES}" && "${current_time}" == "${NEXT_DATE}-${NEXT_TOD}" ]] || return 1
  while IFS= read -r -d '' dir; do
    [[ "${dir}" != "${handoff_dir}" ]] || continue
    base=${dir##*/}
    [[ "${base}" =~ ^cycle_([0-9]+)\.job_([0-9]+)$ ]] || continue
    old_cycle=${BASH_REMATCH[1]}
    old_job=${BASH_REMATCH[2]}
    (( 10#${old_cycle} < 10#${current_cycle} )) || continue
    marker=""
    while IFS= read -r candidate; do
      marker_cycle=$(awk -F= '$1 == "cycle" {print $2; exit}' "${candidate}")
      marker_job=$(awk -F= '$1 == "slurm_job_id" {print $2; exit}' "${candidate}")
      if [[ "${marker_cycle}" == "${old_cycle}" && "${marker_job}" == "${old_job}" ]]; then marker=${candidate}; break; fi
    done < <(find "${my_status_dir}" -maxdepth 1 -type f -name 'cycle_complete.*' -print)
    if [[ -z "${marker}" ]]; then
      echo "WARNING: retaining unverified handoff snapshot: ${dir}" >&2
      continue
    fi
    rm -rf -- "${dir}" || return 1
    echo "Removed older completed handoff snapshot: ${dir}"
  done < <(find "${my_handoff_dir}" -mindepth 1 -maxdepth 1 -type d -name 'cycle_*.job_*' -print0)
}

HANDOFF_TRANSACTION_ACTIVE=FALSE
CYCLE_COUNTER_UPDATED=FALSE
committed_members=()

stage_handoff_member() (
  set -e
  local i="$1" ENSTR CASE_NAME CASE_DIR RUN_DIR member_archive eam_ic run_ic staged_ic
  ENSTR=$(printf "EN%02d" "${i}")
  CASE_NAME="${my_casename}.${ENSTR}"
  CASE_DIR="${CASE_ROOT/EN01/${ENSTR}}"
  RUN_DIR="${RUN_ROOT/EN01/${ENSTR}}"
  member_archive="${my_modeldir}/${ENSTR}/archive"
  eam_ic="${member_archive}/rest/${NEXT_DATE}-${NEXT_TOD}/${CASE_NAME}.eam.i.${NEXT_DATE}-${NEXT_TOD}.nc"
  run_ic="${RUN_DIR}/$(basename "${eam_ic}")"
  staged_ic="${run_ic}.handoff.${handoff_id}"

  [[ -d "${CASE_DIR}" ]]
  [[ -x "${CASE_DIR}/xmlchange" ]]
  [[ -s "${CASE_DIR}/env_run.xml" ]]
  [[ -d "${RUN_DIR}" ]]
  valid_eam_initial_file "${eam_ic}"
  cp -p "${CASE_DIR}/env_run.xml" "${handoff_dir}/${ENSTR}.env_run.xml" || return 1
  printf '%s\n' "${staged_ic}" > "${handoff_dir}/${ENSTR}.staged_ic" || return 1
  cp -p "${eam_ic}" "${staged_ic}" || return 1
  valid_eam_initial_file "${staged_ic}" || return 1
)

wait_for_handoff_batch() {
  local failed=0 idx pid member status
  for idx in "${!handoff_pids[@]}"; do
    pid="${handoff_pids[$idx]}"
    member="${handoff_names[$idx]}"
    if wait "${pid}"; then
      echo "Handoff staging completed for ${member}"
    else
      status=$?
      echo "ERROR: handoff staging failed for ${member} (status ${status})"
      failed=1
    fi
  done
  handoff_pids=()
  handoff_names=()
  (( failed == 0 ))
}

cleanup_staged_handoff() {
  local staged_file
  while IFS= read -r staged_file; do
    [[ -n "${staged_file}" ]] && rm -f "${staged_file}"
  done < <(awk 'NF' "${handoff_dir}"/*.staged_ic 2>/dev/null || true)
}

rollback_committed_members() {
  local ENSTR CASE_DIR restore_tmp
  echo "Rolling back case metadata for committed handoff members"
  for ENSTR in "${committed_members[@]}"; do
    CASE_DIR="${CASE_ROOT/EN01/${ENSTR}}"
    restore_tmp="${CASE_DIR}/env_run.xml.rollback.${handoff_id}"
    if cp -p "${handoff_dir}/${ENSTR}.env_run.xml" "${restore_tmp}"; then
      mv -f "${restore_tmp}" "${CASE_DIR}/env_run.xml" || echo "ERROR: rollback move failed for ${ENSTR}"
    else
      echo "ERROR: rollback copy failed for ${ENSTR}"
    fi
  done
}

restore_cycle_counter() {
  local config_file="${my_wkdir}/create_and_setup_case.sh"
  local restore_tmp="${config_file}.rollback.${handoff_id}"
  [[ "${CYCLE_COUNTER_UPDATED:-FALSE}" == "TRUE" ]] || return 0
  if cp -p "${config_snapshot}" "${restore_tmp}" && mv -f "${restore_tmp}" "${config_file}"; then
    CYCLE_COUNTER_UPDATED=FALSE
    echo "Restored cycle configuration after incomplete handoff"
    return 0
  fi
  echo "ERROR: failed to restore cycle configuration from ${config_snapshot}"
  return 1
}

handoff_abort_on_signal() {
  if [[ "${HANDOFF_TRANSACTION_ACTIVE:-FALSE}" == "TRUE" ]]; then
    if [[ -n "${status_file:-}" && -s "${status_file}" ]]; then
      echo "Handoff completion marker is durable; no signal rollback is needed"
    else
      rollback_committed_members || true
      restore_cycle_counter || true
    fi
  fi
  cleanup_staged_handoff || true
  HANDOFF_TRANSACTION_ACTIVE=FALSE
  return 0
}
commit_handoff_member() (
  local ENSTR="$1" CASE_DIR RUN_DIR staged_ic run_ic
  CASE_DIR="${CASE_ROOT/EN01/${ENSTR}}"
  RUN_DIR="${RUN_ROOT/EN01/${ENSTR}}"
  [[ -s "${handoff_dir}/${ENSTR}.staged_ic" ]] || return 1
  staged_ic=$(<"${handoff_dir}/${ENSTR}.staged_ic") || return 1
  [[ -s "${staged_ic}" ]] || return 1
  run_ic="${RUN_DIR}/$(basename "${staged_ic}" ".handoff.${handoff_id}")"
  cd "${CASE_DIR}" || return 1
  ./xmlchange RUN_STARTDATE="${NEXT_DATE}" || return 1
  ./xmlchange START_TOD="${NEXT_TOD}" || return 1
  ./xmlchange STOP_OPTION="nhours" || return 1
  ./xmlchange STOP_N="${DATA_ASSIMILATION_WINDOW}" || return 1
  ./xmlchange RUN_REFDATE="${NEXT_DATE}" || return 1
  ./xmlchange RUN_REFTOD="${NEXT_TOD}" || return 1
  mv -f "${staged_ic}" "${run_ic}" || return 1
)

echo "Validating and staging all next-cycle members (up to ${MAX_PARALLEL_HANDOFF} concurrently)"
handoff_pids=()
handoff_names=()
for i in $(seq 1 "${my_ensnum}"); do
  ENSTR=$(printf "EN%02d" "${i}")
  stage_log="${handoff_dir}/stage.${ENSTR}.log"
  stage_handoff_member "${i}" > "${stage_log}" 2>&1 &
  handoff_pids+=("$!")
  handoff_names+=("${ENSTR}")
  if (( ${#handoff_pids[@]} >= MAX_PARALLEL_HANDOFF )); then
    wait_for_handoff_batch || { cleanup_staged_handoff; cycling_fail "handoff staging failed; no case state was changed"; }
  fi
done
wait_for_handoff_batch || { cleanup_staged_handoff; cycling_fail "handoff staging failed; no case state was changed"; }
echo "All ${my_ensnum} next-cycle members passed validation and staging"
HANDOFF_TRANSACTION_ACTIVE=TRUE

committed_members=()
for i in $(seq 1 "${my_ensnum}"); do
  ENSTR=$(printf "EN%02d" "${i}")
  echo "Committing ${ENSTR} for ${NEXT_DATE}-${NEXT_TOD}"
  committed_members+=("${ENSTR}")
  if commit_handoff_member "${ENSTR}"; then
    :
  else
    rollback_committed_members
    cleanup_staged_handoff
    cycling_fail "handoff commit failed for ${ENSTR}; previously committed case metadata was rolled back"
  fi
done

cd "${my_wkdir}"
if update_cycle_counter; then
  CYCLE_COUNTER_UPDATED=TRUE
  echo "Updated my_e3sm_completed_cycles to ${DATA_ASSIMILATION_CYCLES}"
else
  rollback_committed_members
  cycling_fail "cycle counter update failed; committed case metadata was rolled back"
fi
status_dir="${my_status_dir}"
status_file="${status_dir}/cycle_complete.${NEXT_DATE}-${NEXT_TOD}"
status_tmp="${status_file}.tmp.${SLURM_JOB_ID:-$$}"
if ! mkdir -p "${status_dir}" ||
   ! printf 'cycle=%s\nvalid_time=%s-%s\ncase=%s\nensemble_size=%s\narchive_layout=%s\ndart_root=%s\nslurm_job_id=%s\ncompleted_at=%s\n' "${DATA_ASSIMILATION_CYCLES}" "${NEXT_DATE}" "${NEXT_TOD}" "${my_casename}" "${my_ensnum}" "per_member" "${my_dart_root}" "${SLURM_JOB_ID:-none}" "$(date '+%Y-%m-%d %H:%M:%S')" > "${status_tmp}" ||
   ! mv -f "${status_tmp}" "${status_file}"; then
  rm -f "${status_tmp}"
  rollback_committed_members
  restore_cycle_counter
  cycling_fail "could not durably record completed handoff; transaction was rolled back"
fi
HANDOFF_TRANSACTION_ACTIVE=FALSE
echo "Wrote cycle completion marker: ${status_file}"
if prune_older_completed_handoffs "${status_file}"; then
  echo "Retained most recent completed handoff snapshot: ${handoff_dir}"
else
  echo "WARNING: older handoff pruning was skipped because safe validation failed" >&2
fi

if [[ "${NEXT_STAMP}" -ge "${END_STAMP}" ]]; then
  echo "Reached configured DA end time: ${NEXT_DATE}-${NEXT_TOD}"
elif [[ "${DEFER_CYCLE_RESUBMIT}" == "TRUE" ]]; then
  echo "Next-cycle submission deferred to the multi-cycle driver"
else
  cycling_fail "cycle worker cannot submit continuations; run it through 4_run_dart_e3sm_cycleda.sh"
fi

echo "==================================="
