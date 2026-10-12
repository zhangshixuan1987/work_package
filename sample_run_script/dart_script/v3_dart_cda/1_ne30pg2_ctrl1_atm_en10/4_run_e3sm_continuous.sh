#!/bin/bash
#------------------------------------------------------------------------------
# Continuous free-run driver (replaces DA cycling for this case).
# Run after Steps 1-3. Every member starts once (hybrid) from its perturbed
# eam.i file at my_e3sm_start_date and runs continuously to my_e3sm_end_date.
# The run is split into my_free_segment_days segments that continue with
# CONTINUE_RUN=TRUE; restarts are written only at segment ends so an
# interrupted job resumes from the last segment. When every member reaches the
# end time, intermediate restart sets are deleted and a completion record is
# written. No eam.i (inithist) files are produced.
#------------------------------------------------------------------------------
#SBATCH --account=e3sm
#SBATCH --time=02:00:00
#SBATCH --partition=short
#SBATCH --job-name=e3sm_ctrl_free
#SBATCH --nodes=20
#SBATCH --output=e3sm_ctrl_free.%j
#SBATCH --exclusive
#SBATCH --no-kill
#SBATCH --requeue

set -Eeuo pipefail

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  [[ -n "${SLURM_SUBMIT_DIR:-}" ]] || fail "SLURM_SUBMIT_DIR is unavailable"
  my_wkdir=$(readlink -f "${SLURM_SUBMIT_DIR}") || fail "could not resolve SLURM_SUBMIT_DIR"
else
  SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || fail "could not resolve script path"
  my_wkdir=$(dirname "${SCRIPT_PATH}")
fi
CONFIG_FILE="${my_wkdir}/create_and_setup_case.sh"
SCRIPT_PATH="${my_wkdir}/4_run_e3sm_continuous.sh"
[[ -r "${SCRIPT_PATH}" ]] || fail "missing persistent driver: ${SCRIPT_PATH}"
[[ -r "${CONFIG_FILE}" ]] || fail "missing configuration: ${CONFIG_FILE}"
source "${CONFIG_FILE}"
[[ -n "${my_dart_env_file:-}" ]] || fail "my_dart_env_file is unset"
[[ -r "${my_dart_env_file}" ]] || fail "configured machine environment is not readable: ${my_dart_env_file}"
source "${my_dart_env_file}"
mkdir -p "${my_log_dir}" "${my_status_dir}" "${my_lock_dir}"
LOG_DIR="${my_log_dir}"

for cmd in awk cksum cp date ex find flock ncdump readlink sbatch squeue srun stat tail; do
  command -v "${cmd}" >/dev/null 2>&1 || fail "required command not found: ${cmd}"
done
exec 9>"${my_lock_dir}/e3sm_free_run.lock"
flock -n 9 || fail "another free-run driver is already running for ${my_wkdir}"

validate_positive_int() {
  [[ "$2" =~ ^[1-9][0-9]*$ ]] || fail "$1 must be a positive integer, got: $2"
}

SEGMENT_DAYS="${my_free_segment_days:-8}"
KEEP_INTERMEDIATE="${my_free_keep_intermediate_restarts:-FALSE}"
KEEP_INTERMEDIATE="${KEEP_INTERMEDIATE^^}"
MIN_SEGMENT_TIME_SEC="${my_free_min_segment_time_sec:-21600}"
SHUTDOWN_MARGIN_SEC="${my_cycle_shutdown_margin_sec:-900}"
NODES_PER_MEMBER="${my_nodes_per_member:-2}"
MAX_PARALLEL_SETUP="${my_max_parallel_setup:-8}"
validate_positive_int "my_free_segment_days" "${SEGMENT_DAYS}"
validate_positive_int "my_free_min_segment_time_sec" "${MIN_SEGMENT_TIME_SEC}"
validate_positive_int "my_cycle_shutdown_margin_sec" "${SHUTDOWN_MARGIN_SEC}"
validate_positive_int "my_nodes_per_member" "${NODES_PER_MEMBER}"
validate_positive_int "my_max_parallel_setup" "${MAX_PARALLEL_SETUP}"
validate_positive_int "my_ensnum" "${my_ensnum}"
validate_positive_int "my_job_nnodes" "${my_job_nnodes}"
[[ "${KEEP_INTERMEDIATE}" == "TRUE" || "${KEEP_INTERMEDIATE}" == "FALSE" ]] || fail "my_free_keep_intermediate_restarts must be TRUE or FALSE"
(( my_job_nnodes % NODES_PER_MEMBER == 0 )) || fail "my_job_nnodes (${my_job_nnodes}) must be divisible by ${NODES_PER_MEMBER}"
BATCH_SIZE=$((my_job_nnodes / NODES_PER_MEMBER))
if [[ -n "${SLURM_JOB_NUM_NODES:-}" && "${SLURM_JOB_NUM_NODES}" != "${my_job_nnodes}" ]]; then
  fail "Slurm allocation has ${SLURM_JOB_NUM_NODES} nodes; configuration requires ${my_job_nnodes}"
fi
case "${my_runtype}" in
  Full-CPL|AMIP) ;;
  *) fail "my_runtype must be either Full-CPL or AMIP, got: ${my_runtype:-unset}" ;;
esac

#------------------------------------------------------------------------------
# Time helpers. Stamps use the E3SM form YYYY-MM-DD-SSSSS.
#------------------------------------------------------------------------------
stamp_to_epoch() {
  local stamp="$1" ymd tod
  [[ "${stamp}" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}-[0-9]{5}$ ]] || return 1
  ymd="${stamp%-*}"; tod=$((10#${stamp##*-}))
  (( tod < 86400 )) || return 1
  echo $(( $(date -u -d "${ymd}" +%s) + tod ))
}
epoch_to_stamp() {
  local ymd h m s
  read -r ymd h m s < <(date -u -d "@$1" '+%Y-%m-%d %H %M %S')
  printf '%s-%05d\n' "${ymd}" "$((10#${h} * 3600 + 10#${m} * 60 + 10#${s}))"
}

START_STAMP="${my_e3sm_start_date}-${my_e3sm_start_tod}"
END_STAMP="${my_e3sm_end_date}-${my_e3sm_end_tod}"
START_EPOCH=$(stamp_to_epoch "${START_STAMP}") || fail "invalid start time: ${START_STAMP}"
END_EPOCH=$(stamp_to_epoch "${END_STAMP}") || fail "invalid end time: ${END_STAMP}"
(( END_EPOCH > START_EPOCH )) || fail "end time ${END_STAMP} must be after start time ${START_STAMP}"
(( START_EPOCH % 86400 == 0 && END_EPOCH % 86400 == 0 )) || fail "start and end times must be at 00Z (whole days)"
[[ "${START_STAMP}" == "${my_refdate}-${my_reftod}" ]] || fail "start time ${START_STAMP} differs from the Step 3 reference time ${my_refdate}-${my_reftod}"

#------------------------------------------------------------------------------
# Upstream (Step 3) completion record.
#------------------------------------------------------------------------------
marker_field() { awk -F= -v key="$2" '$1 == key {sub(/^[^=]*=/, ""); print; found=1; exit} END {if (!found) exit 1}' "$1"; }
perturb_marker="${my_status_dir}/perturb_complete.${my_refdate}-${my_reftod}"
[[ -s "${perturb_marker}" ]] || fail "missing Step 3 completion record: ${perturb_marker}"
[[ "$(marker_field "${perturb_marker}" case)" == "${my_casename}" ]] || fail "Step 3 record is for a different case: ${perturb_marker}"
[[ "$(marker_field "${perturb_marker}" ensemble_size)" == "${my_ensnum}" ]] || fail "Step 3 record has a different ensemble size: ${perturb_marker}"

COMPLETE_MARKER="${my_status_dir}/free_run_complete.${START_STAMP}-${END_STAMP}"

#------------------------------------------------------------------------------
# Member helpers
#------------------------------------------------------------------------------
member_name() { printf 'EN%02d' "$1"; }
member_case_dir() { echo "${my_modeldir}/$1/case_scripts"; }
member_run_dir() { echo "${my_modeldir}/$1/run"; }
member_archive() { echo "${my_modeldir}/$1/archive"; }

# Latest archived time with a real (non-empty) EAM restart; the start time when
# the member has not completed a segment yet. Step 2 leaves a zero-byte eam.r
# placeholder at the start time, which is correctly ignored.
member_state_stamp() {
  local enstr="$1" case_name="${my_casename}.$1" archive dir stamp epoch best_epoch="${START_EPOCH}"
  archive=$(member_archive "${enstr}")
  for dir in "${archive}"/rest/????-??-??-?????; do
    [[ -d "${dir}" ]] || continue
    stamp=$(basename "${dir}")
    epoch=$(stamp_to_epoch "${stamp}") || continue
    (( epoch > best_epoch && epoch <= END_EPOCH )) || continue
    [[ -s "${dir}/${case_name}.eam.r.${stamp}.nc" ]] || continue
    best_epoch=${epoch}
  done
  epoch_to_stamp "${best_epoch}"
}

user_eam_nl() {
  ex user_nl_eam <<EOF
g/^ *ncdata *=/s@=.*@= "$1"@
g/^ *inithist *=/s@=.*@= 'NONE'@
wq
EOF
}

user_elm_nl() {
  ex user_nl_elm <<EOF
g/^ *finidat *=/s@=.*@= "$1"@
g/^ *check_finidat_year_consistency *=/s@=.*@= .false.@
g/^ *check_dynpft_consistency *=/s@=.*@= .false.@
g/^ *check_finidat_fsurdat_consistency *=/s@=.*@= .false.@
g/^ *check_finidat_pct_consistency *=/s@=.*@= .false.@
wq
EOF
}

user_mosart_nl() {
  ex user_nl_mosart <<EOF
g/^ *finidat_rtm *=/s@=.*@= "$1"@
wq
EOF
}

user_mpassi_nl() {
  ex user_nl_mpassi <<EOF
g/^ *config_start_time *=/s@=.*@= '$1_$2'@
g/^ *config_calendar_type *=/s@=.*@= 'gregorian'@
g/^ *config_initial_condition_type *=/s@=.*@= 'restart'@
g/^ *config_do_restart *=/s@=.*@= .true.@
g/^ *config_restart_timestamp_name *=/s@=.*@= 'rpointer.ice'@
wq
EOF
}

user_mpaso_nl() {
  ex user_nl_mpaso <<EOF
g/^ *config_start_time *=/s@=.*@= '$1_$2'@
g/^ *config_calendar_type *=/s@=.*@= 'gregorian'@
g/^ *config_do_restart *=/s@=.*@= .true.@
g/^ *config_output_reference_time *=/s@=.*@= '$3'@
g/^ *config_restart_timestamp_name *=/s@=.*@= 'rpointer.ocn'@
wq
EOF
}

# Point the MPAS initial-condition input stream at the staged restart file
# (same edit as the DA cycle worker). It is read only at a hybrid start.
patch_mpas_streams() {
  local run_dir="$1" streams="$2" src_dir="$3" ic_file="$4" case_dir="$5" line
  local -a lines=("${@:6}")
  pushd "${run_dir}" >/dev/null
  for line in "${lines[@]}"; do
    sed -i "${line}s#filename_template=.*#filename_template=\"${ic_file}\"#" "${streams}"
  done
  cp "${streams}" "${case_dir}/SourceMods/${src_dir}/"
  popd >/dev/null
}

stage_file_atomic() {
  local source_file="$1" destination_file="$2" tmp_file="$2.tmp.${SLURM_JOB_ID:-$$}"
  [[ -s "${source_file}" ]] || fail "restart source is missing or empty: ${source_file}"
  if [[ "${source_file}" == *.nc ]]; then
    ncdump -h "${source_file}" >/dev/null 2>&1 || fail "restart source is invalid NetCDF: ${source_file}"
  fi
  cp -p "${source_file}" "${tmp_file}"
  mv -f "${tmp_file}" "${destination_file}"
}

# Hybrid start from the Step 3 perturbed eam.i file and Step 2 component restarts.
prepare_member_start() (
  local enstr="$1" seg_hours="$2"
  local case_name="${my_casename}.${enstr}" case_dir run_dir archive ref_dir ymd tod hour ocn_time file
  case_dir=$(member_case_dir "${enstr}"); run_dir=$(member_run_dir "${enstr}"); archive=$(member_archive "${enstr}")
  ymd="${my_e3sm_start_date}"; tod="${my_e3sm_start_tod}"
  hour=$(printf '%02d' $((10#${tod} / 3600)))
  ocn_time="${ymd}_${hour}:00:00"
  ref_dir="${archive}/rest/${ymd}-${tod}"
  set -e

  cd "${case_dir}"
  ./xmlchange RUN_TYPE="hybrid"
  ./xmlchange CONTINUE_RUN=FALSE
  ./xmlchange RUN_STARTDATE="${ymd}"
  ./xmlchange START_TOD="${tod}"
  ./xmlchange GET_REFCASE=FALSE
  ./xmlchange RUN_REFCASE="${case_name}"
  ./xmlchange RUN_REFDATE="${ymd}"
  ./xmlchange RUN_REFTOD="${tod}"
  ./xmlchange RUN_REFDIR="${ref_dir}"

  cd "${run_dir}"
  rm -f rpointer*
  echo "${case_name}.eam.r.${ymd}-${tod}.nc" > rpointer.atm
  echo "./${case_name}.elm.r.${ymd}-${tod}.nc" > rpointer.lnd
  echo "./${case_name}.mosart.r.${ymd}-${tod}.nc" > rpointer.rof
  echo "${case_name}.cpl.r.${ymd}-${tod}.nc" > rpointer.drv
  echo "${ocn_time}" > rpointer.ice
  if [[ "${my_runtype}" == "AMIP" ]]; then
    printf '%s\n%s\n' "${case_name}.docn.r.${ymd}_${tod}.nc" "${case_name}.docn.rs1.${ymd}_${tod}.bin" > rpointer.ocn
    for file in "${case_name}.docn.r.${ymd}_${tod}.nc" "${case_name}.docn.rs1.${ymd}_${tod}.bin"; do
      [[ -s "${ref_dir}/${file}" ]] && stage_file_atomic "${ref_dir}/${file}" "${file}"
    done
  else
    echo "${ocn_time}" > rpointer.ocn
    stage_file_atomic "${ref_dir}/${case_name}.mpaso.rst.${ymd}_${tod}.nc" "${case_name}.mpaso.rst.${ymd}_${tod}.nc"
  fi
  # EAM starts from ncdata in hybrid mode; the eam.r placeholder is not staged.
  for file in "${case_name}.elm.r.${ymd}-${tod}.nc" "${case_name}.mosart.r.${ymd}-${tod}.nc" \
              "${case_name}.cpl.r.${ymd}-${tod}.nc" "${case_name}.mpassi.rst.${ymd}_${tod}.nc"; do
    stage_file_atomic "${ref_dir}/${file}" "${file}"
  done

  cd "${case_dir}"
  user_eam_nl "${ref_dir}/${case_name}.eam.i.${ymd}-${tod}.nc"
  user_elm_nl "${ref_dir}/${case_name}.elm.r.${ymd}-${tod}.nc"
  user_mosart_nl "${ref_dir}/${case_name}.mosart.r.${ymd}-${tod}.nc"
  user_mpassi_nl "${ymd}" "${hour}"
  [[ "${my_runtype}" == "Full-CPL" ]] && user_mpaso_nl "${ymd}" "${hour}" "${ocn_time}"
  set_run_length "${seg_hours}" "${archive}"
  ./case.setup
  patch_mpas_streams "${run_dir}" streams.seaice src.mpassi "${ref_dir}/${case_name}.mpassi.rst.${ymd}_${tod}.nc" "${case_dir}" 11 38
  if [[ "${my_runtype}" == "Full-CPL" ]]; then
    patch_mpas_streams "${run_dir}" streams.ocean src.mpaso "${ref_dir}/${case_name}.mpaso.rst.${ymd}_${tod}.nc" "${case_dir}" 12
  fi
  echo "Prepared hybrid start for ${enstr} at ${ymd}-${tod}"
)

# Continue from the member's latest archived restart set.
prepare_member_continue() (
  local enstr="$1" seg_hours="$2" state="$3"
  local case_name="${my_casename}.${enstr}" case_dir run_dir archive rest_dir file atm_rst
  case_dir=$(member_case_dir "${enstr}"); run_dir=$(member_run_dir "${enstr}"); archive=$(member_archive "${enstr}")
  rest_dir="${archive}/rest/${state}"
  atm_rst="${case_name}.eam.r.${state}.nc"
  set -e

  # st_archive leaves the latest restart set and matching rpointers in the run
  # directory. Restage from the archive only when that state is not intact.
  if [[ "$(head -n 1 "${run_dir}/rpointer.atm" 2>/dev/null)" != "${atm_rst}" || ! -s "${run_dir}/${atm_rst}" ]]; then
    echo "Restaging ${enstr} restart set from ${rest_dir}"
    cd "${run_dir}"
    rm -f rpointer*
    for file in "${rest_dir}"/*; do
      case "$(basename "${file}")" in
        *.h[0-9].*|*.h[0-9][0-9].*) continue ;;   # skip in-progress history copies
      esac
      stage_file_atomic "${file}" "$(basename "${file}")"
    done
    echo "${atm_rst}" > rpointer.atm
    echo "./${case_name}.elm.r.${state}.nc" > rpointer.lnd
    echo "./${case_name}.mosart.r.${state}.nc" > rpointer.rof
    echo "${case_name}.cpl.r.${state}.nc" > rpointer.drv
    echo "${state:0:10}_$(printf '%02d' $((10#${state:11:5} / 3600))):00:00" > rpointer.ice
    if [[ "${my_runtype}" == "AMIP" ]]; then
      printf '%s\n%s\n' "${case_name}.docn.r.${state:0:10}_${state:11:5}.nc" "${case_name}.docn.rs1.${state:0:10}_${state:11:5}.bin" > rpointer.ocn
    else
      cp rpointer.ice rpointer.ocn
    fi
  fi

  cd "${case_dir}"
  ./xmlchange CONTINUE_RUN=TRUE
  set_run_length "${seg_hours}" "${archive}"
  ./case.setup
  echo "Prepared continuation for ${enstr} from ${state}"
)

set_run_length() {
  ./xmlchange STOP_OPTION=nhours,STOP_N="$1"
  ./xmlchange REST_OPTION=nhours,REST_N="$1"
  ./xmlchange RESUBMIT=0
  ./xmlchange DOUT_S=TRUE
  ./xmlchange DOUT_S_ROOT="$2"
}

capture_member_failure() {
  local enstr="$1" run_dir log diag
  run_dir=$(member_run_dir "${enstr}")
  diag="${LOG_DIR}/failure.${enstr}.${SLURM_JOB_ID:-$$}.free.log"
  log=$(find "${run_dir}" -maxdepth 1 -type f -name 'e3sm.log.*' -printf '%T@ %p\n' 2>/dev/null | sort -nr | awk 'NR==1 {print $2}')
  {
    echo "Failure diagnostics for ${enstr} at $(date '+%Y-%m-%d %H:%M:%S')"
    tail -n 40 "$(member_case_dir "${enstr}")/CaseStatus" 2>/dev/null || true
    if [[ -n "${log}" ]]; then
      echo "===== ${log} ====="
      grep -nEi 'error|fatal|abort|segmentation|killed|bad termination' "${log}" | tail -n 60 || true
      tail -n 100 "${log}"
    fi
  } > "${diag}" 2>&1
  echo "Captured ${enstr} failure diagnostics: ${diag}"
}

remaining_allocation_seconds() {
  local value days=0 hours minutes seconds
  [[ -n "${SLURM_JOB_ID:-}" ]] || return 1
  value=$(squeue -h -j "${SLURM_JOB_ID}" -o '%L') || return 1
  value=${value//[[:space:]]/}
  [[ -n "${value}" && "${value}" != "UNLIMITED" && "${value}" != "NOT_SET" ]] || return 1
  if [[ "${value}" == *-* ]]; then days=${value%%-*}; value=${value#*-}; fi
  IFS=: read -r hours minutes seconds <<< "${value}"
  if [[ -z "${seconds:-}" ]]; then seconds=${minutes}; minutes=${hours}; hours=0; fi
  echo $((10#${days} * 86400 + 10#${hours} * 3600 + 10#${minutes} * 60 + 10#${seconds}))
}

submit_continuation() {
  local job_name existing next_job_id
  job_name="e3sm-free-$(printf '%s' "${my_wkdir}" | cksum | awk '{print $1}')"
  existing=$(squeue -h -u "${USER}" -n "${job_name}" -o '%A' | awk -v cur="${SLURM_JOB_ID:-none}" 'NF && $1 != cur {print; exit}')
  if [[ -n "${existing}" ]]; then
    echo "Continuation already queued: job ${existing}"
    return 0
  fi
  next_job_id=$(sbatch --parsable --job-name="${job_name}" --dependency="afterany:${SLURM_JOB_ID}" "${SCRIPT_PATH}")
  echo "Submitted continuation job ${next_job_id%%;*}"
}

#------------------------------------------------------------------------------
# Completion: verify history, prune intermediate restarts, write the record.
#------------------------------------------------------------------------------
# With the template's 6-hourly streams (mfilt=4), EAM/ELM open a new history
# file at 00Z each day: file D holds the samples valid at D 00/06/12/18Z (the
# first file also holds the nstep-0 sample), and one last file is opened at the
# end time. Expect one file per day from the start through the end time.
FREE_HIST_STREAMS=("atm:eam.h0" "atm:eam.h1" "atm:eam.h2" "lnd:elm.h0" "lnd:elm.h1")
finalize_run() {
  local i enstr case_name archive epoch stamp missing=0 dir dir_epoch stream file
  echo "All members reached ${END_STAMP}; verifying archived history"
  for ((i=1; i<=my_ensnum; i++)); do
    enstr=$(member_name "${i}"); case_name="${my_casename}.${enstr}"; archive=$(member_archive "${enstr}")
    for ((epoch=START_EPOCH; epoch<=END_EPOCH; epoch+=86400)); do
      stamp=$(epoch_to_stamp "${epoch}")
      for stream in "${FREE_HIST_STREAMS[@]}"; do
        file="${archive}/${stream%%:*}/hist/${case_name}.${stream#*:}.${stamp}.nc"
        [[ -s "${file}" ]] || { echo "ERROR: missing history file: ${file}" >&2; missing=$((missing + 1)); }
      done
    done
  done
  (( missing == 0 )) || fail "${missing} history file(s) missing; restart sets were kept for recovery"

  if [[ "${KEEP_INTERMEDIATE}" == "FALSE" ]]; then
    for ((i=1; i<=my_ensnum; i++)); do
      enstr=$(member_name "${i}"); archive=$(member_archive "${enstr}")
      for dir in "${archive}"/rest/????-??-??-?????; do
        [[ -d "${dir}" ]] || continue
        dir_epoch=$(stamp_to_epoch "$(basename "${dir}")") || continue
        if (( dir_epoch > START_EPOCH && dir_epoch < END_EPOCH )); then
          echo "Removing intermediate restart set: ${dir}"
          rm -rf -- "${dir}"
        fi
      done
    done
  fi

  printf 'valid_time=%s\nstart_time=%s\ncase=%s\nensemble_size=%s\nslurm_job_id=%s\ncompleted_at=%s\n' \
    "${END_STAMP}" "${START_STAMP}" "${my_casename}" "${my_ensnum}" "${SLURM_JOB_ID:-none}" "$(date '+%Y-%m-%d %H:%M:%S')" > "${COMPLETE_MARKER}.tmp"
  mv -f "${COMPLETE_MARKER}.tmp" "${COMPLETE_MARKER}"
  echo "Wrote completion record: ${COMPLETE_MARKER}"
}

#------------------------------------------------------------------------------
# Main loop: one segment per pass, all pending members in parallel waves.
#------------------------------------------------------------------------------
SCRIPT_NAME=$(basename "${SCRIPT_PATH}")
RUN_START_EPOCH=$(date +%s)
trap 'rc=$?; echo "INFO: ${SCRIPT_NAME} exited with status ${rc} after $(( $(date +%s) - RUN_START_EPOCH ))s"' EXIT
on_signal() {
  echo "ERROR: ${SCRIPT_NAME} received a termination signal"
  jobs -pr | while IFS= read -r pid; do kill -TERM "${pid}" 2>/dev/null || true; done
  wait 2>/dev/null || true
  exit 143
}
trap on_signal INT TERM

cd "${my_wkdir}"
echo "== Start of ${SCRIPT_NAME}: ${my_casename}, ${START_STAMP} -> ${END_STAMP}, ${SEGMENT_DAYS}-day segments =="
date
if [[ -s "${COMPLETE_MARKER}" ]]; then
  echo "Free run already complete: ${COMPLETE_MARKER}"
  exit 0
fi

segment_index=0
segments_run=0
while true; do
  segment_index=$((segment_index + 1))
  declare -A member_state=() member_target=()
  pending=()
  for ((i=1; i<=my_ensnum; i++)); do
    enstr=$(member_name "${i}")
    [[ -x "$(member_case_dir "${enstr}")/xmlchange" ]] || fail "missing member case from Step 1: $(member_case_dir "${enstr}")"
    state=$(member_state_stamp "${enstr}")
    state_epoch=$(stamp_to_epoch "${state}")
    if (( state_epoch < END_EPOCH )); then
      target_epoch=$((state_epoch + SEGMENT_DAYS * 86400))
      (( target_epoch > END_EPOCH )) && target_epoch=${END_EPOCH}
      member_state[${enstr}]=${state}
      member_target[${enstr}]=$(epoch_to_stamp "${target_epoch}")
      pending+=("${enstr}")
    fi
  done

  if (( ${#pending[@]} == 0 )); then
    finalize_run
    break
  fi

  required_sec=$((MIN_SEGMENT_TIME_SEC + SHUTDOWN_MARGIN_SEC))
  if remaining_sec=$(remaining_allocation_seconds); then
    echo "Remaining allocation: ${remaining_sec}s; required per segment: ${required_sec}s"
    if (( remaining_sec < required_sec )); then
      # Resubmit only after progress; an allocation too short for even one
      # segment would otherwise resubmit itself forever.
      (( segments_run > 0 )) || fail "job time limit leaves ${remaining_sec}s but one ${SEGMENT_DAYS}-day segment needs ${required_sec}s; increase #SBATCH --time or lower my_free_segment_days and my_free_min_segment_time_sec"
      echo "Insufficient time for another segment; resubmitting"
      submit_continuation
      break
    fi
  elif [[ -n "${SLURM_JOB_ID:-}" ]]; then
    (( segments_run > 0 )) || fail "unable to determine remaining Slurm time before the first segment"
    echo "WARNING: unable to determine remaining Slurm time; resubmitting conservatively"
    submit_continuation
    break
  fi

  echo "== Segment pass ${segment_index}: ${#pending[@]} member(s) =="
  setup_pids=(); setup_names=()
  for enstr in "${pending[@]}"; do
    state=${member_state[${enstr}]}
    seg_hours=$(( ($(stamp_to_epoch "${member_target[${enstr}]}") - $(stamp_to_epoch "${state}")) / 3600 ))
    echo "${enstr}: ${state} -> ${member_target[${enstr}]} (${seg_hours} h)"
    setup_log="${LOG_DIR}/free_setup.${enstr}.${SLURM_JOB_ID:-$$}.${state}.log"
    if [[ "${state}" == "${START_STAMP}" ]]; then
      prepare_member_start "${enstr}" "${seg_hours}" > "${setup_log}" 2>&1 &
    else
      prepare_member_continue "${enstr}" "${seg_hours}" "${state}" > "${setup_log}" 2>&1 &
    fi
    setup_pids+=("$!"); setup_names+=("${enstr}")
    if (( ${#setup_pids[@]} >= MAX_PARALLEL_SETUP )) || [[ "${enstr}" == "${pending[-1]}" ]]; then
      for idx in "${!setup_pids[@]}"; do
        wait "${setup_pids[$idx]}" || fail "setup failed for ${setup_names[$idx]}; see ${LOG_DIR}/free_setup.${setup_names[$idx]}.${SLURM_JOB_ID:-$$}.*.log"
      done
      setup_pids=(); setup_names=()
    fi
  done

  failed=()
  member_pids=(); member_names=()
  for enstr in "${pending[@]}"; do
    (cd "$(member_case_dir "${enstr}")" && ./case.submit --no-batch) \
      > "${LOG_DIR}/e3sm.${enstr}.o${SLURM_JOB_ID:-$$}.free.${member_state[${enstr}]}.log" 2>&1 &
    member_pids+=("$!"); member_names+=("${enstr}")
    if (( ${#member_pids[@]} >= BATCH_SIZE )) || [[ "${enstr}" == "${pending[-1]}" ]]; then
      echo "Waiting for members: ${member_names[*]}"
      for idx in "${!member_pids[@]}"; do
        wait "${member_pids[$idx]}" || echo "WARNING: case.submit returned non-zero for ${member_names[$idx]}"
      done
      member_pids=(); member_names=()
    fi
  done

  # A member succeeded only if its segment-end restart reached the archive.
  for enstr in "${pending[@]}"; do
    target=${member_target[${enstr}]}
    case_name="${my_casename}.${enstr}"
    archived="$(member_archive "${enstr}")/rest/${target}/${case_name}.eam.r.${target}.nc"
    if [[ ! -s "${archived}" && -s "$(member_run_dir "${enstr}")/${case_name}.eam.r.${target}.nc" ]]; then
      echo "Running short-term archive for ${enstr}"
      (cd "$(member_case_dir "${enstr}")" && srun -N 1 -n 1 ./case.st_archive) || true
    fi
    if [[ -s "${archived}" ]]; then
      echo "${enstr} reached ${target}"
    else
      capture_member_failure "${enstr}"
      failed+=("${enstr}")
    fi
  done
  (( ${#failed[@]} == 0 )) || fail "segment failed for: ${failed[*]}; resubmit this script to resume from the last archived restart"
  segments_run=$((segments_run + 1))
done

echo "== End of ${SCRIPT_NAME} =="
date
