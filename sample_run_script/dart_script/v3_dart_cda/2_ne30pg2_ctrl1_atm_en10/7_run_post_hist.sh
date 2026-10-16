#!/bin/bash -l
#SBATCH --account=esmd
#SBATCH --time=2:00:00
#SBATCH --partition=short
#SBATCH --job-name=e3sm_dart_post
#SBATCH --nodes=1
#SBATCH --output=e3sm_dart_post.%j
#SBATCH --exclusive
#SBATCH --no-kill
#SBATCH --requeue

set -o pipefail

# Step 7: regridded history products for every member, per component.
#
# DA-cycled runs (eam_cycled, elm_cycled; workflow_lib/post/cycled_split.py)
# write, under ENxx/archive/post/<atm|lnd>/180x360_aave/:
#   ts/6hourly           forecast: 6-h means/max/min (window middle, with
#                        bounds) and background states (window end)
#   ts/6hourly_analysis  analysis state at every cycle start
#   ts/daily_mean        calendar-day mean: four 6-h forecast means, or the
#                        mean of the 00/06/12/18Z analysis states
#   ts/daily_max, ts/daily_min, ts/daily_00z
#   monthly/             mean of ts/daily_mean per month
# Products of the earlier version (ts/6hourly_da, ts/daily, clim, ...) are
# moved to 180x360_aave/legacy/ on the first run.
# Free runs (my_free_layout set; eam_hfreq, elm_hfreq; free_run_hfreq.py) write
# the same layout without ts/6hourly_analysis (a free run has no analysis):
# ts/6hourly, ts/daily_mean, ts/daily_max, ts/daily_min, ts/daily_00z, monthly/.

# User settings
POST_START="2012-05-01-00000"
POST_END="2012-05-31-00000"     # last day processed (inclusive), midnight stamp
MAX_CONCURRENT_WORKERS=4
OVERWRITE_EXISTING="FALSE"
RUN_EAM="TRUE"
RUN_ELM="TRUE"
# Cycled runs: tapes allowed to have no 6-h records (analysis states only).
ANALYSIS_ONLY_TAPES=""  # none: every tape must hold the 6-h records

fail() { echo "ERROR: $*" >&2; exit 1; }
normalize_bool() {
  local name="$1" value="${2^^}"
  [[ "${value}" == "TRUE" || "${value}" == "FALSE" ]] || fail "${name} must be TRUE or FALSE, got: $2"
  printf '%s' "${value}"
}

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  WORK_DIR=$(readlink -f "${SLURM_SUBMIT_DIR}") || fail "cannot resolve SLURM_SUBMIT_DIR"
else
  SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || fail "cannot resolve script path"
  WORK_DIR=$(dirname "${SCRIPT_PATH}")
fi
cd "${WORK_DIR}"
[[ -r create_and_setup_case.sh ]] || fail "missing workflow configuration"
source ./create_and_setup_case.sh

[[ -n "${my_conda_setup_file:-}" ]] || fail "my_conda_setup_file is unset"
[[ -r "${my_conda_setup_file}" ]] || fail "configured Conda setup is not readable: ${my_conda_setup_file}"
[[ -n "${my_analysis_conda_env:-}" ]] || fail "my_analysis_conda_env is unset"

mkdir -p "${my_log_dir}" "${my_status_dir}" "${my_lock_dir}"

for cmd in awk bash date find flock mkdir mv ncdump readlink sort; do
  command -v "${cmd}" >/dev/null 2>&1 || fail "required command not found: ${cmd}"
done
[[ "${MAX_CONCURRENT_WORKERS}" =~ ^[1-9][0-9]*$ ]] || fail "MAX_CONCURRENT_WORKERS must be positive"
for setting in OVERWRITE_EXISTING RUN_EAM RUN_ELM; do
  printf -v "${setting}" '%s' "$(normalize_bool "${setting}" "${!setting}")"
done

stamp_to_epoch() {
  local stamp="$1" ymd tod h m s
  [[ "${stamp}" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}-[0-9]{5}$ ]] || return 1
  ymd="${stamp%-*}"; tod="${stamp##*-}"
  (( 10#${tod} < 86400 )) || return 1
  h=$((10#${tod}/3600)); m=$(((10#${tod}%3600)/60)); s=$((10#${tod}%60))
  date -u -d "${ymd} $(printf '%02d:%02d:%02d' "${h}" "${m}" "${s}")" +%s
}
epoch_to_stamp() {
  local epoch="$1" ymd h m s tod
  read -r ymd h m s < <(date -u -d "@${epoch}" '+%Y-%m-%d %H %M %S')
  tod=$((10#${h}*3600+10#${m}*60+10#${s}))
  printf '%s-%05d\n' "${ymd}" "${tod}"
}
START_EPOCH=$(stamp_to_epoch "${POST_START}") || fail "invalid POST_START: ${POST_START}"
[[ "${POST_START##*-}" == "00000" && "${POST_END##*-}" == "00000" ]] || fail "POST_START and POST_END must be midnight timestamps; workers process whole calendar days"
END_EPOCH=$(stamp_to_epoch "${POST_END}") || fail "invalid POST_END: ${POST_END}"
(( START_EPOCH <= END_EPOCH )) || fail "POST_START is after POST_END"
STEP_SECONDS=$((my_e3sm_cycle_hours*3600))
(( (END_EPOCH-START_EPOCH)%STEP_SECONDS == 0 )) || fail "post-processing range is not aligned to the DA interval"
POST_START_DATE="${POST_START%-*}"
POST_END_DATE="${POST_END%-*}"
START_YEAR="${POST_START_DATE%%-*}"
END_YEAR="${POST_END_DATE%%-*}"

STATUS_DIR="${my_status_dir}"
mkdir -p "${STATUS_DIR}"
exec 8>"${my_lock_dir}/postprocess.lock"
flock -n 8 || fail "another post-processing driver is active"

validate_history_inputs() {
  local member case_name archive epoch stamp component stream file requirement missing=0
  local -a cycle_streams=() required_stamps=()

  # Post-processing depends on archived history, not on DA cycle status files.
  # Every cycle of the last day (through 18Z) closes that day's forecast windows.
  [[ "${RUN_EAM}" == "TRUE" ]] && cycle_streams+=("atm:eam.h1" "atm:eam.h2")
  [[ "${RUN_ELM}" == "TRUE" ]] && cycle_streams+=("lnd:elm.h1")
  for ((epoch=START_EPOCH; epoch<=END_EPOCH+24*3600-STEP_SECONDS; epoch+=STEP_SECONDS)); do
    required_stamps+=("$(epoch_to_stamp "${epoch}")")
  done
  echo "Checking archived history inputs for ${my_ensnum} members"
  for i in $(seq 1 "${my_ensnum}"); do
    member=$(printf 'EN%02d' "${i}")
    case_name="${my_casename}.${member}"
    archive="${my_modeldir}/${member}/archive"
    for requirement in "${cycle_streams[@]}"; do
      component="${requirement%%:*}"
      stream="${requirement#*:}"
      for stamp in "${required_stamps[@]}"; do
        file="${archive}/${component}/hist/${case_name}.${stream}.${stamp}.nc"
        if [[ ! -s "${file}" ]]; then
          echo "ERROR: missing or empty history input: ${file}" >&2
          missing=$((missing+1))
        fi
      done
    done
  done
  (( missing == 0 )) || fail "history input preflight found ${missing} missing or empty file(s)"
  echo "Archived history input preflight passed"
}
# Free-run cases (Step 4 = 4_run_e3sm_continuous.sh) write every tape 6-hourly,
# four records per daily file, instead of one file per DA cycle; their
# products (eam_hfreq, elm_hfreq) cover every history variable.
FREE_RUN="FALSE"
[[ -n "${my_free_layout:-}" ]] && FREE_RUN="TRUE"
if [[ "${FREE_RUN}" == "TRUE" ]]; then
  # finalize already verified every archived history file of the run.
  free_marker="${my_status_dir}/free_run_complete.${my_e3sm_start_date}-${my_e3sm_start_tod}-${my_e3sm_end_date}-${my_e3sm_end_tod}"
  [[ -s "${free_marker}" ]] || fail "free run not finalized (${free_marker}); run ./4_run_e3sm_continuous.sh finalize"
  [[ ! "${POST_START_DATE}" < "${my_e3sm_start_date}" && "${POST_END_DATE}" < "${my_e3sm_end_date}" ]] ||
    fail "POST_START..POST_END (${POST_START_DATE}..${POST_END_DATE}) must lie within the run, ${my_e3sm_start_date} to the day before ${my_e3sm_end_date}"
  echo "Free run: history verified by ${free_marker}"
else
  validate_history_inputs
fi

product_worker() {
  case "$1" in
    eam_cycled|elm_cycled) echo "${my_workflow_lib}/post/post_process_cycled.sh" ;;
    eam_hfreq|elm_hfreq) echo "${my_workflow_lib}/post/post_process_free_hfreq.sh" ;;
    *) return 1 ;;
  esac
}

manifest_path() { echo "${STATUS_DIR}/postprocess_manifest.$1.$2.${POST_START}-${POST_END}"; }

validate_product_outputs() {
  local product="$1" member="$2" file found=0 manifest
  case "${product}" in
    eam_cycled|elm_cycled|eam_hfreq|elm_hfreq)
      # Every file the worker listed must exist and be readable.
      manifest=$(manifest_path "${product}" "${member}")
      [[ -s "${manifest}" ]] || return 1
      while IFS= read -r file; do
        [[ -s "${file}" ]] || return 1
        ncdump -h "${file}" >/dev/null 2>&1 || return 1
        found=$((found+1))
      done < "${manifest}"
      (( found > 0 )) ;;
    *) return 1 ;;
  esac
}

validate_completion_record() {
  local product="$1" member="$2" record
  record="${STATUS_DIR}/postprocess_complete.${product}.${member}.${POST_START}-${POST_END}"
  [[ -s "${record}" ]] || return 1
  [[ "$(awk -F= '$1 == "product" {print $2; exit}' "${record}")" == "${product}" ]] || return 1
  [[ "$(awk -F= '$1 == "member" {print $2; exit}' "${record}")" == "${member}" ]] || return 1
  [[ "$(awk -F= '$1 == "start" {print $2; exit}' "${record}")" == "${POST_START}" ]] || return 1
  [[ "$(awk -F= '$1 == "end" {print $2; exit}' "${record}")" == "${POST_END}" ]] || return 1
  [[ "$(awk -F= '$1 == "case" {print $2; exit}' "${record}")" == "${my_casename}" ]] || return 1
  [[ "$(awk -F= '$1 == "ensemble_size" {print $2; exit}' "${record}")" == "${my_ensnum}" ]] || return 1
  validate_product_outputs "${product}" "${member}"
}

run_product() {
  local product="$1" member="$2" case_name worker record progress tmp
  case_name="${my_casename}.${member}"
  worker=$(product_worker "${product}") || fail "unknown product: ${product}"
  [[ -x "${worker}" ]] || fail "missing post-processing worker: ${worker}"
  record="${STATUS_DIR}/postprocess_complete.${product}.${member}.${POST_START}-${POST_END}"
  progress="${STATUS_DIR}/postprocess_in_progress.${product}.${member}.${POST_START}-${POST_END}"
  exec 9>"${my_lock_dir}/postprocess.${product}.${member}.lock"
  flock -n 9 || fail "another ${product} worker is active for ${member}"

  if [[ "${OVERWRITE_EXISTING}" == "FALSE" ]] && validate_completion_record "${product}" "${member}"; then
    echo "Skipping verified ${product} for ${member}"
    return 0
  fi

  rm -f -- "${record}"
  printf 'product=%s\nmember=%s\nstart=%s\nend=%s\ncase=%s\nensemble_size=%s\njob_id=%s\nstarted_at=%s\n' \
    "${product}" "${member}" "${POST_START}" "${POST_END}" "${my_casename}" "${my_ensnum}" \
    "${SLURM_JOB_ID:-none}" "$(date '+%Y-%m-%d %H:%M:%S')" > "${progress}"

  echo "Running ${product} for ${member}"
  POSTPROCESS_DRIVER_ACTIVE=TRUE POST_PRODUCT="${product}" OVERWRITE_EXISTING="${OVERWRITE_EXISTING}" \
    POST_ENSTR="${member}" POST_CASE_NAME="${case_name}" POST_ANALYSIS_ONLY_TAPES="${ANALYSIS_ONLY_TAPES}" POST_MANIFEST="$(manifest_path "${product}" "${member}")" \
    POST_START_DATE="${POST_START_DATE}" POST_END_DATE="${POST_END_DATE}" "${worker}" ||
    fail "${product} worker exited with an error for ${member}"
  validate_product_outputs "${product}" "${member}" || fail "invalid or missing ${product} output for ${member}"

  tmp="${record}.tmp.${SLURM_JOB_ID:-$$}"
  printf 'product=%s\nmember=%s\nstart=%s\nend=%s\ncase=%s\nensemble_size=%s\njob_id=%s\ncompleted_at=%s\n' \
    "${product}" "${member}" "${POST_START}" "${POST_END}" "${my_casename}" "${my_ensnum}" \
    "${SLURM_JOB_ID:-none}" "$(date '+%Y-%m-%d %H:%M:%S')" > "${tmp}"
  mv -f "${tmp}" "${record}"
  rm -f -- "${progress}"
}

run_product_batch() {
  local product="$1" failed=0 member index status
  local -a pids=() members=()
  for i in $(seq 1 "${my_ensnum}"); do
    member=$(printf 'EN%02d' "${i}")
    run_product "${product}" "${member}" &
    pids+=("$!"); members+=("${member}")
    if (( ${#pids[@]} >= MAX_CONCURRENT_WORKERS )); then
      for index in "${!pids[@]}"; do wait "${pids[$index]}" || { status=$?; echo "ERROR: ${product} failed for ${members[$index]} (${status})" >&2; failed=1; }; done
      pids=(); members=()
    fi
  done
  for index in "${!pids[@]}"; do wait "${pids[$index]}" || { status=$?; echo "ERROR: ${product} failed for ${members[$index]} (${status})" >&2; failed=1; }; done
  (( failed == 0 )) || fail "${product} failed for one or more members"
}

# Every product writes its own daily and monthly series.
PRODUCTS=()
if [[ "${FREE_RUN}" == "TRUE" ]]; then
  [[ "${RUN_EAM}" == "TRUE" ]] && PRODUCTS+=(eam_hfreq)
  [[ "${RUN_ELM}" == "TRUE" ]] && PRODUCTS+=(elm_hfreq)
else
  [[ "${RUN_EAM}" == "TRUE" ]] && PRODUCTS+=(eam_cycled)
  [[ "${RUN_ELM}" == "TRUE" ]] && PRODUCTS+=(elm_cycled)
fi
(( ${#PRODUCTS[@]} > 0 )) || fail "no products enabled"
for product in "${PRODUCTS[@]}"; do run_product_batch "${product}"; done

echo "Post-processing completed: products=${PRODUCTS[*]}, range=${POST_START} through ${POST_END}"
