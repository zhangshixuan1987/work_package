#!/bin/bash -l 
#SBATCH --account=esmd
#SBATCH --time=24:00:00
#SBATCH --partition=slurm
#SBATCH --job-name=e3sm_dart_post
#SBATCH --nodes=1
#SBATCH --output=e3sm_dart_post.%j
#SBATCH --exclusive
#SBATCH --no-kill
#SBATCH --requeue

set -o pipefail

# User settings: base, monthly, or all.
POST_MODE="all"
POST_START="2011-12-01-00000"
POST_END="2011-12-31-00000"
MAX_CONCURRENT_WORKERS=4
OVERWRITE_EXISTING="FALSE"

RUN_EAM_6HOURLY="TRUE"
RUN_EAM_DAILY="TRUE"
RUN_EAM_CLIM="TRUE"
RUN_ELM_DAILY="TRUE"
RUN_ELM_CLIM="TRUE"
RUN_EAM_MONTHLY="TRUE"
RUN_ELM_MONTHLY="TRUE"

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
[[ "${POST_MODE}" =~ ^(base|monthly|all)$ ]] || fail "POST_MODE must be base, monthly, or all"
[[ "${MAX_CONCURRENT_WORKERS}" =~ ^[1-9][0-9]*$ ]] || fail "MAX_CONCURRENT_WORKERS must be positive"
for setting in OVERWRITE_EXISTING RUN_EAM_6HOURLY RUN_EAM_DAILY RUN_EAM_CLIM RUN_ELM_DAILY RUN_ELM_CLIM RUN_EAM_MONTHLY RUN_ELM_MONTHLY; do
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
  local member case_name archive epoch stamp component stream file requirement
  local month_epoch month_end_epoch month missing=0
  local -a cycle_streams=() monthly_streams=() required_stamps=() required_months=()

  # Post-processing depends on archived history, not on DA cycle status files.
  if [[ "${POST_MODE}" == "base" || "${POST_MODE}" == "all" ]]; then
    [[ "${RUN_EAM_6HOURLY}" == "TRUE" ]] && cycle_streams+=("atm:eam.h2")
    [[ "${RUN_EAM_DAILY}" == "TRUE" ]] && cycle_streams+=("atm:eam.h1")
    [[ "${RUN_ELM_DAILY}" == "TRUE" ]] && cycle_streams+=("lnd:elm.h1")
    [[ "${RUN_EAM_CLIM}" == "TRUE" ]] && monthly_streams+=("atm:eam.h0")
    [[ "${RUN_ELM_CLIM}" == "TRUE" ]] && monthly_streams+=("lnd:elm.h0")
  fi

  month_epoch=$(date -u -d "${POST_START_DATE%-*}-01" +%s)
  month_end_epoch=$(date -u -d "${POST_END_DATE%-*}-01" +%s)
  for ((epoch=START_EPOCH; epoch<=END_EPOCH+18*3600; epoch+=STEP_SECONDS)); do
    required_stamps+=("$(epoch_to_stamp "${epoch}")")
  done
  for ((epoch=month_epoch; epoch<=month_end_epoch; )); do
    month=$(date -u -d "@${epoch}" +%Y-%m)
    required_months+=("${month}")
    epoch=$(date -u -d "${month}-01 +1 month" +%s)
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

    for requirement in "${monthly_streams[@]}"; do
      component="${requirement%%:*}"
      stream="${requirement#*:}"
      for month in "${required_months[@]}"; do
        file="${archive}/${component}/hist/${case_name}.${stream}.${month}.nc"
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
validate_history_inputs

product_worker() {
  case "$1" in
    eam_6hourly) echo "${my_workflow_lib}/post/post_process_eam_6hourly.sh" ;;
    eam_daily) echo "${my_workflow_lib}/post/post_process_eam_daily.sh" ;;
    eam_clim) echo "${my_workflow_lib}/post/post_process_eam_monthly.sh" ;;
    elm_daily) echo "${my_workflow_lib}/post/post_process_elm_daily.sh" ;;
    elm_clim) echo "${my_workflow_lib}/post/post_process_elm_monthly.sh" ;;
    eam_monthly) echo "${my_workflow_lib}/post/hfreq_to_eam_monthly.sh" ;;
    elm_monthly) echo "${my_workflow_lib}/post/hfreq_to_elm_monthly.sh" ;;
    *) return 1 ;;
  esac
}

validate_product_outputs() {
  local product="$1" member="$2" root pattern found=0 found_da=0 file member_archive
  member_archive="${my_modeldir}/${member}/archive"
  case "${product}" in
    eam_6hourly) root="${member_archive}/post/atm/180x360_aave/ts/6hourly"; pattern="*.${member}.*.nc" ;;
    eam_daily) root="${member_archive}/post/atm/180x360_aave/ts/daily"; pattern="*.${member}.*.nc" ;;
    eam_clim) root="${member_archive}/post/atm/180x360_aave/clim"; pattern="*.${member}.*.nc" ;;
    elm_daily) root="${member_archive}/post/lnd/180x360_aave/ts/daily"; pattern="*.${member}.*.nc" ;;
    elm_clim) root="${member_archive}/post/lnd/180x360_aave/clim"; pattern="*.${member}.*.nc" ;;
    eam_monthly) root="${member_archive}/post/atm/180x360_aave/monthly"; pattern="*.${member}.*.nc" ;;
    elm_monthly) root="${member_archive}/post/lnd/180x360_aave/monthly"; pattern="*.${member}.*.nc" ;;
  esac
  [[ -d "${root}" ]] || return 1
  while IFS= read -r file; do
    [[ "${file}" == *"${START_YEAR}"* || "${file}" == *"${END_YEAR}"* ]] || continue
    [[ -s "${file}" ]] || return 1
    ncdump -h "${file}" >/dev/null 2>&1 || return 1
    found=$((found+1))
  done < <(find "${root}" -type f -name "${pattern}" -print | sort)
  if [[ "${product}" == "eam_6hourly" ]]; then
    root="${member_archive}/post/atm/180x360_aave/ts/6hourly_da"
    [[ -d "${root}" ]] || return 1
    while IFS= read -r file; do
      [[ "${file}" == *"${START_YEAR}"* || "${file}" == *"${END_YEAR}"* ]] || continue
      [[ -s "${file}" ]] || return 1
      ncdump -h "${file}" >/dev/null 2>&1 || return 1
      found_da=$((found_da+1))
    done < <(find "${root}" -type f -name "${pattern}" -print | sort)
    (( found_da > 0 )) || return 1
  fi
  (( found > 0 ))
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
  POSTPROCESS_DRIVER_ACTIVE=TRUE POST_ENSTR="${member}" POST_CASE_NAME="${case_name}" \
    POST_START_DATE="${POST_START_DATE}" POST_END_DATE="${POST_END_DATE}" "${worker}"
  validate_product_outputs "${product}" "${member}" || fail "invalid or missing ${product} output for ${member}"

  tmp="${record}.tmp.${SLURM_JOB_ID:-$$}"
  printf 'product=%s\nmember=%s\nstart=%s\nend=%s\ncase=%s\nensemble_size=%s\njob_id=%s\ncompleted_at=%s\n' \
    "${product}" "${member}" "${POST_START}" "${POST_END}" "${my_casename}" "${my_ensnum}" \
    "${SLURM_JOB_ID:-none}" "$(date '+%Y-%m-%d %H:%M:%S')" > "${tmp}"
  mv -f "${tmp}" "${record}"
  rm -f -- "${progress}"
}

run_product_batch() {
  local product="$1" failed=0 member pid index status worker_limit
  local -a pids=() members=()
  worker_limit="${MAX_CONCURRENT_WORKERS}"
  if [[ "${product}" == "elm_clim" ]]; then
    worker_limit=1
    echo "Serializing elm_clim because ncremap -P elm uses internal SGS scratch files"
  fi
  for i in $(seq 1 "${my_ensnum}"); do
    member=$(printf 'EN%02d' "${i}")
    run_product "${product}" "${member}" &
    pids+=("$!"); members+=("${member}")
    if (( ${#pids[@]} >= worker_limit )); then
      for index in "${!pids[@]}"; do wait "${pids[$index]}" || { status=$?; echo "ERROR: ${product} failed for ${members[$index]} (${status})" >&2; failed=1; }; done
      pids=(); members=()
    fi
  done
  for index in "${!pids[@]}"; do wait "${pids[$index]}" || { status=$?; echo "ERROR: ${product} failed for ${members[$index]} (${status})" >&2; failed=1; }; done
  (( failed == 0 )) || fail "${product} failed for one or more members"
}

BASE_PRODUCTS=()
[[ "${RUN_EAM_6HOURLY}" == "TRUE" ]] && BASE_PRODUCTS+=(eam_6hourly)
[[ "${RUN_EAM_DAILY}" == "TRUE" ]] && BASE_PRODUCTS+=(eam_daily)
[[ "${RUN_EAM_CLIM}" == "TRUE" ]] && BASE_PRODUCTS+=(eam_clim)
[[ "${RUN_ELM_DAILY}" == "TRUE" ]] && BASE_PRODUCTS+=(elm_daily)
[[ "${RUN_ELM_CLIM}" == "TRUE" ]] && BASE_PRODUCTS+=(elm_clim)
MONTHLY_PRODUCTS=()
[[ "${RUN_EAM_MONTHLY}" == "TRUE" ]] && MONTHLY_PRODUCTS+=(eam_monthly)
[[ "${RUN_ELM_MONTHLY}" == "TRUE" ]] && MONTHLY_PRODUCTS+=(elm_monthly)

if [[ "${POST_MODE}" == "base" || "${POST_MODE}" == "all" ]]; then
  (( ${#BASE_PRODUCTS[@]} > 0 )) || fail "no base products enabled"
  for product in "${BASE_PRODUCTS[@]}"; do run_product_batch "${product}"; done
fi
if [[ "${POST_MODE}" == "monthly" || "${POST_MODE}" == "all" ]]; then
  [[ "${RUN_EAM_MONTHLY}" != "TRUE" || "${RUN_EAM_DAILY}" == "TRUE" || "${RUN_EAM_6HOURLY}" == "TRUE" || "${POST_MODE}" == "monthly" ]] || fail "EAM monthly processing requires EAM daily or 6-hourly base processing"
  [[ "${RUN_ELM_MONTHLY}" != "TRUE" || "${RUN_ELM_DAILY}" == "TRUE" || "${POST_MODE}" == "monthly" ]] || fail "ELM monthly processing requires ELM daily base processing"
  (( ${#MONTHLY_PRODUCTS[@]} > 0 )) || fail "no monthly products enabled"
  if [[ "${POST_MODE}" == "monthly" ]]; then
    for i in $(seq 1 "${my_ensnum}"); do
      member=$(printf 'EN%02d' "${i}")
      [[ "${RUN_EAM_MONTHLY}" != "TRUE" ]] || validate_completion_record eam_daily "${member}" || validate_completion_record eam_6hourly "${member}" || fail "missing or invalid EAM base completion for ${member}"
      [[ "${RUN_ELM_MONTHLY}" != "TRUE" ]] || validate_completion_record elm_daily "${member}" || fail "missing or invalid ELM base completion for ${member}"
    done
  fi
  for product in "${MONTHLY_PRODUCTS[@]}"; do run_product_batch "${product}"; done
fi

echo "Post-processing completed: mode=${POST_MODE}, range=${POST_START} through ${POST_END}"
