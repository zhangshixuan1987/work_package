#!/bin/bash -l
#SBATCH --account=esmd
#SBATCH --time=24:00:00
#SBATCH --partition=slurm
#SBATCH --job-name=e3sm_dart_post6h
#SBATCH --nodes=1
#SBATCH --output=e3sm_dart_post6h.%j
#SBATCH --exclusive
#SBATCH --no-kill
#SBATCH --requeue

# Step 9: 6-hourly post-processing for DA cycling runs (replaces Steps 7 and 8
# for runs set up with all history streams at nhtfrq=-6).
#
# Every DA cycle is a separate 6-h model run, so the model's own daily and
# monthly history averages only ever cover the last cycle. This step builds
# all products from one 6-h record per cycle and derives daily and monthly
# values from them:
#   history (per member): post6h/<comp>/180x360_aave/{6hourly,daily,monthly}/
#   DART mean/spread on pressure levels: ${my_dart_root}/post6h/180x360_plev/
# A history record is accepted only if its time bounds span exactly one cycle,
# and every month must have all of its cycles; anything else stops the step.

set -o pipefail

# User settings
POST_MONTHS="2011-12"                    # whole calendar months, space separated
# component:stream[:rule]; rule is mean6h (default, 6-h averages) or instant
# (one snapshot per cycle at the analysis time, e.g. the ELM-DART input h1).
HIST_STREAMS="atm:eam.h0 atm:eam.h1 atm:eam.h2 lnd:elm.h0 lnd:elm.h1 rof:mosart.h1"
RUN_DART_STATE="TRUE"                    # DART forecast and analysis mean/spread
MAX_CONCURRENT_WORKERS=8
OVERWRITE_EXISTING="FALSE"

fail() { echo "ERROR: $*" >&2; exit 1; }

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  WORK_DIR=$(readlink -f "${SLURM_SUBMIT_DIR}") || fail "cannot resolve SLURM_SUBMIT_DIR"
else
  WORK_DIR=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")") || fail "cannot resolve script path"
fi
cd "${WORK_DIR}"
[[ -r create_and_setup_case.sh ]] || fail "missing create_and_setup_case.sh"
source ./create_and_setup_case.sh
source "${my_workflow_lib}/post6h/post6h_common.sh"
mkdir -p "${my_log_dir}" "${my_status_dir}" "${my_lock_dir}"
exec 8>"${my_lock_dir}/post6h.lock"
flock -n 8 || fail "another Step 9 driver is active"

[[ "${OVERWRITE_EXISTING}" =~ ^(TRUE|FALSE)$ ]] || fail "OVERWRITE_EXISTING must be TRUE or FALSE"
[[ "${RUN_DART_STATE}" =~ ^(TRUE|FALSE)$ ]] || fail "RUN_DART_STATE must be TRUE or FALSE"
[[ "${MAX_CONCURRENT_WORKERS}" =~ ^[1-9][0-9]*$ ]] || fail "MAX_CONCURRENT_WORKERS must be positive"
(( 24 % my_e3sm_cycle_hours == 0 )) || fail "cycle length ${my_e3sm_cycle_hours} h does not divide a day"
for ym in ${POST_MONTHS}; do
  [[ "${ym}" =~ ^[0-9]{4}-(0[1-9]|1[0-2])$ ]] || fail "invalid month: ${ym}"
done
members=(); for i in $(seq 1 "${my_ensnum}"); do members+=("$(printf 'EN%02d' "${i}")"); done

# Preflight: every cycle file of every requested product must exist before any work starts.
echo "Preflight: months=${POST_MONTHS} streams=${HIST_STREAMS} dart=${RUN_DART_STATE} members=${my_ensnum}"
missing=0
for ym in ${POST_MONTHS}; do
  mapfile -t stamps < <(post6h_month_stamps "${ym}")
  for st in "${stamps[@]}"; do
    for spec in ${HIST_STREAMS}; do
      for m in "${members[@]}"; do
        IFS=: read -r comp stream _ <<<"${spec}"
        f="${my_modeldir}/${m}/archive/${comp}/hist/${my_casename}.${m}.${stream}.${st}.nc"
        [[ -s "${f}" ]] || { (( missing < 20 )) && echo "  missing: ${f}" >&2; missing=$(( missing + 1 )); }
      done
    done
    if [[ "${RUN_DART_STATE}" == "TRUE" ]]; then
      for p in dart.e.eam_forecast_mean dart.e.eam_forecast_sd dart.i.eam_output_mean dart.i.eam_output_sd; do
        f="${my_eam_dart_run_dir}/${st}/${my_casename}.${p}.${st}.nc"
        [[ -s "${f}" ]] || { (( missing < 20 )) && echo "  missing: ${f}" >&2; missing=$(( missing + 1 )); }
      done
    fi
  done
done
(( missing == 0 )) || fail "preflight found ${missing} missing input file(s)"
echo "Preflight passed"

# Run workers with a concurrency limit; each writes its own log.
pids=(); labels=(); failed=0
reap() {
  local k
  for k in "${!pids[@]}"; do
    wait "${pids[$k]}" || { echo "ERROR: ${labels[$k]} failed (see log)" >&2; failed=1; }
  done
  pids=(); labels=()
}
launch() { # label log command...
  local label=$1 log=$2; shift 2
  echo "Running ${label}"
  ( "$@" ) > "${log}" 2>&1 &
  pids+=("$!"); labels+=("${label}")
  (( ${#pids[@]} < MAX_CONCURRENT_WORKERS )) || reap
}
export POSTPROCESS_DRIVER_ACTIVE=TRUE
jid=${SLURM_JOB_ID:-$$}

for ym in ${POST_MONTHS}; do
  if [[ "${RUN_DART_STATE}" == "TRUE" ]]; then
    for kind in forecast analysis; do
      rec="${my_status_dir}/post6h_complete.dart_${kind}.${ym}"
      if [[ "${OVERWRITE_EXISTING}" == "FALSE" && -s "${rec}" ]]; then echo "Skipping verified DART ${kind} ${ym}"; continue; fi
      rm -f -- "${rec}"
      launch "DART ${kind} ${ym}" "${my_log_dir}/post6h.dart_${kind}.${ym}.${jid}.log" \
        env POST6H_MONTH="${ym}" POST6H_KIND="${kind}" bash "${my_workflow_lib}/post6h/post6h_dart_state.sh"
    done
  fi
  for spec in ${HIST_STREAMS}; do
    IFS=: read -r comp stream rule <<<"${spec}"
    for m in "${members[@]}"; do
      rec="${my_status_dir}/post6h_complete.${stream}.${m}.${ym}"
      if [[ "${OVERWRITE_EXISTING}" == "FALSE" && -s "${rec}" ]]; then echo "Skipping verified ${stream} ${m} ${ym}"; continue; fi
      rm -f -- "${rec}"
      launch "${stream} ${m} ${ym}" "${my_log_dir}/post6h.${stream}.${m}.${ym}.${jid}.log" \
        env POST6H_MONTH="${ym}" POST6H_COMP="${comp}" POST6H_STREAM="${stream}" POST6H_RULE="${rule:-mean6h}" \
        POST6H_MEMBER="${m}" bash "${my_workflow_lib}/post6h/post6h_hist_stream.sh"
    done
  done
done
reap
(( failed == 0 )) || fail "one or more workers failed; logs are in ${my_log_dir}/post6h.*.${jid}.log"

# Cross-member check: every member must have produced the same fields per stream.
for ym in ${POST_MONTHS}; do
  for spec in ${HIST_STREAMS}; do
    IFS=: read -r comp stream _ <<<"${spec}"
    ref=""
    for m in "${members[@]}"; do
      rec="${my_status_dir}/post6h_complete.${stream}.${m}.${ym}"
      [[ -s "${rec}" ]] || fail "missing completion record: ${rec}"
      cur=$(awk -F= '$1=="fields"{print $2}' "${rec}")
      [[ -z "${ref}" ]] && ref="${cur}"
      [[ "${cur}" == "${ref}" ]] || fail "${stream} ${ym}: ${m} has a different field list than ${members[0]}"
    done
    echo "Verified ${stream} ${ym}: ${#members[@]} members, $(awk -F, '{print NF}' <<<"${ref}") fields each"
  done
done
echo "Step 9 completed: months=${POST_MONTHS}"
