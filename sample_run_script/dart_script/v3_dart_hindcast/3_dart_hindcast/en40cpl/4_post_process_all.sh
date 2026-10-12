#!/bin/bash
#------------------------------------------------------------------------------
# Step 4 -- submit the post-processing jobs (run on a login node). One batch job per member and per post script in
# workflow_lib/post/ (post_process_eam_daily.sh, ..._eam_6hourly.sh,
# ..._eam_monthly.sh, ..._elm_daily.sh, ..._elm_monthly.sh).
#
#   ./4_post_process_all.sh [members] [products...]
#   ./4_post_process_all.sh all eam_daily elm_daily
#
#   products: eam_6hourly eam_daily eam_clim elm_daily elm_clim (default: all)
#
# Before submitting anything it checks, for every selected member, the
# hindcast_complete record (Step 3) and every history file the products need.
# Each job writes runtmp/status/postprocess_complete.<product>.ENxx.<range>;
# verified products and products already queued are not submitted again.
# After these jobs finish: sbatch 5_post_6hr_to_daily.sh
#------------------------------------------------------------------------------

# User settings (may also be given in the environment)
MEMBER_SELECTION="${1:-${MEMBER_SELECTION:-all}}"   # all | 1,3,7-9 | EN05
POST_START="${POST_START:-}"                        # YYYY-MM-DD; empty -> my_post_start
POST_END="${POST_END:-}"                            # YYYY-MM-DD, exclusive; empty -> my_post_end
OVERWRITE_EXISTING="${OVERWRITE_EXISTING:-FALSE}"   # redo products that already completed
REQUIRE_UPSTREAM_RECORDS="${REQUIRE_UPSTREAM_RECORDS:-TRUE}"   # FALSE only for legacy archives

set -Eeo pipefail
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  WORK_DIR=$(readlink -f "${SLURM_SUBMIT_DIR}")
else
  WORK_DIR=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
fi
[[ -r "${WORK_DIR}/workflow_lib/common.sh" ]] || { echo "ERROR: ${WORK_DIR} is not a hindcast workflow directory; submit from it" >&2; exit 1; }
WF_SCRIPT_NAME=$(basename "${BASH_SOURCE[0]}")
source "${WORK_DIR}/workflow_lib/common.sh"
trap wf_exit_summary EXIT
wf_init "${WORK_DIR}"
source "${my_workflow_lib}/post/post_lib.sh"

shift || true
PRODUCTS=("$@")
(( ${#PRODUCTS[@]} )) || PRODUCTS=($(post_products))
for p in "${PRODUCTS[@]}"; do
  product_info "$p" >/dev/null || die "unknown product '${p}' (have: $(post_products))"
  [[ $p != eam_daily_from_6hourly ]] || die "eam_daily_from_6hourly is run by 5_post_6hr_to_daily.sh"
done
for v in OVERWRITE_EXISTING REQUIRE_UPSTREAM_RECORDS; do printf -v "$v" '%s' "$(normalize_bool "$v" "${!v}")"; done
POST_START=${POST_START:-$my_post_start}
POST_END=${POST_END:-$my_post_end}
for d in "$POST_START" "$POST_END"; do
  [[ $d =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]] && date -d "$d" >/dev/null 2>&1 || die "invalid post date: $d"
done
[[ $POST_START < $POST_END ]] || die "POST_START must be before POST_END"
select_members "$MEMBER_SELECTION"
summary_header
log "post range ${POST_START} .. ${POST_END} (exclusive): ${PRODUCTS[*]}"
stage_lock step4_post_submit
require_own_case
for map in "$my_e3sm_rgdmap" "$my_elm_rgdmap"; do [[ -f $map ]] || die "map file not found: ${map}"; done
command -v ncdump >/dev/null || die "ncdump not on PATH"

#--- preflight: upstream records and every needed history file ---------------------
problems=()
active_jobs >/dev/null
for m in "${MEMBERS[@]}"; do
  t=$(mtag "$m")
  if [[ $REQUIRE_UPSTREAM_RECORDS == TRUE ]] && ! record_valid "hindcast_complete.${t}"; then
    problems+=("${t}: no hindcast_complete record (run ./3_check_hindcast.sh)"); continue
  fi
  case $(member_status "$m") in PENDING|RUNNING) problems+=("${t}: model job still active"); continue ;; esac
  for p in "${PRODUCTS[@]}"; do
    line=$(check_post_inputs "$p" "$m") || problems+=("$line")
  done
done
if (( ${#problems[@]} )); then
  for p in "${problems[@]}"; do printf '  PROBLEM: %s\n' "$p" >&2; done
  die "post-processing preflight failed with ${#problems[@]} problem(s); nothing submitted"
fi
log "preflight passed"

#--- submit one job per member and product -----------------------------------------
export POST_START POST_END
n=0 skipped=0
for m in "${MEMBERS[@]}"; do
  t=$(mtag "$m")
  for p in "${PRODUCTS[@]}"; do
    jobname="hc.${p}.$(mcase "$m")"
    if [[ $OVERWRITE_EXISTING == FALSE ]] && record_valid "$(post_record "$p" "$m")" && validate_post_outputs "$p" "$m"; then
      skipped=$((skipped+1)); continue
    fi
    if active_jobs | awk -v n="$jobname" '$2 == n { found=1 } END { exit !found }'; then
      log "${t} ${p}: already queued/running, not submitted again"; continue
    fi
    worker=${my_workflow_lib}/post/$(product_info "$p" | awk '{print $1}')
    id=$(sbatch --parsable -J "$jobname" -o "${LOG_DIR}/post_${p}_${t}.%j" \
           --export=ALL,POSTPROCESS_DRIVER_ACTIVE=TRUE,WORKFLOW_ROOT="${my_workflow_root}",POST_MEMBER="${m}",POST_START="${POST_START}",POST_END="${POST_END}",FORCE="$([[ $OVERWRITE_EXISTING == TRUE ]] && echo 1)" \
           "$worker") || die "sbatch failed for ${p} ${t}"
    log "${t} ${p}: submitted job ${id}"
    n=$((n+1))
  done
done
log "submitted ${n} job(s), ${skipped} already complete; logs: runtmp/logs/post_<product>_ENxx.<jobid>"
log "when they have finished: sbatch 5_post_6hr_to_daily.sh"
