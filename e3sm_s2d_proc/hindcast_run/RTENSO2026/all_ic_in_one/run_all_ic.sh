#!/bin/bash -l
#------------------------------------------------------------------------------
# Batch system directives (run phase: 3 sets x 10 members x 9 nodes)
#------------------------------------------------------------------------------
#SBATCH  --account=e3sm
#SBATCH  --constraint=cpu
#SBATCH  --qos=regular
#SBATCH  --nodes=270
#SBATCH  --time=16:00:00
#SBATCH  --job-name=rtenso_20261101_all_ic
#SBATCH  --output=slurm-%j.out

# ==============================================================================
# The 1 Nov 2026 RTENSO2026 hindcasts of all three initial-condition sets in
# one job:
#
#   clim/           land, river, coupler and sea ice of the clim year
#   obs/            land, river, coupler and sea ice of the obs year
#   obs_ice_nsidc/  obs/ with the NSIDCv6 sea ice (ice/nsidc/)
#
# Each sub-directory is the matching ../<set>_2026/ run set with
# my_leadymd=("2026-11-01") and my_job_nnodes=90; its run_e3sm_ens_bundle.sh is
# unchanged. This script runs the three bundle drivers side by side: in the
# run phase they share one 270-node allocation (each member is an
# srun --exclusive step of 9 nodes, at most 10 per set).
#
#   ./run_all_ic.sh setup   create, build and set up the 30 cases (login node;
#                           the three sets one after the other, logs in <set>/setup.log)
#   sbatch run_all_ic.sh    short spin-up + 2-year run of all 30 members
#                           (logs in <set>/run.log.<jobid>)
#
# Case names are those of the run sets, e.g.
#   WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_RTENSO2026_obs_ice_nsidc_2026110100.EN00
# ==============================================================================

sets=(clim obs obs_ice_nsidc)
phase=${1:-run}

if [[ -n ${SLURM_SUBMIT_DIR:-} ]]; then
  top_dir=${SLURM_SUBMIT_DIR}
else
  top_dir=$(cd "$(dirname "$0")" && pwd)
fi
cd "${top_dir}" || exit 1

case "${phase}" in
  setup)
    for s in "${sets[@]}"; do
      echo "$(date '+%F %T') setup ${s}"
      ( cd "${s}" && ALL_IC_PHASE=setup ./run_e3sm_ens_bundle.sh > setup.log 2>&1 ) || {
        echo "ERROR: setup of ${s} failed (see ${top_dir}/${s}/setup.log)"; exit 1; }
    done
    echo "$(date '+%F %T') setup done; submit the run with: sbatch run_all_ic.sh"
    ;;
  run)
    [[ -n ${SLURM_JOB_ID:-} ]] || { echo "ERROR: submit the run phase with sbatch run_all_ic.sh"; exit 1; }
    declare -A pid_set=()
    for s in "${sets[@]}"; do
      ( cd "${s}" && ALL_IC_PHASE=run ./run_e3sm_ens_bundle.sh > "run.log.${SLURM_JOB_ID}" 2>&1 ) &
      pid_set[$!]=${s}
      echo "$(date '+%F %T') started ${s} (log ${s}/run.log.${SLURM_JOB_ID})"
      sleep 30     # stagger the drivers' first srun steps
    done
    failed=0
    for pid in "${!pid_set[@]}"; do
      if wait "${pid}"; then
        echo "$(date '+%F %T') ${pid_set[$pid]} finished"
      else
        echo "$(date '+%F %T') ERROR: ${pid_set[$pid]} failed (see ${pid_set[$pid]}/run.log.${SLURM_JOB_ID})"
        failed=$((failed + 1))
      fi
    done
    exit "${failed}"
    ;;
  *)
    echo "usage: ./run_all_ic.sh setup | sbatch run_all_ic.sh" >&2; exit 2 ;;
esac
