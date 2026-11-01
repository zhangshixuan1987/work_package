#!/bin/bash
# Run run_post_<driver>.bash on a login node instead of through sbatch.
# Drivers run one after another; for each driver the IC methods run side by
# side with MAX_JOBS members each and JOB_NBR ncclimo jobs per member. Login
# nodes cap each user at 32 CPUs and 30 GB, so keep MAX_JOBS*JOB_NBR*2 small.
#
#   nohup ./run_login_post.bash atm_mon1 atm_glb > run_login_post.log 2>&1 &
#
# RUN_DIR, IC_METHODS and LEADS (space-separated start dates) override the
# defaults, e.g. after 1_run_mov_archive.bash has moved the 2026 starts:
#   IC_METHODS=obs LEADS="2026020100 2026050100 2026080100" \
#     nohup ./run_login_post.bash > run_login_post.log 2>&1 &

drv_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export RUN_DIR="${RUN_DIR:-/global/cfs/cdirs/e3sm/S2S2D/RTENSO2026/simulation}"   # 1_run_mov_archive.bash out_dir
export MAX_JOBS="${MAX_JOBS:-1}"
export JOB_NBR="${JOB_NBR:-8}"
ic_methods=(${IC_METHODS:-clim obs})
log_dir="${LOG_DIR:-$PSCRATCH/rtenso2026_post_logs}"

drivers=("$@")
if (( ${#drivers[@]} == 0 )); then
  drivers=(atm_mon1 atm_glb ocn_mon1 ice_mon1 lnd_mon1 rof_mon1
           lnd_native_mon1 rof_native_mon1 ocn_rgn)
fi

mkdir -p "${log_dir}" || exit 1
cd "${log_dir}" || exit 1

for drv in "${drivers[@]}"; do
  script="${drv_dir}/run_post_${drv}.bash"
  if [[ ! -f "${script}" ]]; then
    echo "ERROR: missing driver: ${script}" >&2
    exit 1
  fi

  pids=()
  for ic in "${ic_methods[@]}"; do
    log="${log_dir}/${ic}_${drv}.login.log"
    echo "$(date '+%F %T') start ${ic} ${drv} -> ${log}"
    IC_METHOD="${ic}" bash "${script}" > "${log}" 2>&1 &
    pids+=($!)
  done

  for pid in "${pids[@]}"; do
    wait "${pid}" || echo "$(date '+%F %T') WARNING: ${drv} pid ${pid} exited non-zero"
  done
  echo "$(date '+%F %T') done ${drv}"
done
