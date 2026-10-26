#!/bin/bash -el
#SBATCH --job-name=post_process
#SBATCH --account=e3sm
#SBATCH --nodes=1
#SBATCH --output=post_process.o%j
#SBATCH --exclusive
#SBATCH --time=24:00:00
#SBATCH --qos=regular
#SBATCH --constraint=cpu

source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

run_dir="/pscratch/sd/z/zhan391/e3sm_project/E3SMv3_S2D"
out_dir="/global/cfs/cdirs/e3sm/S2S2D/simulation" 
work_dir="/global/cfs/cdirs/e3sm/S2S2D/s2d_diag"

ic_method="${IC_METHOD:-clim}"   # clim | obs | hybrid_ic
exp="WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_RTENSO2026_${ic_method}"
nens=10
max_jobs=12
leads=(2015020100 2015050100 2015080100 2015110100)

log_dir="${work_dir}/post_process_logs"
mkdir -p "${log_dir}"

main_log="${log_dir}/post_process_${SLURM_JOB_ID:-manual}.log"
fail_log="${log_dir}/post_process_${SLURM_JOB_ID:-manual}.failed"

: > "${main_log}"
: > "${fail_log}"

timestamp() {
  date "+%Y-%m-%d %H:%M:%S"
}

log_msg() {
  echo "[$(timestamp)] $*" | tee -a "${main_log}"
}

# Clean shutdown on cancel
trap 'log_msg "Received termination signal. Killing child jobs..."; jobs -rp | xargs -r kill || true; wait || true' SIGINT SIGTERM

move_hist_files() {
  local indir="$1"
  local dstdir="$2"

  if [[ ! -d "${indir}" ]]; then
    echo "[$(timestamp)] SKIP missing source dir: ${indir}" >> "${main_log}"
    return 0
  fi

  mkdir -p "${dstdir}"

  local src_real dst_real
  src_real=$(realpath "${indir}")
  dst_real=$(realpath "${dstdir}")

  if [[ "${src_real}" == "${dst_real}" ]]; then
    echo "[$(timestamp)] SKIP same source/destination: ${indir}" >> "${main_log}"
    return 0
  fi

  shopt -s nullglob
  local files=( "${indir}"/*.nc )
  shopt -u nullglob

  if (( ${#files[@]} == 0 )); then
    echo "[$(timestamp)] SKIP no .nc files: ${indir}" >> "${main_log}"
    return 0
  fi

  echo "[$(timestamp)] START ${indir} -> ${dstdir} (${#files[@]} files)" >> "${main_log}"

  if rsync -a --partial --ignore-times --info=stats1 --remove-source-files \
      "${files[@]}" "${dstdir}/" >> "${main_log}" 2>&1; then
    echo "[$(timestamp)] DONE  ${indir} -> ${dstdir}" >> "${main_log}"
    return 0
  else
    echo "[$(timestamp)] FAIL  ${indir} -> ${dstdir}" >> "${main_log}"
    echo "${indir} -> ${dstdir}" >> "${fail_log}"
    return 1
  fi
}

wait_for_slot() {
  while (( $(jobs -rp | wc -l) >= max_jobs )); do
    wait -n || true
  done
}

log_msg "Starting post-process transfer"
log_msg "run_dir=${run_dir}"
log_msg "out_dir=${out_dir}"
log_msg "max_jobs=${max_jobs}"

for lead in "${leads[@]}"; do
  case_name="${exp}_${lead}"

  for (( ii=0; ii<nens; ii++ )); do
    ENSTR=$(printf "EN%02d" "${ii}")

    for comp_name in atm cpl ice lnd ocn rof; do
      indir="${run_dir}/${case_name}/${ENSTR}/archive/${comp_name}/hist"
      dstdir="${out_dir}/${case_name}/${ENSTR}/archive/${comp_name}/hist"

      wait_for_slot
      move_hist_files "${indir}" "${dstdir}" &
    done
  done
done

# Wait for all jobs safely
while (( $(jobs -rp | wc -l) > 0 )); do
  wait -n || true
done

if [[ -s "${fail_log}" ]]; then
  log_msg "Completed with failures. See ${fail_log}"
  exit 1
else
  log_msg "Completed successfully with no transfer failures."
  rm -f "${fail_log}"
fi

exit 0
