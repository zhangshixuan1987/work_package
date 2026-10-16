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
out_dir="/global/cfs/cdirs/e3smdata/simulations/S2S2D"
work_dir="/global/cfs/cdirs/e3smdata/simulations/S2S2D/s2d_diag"

exp="WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_JRA55_FOSIRL"
nens=10
leads=(1980050100 1980110100 1981050100 1981110100)

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

for lead in "${leads[@]}"; do
  case_name="${exp}_${lead}"

  for (( ii=0; ii<nens; ii++ )); do
    ENSTR=$(printf "EN%02d" "${ii}")

    for comp_name in atm cpl ice lnd ocn rof; do
      indir="${run_dir}/${case_name}/${ENSTR}/archive/${comp_name}/hist"
      dstdir="${out_dir}/${case_name}/${ENSTR}/archive/${comp_name}/hist"

      if [[ ! -d "${indir}" ]]; then
        log_msg "Missing source dir: ${indir}"
        echo "Missing source dir: ${indir}" >> "${fail_log}"
        continue
      fi

      mkdir -p "${dstdir}"

      log_msg "Sync ${indir} -> ${dstdir}"
      if ! rsync -a "${indir}/" "${dstdir}/" >> "${main_log}" 2>&1; then
        log_msg "rsync failed: ${indir} -> ${dstdir}"
        echo "rsync failed: ${indir} -> ${dstdir}" >> "${fail_log}"
        continue
      fi

      find "${indir}" -maxdepth 1 -type f -name "*.nc" -exec bash -c '
        dstdir="$1"
        shift
        for f in "$@"; do
          dst="${dstdir}/$(basename "$f")"
          if [[ -f "$dst" ]] && cmp -s "$f" "$dst"; then
            rm -f "$f"
          fi
        done
      ' bash "${dstdir}" {} +
    done
  done
done

if [[ -s "${fail_log}" ]]; then
  log_msg "Completed with failures. See ${fail_log}"
  exit 1
else
  log_msg "Completed successfully."
  rm -f "${fail_log}"
fi

exit 0
