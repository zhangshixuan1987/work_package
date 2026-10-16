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
work_dir="/global/cfs/cdirs/e3sm/zhan391/s2d_diag"

exp="WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_4DEnVarOcn"
nens=10

max_jobs=12

#leads=(1980050100 1980110100 1981050100 1981110100)
#leads=(1982050100 1982110100 1983050100 1983110100)
#leads=(1984050100 1984110100 1985050100 1985110100)
#leads=(1986050100 1986110100 1987050100 1987110100)
#leads=(1988050100 1988110100 1989050100 1989110100)
#leads=(1990050100 1990110100 1991050100 1991110100)
#leads=(1992050100 1992110100 1993050100 1993110100)
#leads=(1994050100 1994110100 1995050100 1995110100)
#leads=(1996050100 1996110100 1997050100 1997110100)
#leads=(1998050100 1998110100 1999050100 1999110100)
#leads=(2000050100 2000110100 2001050100 2001110100)
#leads=(2002050100 2002110100 2003050100 2003110100)
#leads=(2004050100 2004110100 2005050100 2005110100)
#leads=(2006050100 2006110100 2007050100 2007110100)
#leads=(2008050100 2008110100 2009050100 2009110100)
#leads=(2010050100 2010110100 2011050100 2011110100)
leads=(2012050100 2012110100 2013050100 2013110100)
#leads=(2014050100 2014110100 2015050100 2015110100) 
#leads=(2015050100 2015110100 2016050100 2016110100)
#leads=(2017050100 2017110100 2018050100 2018110100)

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
