#!/bin/bash
#SBATCH --job-name=process_e3sm_post
#SBATCH --account=e3sm
#SBATCH --qos=regular
#SBATCH --constraint=cpu
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --output=process_e3sm_post.o%j
#
# Make the post/ tree of F20TR_ne30pg2_r05_IcoswISC30E3r5_* members, the same products as for the
# existing configurations (CLIM_BASE, HYB_BASE, DA_ATM, DA_LND, DA_COUPLED), with run_process_e3sm_post.sh
# and the templates in e3sm_zppy/:
#   - full-length atm daily/6-hourly and lnd/rof daily series (180x360 and native), read by the
#     fcst_init_eval / hind_lead_eval / hind_ic_test notebooks (regen_ts_<stream>.status);
#   - monthly atm (19 pressure levels), global-mean, lnd and rof series (180x360 and native).
# Climatologies (e3sm_zppy/clim_*) are not made: ncclimo needs whole years of monthly files, so they
# fail for these 2-3 month hindcasts (as they did for every existing configuration).
# Tasks (one product of one member) are spread over the nodes and run PARALLEL at a time per node;
# finished products (status OK) are skipped, so the job can be resubmitted.
#
#   sbatch submit_process_e3sm_post.sbatch                                  # all cases, all products
#   CASES="DA_LND_CLIMFILL_2011121600 ..." sbatch submit_process_e3sm_post.sbatch
#   CASES="DA_COUPLED_2011121600" KINDS="ts_atm_daily_180x360_aave" sbatch submit_process_e3sm_post.sbatch
#   ROOT=/pscratch/sd/z/zhan391/e3sm_dart/simulation CASES="..." sbatch ...   # cases on scratch

ROOT=${ROOT:-/global/cfs/cdirs/m4849/zhan391/eam_hindcast}
PREFIX=F20TR_ne30pg2_r05_IcoswISC30E3r5
KINDS=${KINDS:-"ts_atm_6hourly_180x360_aave ts_atm_daily_180x360_aave ts_lnd_daily_180x360_aave
  ts_rof_daily_180x360_aave ts_lnd_daily_native ts_rof_daily_native ts_atm_monthly_180x360_aave
  ts_atm_monthly_glb ts_lnd_monthly_180x360_aave ts_lnd_monthly_native ts_rof_monthly_180x360_aave
  ts_rof_monthly_native"}   # longest first
PARALLEL=${PARALLEL:-24}
export THREADS=${THREADS:-4}
SCRIPT=${SCRIPT:-$SLURM_SUBMIT_DIR/run_process_e3sm_post.sh}

source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

tasks=$SLURM_SUBMIT_DIR/process_e3sm_post.tasks.$SLURM_JOB_ID
: > "$tasks"
if [[ -n "${CASES:-}" ]]; then
  case_dirs=(); for c in $CASES; do case_dirs+=("$ROOT/${PREFIX}_$c"); done
else
  case_dirs=("$ROOT/${PREFIX}"_*)
fi
for k in $KINDS; do
  for c in "${case_dirs[@]}"; do
    [[ -d $c ]] || { echo "[MISSING] $c"; continue; }
    for m in "$c"/EN[0-9][0-9]; do [[ -d $m ]] && echo "$c $(basename "$m") $k" >> "$tasks"; done
  done
done
echo "$(wc -l < "$tasks") tasks on $SLURM_NNODES nodes, $PARALLEL per node"

# node k takes every SLURM_NNODES-th task starting at k
srun --ntasks="$SLURM_NNODES" --ntasks-per-node=1 bash -c '
  awk -v n="$SLURM_NTASKS" -v k="$SLURM_PROCID" "(NR-1)%n==k" "'"$tasks"'" |
    xargs -P '"$PARALLEL"' -L 1 bash '"$SCRIPT"'
'

# products whose status is not OK (daily products need the full-length status, not zppy's plain OK)
nbad=0
for c in "${case_dirs[@]}"; do
  [[ -d $c ]] || continue
  for m in "$c"/EN[0-9][0-9]; do
    for k in $KINDS; do
      re='^OK'; [[ $k == ts_*_daily_* || $k == ts_atm_6hourly_* ]] && re='^OK files=[0-9]+ records='
      st=("$m"/post/scripts/${k}_[0-9]*.status)
      [[ -f ${st[0]} ]] && grep -qE "$re" "${st[0]}" && continue
      echo "[NOT OK] $(basename "$c") $(basename "$m") $k"; nbad=$((nbad + 1))
    done
  done
done
echo "done: $nbad products not OK"
rm -f "$tasks"
