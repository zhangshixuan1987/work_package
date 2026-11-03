#!/bin/bash
#SBATCH --job-name=regen_f20tr_ts
#SBATCH --account=e3sm
#SBATCH --qos=regular
#SBATCH --constraint=cpu
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --output=regen_f20tr_ts.o%j
#
# Regenerate the full-length 180x360 time series of every F20TR_ne30pg2_r05_IcoswISC30E3r5_* member
# (5 configurations x 8 start dates x 20 members) for the atm daily/6-hourly and lnd/rof daily
# streams with regen_f20tr_ts.sh. Tasks are spread over the nodes (one srun task per node) and run
# PARALLEL at a time per node; finished streams (status OK) are skipped, so the job can be resubmitted.
#
#   sbatch submit_regen_f20tr_ts.sbatch                  # all cases, all streams
#   CASES="DA_COUPLED_2011121600" STREAMS="atm_daily" sbatch submit_regen_f20tr_ts.sbatch
#   OUT_ROOT=/pscratch/sd/z/zhan391/f20tr_post sbatch submit_regen_f20tr_ts.sbatch

ROOT=${ROOT:-/global/cfs/cdirs/m4849/zhan391/eam_hindcast}
PREFIX=F20TR_ne30pg2_r05_IcoswISC30E3r5
STREAMS=${STREAMS:-"atm_6hourly atm_daily lnd_daily rof_daily"}   # longest first
PARALLEL=${PARALLEL:-24}
THREADS=${THREADS:-4}
SCRIPT=${SCRIPT:-$SLURM_SUBMIT_DIR/regen_f20tr_ts.sh}
export OUT_ROOT=${OUT_ROOT:-}

source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

tasks=$SLURM_SUBMIT_DIR/regen_f20tr_ts.tasks.$SLURM_JOB_ID
: > "$tasks"
if [[ -n "${CASES:-}" ]]; then
  case_dirs=(); for c in $CASES; do case_dirs+=("$ROOT/${PREFIX}_$c"); done
else
  case_dirs=("$ROOT/${PREFIX}"_*)
fi
for s in $STREAMS; do
  for c in "${case_dirs[@]}"; do
    for m in "$c"/EN[0-9][0-9]; do echo "$c $(basename "$m") $s $THREADS" >> "$tasks"; done
  done
done
echo "$(wc -l < "$tasks") tasks on $SLURM_NNODES nodes, $PARALLEL per node"

# node k takes every SLURM_NNODES-th task starting at k
srun --ntasks="$SLURM_NNODES" --ntasks-per-node=1 bash -c '
  awk -v n="$SLURM_NTASKS" -v k="$SLURM_PROCID" "(NR-1)%n==k" "'"$tasks"'" |
    xargs -P '"$PARALLEL"' -L 1 bash '"$SCRIPT"'
'

grep -L "^OK" $(for c in "${case_dirs[@]}"; do d=${OUT_ROOT:+$OUT_ROOT/$(basename "$c")}; d=${d:-$c}
  for s in $STREAMS; do echo "$d"/EN*/post/scripts/regen_ts_${s}.status; done; done) 2>/dev/null \
  | sed 's/^/[NOT OK] /'
rm -f "$tasks"
