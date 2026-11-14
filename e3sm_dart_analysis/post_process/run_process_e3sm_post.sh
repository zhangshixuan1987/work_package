#!/bin/bash
# Make one zppy-style post product of one F20TR_ne30pg2_r05_IcoswISC30E3r5_* ensemble member, the same
# way as for the existing configurations (CLIM_BASE, HYB_BASE, DA_ATM, DA_LND, DA_COUPLED).
#
# e3sm_zppy/<kind>.bash are those configurations' own zppy scripts (taken from CLIM_BASE_2011121600/EN01)
# with the run-specific header values replaced by @ROOT@ @CASE@ @ENS@ @YR_SRT@ @YR_END@ @MTH_SRT@ @MTH_END@,
# and MPI parallelism switched to serial (--parallel=serial --jobs=1, ncremap -p serial) so that many
# products run side by side on one node. The daily products (ts_*_daily_*, ts_atm_6hourly_*) make
# full-length series instead of zppy's truncated `ncclimo --split` ones (see e3sm_zppy/README.md);
# the others are unchanged. This script fills in the header, writes the script to
# <member>/post/scripts/<kind>_<member>_<tag>.bash (the name zppy used) and runs it; the template
# writes post/scripts/<kind>_<tag>.status itself (and regen_ts_<stream>.status for 180x360 daily).
#
# Usage: run_process_e3sm_post.sh <case dir> <member> <kind>
#   kind: a template name, e.g. ts_atm_monthly_180x360_aave, ts_atm_daily_180x360_aave
#   THREADS (default 4): ncremap threads of the daily 180x360 products
# Skips the product if its status is already OK, unless FORCE=1. For the daily products, only a
# full-length status (OK files=... records=...) counts: zppy's plain "OK" was for truncated files.
#
# Periods (from the start month S of <case>_<YYYYMMDDHH>, as in the existing configurations):
#   monthly ts and climatologies  S .. S+1
#   lnd daily                     S .. S
#   atm daily/6-hourly            S .. S+2
#   rof daily                     month of the first mosart.h1 file .. S+2
#                                 (S+1 for the starts on the 31st, whose first daily river file is the 1st)

set -u
case_dir=${1:?case directory}
ens=${2:?member, e.g. EN01}
kind=${3:?product, e.g. ts_atm_monthly_180x360_aave}

here=$(cd "$(dirname "$(readlink -f "$0")")" && pwd)
tmpl=$here/e3sm_zppy/$kind.bash
[[ -f $tmpl ]] || { echo "unknown product $kind (no $tmpl)" >&2; exit 1; }

root=$(cd "$(dirname "$case_dir")" && pwd)
case=$(basename "$case_dir")
start=${case##*_}
[[ $start =~ ^[0-9]{10}$ ]] || { echo "cannot read start date from $case" >&2; exit 1; }

# month index helpers (months since year 0)
mi() { echo $(( 10#$1 * 12 + 10#$2 - 1 )); }
ym() { printf '%04d %02d' $(( $1 / 12 )) $(( $1 % 12 + 1 )); }
s=$(mi "${start:0:4}" "${start:4:2}")

case $kind in
  ts_lnd_daily_*)            i0=$s; i1=$s ;;
  ts_atm_daily_*|ts_atm_6hourly_*) i0=$s; i1=$((s + 2)) ;;
  ts_rof_daily_*)
    first=$(ls "$case_dir/$ens/archive/rof/hist/$case.$ens.mosart.h1."*.nc 2>/dev/null | head -1)
    [[ -n $first ]] || { echo "no mosart.h1 files for $case $ens" >&2; exit 1; }
    d=$(basename "$first" | grep -oE '\.[0-9]{4}-[0-9]{2}-[0-9]{2}-' | tr -d '.-')
    i0=$(mi "${d:0:4}" "${d:4:2}"); i1=$((s + 2)) ;;
  *)                         i0=$s; i1=$((s + 1)) ;;   # monthly ts, climatologies
esac
read -r y0 m0 <<< "$(ym "$i0")"
read -r y1 m1 <<< "$(ym "$i1")"
tag="${y0}${m0}-${y1}${m1}-0000"

scr_dir=$case_dir/$ens/post/scripts
status=$scr_dir/${kind}_${tag}.status

# DRYRUN=1 only writes the filled-in script, to DRYRUN_DIR (default: a scratch dir), never to the case
if [[ "${DRYRUN:-0}" == 1 ]]; then
  out=${DRYRUN_DIR:-${SCRATCH:-/tmp}/e3sm_zppy_dryrun}/$case/$ens; mkdir -p "$out"
  script=$out/${kind}_${ens}_${tag}.bash
else
  mkdir -p "$scr_dir"
  done_re='^OK'
  [[ $kind == ts_*_daily_* || $kind == ts_atm_6hourly_* ]] && done_re='^OK files=[0-9]+ records='
  if [[ "${FORCE:-0}" != 1 && -f $status ]] && grep -qE "$done_re" "$status"; then
    echo "[SKIP] $case $ens $kind"; exit 0
  fi
  script=$scr_dir/${kind}_${ens}_${tag}.bash
fi
sed -e "s#@ROOT@#$root#g" -e "s#@CASE@#$case#g" -e "s#@ENS@#$ens#g" \
    -e "s#@YR_SRT@#$y0#g" -e "s#@YR_END@#$y1#g" -e "s#@MTH_SRT@#$m0#g" -e "s#@MTH_END@#$m1#g" \
    "$tmpl" > "$script"

[[ "${DRYRUN:-0}" == 1 ]] && { echo "[DRYRUN] $script"; exit 0; }

# ts/clim scripts expect SLURM_JOBID for the status/workdir names
export SLURM_JOBID=${SLURM_JOBID:-local$$}
bash "$script" > "$scr_dir/${kind}_${ens}_${tag}.log" 2>&1
rc=$?
if (( rc == 0 )) && grep -q '^OK' "$status"; then
  echo "[OK] $case $ens $kind $tag"
else
  echo "[FAIL] $case $ens $kind $tag (rc=$rc, $(cat "$status" 2>/dev/null)); log: $scr_dir/${kind}_${ens}_${tag}.log" >&2
  exit 1
fi
