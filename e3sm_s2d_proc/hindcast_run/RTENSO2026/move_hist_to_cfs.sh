#!/bin/bash
#SBATCH -q xfer
#SBATCH -t 2-00:00:00
#SBATCH -J rtenso_mv_hist
#SBATCH -o /global/cfs/cdirs/e3sm/zhan391/s2d_diag/post_process_logs/move_hist_to_cfs.%x.%j.out
#
# Move the RTENSO2026 history output from pscratch to the shared CFS space, as for the earlier start dates:
#   <SRC_ROOT>/<case>/EN*/archive/{atm,cpl,ice,lnd,ocn,rof}/hist  ->  <DST_ROOT>/<case>/EN*/archive/<comp>/hist
# Restarts (archive/rest), logs, run/, case_scripts/ and build/ stay on pscratch.
#
#   sbatch move_hist_to_cfs.sh <case_dir_name> [<case_dir_name> ...]
# e.g. sbatch move_hist_to_cfs.sh WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_RTENSO2026_obs_2026110100
#
# Per case:
#   1. check the e3sm CFS project space (need the hist size + MARGIN_GB)
#   2. per member and component: rsync the hist directory (group e3sm, world-readable, like the other dates)
#   3. verify with a second rsync dry-run (size + mtime) and equal file counts and byte totals
#   4. only then delete the pscratch hist files of that member and component
# Any failure stops the job; the pscratch files of the failed and later directories are left untouched.

set -u
SRC_ROOT=/pscratch/sd/z/zhan391/e3sm_project/E3SMv3_S2D
DST_ROOT=/global/cfs/cdirs/e3sm/S2S2D/RTENSO2026/simulation
COMPS=(atm cpl ice lnd ocn rof)
MARGIN_GB=2000
RSYNC_OPTS=(-rlt --perms --chmod=Du=rwx,Dgo=rx,Fu=rw,Fgo=r)

bytes() { find "$1" -type f -printf '%s\n' | awk '{s+=$1} END{print s+0}'; }

for name in "$@"; do
  src=$SRC_ROOT/$name; dst=$DST_ROOT/$name
  echo "================ $name  $(date)"
  [[ -d $src ]] || { echo "ERROR: missing source $src"; exit 1; }

  need=$(du -sc --block-size=1G "$src"/EN*/archive/*/hist 2>/dev/null | tail -1 | cut -f1)
  free=$(cfsquota e3sm | awk '$1=="e3sm"{print $3-$2}')
  echo "need=${need}G free=${free}G"
  (( free > need + MARGIN_GB )) || { echo "ERROR: not enough CFS space"; exit 1; }

  for mem in "$src"/EN*; do
    en=$(basename "$mem")
    for c in "${COMPS[@]}"; do
      s=$mem/archive/$c/hist; d=$dst/$en/archive/$c/hist
      [[ -d $s ]] || { echo "ERROR: missing $s"; exit 1; }
      ns=$(find "$s" -type f | wc -l)
      if (( ns == 0 )); then echo "$en/$c: no files on pscratch (already moved?)"; continue; fi
      mkdir -p "$d" || { echo "ERROR: mkdir $d"; exit 1; }
      rsync "${RSYNC_OPTS[@]}" "$s/" "$d/" || { echo "ERROR: rsync $en/$c failed"; exit 1; }
      diffs=$(rsync -rlt --dry-run --itemize-changes "$s/" "$d/" | grep -v '^\.d' | wc -l)
      nd=$(find "$d" -type f | wc -l); bs=$(bytes "$s"); bd=$(bytes "$d")
      if [[ $diffs -ne 0 || $ns -ne $nd || $bs -ne $bd ]]; then
        echo "ERROR: verify $en/$c: diffs=$diffs files=$ns/$nd bytes=$bs/$bd; pscratch kept"; exit 1
      fi
      find "$s" -type f -delete
      echo "$en/$c: moved $ns files ($((bs / 1024**3)) GB)"
    done
  done
  chgrp -R e3sm "$dst" && chmod g-w,o+rx "$dst" "$dst"/EN* "$dst"/EN*/archive "$dst"/EN*/archive/*
  echo "done $name $(date)"
done
echo "ALL DONE $(date)"
