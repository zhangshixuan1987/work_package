#!/bin/bash
# Regenerate the 180x360 time series of one F20TR_ne30pg2_r05_IcoswISC30E3r5_* ensemble member
# for one output stream, over the whole forecast.
#
# Why: the zppy-style ts_* scripts in <case>/<member>/post/scripts call `ncclimo --split`, which for
# these 2-month runs passes only the first 4 history files of a stream to ncrcat (its yearly segment
# logic assumes runs of whole years). The atm daily/6-hourly files thus stop after 4 days and the
# land/river daily files at the end of the start month. This script concatenates every history file
# with ncrcat, remaps once with ncremap, and writes one file per variable to the same place and with
# the same names (<VAR>_<YYYYMM>_<YYYYMM>.nc) as before, replacing the truncated files.
#
# Usage: regen_f20tr_ts.sh <case dir> <member> <stream> [threads]
#   stream: atm_daily | atm_6hourly | lnd_daily | rof_daily
# Writes <member>/post/scripts/regen_ts_<stream>.status (OK when the file count and record count
# check out) and skips the member if that status is already OK, unless FORCE=1.
# OUT_ROOT=<dir> writes post/ under <dir>/<case>/<member>/ instead of next to the archive.

set -u
case_dir=${1:?case directory}
ens=${2:?member, e.g. EN01}
stream=${3:?stream}
thr=${4:-4}

MAP_ATM=/global/cfs/cdirs/e3sm/diagnostics/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc
MAP_LND=/global/cfs/cdirs/e3sm/diagnostics/maps/map_r05_to_cmip6_180x360_aave.20231110.nc

case $stream in
  atm_daily)
    comp=atm; hist=eam.h1; freq=daily; map=$MAP_ATM; prc=eam; xtr=""
    vars=TCO,SCO,FLUT,PRECT,PS,PSL,TS,TBOT,TUQ,TVQ,UBOT,VBOT,ZBOT,U10,TAUX,TAUY,TREFHT,QREFHT,TMQ,CAPE,CIN,TTQ
    for f in U V Z T Q OMEGA; do
      for p in 1000 975 950 925 900 850 800 700 600 500 400 300 200 100 050 010; do
        [[ $f == OMEGA && ( $p == 050 || $p == 010 ) ]] && continue
        vars+=",$f$p"
      done
    done
    vars+=,TREFHTMN,TREFHTMX,PRECTMX ;;
  atm_6hourly)
    comp=atm; hist=eam.h2; freq=6hourly; map=$MAP_ATM; prc=eam; xtr=""
    vars=RHREFHT,TTOP,TOZ,PS,PSL,OMEGA500,PRECT,PRECL,PRECC,PRECSL,TS,CLDLOW,CLDMED,CLDHGH,CLDTOT
    vars+=,TGCLDLWP,TGCLDIWP,TGCLDCWP,TMQ,CAPE,CIN,TREFHT,QREFHT,U10,TAUX,TAUY,SHFLX,LHFLX,TUQ,TVQ
    vars+=,FSNT,FLNT,FSNTOA,FLUT,FSUTOA,FSNTC,FLNTC,FSNTOAC,FLUTC,SWCF,LWCF,FSNS,FLNS,FSDS,FLDS
    vars+=,FSNSC,FLNSC,FSDSC,FLDSC,QFLX,U90M,V90M
    vars+=,CLDTOT_ISCCP,MEANCLDALB_ISCCP,MEANTAU_ISCCP,MEANPTOP_ISCCP,MEANTB_ISCCP
    for c in TOT HGH MED LOW; do vars+=",CLD${c}_CAL,CLD${c}_CAL_LIQ,CLD${c}_CAL_ICE,CLD${c}_CAL_UN"; done ;;
  lnd_daily)
    comp=lnd; hist=elm.h1; freq=daily; map=$MAP_LND; prc=elm; xtr=area,landfrac,landmask,topo
    vars=H2OSNO,H2OSOI,SOILWATER_10CM,TWS,ZWT,TSOI_10CM,TSA ;;
  rof_daily)
    comp=rof; hist=mosart.h1; freq=daily; map=$MAP_LND; prc=""; xtr=mask,area,areatotal,areatotal2  # ncremap has no mosart type
    vars=RIVER_DISCHARGE_OVER_LAND_LIQ ;;
  *) echo "unknown stream $stream" >&2; exit 1 ;;
esac

case_name=$(basename "$case_dir")
mem_dir=$case_dir/$ens
out_dir=${OUT_ROOT:+$OUT_ROOT/$case_name/$ens}
out_dir=${out_dir:-$mem_dir}
scr_dir=$out_dir/post/scripts
dest=$out_dir/post/$comp/180x360_aave/ts/$freq/0yr
status=$scr_dir/regen_ts_${stream}.status
mkdir -p "$scr_dir"

if [[ "${FORCE:-0}" != 1 && -f $status ]] && grep -q '^OK' "$status"; then
  echo "[SKIP] $case_name $ens $stream"; exit 0
fi
echo "RUNNING ${SLURM_JOB_ID:-local}" > "$status"
fail() { echo "ERROR $1" > "$status"; echo "[FAIL] $case_name $ens $stream: $1" >&2; exit 1; }

mapfile -t files < <(ls "$mem_dir"/archive/$comp/hist/${case_name}.${ens}.${hist}.*.nc 2>/dev/null)
(( ${#files[@]} )) || fail "no ${hist} files"

work=$(mktemp -d "${SCRATCH:-/tmp}/regen_ts.${case_name##*_}.${ens}.${stream}.XXXX") || fail "mktemp"
trap 'rm -rf "$work"' EXIT

# 1. concatenate the whole run (native grid), 2. remap once, 3. split per variable
ncrcat -O -7 -L 1 -v "${vars}${xtr:+,$xtr}" "${files[@]}" "$work/native.nc" || fail "ncrcat"
ncremap ${prc:+-P "$prc"} -t "$thr" -m "$map" -i "$work/native.nc" -o "$work/rgr.nc" > "$work/ncremap.log" 2>&1 \
  || { cat "$work/ncremap.log" >&2; fail "ncremap"; }

nrec=$(ncks --trd -M -m "$work/rgr.nc" | awk '/^Record dimension: name = time, size =/ {print $NF; exit}')
[[ -z "$nrec" ]] && nrec=$(ncks -C -v time --cdl -m "$work/rgr.nc" | sed -n 's/.*time = UNLIMITED ; \/\/ (\([0-9]*\) currently).*/\1/p' | head -1)
tag=$(ncks -H -C -v time --cal "$work/rgr.nc" | grep -oE '[0-9]{4}-[0-9]{2}' | sed -n '1p;$p' | tr -d '-' | paste -sd_)
[[ $tag =~ ^[0-9]{6}_[0-9]{6}$ ]] || fail "cannot read time range"

mkdir -p "$work/out" "$dest"
IFS=, read -ra vlist <<< "$vars"
for v in "${vlist[@]}"; do
  ncks -O -7 -L 1 -v "$v${xtr:+,$xtr}" "$work/rgr.nc" "$work/out/${v}_${tag}.nc" || fail "split $v"
done

# replace the truncated files of this stream (old names may carry another month range)
for v in "${vlist[@]}"; do rm -f "$dest/${v}_"[0-9][0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][0-9][0-9].nc; done
mv "$work"/out/*.nc "$dest"/ || fail "move"

nout=$(ls "$dest"/*_"${tag}".nc | wc -l)
(( nout == ${#vlist[@]} )) || fail "wrote $nout of ${#vlist[@]} variables"
echo "OK files=${nout} records=${nrec} range=${tag} inputs=${#files[@]} $(date -Iseconds)" > "$status"
echo "[OK] $case_name $ens $stream records=${nrec} range=${tag}"
