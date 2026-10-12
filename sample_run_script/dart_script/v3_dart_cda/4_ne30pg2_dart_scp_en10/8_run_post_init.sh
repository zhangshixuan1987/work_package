#!/bin/bash
#SBATCH --account=esmd
#SBATCH --time=24:00:00
#SBATCH --partition=slurm
#SBATCH --job-name=e3sm_dart_interp_init
#SBATCH --nodes=1
#SBATCH --output=e3sm_dart_interp_init.%j
#SBATCH --exclusive
#SBATCH --no-kill
#SBATCH --requeue

set -Eeuo pipefail
shopt -s nullglob

INTERP_START="2011-12-01-21600"
INTERP_END="2011-12-01-21600"
MAX_CONCURRENT_WORKERS=4
OVERWRITE_EXISTING="FALSE"
REQUIRE_UPSTREAM_MARKERS="TRUE"
COMPONENT="${1:-BOTH}" # EAM, ELM, or BOTH; first argument overrides this value.
EAM_REGRID="TRUE"
ELM_REGRID="FALSE"
readonly -a EAM_PASSES=(eam_atmosphere eam_land_cross)
readonly -a ELM_PASSES=(elm_land elm_atmosphere_cross)
readonly -a EAM_FILE_PATTERNS=("*.eam_????.e.*" "*.eam_????.i.*" "*.dart.i.eam_output_mean.*" "*.dart.i.eam_output_sd.*" "*.dart.e.eam_forecast_mean.*" "*.dart.e.eam_forecast_sd.*")
readonly -a ELM_FILE_PATTERNS=("preassim_member_*_d??.nc" "preassim_mean_d??.nc" "preassim_sd_d??.nc" "output_mean_d??.nc" "output_sd_d??.nc")
readonly EAM_VARIABLES="lat_d,lon_d,U,V,T,Q,PS"

fail() { echo "ERROR: $*" >&2; exit 1; }
normalize_bool() { local n="$1" v="${2^^}"; [[ "$v" == TRUE || "$v" == FALSE ]] || fail "$n must be TRUE or FALSE"; printf %s "$v"; }
stamp_to_epoch() {
  local s="$1" d t e
  [[ "$s" =~ ^([0-9]{4}-[0-9]{2}-[0-9]{2})-([0-9]{5})$ ]] || return 1
  d=${BASH_REMATCH[1]}; t=${BASH_REMATCH[2]}; (( 10#$t < 86400 )) || return 1
  e=$(date -u -d "$d 00:00:00" +%s) || return 1; echo $((e + 10#$t))
}
epoch_to_stamp() { date -u -d "@$1" '+%Y-%m-%d %H %M %S' | awk '{printf "%s-%05d\n",$1,$2*3600+$3*60+$4}'; }
record_value() { awk -F= -v k="$2" '$1==k {sub(/^[^=]*=/,""); print; exit}' "$1"; }

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then WORK_DIR=$(readlink -f "$SLURM_SUBMIT_DIR"); else SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}"); WORK_DIR=$(dirname "$SCRIPT_PATH"); fi
cd "$WORK_DIR"
[[ -r create_and_setup_case.sh ]] || fail "missing create_and_setup_case.sh"
source ./create_and_setup_case.sh
[[ -r "${my_conda_setup_file:-}" ]] || fail "invalid Conda setup: ${my_conda_setup_file:-unset}"
source "$my_conda_setup_file"
conda activate "$my_analysis_conda_env" || fail "cannot activate $my_analysis_conda_env"
for cmd in awk date dirname flock ln mkdir mktemp mv ncdump ncks ncrename readlink rm; do command -v "$cmd" >/dev/null || fail "missing command: $cmd"; done
(( $# <= 1 )) || fail "usage: $0 [EAM|ELM|BOTH]"
COMPONENT=${COMPONENT^^}
case "$COMPONENT" in EAM) PROCESS_EAM=TRUE; PROCESS_ELM=FALSE;; ELM) PROCESS_EAM=FALSE; PROCESS_ELM=TRUE;; BOTH) PROCESS_EAM=TRUE; PROCESS_ELM=TRUE;; *) fail "invalid component: $COMPONENT";; esac
for n in OVERWRITE_EXISTING REQUIRE_UPSTREAM_MARKERS EAM_REGRID ELM_REGRID; do printf -v "$n" %s "$(normalize_bool "$n" "${!n}")"; done
[[ "$MAX_CONCURRENT_WORKERS" =~ ^[1-9][0-9]*$ ]] || fail "invalid MAX_CONCURRENT_WORKERS"
[[ "$PROCESS_ELM" != TRUE || "$ELM_REGRID" == FALSE ]] || fail "ELM_REGRID=TRUE is unsupported for domain-split files"
if [[ "$PROCESS_EAM" == TRUE && "$EAM_REGRID" == TRUE ]]; then command -v ncremap >/dev/null || fail "missing ncremap"; [[ -s "${my_eam_init_mapping_file:-}" ]] || fail "missing EAM map"; fi
START_EPOCH=$(stamp_to_epoch "$INTERP_START") || fail "invalid INTERP_START"
END_EPOCH=$(stamp_to_epoch "$INTERP_END") || fail "invalid INTERP_END"
(( START_EPOCH <= END_EPOCH )) || fail "start is after end"
OUT_ROOT="$my_dart_root/post/diag_processed"; STATUS_DIR="$my_status_dir"
mkdir -p "$OUT_ROOT" "$STATUS_DIR" "$my_log_dir" "$my_lock_dir"
exec 8>"$my_lock_dir/post_diag_process.lock"; flock -n 8 || fail "another Step 8 job is active"

validate_marker() {
  local s="$1" m
  if [[ "$s" == "$my_e3sm_start_date-$my_e3sm_start_tod" ]]; then m="$STATUS_DIR/perturb_complete.$s"; else m="$STATUS_DIR/cycle_complete.$s"; fi
  if [[ ! -s "$m" ]]; then [[ "$REQUIRE_UPSTREAM_MARKERS" == FALSE ]] || fail "missing marker: $m"; return; fi
  [[ "$(record_value "$m" valid_time)" == "$s" ]] || fail "bad marker: $m"
}
validate_output() { [[ -s "$2" ]] && ncdump -h "$2" >/dev/null 2>&1 && { [[ "$1" != eam ]] || ncks -m -v lat,lon,U,V,T,Q,PS "$2" >/dev/null 2>&1; }; }
process_file() {
  local c="$1" src="$2" dst="$3" d tmp
  d=$(dirname "$dst"); mkdir -p "$d"
  if [[ "$OVERWRITE_EXISTING" == FALSE ]] && validate_output "$c" "$dst"; then echo "Skipping: $dst"; return; fi
  tmp=$(mktemp -d "$d/.diag.XXXXXX") || return 1
  if [[ "$c" == eam ]]; then
    ncks -O -v "$EAM_VARIABLES" "$src" "$tmp/native.nc" && ncrename -O -v lat_d,lat -v lon_d,lon -d ncol_d,ncol "$tmp/native.nc" || { rm -rf "$tmp"; return 1; }
    if [[ "$EAM_REGRID" == TRUE ]]; then ncremap -m "$my_eam_init_mapping_file" -i "$tmp/native.nc" -o "$tmp/final.nc" || { rm -rf "$tmp"; return 1; }; else mv "$tmp/native.nc" "$tmp/final.nc"; fi
  else
    ncks -O "$src" "$tmp/final.nc" || { rm -rf "$tmp"; return 1; }
  fi
  validate_output "$c" "$tmp/final.nc" || { rm -rf "$tmp"; return 1; }
  mv "$tmp/final.nc" "$dst"; rm -f "$tmp/native.nc"; rmdir "$tmp"; echo "Completed: $src"
}

declare -a TASK_C=() TASK_I=() TASK_O=() ALIAS_SRC=() ALIAS_DST=()
declare -A PRIMARY=()
add_tasks() {
  local c="$1" stamp="$2" rg="$3" pass dir pattern src key dst mode
  local -a passes patterns matches
  if [[ "$c" == eam ]]; then passes=("${EAM_PASSES[@]}"); patterns=("${EAM_FILE_PATTERNS[@]}"); else passes=("${ELM_PASSES[@]}"); patterns=("${ELM_FILE_PATTERNS[@]}"); fi
  mode=native; [[ "$rg" == TRUE ]] && mode=regridded
  for pass in "${passes[@]}"; do
    dir="$my_dart_root/$c/$stamp/$pass"; [[ -d "$dir" ]] || fail "missing pass: $dir"
    for pattern in "${patterns[@]}"; do
      matches=("$dir"/$pattern); (( ${#matches[@]} )) || { echo "WARNING: no match: $dir/$pattern" >&2; continue; }
      for src in "${matches[@]}"; do
        [[ -e "$src" ]] || fail "broken link: $src"; ncdump -h "$src" >/dev/null 2>&1 || fail "invalid NetCDF: $src"
        [[ "$c" != eam ]] || ncks -m -v "$EAM_VARIABLES" "$src" >/dev/null 2>&1 || fail "missing EAM variables: $src"
        key="$c:$(readlink -f "$src")"; dst="$OUT_ROOT/$mode/$c/$stamp/$pass/$(basename "$src")"
        if [[ -n "${PRIMARY[$key]:-}" ]]; then ALIAS_SRC+=("${PRIMARY[$key]}"); ALIAS_DST+=("$dst"); else PRIMARY[$key]="$dst"; TASK_C+=("$c"); TASK_I+=("$src"); TASK_O+=("$dst"); fi
      done
    done
  done
}
run_tasks() {
  local j k rc bad=0; local -a pids=() labels=()
  for j in "${!TASK_I[@]}"; do
    process_file "${TASK_C[$j]}" "${TASK_I[$j]}" "${TASK_O[$j]}" & pids+=("$!"); labels+=("${TASK_I[$j]}")
    if (( ${#pids[@]} >= MAX_CONCURRENT_WORKERS )); then for k in "${!pids[@]}"; do wait "${pids[$k]}" || { rc=$?; echo "ERROR ($rc): ${labels[$k]}" >&2; bad=1; }; done; pids=(); labels=(); fi
  done
  for k in "${!pids[@]}"; do wait "${pids[$k]}" || { rc=$?; echo "ERROR ($rc): ${labels[$k]}" >&2; bad=1; }; done
  (( bad == 0 )) || fail "diagnostic processing failed"
}

for ((epoch=START_EPOCH; epoch<=END_EPOCH; epoch+=86400)); do
  stamp=$(epoch_to_stamp "$epoch"); validate_marker "$stamp"
  TASK_C=(); TASK_I=(); TASK_O=(); ALIAS_SRC=(); ALIAS_DST=(); PRIMARY=()
  [[ "$PROCESS_EAM" != TRUE ]] || add_tasks eam "$stamp" "$EAM_REGRID"
  [[ "$PROCESS_ELM" != TRUE ]] || add_tasks elm "$stamp" "$ELM_REGRID"
  (( ${#TASK_I[@]} )) || fail "no requested files at $stamp"
  echo "Processing ${#TASK_I[@]} unique files; ${#ALIAS_DST[@]} duplicate targets"
  run_tasks
  for j in "${!ALIAS_DST[@]}"; do mkdir -p "$(dirname "${ALIAS_DST[$j]}")"; ln -f "${ALIAS_SRC[$j]}" "${ALIAS_DST[$j]}" || fail "alias failed: ${ALIAS_DST[$j]}"; done
  marker="$STATUS_DIR/post_diag_complete.$COMPONENT.$stamp"
  printf 'valid_time=%s\ncomponent=%s\nunique_files=%s\nduplicate_targets=%s\ncompleted_at=%s\n' "$stamp" "$COMPONENT" "${#TASK_I[@]}" "${#ALIAS_DST[@]}" "$(date '+%F %T')" > "$marker.tmp"
  mv "$marker.tmp" "$marker"
done
echo "Step 8 completed: $INTERP_START through $INTERP_END"
echo "Outputs: $OUT_ROOT"
