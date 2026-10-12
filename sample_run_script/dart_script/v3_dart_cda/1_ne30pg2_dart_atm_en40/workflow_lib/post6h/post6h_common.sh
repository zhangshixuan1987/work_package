#!/bin/bash
# Shared helpers for the 6-hourly post-processing step (9_run_post_6hourly.sh).
#
# Every DA cycle is a separate ${my_e3sm_cycle_hours}-h model run, so history
# averages longer than one cycle are never valid. All products here are built
# from one record per cycle, and a record is accepted only if its time bounds
# span exactly one cycle. Daily and monthly values are derived from those.

post6h_fail() { echo "ERROR: $*" >&2; exit 1; }

post6h_activate_env() {
  if [[ "${my_analysis_env_type:-}" == "custom" ]]; then
    source "${my_conda_setup_file}" || post6h_fail "could not source ${my_conda_setup_file}"
    conda activate "${my_analysis_conda_env}" || post6h_fail "could not activate ${my_analysis_conda_env}"
  else
    local pixi_dir="/share/apps/E3SM/conda_envs/e3smu_1_13_0/compy/pixi_login"
    [[ -x "${pixi_dir}/bin/pixi" ]] || post6h_fail "Pixi executable not found: ${pixi_dir}/bin/pixi"
    eval "$("${pixi_dir}/bin/pixi" shell-hook --as-is --shell bash --manifest-path "${pixi_dir}/pixi.toml")"
    hash -r
  fi
  local cmd
  for cmd in ncks ncrcat ncra ncap2 ncremap ncrename ncdump; do
    command -v "${cmd}" >/dev/null 2>&1 || post6h_fail "required command not found: ${cmd}"
  done
}

post6h_cycles_per_day() { echo $(( 24 / my_e3sm_cycle_hours )); }

post6h_days_in_month() { date -u -d "$1-01 +1 month -1 day" +%-d; }

# One stamp (YYYY-MM-DD-SSSSS) per DA cycle that starts in month $1 (YYYY-MM).
post6h_month_stamps() {
  local step=$(( my_e3sm_cycle_hours * 3600 )) t end
  t=$(date -u -d "$1-01 00:00:00" +%s)
  end=$(date -u -d "$1-01 +1 month" +%s)
  while (( t < end )); do
    printf '%s-%05d\n' "$(date -u -d "@${t}" +%Y-%m-%d)" $(( t % 86400 ))
    t=$(( t + step ))
  done
}

# Name of the time-bounds variable in $1 (EAM: time_bnds, ELM/MOSART: time_bounds).
post6h_bounds_name() {
  local b
  for b in time_bnds time_bounds; do
    ncks -m -v "${b}" "$1" >/dev/null 2>&1 && { echo "${b}"; return 0; }
  done
  return 1
}

# Index of the single record in $1 whose time bounds span exactly one cycle.
# Fails if there is none or more than one, or if time is not in days.
post6h_cycle_record() {
  local f=$1 b width hdr
  # Read the header first: grep -q on a pipe exits early and, under pipefail,
  # the resulting SIGPIPE in ncdump would make a passing check fail.
  hdr=$(ncdump -h "${f}") || return 1
  grep -qE 'time:units = "days since' <<<"${hdr}" || return 1
  b=$(post6h_bounds_name "${f}") || return 1
  width=$(awk -v h="${my_e3sm_cycle_hours}" 'BEGIN{printf "%.10f", h/24}')
  ncdump -v "${b}" "${f}" | sed -n "/^ ${b} =/,/;/p" | sed "s/^ ${b} =//" | tr -d ';\n' | tr ',' ' ' |
    awk -v w="${width}" '{
      n = 0
      for (i = 1; i < NF; i += 2) { d = ($(i+1) - $i) - w; if (d < 1e-6 && d > -1e-6) { idx = (i - 1) / 2; n++ } }
      if (n == 1) print idx; else exit 1
    }'
}

# Index of the single record in $1 valid at the file's own start time, for
# instantaneous streams (e.g. the ELM-DART forecast input written at the
# analysis time). The file's time units must be "days since <its stamp>".
post6h_instant_record() {
  local f=$1 st=$2 hdr ref
  hdr=$(ncdump -h "${f}") || return 1
  ref="$(date -u -d "${st%-*} 00:00:00 $((10#${st##*-})) seconds" '+%Y-%m-%d %H:%M:%S')"
  grep -qE "time:units = \"days since ${ref}(\.0)?\"" <<<"${hdr}" || return 1
  ncdump -v time "${f}" | sed -n '/^ time =/,/;/p' | sed 's/^ time =//' | tr -d ';\n' | tr ',' ' ' |
    awk '{ n = 0; for (i = 1; i <= NF; i++) if ($i < 1e-6 && $i > -1e-6) { idx = i - 1; n++ }
           if (n == 1) print idx; else exit 1 }'
}

# Number of records along the time dimension of $1.
post6h_ntime() {
  local hdr
  hdr=$(ncdump -h "$1") || return 1
  grep -m1 -oE 'time = UNLIMITED ; // \([0-9]+ currently|time = [0-9]+ ;' <<<"${hdr}" | grep -oE '[0-9]+'
}

# Publish $1 as $2: compress, verify it is readable with $3 records, then move into place.
post6h_publish() {
  local src=$1 dst=$2 nt=$3 tmp
  tmp="${src}.pub"
  ncks -O -4 -L 1 "${src}" "${tmp}" || return 1
  [[ "$(post6h_ntime "${tmp}")" == "${nt}" ]] || { echo "ERROR: ${dst##*/} has $(post6h_ntime "${tmp}") records, expected ${nt}" >&2; return 1; }
  mkdir -p "${dst%/*}" && mv -f "${tmp}" "${dst}"
}

# Write a completion record atomically: post6h_record <file> key=value ...
post6h_record() {
  local rec=$1; shift
  printf '%s\n' "$@" "job_id=${SLURM_JOB_ID:-none}" "completed_at=$(date '+%Y-%m-%d %H:%M:%S')" > "${rec}.tmp.$$" &&
    mv -f "${rec}.tmp.$$" "${rec}"
}
