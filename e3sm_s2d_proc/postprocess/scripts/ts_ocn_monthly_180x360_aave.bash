#!/bin/bash

# Running on pm-cpu

#SBATCH --job-name=ts_ocn_monthly_180x360_aave
#SBATCH --account=e3sm
#SBATCH --nodes=1
#SBATCH --output=ts_ocn_monthly_180x360_aave.o%j
#SBATCH --exclusive
#SBATCH --time=0:30:00
#SBATCH --qos=regular
#SBATCH --constraint=cpu

package="/global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh"
map="/global/cfs/cdirs/e3sm/diagnostics/maps/map_IcoswISC30E3r5_to_cmip6_180x360_traave.20240221.nc"

diri="/pscratch/sd/z/zhan391/e3sm_project/E3SMv3_S2D"
diro="/global/cfs/cdirs/e3smdata/simulations/S2S2D"

case="WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_JRA55_FOSIRL_1980050100"
ens="EN00"

comp="ocn"
ckey="mpasocean"
group="ts"
grid="180x360_aave"
freq="monthly"
hist="mpaso.hist.am.timeSeriesStatsMonthly"

yr_srt=1980
yr_end=1982
mth_srt=05
mth_end=04

ypf=$((yr_end-yr_srt+1))

ref_yr=1850

start_total=$(( yr_srt * 12 + 10#${mth_srt} ))
end_total=$(( yr_end * 12 + 10#${mth_end} ))
n_months=$(( end_total - start_total + 1 ))
yint=$(( n_months / 12 ))

if [[ -n "${ens}" ]]; then
  case_name="${case}.${ens}"
else
  case_name="${case}"
fi

ttag="$(printf "%04d" ${yr_srt})$(printf "%02d" ${mth_srt})-$(printf "%04d" ${yr_end})$(printf "%02d" ${mth_end})-$(printf "%04d" ${yint})"

debug=False
if [[ "${debug,,}" == "true" ]]; then
  set -x
fi

if [[ -n "${ens}" ]]; then
  scrip_dir="${diro}/${case}/${ens}/post/scripts"
else
  scrip_dir="${diro}/${case_name}/post/scripts"
fi

mkdir -p "${scrip_dir}"
cd "${scrip_dir}" || exit 1

id=${SLURM_JOBID}
STARTTIME=$(date +%s)

status_file="${group}_${comp}_${freq}_${grid}_${ttag}.status"
echo "RUNNING ${id}" > "${status_file}"

write_status() {
  echo "$1" > "${scrip_dir}/${status_file}"
}

fail() {
  local code=$1
  local message=$2

  echo "${message}"
  write_status "ERROR (${code})"
  exit "${code}"
}

on_error() {
  local rc=$?
  local line_no=$1

  if (( rc != 0 )); then
    echo "ERROR: Command failed at line ${line_no} with exit code ${rc}"
    write_status "ERROR (${rc})"
  fi

  exit "${rc}"
}

require_file() {
  local path=$1
  local description=$2

  if [[ ! -f "${path}" ]]; then
    fail 1 "ERROR: Missing ${description}: ${path}"
  fi
}

require_cmd() {
  local cmd=$1

  if ! command -v "${cmd}" >/dev/null 2>&1; then
    fail 1 "ERROR: Required command not found: ${cmd}"
  fi
}

select_input_dir() {
  local candidate
  local file

  shopt -s nullglob
  for candidate in "$@"; do
    [[ -d "${candidate}" ]] || continue

    for file in "${candidate}/${case_name}.${hist}."????-*.nc; do
      input="${candidate}"
      shopt -u nullglob
      return 0
    done
  done
  shopt -u nullglob

  return 1
}

set_pkg_manager() {
  if echo "source ${package}" | grep -Eq 'load_latest_e3sm_unified'; then
    pkg_manager="pixi"
  else
    pkg_manager="conda"
  fi

  echo "${pkg_manager} list python:"
  ${pkg_manager} list python || true
}

run_nco() {
  local environment_path=""
  local nco_cmd=$1
  shift

  if [[ -n "${environment_path}" && ( "${nco_cmd}" == "ncclimo" || "${nco_cmd}" == "ncremap" ) ]]; then
    "${environment_path%/}/${nco_cmd}" --npo "$@"
  else
    "${nco_cmd}" "$@"
  fi
}

read_first_nco_value() {
  local var_name=$1
  local nc_file=$2

  run_nco ncks -H -C -v "${var_name}" "${nc_file}" | \
    awk -v var="${var_name}" '
      index($0, "=") == 0 { next }
      first != "" { next }
      $0 ~ "^[[:space:]]*" var "([[:space:]]|\\[|=)" {
        value = $0
        sub(/^[^=]*=/, "", value)
        sub(/^[[:space:]]*/, "", value)
        split(value, values, /[,[:space:];]+/)
        first = values[1]
      }
      END { print first }
    '
}

require_cmd awk
require_cmd basename
require_cmd find
require_cmd sort
require_cmd tr

require_file "${package}" "environment package"
require_file "${map}" "mapping file"

source "${package}" || fail 1 "ERROR: Failed to source environment package: ${package}"

set_pkg_manager
echo "${pkg_manager} list nco:"
${pkg_manager} list nco || true

require_cmd ncclimo
require_cmd ncap2
require_cmd ncatted
require_cmd ncks
require_cmd ncrename
set -e
set -E
set -o pipefail
trap 'on_error ${LINENO}' ERR

workdir=$(mktemp -d "tmp.${group}_${comp}_${freq}_${grid}_${ttag}.${id:-nojob}.XXXXXX") || \
  fail 1 "ERROR: Failed to create workdir"

cd "${workdir}" || fail 1 "ERROR: Failed to cd into workdir"

vars=timeMonthly_avg_tThreshMLD,timeMonthly_avg_dThreshMLD,timeMonthly_avg_atmosphericPressure,timeMonthly_avg_ssh,timeMonthly_avg_pressureAdjustedSSH,timeMonthly_avg_oceanHeatContentSfcTo700m,timeMonthly_avg_oceanHeatContent700mTo2000m,timeMonthly_avg_oceanHeatContentSfcToBot

vext=xtime_startMonthly,xtime_endMonthly,time

if [[ -n "${ens}" ]]; then
  input_candidates=(
    "${diri}/${case}/${ens}/archive/${comp}/hist"
    "${diri}/${case}/${ens}/run"
  )
else
  input_candidates=(
    "${diri}/${case_name}/archive/${comp}/hist"
    "${diri}/${case_name}/run"
  )
fi

if ! select_input_dir "${input_candidates[@]}"; then
  fail 1 "ERROR: Missing input files matching ${case_name}.${hist}.????-*.nc in: ${input_candidates[*]}"
fi

echo "Input directory: ${input}"

# ------------------------------------------------------------
# Compute noleap offset for the unified time origin.
#
# Desired first timestamp:
#   ${yr_srt}-${mth_srt}-01 00:00:00
#
# Output time coordinate:
#   time = timeMonthly_avg_daysSinceStartOfSim + (start_offset - first_input_day)
#   time@units    = "days since 1850-01-01 00:00:00.0"
#   time@calendar = "noleap"
# ------------------------------------------------------------

sim_yr=$((10#${yr_srt}))
sim_mo=$((10#${mth_srt}))
sim_dy=1

month_days=(0 31 28 31 30 31 30 31 31 30 31 30 31)

doy=0
for (( m=1; m<sim_mo; m++ )); do
  doy=$(( doy + ${month_days[$m]} ))
done
doy=$(( doy + sim_dy - 1 ))

start_offset=$(( (sim_yr - ref_yr) * 365 + doy ))

echo "Time anchor: ${yr_srt}-${mth_srt}-01 00:00:00"
echo "Unified time units: days since 1850-01-01 00:00:00.0"
echo "Noleap start offset: ${start_offset}"

# Build a single global shift so all files share one continuous time axis.
shopt -s nullglob
first_input_file=""

anchor_ymd="$(printf "%04d-%02d-01" "${sim_yr}" "${sim_mo}")"
for candidate in "${input}/${case_name}.${hist}.${anchor_ymd}.nc" "${input}/${case_name}.${hist}.${anchor_ymd}_"*.nc; do
  [[ -f "${candidate}" ]] || continue
  first_input_file="${candidate}"
  break
done

# Fallback: use first file in requested year range.
if [[ -z "${first_input_file}" ]]; then
for (( year=${yr_srt}; year<=${yr_end}; year++ )); do
  YYYY=$(printf "%04d" ${year})
  for candidate in "${input}/${case_name}.${hist}.${YYYY}"-*.nc; do
    first_input_file="${candidate}"
    break
  done
  if [[ -n "${first_input_file}" ]]; then
    break
  fi
done
fi

if [[ -z "${first_input_file}" ]]; then
  fail 1 "ERROR: Missing input files in ${input}"
fi

first_sim_day=$(read_first_nco_value timeMonthly_avg_daysSinceStartOfSim "${first_input_file}")

if ! [[ "${first_sim_day}" =~ ^-?[0-9]+([.][0-9]+)?$ ]]; then
  fail 3 "ERROR: Failed to read numeric timeMonthly_avg_daysSinceStartOfSim from ${first_input_file}"
fi

if awk -v fsd="${first_sim_day}" 'BEGIN { exit !(fsd > 10000.0) }'; then
  time_shift=0
else
  time_shift=$(awk -v so="${start_offset}" -v fsd="${first_sim_day}" 'BEGIN{printf "%.15g", so-fsd}')
fi

echo "First simulation day in input: ${first_sim_day}"
echo "Time anchor file: ${first_input_file}"
echo "Global time shift applied: ${time_shift}"

input_end_yr=${yr_end}
input_end_mo=$((10#${mth_end} + 1))
if (( input_end_mo > 12 )); then
  input_end_mo=1
  input_end_yr=$((input_end_yr + 1))
fi

input_start_ymd="$(printf "%04d-%02d-01" "${yr_srt}" "$((10#${mth_srt}))")"
input_end_ymd="$(printf "%04d-%02d-01" "${input_end_yr}" "${input_end_mo}")"

echo "Preprocessing input dates: ${input_start_ymd} through ${input_end_ymd} (inclusive)"

mkdir -p preproc

for (( year=${yr_srt}; year<=${input_end_yr}; year++ )); do
  YYYY=$(printf "%04d" ${year})

  for file in "${input}/${case_name}.${hist}.${YYYY}"-*.nc; do
    bname=$(basename "${file}")

    if [[ "${bname}" =~ ([0-9]{4}-[0-9]{2}-[0-9]{2}) ]]; then
      file_ymd="${BASH_REMATCH[1]}"
    else
      fail 1 "ERROR: Could not parse date from input filename: ${bname}"
    fi

    if [[ "${file_ymd}" < "${input_start_ymd}" || "${file_ymd}" > "${input_end_ymd}" ]]; then
      continue
    fi

    outfile="preproc/${bname}"
    tmp_with_time="preproc/tmp_with_time_${bname}"

    echo "Preprocessing ${bname}"

    run_nco ncks -O \
      -v "${vars},xtime_startMonthly,xtime_endMonthly,timeMonthly_avg_daysSinceStartOfSim" \
      "${file}" "${outfile}"

    run_nco ncrename -O \
      -d Time,time \
      "${outfile}" "${outfile}"

    run_nco ncap2 -O \
      -s "time=timeMonthly_avg_daysSinceStartOfSim+${time_shift}; time@units=\"days since 1850-01-01 00:00:00.0\"; time@calendar=\"noleap\"; time@standard_name=\"time\"; time@long_name=\"time\"; time@axis=\"T\"" \
      "${outfile}" "${tmp_with_time}"

    run_nco ncks -O -x -v timeMonthly_avg_daysSinceStartOfSim \
      "${tmp_with_time}" "${outfile}"

    run_nco ncatted -O -a _FillValue,time,d,, \
      "${outfile}" "${outfile}"

    rm -f "${tmp_with_time}"
  done
done

shopt -u nullglob

find preproc -name "${case_name}.${hist}.????-*.nc" | sort > input.txt

if [[ ! -s input.txt ]]; then
  fail 1 "ERROR: Missing preprocessed input files"
fi

mkdir -p trash output

echo "ncclimo settings: vext=${vext}"

if ! cat input.txt | run_nco ncclimo \
  -c "${case_name}" \
  -v "${vars}" \
  --var_xtr="${vext}" \
  --split \
  --yr_srt="${yr_srt}" \
  --yr_end="${yr_end}" \
  --mth_srt="${mth_srt}" \
  --mth_end="${mth_end}" \
  --ypf="${ypf}" \
  --map="${map}" \
  -o trash \
  -O output \
  --prc_typ="${ckey}"
then
  fail 2 "ERROR: ncclimo failed"
fi

if [[ -n "${ens}" ]]; then
  dest="${diro}/${case}/${ens}/post/${comp}/${grid}/${group}/${freq}/${yint}yr"
else
  dest="${diro}/${case_name}/post/${comp}/${grid}/${group}/${freq}/${yint}yr"
fi

mkdir -p "${dest}"

shopt -s nullglob
output_files=(output/*.nc)
shopt -u nullglob

if (( ${#output_files[@]} == 0 )); then
  fail 6 "ERROR: ncclimo produced no output files"
fi

for file in "${output_files[@]}"; do
  run_nco ncatted -O -a _FillValue,time,d,, \
    "${file}" "${file}"
done

if ! mv "${output_files[@]}" "${dest}"; then
  fail 6 "ERROR: Failed to move output files to ${dest}"
fi

cd ..

if [[ "${debug,,}" != "true" ]]; then
  rm -rf "${workdir}"
fi

ENDTIME=$(date +%s)
ELAPSEDTIME=$((ENDTIME - STARTTIME))

echo "=============================================="
echo "Elapsed time: ${ELAPSEDTIME} seconds"
echo "=============================================="

cd "${scrip_dir}"
rm -f "${status_file}"
write_status "OK"

exit 0
