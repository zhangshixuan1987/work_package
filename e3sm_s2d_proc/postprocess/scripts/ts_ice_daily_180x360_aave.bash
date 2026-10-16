#!/bin/bash

# Running on pm-cpu

#SBATCH --job-name=ts_ice_daily_180x360_aave
#SBATCH --account=e3sm
#SBATCH --nodes=1
#SBATCH --output=ts_ice_daily_180x360_aave.o%j
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

comp="ice"
ckey="mpasseaice"
group="ts"
grid="180x360_aave"
freq="daily"
hist="mpassi.hist.am.timeSeriesStatsDaily"

yr_srt=1980
yr_end=1982
mth_srt=05
mth_end=04

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
expected_span="$(printf "%04d%02d_%04d%02d" "${yr_srt}" "$((10#${mth_srt}))" "${yr_end}" "$((10#${mth_end}))")"

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

xtime_ncap_array() {
  local var_name=$1
  local nc_file=$2

  run_nco ncks -H -C -v "${var_name}" "${nc_file}" | \
    awk -v ref_yr="${ref_yr}" '
      BEGIN {
        md[1]=31; md[2]=28; md[3]=31; md[4]=30; md[5]=31; md[6]=30
        md[7]=31; md[8]=31; md[9]=30; md[10]=31; md[11]=30; md[12]=31
      }
      {
        line = $0
        while (match(line, /"[0-9]{4}-[0-9]{2}-[0-9]{2}_[0-9]{2}:[0-9]{2}:[0-9]{2}"/)) {
          stamp = substr(line, RSTART + 1, RLENGTH - 2)

          yr = substr(stamp, 1, 4) + 0
          mo = substr(stamp, 6, 2) + 0
          dy = substr(stamp, 9, 2) + 0
          hr = substr(stamp, 12, 2) + 0
          mn = substr(stamp, 15, 2) + 0
          sc = substr(stamp, 18, 2) + 0

          doy = dy - 1
          for (m = 1; m < mo; m++) {
            doy += md[m]
          }

          val = (yr - ref_yr) * 365 + doy + (hr * 3600 + mn * 60 + sc) / 86400.0
          values[++n] = sprintf("%.15g", val) "d"

          line = substr(line, RSTART + RLENGTH)
        }
      }
      END {
        if (n == 0) {
          exit 1
        }

        printf "{"
        for (i = 1; i <= n; i++) {
          printf "%s%s", (i == 1 ? "" : ","), values[i]
        }
        print "}"
      }
    '
}

require_cmd awk
require_cmd basename
require_cmd find
require_cmd grep
require_cmd sort
require_cmd tr

require_file "${package}" "environment package"
require_file "${map}" "mapping file"

source "${package}" || fail 1 "ERROR: Failed to source environment package: ${package}"

set_pkg_manager
echo "${pkg_manager} list nco:"
${pkg_manager} list nco || true

require_cmd ncap2
require_cmd ncatted
require_cmd ncks
require_cmd ncrcat
require_cmd ncrename
require_cmd ncremap

set -e
set -E
set -o pipefail
trap 'on_error ${LINENO}' ERR

workdir=$(mktemp -d "tmp.${group}_${comp}_${freq}_${grid}_${ttag}.${id:-nojob}.XXXXXX") || \
  fail 1 "ERROR: Failed to create workdir"

cd "${workdir}" || fail 1 "ERROR: Failed to cd into workdir"

vars=timeDaily_avg_iceAreaCell,timeDaily_avg_iceVolumeCell,timeDaily_avg_snowVolumeCell

vext=xtime_startDaily,xtime_endDaily
remap_flags=(--prc_typ=mpasocean --d2f)

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
# Compute a CF-style noleap time coordinate from MPAS-Seaice timestamps.
#
# Output time coordinate:
#   time = noleap day offset decoded from xtime_endDaily
#   time@units    = "days since 1850-01-01 00:00:00.0"
#   time@calendar = "noleap"
# ------------------------------------------------------------

echo "Time source: xtime_endDaily"
echo "Unified time units: days since 1850-01-01 00:00:00.0"

shopt -s nullglob

input_end_yr=${yr_end}
input_end_mo=$((10#${mth_end}))

input_start_ymd="$(printf "%04d-%02d-01" "${yr_srt}" "$((10#${mth_srt}))")"
input_end_ymd="$(printf "%04d-%02d-01" "${input_end_yr}" "${input_end_mo}")"

echo "Preprocessing input dates: ${input_start_ymd} through ${input_end_ymd} (inclusive)"

mkdir -p preproc

for (( year=${yr_srt}; year<=${input_end_yr}; year++ )); do
  YYYY=$(printf "%04d" ${year})

  for file in "${input}/${case_name}.${hist}.${YYYY}"-*.nc; do
    bname=$(basename "${file}")

    if [[ "${bname}" =~ ([0-9]{4})-([0-9]{2})-([0-9]{2}) ]]; then
      file_yr="${BASH_REMATCH[1]}"
      file_mo="${BASH_REMATCH[2]}"
      file_dy="${BASH_REMATCH[3]}"
      file_ymd="${file_yr}-${file_mo}-${file_dy}"
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
      -v "${vars},xtime_startDaily,xtime_endDaily" \
      "${file}" "${outfile}"

    run_nco ncrename -O \
      -d Time,time \
      "${outfile}" "${outfile}"

    time_values=$(xtime_ncap_array xtime_endDaily "${file}") || \
      fail 3 "ERROR: Failed to decode xtime_endDaily from ${file}"

    run_nco ncap2 -O \
      -s "time[\$time]=${time_values}; time@units=\"days since 1850-01-01 00:00:00.0\"; time@calendar=\"noleap\"; time@standard_name=\"time\"; time@long_name=\"time\"; time@axis=\"T\"" \
      "${outfile}" "${tmp_with_time}"

    run_nco ncks -O "${tmp_with_time}" "${outfile}"

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

echo "Native concat settings: vext=${vext}, expected_span=${expected_span}"
echo "Manual regrid settings: ${remap_flags[*]}, map=${map}"

mapfile -t preproc_files < input.txt

if (( ${#preproc_files[@]} == 0 )); then
  fail 1 "ERROR: Missing preprocessed input files"
fi

for var in ${vars//,/ }; do
  native_f="trash/${var}_${expected_span}.nc"

  echo "Concatenating ${var} to ${native_f}"
  if ! run_nco ncrcat -O \
    -v "${var},time,${vext}" \
    "${preproc_files[@]}" \
    "${native_f}"
  then
    fail 2 "ERROR: ncrcat failed for ${var}"
  fi

  run_nco ncatted -O -a _FillValue,time,d,, \
    "${native_f}" "${native_f}"
done

shopt -s nullglob
native_files=(trash/*.nc)
shopt -u nullglob

if (( ${#native_files[@]} == 0 )); then
  fail 6 "ERROR: ncrcat produced no native split files"
fi

for file in "${native_files[@]}"; do
  bname=$(basename "${file}")
  tmp_f="tmp_regrid_${bname}"
  out_f="output/${bname}"

  echo "Regridding ${bname}"
  cp "${file}" "${tmp_f}"

  if ! run_nco ncremap "${remap_flags[@]}" --no_stagger -t 1 \
    --nco_opt='--no_tmp_fl --hdr_pad=10000' \
    --map="${map}" \
    "${tmp_f}" "${out_f}"
  then
    fail 4 "ERROR: ncremap failed for ${bname}"
  fi

  for tvar in time xtime_startDaily xtime_endDaily; do
    if run_nco ncks -m "${tmp_f}" | grep -qE "^[[:space:]]*(char|float|double|int)[[:space:]]+${tvar}(\\(|\\[)" && \
       ! run_nco ncks -m "${out_f}" | grep -qE "^[[:space:]]*(char|float|double|int)[[:space:]]+${tvar}(\\(|\\[)"; then
      run_nco ncks -A -v "${tvar}" "${tmp_f}" "${out_f}"
    fi
  done

  if run_nco ncks -m "${out_f}" | grep -qE "^[[:space:]]*(float|double)[[:space:]]+time(\\(|\\[)"; then
    run_nco ncatted -O -a _FillValue,time,d,, \
      "${out_f}" "${out_f}"
  fi

  rm -f "${tmp_f}"
done

shopt -s nullglob
output_files=(output/*.nc)
shopt -u nullglob

if (( ${#output_files[@]} == 0 )); then
  fail 6 "ERROR: ncremap produced no output files"
fi

missing_outputs=()
for var in ${vars//,/ }; do
  if [[ ! -f "output/${var}_${expected_span}.nc" ]]; then
    missing_outputs+=("${var}")
  fi
done

if (( ${#missing_outputs[@]} > 0 )); then
  fail 6 "ERROR: Missing expected ${expected_span} output files for variables: ${missing_outputs[*]}"
fi

echo "Output files ready:"
printf '  %s\n' "${output_files[@]}"

# Remove _FillValue from time after ncrcat/ncremap,
# because ncremap may add it back to all float/double variables.
for file in "${output_files[@]}"; do
  if run_nco ncks -m "${file}" | grep -qE "^[[:space:]]*(float|double)[[:space:]]+time(\(|\[)"; then
    run_nco ncatted -O -a _FillValue,time,d,, \
      "${file}" "${file}"
  fi
done

if [[ -n "${ens}" ]]; then
  dest="${diro}/${case}/${ens}/post/${comp}/${grid}/${group}/${freq}/${yint}yr"
else
  dest="${diro}/${case_name}/post/${comp}/${grid}/${group}/${freq}/${yint}yr"
fi

mkdir -p "${dest}"

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
