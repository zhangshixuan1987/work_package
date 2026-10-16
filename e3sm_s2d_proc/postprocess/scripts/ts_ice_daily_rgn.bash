#!/bin/bash

# Running on pm-cpu

#SBATCH  --job-name=ts_ice_daily_rgn
#SBATCH  --account=e3sm
#SBATCH  --nodes=1
#SBATCH  --output=ts_ice_daily_rgn.o%j
#SBATCH  --exclusive
#SBATCH  --time=0:30:00


#SBATCH  --qos=regular
#SBATCH  --constraint=cpu


######################
package="/global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh"
plevdata="/global/cfs/cdirs/e3sm/diagnostics/e3sm_to_cmip_data/grids/vrt_remap_plev19.nc"
map="/global/cfs/cdirs/e3sm/diagnostics/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc"
diri="/pscratch/sd/z/zhan391/e3sm_project/E3SMv3_S2D"
diro="/global/cfs/cdirs/e3smdata/simulations/S2S2D"
case="WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_JRA55_FOSIRL_1980050100"
ens="EN00"
comp="ice"
ckey="mpassi"
group="ts"
grid="rgn"
freq="daily"
hist="mpassi.hist.am.regionalStatistics"

yr_srt=1980
yr_end=1982

mth_srt=05
mth_end=04

ref_yr=1850

# number of complete 12-month periods
start_total=$(( yr_srt * 12 + 10#${mth_srt} ))
end_total=$(( yr_end * 12 + 10#${mth_end} ))
n_months=$(( end_total - start_total +1 ))
yint=$(( n_months / 12 ))

if [[ -n "$ens" ]]; then
    case_name="${case}.${ens}"
else
    case_name="${case}"
fi
ttag="`printf "%04d" ${yr_srt}``printf "%02d" ${mth_srt}`-`printf "%04d" ${yr_end}``printf "%02d" ${mth_end}`-`printf "%04d" ${yint}`"

#######################

# Turn on debug output if needed
debug=False
if [[ "${debug,,}" == "true" ]]; then
  set -x
fi

# Script dir
if [[ -n "$ens" ]]; then
  scrip_dir="${diro}/${case}/${ens}/post/scripts"
else
  scrip_dir="${diro}/${case_name}/post/scripts"
fi
if [ ! -d "${scrip_dir}" ];then 
   mkdir -p ${scrip_dir}
fi 
cd ${scrip_dir}

# Get jobid
id=${SLURM_JOBID}

# Update status file
STARTTIME=$(date +%s)
echo "RUNNING ${id}" > ${group}_${comp}_${freq}_${grid}_${ttag}.status

set_pkg_manager() {
  # Detect whether to use pixi or conda based on environment_commands
  # Example environment commands:

  # These would indicate pixi is the package manager:
  # source /share/apps/E3SM/conda_envs/test_e3sm_unified_1.13.0rc7_compy.sh
  # source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh

  # This would indicate conda is the package manager:
  # source /qfs/people/fors729/miniforge3/etc/profile.d/conda.sh; conda activate e3sm_to_cmip_20260415

  # Important: we need the quotes since environment_commands might
  # include multiple commands, separated by semi-colon.
  if echo "source ${package}" | grep -Eq 'test_e3sm_unified|load_latest_e3sm_unified'; then
    pkg_manager="pixi"
  else
    pkg_manager="conda"
  fi

  # We always want to know what the python version is.
  echo "${pkg_manager} list python:"
  ${pkg_manager} list python || true # If we can't print this, just continue on.
}

run_nco() {
  # Example:
  # environment_path=/global/cfs/cdirs/e3sm/zender/bin_perlmutter/ (set with jinja2 substitution from config file)
  # nco_cmd=ncclimo (passed in directly)
  # Runs `/global/cfs/cdirs/e3sm/zender/bin_perlmutter/ncclimo --npo` on the original arguments

  # Full list of NCO commands: https://nco.sourceforge.net/nco.html#Reference-Manual

  local environment_path=""
  local nco_cmd=$1
  shift # Remove nco_cmd from the argument list

  if [[ -n "$environment_path" && ( "$nco_cmd" == "ncclimo" || "$nco_cmd" == "ncremap" ) ]]; then
    # nco_path is non-empty AND nco_cmd is a script (either ncclimo or ncremap),
    # so use the full development setup:
    "${environment_path%/}/$nco_cmd" --npo "$@"
  else
    # nco_path is empty OR nco_cmd is a binary executable (not a script),
    # so use just the nco_cmd itself:
    "$nco_cmd" "$@"
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

noleap_day_offset() {
  local yr=$1
  local mo=$2
  local dy=$3

  awk -v ref_yr="${ref_yr}" -v yr="${yr}" -v mo="${mo}" -v dy="${dy}" '
    BEGIN {
      md[1]=31; md[2]=28; md[3]=31; md[4]=30; md[5]=31; md[6]=30
      md[7]=31; md[8]=31; md[9]=30; md[10]=31; md[11]=30; md[12]=31

      doy = dy - 1
      for (m = 1; m < mo; m++) {
        doy += md[m]
      }

      printf "%.15g", (yr - ref_yr) * 365 + doy
    }
  '
}

set -e
source ${package}
set +e

set_pkg_manager
echo "${pkg_manager} list nco:"
${pkg_manager} list nco || true # If we can't print this, just continue on.

# Create temporary workdir
hash=`mktemp --dry-run -d XXXX`
workdir=tmp.${group}_${comp}_${freq}_${grid}_${ttag}.${id}.${hash}
mkdir ${workdir} || { echo 'ERROR: Failed to create workdir'; echo 'ERROR (1)' > ${scrip_dir}/${group}_${comp}_${freq}_${grid}_${ttag}.status; exit 1; }
cd ${workdir} || { echo 'ERROR: Failed to cd into workdir'; echo 'ERROR (1)' > ${scrip_dir}/${group}_${comp}_${freq}_${grid}_${ttag}.status; exit 1; }

if [[ -n "$ens" ]]; then
  input="${diri}/${case}/${ens}/archive/${comp}/hist"
else
  input="${diri}/${case_name}/archive/${comp}/hist"
fi

vars=totalIceArea,totalIceExtent,totalIceVolume,totalSnowVolume,totalKineticEnergy,rmsIceSpeed,averageAlbedo,maximumIceVolume,maximumIcePressure

vext=xtime

# Collect and time-fix files month by month over the requested date range.
tmp_fix_dir="tmp_timefix_${ttag}"
rm -rf "${tmp_fix_dir}"
mkdir -p "${tmp_fix_dir}"

echo "Time source: ${vext}"
echo "Unified time units: days since 1850-01-01 00:00:00.0"

input_end_yr=${yr_end}
input_end_mo=$((10#${mth_end} + 1))
if (( input_end_mo > 12 )); then
  input_end_mo=1
  input_end_yr=$((input_end_yr + 1))
fi

echo "Preprocessing input months: $(printf "%04d.%02d" "${yr_srt}" "$((10#${mth_srt}))") through $(printf "%04d.%02d" "${input_end_yr}" "${input_end_mo}") (inclusive)"
echo "Keeping boundary record: $(printf "%04d-%02d-01_00:00:00" "${input_end_yr}" "${input_end_mo}")"

fixed_files=()
cur_year=${yr_srt}
cur_month=${mth_srt}

while (( cur_year < input_end_yr || (cur_year == input_end_yr && cur_month <= input_end_mo) )); do
  YYYY=$(printf "%04d" ${cur_year})
  MM=$(printf "%02d" ${cur_month})

  f="${input}/${case_name}.${hist}.${YYYY}.${MM}.nc"

  if [[ -f "$f" ]]; then
    bname=$(basename "$f")
    fout="${tmp_fix_dir}/${bname}"
    tmp_with_time="${tmp_fix_dir}/tmp_with_time_${bname}"
    tmp_trimmed="${tmp_fix_dir}/tmp_trimmed_${bname}"

    echo "Processing: ${bname}"

    run_nco ncks -O \
      -v "${vars},${vext}" \
      "$f" "$fout"

    run_nco ncrename -O -d Time,time "$fout" "$fout" 2>/dev/null || true

    time_values=$(xtime_ncap_array "${vext}" "$f") || {
      echo "ERROR: Failed to decode ${vext} from ${f}"
      echo 'ERROR (3)' > ${scrip_dir}/${group}_${comp}_${freq}_${grid}_${ttag}.status
      exit 3
    }

    run_nco ncap2 -O \
      -s "time[\$time]=${time_values}; time@units=\"days since 1850-01-01 00:00:00.0\"; time@calendar=\"noleap\"; time@standard_name=\"time\"; time@long_name=\"time\"; time@axis=\"T\"" \
      "$fout" "$tmp_with_time"

    if (( cur_year == input_end_yr && cur_month == input_end_mo )); then
      run_nco ncks -O -d time,0,0 \
        "$tmp_with_time" "$tmp_trimmed"
      run_nco ncks -O "$tmp_trimmed" "$fout"
    else
      run_nco ncks -O "$tmp_with_time" "$fout"
    fi

    run_nco ncatted -O -a _FillValue,time,d,, \
      "$fout" "$fout"

    rm -f "$tmp_with_time" "$tmp_trimmed"

    fixed_files+=("$fout")
  else
    echo "WARNING: Missing file: ${f}"
  fi

  # Advance to next month
  cur_month=$(( cur_month + 1 ))
  if (( cur_month > 12 )); then
    cur_month=1
    cur_year=$(( cur_year + 1 ))
  fi
done

if [[ ${#fixed_files[@]} -eq 0 ]]; then
  cd ${scrip_dir}
  echo 'ERROR: No input files found'
  echo 'ERROR (1)' > ${group}_${comp}_${freq}_${grid}_${ttag}.status
  exit 1
fi

# Generate per-variable time-series files using ncrcat
mkdir -p output
IFS=',' read -ra var_list <<< "${vars}"
ts_err=0
for var in "${var_list[@]}"; do
  outfile="output/${var}_${ttag}.nc"
  # Include extra time strings when present; fall back without them if missing.
  run_nco ncrcat -O -v "${var},time,${vext}" "${fixed_files[@]}" "${outfile}" 2>/dev/null || \
  run_nco ncrcat -O -v "${var},time" "${fixed_files[@]}" "${outfile}"
  if [ $? != 0 ]; then
    echo "ERROR: ncrcat failed for variable ${var}"
    ts_err=1
  else
    run_nco ncatted -O -a _FillValue,time,d,, "${outfile}" "${outfile}"
  fi
done

if [ "${ts_err}" == "1" ]; then
  cd ${scrip_dir}
  echo 'ERROR (2)' > ${group}_${comp}_${freq}_${grid}_${ttag}.status
  exit 2
fi

missing_outputs=()
for var in "${var_list[@]}"; do
  if [[ ! -f "output/${var}_${ttag}.nc" ]]; then
    missing_outputs+=("${var}")
  fi
done

if (( ${#missing_outputs[@]} > 0 )); then
  cd ${scrip_dir}
  echo "ERROR: Missing output files for variables: ${missing_outputs[*]}"
  echo 'ERROR (3)' > ${group}_${comp}_${freq}_${grid}_${ttag}.status
  exit 3
fi

echo "Output files ready:"
printf '  %s\n' output/*.nc

# Move output ts files to final destination
{
  if [[ -n "$ens" ]]; then
    dest="${diro}/${case}/${ens}/post/${comp}/${grid}/${group}/${freq}/${yint}yr"
  else
    dest="${diro}/${case_name}/post/${comp}/${grid}/${group}/${freq}/${yint}yr"
  fi
  mkdir -p ${dest}
  mv output/*.nc ${dest}
}
if [ $? != 0 ]; then
  cd ${scrip_dir}
  echo 'ERROR (3)' > ${group}_${comp}_${freq}_${grid}_${ttag}.status
  exit 3
fi

# Delete temporary workdir
cd ..
if [[ "${debug,,}" != "true" ]]; then
  rm -rf ${workdir}
fi

# Update status file and exit

ENDTIME=$(date +%s)
ELAPSEDTIME=$(($ENDTIME - $STARTTIME))

echo ==============================================
echo "Elapsed time: $ELAPSEDTIME seconds"
echo ==============================================
cd ${scrip_dir}
rm -f ${group}_${comp}_${freq}_${grid}_${ttag}.status
echo 'OK' > ${group}_${comp}_${freq}_${grid}_${ttag}.status
exit 0
