#!/bin/bash

# Running on pm-cpu

#SBATCH  --job-name=ts_ocn_monthly_rgn
#SBATCH  --account=e3sm
#SBATCH  --nodes=1
#SBATCH  --output=ts_ocn_monthly_rgn.o%j
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
comp="ocn"
ckey="mpasocean"
group="ts"
grid="rgn"
freq="monthly"
hist="mpaso.hist.am.timeSeriesStatsMonthly"
yr_srt=1980
yr_end=1982
mth_srt=05
mth_end=04
ypf=$((yr_end-yr_srt+1))
ref_yr=1850

# number of complete 12-month periods
start_total=$(( yr_srt * 12 + mth_srt ))
end_total=$(( yr_end * 12 + mth_end ))
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


vars=\
timeMonthly_avg_mocStreamvalLatAndDepth,\
timeMonthly_avg_mocStreamvalLatAndDepthRegion,\
timeMonthly_avg_avgValueWithinOceanRegion_avgSurfaceTemperature,\
timeMonthly_avg_avgValueWithinOceanRegion_avgSurfaceSalinity,\
timeMonthly_avg_avgValueWithinOceanRegion_avgSurfaceSSH,\
timeMonthly_avg_avgValueWithinOceanRegion_avgSurfaceNetHeatFlux,\
timeMonthly_avg_avgValueWithinOceanRegion_avgSurfaceNetFreshWaterFlux,\
timeMonthly_avg_avgValueWithinOceanRegion_avgSurfaceWindStressMagnitude,\
timeMonthly_avg_avgValueWithinOceanRegion_avgWindStressZonal,\
timeMonthly_avg_avgValueWithinOceanRegion_avgWindStressMeridional,\
timeMonthly_avg_avgValueWithinOceanRegion_avgBoundaryLayerDepth,\
timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature,\
timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerSalinity,\
timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerThickness,\
timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerPotentialDensity,\
timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerVelocityZonal,\
timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerVelocityMeridional

vext=binBoundaryMocStreamfunction,xtime_startMonthly,xtime_endMonthly,time

if [[ -n "$ens" ]]; then
  input="${diri}/${case}/${ens}/archive/${comp}/hist"
else
  input="${diri}/${case_name}/archive/${comp}/hist"
fi

echo "Input directory: ${input}"

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

shopt -s nullglob
first_input_file=""

anchor_ymd="$(printf "%04d-%02d-01" "${sim_yr}" "${sim_mo}")"
for candidate in "${input}/${case_name}.${hist}.${anchor_ymd}.nc" "${input}/${case_name}.${hist}.${anchor_ymd}_"*.nc; do
  [[ -f "${candidate}" ]] || continue
  first_input_file="${candidate}"
  break
done

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
  cd ${scrip_dir}
  echo "Missing input files in ${input}"
  echo 'ERROR (1)' > ${group}_${comp}_${freq}_${grid}_${ttag}.status
  exit 1
fi

first_sim_day=$(read_first_nco_value timeMonthly_avg_daysSinceStartOfSim "${first_input_file}")

if ! [[ "${first_sim_day}" =~ ^-?[0-9]+([.][0-9]+)?$ ]]; then
  cd ${scrip_dir}
  echo "ERROR: Failed to read numeric timeMonthly_avg_daysSinceStartOfSim from ${first_input_file}"
  echo 'ERROR (3)' > ${group}_${comp}_${freq}_${grid}_${ttag}.status
  exit 3
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
      cd ${scrip_dir}
      echo "ERROR: Could not parse date from input filename: ${bname}"
      echo 'ERROR (1)' > ${group}_${comp}_${freq}_${grid}_${ttag}.status
      exit 1
    fi

    if [[ "${file_ymd}" < "${input_start_ymd}" || "${file_ymd}" > "${input_end_ymd}" ]]; then
      continue
    fi

    outfile="preproc/${bname}"
    tmp_with_time="preproc/tmp_with_time_${bname}"

    echo "Preprocessing ${bname}"

    run_nco ncks -O \
      -v "${vars},binBoundaryMocStreamfunction,xtime_startMonthly,xtime_endMonthly,timeMonthly_avg_daysSinceStartOfSim" \
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
  cd ${scrip_dir}
  echo 'Missing preprocessed input files'
  echo 'ERROR (1)' > ${group}_${comp}_${freq}_${grid}_${ttag}.status
  exit 1
fi

mkdir -p output

echo "ncclimo settings: vext=${vext}"

# Generate time series files
# If the user-defined parameter "vars" is "", then ${vars}, defined above, will be too.
cat input.txt | run_nco ncclimo \
-c "${case_name}" \
-v "${vars}" \
--split \
--yr_srt="${yr_srt}" \
--yr_end="${yr_end}" \
--mth_srt="${mth_srt}" \
--mth_end="${mth_end}" \
--ypf="${ypf}" \
--var_xtr="${vext}" \
-o output \
--prc_typ="${ckey}"

if [ $? != 0 ]; then
  cd ${scrip_dir}
  echo 'ERROR (2)' > ${group}_${comp}_${freq}_${grid}_${ttag}.status
  exit 2
fi

for file in output/*.nc; do
  [[ -f "${file}" ]] || continue
  run_nco ncatted -O -a _FillValue,time,d,, \
    "${file}" "${file}"
done

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
  echo 'ERROR (4)' > ${group}_${comp}_${freq}_${grid}_${ttag}.status
  exit 4
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
