#!/bin/bash

# Running on pm-cpu

#SBATCH  --job-name=ts_atm_daily_180x360_aave
#SBATCH  --account=e3sm
#SBATCH  --nodes=1
#SBATCH  --output=ts_atm_daily_180x360_aave.o%j
#SBATCH  --exclusive
#SBATCH  --time=0:30:00


#SBATCH  --qos=regular
#SBATCH  --constraint=cpu


######################
package="/share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh"
map="/global/cfs/cdirs/e3sm/diagnostics/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc"
diri="/pscratch/sd/z/zhan391/e3sm_project/E3SMv3_S2D"
diro="/global/cfs/cdirs/e3smdata/simulations/S2S2D"
case="WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_JRA55_FOSIRL_1980050100"
ens="EN00"
comp="atm"
ckey="eam"
group="ts"
grid="180x360_aave"
freq="daily"
hist="${ckey}.h1"
yr_srt=1980
yr_end=1982
mth_srt=05
mth_end=04
ypf=$((yr_end-yr_srt+1))
dpf=30
tpd=1

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

# Create symbolic links to input files
if [[ -n "$ens" ]]; then
  input="${diri}/${case}/${ens}/archive/${comp}/hist"
else
  input="${diri}/${case_name}/archive/${comp}/hist"
fi
for (( year=${yr_srt}; year<=${yr_end}; year++ ))
do
  YYYY=`printf "%04d" ${year}`
  for file in ${input}/${case_name}.${hist}.${YYYY}-*.nc
  do
    if [[ -f "${file}" ]]; then
      ln -s ${file} .
    fi 
  done
done

# For non-monthly input files, need to add the last file of the previous year
year=$((yr_srt - 1))
YYYY=`printf "%04d" ${year}`
mapfile -t files < <( ls ${input}/${case_name}.${hist}.${YYYY}-*.nc 2> /dev/null )
if [ ${#files[@]} -ne 0 ]
then
  ln -s ${files[-1]} .
fi
# as well as first file of next year to ensure that first and last years are complete
year=$((yr_end + 1))
YYYY=`printf "%04d" ${year}`
mapfile -t files < <( ls ${input}/${case_name}.${hist}.${YYYY}-*.nc 2> /dev/null )
if [ ${#files[@]} -ne 0 ]
then
  ln -s ${files[0]} .
fi


vars=FLUT,FLDS,FLNS,FSDS,FSNS,SHFLX,LHFLX,FSNT,FLNT,CLDTOT,TGCLDCWP,TGCLDIWP,TGCLDLWP,TAUX,TAUY,TH7001000,PSL,PS,U10,TS,TMQ,TUQ,TVQ,TOZ,PRECSC,PRECSL,PRECT,PRECTMX,PRECC,TREFHT,TREFHTMN,TREFHTMX,QREFHT,RHREFHT,UBOT,VBOT,TBOT,QBOT,ZBOT,UTOP,TTOP,V250,U250,U150,U025,U005,U002,U001,T250,T150,T025,T005,T002,T001,U1000,U975,U950,U925,U900,U850,U800,U700,U600,U500,U400,U300,U200,U100,U050,U010,V1000,V975,V950,V925,V900,V850,V800,V700,V600,V500,V400,V300,V200,V100,V050,V010,T1000,T975,T950,T925,T900,T850,T800,T700,T600,T500,T400,T300,T200,T100,T050,T010,Q1000,Q975,Q950,Q925,Q900,Q850,Q800,Q700,Q600,Q500,Q400,Q300,Q200,Q100,Q050,Q010,Z1000,Z975,Z950,Z925,Z900,Z850,Z800,Z700,Z600,Z500,Z400,Z300,Z200,Z100,Z050,Z010,OMEGA500,TROPF_P,TROPF_T,TROPF_Z,TCO,SCO

var_xtr=LANDFRAC,OCNFRAC,ICEFRAC

if [[ -n "${var_xtr}" ]]; then
  var_xtr_opt="--var_xtr=${var_xtr}"
fi

ls ${case_name}.${hist}.????-*.nc > input.txt
if grep -q "*" input.txt; then
  cd ${scrip_dir}
  echo 'Missing input files'
  echo 'ERROR (1)' > ${group}_${comp}_${freq}_${grid}_${ttag}.status
  exit 1
fi
# Generate time series files
# If the user-defined parameter "vars" is "", then ${vars}, defined above, will be too.
# Process complete calendar years first. This avoids an ncclimo 5.3.9 bug
# that constructs an invalid month 13 endpoint for partial ranges ending in December.
mth_srt1=1
mth_end1=12
ypf1=$((yr_end-yr_srt+1))
cat input.txt | run_nco ncclimo \
-c ${case_name} \
-v ${vars} \
--split \
${var_xtr_opt} \
--yr_srt=${yr_srt} \
--yr_end=${yr_end} \
--mth_srt=${mth_srt1} \
--mth_end=${mth_end1} \
--ypf=${ypf1} \
--map=${map} \
-o trash \
-O output1 \
--clm_md=hfs \
--dpf=${dpf} \
--tpd=${tpd} \
--prc_typ=${ckey}



if [ $? != 0 ]; then
  cd ${scrip_dir}
  echo 'ERROR (2)' > ${group}_${comp}_${freq}_${grid}_${ttag}.status
  exit 2
fi

# Trim the complete-year files to the requested hindcast range.
tag_srt=$(printf "%04d%02d" "${yr_srt}" "${mth_srt}")
tag_end=$(printf "%04d%02d" "${yr_end}" "${mth_end}")

srt_ym="${tag_srt:0:4}-${tag_srt:4:2}"
end_ym="${tag_end:0:4}-${tag_end:4:2}"

srt_date=$(date -d "${srt_ym}-01" +"%Y-%m-%dT00:00:00")
# Include midnight at the start of the following month for downstream boundaries.
end_date=$(date -d "${end_ym}-01 +1 month" +"%Y-%m-%dT00:00:00")

mkdir -p output

IFS=',' read -r -a var_list <<< "${vars}"

for var in "${var_list[@]}"; do
    var="${var//[[:space:]]/}"
    [[ -z "${var}" ]] && continue

    fout="${var}_${tag_srt}_${tag_end}.nc"

    shopt -s nullglob
    input_files=(output1/"${var}"_*.nc)
    shopt -u nullglob

    if (( ${#input_files[@]} == 0 )); then
        echo "WARNING: No input files found for variable ${var}" >&2
        continue
    fi

    echo "Trimming ${#input_files[@]} file(s) for ${var} to ${tag_srt}-${tag_end}..."
    if ! run_nco ncrcat -O \
        -d time,"${srt_date}","${end_date}" \
        "${input_files[@]}" \
        "output/${fout}"; then
      cd ${scrip_dir}
      echo 'ERROR (3)' > ${group}_${comp}_${freq}_${grid}_${ttag}.status
      exit 3
    fi
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
