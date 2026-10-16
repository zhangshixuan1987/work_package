#!/bin/bash

# Running on pm-cpu

#SBATCH  --job-name=ts_rof_daily
#SBATCH  --account=e3sm
#SBATCH  --nodes=1
#SBATCH  --output=ts_rof_daily.o%j
#SBATCH  --exclusive
#SBATCH  --time=0:30:00


#SBATCH  --qos=regular
#SBATCH  --constraint=cpu


######################
package="/share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh"
map="/global/cfs/cdirs/e3sm/diagnostics/maps/map_r05_to_cmip6_180x360_aave.20231110.nc"
diri="/pscratch/sd/z/zhan391/e3sm_project/E3SMv3_S2D"
diro="/global/cfs/cdirs/e3smdata/simulations/S2S2D"
case="WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_JRA55_FOSIRL_1980050100"
ens="EN00"
comp="lnd"
ckey="elm"
group="ts"
grid="native"
freq="daily"
hist="${ckey}.h1"
yr_srt=1980
yr_end=1982
mth_srt=05
mth_end=04
ypf=$((yr_end-yr_srt+1))
dpf=365
tpd=1

# number of complete 12-month periods
start_total=$(( yr_srt * 12 + mth_srt ))
end_total=$(( yr_end * 12 + mth_end ))
n_months=$(( end_total - start_total + 1 ))
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

# For non-daily input files, need to add the last file of the previous year
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

vars=H2OSNO,FSNO,QRUNOFF,QSNOMELT,FSNO_EFF,SNORDSL,SNOW,FSDS,FSR,FLDS,FIRE,FIRA
var_xtr=area,landfrac,landmask,topo

ls ${case_name}.${hist}.????-*.nc > input.txt
if grep -q "*" input.txt; then
  cd ${scrip_dir}
  echo 'Missing input files'
  echo 'ERROR (1)' > ${group}_${comp}_${freq}_${grid}_${ttag}.status
  exit 1
fi

# Generate time series files
# If the user-defined parameter "vars" is "", then ${vars}, defined above, will be too.
mth_srt1=1
mth_end1=12
ypf1=$((yr_end-yr_srt+1))
cat input.txt | run_nco ncclimo \
-c ${case_name} \
-v ${vars} \
--split \
--var_xtr=${var_xtr} \
--yr_srt=${yr_srt} \
--yr_end=${yr_end} \
--mth_srt=${mth_srt1} \
--mth_end=${mth_end1} \
--ypf=${ypf1} \
-o output1 \
--clm_md=hfs \
--dpf=${dpf} \
--tpd=${tpd} \
--prc_typ=${ckey}

if [ $? != 0 ]; then
  cd ${scrip_dir}
  echo 'ERROR (2)' > ${group}_${comp}_${freq}_${grid}_${ttag}.status
  exit 2
fi

# Trim to the expected hindcast range
mkdir -p output

tag_srt=$(printf "%04d%02d" ${yr_srt} ${mth_srt})
tag_end=$(printf "%04d%02d" ${yr_end} ${mth_end})
fout="${vars}_${tag_srt}_${tag_end}.nc"

srt_ym="${tag_srt:0:4}-${tag_srt:4:2}"
end_ym="${tag_end:0:4}-${tag_end:4:2}"
srt_date=$(date -d "${srt_ym}-01" +"%Y-%m-%dT00:00:00")
end_date=$(date -d "${end_ym}-01 +1 month -1 day" +"%Y-%m-%dT23:59:59")

run_nco ncrcat -O \
  -d time,"${srt_date}","${end_date}" \
  output1/${vars}_*.nc \
  "output/${fout}"

if [ $? != 0 ]; then
  cd ${scrip_dir}
  echo 'ERROR (3)' > ${group}_${comp}_${freq}_${grid}_${ttag}.status
  exit 3
fi

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
