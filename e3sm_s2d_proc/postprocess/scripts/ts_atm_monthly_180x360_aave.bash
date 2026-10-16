#!/bin/bash

# Running on pm-cpu

#SBATCH  --job-name=ts_atm_monthly_180x360_aave
#SBATCH  --account=e3sm
#SBATCH  --nodes=1
#SBATCH  --output=ts_atm_monthly_180x360_aave.o%j
#SBATCH  --exclusive
#SBATCH  --time=0:30:00


#SBATCH  --qos=regular
#SBATCH  --constraint=cpu


######################
package="/share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh"
plevdata="/global/cfs/cdirs/e3sm/diagnostics/e3sm_to_cmip_data/grids/vrt_remap_plev19.nc"
map="/global/cfs/cdirs/e3sm/diagnostics/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc"
diri="/pscratch/sd/z/zhan391/e3sm_project/E3SMv3_S2D"
diro="/global/cfs/cdirs/e3smdata/simulations/S2S2D"
case="WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_JRA55_FOSIRL_1980050100"
ens="EN00"
comp="atm"
ckey="eam"
group="ts"
grid="180x360_aave"
freq="monthly"
hist="${ckey}.h0"
yr_srt=1980
yr_end=1982
mth_srt=05
mth_end=04
ypf=$((yr_end-yr_srt+1))

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
# input data dir
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

vars=ICEFRAC,LANDFRAC,OCNFRAC,PSL,FSNTC,FSNTOAC,SWCF,LWCF,FLUT,FSNT,FSNTOA,FLNT,FLNTC,FSNS,FLNS,SHFLX,QFLX,LHFLX,TAUX,TAUY,PRECC,PRECL,PRECSC,PRECSL,TS,TREFHT,U10,QREFHT,TMQ,CLDTOT,CLDHGH,CLDMED,CLDLOW,FLDS,FSDS,TGCLDIWP,TGCLDCWP,TGCLDLWP,FLNSC,FLUTC,FSDSC,SOLIN,FSNSC,AODBC,AODDUST,AODSO4,AODSOA,AODSS,AODVIS,TREFMNAV,TREFMXAV,TCO,PS,PHIS,AREL,U,V,T,Q,Z3,OMEGA,RELHUM
#TROP_P,CO_SRF,AODPOM,AODALL

interp_vars=AREL,U,V,T,Q,Z3,OMEGA,RELHUM
derive_vars=(PRECT  PRECST)
derive_var1=(PRECL  PRECSL)
derive_var2=(PRECC  PRECSC)

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
cat input.txt | run_nco ncclimo \
-c ${case_name} \
-v ${vars} \
--split \
${var_xtr_opt} \
--yr_srt=${yr_srt} \
--yr_end=${yr_end} \
--mth_srt=${mth_srt} \
--mth_end=${mth_end} \
--ypf=${ypf} \
--map=${map} \
-o trash \
-O output \
--prc_typ=${ckey}

if [ $? != 0 ]; then
  cd ${scrip_dir}
  echo 'ERROR (2)' > ${group}_${comp}_${freq}_${grid}_${ttag}.status
  exit 2
fi

# process 3D model level to pressure level data
y1=$(printf "%04d" "${yr_srt}")
y2=$(printf "%04d" "${yr_end}")
IFS=',' read -ra mlvars <<< "${interp_vars}"
for var in "${mlvars[@]}"; do
  for file in output/${var}_${y1}??_${y2}??.nc; do
    [[ -e "$file" ]] || continue
    tmp_file="${file}.plev"
    run_nco ncremap -p mpi --vrt_ntp=log --vrt_xtr=mss_val --vrt_out="${plevdata}" -i "${file}" -o "${tmp_file}"
    mv "${tmp_file}" "${file}"
  done
done

# process derived variables
derive_err=0
for i in "${!derive_vars[@]}"; do
  target="${derive_vars[$i]}"
  var1="${derive_var1[$i]}"
  var2="${derive_var2[$i]}"
    
  for file1 in output/${var1}_${y1}??_${y2}??.nc; do
    [[ -e "$file1" ]] || continue

    suffix="${file1#output/${var1}_}"
    file2="output/${var2}_${suffix}"
    outfile="output/${target}_${suffix}"
    tmpfile="output/tmp_${target}_${suffix}"

    # skip if already exists
    [[ -e "$outfile" ]] && continue

    # skip if second file missing
    [[ -e "$file2" ]] || continue

    echo "Deriving $target for $suffix"

    cp "$file1" "$tmpfile" || { derive_err=1; rm -f "$tmpfile"; continue; }
    run_nco ncks -A "$file2" "$tmpfile" || { derive_err=1; rm -f "$tmpfile"; continue; }
    echo $target $var1 $var2
    run_nco ncap2 -O -s "${target}=${var1}+${var2}" "$tmpfile" "$outfile" || { derive_err=1; rm -f "$tmpfile"; continue; }
    rm -f "$tmpfile"
  done
done

if [ "${derive_err}" == "1" ]; then
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
