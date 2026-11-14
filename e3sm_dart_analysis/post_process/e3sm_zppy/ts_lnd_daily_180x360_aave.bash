#!/bin/bash

# Running on pm-cpu

#SBATCH  --job-name=ts_lnd_daily_180x360_aave
#SBATCH  --account=e3sm
#SBATCH  --nodes=1
#SBATCH  --output=ts_lnd_daily_180x360_aave.o%j
#SBATCH  --exclusive
#SBATCH  --time=0:30:00
#SBATCH  --qos=regular
#SBATCH  --constraint=cpu


######################
package=/global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
map=/global/cfs/cdirs/e3sm/diagnostics/maps/map_r05_to_cmip6_180x360_aave.20231110.nc
diri=@ROOT@
diro=@ROOT@
case=@CASE@
ens=@ENS@
comp=lnd
ckey=elm
group="ts"
grid="180x360_aave"
freq="daily"
hist=elm.h1
yr_srt=@YR_SRT@
yr_end=@YR_END@
mth_srt=@MTH_SRT@
mth_end=@MTH_END@
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

# Full-length time series. These 2-3 month hindcasts are shorter than a year, and zppy's
# `ncclimo --split` then reads only the first few history files (atm stops after 4 days, lnd/rof at the
# end of the start month). Instead, concatenate every history file of the stream, remap it once
# (180x360 products), and write one file per variable named by the data's time range
# (<VAR>_<YYYYMM>_<YYYYMM>.nc), replacing earlier files of these variables.
status=${scrip_dir}/${group}_${comp}_${freq}_${grid}_${ttag}.status
# 180x360 daily products also record post/scripts/regen_ts_<stream>.status, which the hindcast
# notebooks check (fcst_init_experiment_config.stream_ready, hind_ic_test_eval)
regen_stream=lnd_daily
[[ -n "${regen_stream}" ]] && echo "RUNNING ${id}" > ${scrip_dir}/regen_ts_${regen_stream}.status
fail() {
  cd ${scrip_dir}
  echo "ERROR ($1)" > ${status}
  [[ -n "${regen_stream}" ]] && echo "ERROR $2" > ${scrip_dir}/regen_ts_${regen_stream}.status
  echo "ERROR: $2"
  exit $1
}

if [[ -n "$ens" ]]; then
  input="${diri}/${case}/${ens}/archive/${comp}/hist"
else
  input="${diri}/${case_name}/archive/${comp}/hist"
fi
mapfile -t files < <( ls ${input}/${case_name}.${hist}.*.nc 2> /dev/null )
(( ${#files[@]} )) || fail 1 "no ${hist} files in ${input}"

vars=H2OSNO,H2OSOI,SOILWATER_10CM,TWS,ZWT,TSOI_10CM,TSA
var_xtr=area,landfrac,landmask,topo
xtr_opt=${var_xtr:+,${var_xtr}}

# 1. concatenate the whole run (native grid)
run_nco ncrcat -O -7 -L 1 -v "${vars}${xtr_opt}" "${files[@]}" native.nc || fail 2 "ncrcat"
# 2. remap once
if [[ "${grid}" == "native" ]]; then
  mv native.nc rgr.nc
else
  prc=${ckey}; [[ "${ckey}" == "mosart" ]] && prc=""   # ncremap has no mosart type
  run_nco ncremap ${prc:+-P ${prc}} -t ${THREADS:-4} -m ${map} -i native.nc -o rgr.nc > ncremap.log 2>&1 \
    || { cat ncremap.log; fail 2 "ncremap"; }
fi
nrec=$(ncks --trd -M -m rgr.nc | awk '/^Record dimension: name = time, size =/ {print $NF; exit}')
[[ -z "${nrec}" ]] && nrec=$(ncks -C -v time --cdl -m rgr.nc | sed -n 's/.*time = UNLIMITED ; \/\/ (\([0-9]*\) currently).*/\1/p' | head -1)
tag=$(ncks -H -C -v time --cal rgr.nc | grep -oE '[0-9]{4}-[0-9]{2}' | sed -n '1p;$p' | tr -d '-' | paste -sd_)
[[ ${tag} =~ ^[0-9]{6}_[0-9]{6}$ ]] || fail 3 "cannot read time range"

# 3. one file per variable
mkdir -p output
IFS=',' read -ra var_list <<< "${vars}"
mapfile -t var_list < <( printf '%s\n' "${var_list[@]}" | awk 'NF && !seen[$0]++' )
for var in "${var_list[@]}"; do
  run_nco ncks -O -7 -L 1 -v "${var}${xtr_opt}" rgr.nc output/${var}_${tag}.nc || fail 3 "split ${var}"
done

# Move output ts files to final destination, replacing older files of these variables
# (zppy's truncated files may carry another month range)
if [[ -n "$ens" ]]; then
  dest="${diro}/${case}/${ens}/post/${comp}/${grid}/${group}/${freq}/${yint}yr"
else
  dest="${diro}/${case_name}/post/${comp}/${grid}/${group}/${freq}/${yint}yr"
fi
mkdir -p ${dest} || fail 4 "mkdir ${dest}"
for var in "${var_list[@]}"; do
  rm -f ${dest}/${var}_[0-9][0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][0-9][0-9].nc
done
mv output/*.nc ${dest} || fail 4 "move"
nout=$(ls ${dest}/*_${tag}.nc | wc -l)
(( nout == ${#var_list[@]} )) || fail 4 "wrote ${nout} of ${#var_list[@]} variables"
okmsg="files=${nout} records=${nrec} range=${tag} inputs=${#files[@]} $(date -Iseconds)"
[[ -n "${regen_stream}" ]] && echo "OK ${okmsg}" > ${scrip_dir}/regen_ts_${regen_stream}.status

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
echo "OK ${okmsg}" > ${group}_${comp}_${freq}_${grid}_${ttag}.status
exit 0
