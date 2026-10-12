#!/bin/bash -el
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --account=esmd
#SBATCH  --time=2:00:00
#SBATCH  --partition=short
#SBATCH  --job-name=eam_darts1_compare
#SBATCH  --nodes=1
#SBATCH  --output=eam_darts1_compare.%j
#SBATCH  --exclusive
#SBATCH  --no-kill
#SBATCH  --requeue

set -Eeuo pipefail
shopt -s nullglob

fail() { echo "ERROR: $*" >&2; exit 1; }
trap 'echo "ERROR: comparison failed at line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

echo == Start of e3sm dart diagnostic ==
date
echo ============================================

#For cshell:
#limit stacksize unlimited
#limit datasize unlimited

#For bash 
#ulimit -s unlimited
#ulimit -d unlimited

#export SLURM_NNODES=20
#export SLURM_NTASKS=800

VERBOSE='-v'
MOVE='/usr/bin/mv'
COPY='/usr/bin/cp --preserve=timestamps'
LINK='/usr/bin/ln -fs'
LINKV=TRUE
LIST='/usr/bin/ls'
REMOVE='/usr/bin/rm'
LAUNCHCMD='srun -N 1 -n 24'
OBSDIAG_LAUNCHCMD='srun -N 1 -n 1'

# Select a resumable workflow stage: all, common, obs2nc, or obs_diag.
COMPARE_STAGE="obs_diag"
case "${COMPARE_STAGE}" in
  all|common|obs2nc|obs_diag) ;;
  *) fail "invalid COMPARE_STAGE=${COMPARE_STAGE}; use all, common, obs2nc, or obs_diag" ;;
esac

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  my_wkdir=$(readlink -f "${SLURM_SUBMIT_DIR}") || fail "cannot resolve SLURM_SUBMIT_DIR"
else
  SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || fail "cannot resolve script path"
  my_wkdir=$(dirname "${SCRIPT_PATH}")
fi
scomp="eam"
cd "${my_wkdir}"

source ./create_and_setup_case.sh

# This driver compares the May 2012 CAM6 assimilation with DARTEN40S1.
# Both endpoints are inclusive and span 124 six-hour observation sequences.
DART_DIAG_YMDS="2012-05-01-00000"
DART_DIAG_YMDE="2012-05-31-64800"
# Use the final per-experiment cam6_common directory directly as the output workspace.
COMPARE_WORKDIR="${my_wkdir}/runtmp/compare_2012_cam6_darten40s1"
mkdir -p "${COMPARE_WORKDIR}"
exec 8>"${my_wkdir}/runtmp/compare_2012_cam6_darten40s1.lock"
flock -n 8 || fail "another CAM6-DARTEN40S1 comparison is active"

# Set paths
DART_ROOT=${my_dart_code}
DART_MODEL=${my_dart_eam}
DART_SCPTDIR=${DART_ROOT}/models/${DART_MODEL}/shell_scripts
DART_WORKDIR=${DART_ROOT}/models/${DART_MODEL}/work
BASE_OBSDIR=${my_dart_obsdir}
BASE_PHIS=${my_e3sm_topo}
BASE_SEMAPS=${my_e3sm_semap}
BASE_CSGRID=${my_e3sm_csgrid}
# Explicit two-case comparison. CAM6 has an archived DART tree but no workflow
# configuration here, so experiment 2 supplies the compatible environment and
# diagnostics namelist for both inputs.
dart_experiment_root="/qfs/people/zhan391/e3sm_dart_work/run_script/v3_dart_cda/2_ne30pg2_dart_atm_en40"
dart_config="${dart_experiment_root}/create_and_setup_case.sh"
[[ -r "${dart_config}" ]] || fail "missing DARTEN40S1 configuration: ${dart_config}"
dart_metadata=$(bash -c 'source "$1" >/dev/null 2>&1; printf "%s|%s|%s|%s|%s" "$my_casename" "$my_ensnum" "$my_eam_dart_run_dir" "$my_dart_env_file" "$my_eam_diag_nml"' _ "${dart_config}")
IFS='|' read -r dart_case dart_ens dart_rundir dart_env_file dart_diag_nml <<< "${dart_metadata}"

cam6_case="f.e22.FHIST_BGC.f09_025.CAM6assim.011"
cam6_rundir="/compyfs/zhan391/v3_dart_cda_scratch/${cam6_case}/archive/dart_en80"
[[ "${dart_case}" == "DARTEN40S1_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy" ]] ||
  fail "unexpected experiment-4 case name: ${dart_case}"
[[ "${dart_ens}" == "40" && -d "${dart_rundir}" && -d "${cam6_rundir}" &&
   -r "${dart_env_file}" && -r "${dart_diag_nml}" ]] ||
  fail "invalid CAM6-DARTEN40S1 comparison metadata"

declare -a refcases=("${cam6_case}" "${dart_case}")
declare -a refnens=(80 "${dart_ens}")
declare -a rundirs=("${cam6_rundir}" "${dart_rundir}")
declare -a dart_env_files=("${dart_env_file}" "${dart_env_file}")
declare -a diag_nml_files=("${dart_diag_nml}" "${dart_diag_nml}")
num_to_compare_at_once=${#rundirs[@]}
echo "Using configured DART machine environment: ${dart_env_files[0]}"
source "${dart_env_files[0]}"
DART_DIAG_NML="${diag_nml_files[0]}"

# Retain the final experiment as the template/grid owner for legacy code below.
DART_ENSNUM=${refnens[num_to_compare_at_once-1]}
DART_CASE=${refcases[num_to_compare_at_once-1]}
DART_RUNDIR="${rundirs[num_to_compare_at_once-1]}/dart_diagnostics/cam6_common"
mkdir -p "${DART_RUNDIR}"
DART_NTASKS=${my_dart_ntask}
DART_ON_PGRID=${my_dart_pgrid}

DATA_ASSIMILATION_CYCLES=124
DATA_ASSIMILATION_WINDOW=${my_dart_window}
DATA_ASSIMILATION_ATM=TRUE

homme_map_file="SEMapping.nc"
cs_grid_file="SEMapping_cs_grid.nc"

CUR_YMD="${DART_DIAG_YMDS%-*}"
CUR_TOD="${DART_DIAG_YMDS##*-}"
CUR_DATE=( `echo ${CUR_YMD}-${CUR_TOD} | sed -e "s#-# #g"` )
CUR_YEAR=`echo "${CUR_DATE[0]}" | bc`
CUR_MONTH=`echo "${CUR_DATE[1]}" | bc`
CUR_DAY=`echo "${CUR_DATE[2]}" | bc`
CUR_SECONDS=`echo "${CUR_DATE[3]}" | bc`
CUR_HOUR=`echo "${CUR_DATE[3]}" / 3600 | bc`
echo "valid time for eam forecast cycle is $CUR_YEAR $CUR_MONTH $CUR_DAY $CUR_SECONDS (seconds)"
echo "valid time for eam forecast cycle is $CUR_YEAR $CUR_MONTH $CUR_DAY $CUR_HOUR (hours)"

CURRENT_DADIR="${DART_RUNDIR}"
if [ ! -d ${CURRENT_DADIR} ]; then
  mkdir -p ${CURRENT_DADIR}
fi
cd ${CURRENT_DADIR}

user_dart_nl() {
  if [[ -r "${DART_DIAG_NML}" ]]; then
    ${2} "${DART_DIAG_NML}" input.nml  || exit 123
    sed -i "/#/d;/^\!/d;/^[ ]*\!/d" input.nml
    sed -i '1,1i\WARNING: Changes to this file will be ignored. \n Edit \$DART_WORKDIR/diag_dart_input.nml instead.\n\n\n'  input.nml
  else
    echo "ERROR ... configured DART diagnostics namelist ${DART_DIAG_NML} not found ... ERROR"
    exit 124
  fi

  xlist=`grep '^[ ]*vertical_localization_coord' input.nml`
  xlist=( `echo $xlist | sed -e "s#[=,']# #g"` )
  if [ "${xlist[1]}" == "SCALEHEIGHT" ]; then
     list1=`grep '^[ ]*vert_normalization_scale_height' input.nml `
     list1=( `echo $list1 | sed -e "s#[=,]##g"` )
     if [ "${list1[1]}" != "1.5" ]; then
        echo "WARNING!  input.nml is not using 1.5 for vert_normalization_scale_height."
        echo "          Use a different value only if you definitely want to. "
     fi
  else
     echo "WARNING!  input.nml is not using SCALEHEIGHT for vertical_localization_coord."
     echo "          SCALEHEIGHT is highly recommended for EAM"
  fi
  
  # If possible, use the round-robin approach to deal out the tasks.
  # This facilitates using multiple nodes for the simultaneous I/O operations.
  if [ -v ${3} ]; then
     if [ ${#3[@]} -gt 0 ]; then
        sed -i "s#layout.*#layout = 2#"  input.nml
        sed -i "s#tasks_per_node.*#tasks_per_node = ${3}#" input.nml
     fi
  fi
  
}

user_obs_diag_nl() {
  local ens_size="$1"
  local use_pgrid="$2"
  local obs_list="$3"
  local bin_start="$4"
  local bin_end="$5"
  local bin_sep="$6"
  local bin_width="$7"
  local skip_time="$8"
  local trusted="$9"

  [[ -f input.nml ]] || { echo "input.nml not found"; return 1; }

  cat <<EOF | ex -s input.nml
g#^ *ens_size *=#s#=.*#= ${ens_size}#
g#^ *num_output_state_members *=#s#=.*#= ${ens_size}#
g#^ *num_output_obs_members *=#s#=.*#= ${ens_size}#
g#^ *eam_use_pgrid *=#s#=.*#= ${use_pgrid}#
g#^ *obs_sequence_name *=#s#=.*#= ''#
g#^ *obs_sequence_list *=#s#=.*#= '${obs_list}'#
g#^ *first_bin_center *=#s#=.*#= ${bin_start}#
g#^ *last_bin_center *=#s#=.*#= ${bin_end}#
g#^ *bin_separation *=#s#=.*#= ${bin_sep}#
g#^ *bin_width *=#s#=.*#= ${bin_width}#
g#^ *time_to_skip *=#s#=.*#= ${skip_time}#
g#^ *trusted_obs *=#s#=.*#= ${trusted}#
wq
EOF
}

user_obs2nc_nl() {
  local ens_size="$1"
  local use_pgrid="$2"
  local obs_list="$3"
  local first_bin_start="$4"
  local first_bin_end="$5"
  local last_bin_end="$6"
  local bin_days="$7"
  local bin_seconds="$8"

  cat <<EOF | ex -s input.nml
g#^ *ens_size *=#s#=.*#= ${ens_size}#
g#^ *num_output_state_members *=#s#=.*#= ${ens_size}#
g#^ *num_output_obs_members *=#s#=.*#= ${ens_size}#
g#^ *eam_use_pgrid *=#s#=.*#= ${use_pgrid}#
g#^ *obs_sequence_name *=#s#=.*#= ''#
g#^ *obs_sequence_list *=#s#=.*#= '${obs_list}'#
g#^ *first_bin_start *=#s#=.*#= ${first_bin_start}#
g#^ *first_bin_end *=#s#=.*#= ${first_bin_end}#
g#^ *last_bin_end *=#s#=.*#= ${last_bin_end}#
g#^ *bin_interval_days *=#s#=.*#= ${bin_days}#
g#^ *bin_interval_seconds *=#s#=.*#= ${bin_seconds}#
wq
EOF
}

user_obs_sequence_nl() {
  local filename_seq_list="$1"
  local filename_out="$2"
  local qc_metadata="$3"
  local min_qc="$4"
  local max_qc="$5"
  local calendar="$6"
  local print_only="$7"
  local edit_copies="$8"
  local new_copy_index="$9"

  [[ -f input.nml ]] || { echo "input.nml not found"; return 1; }

  cat <<EOF | ex -s input.nml
g#^ *filename_seq *=#s#=.*#= '',#
g#^ *filename_seq_list *=#s#=.*#= ${filename_seq_list},#
g#^ *filename_out *=#s#=.*#= ${filename_out},#
g#^ *qc_metadata *=#s#=.*#= ${qc_metadata},#
g#^ *min_qc *=#s#=.*#= ${min_qc},#
g#^ *max_qc *=#s#=.*#= ${max_qc},#
g#^ *calendar *=#s#=.*#= '${calendar}',#
g#^ *print_only *=#s#=.*#= ${print_only},#
g#^ *edit_copies *=#s#=.*#= ${edit_copies},#
g#^ *new_copy_index *=#s#=.*#= ${new_copy_index},#
wq
EOF
}

user_obs_common_nl() {
  local num_to_compare_at_once="$1"
  local filename_seq_list="$2"
  local filename_out_suffix="$3"
  local print_every="$4"
  local dart_qc_threshold="$5"
  local calendar="$6"
  local print_only="$7"
  local eval_and_assim_can_match="$8"

  [[ -f input.nml ]] || { echo "input.nml not found"; return 1; }

  cat <<EOF | ex -s input.nml
g#^ *num_to_compare_at_once *=#s#=.*#= ${num_to_compare_at_once},#
g#^ *filename_seq *=#s#=.*#= '',#
g#^ *filename_seq_list *=#s#=.*#= ${filename_seq_list},#
g#^ *filename_out_suffix *=#s#=.*#= ${filename_out_suffix},#
g#^ *print_every *=#s#=.*#= ${print_every},#
g#^ *dart_qc_threshold *=#s#=.*#= ${dart_qc_threshold},#
g#^ *calendar *=#s#=.*#= '${calendar}',#
g#^ *print_only *=#s#=.*#= ${print_only},#
g#^ *eval_and_assim_can_match *=#s#=.*#= ${eval_and_assim_can_match},#
wq
EOF
}

cd ${CURRENT_DADIR}

stage_common_files() {
  local target_dir="$1"
  local stage_dir="$2"
  mkdir -p "${target_dir}"
  mkdir -p "${stage_dir}"
  find "${stage_dir}" -mindepth 1 -maxdepth 1 -exec rm -rf {} +
  find "${target_dir}" -mindepth 1 -maxdepth 1 -exec rm -rf {} +
}

############################################################
# run obs_diag (observation space diagnostics) 
############################################################
OBSCOMMON_DIR="${CURRENT_DADIR}"
mkdir -p "${OBSCOMMON_DIR}"
cd ${OBSCOMMON_DIR}

user_dart_nl ${DART_WORKDIR} ${COPY} ${my_task_per_node}

#delete line that is deprecated by current version of DART
MYSTRING=`grep eam_template_filename input.nml`
MYSTRING=( `echo $MYSTRING | sed -e "s#[=,']# #g"` )
EAMINPUT=${MYSTRING[1]}

MYSTRING=`grep eam_phis_filename input.nml`
MYSTRING=( `echo $MYSTRING | sed -e "s#[=,']# #g"` )
EAM_PHIS=${MYSTRING[1]}
${LINK} ${BASE_PHIS} ${EAM_PHIS} || exit 99

#Now, Link the grid information files 
if [ ! -f ${BASE_SEMAPS} ] && [ ! -f ${BASE_CSGRID} ]; then
  echo "ERROR ... no mapping file ${homme_map_file}"
  echo "ERROR ... no gridinfo file ${cs_grid_file}"
  echo "ERROR ... must provide either of them"
  exit 100
else
  if [ -f ${BASE_SEMAPS} ]; then
    ${COPY} -r ${BASE_SEMAPS} ${homme_map_file} || exit 101
  fi
  if [ -f ${BASE_CSGRID} ]; then
    ${COPY} -r ${BASE_CSGRID} ${cs_grid_file} || exit 102
  fi
fi

[[ -f input.nml ]] || { echo "input.nml not found"; return 1; }

#use obs_common_subset to generate diagnostics 
${COPY} -f ${DART_WORKDIR}/obs_common_subset  ${OBSCOMMON_DIR} || exit 103
declare -a outlists=()
filename_seq_list=""
for ((i=0; i<num_to_compare_at_once; i++)); do
  outlists+=("${OBSCOMMON_DIR}/obs_seq${i}")
  [[ -z "${filename_seq_list}" ]] || filename_seq_list+=","
  filename_seq_list+="'obs_seq${i}'"
done
filename_out_suffix="'.common'"
print_every=10000
dart_qc_threshold=8
calendar='Gregorian'
print_only=".false."
eval_and_assim_can_match=".true."
user_obs_common_nl \
  "${num_to_compare_at_once}" \
  "${filename_seq_list}" \
  "${filename_out_suffix}" \
  "${print_every}" \
  "${dart_qc_threshold}" \
  "${calendar}" \
  "${print_only}" \
  "${eval_and_assim_can_match}"

# Build the exact shared six-hour timestamp sequence.
stamp_to_epoch() {
  local stamp="$1" ymd tod hour
  ymd="${stamp%-*}"
  tod="${stamp##*-}"
  hour=$((10#${tod}/3600))
  date -u -d "${ymd} $(printf '%02d:00:00' "${hour}")" +%s
}
epoch_to_stamp() {
  local epoch="$1" ymd hour
  read -r ymd hour < <(date -u -d "@${epoch}" '+%Y-%m-%d %H')
  printf '%s-%05d\n' "${ymd}" "$((10#${hour}*3600))"
}
start_epoch=$(stamp_to_epoch "${DART_DIAG_YMDS}") || fail "invalid comparison start: ${DART_DIAG_YMDS}"
end_epoch=$(stamp_to_epoch "${DART_DIAG_YMDE}") || fail "invalid comparison end: ${DART_DIAG_YMDE}"
step_seconds=$((DATA_ASSIMILATION_WINDOW*3600))
(( start_epoch <= end_epoch && (end_epoch-start_epoch)%step_seconds == 0 )) ||
  fail "comparison range is invalid or not aligned to the DA interval"
declare -a timestamps=()
for ((epoch=start_epoch; epoch<=end_epoch; epoch+=step_seconds)); do
  timestamps+=("$(epoch_to_stamp "${epoch}")")
done

# Preflight and write deterministic, identically aligned input lists.
for (( i=0; i< num_to_compare_at_once; i++ )); do
    rundir="${rundirs[$i]}"
    refcase="${refcases[$i]}"
    outlist="${outlists[$i]}"
    src_dir="${rundir}/dart_diagnostics/obs_seq"
    prefix="${refcase}.dart.e.${scomp}_obs_seq.processed"
  
    [[ -d "${src_dir}" ]] || fail "missing processed-observation directory: ${src_dir}"
    : > "${outlist}"
    for ts in "${timestamps[@]}"; do
      f="${src_dir}/${prefix}.${ts}"
      [[ -s "${f}" ]] || fail "missing or empty processed observation sequence: ${f}"
      printf '%s\n' "${f}" >> "${outlist}"
    done
    [[ $(wc -l < "${outlist}") -eq ${#timestamps[@]} ]] ||
      fail "incorrect input count in ${outlist}"
    echo "[OK] Preflighted ${#timestamps[@]} inputs: ${outlist}"
done

if [[ "${COMPARE_PREFLIGHT_ONLY:-FALSE}" == "TRUE" ]]; then
  echo "Comparison preflight completed successfully for ${num_to_compare_at_once} experiments and ${#timestamps[@]} timestamps"
  exit 0
fi

# Remove stale .common files from any legacy obs_seq directories before the next run.
for (( i=0; i<num_to_compare_at_once; i++ )); do
  legacy_dir="${rundirs[$i]}/dart_diagnostics/obs_seq"
  if [[ -d "${legacy_dir}" ]]; then
    find "${legacy_dir}" -maxdepth 1 -type f -name '*.common' -delete
  fi
done

common_outputs_complete() {
  local i input_file common_file tgt_dir
  for (( i=0; i<num_to_compare_at_once; i++ )); do
    tgt_dir="${rundirs[$i]}/dart_diagnostics/cam6_common"
    [[ -d "${tgt_dir}" ]] || return 1
    while IFS= read -r input_file; do
      common_file="${tgt_dir}/$(basename "${input_file}").common"
      [[ -s "${common_file}" ]] || return 1
    done < "${outlists[$i]}"
  done
}

run_common=FALSE
if [[ "${my_dartdiag_common}" == ".true." &&
      ( "${COMPARE_STAGE}" == "all" || "${COMPARE_STAGE}" == "common" ) ]]; then
  run_common=TRUE
  if [[ "${REUSE_COMMON_OUTPUTS:-TRUE}" == "TRUE" ]] && common_outputs_complete; then
    echo "[OK] Reusing complete .common outputs for all ${num_to_compare_at_once} experiments"
    run_common=FALSE
  fi
fi

if [[ "${COMPARE_STAGE}" == "obs2nc" || "${COMPARE_STAGE}" == "obs_diag" ]]; then
  common_outputs_complete ||
    fail "COMPARE_STAGE=${COMPARE_STAGE} requires a complete existing set of .common files"
  echo "[OK] Validated existing .common outputs; skipping obs_common_subset"
fi

if [[ "${run_common}" == "TRUE" ]]; then
  # obs_common_subset can terminate worker ranks after their assigned files are
  # complete. Slurm reports those exits as a nonzero step status even when all
  # requested output was written. Let the exhaustive checks below decide.
  if ! ${LAUNCHCMD} "${OBSCOMMON_DIR}/obs_common_subset"; then
    echo "[WARN] obs_common_subset launcher returned nonzero; validating all expected outputs"
  fi
fi

if [[ "${run_common}" == "TRUE" ]]; then
  # Directly generate each set of common files in the final cam6_common output directory.
  for (( i=0; i<num_to_compare_at_once; i++ )); do
    tgt_dir="${rundirs[$i]}/dart_diagnostics/cam6_common"
    mkdir -p "${tgt_dir}"
    find "${tgt_dir}" -maxdepth 1 -type f -name '*.common' -delete
    moved=0
    while IFS= read -r input_file; do
      common_file="${input_file}.common"
      [[ -s "${common_file}" ]] || fail "obs_common_subset did not create ${common_file}"
      ${MOVE} "${common_file}" "${tgt_dir}/"
      moved=$((moved+1))
    done < "${outlists[$i]}"
    (( moved == ${#timestamps[@]} )) || fail "generated only ${moved} common files for ${refcases[$i]} in ${tgt_dir}"
    printf '[OK] Generated %d common files for %s in %s\n' "${moved}" "${refcases[$i]}" "${tgt_dir}"
  done
fi

#now process the obs diagnostics for obs_common
if [[ "${my_dartdiag_obs2nc}" == ".true." &&
      ( "${COMPARE_STAGE}" == "all" || "${COMPARE_STAGE}" == "obs2nc" ) ]]; then

  # normalize suffix once (handles the case filename_out_suffix="'.common'")
  filename_out_suffix="${filename_out_suffix//\'/}"

  for (( i=0; i<num_to_compare_at_once; i++ )); do
    rundir="${rundirs[$i]}"
    refnen="${refnens[$i]}"        # <- use a new var; don't overwrite the array name
    refcase="${refcases[$i]}"
    outlist="${outlists[$i]}"
    prefix="${refcase}.dart.e.${scomp}_obs_seq.processed"
    OBSSEQ_DIR="${rundir}/dart_diagnostics/cam6_common"

    if [[ ! -d "${OBSSEQ_DIR}" ]]; then
       echo "ERROR: ${OBSSEQ_DIR} does not exist..."
       exit 105
    fi

    cd "${OBSSEQ_DIR}"
    ${COPY} -f "${DART_WORKDIR}/obs_seq_to_netcdf" "${OBSSEQ_DIR}" || exit 106
    user_dart_nl "${DART_WORKDIR}" "${COPY}" "${my_task_per_node}"

    # delete line that is deprecated by current version of DART
    MYSTRING="$(grep eam_template_filename input.nml || true)"
    MYSTRING=( $(echo "$MYSTRING" | sed -e "s#[=,']# #g") )
    EAMINPUT="${MYSTRING[1]}"

    MYSTRING="$(grep eam_phis_filename input.nml || true)"
    MYSTRING=( $(echo "$MYSTRING" | sed -e "s#[=,']# #g") )
    EAM_PHIS="${MYSTRING[1]}"
    ${LINK} "${BASE_PHIS}" "${EAM_PHIS}" || exit 107

    # Now, Link the grid information files
    if [[ ! -f ${BASE_SEMAPS} && ! -f ${BASE_CSGRID} ]]; then
      echo "ERROR ... no mapping file ${homme_map_file}"
      echo "ERROR ... no gridinfo file ${cs_grid_file}"
      echo "ERROR ... must provide either of them"
      exit 108
    else
      if [[ -f ${BASE_SEMAPS} ]]; then
        ${COPY} -r "${BASE_SEMAPS}" "${homme_map_file}" || exit 109
      fi
      if [[ -f ${BASE_CSGRID} ]]; then
        ${COPY} -r "${BASE_CSGRID}" "${cs_grid_file}" || exit 110
      fi
    fi

    [[ -f input.nml ]] || { echo "input.nml not found"; exit 111; }

    obs_sequence_list="${OBSSEQ_DIR}/obs_sequence_list.txt"
    ${LIST} "${OBSSEQ_DIR}/${prefix}."*"${filename_out_suffix}" > "${obs_sequence_list}"
    cat "${obs_sequence_list}"

    # Extract timestamps like YYYY-MM-DD-SSSSS (unique + sorted)
    mapfile -t timestamps < <(grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2}-[0-9]{5}' "${obs_sequence_list}" | sort -u)

    if (( ${#timestamps[@]} == 0 )); then
        echo "ERROR: No timestamps found in ${obs_sequence_list}"
        exit 112
    fi

    first_ts="${timestamps[0]}"
    last_ts="${timestamps[$((${#timestamps[@]}-1))]}"

    # Seconds per DA window (e.g., 6h => 21600)
    DART_SECINC=$((3600 * DATA_ASSIMILATION_WINDOW))

    # Split components for first timestamp
    IFS='-' read -r DART_YEAR1 DART_MONTH1 DART_DAY1 DART_SEC1 <<< "${first_ts}"
    DART_HOUR1=$((10#${DART_SEC1} / 3600))                 # 21600 → 6
    DART_HOUR1e=$((DART_HOUR1 + DART_SECINC / 3600))       # window end hour (same-day assumption)

    # Split components for last timestamp
    IFS='-' read -r DART_YEAR2 DART_MONTH2 DART_DAY2 DART_SEC2 <<< "${last_ts}"
    DART_HOUR2=$((10#${DART_SEC2} / 3600))

    echo "[INFO] First timestamp: ${first_ts} → ${DART_YEAR1}-${DART_MONTH1}-${DART_DAY1} hour=${DART_HOUR1}"
    echo "[INFO] Last  timestamp: ${last_ts} → ${DART_YEAR2}-${DART_MONTH2}-${DART_DAY2} hour=${DART_HOUR2}"

    # Assemble bin times
    first_bin_start="${DART_YEAR1},${DART_MONTH1},${DART_DAY1},${DART_HOUR1},0,0"
    first_bin_end="${DART_YEAR1},${DART_MONTH1},${DART_DAY1},${DART_HOUR1e},0,0"
    last_bin_end="${DART_YEAR2},${DART_MONTH2},${DART_DAY2},${DART_HOUR2},0,0"
    bin_interval_days="0"
    bin_interval_seconds="${DART_SECINC}"

    echo "${first_bin_start} ${first_bin_end} ${last_bin_end}"

    user_obs2nc_nl "${refnen}" "${DART_ON_PGRID}" "${obs_sequence_list}" \
                   "${first_bin_start}" "${first_bin_end}" "${last_bin_end}" \
                   "${bin_interval_days}" "${bin_interval_seconds}"

    # convert obs sequence to netcdf file
    ${LAUNCHCMD} "${OBSSEQ_DIR}/obs_seq_to_netcdf" || exit 113

    # rename the files (use j to avoid clobbering i)
    for (( j=0; j<${#timestamps[@]}; j++ )); do
        fepoch="${OBSSEQ_DIR}/obs_epoch_$(printf "%03d" ${j}).nc"
        OUT_DATE="${timestamps[$j]}"
        if [[ -f "${fepoch}" ]]; then
            ${MOVE} "${fepoch}" "${OBSSEQ_DIR}/${refcase}.dart.e.${scomp}_obs_seq_common.${OUT_DATE}.nc"
        else
            echo "[WARN] Missing file: ${fepoch}"
        fi
    done
  done
fi

if [[ "${my_dartdiag_obdiag}" == ".true." &&
      ( "${COMPARE_STAGE}" == "all" || "${COMPARE_STAGE}" == "obs_diag" ) ]]; then

  # normalize suffix once (handles the case filename_out_suffix="'.common'")
  filename_out_suffix="${filename_out_suffix//\'/}"

  for (( i=0; i<num_to_compare_at_once; i++ )); do
    rundir="${rundirs[$i]}"
    refnen="${refnens[$i]}"
    refcase="${refcases[$i]}"
    outlist="${outlists[$i]}"
    prefix="${refcase}.dart.e.${scomp}_obs_seq.processed"
    OBSDIAG_DIR="${rundir}/dart_diagnostics/cam6_common"

    if [[ ! -d "${OBSDIAG_DIR}" ]]; then
       echo "ERROR: ${OBSDIAG_DIR} does not exist..."
       exit 114
    fi

    cd "${OBSDIAG_DIR}"
    ${COPY} -f "${DART_WORKDIR}/obs_diag" "${OBSDIAG_DIR}" || exit 115
    user_dart_nl "${DART_WORKDIR}" "${COPY}" "${my_task_per_node}"

    # delete line that is deprecated by current version of DART
    MYSTRING="$(grep eam_template_filename input.nml || true)"
    MYSTRING=( $(echo "$MYSTRING" | sed -e "s#[=,']# #g") )
    EAMINPUT="${MYSTRING[1]}"

    MYSTRING="$(grep eam_phis_filename input.nml || true)"
    MYSTRING=( $(echo "$MYSTRING" | sed -e "s#[=,']# #g") )
    EAM_PHIS="${MYSTRING[1]}"
    ${LINK} "${BASE_PHIS}" "${EAM_PHIS}" || exit 116

    # Now, Link the grid information files
    if [[ ! -f ${BASE_SEMAPS} && ! -f ${BASE_CSGRID} ]]; then
      echo "ERROR ... no mapping file ${homme_map_file}"
      echo "ERROR ... no gridinfo file ${cs_grid_file}"
      echo "ERROR ... must provide either of them"
      exit 117
    else
      if [[ -f ${BASE_SEMAPS} ]]; then
        ${COPY} -r "${BASE_SEMAPS}" "${homme_map_file}" || exit 118
      fi
      if [[ -f ${BASE_CSGRID} ]]; then
        ${COPY} -r "${BASE_CSGRID}" "${cs_grid_file}" || exit 119
      fi
    fi

    [[ -f input.nml ]] || { echo "input.nml not found"; exit 120; }

    obs_sequence_list="${OBSDIAG_DIR}/obs_sequence_list.txt"
    ${LIST} "${OBSDIAG_DIR}/${prefix}."*"${filename_out_suffix}" > "${obs_sequence_list}"
    cat "${obs_sequence_list}"

    # Extract timestamps like YYYY-MM-DD-SSSSS (unique + sorted)
    mapfile -t timestamps < <(grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2}-[0-9]{5}' "${obs_sequence_list}" | sort -u)

    if (( ${#timestamps[@]} == 0 )); then
        echo "ERROR: No timestamps found in ${obs_sequence_list}"
        exit 121
    fi

    first_ts="${timestamps[0]}"
    last_ts="${timestamps[$((${#timestamps[@]}-1))]}"

    # Seconds per DA window (e.g., 6h => 21600)
    DART_SECINC=$((3600 * DATA_ASSIMILATION_WINDOW))

    # Split components for first timestamp
    IFS='-' read -r DART_YEAR1 DART_MONTH1 DART_DAY1 DART_SEC1 <<< "${first_ts}"
    DART_HOUR1=$((10#${DART_SEC1} / 3600))
    DART_HOUR1e=$((DART_HOUR1 + DART_SECINC / 3600))   # window end hour

    # Split components for last timestamp
    IFS='-' read -r DART_YEAR2 DART_MONTH2 DART_DAY2 DART_SEC2 <<< "${last_ts}"
    DART_HOUR2=$((10#${DART_SEC2} / 3600))

    echo "[INFO] First timestamp: ${first_ts} → ${DART_YEAR1}-${DART_MONTH1}-${DART_DAY1} hour=${DART_HOUR1}"
    echo "[INFO] Last  timestamp: ${last_ts} → ${DART_YEAR2}-${DART_MONTH2}-${DART_DAY2} hour=${DART_HOUR2}"

    # Assemble bin times
    first_bin_center="${DART_YEAR1},${DART_MONTH1},${DART_DAY1},${DART_HOUR1},0,0"
    last_bin_center="${DART_YEAR2},${DART_MONTH2},${DART_DAY2},${DART_HOUR2},0,0"
    bin_separation="0, 0, 0, 6, 0, 0"
    bin_width="0, 0, 0, 6, 0, 0"
    time_to_skip="0, 0, 0, 6, 0, 0"
    #trusted_obs="'RADIOSONDE_TEMPERATURE', 'RADIOSONDE_SPECIFIC_HUMIDITY', 'RADIOSONDE_U_WIND_COMPONENT', 'RADIOSONDE_V_WIND_COMPONENT'"
    trusted_obs="'null'"
    obs_sequence_list="${OBSDIAG_DIR}/obs_sequence_list.txt"

    user_obs_diag_nl \
      "${refnen}" \
      "${DART_ON_PGRID}" \
      "${obs_sequence_list}" \
      "${first_bin_center}" \
      "${last_bin_center}" \
      "${bin_separation}" \
      "${bin_width}" \
      "${time_to_skip}" \
      "${trusted_obs}"

    echo "${first_bin_center} ${last_bin_center}"

    # obs_diag writes one shared NetCDF file and must not be launched as 24
    # independent writers. Preserve a partial file left by an interrupted run.
    [[ ! -e "${OBSDIAG_DIR}/obs_diag_output.nc" ]] ||
      ${MOVE} "${OBSDIAG_DIR}/obs_diag_output.nc" \
              "${OBSDIAG_DIR}/obs_diag_output.nc.incomplete.$(date +%s)"
    ${OBSDIAG_LAUNCHCMD} "${OBSDIAG_DIR}/obs_diag" || exit 122

    # rename the files
    fout="${OBSDIAG_DIR}/obs_diag_output.nc"
    DART_DATE1=`printf "%04d" ${DART_YEAR1}``printf "%02d" ${DART_MONTH1}``printf "%02d" ${DART_DAY1}``printf "%02d" ${DART_HOUR1}`
    DART_DATE2=`printf "%04d" ${DART_YEAR2}``printf "%02d" ${DART_MONTH2}``printf "%02d" ${DART_DAY2}``printf "%02d" ${DART_HOUR2}`
    [[ -s "${fout}" ]] || fail "obs_diag did not create ${fout}"
    ncdump -h "${fout}" >/dev/null 2>&1 || fail "invalid obs_diag output: ${fout}"
    ${MOVE} "${fout}" "${OBSDIAG_DIR}/${refcase}.dart.e.${scomp}_obs_diag_output.${DART_DATE1}-${DART_DATE2}.nc"
  done
fi

echo ===== End of e3sm dart diagnostic =====
date
echo =====================================
