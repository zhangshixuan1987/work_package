#!/bin/bash -el 
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --account=esmd
#SBATCH  --time=00:30:00
#SBATCH  --partition=short
#SBATCH  --job-name=e3sm_dart_ensda_pert 
#SBATCH  --nodes=20
#SBATCH  --output=e3sm_dart_ensda_pert.%j 

fail() {
   echo "ERROR: $*"
   exit 1
}

check_command() {
   local cmd="$1"
   command -v "${cmd}" >/dev/null 2>&1 || fail "required command not found: ${cmd}"
}

validate_positive_int() {
   local name="$1"
   local value="$2"
   if ! [[ "${value}" =~ ^[0-9]+$ ]] || (( value < 1 )); then
      fail "${name} must be a positive integer, got: ${value}"
   fi
}

SCRIPT_NAME=$(basename "$0")
RUN_START_EPOCH=$(date +%s)
on_exit_summary() {
   local exit_code=$?
   local run_end_epoch
   local elapsed
   run_end_epoch=$(date +%s)
   elapsed=$((run_end_epoch - RUN_START_EPOCH))
   if (( exit_code == 0 )); then
      echo "INFO: ${SCRIPT_NAME} completed successfully in ${elapsed}s at $(date '+%Y-%m-%d %H:%M:%S')"
   else
      echo "ERROR: ${SCRIPT_NAME} failed with exit ${exit_code} after ${elapsed}s at $(date '+%Y-%m-%d %H:%M:%S')"
   fi
}
trap on_exit_summary EXIT

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
   WORK_DIR=$(readlink -f "${SLURM_SUBMIT_DIR}") || fail "cannot resolve SLURM_SUBMIT_DIR"
else
   SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || fail "cannot resolve script path"
   WORK_DIR=$(dirname "${SCRIPT_PATH}")
fi
cd "${WORK_DIR}"
source "${WORK_DIR}/create_and_setup_case.sh"
[[ -n "${my_conda_setup_file:-}" ]] || fail "my_conda_setup_file is unset"
[[ -r "${my_conda_setup_file}" ]] || fail "configured Conda setup is not readable: ${my_conda_setup_file}"
[[ -n "${my_analysis_conda_env:-}" ]] || fail "my_analysis_conda_env is unset"
echo "Activating configured analysis environment: ${my_analysis_conda_env}"
source "${my_conda_setup_file}"
conda activate "${my_analysis_conda_env}" || fail "could not activate Conda environment: ${my_analysis_conda_env}"
mkdir -p "${my_log_dir}" "${my_status_dir}" "${my_lock_dir}"

for cmd in awk bc ex sed grep date dirname find flock mkdir cp ln ncdump ncks readlink rm srun wc; do
   check_command "${cmd}"
done

validate_positive_int "my_ensnum" "${my_ensnum}"
validate_positive_int "my_task_per_node" "${my_task_per_node}"
exec 8>"${my_lock_dir}/step3_perturb.lock"
flock -n 8 || fail "another step-3 ensemble perturbation is already running"

marker_field() { awk -F= -v key="$2" '$1 == key {sub(/^[^=]*=/, ""); print; found=1; exit} END {if (!found) exit 1}' "$1"; }
validate_marker() {
 local marker="$1" expected_time="$2" expected_case="$3" expected_size="$4"
 local actual_time actual_case actual_size actual_layout actual_dart_root
 [[ -s "${marker}" ]] || fail "missing upstream completion record: ${marker}"
 actual_time=$(marker_field "${marker}" valid_time) || fail "completion record lacks valid_time: ${marker}"
 actual_case=$(marker_field "${marker}" case) || fail "completion record lacks case: ${marker}"
 actual_size=$(marker_field "${marker}" ensemble_size) || fail "completion record lacks ensemble_size: ${marker}"
 actual_layout=$(marker_field "${marker}" archive_layout) || fail "completion record lacks archive_layout; rerun Step 2: ${marker}"
 actual_dart_root=$(marker_field "${marker}" dart_root) || fail "completion record lacks dart_root; rerun Step 2: ${marker}"
 [[ "${actual_time}" == "${expected_time}" && "${actual_case}" == "${expected_case}" && "${actual_size}" == "${expected_size}" ]] || fail "upstream completion record does not match time, case, or ensemble size: ${marker}"
 [[ "${actual_layout}" == "per_member" ]] || fail "Step 2 archive layout mismatch: expected per_member, found ${actual_layout}"
 [[ "${actual_dart_root}" == "${my_dart_root}" ]] || fail "Step 2 DART root mismatch: expected ${my_dart_root}, found ${actual_dart_root}"
}
validate_marker "${my_status_dir}/icbc_complete.${my_refdate}-${my_reftod}" "${my_refdate}-${my_reftod}" "${my_casename}" "${my_ensnum}"
PERTURB_MARKERS=()
for i in $(seq 1 "${my_ensnum}"); do
   ENSTR=$(printf 'EN%02d' "${i}")
   MEMBER_ARCHIVE_DIR="${my_modeldir}/${ENSTR}/archive"
   member_marker="${MEMBER_ARCHIVE_DIR}/rest/${my_refdate}-${my_reftod}/.dart_perturb_in_progress"
   [[ ! -e "${member_marker}" ]] || fail "previous perturbation did not finish cleanly for ${ENSTR}; rerun Step 2: ${member_marker}"
   PERTURB_MARKERS+=("${member_marker}")
done
perturb_status="${my_status_dir}/perturb_complete.${my_refdate}-${my_reftod}"
rm -f -- "${perturb_status}"

#For cshell:
#limit stacksize unlimited
#limit datasize unlimited

#For bash 
#ulimit -s unlimited
#ulimit -d unlimited

# ---------------------
# Purpose
# ---------------------
#
# This script is used to generate initial ensembles with DART 
# The perturbations were added to Temperature fields and then the model 
# was run for 15-days to spin up the ensemble 
#
#*******************************************************************************
# Set paths
E3SM_ROOT=${my_e3sm_code}
DART_ROOT=${my_eam_dart_code}
DART_MODEL=${my_eam_dart_model}
DART_SCPTDIR=${DART_ROOT}/models/${DART_MODEL}/shell_scripts
DART_WORKDIR=${DART_ROOT}/models/${DART_MODEL}/work
BASE_OBSDIR=${my_eam_dart_obsdir}
BASE_PHIS=${my_eam_topography_file} 
BASE_SEMAPS=${my_eam_se_mapping_file}
BASE_CSGRID=${my_eam_cs_grid_file}

# Run options
START_DATE=${my_refdate}
START_TOD=${my_reftod}
DART_CASE=${my_casename}
DART_ENSNUM=${my_ensnum}
DART_RUNDIR=${my_eam_dart_run_dir}
REQUEST_NODES=${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-${REQUEST_NODES:-}}}
validate_positive_int "REQUEST_NODES" "${REQUEST_NODES}"
DART_NTASKS=$((REQUEST_NODES * my_task_per_node))
DART_ON_PGRID=${my_eam_dart_pgrid}
DATA_ASSIMILATION_CYCLES=0
DATA_ASSIMILATION_WINDOW=${my_eam_dart_cycle_hours}

homme_map_file="SEMapping.nc"
cs_grid_file="SEMapping_cs_grid.nc"

# Assuming perturbations always applied to atmosphere model 
scomp="eam"

# ==============================================================================
# standard commands:
# Make sure that this script is using standard system commands
# instead of aliases defined by the user.
# If the standard commands are not in the location listed below,
# The 'force' (-f) options listed are added to commands where they are used.
# The verbose (-v) argument has been separated from these command definitions
# because these commands may not accept it on some systems.  On those systems
# set VERBOSE = ''
# ==============================================================================
# machine-specific dereferencing
# suppress "rm" warnings if wildcard does not match anything
nonomatch=1
case ${my_machine} in
        "compy")
                VERBOSE='-v'
                MOVE='/usr/bin/mv'
                COPY='/usr/bin/cp --preserve=timestamps'
                LINK='/usr/bin/ln -fs'
                LINKV=TRUE
                LIST='/usr/bin/ls'
                REMOVE='/usr/bin/rm -fr'
                LAUNCHCMD="srun --mpi=pmi2 --ntasks=${DART_NTASKS} --kill-on-bad-exit -l --cpu_bind=cores -c 1 -m plane=${my_task_per_node}"
                ;;
        "pm-cpu")
                VERBOSE='-v'
                MOVE='/usr/bin/mv'
                COPY='/usr/bin/cp --preserve=timestamps'
                LINK='/usr/bin/ln -fs'
                LINKV=TRUE
                LIST='/usr/bin/ls'
                REMOVE='/usr/bin/rm'
                LAUNCHCMD=mpirun.lsf 
                ;;
         *)
                VERBOSE='-v'
                MOVE='/usr/bin/mv'
                COPY='/usr/bin/cp --preserve=timestamps'
                LINK='/usr/bin/ln -fs'
                LINKV=TRUE
                LIST='/usr/bin/ls'
                REMOVE='/usr/bin/rm -fr'
                LAUNCHCMD="srun --mpi=pmi2 --ntasks=${DART_NTASKS} "
                ;;

esac

echo "`date` -- START EAM_DART_PERTURBATION"

# ==============================================================================
# Block 0: Set command environment
# ==============================================================================
# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.

echo "`date` -- BEGIN EAM_ASSIMILATE"

# ==============================================================================
# Make sure the DART executables exist or build them if we can't find them.
# The DART input.nml in the model directory IS IMPORTANT during this part
# because it defines what observation types are supported.
# ==============================================================================
targetdir=${DART_ROOT}/models/${DART_MODEL}/work
if [ ! -x ${targetdir}/filter ]; then
   if [[ "${ALLOW_DART_REBUILD:-FALSE}" != "TRUE" ]]; then
      echo "ERROR: DART filter is missing: ${targetdir}/filter"
      echo "ERROR: build DART before step 3 or explicitly set ALLOW_DART_REBUILD=TRUE"
      exit 41
   fi
   echo ""
   echo "WARNING: executable file 'filter' not found."
   echo "         Looking for: $targetdir/filter "
   echo "         Trying to rebuild all executables for $DART_MODEL now ..."
   echo "         This will be incorrect, if input.nml:preprocess_nml is not correct."
   cd $targetdir 
   ./quickbuild.sh mpi
   if [ ! -x ${targetdir}/filter ]; then
      echo "ERROR: executable file 'filter' not found."
      echo "       Unsuccessfully tried to rebuild: $targetdir/filter "
      echo "       Required DART assimilation executables are not found."
      echo "       Stopping prematurely."
      exit 01
   fi
fi

# ==============================================================================
# Make a place to perform DART data assimilation 
# st_archive can make a home for them.
# ==============================================================================
ATM_DATE_EXT=${START_DATE}-${START_TOD}
member_eam_initial_file() {
  local enstr="$1" member_archive
  member_archive="${my_modeldir}/${enstr}/archive"
  printf '%s/rest/%s-%s/%s.%s.eam.i.%s.nc\n' "${member_archive}" "${START_DATE}" "${START_TOD}" "${DART_CASE}" "${enstr}" "${ATM_DATE_EXT}"
}
validate_step3_input_ensemble() {
  local i enstr input_file
  for i in `seq 1 ${DART_ENSNUM}`; do
    enstr=EN`printf "%02d" ${i}`
    input_file=$(member_eam_initial_file "${enstr}") || fail "cannot resolve EAM state for ${enstr}"
    if [ ! -s "${input_file}" ] || ! ncdump -h "${input_file}" >/dev/null 2>&1; then
      fail "missing or invalid step-2 EAM state for ${enstr}: ${input_file}"
    fi
    if ! ncks -m -v PS,U,V,T,Q "${input_file}" >/dev/null 2>&1; then
      fail "step-2 EAM state lacks core variables for ${enstr}: ${input_file}"
    fi
  done
}

validate_step3_dependencies() {
  local required_file
  for required_file in "${DART_WORKDIR}/filter" "${DART_WORKDIR}/perfect_model_obs" "${DART_WORKDIR}/fill_inflation_restart" "${my_eam_perturb_nml}" "${DART_WORKDIR}/qceff_table.csv" "${BASE_PHIS}"; do
    [ -s "${required_file}" ] || fail "missing required step-3 dependency: ${required_file}"
  done
  if [ ! -s "${BASE_SEMAPS}" ] && [ ! -s "${BASE_CSGRID}" ]; then
    fail "neither EAM mapping file is available: ${BASE_SEMAPS} or ${BASE_CSGRID}"
  fi
}

safe_reset_dart_workdir() {
  local target="$1"
  case "${target}" in
    "${DART_RUNDIR}/${START_DATE}-${START_TOD}") ;;
    *) fail "refusing to clean unexpected DART work directory: ${target}" ;;
  esac
  mkdir -p "${target}"
  find "${target}" -mindepth 1 -maxdepth 1 -exec rm -rf -- {} +
}

validate_step3_dependencies
validate_step3_input_ensemble
for member_marker in "${PERTURB_MARKERS[@]}"; do
   printf 'valid_time=%s-%s\ncase=%s\nensemble_size=%s\narchive_layout=%s\ndart_root=%s\nslurm_job_id=%s\nstarted_at=%s\n' "${my_refdate}" "${my_reftod}" "${my_casename}" "${my_ensnum}" "per_member" "${my_dart_root}" "${SLURM_JOB_ID:-none}" "$(date '+%Y-%m-%d %H:%M:%S')" > "${member_marker}"
done
CURRENT_DADIR="${DART_RUNDIR}/${START_DATE}-${START_TOD}"
safe_reset_dart_workdir "${CURRENT_DADIR}"

#=========================================================================
# Loop over members and link eam files for data assimilation
# As implemented, the input filenames are static in the DART namelists.
# We must link the new uniquely-named files to static names.
#=========================================================================
for i in `seq 1 ${DART_ENSNUM}`;do
  cd ${CURRENT_DADIR}
  ENSTR=EN`printf "%02d" ${i}`
  ATM_INITIAL_FILENAME=$(member_eam_initial_file "${ENSTR}") || fail "cannot resolve EAM state for ${ENSTR}"
  if [ ! -s "${ATM_INITIAL_FILENAME}" ] || ! ncdump -h "${ATM_INITIAL_FILENAME}" >/dev/null 2>&1; then
    echo "ERROR: required file missing ${ATM_INITIAL_FILENAME}"
    exit 1
  fi 
  inst_string=`printf _%04d ${i}`
  ATM_DART_FILENAME=${DART_CASE}.eam${inst_string}.i.${ATM_DATE_EXT}.nc
  #echo $ATM_INITIAL_FILENAME 
  #echo $ATM_DART_FILENAME
  if [ $LINKV == TRUE ]; then
    echo "Linking $ATM_INITIAL_FILENAME  $ATM_DART_FILENAME"
  fi
  $LINK ${ATM_INITIAL_FILENAME} ${ATM_DART_FILENAME} || exit 03
  if [ ! -L "${ATM_DART_FILENAME}" ] || [ "$(readlink -f "${ATM_DART_FILENAME}")" != "$(readlink -f "${ATM_INITIAL_FILENAME}")" ]; then
    fail "incorrect DART member link for ${ENSTR}: ${ATM_DART_FILENAME}"
  fi
done

# ==============================================================================
# Block 1: Determine time of current model state from file name of member 1
# These are of the form "${CASE}.eam_${ensemble_member}.i.2000-01-06-00000.nc"
# ==============================================================================
cd ${CURRENT_DADIR}
ATM_DATE=( `echo $ATM_DATE_EXT | sed -e "s#-# #g"` )
ATM_YEAR=`echo "${ATM_DATE[0]}" | bc`
ATM_MONTH=`echo "${ATM_DATE[1]}" | bc`
ATM_DAY=`echo "${ATM_DATE[2]}" | bc`
ATM_SECONDS=`echo "${ATM_DATE[3]}" | bc`
ATM_HOUR=`echo "${ATM_DATE[3]}" / 3600 | bc`
echo "valid time of model is $ATM_YEAR $ATM_MONTH $ATM_DAY $ATM_SECONDS (seconds)"
echo "valid time of model is $ATM_YEAR $ATM_MONTH $ATM_DAY $ATM_HOUR (hours)"

#determine the time for previous DA cycle 
OLD_DATE=`date -d "$ATM_HOUR:00 ${ATM_YEAR}-${ATM_MONTH}-${ATM_DAY} -${DATA_ASSIMILATION_WINDOW} hours" +"%Y-%m-%d %H"`
OLD_DATE=( `echo $OLD_DATE | sed -e "s#-# #g"` )
OLD_YEAR=`echo "${OLD_DATE[0]}" | bc`
OLD_MONTH=`echo "${OLD_DATE[1]}" | bc`
OLD_DAY=`echo "${OLD_DATE[2]}" | bc`
OLD_HOUR=`echo "${OLD_DATE[3]}" | bc`
OLD_SECONDS=`echo "${OLD_DATE[3]}" \* 3600 | bc`
PRE_DATE=`printf "%04d" ${OLD_YEAR}`-`printf "%02d" ${OLD_MONTH}`-`printf "%02d" ${OLD_DAY}`
PRE_TOD=`printf "%05d" ${OLD_SECONDS}`
echo "valid time for previous DA cycle is $OLD_YEAR $OLD_MONTH $OLD_DAY $OLD_SECONDS (seconds)"
echo "valid time for previous DA cycle is $OLD_YEAR $OLD_MONTH $OLD_DAY $OLD_HOUR (hours)"

# ==============================================================================
# Block 2: Populate a run-time directory with the input needed to run DART.
# ==============================================================================
echo "`date` -- BEGIN COPY BLOCK"
cd ${CURRENT_DADIR}
# Put a pared down copy (no comments) of input.nml in this assimilate_eam directory.
# The contents may change from one cycle to the next, so always start from
# the known configuration in the CASE_ROOT directory.
if [ -e "${my_eam_perturb_nml}" ]; then
  ${COPY} "${my_eam_perturb_nml}" input.nml  || exit 04
  sed -i "/#/d;/^\!/d;/^[ ]*\!/d" input.nml
  sed -i '1,1i\WARNING: Changes to this runtime file will be ignored. \n Edit workflow_lib/namelists/eam/perturb.nml instead.\n\n\n' input.nml
else
  echo "ERROR ... workflow EAM perturbation namelist ${my_eam_perturb_nml} not found ... ERROR"
  exit 05
fi

#This file is needed for DART code after 2025/02/01
if [ -e "${DART_WORKDIR}/qceff_table.csv" ]; then
  ${COPY} ${DART_WORKDIR}/qceff_table.csv qceff_table.csv  || exit 04
else
  echo "ERROR ... DART required file ${DART_WORKDIR}/qceff_table.csv not found ... ERROR"
  exit 05
fi

# Ensure that the input.nml ensemble size matches the number of instances.
# WARNING: the output files contain ALL ensemble members ==> BIG
ex input.nml <<ex_end
g;ens_size ;s;= .*;= ${DART_ENSNUM};
g;num_output_state_members ;s;= .*;= ${DART_ENSNUM};
g;num_output_obs_members ;s;= .*;= ${DART_ENSNUM};
g;eam_use_pgrid ;s;= .*;= ${DART_ON_PGRID};
wq
ex_end

if ! grep -Eiq '^[[:space:]]*perturb_from_single_instance[[:space:]]*=[[:space:]]*\.true\.' input.nml; then
   fail "the EAM perturbation template must enable perturb_from_single_instance"
fi
if ! grep -Eq "fields_to_perturb[[:space:]]*=.*QTY_TEMPERATURE" input.nml; then
   fail "the EAM perturbation template must perturb QTY_TEMPERATURE"
fi

list=`grep '^[ ]*vertical_localization_coord' input.nml`
list=( `echo $list | sed -e "s#[=,']# #g"` )
if [ "${list[1]}" == "SCALEHEIGHT" ]; then
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
if [[ -n "${my_task_per_node:-}" ]]; then
   sed -i "s#layout.*#layout = 2#"  input.nml
   sed -i "s#tasks_per_node.*#tasks_per_node = $my_task_per_node#" input.nml
fi
#echo ${my_task_per_node}

#stage_dart_files 
${COPY} -f ${DART_WORKDIR}/filter                 ${CURRENT_DADIR} || exit 06
${COPY} -f ${DART_WORKDIR}/perfect_model_obs      ${CURRENT_DADIR} || exit 07
${COPY} -f ${DART_WORKDIR}/fill_inflation_restart ${CURRENT_DADIR} || exit 08

if [ $DART_ENSNUM -gt 1 ] ; then
   SAMP_ERR_DIR=assimilation_code/programs/gen_sampling_err_table/work
   SAMP_ERR_FILE=${DART_ROOT}/${SAMP_ERR_DIR}/sampling_error_correction_table.nc
   if [ -e ${SAMP_ERR_FILE} ]; then
      ${COPY} -f ${VERBOSE} ${SAMP_ERR_FILE} ${CURRENT_DADIR}  || exit 09
      if [ ${DART_ENSNUM} -lt 3 ] || [ ${DART_ENSNUM} -gt 200 ]; then
         echo ""
         echo "ERROR: sampling_error_correction_table.nc handles ensemble sizes 3...200."
         echo "ERROR: Yours is $DART_ENSNUM"
         echo ""
         exit 10
      fi
   else
      list=`grep sampling_error_correction input.nml` 
      list=( `echo $list | sed -e "s/[=\.,]//g"` ) 
      if [ ${list[1]} == "true" ]; then
         echo ""
         echo "ERROR: No sampling_error_correction_table.nc file found ..."
         echo "ERROR: the input.nml:assim_tool_nml:sampling_error_correction"
         echo "ERROR: is 'true' so this file must exist."
         echo ""
         exit 11
      fi
   fi

fi

echo "`date` -- END COPY BLOCK"

# ==============================================================================
# Block 3: Identify requested output stages, warn about redundant output.
# ==============================================================================
cd ${CURRENT_DADIR}

MYSTRING=`grep stages_to_write input.nml`
MYSTRING=( `echo $MYSTRING | sed -e "s#[=,'\.]# #g"` )
STAGE_input=FALSE
STAGE_forecast=FALSE
STAGE_preassim=FALSE
STAGE_postassim=FALSE
STAGE_analysis=FALSE
STAGE_output=FALSE

# Assemble lists of stages to write out, which are not the 'output' stage.
stages_except_output="{"
stage=1
nstage=`expr ${#MYSTRING[@]} - 1`
while [ $stage -le ${nstage} ];do
  if [ ${MYSTRING[$stage]} == 'input' ]; then
      STAGE_input=TRUE
      if [ $stage -gt 1 ]; then 
        stages_except_output="${stages_except_output},"
      fi 
      stages_except_output="${stages_except_output}input"
   fi
   if [ ${MYSTRING[$stage]} == 'forecast' ]; then
      STAGE_forecast=TRUE
      if [ $stage -gt 1 ]; then 
        stages_except_output="${stages_except_output},"
      fi 
      stages_except_output="${stages_except_output}forecast"
   fi
   if [ ${MYSTRING[$stage]} == 'preassim' ]; then
      STAGE_preassim=TRUE
      if [ $stage -gt 1 ]; then 
        stages_except_output="${stages_except_output},"
      fi
      stages_except_output="${stages_except_output}preassim"
   fi
   if [ ${MYSTRING[$stage]} == 'postassim' ]; then
      STAGE_postassim=TRUE
      if [ $stage -gt 1 ]; then 
        stages_except_output="${stages_except_output},"
      fi
      stages_except_output="${stages_except_output}postassim"
   fi
   if [ ${MYSTRING[$stage]} == 'analysis' ]; then
      STAGE_analysis=TRUE
      if [ $stage -gt 1 ]; then 
        stages_except_output="${stages_except_output},"
      fi
      stages_except_output="${stages_except_output}analysis"
   fi
   if [ $stage == ${nstage} ]; then
      stages_all="${stages_except_output}"
      if [ ${MYSTRING[$stage]} == 'output' ]; then
        STAGE_output=TRUE
        stages_all="${stages_all},output"
      fi
   fi
   stage=$((stage + 1))
done 

# Add the closing }
stages_all="${stages_all}}"
stages_except_output="${stages_except_output}}"

# Checking
echo "stages_except_output = $stages_except_output"
echo "stages_all = $stages_all"
if [ ${STAGE_output} != TRUE ];  then
   echo "ERROR: assimilate.csh requires that input.nml:filter_nml:stages_to_write includes stage 'output'"
   exit 12
fi

# ==============================================================================
# Block 5: Get observation sequence file ... or die right away.
# The observation file names have a time that matches the stopping time of EAM.
#
# Make sure the file name structure matches the obs you will be using.
# PERFECT model obs output appends .perfect to the filenames
# ==============================================================================
cd ${CURRENT_DADIR}
YYYYMM=`printf %04d%02d ${ATM_YEAR} ${ATM_MONTH}`
if [ ! -d ${BASE_OBSDIR}/${YYYYMM}_6H_CESM ]; then
   echo "E3SM+DART requires 6 hourly obs_seq files in directories of the form YYYYMM_6H_CESM"
   echo "The directory ${BASE_OBSDIR}/${YYYYMM}_6H_CESM is not found.  Exiting"
   exit 13
fi

OBSFNAME=`printf obs_seq.%04d-%02d-%02d-%05d ${ATM_YEAR} ${ATM_MONTH} ${ATM_DAY} ${ATM_SECONDS}`
OBS_FILE=${BASE_OBSDIR}/${YYYYMM}_6H_CESM/${OBSFNAME}
#echo "OBS_FILE = $OBS_FILE"

${REMOVE} obs_seq.out
if [ -e ${OBS_FILE} ]; then
   ${LINK} ${OBS_FILE} obs_seq.out || exit 14
else
   echo "ERROR ... no observation file ${OBS_FILE}"
   echo "ERROR ... no observation file ${OBS_FILE}"
   exit 15
fi

# ==============================================================================
# Block 6: DART INFLATION
# This block is only relevant if 'inflation' is turned on AND
# inflation values change through time:
# filter_nml
#    inf_flavor(:)  = 2  (or 3 (or 4 for posterior))
#    inf_initial_from_restart    = .TRUE.
#    inf_sd_initial_from_restart = .TRUE.
#
#
# This block stages the files that contain the inflation values.
# The inflation files are essentially duplicates of the DART model state,
# which have names in the E3SM style, something like
#    ${case}.dart.rh.${scomp}_output_priorinf_{mean,sd}.YYYY-MM-DD-SSSSS.nc
# The strategy is to use the latest such files in ${RUNDIR}.
# If those don't exist at the start of an assimilation,
# this block creates them with 'fill_inflation_restart'.
# If they don't exist AFTER the first cycle, the script will exit
# because they should have been available from a previous cycle.
# The script does NOT check the model date of the files for consistency
# with the current forecast time, so check that the inflation mean
# files are evolving as expected.
#
# E3SM's st_archive should archive the inflation restart files
# like any other "restart history" (.rh.) files; copying the latest files
# to the archive directory, and moving all of the older ones.
# ==============================================================================

cd ${CURRENT_DADIR}

# If we need to run fill_inflation_restart, EAM:static_init_model()
# always needs a eaminput.nc and a eam_phis.nc for geometry information, etc.
MYSTRING=`grep eam_template_filename input.nml`
MYSTRING=( `echo $MYSTRING | sed -e "s#[=,']# #g"` )
EAMINPUT=${MYSTRING[1]}
${REMOVE} ${EAMINPUT}
${LINK} ${DART_CASE}.eam_0001.i.${ATM_DATE_EXT}.nc ${EAMINPUT} || exit 16

MYSTRING=`grep eam_phis_filename input.nml`
MYSTRING=( `echo $MYSTRING | sed -e "s#[=,']# #g"` )
EAM_PHIS=${MYSTRING[1]}
${REMOVE} ${EAM_PHIS}
${LINK} ${BASE_PHIS} ${EAM_PHIS} || exit 17

#Now, Link the grid information files 
if [ ! -f ${BASE_SEMAPS} ] && [ ! -f ${BASE_CSGRID} ]; then 
  echo "ERROR ... no mapping file ${homme_map_file}"
  echo "ERROR ... no gridinfo file ${cs_grid_file}"
  echo "ERROR ... must provide either of them"
  exit 18    
else 
  if [ -f ${BASE_SEMAPS} ]; then
    ${REMOVE} ${homme_map_file}
    ${COPY} -r ${BASE_SEMAPS} ${homme_map_file} || exit 19
  fi
  if [ -f ${BASE_CSGRID} ]; then  
    ${REMOVE} ${cs_grid_file}
    ${COPY} -r ${BASE_CSGRID} ${cs_grid_file} || exit 20
  fi
fi 

# Now, actually check the inflation settings
MYSTRING=`grep inf_flavor input.nml`
MYSTRING=( `echo $MYSTRING | sed -e "s#[=,'\.]# #g"` )
PRIOR_INF=${MYSTRING[1]}
POSTE_INF=${MYSTRING[2]}
#echo $PRIOR_INF $POSTE_INF

MYSTRING=`grep inf_initial_from_restart input.nml`
MYSTRING=( `echo $MYSTRING | sed -e "s#[=,'\.]# #g"` )
#echo ${MYSTRING[@]}

# If no inflation is requested, the inflation restart source is ignored
if [ ${PRIOR_INF} -eq 0 ];  then
  PRIOR_INFLATION_FROM_RESTART=ignored
  USING_PRIOR_INFLATION=false
else
  PRIOR_INFLATION_FROM_RESTART=`echo ${MYSTRING[1]} | tr '[:upper:]' '[:lower:]'`
  USING_PRIOR_INFLATION=true
fi

if [ ${POSTE_INF} -eq 0 ]; then
  POSTE_INFLATION_FROM_RESTART=ignored
  USING_POSTE_INFLATION=false
else
  POSTE_INFLATION_FROM_RESTART=`echo ${MYSTRING[2]} | tr '[:upper:]' '[:lower:]'`
  USING_POSTE_INFLATION=true
fi

#echo $PRIOR_INF $PRIOR_INFLATION_FROM_RESTART $USING_PRIOR_INFLATION 
#echo $POSTE_INF $POSTE_INFLATION_FROM_RESTART $USING_POSTE_INFLATION 

if [ ${USING_PRIOR_INFLATION} == false ]; then
   stages_requested=0
   if [ ${STAGE_input}  == TRUE ]; then 
     stages_requested=$((stages_requested+1)) 
   fi
   if [ ${STAGE_forecast} == TRUE ]; then 
     stages_requested=$((stages_requested+1))
   fi 
   if [ ${STAGE_preassim} == TRUE ]; then
     stages_requested=$((stages_requested+1))
   fi
   if [ ${stages_requested} -gt 1 ]; then
      echo " "
      echo "WARNING ! ! Redundant output is requested at multiple stages before assimilation."
      echo "            Stages 'input' and 'forecast' are always redundant."
      echo "            Prior inflation is OFF, so stage 'preassim' is also redundant. "
      echo "            We recommend requesting just 'preassim'."
      echo " "
   fi
fi

#echo $STAGE_input $STAGE_forecast $STAGE_preassim $stages_requested 

if [ ${USING_POSTE_INFLATION} == false ]; then
   stages_requested=0
   if [ ${STAGE_postassim} == TRUE ]; then 
      stages_requested=$((stages_requested+1))
   fi
   if [ ${STAGE_analysis}  == TRUE ]; then 
     stages_requested=$((stages_requested+1))
   fi
   if [ ${STAGE_output}    == TRUE ]; then 
     stages_requested=$((stages_requested+1))
   fi
   if [ ${stages_requested} -gt 1 ];  then
      echo " "
      echo "WARNING ! ! Redundant output is requested at multiple stages after assimilation."
      echo "            Stages 'output' and 'analysis' are always redundant."
      echo "            Posterior inflation is OFF, so stage 'postassim' is also redundant. "
      echo "            We recommend requesting just 'output'."
      echo " "
   fi
fi

#echo $STAGE_postassim $STAGE_analysis $STAGE_output $stages_requested

# IF we want PRIOR inflation:
if [ ${USING_PRIOR_INFLATION} == true ]; then
   if [ ${PRIOR_INFLATION_FROM_RESTART} == false ]; then
      echo "inf_flavor(1) = $PRIOR_INF, using namelist values."
   else
      # Look for the output from the previous assimilation (or fill_inflation_restart)
      # If inflation files exists, use them as input for this assimilation
      OLD_DADIR="${DART_RUNDIR}/${PRE_DATE}-${PRE_TOD}"
      (${LIST} -rt1 ${OLD_DADIR}/*.dart.rh.${scomp}_output_priorinf_mean* | tail -n 1 >  latestfile) >& /dev/null
      (${LIST} -rt1 ${OLD_DADIR}/*.dart.rh.${scomp}_output_priorinf_sd*   | tail -n 1 >> latestfile) >& /dev/null
      nfiles=`cat latestfile | wc -l`
      # If one exists, use it as input for this assimilation
      if [ ${nfiles} -gt 0 ]; then
         latest_mean=`head -n 1 latestfile`
         latest_sd=`tail -n 1 latestfile`
         # Need to COPY instead of link because of short-term archiver and disk management.
         ${COPY} $latest_mean input_priorinf_mean.nc
         ${COPY} $latest_sd   input_priorinf_sd.nc
      elif [ ${DATA_ASSIMILATION_CYCLES} -eq 0 ]; then
         # It's the first assimilation; try to find some inflation restart files
         # or make them using fill_inflation_restart.
         # Fill_inflation_restart needs eaminput.nc and eam_phis.nc for static_model_init,
         # so this staging is done in assimilate.csh (after a forecast) instead of stage_e3sm_files.
         if [ -x ${CURRENT_DADIR}/fill_inflation_restart ]; then
            ${CURRENT_DADIR}/fill_inflation_restart
         else
            echo "ERROR: Requested PRIOR inflation restart for the first cycle."
            echo "       There are no existing inflation files available "
            echo "       and ${CURRENT_DADIR}/fill_inflation_restart is missing."
            echo "EXITING"
            exit 21
         fi
      else
         echo "ERROR: Requested PRIOR inflation restart, "
         echo "       but files *.dart.rh.${scomp}_output_priorinf_* do not exist in the ${CURRENT_DADIR}."
         echo "       If you are changing from eam_no_assimilate.csh to assimilate.csh,"
         echo "       you might be able to continue by changing CONTINUE_RUN = FALSE for this cycle,"
         echo "       and restaging the initial ensemble."
         ${LIST} -l *inf*
         echo "EXITING"
         exit 22
      fi
   fi
else
   echo "Prior Inflation not requested for this assimilation."
fi
      
# POSTERIOR: We look for the 'newest' and use it - IFF we need it.
if [ ${USING_POSTE_INFLATION} == true ] ; then
   if [ ${POSTE_INFLATION_FROM_RESTART} == false ]; then
      # we are not using an existing inflation file.
      echo "inf_flavor(2) = $POSTE_INF, using namelist values."
   else
      # Look for the output from the previous assimilation (or fill_inflation_restart).
      # (The only stage after posterior inflation.)
      OLD_DADIR="${DART_RUNDIR}/${PRE_DATE}-${PRE_TOD}"
      (${LIST} -rt1 ${OLD_DADIR}/*.dart.rh.${scomp}_output_postinf_mean* | tail -n 1 >  latestfile) >& /dev/null
      (${LIST} -rt1 ${OLD_DADIR}/*.dart.rh.${scomp}_output_postinf_sd*   | tail -n 1 >> latestfile) >& /dev/null
      nfiles=`cat latestfile | wc -l`
      # If one exists, use it as input for this assimilation
      if [ $nfiles -gt 0 ] ; then
         latest_mean=`head -n 1 latestfile`
         latest_sd=`tail -n 1 latestfile`
         ${LINK} $latest_mean input_postinf_mean.nc || exit 23
         ${LINK} $latest_sd   input_postinf_sd.nc   || exit 24
      elif [ ${DATA_ASSIMILATION_CYCLES} -eq 0 ]; then
         # It's the first assimilation; try to find some inflation restart files
         # or make them using fill_inflation_restart.
         # Fill_inflation_restart needs eaminput.nc and eam_phis.nc for static_model_init,
         # so this staging is done in assimilate.csh (after a forecast).
         if [ -x ${CURRENT_DADIR}/fill_inflation_restart ]; then
            ${CURRENT_DADIR}/fill_inflation_restart
            ${MOVE} prior_inflation_mean.nc input_postinf_mean.nc || exit 25
            ${MOVE} prior_inflation_sd.nc   input_postinf_sd.nc   || exit 26
         else
            echo "ERROR: Requested POSTERIOR inflation restart for the first cycle."
            echo "       There are no existing inflation files available "
            echo "       and ${CURRENT_DADIR}/fill_inflation_restart is missing."
            echo "EXITING"
            exit 27
         fi
      else
         echo "ERROR: Requested POSTERIOR inflation restart, "
         echo "       but files *.dart.rh.${scomp}_output_postinf_* do not exist in the ${CURRENT_DADIR}."
         ${LIST} -l *inf*
         echo "EXITING"
         exit 28
      fi
   fi
else
   echo "Posterior Inflation not requested for this assimilation."
fi

# ==============================================================================
# Block 7: Actually run the assimilation.
#
# DART namelist settings required:
# &filter_nml
#    adv_ens_command         = "no_eam-se_advance_script",
#    obs_sequence_in_name    = 'obs_seq.out'
#    obs_sequence_out_name   = 'obs_seq.final'
#    single_file_in          = .false.,
#    single_file_out         = .false.,
#    stages_to_write         = stages you want + ,'output'
#    input_state_file_list   = 'eam_init_files'
#    output_state_file_list  = 'eam_init_files',
#
# WARNING: the default mode of this script assumes that
#          input_state_file_list = output_state_file_list, so that
#          the EAM initial files used as input to filter will be overwritten.
#          The input model states can be preserved by requesting that stage
#          'forecast' be output.
#
# ==============================================================================

# In the default mode of EAM assimilations, filter gets the model state(s)
# from EAM initial files.  This section puts the names of those files into a text file.
# The name of the text file is provided to filter in filter_nml:input_state_file_list.

# NOTE:
# If the files in input_state_file_list are eam-se initial files (all vars and
# all meta data), then they will end up with a different structure than
# the non-'output', stage output written by filter ('preassim', 'postassim', etc.).
# This can be prevented (at the cost of more disk space) by copying
# the eam-se format initial files into the names filter will use for preassim, etc.:
#    > cp $case.eam_0001.i.$date.nc  preassim_member_0001.nc.
#    > ... for all members
# Filter will replace the state variables in preassim_member* with updated versions,
# but leave the other variables and all metadata unchanged.

# If filter will create an ensemble from a single state,
#    filter_nml: perturb_from_single_instance = .true.
# it's fine (and convenient) to put the whole list of files in input_state_file_list.
# Filter will just use the first as the base to perturb.

cd ${CURRENT_DADIR} 

line=( `grep input_state_file_list input.nml | sed -e "s#[=,'\.]# #g"` )
input_file_list_name=${line[1]}

${LIST} -1 ${DART_CASE}.eam_[0-9][0-9][0-9][0-9].i.${ATM_DATE_EXT}.nc > $input_file_list_name
member_file_count=$(wc -l < "${input_file_list_name}")
if [ "${member_file_count}" -ne "${DART_ENSNUM}" ]; then
   fail "DART member list has ${member_file_count} entries; expected ${DART_ENSNUM}"
fi

# If the file names in $output_state_file_list = names in $input_state_file_list,
# then the restart file contents will be overwritten with the states updated by DART.
line=( `grep output_state_file_list input.nml | sed -e "s#[=,'\.]# #g"` )
output_file_list_name=${line[1]}

if [ $input_file_list_name != $output_file_list_name ]; then
   echo "ERROR: assimilate.csh requires that input_file_list = output_file_list"
   echo "       You can probably find the data you want in stage 'forecast'."
   echo "       If you truly require separate copies of EAM's initial files"
   echo "       before and after the assimilation, see revision 12603, and note that"
   echo "       it requires changing the linking to eam_initial_####.nc, below."
   exit 29
fi

#TEMPSZ: comment out for production run
echo "`date` -- BEGIN FILTER"
${LAUNCHCMD} ${CURRENT_DADIR}/filter || exit 30
echo "`date` -- END FILTER"

for i in `seq 1 ${DART_ENSNUM}`; do
   inst_string=`printf _%04d ${i}`
   perturbed_file="${CURRENT_DADIR}/${DART_CASE}.eam${inst_string}.i.${ATM_DATE_EXT}.nc"
   if [ ! -s "${perturbed_file}" ] || ! ncdump -h "${perturbed_file}" >/dev/null 2>&1; then
      echo "ERROR: invalid perturbed member ${i}: ${perturbed_file}"
      exit 39
   fi
   if ! ncks -m -v PS,U,V,T,Q "${perturbed_file}" >/dev/null 2>&1; then
      echo "ERROR: perturbed member ${i} lacks core EAM variables: ${perturbed_file}"
      exit 42
   fi
done

# ==============================================================================
# Block 8: Rename the output using the eam-se file-naming convention.
# ==============================================================================
# If output_state_file_list is filled with custom (eam-se) filenames,
# then 'output' ensemble members will not appear with filter's default,
# hard-wired names.  But file types output_{mean,sd} will appear and be
# renamed here.
#
# We don't know the exact set of files which will be written,
# so loop over all possibilities: use LIST in the foreach.
# LIST will expand the variables and wildcards, only existing files will be
# in the foreach loop. (If the input.nml has num_output_state_members = 0,
# there will be no output_member_xxxx.nc even though the 'output' stage
# may be requested - for the mean and sd) 
#
# Handle files with instance numbers first.
#    split off the .nc
#    separate the pieces of the remainder
#    grab all but the trailing 'member' and #### parts.
#    and join them back together

echo "`date` -- BEGIN FILE RENAMING"

# The short-term archiver archives files depending on pieces of their names.
# '_####.i.' files are eam-se initial files.
# '.dart.i.' files are ensemble statistics (mean, sd) of just the state variables 
#            in the initial files.
# '.e.'      designates a file as something from the 'external system processing ESP', e.g. DART.

stages_all=`echo $stages_all | sed -e "s#"}"##g"`
stages_all=`echo $stages_all | sed -e "s#"{"##g"`
stages_all=(`echo $stages_all | sed -e "s#\,# #g"`)
for stage in ${stages_all[@]}; do 
  for FILE in `${LIST} ${stage}_member_*.nc` ; do
    parts=( `echo $FILE | sed -e "s#\.# #g"` )
    list=( `echo ${parts[0]}  | sed -e "s#_# #g"` )
    last=`expr ${#list[@]} - 2`
    dart_file=( `echo ${list[$last-1]} | sed -e "s# #_#g"` )
    # DART 'output_member_****.nc' files are actually linked to eam input files
    echo $FILE > dart_tmp.txt 
    if [ `grep -c "put" dart_tmp.txt` -gt 0 ] ; then 
      type="i"
    else 
      type="e"
    fi
    #echo ${DART_CASE}.${scomp}_${list[${#list[@]}-1]}.${type}.${dart_file}.${ATM_DATE_EXT} ; exit
    ${MOVE} $FILE \
        ${DART_CASE}.${scomp}_${list[${#list[@]}-1]}.${type}.${dart_file}.${ATM_DATE_EXT} || exit 43
    ${REMOVE} dart_tmp.txt
  done
done 

# Files without instance numbers need to have the scomp part of their names = "dart".
# This is because in st_archive, all files with  scomp = "eam"
# (= compname in env_archive.xml) will be st_archived using a pattern
# which has the instance number added onto it.  {mean,sd} files don't have 
# instance numbers, so they need to be archived by the "dart" section of env_archive.xml.
# But they still need to be different for each component, so include $scomp in the
# ".dart_file" part of the file name.  Somewhat awkward and inconsistent, but effective.

# Means and standard deviation files (except for inflation).
for stage in ${stages_all[@]}; do
  for FILE in `${LIST} ${stage}_{mean,sd}*.nc`; do 
     parts=( `echo $FILE | sed -e "s#\.# #g"` )
     list=( `echo ${parts[0]}  | sed -e "s#_# #g"` )
     #last=`expr ${#list[@]} - 2`
     #dart_file=( `echo ${list[$last-1]} | sed -e "s# #_#g"` )
     echo $FILE > dart_tmp.txt
     if [ `grep -c "put" dart_tmp.txt` -gt 0 ] ; then
       type="i"
     else
       type="e"
     fi
     echo ${DART_CASE}.dart.${type}.${scomp}_${parts[0]}.${ATM_DATE_EXT}.nc 
     ${MOVE} $FILE ${DART_CASE}.dart.${type}.${scomp}_${parts[0]}.${ATM_DATE_EXT}.nc || exit 31
  done 
done 

# Rename the observation file and run-time output
${MOVE} obs_seq.final ${DART_CASE}.dart.e.${scomp}_obs_seq_final.${ATM_DATE_EXT} || exit 32
${MOVE} dart_log.out  ${scomp}_dart_log.${ATM_DATE_EXT}.out || exit 33

# Rename the inflation files and designate them as 'rh' files - which get
# reinstated in the run directory by the short-term archiver and are then
# available for the next assimilation cycle.
#
# Accommodate any possible inflation files.
# The .${scomp}_ part is needed by DART to distinguish
# between inflation files from separate components in coupled assims.
for stage in ${stages_all[@]}; do
  for FILE in `${LIST} ${stage}_{prior,post}inf_*`; do
   parts=( `echo $FILE | sed -e "s#\.# #g"` ) 
   ${MOVE} $FILE  ${DART_CASE}.dart.rh.${scomp}_${parts[0]}.${ATM_DATE_EXT}.nc || exit 34
  done 
done 

# Handle localization_diagnostics_files
MYSTRING=`grep 'localization_diagnostics_file' input.nml`
MYSTRING=`echo $MYSTRING | sed -e "s#[=,']# #g"`
MYSTRING=( `echo $MYSTRING | sed -e 's#"# #g'` )
loc_diag=$MYSTRING[1]
if [ -f $loc_diag ]; then
   ${MOVE} $loc_diag  ${scomp}_${loc_diag}.dart.e.${ATM_DATE_EXT} || exit 35
fi

# Handle regression diagnostics
MYSTRING=`grep 'reg_diagnostics_file' input.nml`
MYSTRING=`echo $MYSTRING | sed -e "s#[=,']# #g"`
MYSTRING=( `echo $MYSTRING | sed -e 's#"# #g'` ) 
reg_diag=$MYSTRING[1]
if [ -f $reg_diag ] ; then
   ${MOVE} $reg_diag  ${scomp}_${reg_diag}.dart.e.${ATM_DATE_EXT} || exit 36
fi 

# Then this script will need to feed the files in output_restart_list_file
# to the next model advance.
# This gets the .i. or .r. piece from the EAM-SE format file name.

# Step 3 is complete only after filter output validation and all diagnostic
# post-processing have succeeded.
perturb_tmp="${perturb_status}.tmp.${SLURM_JOB_ID:-$$}"
printf 'valid_time=%s-%s\ncase=%s\nensemble_size=%s\narchive_layout=%s\ndart_root=%s\nslurm_job_id=%s\ncompleted_at=%s\n' "${START_DATE}" "${START_TOD}" "${DART_CASE}" "${DART_ENSNUM}" "per_member" "${my_dart_root}" "${SLURM_JOB_ID:-none}" "$(date '+%Y-%m-%d %H:%M:%S')" > "${perturb_tmp}"
mv -f "${perturb_tmp}" "${perturb_status}" || exit 44
for member_marker in "${PERTURB_MARKERS[@]}"; do
   rm -f -- "${member_marker}" || exit 40
done
echo "Validated all ${DART_ENSNUM} perturbed ensemble members"
echo "Wrote step-3 completion record: ${perturb_status}"
echo "`date` -- END EAM_DART_PERTURBATION"

# Ensure the removal of unneeded restart sets and copy of obs_seq.final are finished.
wait

exit 0
