#!/bin/bash -el
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --account=esmd
#SBATCH  --time=04:00:00
#SBATCH  --partition=slurm
#SBATCH  --job-name=e3sm_dart_diag 
#SBATCH  --nodes=1
#SBATCH  --output=runtmp/logs/e3sm_dart_diag.%j 
#SBATCH  --exclusive
#SBATCH  --no-kill
#SBATCH  --requeue

set -Eeo pipefail
echo == Start of e3sm dart diagnostic ==
date
echo ============================================

#source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh
#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.sh
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
# obs_sequence_tool is a serial utility. Multiple srun tasks would start
# independent copies in this shared directory, racing on its inputs and logs.
SERIAL_LAUNCHCMD="srun -N 1 -n 1"

my_wkdir=${PWD}
scomp="eam"
cd ${my_wkdir}

source ./create_and_setup_case.sh
[[ -n "${my_dart_env_file:-}" ]] || { echo "ERROR: my_dart_env_file is unset" >&2; exit 1; }
[[ -r "${my_dart_env_file}" ]] || { echo "ERROR: configured DART environment is not readable: ${my_dart_env_file}" >&2; exit 1; }
echo "Using configured DART machine environment: ${my_dart_env_file}"
source "${my_dart_env_file}"
my_eam_dart_diag_start="${MY_DART_DIAG_START:-${my_eam_dart_diag_start}}"
my_eam_dart_diag_end="${MY_DART_DIAG_END:-${my_eam_dart_diag_end}}"
if [[ "${STEP6_DRIVER_ACTIVE:-FALSE}" != "TRUE" ]]; then
  echo "ERROR: internal diagnostic worker; run 6_run_dart_diag.sh" >&2
  exit 1
fi

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
CASE_ROOT=${my_modelcase}

# Run options
DART_ENSNUM=${my_ensnum}
DART_CASE=${my_casename}
DART_RUNDIR=${my_eam_dart_run_dir}
DART_ON_PGRID=${my_eam_dart_pgrid}

DATA_ASSIMILATION_CYCLES=${EAM_DA_COMPLETED_CYCLES:?EAM_DA_COMPLETED_CYCLES is unset}
DATA_ASSIMILATION_WINDOW=${my_eam_dart_cycle_hours}
DATA_ASSIMILATION_ATM=TRUE

homme_map_file="SEMapping.nc"
cs_grid_file="SEMapping_cs_grid.nc"

CURRENT_DADIR="${DART_RUNDIR}/dart_diagnostics"
if [ ! -d ${CURRENT_DADIR} ]; then
  mkdir -p ${CURRENT_DADIR}
fi
cd ${CURRENT_DADIR}

dart_to_epoch () {
  local dart="$1"                # YYYY-MM-DD-SSSSS (seconds since midnight)
  local date_part=${dart%-*}
  local sod=${dart##*-}
  local hh=$(( sod / 3600 ))
  local mm=$(( (sod % 3600) / 60 ))
  local ss=$(( sod % 60 ))
  date -u -d "${date_part} $(printf '%02d:%02d:%02d' "$hh" "$mm" "$ss")" +%s
}

epoch_to_dart () {
  local epoch="$1"
  local date_part
  date_part=$(date -u -d "@$epoch" +%Y-%m-%d)
  local midnight
  midnight=$(date -u -d "${date_part} 00:00:00" +%s)
  local sod=$(( epoch - midnight ))
  printf "%s-%05d" "$date_part" "$sod"
}

user_dart_nl() {
  if [ -e "${my_eam_diag_nml}" ]; then
    ${2} "${my_eam_diag_nml}" input.nml  || exit 10
    sed -i "/#/d;/^\!/d;/^[ ]*\!/d" input.nml
    sed -i '1,1i\WARNING: Changes to this runtime file will be ignored. \n Edit workflow_lib/namelists/eam/diagnostics.nml instead.\n\n\n' input.nml
  else
    echo "ERROR ... workflow EAM diagnostics namelist ${my_eam_diag_nml} not found ... ERROR"
    exit 11
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

cd ${CURRENT_DADIR}

############################################################
# run obs_diag (observation space diagnostics) 
############################################################
OBSSEQ_DIR="${CURRENT_DADIR}/obs_seq"
if [ ! -d "${OBSSEQ_DIR}" ];then
   mkdir ${OBSSEQ_DIR}
fi
cd ${OBSSEQ_DIR}

user_dart_nl ${DART_WORKDIR} ${COPY} ${my_task_per_node}

#delete line that is deprecated by current version of DART
MYSTRING=`grep eam_template_filename input.nml`
MYSTRING=( `echo $MYSTRING | sed -e "s#[=,']# #g"` )
EAMINPUT=${MYSTRING[1]}

MYSTRING=`grep eam_phis_filename input.nml`
MYSTRING=( `echo $MYSTRING | sed -e "s#[=,']# #g"` )
EAM_PHIS=${MYSTRING[1]}
${LINK} ${BASE_PHIS} ${EAM_PHIS} || exit 100

#Now, Link the grid information files 
if [ ! -f ${BASE_SEMAPS} ] && [ ! -f ${BASE_CSGRID} ]; then
  echo "ERROR ... no mapping file ${homme_map_file}"
  echo "ERROR ... no gridinfo file ${cs_grid_file}"
  echo "ERROR ... must provide either of them"
  exit 91
else
  if [ -f ${BASE_SEMAPS} ]; then
    ${COPY} -r ${BASE_SEMAPS} ${homme_map_file} || exit 101
  fi
  if [ -f ${BASE_CSGRID} ]; then
    ${COPY} -r ${BASE_CSGRID} ${cs_grid_file} || exit 102
  fi
fi

[[ -f input.nml ]] || { echo "input.nml not found"; return 1; }

#first step: use obs_sequence_tool to convert the file 
${COPY} -f ${DART_WORKDIR}/obs_sequence_tool  ${OBSSEQ_DIR} || exit 103
${COPY} -f ${DART_WORKDIR}/obs_seq_to_netcdf  ${OBSSEQ_DIR} || exit 104

#observation sequence file to process
filename_seq_list="obs_seq"
filename_out="obs_seq.processed"
qc_metadata="''" #"'Data QC','DART quality control'"
print_every=10000
min_qc=-888888.0 #0
max_qc=-888888.0 #5
calendar='Gregorian'
print_only=".false."
edit_copies=".true."
new_copy_index="1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13"

user_obs_sequence_nl \
  "'${filename_seq_list}'" \
  "'${filename_out}'" \
  "${qc_metadata}" \
  "${min_qc}" \
  "${max_qc}" \
  "${calendar}" \
  "${print_only}" \
  "${edit_copies}" \
  "${new_copy_index}"

start_epoch=$(dart_to_epoch "$my_eam_dart_diag_start")
end_epoch=$(dart_to_epoch "$my_eam_dart_diag_end")
window_sec=$(( DATA_ASSIMILATION_WINDOW * 3600 ))

current_epoch=$start_epoch
while [[ $current_epoch -lt $end_epoch ]]; do
  current_dart=$(epoch_to_dart "$current_epoch")
  current_dart=$(date -u -d "@${current_epoch}" +%Y-%m-%d)-$(date -u -d "@${current_epoch}" +%s | awk '{printf "%05d", $1 % 86400}')
  echo "Processing window starting at $current_dart"
  obs_sequence_list="${OBSSEQ_DIR}/obs_seq"
  for file in ${DART_RUNDIR}/*/${DART_CASE}.dart.e.${scomp}_obs_seq_final.${current_dart}* ;do 
     ${LIST} ${file} > ${obs_sequence_list}
     fout=`basename ${file}`
     fout=`echo ${fout} | sed "s/obs_seq_final/obs_seq.processed/g"` 
     ${SERIAL_LAUNCHCMD} ${OBSSEQ_DIR}/obs_sequence_tool || exit 56
     mv ${filename_out} ${fout}
  done 
  current_epoch=$(( current_epoch + window_sec ))
done

echo ===== End of e3sm dart diagnostic =====
date
echo =====================================
