#!/bin/bash -el 
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --account=esmd
#SBATCH  --time=2:00:00
#SBATCH  --partition=short
#SBATCH  --job-name=e3sm_dart_diag 
#SBATCH  --nodes=1
#SBATCH  --output=e3sm_dart_diag.%j 
#SBATCH  --exclusive
#SBATCH  --no-kill
#SBATCH  --requeue

echo == Start of e3sm dart diagnostic ==
date
echo ============================================

#source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh
#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.sh
source /qfs/people/zhan391/e3sm_dart_work/code/DART/models/eam-se/work/env_mach_specific.sh

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
LAUNCHCMD='srun -N 1 -n 1'

my_wkdir=${PWD}
scomp="eam"
cd ${my_wkdir}

source ./create_and_setup_case.sh

run_closest_member=.false. # ".false."
run_obs2netcdf=.true.    # ".false."
run_obs_diag=.true.      # ".false."

# Set paths
E3SM_ROOT=${my_e3sm_code}
DART_ROOT=${my_dart_code}
DART_MODEL=${my_dart_eam}
DART_SCPTDIR=${DART_ROOT}/models/${DART_MODEL}/shell_scripts
DART_WORKDIR=${DART_ROOT}/models/${DART_MODEL}/work
BASE_OBSDIR=${my_dart_obsdir}
BASE_PHIS=${my_e3sm_topo}
BASE_SEMAPS=${my_e3sm_semap}
BASE_CSGRID=${my_e3sm_csgrid}
CASE_ROOT=${my_modelcase}

# Run options
DART_ENSNUM=${my_ensnum}
DART_CASE=${my_casename}
DART_RUNDIR="${my_dart_runpath}/dart_en"`printf "%02d" ${DART_ENSNUM}`
DART_NTASKS=${my_dart_ntask}
DART_ON_PGRID=${my_dart_pgrid}

DATA_ASSIMILATION_CYCLES=${my_dart_cycle}
DATA_ASSIMILATION_WINDOW=${my_dart_window}
DATA_ASSIMILATION_ATM=TRUE

homme_map_file="SEMapping.nc"
cs_grid_file="SEMapping_cs_grid.nc"

CUR_YMD=${my_dartymds}
CUR_TOD=${my_darttods}
CUR_DATE=( `echo ${CUR_YMD}-${CUR_TOD} | sed -e "s#-# #g"` )
CUR_YEAR=`echo "${CUR_DATE[0]}" | bc`
CUR_MONTH=`echo "${CUR_DATE[1]}" | bc`
CUR_DAY=`echo "${CUR_DATE[2]}" | bc`
CUR_SECONDS=`echo "${CUR_DATE[3]}" | bc`
CUR_HOUR=`echo "${CUR_DATE[3]}" / 3600 | bc`
echo "valid time for eam forecast cycle is $CUR_YEAR $CUR_MONTH $CUR_DAY $CUR_SECONDS (seconds)"
echo "valid time for eam forecast cycle is $CUR_YEAR $CUR_MONTH $CUR_DAY $CUR_HOUR (hours)"

CURRENT_DADIR="${DART_RUNDIR}/dart_diagnostics"
if [ ! -d ${CURRENT_DADIR} ]; then
  mkdir -p ${CURRENT_DADIR}
fi
cd ${CURRENT_DADIR}

user_dart_nl() {
  if [ -e "${1}/diag_dart_input.nml" ]; then
    ${2} ${1}/diag_dart_input.nml input.nml  || exit 10
    sed -i "/#/d;/^\!/d;/^[ ]*\!/d" input.nml
    sed -i '1,1i\WARNING: Changes to this file will be ignored. \n Edit \$DART_WORKDIR/diag_dart_input.nml instead.\n\n\n'  input.nml
  else
    echo "ERROR ... DART required file ${1}/eam_dart_input.nml not found ... ERROR"
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

user_closest_member_nl() {
# Ensure that the input.nml ensemble size matches the number of instances.
# WARNING: the output files contain ALL ensemble members ==> BIG
ex input.nml <<ex_end
  g;ens_size ;s;= .*;= ${1};
  g;num_output_state_members ;s;= .*;= ${1};
  g;num_output_obs_members ;s;= .*;= ${1};
  g;eam_use_pgrid ;s;= .*;= ${2};
  g;input_restart_file_list ;s;= .*;= \'${input_restart_file_list}\';
  g;output_file_name        ;s;= .*;= \'${output_file_name}\';
  g;difference_method       ;s;= .*;= ${difference_method};
  g;single_restart_file_in  ;s;= .*;= .false.;
  g;use_only_qtys           ;s;= .*;= ${use_only_qtys};
  wq
ex_end
}

user_obs2nc_nl() {
  ex input.nml <<ex_end
  g;ens_size ;s;= .*;= ${1};
  g;num_output_state_members ;s;= .*;= ${1};
  g;num_output_obs_members ;s;= .*;= ${1};
  g;eam_use_pgrid ;s;= .*;= ${2};
  g;obs_sequence_name ;s;= .*;= \'\';
  g;obs_sequence_list ;s;= .*;= \'${3}\';
  g;first_bin_start      ;s;= .*;= ${4};
  g;first_bin_end        ;s;= .*;= ${5};
  g;last_bin_end         ;s;= .*;= ${6};
  g;bin_interval_days    ;s;= .*;= ${7};
  g;bin_interval_seconds ;s;= .*;= ${8};
  wq
ex_end

}

user_obs_diag_nl() {
  #run obs diagnostics 
  ex input.nml <<ex_end
  g;ens_size ;s;= .*;= ${1};
  g;num_output_state_members ;s;= .*;= ${1};
  g;num_output_obs_members ;s;= .*;= ${1};
  g;eam_use_pgrid ;s;= .*;= ${2};
  g;obs_sequence_name ;s;= .*;= \'\';
  g;obs_sequence_list ;s;= .*;= \'${3}\';
  g;first_bin_center  ;s;= .*;= ${4};
  g;last_bin_center   ;s;= .*;= ${5};
  g;bin_separation    ;s;= .*;= ${6};
  g;bin_width         ;s;= .*;= ${7};
  g;time_to_skip      ;s;= .*;= ${8};
  g;trusted_obs       ;s;= .*;= ${9};
  wq
ex_end
}

cd ${CURRENT_DADIR}
############################################################
# run closest_member_tool (closest member to ensemble mean) 
############################################################
if [[ ${run_closest_member}  == '.true.' ]];then 
   
  WORKDIR="${CURRENT_DADIR}/closest_member"
  if [ ! -d "${WORKDIR}" ];then 
     mkdir ${WORKDIR}
  fi

  cd ${WORKDIR}
  echo ${BASE_PHIS}

  user_dart_nl ${DART_WORKDIR} ${COPY} ${my_task_per_node}
 
  #delete line that is deprecated by current version of DART
  MYSTRING=`grep eam_template_filename input.nml`
  MYSTRING=( `echo $MYSTRING | sed -e "s#[=,']# #g"` )
  EAMINPUT=${MYSTRING[1]}
  echo $EAMINPUT 
 
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

  ${COPY} -f ${DART_WORKDIR}/closest_member_tool    ${WORKDIR} || exit 59

  DART_DATE1=`date -d "$CUR_HOUR:00 ${CUR_YEAR}-${CUR_MONTH}-${CUR_DAY} +0 hours" +"%Y-%m-%d %H"`
  DART_DATE1=( `echo $DART_DATE1 | sed -e "s#-# #g"` )
  DART_YEAR1=`echo "${DART_DATE1[0]}" | bc`
  DART_MONTH1=`echo "${DART_DATE1[1]}" | bc`
  DART_DAY1=`echo "${DART_DATE1[2]}" | bc`
  DART_HOUR1=`echo "${DART_DATE1[3]}" | bc`
  DART_SECONDS1=`echo "${DART_DATE1[3]}" \* 3600 | bc`
  DART_DATE1=`printf "%04d" ${DART_YEAR1}``printf "%02d" ${DART_MONTH1}``printf "%02d" ${DART_DAY1}``printf "%02d" ${DART_HOUR1}`

  input_restart_file_list="eam_in.txt"
  output_file_name="closest_restart"
  for use_only_qtys in 'QTY_U_WIND_COMPONENT' 'QTY_V_WIND_COMPONENT' 'QTY_TEMPERATURE' 'QTY_SPECIFIC_HUMIDITY' 'QTY_CLOUD_LIQUID_WATER' 'QTY_CLOUD_ICE' 'QTY_SURFACE_PRESSURE';do 

    # namelist 
    user_closest_member_nl ${DART_ENSNUM} ${DART_ON_PGRID} ${input_restart_file_list} ${output_file_name} ${difference_method} ${use_only_qtys}

    k=0
    while [ ${k} -lt  ${DATA_ASSIMILATION_CYCLES} ] ; do 
      total_time=$((DATA_ASSIMILATION_WINDOW * k))
      DART_DATE2=`date -d "$CUR_HOUR:00 ${CUR_YEAR}-${CUR_MONTH}-${CUR_DAY} +${total_time} hours" +"%Y-%m-%d %H"`
      DART_DATE2=( `echo $DART_DATE2 | sed -e "s#-# #g"` )
      DART_YEAR2=`echo "${DART_DATE2[0]}" | bc`
      DART_MONTH2=`echo "${DART_DATE2[1]}" | bc`
      DART_DAY2=`echo "${DART_DATE2[2]}" | bc`
      DART_HOUR2=`echo "${DART_DATE2[3]}" | bc`
      DART_SECONDS2=`echo "${DART_DATE2[3]}" \* 3600 | bc`
      my_dartdate=`printf "%04d" ${DART_YEAR2}`-`printf "%02d" ${DART_MONTH2}`-`printf "%02d" ${DART_DAY2}`
      my_darttod=`printf "%05d" ${DART_SECONDS2}`
      # Loop over members to generate file list 
      rm -rvf ${input_restart_file_list}
      for i in `seq 1 $my_ensnum`;do
        echo === Starting member ${i} ===
        ENSTR=EN`printf "%02d" ${i}`
        CASE_NAME=${my_casename}.${ENSTR}
        ARC_DIR="${my_modeldir}/archive/rest/${my_dartdate}-${my_darttod}"
        ${LIST} ${ARC_DIR}/${CASE_NAME}.*eam.i*${my_dartdate}-${my_darttod}.nc >> ${input_restart_file_list}
        if [ $i == 1 ]; then   
          if [ -f ${EAMINPUT} ]; then 
            ${REMOVE} ${EAMINPUT}
          fi 
          echo ${ARC_DIR}/*eam.i*${my_dartdate}-${my_darttod}.nc
          ${LINK} ${ARC_DIR}/${CASE_NAME}.*eam.i*${my_dartdate}-${my_darttod}.nc ${EAMINPUT} || exit 90
        fi 
      done
      #run closest_member 
      ${LAUNCHCMD} ${WORKDIR}/closest_member_tool || exit 140

      fobs="${WORKDIR}/closest_restart"
      if [ -f ${fobs} ]; then
        if [ $k == 0 ]; then 
          cat ${fobs}  > ${WORKDIR}/${DART_CASE}.dart.e.${scomp}.${use_only_qtys}.closest_member.txt
        else      
          cat ${fobs}  >> ${WORKDIR}/${DART_CASE}.dart.e.${scomp}.${use_only_qtys}.closest_member.txt
        fi 
      fi
      k=$((k+1))
    done 
  done 
fi 

############################################################
# run obs2netcdf (observation file to netcdf) 
############################################################
if [ ${run_obs2netcdf}  = '.true.' ];then

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
  echo $EAMINPUT 

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

  ${COPY} -f ${DART_WORKDIR}/obs_seq_to_netcdf      ${OBSSEQ_DIR} || exit 55

  DART_DATE1=`date -d "$CUR_HOUR:00 ${CUR_YEAR}-${CUR_MONTH}-${CUR_DAY} +0 hours" +"%Y-%m-%d %H"`
  DART_DATE1=( `echo $DART_DATE1 | sed -e "s#-# #g"` )
  DART_YEAR1=`echo "${DART_DATE1[0]}" | bc`
  DART_MONTH1=`echo "${DART_DATE1[1]}" | bc`
  DART_DAY1=`echo "${DART_DATE1[2]}" | bc`
  DART_HOUR1=`echo "${DART_DATE1[3]}" | bc`
  DART_SECONDS1=`echo "${DART_DATE1[3]}" \* 3600 | bc`

  total_time=$((DATA_ASSIMILATION_WINDOW * DATA_ASSIMILATION_CYCLES))
  DART_DATE2=`date -d "$CUR_HOUR:00 ${CUR_YEAR}-${CUR_MONTH}-${CUR_DAY} +${total_time} hours" +"%Y-%m-%d %H"`
  DART_DATE2=( `echo $DART_DATE2 | sed -e "s#-# #g"` )
  DART_YEAR2=`echo "${DART_DATE2[0]}" | bc`
  DART_MONTH2=`echo "${DART_DATE2[1]}" | bc`
  DART_DAY2=`echo "${DART_DATE2[2]}" | bc`
  DART_HOUR2=`echo "${DART_DATE2[3]}" | bc`
  DART_SECONDS2=`echo "${DART_DATE2[3]}" \* 3600 | bc`

  DART_HOUR1e=$((DART_HOUR1 + DATA_ASSIMILATION_WINDOW))
  DART_SECINC=$((3600 * DATA_ASSIMILATION_WINDOW))
  first_bin_start="`echo ${DART_YEAR1},${DART_MONTH1},${DART_DAY1},${DART_HOUR1},0,0`"
  first_bin_end="`echo ${DART_YEAR1},${DART_MONTH1},${DART_DAY1},${DART_HOUR1e},0,0`"
  last_bin_end="`echo ${DART_YEAR2},${DART_MONTH2},${DART_DAY2},${DART_HOUR2},0,0`"
  bin_interval_days="0"
  bin_interval_seconds="${DART_SECINC}"
  obs_sequence_list="${OBSSEQ_DIR}/obs_sequence_list.txt"

  echo ${first_bin_start} # ${first_bin_start[5]} #$first_bin_end
  echo $last_bin_end 
  user_obs2nc_nl ${DART_ENSNUM} ${DART_ON_PGRID} ${obs_sequence_list} ${first_bin_start} ${first_bin_end} ${last_bin_end} ${bin_interval_days} ${bin_interval_seconds}
  
  ${LIST} ${DART_RUNDIR}/*/${DART_CASE}.dart.e.${scomp}_obs_seq_final.* > ${obs_sequence_list}

  #convert obs sequence to netcdf file 
  ${LAUNCHCMD} ${OBSSEQ_DIR}/obs_seq_to_netcdf || exit 130

  i=0
  while [ $i -le $DATA_ASSIMILATION_CYCLES ]; do
    hour=$((DATA_ASSIMILATION_WINDOW * i))
    DART_DATE=`date -d "$CUR_HOUR:00 ${CUR_YEAR}-${CUR_MONTH}-${CUR_DAY} +${hour} hours" +"%Y-%m-%d %H"`
    DART_DATE=( `echo $DART_DATE | sed -e "s#-# #g"` )
    DART_YEAR=`echo "${DART_DATE[0]}" | bc`
    DART_MONTH=`echo "${DART_DATE[1]}" | bc`
    DART_DAY=`echo "${DART_DATE[2]}" | bc`
    DART_HOUR=`echo "${DART_DATE[3]}" | bc`
    DART_SECONDS=`echo "${DART_DATE[3]}" \* 3600 | bc`
    echo "valid time of current DA cycle is $DART_YEAR $DART_MONTH $DART_DAY $DART_SECONDS (seconds)"
    echo "valid time of current DA cycle is $DART_YEAR $DART_MONTH $DART_DAY $DART_HOUR (hours)"
    START_DATE=`printf "%04d" ${DART_YEAR}`-`printf "%02d" ${DART_MONTH}`-`printf "%02d" ${DART_DAY}`
    START_TOD=`printf "%05d" ${DART_SECONDS}`

    fepoch="${OBSSEQ_DIR}/obs_epoch_`printf "%03d" ${i}`.nc"
    OUT_DATE=${START_DATE}-${START_TOD}
    if [ -f ${fepoch} ]; then
      ${MOVE} ${fepoch}  ${OBSSEQ_DIR}/${DART_CASE}.dart.e.${scomp}_obs_seq_final.${OUT_DATE}.nc
    fi
    i=$((i+1))
  done
fi 

############################################################
# run obs_diag (observation space diagnostics) 
############################################################
if [[ ${run_obs_diag}  == '.true.' ]];then

  OBSDIAG_DIR="${CURRENT_DADIR}/obs_diag"
  if [ ! -d "${OBSDIAG_DIR}" ];then
     mkdir ${OBSDIAG_DIR}
  fi

  cd ${OBSDIAG_DIR}
  user_dart_nl ${DART_WORKDIR} ${COPY} ${my_task_per_node}

  #delete line that is deprecated by current version of DART
  MYSTRING=`grep eam_template_filename input.nml`
  MYSTRING=( `echo $MYSTRING | sed -e "s#[=,']# #g"` )
  EAMINPUT=${MYSTRING[1]}
  echo $EAMINPUT 

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

  ${COPY} -f ${DART_WORKDIR}/obs_diag               ${OBSDIAG_DIR} || exit 56

  #run obs diagnostics 
  DART_DATE1=`date -d "$CUR_HOUR:00 ${CUR_YEAR}-${CUR_MONTH}-${CUR_DAY} +0 hours" +"%Y-%m-%d %H"`
  DART_DATE1=( `echo $DART_DATE1 | sed -e "s#-# #g"` )
  DART_YEAR1=`echo "${DART_DATE1[0]}" | bc`
  DART_MONTH1=`echo "${DART_DATE1[1]}" | bc`
  DART_DAY1=`echo "${DART_DATE1[2]}" | bc`
  DART_HOUR1=`echo "${DART_DATE1[3]}" | bc`
  DART_SECONDS1=`echo "${DART_DATE1[3]}" \* 3600 | bc`
  DART_DATE1=`printf "%04d" ${DART_YEAR1}``printf "%02d" ${DART_MONTH1}``printf "%02d" ${DART_DAY1}``printf "%02d" ${DART_HOUR1}`

  total_time=$((DATA_ASSIMILATION_WINDOW * DATA_ASSIMILATION_CYCLES))
  DART_DATE2=`date -d "$CUR_HOUR:00 ${CUR_YEAR}-${CUR_MONTH}-${CUR_DAY} +${total_time} hours" +"%Y-%m-%d %H"`
  DART_DATE2=( `echo $DART_DATE2 | sed -e "s#-# #g"` )
  DART_YEAR2=`echo "${DART_DATE2[0]}" | bc`
  DART_MONTH2=`echo "${DART_DATE2[1]}" | bc`
  DART_DAY2=`echo "${DART_DATE2[2]}" | bc`
  DART_HOUR2=`echo "${DART_DATE2[3]}" | bc`
  DART_SECONDS2=`echo "${DART_DATE2[3]}" \* 3600 | bc`
  DART_DATE2=`printf "%04d" ${DART_YEAR2}``printf "%02d" ${DART_MONTH2}``printf "%02d" ${DART_DAY2}``printf "%02d" ${DART_HOUR2}`

  first_bin_center="`echo ${DART_YEAR1},${DART_MONTH1},${DART_DAY1},${DART_HOUR1},0,0`"
  last_bin_center="`echo ${DART_YEAR2},${DART_MONTH2},${DART_DAY2},${DART_HOUR2},0,0`"
  bin_separation="0, 0, 0, 6, 0, 0"
  bin_width="0, 0, 0, 6, 0, 0"
  time_to_skip="0, 0, 1, 0, 0, 0"
  trusted_obs="'RADIOSONDE_TEMPERATURE', 'RADIOSONDE_Q_WIND_COMPONENT', 'RADIOSONDE_U_WIND_COMPONENT', 'RADIOSONDE_V_WIND_COMPONENT'"
  obs_sequence_list="${OBSDIAG_DIR}/obs_sequence_list.txt"

  user_obs_diag_nl ${DART_ENSNUM} ${DART_ON_PGRID} ${obs_sequence_list} ${first_bin_center} ${last_bin_center} ${bin_separation} ${bin_width} ${time_to_skip} ${trusted_obs}

  ${LIST} ${DART_RUNDIR}/*/${DART_CASE}.dart.e.${scomp}_obs_seq_final.* > ${obs_sequence_list}
  #cat $obs_sequence_list

  ${LAUNCHCMD} ${OBSDIAG_DIR}/obs_diag || exit 140
  
  fout="${OBSDIAG_DIR}/obs_diag_output.nc"
  if [ -f ${fout} ]; then
    ${MOVE} ${fout}  ${OBSDIAG_DIR}/${DART_CASE}.dart.e.${scomp}_obs_diag_output.${DART_DATE1}-${DART_DATE2}.nc
  fi

fi 

echo ===== End of e3sm dart diagnostic =====
date
echo =====================================
