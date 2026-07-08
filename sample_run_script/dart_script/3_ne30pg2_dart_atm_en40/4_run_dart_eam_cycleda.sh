#!/bin/bash -el 
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --account=esmd
#SBATCH  --time=6:00:00
#SBATCH  --partition=slurm
#SBATCH  --job-name=e3sm_dart_ensda_cyc 
#SBATCH  --nodes=40
#SBATCH  --output=e3sm_dart_ensda_cyc.%j 
#SBATCH  --exclusive
#SBATCH  --no-kill
#SBATCH  --requeue

#For cshell:
#limit stacksize unlimited
#limit datasize unlimited

#For bash 
#ulimit -s unlimited
#ulimit -d unlimited

#export SLURM_NNODES=20
#export SLURM_NTASKS=800

echo == Start of 4_run_dart_eam_cycleda.sh ==
date
echo ============================================

#source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh
#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.sh
source /qfs/people/zhan391/e3sm_dart_work/code/DART/models/eam-se/work/env_mach_specific.sh

my_wkdir=${PWD}

cd ${my_wkdir}
source ./create_and_setup_case.sh

DATA_ASSIMILATION_ATM=TRUE
DATA_ASSIMILATION_CYCLES=${my_dart_cycle}
DATA_ASSIMILATION_WINDOW=${my_dart_window}

CASE_ROOT=${my_modeldir}/EN01/case_scripts
RUN_ROOT=${my_modeldir}/EN01/run
ARCHIVE_DIR="${my_modeldir}/archive"

if [ -d "${CASE_ROOT}" ]; then
  cd ${CASE_ROOT}
  if [[ ${DATA_ASSIMILATION_CYCLES} -eq 0 ]]; then 
    CUR_YMD=${my_casedate}
    CUR_TOD=${my_casetod}
  else 
    CUR_YMD=`./xmlquery RUN_STARTDATE      --value`
    CUR_TOD=`./xmlquery START_TOD          --value`
    CUR_TOD=`printf %05d ${CUR_TOD}`
  fi 
else
  echo "ERROR: Case directory does not exist: ${my_modelcase}"
  exit
fi

#determine the time for previous DA cycle 
CUR_DATE=( `echo ${CUR_YMD}-${CUR_TOD} | sed -e "s#-# #g"` )
CUR_YEAR=`echo "${CUR_DATE[0]}" | bc`
CUR_MONTH=`echo "${CUR_DATE[1]}" | bc`
CUR_DAY=`echo "${CUR_DATE[2]}" | bc`
CUR_HOUR=`echo "${CUR_DATE[3]}" / 3600 | bc`
CUR_SECONDS=`echo "${CUR_DATE[3]}" | bc`
echo "valid time for eam forecast cycle is $CUR_YEAR $CUR_MONTH $CUR_DAY $CUR_SECONDS (seconds)"

#function to modify namlist 
user_eam_nl() {
  local file="user_nl_eam"
  local ncdata_path="$1"
  local hist_freq="$2"
  local inithist_freq="${hist_freq}-HOURLY"

  if [[ ! -f "$file" ]]; then
    echo "[ERROR] Namelist file '$file' not found!"
    return 1
  fi

  ex "$file" <<EOF
g/^ *ncdata *=/s@=.*@= "${ncdata_path}"@
g/^ *inithist *=/s@=.*@= '${inithist_freq}'@
g/^ *inithist_all *=/s@=.*@= .true.@
wq
EOF
}

user_elm_nl() {
  local file="user_nl_elm"
  local finidat="$1"
  local yr_check="$2"
  local dynpft_check="$3"
  local fsurdat_check="$4"
  local pct_check="$5"

  ex "$file" <<EOF
g/^ *finidat *=/s@=.*@= "${finidat}"@
g/^ *check_finidat_year_consistency *=/s@=.*@= ${yr_check}@
g/^ *check_dynpft_consistency *=/s@=.*@= ${dynpft_check}@
g/^ *check_finidat_fsurdat_consistency *=/s@=.*@= ${fsurdat_check}@
g/^ *check_finidat_pct_consistency *=/s@=.*@= ${pct_check}@
wq
EOF
}

user_mosart_nl() {
  local file="user_nl_mosart"
  local finidat_rtm="$1"

  ex "$file" <<EOF
g/^ *finidat_rtm *=/s@=.*@= "${finidat_rtm}"@
wq
EOF
}

user_mpassi_nl() {
  local ymd="$1"
  local hour="$2"
  local ref_time="$3"
  local file="user_nl_mpassi"

  ex "$file" <<EOF
g/^ *config_start_time *=/s@=.*@= '${ymd}_${hour}'@
g/^ *config_calendar_type *=/s@=.*@= 'gregorian'@
g/^ *config_initial_condition_type *=/s@=.*@= 'restart'@
g/^ *config_do_restart *=/s@=.*@= .true.@
g/^ *config_restart_timestamp_name *=/s@=.*@= 'rpointer.ice'@
wq
EOF
}

user_mpaso_nl() {
  local ymd="$1"
  local hour="$2"
  local ref_time="$3"
  local file="user_nl_mpaso"

  ex "$file" <<EOF
g/^ *config_start_time *=/s@=.*@= '${ymd}_${hour}'@
g/^ *config_calendar_type *=/s@=.*@= 'gregorian'@
g/^ *config_do_restart *=/s@=.*@= .true.@
g/^ *config_output_reference_time *=/s@=.*@= '${ref_time}'@
g/^ *config_restart_timestamp_name *=/s@=.*@= 'rpointer.ocn'@
wq
EOF
}

# =====================================
# Customize MPAS stream files if needed
# =====================================
patch_mpaso_streams() {
echo
echo 'Modifying MPAS (OCEAN) streams files'
pushd ${1}

rline=`sed -n -e 12p streams.ocean`
rline=`echo ${rline}| sed "s/filename_template=//g"`

patch streams.ocean << EOF
--- streams.ocean
+++ streams.ocean
@@ -12,1 +12,1 @@
-                  filename_template=${rline}
+                  filename_template="${2}"
EOF

# copy to SourceMods
cp streams.ocean  ${3}/SourceMods/src.mpaso/

popd

}

patch_mpassi_streams() {
echo
echo 'Modifying MPAS streams files'
pushd ${1}

rlin1=`sed -n -e 11p streams.seaice`
rlin1=`echo ${rlin1}| sed "s/filename_template=//g"`

rlin2=`sed -n -e 38p streams.seaice`
rlin2=`echo ${rlin2}| sed "s/filename_template=//g"`

patch streams.seaice << EOF
--- streams.seaice    
+++ streams.seaice    
@@ -11,1 +11,1 @@
-                  filename_template=${rlin1}
+                  filename_template="${2}"
@@ -38,1 +38,1 @@
-                  filename_template=${rlin2}
+                  filename_template="${2}"
EOF

# copy to SourceMods
cp streams.seaice ${3}/SourceMods/src.mpassi/

popd

}

# First Step: Loop over members and run 6-hour ensembel forecast
for i in `seq 1 ${my_ensnum}`;do
  echo === Starting member ${i} ===
  ENSTR=EN`printf "%02d" ${i}`
  CASE_NAME=${my_casename}.${ENSTR}
  CASE_DIR=`echo ${CASE_ROOT} | sed "s/EN01/${ENSTR}/g"`
  RUN_DIR=`echo ${RUN_ROOT} | sed "s/EN01/${ENSTR}/g"`
  REF_DIR="${ARCHIVE_DIR}/rest/${CUR_YMD}-${CUR_TOD}"
  cd ${CASE_DIR}
  ./xmlchange run_exe="--kill-on-bad-exit=1 --job-name=${CASE_NAME} \${EXEROOT}/e3sm.exe "
  ./xmlchange RUN_TYPE="hybrid"
  ./xmlchange CONTINUE_RUN=FALSE
  ./xmlchange RUN_STARTDATE="${CUR_YMD}"
  ./xmlchange START_TOD="${CUR_TOD}"
  ./xmlchange REST_OPTION="nhours"
  ./xmlchange REST_N="${DATA_ASSIMILATION_WINDOW}"
  ./xmlchange STOP_OPTION="nhours"
  ./xmlchange STOP_N="${DATA_ASSIMILATION_WINDOW}"
  ./xmlchange GET_REFCASE=FALSE
  ./xmlchange RUN_REFCASE="${CASE_NAME}"
  ./xmlchange RUN_REFDATE="${CUR_YMD}"
  ./xmlchange RUN_REFTOD="${CUR_TOD}"
  ./xmlchange RUN_REFDIR="${REF_DIR}"
  ./xmlchange DOUT_S=True  #FALSE
  ./xmlchange DOUT_S_ROOT="${ARCHIVE_DIR}"

  atm_in="${REF_DIR}/${CASE_NAME}.eam.i.${CUR_YMD}-${CUR_TOD}.nc"
  lnd_in="${REF_DIR}/${CASE_NAME}.elm.r.${CUR_YMD}-${CUR_TOD}.nc"
  rof_in="${REF_DIR}/${CASE_NAME}.mosart.r.${CUR_YMD}-${CUR_TOD}.nc"
  ocn_in="${REF_DIR}/${CASE_NAME}.mpaso.rst.${CUR_YMD}_${CUR_TOD}.nc"
  ice_in="${REF_DIR}/${CASE_NAME}.mpassi.rst.${CUR_YMD}_${CUR_TOD}.nc"
  drv_in="${REF_DIR}/${CASE_NAME}.cpl.r.${CUR_YMD}-${CUR_TOD}.nc"

  #check and correct the pointer file 
  cd ${RUN_DIR}
  rm -rvf rpointer*
  atm_rst="${CASE_NAME}.eam.r.${CUR_YMD}-${CUR_TOD}.nc"
  lnd_rst="./${CASE_NAME}.elm.r.${CUR_YMD}-${CUR_TOD}.nc"
  rof_rst="./${CASE_NAME}.mosart.r.${CUR_YMD}-${CUR_TOD}.nc"
  drv_rst="${CASE_NAME}.cpl.r.${CUR_YMD}-${CUR_TOD}.nc"
  ice_rst="${CASE_NAME}.mpassi.rst.${CUR_YMD}_${CUR_TOD}.nc"
  echo "${atm_rst}"   >  rpointer.atm
  echo "${lnd_rst}"   >  rpointer.lnd
  echo "${rof_rst}"   >  rpointer.rof
  echo "${drv_rst}"   >  rpointer.drv
  echo "${CUR_YMD}_`printf "%02d" ${CUR_HOUR}`:00:00"  > rpointer.ice
  ######################################
  if [[ ${my_runtype} == "AMIP" ]];then 
    ocn1_rst="${CASE_NAME}.docn.r.${CUR_YMD}_${CUR_TOD}.nc"
    ocn2_rst="${CASE_NAME}.docn.rs1.${CUR_YMD}_${CUR_TOD}.bin"
    echo "${ocn1_rst}"  >  rpointer.ocn
    echo "${ocn2_rst}"  >> rpointer.ocn
    if [ ! -f ${ocn1_rst} ] && [ -f ${REF_DIR}/${ocn1_rst} ] ;then 
      cp -rp ${REF_DIR}/${ocn1_rst} ${ocn1_rst}
    fi 
    if [ ! -f ${ocn2_rst} ] && [ -f ${REF_DIR}/${ocn2_rst} ] ;then
      cp -rp ${REF_DIR}/${ocn2_rst} ${ocn2_rst}
    fi 
  else 
    ocn_rst="${CASE_NAME}.mpaso.rst.${CUR_YMD}_${CUR_TOD}.nc"
    echo "${CUR_YMD}_`printf "%02d" ${CUR_HOUR}`:00:00"  > rpointer.ocn
    if [ ! -f ${ocn_rst} ]; then 
      cp -rp ${REF_DIR}/${ocn_rst} ${ocn_rst}
    fi 
  fi
  for file in ${atm_rst} ${lnd_rst} ${rof_rst} ${drv_rst} ${ice_rst};do 
    fname=`basename ${file}`
    if [ ! -f ${file} ]; then 
      cp -rp ${REF_DIR}/${fname} ${file}
    fi 
  done 

  OCN_DATE=${CUR_YMD}
  OCN_HOUR=`echo "${CUR_TOD}" / 3600 | bc`
  OCN_TIME="${OCN_DATE}_"`printf "%02d" ${OCN_HOUR}`":00:00"

  #############################
  # revise namelist 
  #############################
  cd ${CASE_DIR}
  user_eam_nl ${atm_in} ${DATA_ASSIMILATION_WINDOW}  
  if [[ ${DATA_ASSIMILATION_CYCLES} -eq 0 ]];then
    user_elm_nl ${lnd_in} .false. .false. .false. .false.
  else 
    user_elm_nl ${lnd_in} .false. .false. .true. .true.
  fi 
  user_mosart_nl ${rof_in} 
  user_mpassi_nl ${OCN_DATE} ${OCN_HOUR} ${OCN_TIME}
  if [[ ${my_runtype} == "Full-CPL" ]];then
    user_mpaso_nl ${OCN_DATE} ${OCN_HOUR} ${OCN_TIME}
  fi
 
  # Finally, run CIME case.setup
  ./case.setup

  ##############################
  # Patch mpas streams files
  ##############################
  cd ${CASE_DIR}
  patch_mpassi_streams ${RUN_DIR} ${ice_in} ${CASE_DIR}
  if [[ ${my_runtype} == "Full-CPL" ]];then
    patch_mpaso_streams ${RUN_DIR} ${ocn_in} ${CASE_DIR}
  fi

  ##################
  #run model 
  ##################
  cd ${CASE_DIR}
  ./case.submit --no-batch 2>&1 > e3sm.log.o${SLURM_JOB_ID} &
  PID=$!
  kpaulse="$(echo "scale=1; ( $i * 2 ) / ${my_job_nnodes}" | bc)"
  if [[ $kpaulse == "1.0" || $kpaulse == "2.0" ||  \
        $kpaulse == "3.0" || $kpaulse == "4.0" ||  \
        $kpaulse == "5.0" || $kpaulse == "6.0" ||  \
        $kpaulse == "7.0" || $kpaulse == "8.0" ]] ; then
    wait ${PID}
  fi
  echo ${PID}

  echo ============================

done

# Wait loop with external hook
while true
do

  sleep 60

  # Execute extra instructions
  cd ${my_wkdir}
  . ./eam_boundle_extra.sh

  # List running background processes.
  # (Needed for the stop clause below to work)
  k=$((k+1))
  if (( k % 5 == 0 ))
  then
    echo ============================
    date
    jobs -l
    echo ----------------------------
    squeue --job=${SLURM_JOBID} --steps
    echo ============================
  fi

  # Stop when all processes are done
  n=`jobs -l | wc -l`
  if (( n == 0 ))
  then
     echo ============================
     date
     echo No running jobs left
     echo ============================
     break
  fi
done

# Second Step: run eam-dart data assimilation 
#-----------------------------------------------------------------------------------#
# Alway assume we perform forecast first, then DART data assimilation. At begining 
# of each forecast, we link the DART modified IC/BC files to the run directory to 
# enable the forecast for cycling data assimilation 
#-----------------------------------------------------------------------------------# 
DART_DATE=`date -d "$CUR_HOUR:00 ${CUR_YEAR}-${CUR_MONTH}-${CUR_DAY} + ${DATA_ASSIMILATION_WINDOW} hours" +"%Y-%m-%d %H"`
DART_DATE=( `echo $DART_DATE | sed -e "s#-# #g"` )
DART_YEAR=`echo "${DART_DATE[0]}" | bc`
DART_MONTH=`echo "${DART_DATE[1]}" | bc`
DART_DAY=`echo "${DART_DATE[2]}" | bc`
DART_HOUR=`echo "${DART_DATE[3]}" | bc`
DART_SECONDS=`echo "${DART_DATE[3]}" \* 3600 | bc`
echo "valid time of current DA cycle is $DART_YEAR $DART_MONTH $DART_DAY $DART_SECONDS (seconds)"

# check if needed data is ready for data assimilation 
for i in `seq 1 $my_ensnum`;do
  ENSTR=EN`printf "%02d" ${i}`
  CASE_NAME=${my_casename}.${ENSTR}
  TMP_DATE=`printf "%04d" ${DART_YEAR}`-`printf "%02d" ${DART_MONTH}`-`printf "%02d" ${DART_DAY}`
  TMP_TOD=`printf "%05d" ${DART_SECONDS}`
  atm_in="${REF_DIR}/${CASE_NAME}.eam.i.${TMP_DATE}-${TMP_TOD}.nc"
  if [ ! -f "${atm_in}" ];then
     CASE_DIR=`echo ${CASE_ROOT} | sed "s/EN01/${ENSTR}/g"`
     RUN_DIR=`echo ${RUN_ROOT} | sed "s/EN01/${ENSTR}/g"`
     atm_in1="${RUN_DIR}/${CASE_NAME}.eam.i.${TMP_DATE}-${TMP_TOD}.nc"
     if [ -f "${atm_in1}" ]; then 
       srun -N 1 -n 1 ${CASE_DIR}/case.st_archive
     else
       cd ${CASE_DIR}
       ./case.submit --no-batch 2>&1 > e3sm.log.o${SLURM_JOB_ID} &
       PID=$!
       echo "run model again to obtain ${atm_in1}"
     fi 
  fi
done
wait

# call eam-dart data assimilation 
cd ${my_wkdir}
. ./eam_dart_assimilation.sh
# Wait for dart data assimilation to complete 
wait

DATA_ASSIMILATION_CYCLES=$((DATA_ASSIMILATION_CYCLES+1))

# Post steps
cd ${my_wkdir}
. ./eam_boundle_cycling.sh 

# That's all folks!
sleep 10

echo ===== End of 4_run_dart_eam_cycleda.sh =====
date
echo =====================================
