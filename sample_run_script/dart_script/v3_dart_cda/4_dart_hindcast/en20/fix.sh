#!/bin/bash -el 
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --account=esmd
#SBATCH  --time=02:00:00
#SBATCH  --partition=short
#SBATCH  --job-name=e3sm_hindcast_run 
#SBATCH  --nodes=1
#SBATCH  --output=e3sm_hindcast_run.%j 
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

echo == Start of 2_run_e3sm_hindcast.sh ==
date
echo ============================================

#source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh
#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.sh
source /qfs/people/zhan391/e3sm_dart_work/code/DART/models/eam-se/work/env_mach_specific.sh

my_wkdir=${PWD}
my_enstart=${my_enstart}

cd ${my_wkdir}
source ./create_and_setup_case.sh

CASE_ROOT=${my_modeldir}/EN01/case_scripts
RUN_ROOT=${my_modeldir}/EN01/run
ARCHIVE_DIR="${my_modeldir}/archive"

if [ -d "${CASE_ROOT}" ]; then
  cd ${CASE_ROOT}
  CUR_YMD=${my_casedate}
  CUR_TOD=${my_casetod}
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

  ex "$file" <<EOF
g/^ *ncdata *=/s@=.*@= "${ncdata_path}"@
g/^ *inithist *=/s@=.*@= 'MONTHLY'@
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
  local year="$1"
  local month=$(printf "%02d" "$2")
  local file="user_nl_mpassi"

  ex "$file" <<EOF
g/^ *config_calendar_type *=/s@=.*@= 'gregorian'@
g/^ *config_start_time *=/s@=.*@= '${year}_${month}'@
g/^ *config_initial_condition_type *=/s@=.*@= 'restart'@
g/^ *config_do_restart *=/s@=.*@= .true.@
g/^ *config_restart_timestamp_name *=/s@=.*@= 'rpointer.ice'@
wq
EOF
}

user_mpaso_nl() {
  local year="$1"
  local month=$(printf "%02d" "$2")
  local ref_time="$3"
  local file="user_nl_mpaso"

  ex "$file" <<EOF
g/^ *config_calendar_type *=/s@=.*@= 'gregorian'@
g/^ *config_do_restart *=/s@=.*@= .true.@
g/^ *config_start_time *=/s@=.*@= '${year}_${month}'@
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

# Loop over members and run ensembel forecast
for i in `seq 1 8`;do
  echo === Starting member ${i} ===
  ENSTR=EN`printf "%02d" ${i}`
  CASE_NAME=${my_casename}.${ENSTR}
  CASE_DIR=`echo ${CASE_ROOT} | sed "s/EN01/${ENSTR}/g"`
  RUN_DIR=`echo ${RUN_ROOT} | sed "s/EN01/${ENSTR}/g"`
  REF_NAME=${my_refcase}.${ENSTR}
  REF_DIR="${my_refdir}/rest/${CUR_YMD}-${CUR_TOD}"
  cd ${CASE_DIR}
  ./xmlchange run_exe="--kill-on-bad-exit=1 --job-name=${CASE_NAME} \${EXEROOT}/e3sm.exe "
  ./xmlchange RUN_TYPE="hybrid"
  ./xmlchange CONTINUE_RUN=FALSE
  ./xmlchange RUN_STARTDATE="${CUR_YMD}"
  ./xmlchange START_TOD="${CUR_TOD}"
  ./xmlchange REST_OPTION="nmonths"
  ./xmlchange REST_N="1"
  ./xmlchange STOP_OPTION="nmonths"
  ./xmlchange STOP_N="2"
  ./xmlchange GET_REFCASE=FALSE
  ./xmlchange RUN_REFCASE="${REF_NAME}"
  ./xmlchange RUN_REFDATE="${CUR_YMD}"
  ./xmlchange RUN_REFTOD="${CUR_TOD}"
  ./xmlchange RUN_REFDIR="${REF_DIR}"
  ./xmlchange DOUT_S=True  #FALSE
  ./xmlchange DOUT_S_ROOT="${ARCHIVE_DIR}"

  atm_in="${REF_DIR}/${REF_NAME}.eam.i.${CUR_YMD}-${CUR_TOD}.nc"
  lnd_in="${REF_DIR}/${REF_NAME}.elm.r.${CUR_YMD}-${CUR_TOD}.nc"
  rof_in="${REF_DIR}/${REF_NAME}.mosart.r.${CUR_YMD}-${CUR_TOD}.nc"
  ocn_in="${REF_DIR}/${REF_NAME}.mpaso.rst.${CUR_YMD}_${CUR_TOD}.nc"
  ice_in="${REF_DIR}/${REF_NAME}.mpassi.rst.${CUR_YMD}_${CUR_TOD}.nc"
  drv_in="${REF_DIR}/${REF_NAME}.cpl.r.${CUR_YMD}-${CUR_TOD}.nc"

  #check and correct the pointer file 
  cd ${RUN_DIR}
  rm -rvf rpointer*
  atm_rst="${REF_NAME}.eam.r.${CUR_YMD}-${CUR_TOD}.nc"
  lnd_rst="./${REF_NAME}.elm.r.${CUR_YMD}-${CUR_TOD}.nc"
  rof_rst="./${REF_NAME}.mosart.r.${CUR_YMD}-${CUR_TOD}.nc"
  drv_rst="${REF_NAME}.cpl.r.${CUR_YMD}-${CUR_TOD}.nc"
  ice_rst="${REF_NAME}.mpassi.rst.${CUR_YMD}_${CUR_TOD}.nc"
  echo "${atm_rst}"   >  rpointer.atm
  echo "${lnd_rst}"   >  rpointer.lnd
  echo "${rof_rst}"   >  rpointer.rof
  echo "${drv_rst}"   >  rpointer.drv
  echo "${CUR_YMD}_`printf "%02d" ${CUR_HOUR}`:00:00"  > rpointer.ice
  ######################################
  if [[ ${my_runtype} == "AMIP" ]];then 
    ocn1_rst="${REF_NAME}.docn.r.${CUR_YMD}_${CUR_TOD}.nc"
    ocn2_rst="${REF_NAME}.docn.rs1.${CUR_YMD}_${CUR_TOD}.bin"
    echo "${ocn1_rst}"  >  rpointer.ocn
    echo "${ocn2_rst}"  >> rpointer.ocn
    if [ ! -f ${ocn1_rst} ] && [ -f ${REF_DIR}/${ocn1_rst} ] ;then 
      ln -sf ${REF_DIR}/${ocn1_rst} ${ocn1_rst}
    fi 
    if [ ! -f ${ocn2_rst} ] && [ -f ${REF_DIR}/${ocn2_rst} ] ;then
      ln -sf ${REF_DIR}/${ocn2_rst} ${ocn2_rst}
    fi 
  else 
    ocn_rst="${REF_NAME}.mpaso.rst.${CUR_YMD}_${CUR_TOD}.nc"
    echo "${CUR_YMD}_`printf "%02d" ${CUR_HOUR}`:00:00"  > rpointer.ocn
    if [ ! -f ${ocn_rst} ]; then 
      ln -sf ${REF_DIR}/${ocn_rst} ${ocn_rst}
    fi 
  fi
  for file in ${atm_rst} ${lnd_rst} ${rof_rst} ${drv_rst} ${ice_rst};do 
    fname=`basename ${file}`
    if [ ! -f ${file} ]; then 
      ln -sf ${REF_DIR}/${fname} ${file}
    fi 
  done 

  OCN_DATE=${CUR_YMD}
  OCN_HOUR=`echo "${CUR_TOD}" / 3600 | bc`
  OCN_TIME="${OCN_DATE}_"`printf "%02d" ${OCN_HOUR}`":00:00"

  #############################
  # revise namelist 
  #############################
  cd ${CASE_DIR}

  user_eam_nl ${atm_in} 

  user_elm_nl ${lnd_in} .false. .false. .false. .false.

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
  ./case.submit 

  PID=$!
  echo ${PID}

  echo ============================

done

echo ===== End of 2_run_e3sm_hindcast.sh =====
date
echo =====================================
