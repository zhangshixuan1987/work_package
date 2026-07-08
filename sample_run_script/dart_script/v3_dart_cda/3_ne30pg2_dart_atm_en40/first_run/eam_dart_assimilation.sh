#!/bin/bash -el 
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
# This script is used to perform a single-cycle eam-dart data assimilation 
#
#*******************************************************************************
# Set paths
WORK_DIR=`pwd`
DART_ROOT=${my_dart_code}
DART_MODEL=${my_dart_eam}
DART_SCPTDIR=${DART_ROOT}/models/${DART_MODEL}/shell_scripts
DART_WORKDIR=${DART_ROOT}/models/${DART_MODEL}/work
BASE_OBSDIR=${my_dart_obsdir}
BASE_PHIS=${my_e3sm_topo} 
BASE_SEMAPS=${my_e3sm_semap}
BASE_CSGRID=${my_e3sm_csgrid}

# Run options
DART_CASE=${my_casename}
DART_ENSNUM=${my_ensnum}
DART_RUNDIR="${my_dart_runpath}/dart_en"`printf "%02d" ${DART_ENSNUM}`
DART_NTASKS=${my_dart_ntask}
DART_ON_PGRID=${my_dart_pgrid}

homme_map_file="SEMapping.nc"
cs_grid_file="SEMapping_cs_grid.nc"

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
                LAUNCHCMD="srun --mpi=pmi2 --ntasks=${DART_NTASKS} --kill-on-bad-exit -l --cpu_bind=cores -c 1 -m plane=40"
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

# ==============================================================================
# Block 0: Set command environment
# ==============================================================================
# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.

echo "`date` -- BEGIN EAM_DART_ASSIMILATION"

cd ${CASE_ROOT}
scomp=`./xmlquery COMP_ATM             --value`
#EAM_DYCORE=`./xmlquery NINST_ATM       --value`
#TOTALPES=`./xmlquery TOTALPES          --value`
#CONT_RUN=`./xmlquery CONTINUE_RUN      --value`
#CHECK_TIMING=`./xmlquery CHECK_TIMING  --value`

# ==============================================================================
# Make sure the DART executables exist or build them if we can't find them.
# The DART input.nml in the model directory IS IMPORTANT during this part
# because it defines what observation types are supported.
# ==============================================================================
targetdir=${DART_ROOT}/models/${DART_MODEL}/work
if [ ! -x ${targetdir}/filter ]; then
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
TMP_DATE=`printf "%04d" ${DART_YEAR}`-`printf "%02d" ${DART_MONTH}`-`printf "%02d" ${DART_DAY}`
TMP_TOD=`printf "%05d" ${DART_SECONDS}`
ATM_DATE_EXT=${TMP_DATE}-${TMP_TOD}
CURRENT_DADIR=${DART_RUNDIR}/${ATM_DATE_EXT} 
if [ ! -d ${CURRENT_DADIR} ]; then
  mkdir -p ${CURRENT_DADIR}
else
  ${REMOVE} ${CURRENT_DADIR}/*
fi
#=========================================================================
# Loop over members and link eam files for data assimilation
# As implemented, the input filenames are static in the DART namelists.
# We must link the new uniquely-named files to static names.
#=========================================================================
if [ ${my_ensnum} -eq 1 ]; then
  cd ${CURRENT_DADIR}
  ATM_INITIAL_FILENAME="${ARCHIVE_DIR}/rest/${ATM_DATE_EXT}/${DART_CASE}.eam.i.${ATM_DATE_EXT}.nc"
  if [ ! -f "${ATM_INITIAL_FILENAME}" ];then
    echo "ERROR: required file missing ${ATM_INITIAL_FILENAME}"
    exit
  fi
  ATM_DART_FILENAME="${DART_CASE}.eam.i.${ATM_DATE_EXT}.nc"
  if [ $LINKV == TRUE ]; then 
    echo "Linking $ATM_INITIAL_FILENAME  $ATM_DART_FILENAME"
  fi
  #TEMPSZ: comment out for production run 
  $LINK $ATM_INITIAL_FILENAME  $ATM_DART_FILENAME || exit  02
  #${COPY} $ATM_INITIAL_FILENAME  $ATM_DART_FILENAME || exit  02
else 
  for i in `seq 1 ${my_ensnum}`;do
    cd ${CURRENT_DADIR}
    ENSTR=EN`printf "%02d" ${i}`
    ATM_INITIAL_FILENAME="${ARCHIVE_DIR}/rest/${ATM_DATE_EXT}/${DART_CASE}.${ENSTR}.eam.i.${ATM_DATE_EXT}.nc"
    if [ ! -f "${ATM_INITIAL_FILENAME}" ];then
      echo "ERROR: $ATM_INITIAL_FILENAME not exist!"
      exit
    fi 
    inst_string=`printf _%04d ${i}`
    ATM_DART_FILENAME=${DART_CASE}.eam${inst_string}.i.${ATM_DATE_EXT}.nc
    #echo $ATM_INITIAL_FILENAME 
    #echo $ATM_DART_FILENAME
    if [ $LINKV == TRUE ]; then 
      echo "Linking $ATM_INITIAL_FILENAME  $ATM_DART_FILENAME"
    fi 
    #TEMPSZ: comment out for production run 
    $LINK ${ATM_INITIAL_FILENAME} ${ATM_DART_FILENAME} || exit 03
    #${COPY} ${ATM_INITIAL_FILENAME} ${ATM_DART_FILENAME} || exit 03
  done
fi 

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
#echo "valid time of DA is $ATM_YEAR $ATM_MONTH $ATM_DAY $ATM_SECONDS (seconds)"
#echo "valid time of DA is $ATM_YEAR $ATM_MONTH $ATM_DAY $ATM_HOUR (hours)"

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
#echo "valid time for previous DA cycle is $OLD_YEAR $OLD_MONTH $OLD_DAY $OLD_SECONDS (seconds)"
#echo "valid time for previous DA cycle is $OLD_YEAR $OLD_MONTH $OLD_DAY $OLD_HOUR (hours)"

#=============================================
# setup ensemble size dependent parameter
#============================================
if [[ ${DART_ENSNUM} -le 20 ]];then 
  cutoff=0.1
  inf_damping=0.6 
elif [[ ${DART_ENSNUM} -le 40 ]];then
  cutoff=0.15
  inf_damping=0.6
else 
  cutoff=0.2
  inf_damping=0.6
fi 

# ==============================================================================
# Block 2: Populate a run-time directory with the input needed to run DART.
# ==============================================================================
echo "`date` -- BEGIN COPY BLOCK"
cd ${CURRENT_DADIR}
# Put a pared down copy (no comments) of input.nml in this assimilate_eam directory.
# The contents may change from one cycle to the next, so always start from
# the known configuration in the CASE_ROOT directory.
if [ -e "${DART_WORKDIR}/eam_dart_input.nml" ]; then
  ${COPY} ${DART_WORKDIR}/eam_dart_input.nml input.nml  || exit 04
  sed -i "/#/d;/^\!/d;/^[ ]*\!/d" input.nml
  sed -i '1,1i\WARNING: Changes to this file will be ignored. \n Edit \$DART_WORKDIR/eam_dart_input.nml instead.\n\n\n'  input.nml
else
  echo "ERROR ... DART required file ${DART_WORKDIR}/eam_dart_input.nml not found ... ERROR"
  echo "ERROR ... DART required file ${DART_WORKDIR}/eam_dart_input.nml not found ... ERROR"
  exit 05
fi

#This file is possibly needed for DART code after 2025/02/01
if [ -e "${DART_WORKDIR}/qceff_table.csv" ]; then
  ${COPY} ${DART_WORKDIR}/qceff_table.csv qceff_table.csv  || exit 06
else
  echo "ERROR ... DART required file ${DART_WORKDIR}/qceff_table.csv not found ... ERROR"
  exit 07
fi

# Ensure that the input.nml ensemble size matches the number of instances.
# WARNING: the output files contain ALL ensemble members ==> BIG
ex input.nml <<ex_end
g;ens_size ;s;= .*;= ${DART_ENSNUM};
g;num_output_state_members ;s;= .*;= ${DART_ENSNUM};
g;num_output_obs_members ;s;= .*;= ${DART_ENSNUM};
g;eam_use_pgrid ;s;= .*;= ${DART_ON_PGRID};
g;cutoff ;s;= .*;= ${cutoff};
g;inf_damping ;s;= .*;= ${inf_damping}, ${inf_damping};
wq
ex_end

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
if [ -v ${my_task_per_node} ]; then
   if [ ${#my_task_per_node[@]} -gt 0 ]; then
      sed -i "s#layout.*#layout = 2#"  input.nml
      sed -i "s#tasks_per_node.*#tasks_per_node = $my_task_per_node#" input.nml
   fi
fi
#echo ${#my_task_per_node[@]}
#echo ${my_task_per_node}

#stage_dart_files 
${COPY} -f ${DART_WORKDIR}/filter                 ${CURRENT_DADIR} || exit 08
${COPY} -f ${DART_WORKDIR}/perfect_model_obs      ${CURRENT_DADIR} || exit 09
${COPY} -f ${DART_WORKDIR}/fill_inflation_restart ${CURRENT_DADIR} || exit 10

if [ ${DART_ENSNUM} -gt 1 ] ; then

   SAMP_ERR_DIR=assimilation_code/programs/gen_sampling_err_table/work
   SAMP_ERR_FILE=${DART_ROOT}/${SAMP_ERR_DIR}/sampling_error_correction_table.nc

   if [ -e ${SAMP_ERR_FILE} ]; then
      ${COPY} -f ${VERBOSE} ${SAMP_ERR_FILE} ${CURRENT_DADIR}  || exit 11
      if [ ${DART_ENSNUM} -lt 3 ] || [ ${DART_ENSNUM} -gt 200 ]; then
         echo ""
         echo "ERROR: sampling_error_correction_table.nc handles ensemble sizes 3...200."
         echo "ERROR: Yours is $DART_ENSNUM"
         echo ""
         exit 12
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
         exit 13
      fi
   fi

fi

#if [ -f ${DART_SCPTDIR}/compress.csh ]; then
#   $COPY -f ${VERBOSE} ${DART_SCPTDIR}/compress.csh ${CURRENT_DADIR} || exit 14
#else
#   echo "ERROR: no compress.csh in  ${DART_SCPTDIR}"
#   exit 15
#fi

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
   exit 16
fi

# ==============================================================================
# Block 4: Preliminary clean up, which can run in the background.
# ==============================================================================
# E3SM2_0's new archiver has a mechanism for removing restart file sets,
# which we don't need, but it runs only after the (multicycle) job finishes.
# We'd like to remove unneeded restarts as the job progresses, allowing more
# cycles to run before needing to stop to archive data.  So clean them out of
# RUNDIR, and st_archive will never have to deal with them.
# ------------------------------------------------------------------------------


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
   echo "E3SM+DART requires 6 hourly obs_seq files in directories of the form YYYYMM_6H_E3SM"
   echo "The directory ${BASE_OBSDIR}/${YYYYMM}_6H_E3SM is not found.  Exiting"
   exit 17
fi

OBSFNAME=`printf obs_seq.%04d-%02d-%02d-%05d ${ATM_YEAR} ${ATM_MONTH} ${ATM_DAY} ${ATM_SECONDS}`
OBS_FILE=${BASE_OBSDIR}/${YYYYMM}_6H_CESM/${OBSFNAME}
#echo "OBS_FILE = $OBS_FILE"

${REMOVE} obs_seq.out
if [ -e ${OBS_FILE} ]; then
   ${LINK} ${OBS_FILE} obs_seq.out || exit 18
else
   echo "ERROR ... no observation file ${OBS_FILE}"
   echo "ERROR ... no observation file ${OBS_FILE}"
   exit 19
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
${LINK} ${DART_CASE}.eam_0001.i.${ATM_DATE_EXT}.nc ${EAMINPUT} || exit 20

MYSTRING=`grep eam_phis_filename input.nml`
MYSTRING=( `echo $MYSTRING | sed -e "s#[=,']# #g"` )
EAM_PHIS=${MYSTRING[1]}
${REMOVE} ${EAM_PHIS}
${LINK} ${BASE_PHIS} ${EAM_PHIS} || exit 21

#Now, Link the grid information files 
if [ ! -f ${BASE_SEMAPS} ] && [ ! -f ${BASE_CSGRID} ]; then 
  echo "ERROR ... no mapping file ${homme_map_file}"
  echo "ERROR ... no gridinfo file ${cs_grid_file}"
  echo "ERROR ... must provide either of them"
  exit 22    
else 
  if [ -f ${BASE_SEMAPS} ]; then
    ${REMOVE} ${homme_map_file}
    ${COPY} -r ${BASE_SEMAPS} ${homme_map_file} || exit 23
  fi
  if [ -f ${BASE_CSGRID} ]; then  
    ${REMOVE} ${cs_grid_file}
    ${COPY} -r ${BASE_CSGRID} ${cs_grid_file} || exit 24
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
            exit 25
         fi
      else
         echo "ERROR: Requested PRIOR inflation restart, "
         echo "       but files *.dart.rh.${scomp}_output_priorinf_* do not exist in the ${CURRENT_DADIR}."
         echo "       If you are changing from eam_no_assimilate.csh to assimilate.csh,"
         echo "       you might be able to continue by changing CONTINUE_RUN = FALSE for this cycle,"
         echo "       and restaging the initial ensemble."
         ${LIST} -l *inf*
         echo "EXITING"
         exit 26
      fi
   fi
else
   echo "Prior Inflation not requested for this assimilation."
fi
      
# POSTERIOR: We look for the 'newest' and use it - IFF we need it.
if [ ${USING_POSTE_INFLATION} == true ] ; then
   if [ ${POSTE_INFLATION_FROM}_RESTART == false ]; then
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
         ${LINK} $latest_mean input_postinf_mean.nc || exit 27
         ${LINK} $latest_sd   input_postinf_sd.nc   || exit 28
      elif [ ${DATA_ASSIMILATION_CYCLES} -eq 0 ]; then
         # It's the first assimilation; try to find some inflation restart files
         # or make them using fill_inflation_restart.
         # Fill_inflation_restart needs eaminput.nc and eam_phis.nc for static_model_init,
         # so this staging is done in assimilate.csh (after a forecast).
         if [ -x ${CURRENT_DADIR}/fill_inflation_restart ]; then
            ${CURRENT_DADIR}/fill_inflation_restart
            ${MOVE} prior_inflation_mean.nc input_postinf_mean.nc || exit 29
            ${MOVE} prior_inflation_sd.nc   input_postinf_sd.nc   || exit 30
         else
            echo "ERROR: Requested POSTERIOR inflation restart for the first cycle."
            echo "       There are no existing inflation files available "
            echo "       and ${CURRENT_DADIR}/fill_inflation_restart is missing."
            echo "EXITING"
            exit 31
         fi
      else
         echo "ERROR: Requested POSTERIOR inflation restart, "
         echo "       but files *.dart.rh.${scomp}_output_postinf_* do not exist in the ${CURRENT_DADIR}."
         ${LIST} -l *inf*
         echo "EXITING"
         exit 32
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
   exit 33
fi

#TEMPSZ: comment out for production run 
echo "`date` -- BEGIN FILTER"
${LAUNCHCMD} ${CURRENT_DADIR}/filter || exit 34
echo "`date` -- END FILTER"

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
        ${DART_CASE}.${scomp}_${list[${#list[@]}-1]}.${type}.${dart_file}.${ATM_DATE_EXT}
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
     ${MOVE} $FILE ${DART_CASE}.dart.${type}.${scomp}_${parts[0]}.${ATM_DATE_EXT}.nc || exit 35
  done 
done 

# Rename the observation file and run-time output
${MOVE} obs_seq.final ${DART_CASE}.dart.e.${scomp}_obs_seq_final.${ATM_DATE_EXT} || exit 36
${MOVE} dart_log.out  ${scomp}_dart_log.${ATM_DATE_EXT}.out || exit 37

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
   ${MOVE} $FILE  ${DART_CASE}.dart.rh.${scomp}_${parts[0]}.${ATM_DATE_EXT}.nc || exit 38
  done 
done 

# Handle localization_diagnostics_files
MYSTRING=`grep 'localization_diagnostics_file' input.nml`
MYSTRING=`echo $MYSTRING | sed -e "s#[=,']# #g"`
MYSTRING=( `echo $MYSTRING | sed -e 's#"# #g'` )
loc_diag=$MYSTRING[1]
if [ -f $loc_diag ]; then
   ${MOVE} $loc_diag  ${scomp}_${loc_diag}.dart.e.${ATM_DATE_EXT} || exit 39
fi

# Handle regression diagnostics
MYSTRING=`grep 'reg_diagnostics_file' input.nml`
MYSTRING=`echo $MYSTRING | sed -e "s#[=,']# #g"`
MYSTRING=( `echo $MYSTRING | sed -e 's#"# #g'` ) 
reg_diag=$MYSTRING[1]
if [ -f $reg_diag ] ; then
   ${MOVE} $reg_diag  ${scomp}_${reg_diag}.dart.e.${ATM_DATE_EXT} || exit 40
fi 

# Then this script will need to feed the files in output_restart_list_file
# to the next model advance.
# This gets the .i. or .r. piece from the EAM-SE format file name.

echo "`date` -- END EAM_DART_ASSIMILATION"

# Ensure the removal of unneeded restart sets and copy of obs_seq.final are finished.
wait
#exit 41
