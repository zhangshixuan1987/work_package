#!/bin/csh
#SBATCH -q slurm
#SBATCH -N 10
#SBATCH -A ESMD
#SBATCH --time=1:00:00 
#SBATCH --error=dart-%j.err
#SBATCH --output=dart-%j.out
##SBATCH --mail-user=elvis@nersc.gov


#module purge
#module load intel/19.0.3
#module load netcdf/4.6.3
#module load openmpi/4.0.1
#module load pnetcdf/1.9.0
#module load mkl/2019u5
#module load cmake/3.11.4
#module load nco/4.7.9
#module load cdo/1.9.6

limit coredumpsize unlimited

limit stacksize unlimited

#setenv OMP_STACKSIZE      256M
#setenv OMP_NUM_THREADS    1
#setenv MPI_TYPE_DEPTH     16
#setenv MPI_IB_CONGESTED   1
#setenv MPIEXEC_MPT_DEBUG  0
#setenv MP_DEBUG_NOTIMEOUT yes

setenv OMP_PROC_BIND true
setenv OMP_PLACES threads
setenv OMP_NUM_THREADS 4

setenv  LAUNCHCMD  "srun"

#=========================================================================
# Block 0: Set command environment
#=========================================================================
# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.

setenv  MOVEV "FALSE"

echo "`date` -- BEGIN CAM_ASSIMILATE"
pwd

set ncycle        = 0
set num_instances = 10
set case          = DENS_F20TRC5-AMIP_f19_g16_DT1800
setenv ensemble_size  $num_instances

setenv CASEDIR            /compyfs/zhan391/fens_e3sm_cam5/cases/$case
setenv SAVEDIR            /compyfs/zhan391/fens_e3sm_cam5/archive/$case
setenv DARTROOT           /compyfs/zhan391/code/DART_E3SM
setenv DART_SCRIPTS_DIR   $DARTROOT/models/E3SM/shell_scripts/e3sm_cam5
setenv BASEOBSDIR         /compyfs/zhan391/acme_init/Observations/PREBUFR
setenv SAVE_EVERY_MTH_RESTART 1

###source the environment for compy and DART run
#source /compyfs/zhan391/code/DART_E3SM/models/E3SM/environment.sh

echo "DART_SCRIPTS_DIR = $DART_SCRIPTS_DIR"
cd $CASEDIR

if ( ! -e ./xmlquery ) then
   echo "ERROR: $0 must be run from a CASEDIR directory".
   exit 1
endif

setenv CASE                      `./xmlquery CASE          --value`
setenv CASEROOT                  `./xmlquery CASEROOT      --value`
setenv EXEROOT                   `./xmlquery EXEROOT       --value`
setenv RUNDIR                    `./xmlquery RUNDIR        --value`
setenv COMPSET                   `./xmlquery COMPSET       --value`
setenv scomp                     `./xmlquery COMP_ATM      --value`
#setenv ensemble_size             `./xmlquery NINST_ATM     --value`
setenv CAM_DYCORE                `./xmlquery CAM_DYCORE    --value`
#setenv archive                   `./xmlquery DOUT_S_ROOT   --value`
setenv TOTALPES                  `./xmlquery TOTALPES      --value`
setenv CONT_RUN                  `./xmlquery CONTINUE_RUN  --value`
setenv TASKS_PER_NODE            `./xmlquery NTASKS_ESP   --value`
setenv DATA_ASSIMILATION_CYCLES  `./xmlquery DATA_ASSIMILATION_CYCLES --value`

# ==============================================================================
# Turn on the assimilation in CESM
./xmlchange DATA_ASSIMILATION_ATM=TRUE
./xmlchange DATA_ASSIMILATION_CYCLES=1

# ==============================================================================

# ==============================================================================
# Set the system commands to avoid user's aliases.
# ==============================================================================
set nonomatch       # suppress "rm" warnings if wildcard does not match anything

# The FORCE options are not optional.
# The VERBOSE options are useful for debugging though
# some systems don't like the -v option to any of the following.

set   MOVE = '/usr/bin/mv -f'
set   COPY = '/usr/bin/cp -f --preserve=timestamps'
set   LINK = '/usr/bin/ln -fs'
set   LIST = '/usr/bin/ls '
set REMOVE = '/usr/bin/rm -fvr'

echo ""

# Make a place to store inflation restarts to protect from purging until
# st_archive can make a home for them.
if (! -d ${EXEROOT}/archive/esp/hist) mkdir -p ${EXEROOT}/archive/esp/hist

foreach DIR ( $CASEROOT $DART_SCRIPTS_DIR ${EXEROOT}/archive/esp/hist)
   if ( ! -d $DIR ) then
      echo "ERROR: directory '$DIR' not found"
      echo "       In the setup script check the setting of: $VAR"
      exit 10
   endif
end

# ==============================================================================
# Make sure the DART executables exist or build them if we can't find them.
# The DART input.nml in the model directory IS IMPORTANT during this part
# because it defines what observation types are supported.
# ==============================================================================
foreach MODEL ( E3SM )
   set targetdir = $DARTROOT/models/$MODEL/work
   if ( ! -x $targetdir/filter ) then
      echo ""
      echo "WARNING: executable file 'filter' not found."
      echo "         Looking for: $targetdir/filter "
      echo "         Trying to rebuild all executables for $MODEL now ..."
      echo "         This will be incorrect, if input.nml:preprocess_nml is not correct."
      (cd $targetdir; source environment.sh)
      (cd $targetdir; ./quickbuild.csh -mpi)
      if ( ! -x $targetdir/filter ) then
         echo "ERROR: executable file 'filter' not found."
         echo "       Unsuccessfully tried to rebuild: $targetdir/filter "
         echo "       Required DART assimilation executables are not found."
         echo "       Stopping prematurely."
         exit 20
      endif
   endif
end

cd $CASEROOT
# ==============================================================================
# Stage the required parts of DART in the CASEROOT directory.
# ==============================================================================

sed -e "s#BOGUSNUMINST#$num_instances#" \
    ${DART_SCRIPTS_DIR}/no_assimilate.csh.template > no_assimilate.csh || exit 30

sed -e "s#BOGUSBASEOBSDIR#$BASEOBSDIR#"  \
    -e "s#BOGUS_save_every_Mth#$SAVE_EVERY_MTH_RESTART#" \
    ${DART_SCRIPTS_DIR}/assimilate.csh.template  > assimilate.csh  || exit 40

chmod 755 assimilate.csh
chmod 755 no_assimilate.csh

# ==============================================================================
# Stage the DART executables in the CESM execution root directory: EXEROOT
# If you recompile the DART code (maybe to support more observation types)
# we're making a script to make it easy to install new DART executables.
# ==============================================================================
cat << EndOfText >! stage_dart_files
#!/bin/sh

# Run this script in the ${CASEROOT} directory.
# This script copies over the dart executables and POSSIBLY a namelist
# to the proper directory.  If you have to update any dart executables,
# do it in the ${DARTROOT} directory and then rerun stage_dart_files.
# If an input.nml does not exist in the ${CASEROOT} directory,
# a default one will be copied into place.
#
# This script was autogenerated by $0 using the variables set in that script.

if [[ -e input.nml ]]; then
   echo "stage_dart_files: Using existing ${CASEROOT}/input.nml"
   if [[ -e input.nml.original ]]; then
      echo "input.nml.original already exists - not making another"
   else
      ${COPY} input.nml input.nml.original
   fi

elif [[ -e ${DARTROOT}/models/E3SM/work/input.nml ]]; then
   ${COPY} ${DARTROOT}/models/E3SM/work/input.nml  input.nml
   if [[ -x update_dart_namelists ]]; then
          ./update_dart_namelists
   fi
else
   echo "ERROR: stage_dart_files could not find an input.nml.  Aborting"
   exit 50
fi

${COPY} ${DARTROOT}/models/E3SM/work/filter                 ${EXEROOT}
${COPY} ${DARTROOT}/models/E3SM/work/perfect_model_obs      ${EXEROOT}
${COPY} ${DARTROOT}/models/E3SM/work/fill_inflation_restart ${EXEROOT}

exit 0

EndOfText
chmod 0755 stage_dart_files

./stage_dart_files  || exit 60

# ==============================================================================
# Ensure the DART namelists are consistent with the ensemble size,
# suggest settings for num members in the output diagnostics files, etc.
# The user is free to update these after setup and before running.
# ==============================================================================

# If we are using WACCM{-X} (i.e. WCxx or WXxx) we have preferred namelist values.
# Extract pieces of the COMPSET for choosing correct setup parameters.
# E.g. "AMIP_CAM5_CLM50%BGC_CICE%PRES_DOCN%DOM_MOSART_CISM1%NOEVOLVE_SWAV"

set comp_list = `echo $COMPSET   | sed -e "s/_/ /g"`
set waccm = "false"
set atm = `echo $comp_list[2] | sed -e "s#%# #"`
if ($#atm > 1) then
   echo $atm[2] | grep WC
   if ($status == 0) set waccm = "true"
endif

cat << EndOfText >! update_dart_namelists
#!/bin/sh

# This script makes certain namelist settings consistent with the number
# of ensemble members built by the setup script.
# This script was autogenerated by $0 using the variables set in that script.

# Ensure that the input.nml ensemble size matches the number of instances.
# WARNING: the output files contain ALL ensemble members ==> BIG

ex input.nml <<ex_end
g;ens_size ;s;= .*;= ${num_instances};
g;num_output_state_members ;s;= .*;= ${num_instances};
g;num_output_obs_members ;s;= .*;= ${num_instances};
wq
ex_end

if [[ "$waccm" = "true" ]]; then 
   list=\`grep '^[ ]*vertical_localization_coord' input.nml | sed -e "s#[=,]##g"\`
   if [[ "\$list[3]" = "SCALEHEIGHT" ]]; then
      list=`grep '^[ ]*vert_normalization_scale_height' input.nml | sed -e "s#[=,]##g"`
      if [[ "\$list[2]" != "1.5" ]]; then
         echo "WARNING!  input.nml is not using 1.5 for vert_normalization_scale_height."
         echo "          Use a different value only if you definitely want to. "
      fi
   else
      echo "WARNING!  input.nml is not using SCALEHEIGHT for vertical_localization_coord."
      echo "          SCALEHEIGHT is highly recommended for WACCM{-X}"
   fi
else
   echo "This model is not configured for WACCM"
   echo "COMPSET is $COMPSET"
fi

exit 0

EndOfText
 
chmod 0755 update_dart_namelists

./update_dart_namelists || exit 70

#=========================================================================
# Stage the files needed for SAMPLING ERROR CORRECTION - even if not
# initially requested. The file is static, small, and may be needed later.
#
# If it is requested and is not present ... it is an error.
#
# The sampling error correction is a lookup table.  A selection of common
# ensemble sizes should be found in the file named below.
# It is only needed if
# input.nml:&assim_tools_nml:sampling_error_correction = .true.,
#=========================================================================

if ( $num_instances > 1 ) then
   set SAMP_ERR_FILE = \
       ${DARTROOT}/assimilation_code/programs/gen_sampling_err_table/work/sampling_error_correction_table.nc

   if (  -e   ${SAMP_ERR_FILE} ) then
      ${COPY} ${SAMP_ERR_FILE} .
   else
      echo ""
      echo "WARNING: no sampling_error_correction_table.nc file found."
      echo "         This file is NOT needed unless you want to turn on the"
      echo "         sampling_error_correction feature in any of the models."
      echo ""

      set list = `grep sampling_error_correction input.nml | sed -e "s/[=\.,]//g`
      if ($list[2] == "true") exit 80

   endif

endif

# Python uses C indexing on loops; cycle = [0,....,$DATA_ASSIMILATION_CYCLES - 1]
# "Fix" that here, so the rest of the script isn't confusing.
@ cycle = $ncycle + 1

cd $RUNDIR

# A switch to save all the inflation files
setenv save_all_inf TRUE
if (! -d ${EXEROOT}/archive/esp/hist) mkdir -p ${EXEROOT}/archive/esp/hist

# A switch to signal how often to save the stages' ensemble members: NONE, RESTART_TIMES, ALL
# Mean and sd will always be saved.
setenv save_stages_freq RESTART_TIMES

#=========================================================================
# Block 1: Populate a run-time directory with the input needed to run DART.
#=========================================================================

echo "`date` -- BEGIN COPY BLOCK"

# Put a pared down copy (no comments) of input.nml in this assimilate_cam directory.
# The contents may change from one cycle to the next, so always start from 
# the known configuration in the CASEROOT directory.

if (  -e   ${CASEROOT}/input.nml ) then
   
   sed -e "/#/d;/^\!/d;/^[ ]*\!/d" \
       -e '1,1i\WARNING: Changes to this file will be ignored. \n Edit \$CASEROOT/input.nml instead.\n\n\n' \
       ${CASEROOT}/input.nml >! input.nml  || exit 20
else
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   exit 21
endif

echo "`date` -- END COPY BLOCK"

# If possible, use the round-robin approach to deal out the tasks.
# This facilitates using multiple nodes for the simultaneous I/O operations. 

if ($?TASKS_PER_NODE) then
   if ($#TASKS_PER_NODE > 0) then
      ${MOVE} input.nml input.nml.$$
      sed -e "s#layout.*#layout = 2#" \
          -e "s#tasks_per_node.*#tasks_per_node = $TASKS_PER_NODE#" \
          input.nml.$$ >! input.nml || exit 30
      ${REMOVE} input.nml.$$
   endif
endif

#=========================================================================
# Block 2: Identify requested output stages, to warn about redundant output.
#=========================================================================

set MYSTRING = `grep stages_to_write input.nml`
set MYSTRING = (`echo $MYSTRING | sed -e "s#[=,'\.]# #g"`)
set STAGE_input     = FALSE
set STAGE_forecast  = FALSE
set STAGE_preassim  = FALSE
set STAGE_postassim = FALSE
set STAGE_analysis  = FALSE
set STAGE_output    = FALSE

# Assemble lists of stages to write out, which are not the 'output' stage.

set stages_except_output = "{"
@ stage = 2
while ($stage <= $#MYSTRING) 
   if ($MYSTRING[$stage] == 'input')  then
      set STAGE_input = TRUE
      if ($stage > 2) set stages_except_output = "${stages_except_output},"
      set stages_except_output = "${stages_except_output}input"
   endif
   if ($MYSTRING[$stage] == 'forecast')  then
      set STAGE_forecast = TRUE
      if ($stage > 2) set stages_except_output = "${stages_except_output},"
      set stages_except_output = "${stages_except_output}forecast"
   endif
   if ($MYSTRING[$stage] == 'preassim')  then
      set STAGE_preassim = TRUE
      if ($stage > 2) set stages_except_output = "${stages_except_output},"
      set stages_except_output = "${stages_except_output}preassim"
   endif
   if ($MYSTRING[$stage] == 'postassim') then
      set STAGE_postassim = TRUE
      if ($stage > 2) set stages_except_output = "${stages_except_output},"
      set stages_except_output = "${stages_except_output}postassim"
   endif
   if ($MYSTRING[$stage] == 'analysis')  then
      set STAGE_analysis = TRUE
      if ($stage > 2) set stages_except_output = "${stages_except_output},"
      set stages_except_output = "${stages_except_output}analysis"
   endif
   if ($stage == $#MYSTRING) then
      set stages_all = "${stages_except_output}"
      if ($MYSTRING[$stage] == 'output')  then
         set STAGE_output = TRUE
         set stages_all = "${stages_all},output"
      endif
   endif
   @ stage++
end

# Add the closing }
set stages_all = "${stages_all}}"
set stages_except_output = "${stages_except_output}}"

# Checking
echo "stages_except_output = $stages_except_output"
echo "stages_all = $stages_all"
if ($STAGE_output != TRUE) then
   echo "ERROR: assimilate.csh requires that input.nml:filter_nml:stages_to_write includes stage 'output'"
   exit 40
endif

#=========================================================================
# Block 3: Preliminary clean up, which can run in the background.
#=========================================================================

#=========================================================================
# Block 4: Determine time of model state 
#=========================================================================
# ... from file name of first member
# of the form "./${CASE}.cam_${ensemble_member}.i.2000-01-06-00000.nc"
#
# Piping stuff through 'bc' strips off any preceeding zeros.
#-------------------------------------------------------------------------

set FILE = `head -n 1 rpointer.atm_0001`
set FILE = $FILE:r
set ATM_DATE_EXT = `echo $FILE:e`
set ATM_DATE     = `echo $FILE:e | sed -e "s#-# #g"`
set ATM_YEAR     = `echo $ATM_DATE[1] | bc`
set ATM_MONTH    = `echo $ATM_DATE[2] | bc`
set ATM_DAY      = `echo $ATM_DATE[3] | bc`
set ATM_SECONDS  = `echo $ATM_DATE[4] | bc`
set ATM_HOUR     = `echo $ATM_DATE[4] / 3600 | bc`

echo "valid time of model is $ATM_YEAR $ATM_MONTH $ATM_DAY $ATM_SECONDS (seconds)"
echo "valid time of model is $ATM_YEAR $ATM_MONTH $ATM_DAY $ATM_HOUR (hours)"

#-----------------------------------------------------------------------------
# Get observation sequence file ... or die right away.
# The observation file names have a time that matches the stopping time of CAM.
#-----------------------------------------------------------------------------
# Make sure the file name structure matches the obs you will be using.
# PERFECT model obs output appends .perfect to the filenames

set YYYYMM   = `printf %04d%02d                ${ATM_YEAR} ${ATM_MONTH}`
if (! -d ${BASEOBSDIR}/${YYYYMM}_6H_E3SM) then
   echo "E3SM+DART requires 6 hourly obs_seq files in directories of the form YYYYMM_6H_E3SM"
   echo "The directory ${BASEOBSDIR}/${YYYYMM}_6H_E3SM is not found.  Exiting"
   exit 60
endif

#set OBSFNAME = `printf obs_seq.%04d-%02d-%02d-%05d ${ATM_YEAR} ${ATM_MONTH} ${ATM_DAY} ${ATM_SECONDS}`
set OBSFNAME = `printf obs_seq%04d%02d%02d%02d ${ATM_YEAR} ${ATM_MONTH} ${ATM_DAY} ${ATM_HOUR}`

set OBS_FILE = ${BASEOBSDIR}/${YYYYMM}_6H_E3SM/${OBSFNAME}
echo "OBS_FILE = $OBS_FILE"

if (  -e ${OBS_FILE} ) then
   ${LINK} ${OBS_FILE} obs_seq.out
else
   echo "ERROR ... no observation file ${OBS_FILE}"
   echo "ERROR ... no observation file ${OBS_FILE}"
   exit 70
endif

#=========================================================================
# Block 5: DART INFLATION
# This block is only relevant if 'inflation' is turned on AND 
# inflation values change through time:
# filter_nml
#    inf_flavor(:)  = 2  (or 3 (or 4 for posterior))
#    inf_initial_from_restart    = .TRUE.
#    inf_sd_initial_from_restart = .TRUE.
#
# This block stages the files that contain the inflation values.
# The inflation files are essentially duplicates of the DART model state,
# which have names in the E3SM style, something like 
#    ${case}.dart.rh.${scomp}_output_priorinf_{mean,sd}.YYYY-MM-DD-SSSSS.nc
# The strategy is to use the latest such files in $rundir.
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
#=========================================================================

# If we need to run fill_inflation_restart, CAM:static_init_model() 
# always needs a caminput.nc and a cam_phis.nc for geometry information, etc.

set MYSTRING = `grep cam_template_filename input.nml`
set MYSTRING = `echo $MYSTRING | sed -e "s#[=,']# #g"`
set CAMINPUT = $MYSTRING[2]
${LINK} ${CASE}.cam_0001.i.${ATM_DATE_EXT}.nc $CAMINPUT
#${LINK} ${CASE}.cam.i.${ATM_DATE_EXT}.nc $CAMINPUT

# All of the .h0. files contain the same PHIS field, so we can link to any of them.

set hists = `${LIST} ${CASE}.cam_0001.h1.*.nc`
#set hists = `${LIST} ${CASE}.cam.h1.*.nc`
set MYSTRING = `grep cam_phis_filename input.nml`
set MYSTRING = `echo $MYSTRING | sed -e "s#[=,']# #g"`
${LINK} $hists[1] $MYSTRING[2]

# Now, actually check the inflation settings

set  MYSTRING = `grep inf_flavor input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  PRIOR_INF = $MYSTRING[2]
set  POSTE_INF = $MYSTRING[3]

set  MYSTRING = `grep inf_initial_from_restart input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`

# If no inflation is requested, the inflation restart source is ignored

if ( $PRIOR_INF == 0 ) then
   set  PRIOR_INFLATION_FROM_RESTART = ignored
   set  USING_PRIOR_INFLATION = false
else
   set  PRIOR_INFLATION_FROM_RESTART = `echo $MYSTRING[2] | tr '[:upper:]' '[:lower:]'`
   set  USING_PRIOR_INFLATION = true
endif

if ( $POSTE_INF == 0 ) then
   set  POSTE_INFLATION_FROM_RESTART = ignored
   set  USING_POSTE_INFLATION = false
else
   set  POSTE_INFLATION_FROM_RESTART = `echo $MYSTRING[3] | tr '[:upper:]' '[:lower:]'`
   set  USING_POSTE_INFLATION = true
endif

if ($USING_PRIOR_INFLATION == false ) then
   set stages_requested = 0
   if ( $STAGE_input    == TRUE ) @ stages_requested++
   if ( $STAGE_forecast == TRUE ) @ stages_requested++ 
   if ( $STAGE_preassim == TRUE ) @ stages_requested++
   if ( $stages_requested > 1 ) then
      echo " "
      echo "WARNING ! ! Redundant output is requested at multiple stages before assimilation."
      echo "            Stages 'input' and 'forecast' are always redundant."
      echo "            Prior inflation is OFF, so stage 'preassim' is also redundant. "
      echo "            We recommend requesting just 'preassim'."
      echo " "
   endif
endif

if ($USING_POSTE_INFLATION == false ) then
   set stages_requested = 0
   if ( $STAGE_postassim == TRUE ) @ stages_requested++
   if ( $STAGE_analysis  == TRUE ) @ stages_requested++ 
   if ( $STAGE_output    == TRUE ) @ stages_requested++
   if ( $stages_requested > 1 ) then
      echo " "
      echo "WARNING ! ! Redundant output is requested at multiple stages after assimilation."
      echo "            Stages 'output' and 'analysis' are always redundant."
      echo "            Posterior inflation is OFF, so stage 'postassim' is also redundant. "
      echo "            We recommend requesting just 'output'."
      echo " "
   endif
endif

# IFF we want PRIOR inflation:

if ($USING_PRIOR_INFLATION == true) then
   if ($PRIOR_INFLATION_FROM_RESTART == false) then
      
      echo "inf_flavor(1) = $PRIOR_INF, using namelist values."

   else
      # Look for the output from the previous assimilation (or fill_inflation_restart)
      # If inflation files exists, use them as input for this assimilation
      (${LIST} -rt1 *.dart.rh.${scomp}_output_priorinf_mean* | tail -n 1 >! latestfile) > & /dev/null
      (${LIST} -rt1 *.dart.rh.${scomp}_output_priorinf_sd*   | tail -n 1 >> latestfile) > & /dev/null
      set nfiles = `cat latestfile | wc -l`

      if ( $nfiles > 0 ) then

         set latest_mean = `head -n 1 latestfile`
         set latest_sd   = `tail -n 1 latestfile`
         # Need to COPY instead of link because of short-term archiver and disk management.
         ${COPY} $latest_mean input_priorinf_mean.nc
         ${COPY} $latest_sd   input_priorinf_sd.nc

      else if ($CONT_RUN == FALSE) then

         # It's the first assimilation; try to find some inflation restart files
         # or make them using fill_inflation_restart.
         # Fill_inflation_restart needs caminput.nc and cam_phis.nc for static_model_init,
         # so this staging is done in assimilate.csh (after a forecast) instead of stage_cesm_files.
   
         if (-x ${EXEROOT}/fill_inflation_restart) then

            ${EXEROOT}/fill_inflation_restart
            ${MOVE} prior_inflation_mean.nc input_priorinf_mean.nc
            ${MOVE} prior_inflation_sd.nc   input_priorinf_sd.nc

         else
            echo "ERROR: Requested PRIOR inflation restart for the first cycle."
            echo "       There are no existing inflation files available "
            echo "       and ${EXEROOT}/fill_inflation_restart is missing."
            echo "EXITING"
            exit 85
         endif

      else
         echo "ERROR: Requested PRIOR inflation restart, "
         echo '       but files *.dart.rh.${scomp}_output_priorinf_* do not exist in the $rundir.'
         echo '       If you are changing from cam_no_assimilate.csh to assimilate.csh,'
         echo '       you might be able to continue by changing CONTINUE_RUN = FALSE for this cycle,'
         echo '       and restaging the initial ensemble.'
         ${LIST} -l *inf*
         echo "EXITING"
         exit 90
      endif
   endif
else
   echo "Prior Inflation not requested for this assimilation."
endif

# POSTERIOR: We look for the 'newest' and use it - IFF we need it.

if ($USING_POSTE_INFLATION == true) then
   if ($POSTE_INFLATION_FROM_RESTART == false) then

      # we are not using an existing inflation file.
      echo "inf_flavor(2) = $POSTE_INF, using namelist values."

   else
      # Look for the output from the previous assimilation (or fill_inflation_restart).
      # (The only stage after posterior inflation.)
      (${LIST} -rt1 *.dart.rh.${scomp}_output_postinf_mean* | tail -n 1 >! latestfile) > & /dev/null
      (${LIST} -rt1 *.dart.rh.${scomp}_output_postinf_sd*   | tail -n 1 >> latestfile) > & /dev/null
      set nfiles = `cat latestfile | wc -l`

      # If one exists, use it as input for this assimilation
      if ( $nfiles > 0 ) then

         set latest_mean = `head -n 1 latestfile`
         set latest_sd   = `tail -n 1 latestfile`
         ${LINK} $latest_mean input_postinf_mean.nc
         ${LINK} $latest_sd   input_postinf_sd.nc

      else if ($CONT_RUN == FALSE) then
         # It's the first assimilation; try to find some inflation restart files
         # or make them using fill_inflation_restart.
         # Fill_inflation_restart needs caminput.nc and cam_phis.nc for static_model_init,
         # so this staging is done in assimilate.csh (after a forecast) instead of stage_cesm_files.
   
         if (-x ${EXEROOT}/fill_inflation_restart) then
            ${EXEROOT}/fill_inflation_restart
            ${MOVE} prior_inflation_mean.nc input_postinf_mean.nc
            ${MOVE} prior_inflation_sd.nc   input_postinf_sd.nc

         else
            echo "ERROR: Requested POSTERIOR inflation restart for the first cycle."
            echo "       There are no existing inflation files available "
            echo "       and ${EXEROOT}/fill_inflation_restart is missing."
            echo "EXITING"
            exit 95
         endif

      else
         echo "ERROR: Requested POSTERIOR inflation restart, " 
         echo '       but files *.dart.rh.${scomp}_output_postinf_* do not exist in the $rundir.'
         ${LIST} -l *inf*
         echo "EXITING"
         exit 100
      endif
   endif
else
   echo "Posterior Inflation not requested for this assimilation."
endif

#=========================================================================
# Block 6: Actually run the assimilation. 

# DART namelist settings required:
# &filter_nml           
#    adv_ens_command         = "no_E3SM_advance_script",
#    obs_sequence_in_name    = 'obs_seq.out'
#    obs_sequence_out_name   = 'obs_seq.final'
#    single_file_in          = .false.,
#    single_file_out         = .false.,
#    stages_to_write         = stages you want + ,'output'
#    input_state_file_list   = 'cam_init_files'
#    output_state_file_list  = 'cam_init_files',

# WARNING: the default mode of this script assumes that 
#          input_state_file_list = output_state_file_list,
#          so the CAM initial files used as input to filter will be overwritten.
#          The input model states can be preserved by requesting that stage 'forecast'
#          be output.

#=========================================================================

# In the default mode of CAM assimilations, filter gets the model state(s) 
# from CAM initial files.  This section puts the names of those files into a text file.
# The name of the text file is provided to filter in filter_nml:input_state_file_list.

# NOTE: 
# If the files in input_state_file_list are E3SM initial files (all vars and 
# all meta data), then they will end up with a different structure than 
# the non-'output', stage output written by filter ('preassim', 'postassim', etc.).  
# This can be prevented (at the cost of more disk space) by copying 
# the E3SM format initial files into the names filter will use for preassim, etc.:
#    > cp $case.cam_0001.i.$date.nc  preassim_member_0001.nc.  
#    > ... for all members
# Filter will replace the state variables in preassim_member* with updated versions, 
# but leave the other variables and all metadata unchanged.

# If filter will create an ensemble from a single state,
#    filter_nml: perturb_from_single_instance = .true.
# it's fine (and convenient) to put the whole list of files in input_state_file_list.  
# Filter will just use the first as the base to perturb.

set line = `grep input_state_file_list input.nml | sed -e "s#[=,'\.]# #g"`
echo "$line"
set input_file_list = $line[2]

#${LIST} -1 ${CASE}.cam.i.${ATM_DATE_EXT}.nc >! $input_file_list
${LIST} -1 ${CASE}.cam_[0-9][0-9][0-9][0-9].i.${ATM_DATE_EXT}.nc >! $input_file_list

# If the file names in $output_state_file_list = names in $input_state_file_list,
# then the restart file contents will be overwritten with the states updated by DART.

set line = `grep output_state_file_list input.nml | sed -e "s#[=,'\.]# #g"` 
set output_file_list = $line[2]

if ($input_file_list != $output_file_list) then
   echo "ERROR: assimilate.csh requires that input_file_list = output_file_list"
   echo "       You can probably find the data you want in stage 'forecast'."
   echo "       If you truly require separate copies of CAM's initial files"
   echo "       before and after the assimilation, see revision 12603, and note that"
   echo "       it requires changing the linking to cam_initial_####.nc, below."
   exit 105
endif

echo "`date` -- BEGIN FILTER"
${LAUNCHCMD} ${EXEROOT}/filter || exit 110
echo "`date` -- END FILTER"

#========================================================================
# Block 7: Rename the output using the E3SM file-naming convention.
#=========================================================================

# If output_state_file_list is filled with custom (E3SM) filenames,
# then 'output' ensemble members will not appear with filter's default,
# hard-wired names.  But file types output_{mean,sd} will appear and be
# renamed here.

# RMA; we don't know the exact set of files which will be written,
# so loop over all possibilities.

# Handle files with instance numbers first.
foreach FILE (`$LIST ${stages_all}_member_*.nc`)
   # split off the .nc
   set parts = `echo $FILE | sed -e "s#\.# #g"`
   # separate the pieces of the remainder
   set list = `echo $parts[1]  | sed -e "s#_# #g"`
   # grab all but the trailing 'member' and #### parts.
   @ last = $#list - 2
   # and join them back together
   set dart_file = `echo $list[1-$last] | sed -e "s# #_#g"`

   set type = "e"
   echo $FILE | grep "put"
   if ($status == 0) set type = "i"

   if ($MOVEV == FALSE) \
      echo "moving $FILE ${CASE}.${scomp}_$list[$#list].${type}.${dart_file}.${ATM_DATE_EXT}.nc"
   $MOVE           $FILE ${CASE}.${scomp}_$list[$#list].${type}.${dart_file}.${ATM_DATE_EXT}.nc
end

# Files without instance numbers need to have the scomp part of their names = "dart".
# This is because in st_archive, all files with  scomp = "cam"
# (= compname in env_archive.xml) will be st_archived using a pattern 
# which has the instance number added onto it.  {mean,sd} files don't instance numbers, 
# so they need to be archived by the "dart" section of env_archive.xml.
# But they still need to be different for each component, so include $scomp in the 
# ".dart_file" part of the file name.  Somewhat awkward and inconsistent, but effective.

# Means and standard deviation files (except for inflation).
foreach FILE (`${LIST} ${stages_all}_{mean,sd}*.nc`)
   set parts = `echo $FILE | sed -e "s#\.# #g"`

   set type = "e"
   echo $FILE | grep "put"
   if ($status == 0) set type = "i"

   ${MOVE}  $FILE ${CASE}.dart.${type}.${scomp}_$parts[1].${ATM_DATE_EXT}.nc
end

# Rename the observation file and run-time output

${MOVE} obs_seq.final ${CASE}.dart.e.${scomp}_obs_seq_final.${ATM_DATE_EXT}
${MOVE} dart_log.out                 ${scomp}_dart_log.${ATM_DATE_EXT}.out

# Rename the inflation files

# Accommodate any possible inflation files.
# The .${scomp}_ part is needed by DART to distinguish
# between inflation files from separate components in coupled assims.

foreach FILE ( `$LIST ${stages_all}_{prior,post}inf_*`)
   set parts = `echo $FILE | sed -e "s#\.# #g"`
   if ($MOVEV == FALSE ) \
      echo "Moved $FILE  $CASE.dart.rh.${scomp}_$parts[1].${ATM_DATE_EXT}.nc"
   ${MOVE}        $FILE  $CASE.dart.rh.${scomp}_$parts[1].${ATM_DATE_EXT}.nc
end

# RMA; do these files have new names?
# Handle localization_diagnostics_files
set MYSTRING = `grep 'localization_diagnostics_file' input.nml`
set MYSTRING = `echo $MYSTRING | sed -e "s#[=,']# #g"`
set MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set loc_diag = $MYSTRING[2]
if (-f $loc_diag) then
   ${MOVE}           $loc_diag  ${scomp}_${loc_diag}.dart.e.${ATM_DATE_EXT}
endif

# Handle regression diagnostics
set MYSTRING = `grep 'reg_diagnostics_file' input.nml`
set MYSTRING = `echo $MYSTRING | sed -e "s#[=,']# #g"`
set MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set reg_diag = $MYSTRING[2]
if (-f $reg_diag) then
   ${MOVE}           $reg_diag  ${scomp}_${reg_diag}.dart.e.${ATM_DATE_EXT}
endif

# RMA
# Then this script will need to feed the files in output_restart_list_file
# to the next model advance.  
# This gets the .i. or .r. piece from the E3SM format file name.
set line = `grep 0001 $output_file_list | sed -e "s#[\.]# #g"` 
set l = 1
while ($l < $#line)
   if ($line[$l] =~ ${scomp}_0001) then
      @ l++
      set file_type = $line[$l]
      break
   endif
   @ l++
end

set member = 1
while ( ${member} <= ${ensemble_size} )

   set inst_string = `printf _%04d $member`
   set ATM_INITIAL_FILENAME = ${CASE}.${scomp}${inst_string}.${file_type}.${ATM_DATE_EXT}.nc

   # TJH I dont think link -f can fail ... this may never exit as intended

   ${LINK} $ATM_INITIAL_FILENAME ${scomp}_initial${inst_string}.nc || exit 120

   @ member++

end

date
echo "`date` -- END CAM_ASSIMILATE"

# Be sure that the removal of unneeded restart sets and copy of obs_seq.final are finished.
wait

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

