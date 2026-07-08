#!/bin/csh 
date

#load the required module if needed
#module purge
#module load python/2.7.9

set my_branch_shortname = "master"
#=====================================
# HOW TO USE THIS SCRIPT:
#
# 1. If the model source code has already been checked out 
#    to ~/codes/$my_branch_shortname, set fetch_code to 0 (off);
#    otherwise, set it to 1 (on) to obtain code.
#
set fetch_code    = 0   # 0 = No, >0 = Yes

# 2. If the executable has already been created under
#    $PTMP/exe/, set compile_model to 0 (off);
#    otherwise set it to 1 to build the model.
#
set compile_model = 0  # 0 = No, >0 = Yes

# 3. 
#
set run_model  = 1   # 0 = No, >0 = Yes

#-----------------------------
set debug = 'TRUE'
set debug = 'FALSE'

set git_dir = `pwd`

####################################################################
# Fetch code
####################################################################
setenv CCSMTAG   maint-1.0_radcpl_dribbling
setenv CCSMROOT  /compyfs/zhan391/code/$CCSMTAG

if ($fetch_code > 0) then
   mkdir -p $CCSMROOT:h
   cd $CCSMROOT:h
   git clone -b zhan391/atm/$my_branch_shortname git@github.com:E3SM-Project/E3SM.git $CCSMTAG
  
endif

####################################################################
# Machine, compset, PE layout etc.
####################################################################

setenv CESM_EMAIL shixuan.zhang@pnnl.gov
setenv PROJECT   ESMD
setenv CESM_PROJ $PROJECT

setenv postCIME 1
setenv COMPSET FC5AV1C-L
setenv RESOLUTION ne30_ne30 #ne30_oECv3
setenv MACH compy
setenv COMPILER intel
setenv jobQueue 'slurm'

setenv NTASKS_PER_INST  2016
setenv NINST  1
setenv NTHRDS 1
setenv NCORES_PER_NODE  64

### BASIC INFO ABOUT RUN
set start_date = 0000-01-01
set dtime      = 1800  #model time step size
set walltime   = '2:00:00'

set taskname   = Maint1.0_Rdiation_Coupling
set job_name   = v1_CTRL_5day_DT`printf "%04d" ${dtime}`
set case_name  = ${job_name}.${RESOLUTION}.${MACH}

### LENGTH OF SIMULATION, RESTARTS, AND ARCHIVING
set stop_units       = ndays
set stop_num         = 5
set restart_units    = ndays
set restart_num      = 5
set num_resubmits    = 0 #2

set do_short_term_archiving      = false

### STARTUP TYPE
set model_start_type  = initial #branch

set restart_case_name = E3SMv1_FC5_NDGDAT_DT1800.ne30_ne30.cori-knl
set restart_fileymd   = "0000-10-01"
set restart_filetod   = "00000"

### Nudging data information
set ndgdata_case_name = E3SMv1_FC5_NDGDAT_DT1800.ne30_ne30.cori-knl

### SIMULATION OPTIONS
if ( $MACH == 'cori-knl' ) then
 setenv WORKDIR         /global/cscratch1/sd/zhan391
 setenv CSMDATA         /project/projectdirs/e3sm/inputdata
 setenv INPUT_NUDGING   /global/cscratch1/sd/zhan391/$taskname/run
 setenv restart_dataDir /global/cscratch1/sd/zhan391/$taskname/run/$restart_case_name
                        #/global/cscratch1/sd/zhan391/acme_init/$restart_case_name
 setenv EMIDATA         /project/projectdirs/e3sm/inputdata/atm/cam/chem
else if ( $MACH == 'compy' ) then
  setenv WORKDIR         /compyfs/zhan391 
  setenv CSMDATA         /compyfs/inputdata
  setenv INPUT_NUDGING   /compyfs/zhan391/202010_SciDAC_simulation/run/$ndgdata_case_name
  setenv restart_dataDir /compyfs/zhan391/202010_SciDAC_simulation/run/$restart_case_name
  setenv EMIDATA          $CSMDATA//atm/cam/chem
endif

setenv PTMP     $WORKDIR/$taskname
setenv EXELOC   $PTMP/exe

#-----------------
mkdir -p $PTMP

if ($NINST > 1) then
   set execase = compile_${CCSMTAG}_${COMPSET}_${RESOLUTION}_${MACH}_${COMPILER}_${NINST}x${NTASKS_PER_INST}x${NTHRDS}bundle
else
   set execase = compile_${CCSMTAG}_${COMPSET}_${RESOLUTION}_${MACH}_${COMPILER}_${NTASKS_PER_INST}x${NTHRDS}threads
endif

if ($debug == 'TRUE') then
   set execase = ${execase}_debug
endif

setenv EXEDIR ${EXELOC}/$execase

set script_output_dir = $PTMP/scripts
if ( -e $script_output_dir )  then
   mv $script_output_dir ${script_output_dir}_bak_`date +%F-%H%M%S-%N`
endif
mkdir -p $script_output_dir

set case_setup_script = "create_and_setup_bundled_case.csh"
cp $git_dir/$case_setup_script $script_output_dir/$case_setup_script

####################################################################
# Compile model
####################################################################
if ($compile_model > 0) then

   setenv CASE     $execase
   setenv CASEROOT ${PTMP}/cases/$CASE
   setenv RUNDIR   ${PTMP}/run/$CASE

   # Create and set up new case

   source ${case_setup_script}

   # Build the model

   cd $CASEROOT

  # Build with COSP
   ./xmlchange --id CAM_CONFIG_OPTS --append --val='-cosp'
  #./xmlchange --id GMAKE_J --val '16'
   ./xmlchange --id DEBUG   --val $debug
   ./case.build

endif

#####################################################################
# Conduct simulation
#####################################################################
if ($run_model > 0) then

   cd $git_dir

   setenv CASE     $case_name
   setenv CASEROOT ${PTMP}/cases/$CASE
   setenv RUNDIR   ${PTMP}/run/$CASE

   # Create and set up new case; no need to build model

   source ${case_setup_script}

   cd $CASEROOT

   ./xmlchange --id BUILD_COMPLETE  --val 'TRUE'
   ./xmlchange --id CAM_CONFIG_OPTS --append --val='-cosp'

   # Runtime options: edit env_run.xml
   ./xmlchange  --id  RUN_STARTDATE   --val $start_date
   ./xmlchange  --id  REST_N          --val $restart_num
   ./xmlchange  --id  REST_OPTION     --val $restart_units

   ./xmlchange  --id  STOP_N          --val $stop_num
   ./xmlchange  --id  STOP_OPTION     --val $stop_units

   ./xmlchange  --id  RESUBMIT        --val $num_resubmits

   ./xmlchange  --id  SAVE_TIMING_DIR --val $PTMP
   ./xmlchange  --id  DOUT_S          --val 'FALSE'
   ./xmlchange  --id  DOUT_L_MS       --val 'FALSE'

   #================================
   # SET WALLTIME FOR CREATE_NEWCASE
   #================================
   ./xmlchange JOB_WALLCLOCK_TIME=$walltime

   #================================
   # SET BATCH SYSTEM FOR THE SIMULATION
   #================================
   ./xmlchange JOB_QUEUE=$jobQueue

   #============================================
   # RUN CONFIGURATION OPTIONS
   #============================================

   #NOTE:  This section is for making specific changes to the run options (ie env_run.xml).
   ./xmlchange  --id PIO_TYPENAME  --val "netcdf"
   ### SST files :
   ./xmlchange --id SSTICE_DATA_FILENAME --val "$CSMDATA/atm/cam/sst/sst_HadOIBl_bc_1x1_clim_c101029.nc"
   #$xmlchange_exe --id SSTICE_YEAR_ALIGN --val "1849"
   #$xmlchange_exe --id SSTICE_YEAR_START --val "1849"
   #$xmlchange_exe --id SSTICE_YEAR_END --val "2016"

   # Properly change time step for both atm and lnd.

   @ ncpl = 86400 / $dtime

   ./xmlchange  --id  ATM_NCPL          --val $ncpl
   ./xmlchange  --id  CAM_NAMELIST_OPTS --val dtime=$dtime
   ./xmlchange  --id  CLM_NAMELIST_OPTS --val dtime=$dtime

 ###copy the required restart files to the directory for branch run#####
   if ( ! -d $RUNDIR ) then
      mkdir -p $RUNDIR
   endif

 if ( $model_start_type == 'initial' ) then
  ### 'initial' run: cobble together files, with RUN_TYPE= 'startup' or 'hybrid'.
  ./xmlchange --id RUN_TYPE --val "startup"
  ./xmlchange --id CONTINUE_RUN --val "FALSE"

  set camiin  = "$restart_case_name.cam.i.$restart_fileymd-$restart_filetod.nc"
  set lndiin  = "$restart_case_name.clm2.r.$restart_fileymd-$restart_filetod.nc"
  set cicerin = "$restart_case_name.cice.r.$restart_fileymd-$restart_filetod.nc"

  set atm_init = ${restart_dataDir}/${camiin}
  set lnd_init = ${restart_dataDir}/${lndiin}

 else if ( $model_start_type == 'continue' ) then

  ### This is a standard restart.

  ./xmlchange --id CONTINUE_RUN --val "TRUE"

 else if ( $model_start_type == 'branch' || $model_start_type == 'hybrid' ) then

   set camiin  = "$restart_case_name.cam.i.$restart_fileymd-$restart_filetod.nc"
   set camrin  = "$restart_case_name.cam.r.$restart_fileymd-$restart_filetod.nc"
   set camrsin = "$restart_case_name.cam.rs.$restart_fileymd-$restart_filetod.nc"
  #set camrhin = "$restart_case_name.cam.rh0.$restart_fileymd-$restart_filetod.nc"
   set lndiin  = "$restart_case_name.clm2.r.$restart_fileymd-$restart_filetod.nc"
  #set lndrin  = "$restart_case_name.clm2.rh0.$restart_fileymd-$restart_filetod.nc"
   set cicerin = "$restart_case_name.cice.r.$restart_fileymd-$restart_filetod.nc"
   set cplrin  = "$restart_case_name.cpl.r.$restart_fileymd-$restart_filetod.nc"
   set docnrin  = "$restart_case_name.docn.rs1.$restart_fileymd-$restart_filetod.nc"
   set docnrsin = "$restart_case_name.docn.rs1.$restart_fileymd-$restart_filetod.bin"
   ln -sf ${restart_dataDir}/${camiin}  $RUNDIR/${camiin}
   ln -sf ${restart_dataDir}/${camrin}  $RUNDIR/${camrin}
   ln -sf ${restart_dataDir}/${camrsin} $RUNDIR/${camrsin}
  #ln -sf ${restart_dataDir}/${camrhin} $RUNDIR/${camrhin}
   ln -sf ${restart_dataDir}/${lndiin}  $RUNDIR/${lndiin}
  #ln -sf ${restart_dataDir}/${lndrin}  $RUNDIR/${lndrin}
   ln -sf ${restart_dataDir}/${cicerin} $RUNDIR/${cicerin}
   ln -sf ${restart_dataDir}/${cplrin}  $RUNDIR/${cplrin}
   ln -sf ${restart_dataDir}/${docnrsin} $RUNDIR/${docnrsin}

   echo "${camrin}"     >! $RUNDIR/rpointer.atm
   echo "${lndiin}"     >! $RUNDIR/rpointer.lnd
   echo "${cicerin}"    >! $RUNDIR/rpointer.ice
   echo "${cplrin}"     >! $RUNDIR/rpointer.drv
   echo "${docnrin}"    >! $RUNDIR/rpointer.ocn
   echo "${docnrsin}"   >> $RUNDIR/rpointer.ocn

   set atm_init   = "$RUNDIR/${camiin}"
   set lnd_init   = "$RUNDIR/${lndiin}"
   set ice_init   = "${cicerin}"

  ./xmlchange --id CONTINUE_RUN  --val   "FALSE"
  ./xmlchange --id RUN_TYPE      --val   $model_start_type
  ./xmlchange --id RUN_REFCASE   --val   $restart_case_name
  ./xmlchange --id RUN_REFDATE   --val   $restart_fileymd
  ./xmlchange --id RUN_REFTOD    --val   $restart_filetod
  ./xmlchange --id RUN_STARTDATE --val   $restart_fileymd
  ./xmlchange --id START_TOD     --val   $restart_filetod

 else

  echo 'ERROR: $model_start_type = '${model_start_type}' is unrecognized.   Exiting.'
  exit 380

 endif

# Namelist variables


cat <<EOF >> user_nl_cam
 radheat_cpl_opt      = 0
 srf_emis_cycle_yr    = 2000
 srf_emis_type        = 'CYCLICAL'
 ext_frc_cycle_yr     = 2000
 ext_frc_type         = 'CYCLICAL'
 rad_diag_1           = 'A:Q:H2O', 'N:O2:O2', 'N:CO2:CO2', 'A:O3:O3', 'N:N2O:N2O', 'N:CH4:CH4', 'N:CFC11:CFC11', 'N:CFC12:CFC12',
 ext_frc_specifier              = 'SO2         -> ${EMIDATA}/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so2_elev_1850-2014_c180205.nc',
         'SOAG        -> ${EMIDATA}/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_soag_elev_1850-2014_c180205.nc',
         'bc_a4       -> ${EMIDATA}/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_bc_a4_elev_1850-2014_c180205.nc',
         'num_a1      -> ${EMIDATA}/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a1_elev_1850-2014_c180205.nc',
         'num_a2      -> ${EMIDATA}/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a2_elev_1850-2014_c180205.nc',
         'num_a4      -> ${EMIDATA}/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a4_elev_1850-2014_c180205.nc',
         'pom_a4      -> ${EMIDATA}/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_pom_a4_elev_1850-2014_c180205.nc',
         'so4_a1      -> ${EMIDATA}/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a1_elev_1850-2014_c180205.nc',
         'so4_a2      -> ${EMIDATA}/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a2_elev_1850-2014_c180205.nc'
 srf_emis_specifier             = 'DMS       -> ${EMIDATA}/trop_mozart_aero/emis/DMSflux.1850-2100.1deg_latlon_conserv.POPmonthlyClimFromACES4BGC_c20160727.nc',
         'SO2       -> ${EMIDATA}/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so2_surf_1850-2014_c180205.nc',
         'bc_a4     -> ${EMIDATA}/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_bc_a4_surf_1850-2014_c180205.nc',
         'num_a1    -> ${EMIDATA}/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a1_surf_1850-2014_c180205.nc',
         'num_a2    -> ${EMIDATA}/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a2_surf_1850-2014_c180205.nc',
         'num_a4    -> ${EMIDATA}/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_num_a4_surf_1850-2014_c180205.nc',
         'pom_a4    -> ${EMIDATA}/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_pom_a4_surf_1850-2014_c180205.nc',
         'so4_a1    -> ${EMIDATA}/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a1_surf_1850-2014_c180205.nc',
         'so4_a2    -> ${EMIDATA}/trop_mozart_aero/emis/DECK_ne30/cmip6_mam4_so4_a2_surf_1850-2014_c180205.nc'
 ncdata               = '$atm_init'
 inithist             = 'ENDOFRUN'
 inithist_all         = .true.
 history_amwg         = .true.
 history_verbose      = .true.
 history_aero_optics  = .false.
 history_aerosol      = .false.
 history_clubb        = .true.
 history_budget       = .true.
!history_microphysics = .true.
 empty_htapes         = .false.
 do_aerocom_ind3      = .false.
 cosp_lite            = .false.
 docosp               = .false.
 iradsw               =   -1,
 iradlw               =   -1,
 nhtfrq               =    1,
 mfilt                =   48,
 avgflag_pertape      = 'A',
EOF

# ==============
# CLM Namelist 
# With a RUN_TYPE=hybrid the finidat is automatically specified
# ==============
cat > user_nl_clm <<EOF
 finidat = '$lnd_init'
EOF


   # Run model 
   cd $CASEROOT
   ./case.submit

date
endif


