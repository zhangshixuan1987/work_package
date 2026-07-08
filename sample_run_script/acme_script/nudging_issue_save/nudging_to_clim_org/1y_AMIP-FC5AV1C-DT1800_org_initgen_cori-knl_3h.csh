#!/bin/csh 
date

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
set compile_model = 1  # 0 = No, >0 = Yes

# 3. 
#
set run_model  = 1   # 0 = No, >0 = Yes

#-----------------------------
set debug = 'TRUE'
set debug = 'FALSE'

set git_dir = `pwd`

set dtime = 1800

####################################################################
# Fetch code
####################################################################
setenv CCSMTAG E3SM_MAINT1.0_NDG
setenv CCSMROOT /compyfs/zhan391/code/$CCSMTAG

if ($fetch_code > 0) then
   mkdir -p $CCSMROOT:h
   cd $CCSMROOT:h
   git clone -b huiwanpnnl/atm/$my_branch_shortname git@github.com:ACME-Climate/ACME.git $CCSMTAG

endif

####################################################################
# Machine, compset, PE layout etc.
####################################################################

setenv CESM_EMAIL shixuan.zhang@pnnl.gov
setenv PROJECT   ESMD
setenv CESM_PROJ $PROJECT

setenv postCIME 1
setenv COMPSET F20TRC5-CMIP6
setenv RESOLUTION ne30_ne30 #ne30_oECv3
setenv MACH compy
setenv COMPILER intel
setenv scheduler 'PBS'

setenv NTASKS_PER_INST  1024 #2016
setenv NINST  1
setenv NTHRDS 1
setenv NCORES_PER_NODE  24

setenv WORKDIR /compyfs/zhan391

setenv CSMDATA /compyfs/inputdata/

set taskname = F20TRC5-CMIP6_NUDG_ISSUE_FNL

setenv INPUT_NUDGING /compyfs/zhan391/F20TRC5-CMIP6_NUDG_ISSUE_FNL/run
setenv NUDG_NAME  AMIP_NDG_DATA_F20TRC5-CMIP6_DT1800

set casename = AMIP_ORG_NDG_3h_${COMPSET}_DT`printf "%04d" ${dtime}`

#-----------------
#setup run information (initial conditions)
set dtime1      = 1800 
set run_refcase = AMIP_CTRL_F20TRC5-CMIP6_DT`printf "%04d" ${dtime1}`
set run_refdate = "2009-11-01"  
set run_reftod  = "00000"
set initDir     = ${INPUT_NUDGING}/${run_refcase}/
#set atm_file    = "${initDir}/${run_refcase}.cam.i.2010-01-01-00000.nc"
#set lnd_file    = "${initDir}/${run_refcase}.clm2.r.2010-01-01-00000.nc"

setenv PTMP     $WORKDIR/$taskname
setenv EXELOC   $PTMP/exe/

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
   ./xmlchange -file env_build.xml -id GMAKE_J -val '16'
   ./xmlchange -file env_build.xml -id DEBUG   -val $debug
   ./case.build

endif

#####################################################################
# Conduct simulation
#####################################################################
if ($run_model > 0) then

   cd $git_dir

   setenv CASE     $casename
   setenv CASEROOT ${PTMP}/cases/$CASE
   setenv RUNDIR   ${PTMP}/run/$CASE

###copy the required restart files to the directory for branch run#####
   if ( ! -d $RUNDIR ) then
      mkdir -p $RUNDIR
   endif
   set camiin   = "$run_refcase.cam.i.$run_refdate-$run_reftod.nc"
   set camrin   = "$run_refcase.cam.r.$run_refdate-$run_reftod.nc"
   set camrsin  = "$run_refcase.cam.rs.$run_refdate-$run_reftod.nc"
  #set camrhin  = "$run_refcase.cam.rh0.$run_refdate-$run_reftod.nc"
   set lndiin   = "$run_refcase.clm2.r.$run_refdate-$run_reftod.nc"
  #set lndrhin  = "$run_refcase.clm2.rh0.$run_refdate-$run_reftod.nc"
   set cicerin  = "$run_refcase.cice.r.$run_refdate-$run_reftod.nc"
   set cplrin   = "$run_refcase.cpl.r.$run_refdate-$run_reftod.nc"
   set docnrin  = "$run_refcase.docn.r.$run_refdate-$run_reftod.nc"
   set docnrsin = "$run_refcase.docn.rs1.$run_refdate-$run_reftod.bin"

   cp -r ${initDir}${camiin}   $RUNDIR/${camiin}
   cp -r ${initDir}${camrin}   $RUNDIR/${camrin}
   cp -r ${initDir}${camrsin}  $RUNDIR/${camrsin}
  #cp -r ${initDir}${camrhin}  $RUNDIR/${camrhin}
   cp -r ${initDir}${lndiin}   $RUNDIR/${lndiin}
  #cp -r ${initDir}${lndrhin}   $RUNDIR/${lndrhin}
   cp -r ${initDir}${cicerin}  $RUNDIR/${cicerin}
   cp -r ${initDir}${cplrin}   $RUNDIR/${cplrin}
   cp -r ${initDir}${docnrin}  $RUNDIR/${docnrin}
   cp -r ${initDir}${docnrsin} $RUNDIR/${docnrsin}

   echo "${camrin}"     >! $RUNDIR/rpointer.atm
   echo "${lndiin}"     >! $RUNDIR/rpointer.lnd
   echo "${cicerin}"    >! $RUNDIR/rpointer.ice
   echo "${cplrin}"     >! $RUNDIR/rpointer.drv
   echo "${docnrin}"    >! $RUNDIR/rpointer.ocn
   echo "${docnrsin}"   >> $RUNDIR/rpointer.ocn

   set atm_init   = "$RUNDIR/${camiin}"
   set lnd_init   = "$RUNDIR/${lndiin}"
   #set atm_init   = ${atm_file}
   #set lnd_init   = ${lnd_file}

   # Create and set up new case; no need to build model

   source ${case_setup_script}

   cd $CASEROOT

   ./xmlchange -file env_build.xml -id BUILD_COMPLETE  -val 'TRUE'
   #./xmlchange -file env_build.xml -id CAM_CONFIG_OPTS -append --val='-cosp'

   # Runtime options: edit env_run.xml
   #
   # Note: with 2048 cores (128 nodes) the L72 model integrates at 
   # the speed of roughly 1 hour wall clock time per model month.
   # Set restart cycle to 6 months to make full use of the max.
   # run time for the job size bin.
   #./xmlchange  -file env_run.xml -id SSTICE_DATA_FILENAME -val "$CSMDATA/ocn/docn7/SSTDATA/sst_ice_CMIP6_DECK_E3SM_1x1_c20180213.nc"
   ./xmlchange  -file env_run.xml -id SSTICE_DATA_FILENAME -val "/compyfs/zhan391/acme_init/csmdata/ocn/docn7/SSTDATA/sst_weekly_cdcunits_1x1_1990-01-01-2013-08-11_n3.nc"
   ./xmlchange  -file env_run.xml -id SSTICE_YEAR_ALIGN -val "1990"
   ./xmlchange  -file env_run.xml -id SSTICE_YEAR_START -val "1990"
   ./xmlchange  -file env_run.xml -id SSTICE_YEAR_END   -val "2013"

  # ./xmlchange  -file env_run.xml -id  RUN_STARTDATE  -val $start_date
   ./xmlchange  -file env_run.xml -id  REST_N         -val '3'
   ./xmlchange  -file env_run.xml -id  REST_OPTION    -val 'nmonths'

  # Set sim length so that a cycle can be finished in 8 hours
  #
   @ nmon = $dtime * 16 / 2000
   ./xmlchange  -file env_run.xml -id  STOP_N       -val '3'
   ./xmlchange  -file env_run.xml -id  STOP_OPTION  -val 'nmonths'

   @ nresub = 68 / $nmon
   ./xmlchange  -file env_run.xml -id  RESUBMIT     -val '0'

   ./xmlchange  -file env_run.xml -id  SAVE_TIMING_DIR -val $PTMP
   ./xmlchange  -file env_run.xml -id  DOUT_S       -val 'FALSE'
   ./xmlchange  -file env_run.xml -id  DOUT_L_MS    -val 'FALSE'

   ./xmlchange  -file env_batch.xml -id JOB_WALLCLOCK_TIME -val "6:00:00"
   # Properly change time step for both atm and lnd.
   
   @ ncpl = 86400 / $dtime

   ./xmlchange  -file env_run.xml -id  ATM_NCPL          -val $ncpl
   ./xmlchange  -file env_run.xml -id  CAM_NAMELIST_OPTS -val dtime=$dtime
   ./xmlchange  -file env_run.xml -id  CLM_NAMELIST_OPTS -val dtime=$dtime

   ./xmlchange  -file env_run.xml -id  PIO_TYPENAME_ATM  -val netcdf
   ./xmlchange  -file env_run.xml -id  PIO_TYPENAME_CPL  -val netcdf
   ./xmlchange  -file env_run.xml -id  PIO_TYPENAME_OCN  -val netcdf
   ./xmlchange  -file env_run.xml -id  PIO_TYPENAME_WAV  -val netcdf
   ./xmlchange  -file env_run.xml -id  PIO_TYPENAME_GLC  -val netcdf
   ./xmlchange  -file env_run.xml -id  PIO_TYPENAME_ICE  -val netcdf
   ./xmlchange  -file env_run.xml -id  PIO_TYPENAME_ROF  -val netcdf
   ./xmlchange  -file env_run.xml -id  PIO_TYPENAME_LND  -val netcdf
   ./xmlchange  -file env_run.xml -id  PIO_TYPENAME_ESP  -val netcdf

   # Namelist variables

cat <<EOF >> user_nl_cam
 ncdata               = '$atm_init'
 history_amwg         = .true.
 docosp               = .false.
 cosp_amwg            = .false.
 history_verbose      = .false.
 history_aero_optics  = .true.
 history_aerosol      = .true.
 history_clubb        = .false.
 history_budget       = .true.
!history_microphysics = .true.
 iradsw               = 2,
 iradlw               = 2,
 inithist             = 'YEARLY'
 inithist_all         = .true.
 empty_htapes         = .false.
 do_aerocom_ind3      = .true.
 nhtfrq               = 0,-24,-6,-3,-1
 mfilt                = 1,1,4,8,24

 fincl1               = 'OMEGA500','RH700',
                        'Nudge_U', 'Nudge_V','Nudge_T',

 fincl2               = 'FLNT','FLNT_d1','FSNT','FSNT_d1',
                        'FLNTC','FLNTC_d1','FSNTC','FSNTC_d1',
                        'FSNTOA','FSNTOAC','FSUTOA','FSUTOAC',
                        'SWCF','SWCF_d1','LWCF','LWCF_d1',
                        'FSNS','FSNS_d1','FSDS','FSDS_d1',
                        'QRL','QRL_d1','QRS','QRS_d1',
                        'FLUT','FLUT_d1','FLDS','FLDS_d1',
                        'PRECC','PRECL','CLDHGH','CLDLOW',
                        'CLDMED','CLDTOT','AODVIS','PS',
                        'LHFLX','SHFLX','QFLX','FLNS','FLUTC',
                        'Nudge_U', 'Nudge_V','Nudge_T',

 fincl3               = 'PS','U','V','T','Q','Z3','OMEGA',

 fincl3               = 'PS','U','V','T','Q','Z3','OMEGA',

 fincl4               = 'AODVIS','angstrm','cod','cdr','cdnc',
                        'cdnum','icnum','clt','lcc','lwp',
                        'iwp','icc','icnc','icr','LHFLX',
                        'SHFLX','OMEGA500','rh700','colrv',
                        'ccn','ccn.1bl','ccn.3bl','ptop',
                        'ttop','rwp','lwp2','iwp2','autoconv',
                        'accretn','FSUTOA_d1','FSUTOAC_d1',
                        'FSUTOA','FSUTOAC','FLUT','FLUTC',
                        'TH7001000','THE7001000','PRECL','PRECC',

 fincl5               = 'QFLX','PRECC','PRECL','LHFLX','SHFLX','FLNT',
                        'FSNT','FLNS','FSNS','TREFHT','CLDHGH','CLDLOW',
                        'CLDMED','CLDTOT','TGCLDLWP','TGCLDIWP',
                        'RELHUM','RHREFHT','QREFHT','PBLH','U10','PSL',

 avgflag_pertape(2)   = 'A',
 avgflag_pertape(3)   = 'I',
 avgflag_pertape(4)   = 'I',
 avgflag_pertape(5)   = 'I',
 rad_diag_1           = 'A:Q:H2O', 'N:O2:O2', 'N:CO2:CO2', 'A:O3:O3', 'N:N2O:N2O', 'N:CH4:CH4', 'N:CFC11:CFC11', 'N:CFC12:CFC12',
EOF

cat <<EOF >> user_nl_cam
 !.......................................................
 ! nudging
 !.......................................................
  Nudge_Model  = .True.
  Nudge_Method = 'Linear'
  Nudge_Tau    = 6.0  !nudging tau in h 
  Nudge_Loc    = .False.
  Nudge_Curr   = .False.
  Num_Slice    = 1
  Nudge_Path   = '${INPUT_NUDGING}/${NUDG_NAME}/'
  Nudge_File_Template = '${NUDG_NAME}.cam.h0.%y-%m-%d-%s.nc'
  Nudge_Times_Per_Day = 8  !! nudging input data frequency
  Model_Times_Per_Day = 48 !! should not be larger than 48 if dtime = 1800s
  Nudge_Uprof = 1
  Nudge_Ucoef = 1.0
  Nudge_Vprof = 1
  Nudge_Vcoef = 1.0
  Nudge_Tprof = 1
  Nudge_Tcoef = 1.0
  Nudge_Qprof = 0
  Nudge_Qcoef = 0.0
  Nudge_PSprof = 0
  Nudge_PScoef = 0.0
  Nudge_Beg_Year  = 2009
  Nudge_Beg_Month = 11
  Nudge_Beg_Day   = 1
  Nudge_End_Year  = 9999
  Nudge_End_Month = 1
  Nudge_End_Day   = 1
EOF

cat > user_nl_clm <<EOF
 finidat = '$lnd_init'
EOF

   # Run model 
   cd $CASEROOT
   ./case.submit

date
endif


