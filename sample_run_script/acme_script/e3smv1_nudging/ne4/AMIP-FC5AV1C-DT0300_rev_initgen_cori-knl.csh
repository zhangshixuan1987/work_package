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
set compile_model = 0  # 0 = No, >0 = Yes

# 3. 
#
set run_model  = 1   # 0 = No, >0 = Yes

#-----------------------------
set debug = 'TRUE'
set debug = 'FALSE'

set git_dir = `pwd`

set dtime = 300

####################################################################
# Fetch code
####################################################################
setenv CCSMTAG ACME_master_20180611
setenv CCSMROOT /project/projectdirs/m3089/shixuan/codes/$CCSMTAG

if ($fetch_code > 0) then
   mkdir -p $CCSMROOT:h
   cd $CCSMROOT:h
   git clone -b huiwanpnnl/atm/$my_branch_shortname git@github.com:ACME-Climate/ACME.git $CCSMTAG

endif

####################################################################
# Machine, compset, PE layout etc.
####################################################################

setenv CESM_EMAIL shixuan.zhang@pnnl.gov
setenv PROJECT   m3089 
setenv CESM_PROJ $PROJECT

setenv postCIME 1
setenv COMPSET F20TRC5-CMIP6 
setenv RESOLUTION ne30_ne30 #ne30_oECv3
setenv MACH cori-knl
setenv COMPILER intel
setenv scheduler 'PBS'

setenv NTASKS_PER_INST  384
setenv NINST  1
setenv NTHRDS 4
setenv NCORES_PER_NODE  24

setenv WORKDIR /global/cscratch1/sd/zhan391

setenv CSMDATA /project/projectdirs/acme/inputdata

setenv INPUT_NUDGING /global/cscratch1/sd/yanhp/F20TRC5-CMIP6_201806_initgen_84nodes/run

set taskname = PH7_F20TRC5-CMIP6_201806_initgen_84nodes
set casename = AMIPNUDG_${COMPSET}_rev_DT`printf "%04d" ${dtime}`

#-----------------
#setup branch run information
set dtime1      = 1800 
set run_refcase = AMIPNUDG_${COMPSET}_DT`printf "%04d" ${dtime1}`
set run_refdate = "2011-03-01"
set run_reftod  = "00000"
set initDir     = ${INPUT_NUDGING}/${run_refcase}/

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
   set camiin  = "$run_refcase.cam.i.$run_refdate-$run_reftod.nc"
   set camrin  = "$run_refcase.cam.r.$run_refdate-$run_reftod.nc"
   set camrsin = "$run_refcase.cam.rs.$run_refdate-$run_reftod.nc"
  #set camrhin = "$run_refcase.cam.rh0.$run_refdate-$run_reftod.nc"
   set lndiin  = "$run_refcase.clm2.r.$run_refdate-$run_reftod.nc"
  #set lndrin  = "$run_refcase.clm2.rh0.$run_refdate-$run_reftod.nc"
   set cicerin = "$run_refcase.cice.r.$run_refdate-$run_reftod.nc"
   set cplrin  = "$run_refcase.cpl.r.$run_refdate-$run_reftod.nc"
   set docnrin = "$run_refcase.docn.rs1.$run_refdate-$run_reftod.bin"
   cp -r ${initDir}${camiin}  $RUNDIR/${camiin}
   cp -r ${initDir}${camrin}  $RUNDIR/${camrin}
   cp -r ${initDir}${camrsin} $RUNDIR/${camrsin}
  #cp -r ${initDir}${camrhin} $RUNDIR/${camrhin}
   cp -r ${initDir}${lndiin}  $RUNDIR/${lndiin}
  #cp -r ${initDir}${lndrin}  $RUNDIR/${lndrin}
   cp -r ${initDir}${cicerin} $RUNDIR/${cicerin}
   cp -r ${initDir}${cplrin}  $RUNDIR/${cplrin}
   cp -r ${initDir}${docnrin} $RUNDIR/${docnrin}
   cp -r ${initDir}rpoint*    $RUNDIR/
   foreach  fil ( $RUNDIR/rpoint* )
      sed -i "s/2011-03-01-00000/${run_refdate}-${run_reftod}/g" $fil
   end

   set atm_init   = "$RUNDIR/${camiin}"
   set lnd_init   = "$RUNDIR/${lndiin}"

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

   ./xmlchange  -file env_run.xml -id SSTICE_DATA_FILENAME -val "$initDir/sst_ice_CMIP6_DECK_E3SM_1x1_c20180213.nc"
   ./xmlchange  -file env_run.xml -id SSTICE_YEAR_ALIGN -val "1849"
   ./xmlchange  -file env_run.xml -id SSTICE_YEAR_START -val "1849"
   ./xmlchange  -file env_run.xml -id SSTICE_YEAR_END   -val "2016"

  # ./xmlchange  -file env_run.xml -id  RUN_STARTDATE  -val $start_date
   ./xmlchange  -file env_run.xml -id  REST_N         -val '1'
   ./xmlchange  -file env_run.xml -id  REST_OPTION    -val 'nmonths'

  # Set sim length so that a cycle can be finished in 8 hours
  #
   @ nmon = $dtime * 16 / 2000
   ./xmlchange  -file env_run.xml -id  STOP_N       -val '48'
   ./xmlchange  -file env_run.xml -id  STOP_OPTION  -val 'nsteps'

   @ nresub = 68 / $nmon
   ./xmlchange  -file env_run.xml -id  RESUBMIT     -val '0'

   ./xmlchange  -file env_run.xml -id  SAVE_TIMING_DIR -val $PTMP
   ./xmlchange  -file env_run.xml -id  DOUT_S       -val 'FALSE'
   ./xmlchange  -file env_run.xml -id  DOUT_L_MS    -val 'FALSE'

   ./xmlchange  -file env_batch.xml -id JOB_WALLCLOCK_TIME -val "2:00:00"
   # Properly change time step for both atm and lnd.
   
   @ ncpl = 86400 / $dtime

   ./xmlchange  -file env_run.xml -id  ATM_NCPL          -val $ncpl
   ./xmlchange  -file env_run.xml -id  CAM_NAMELIST_OPTS -val dtime=$dtime
   ./xmlchange  -file env_run.xml -id  CLM_NAMELIST_OPTS -val dtime=$dtime


   # Namelist variables


cat <<EOF >> user_nl_cam
 cld_macmic_num_steps = 1
 se_nsplit            = 1,
 rsplit               = 1,
 ncdata               = '$atm_init'
 l_dribling_tend      = .true.
 l_dribling_uv        = .false.
 l_dribling_w         = .false.
 macmic_extra_diag    = .true.
 history_amwg         = .true.
 docosp               = .false.
 cosp_amwg            = .false.
 history_verbose      = .true.
 history_aero_optics  = .false.
 history_aerosol      = .false.
 history_clubb        = .true.
 history_budget       = .true.
!history_microphysics = .true.
 iradsw               = 2,
 iradlw               = 2,
 inithist             = 'MONTHLY'
 inithist_all         = .true.
 clubb_c14            = 1.3D0
 l_ieflx_fix          = .false.
 empty_htapes         = .false.
 pergro_test_active   = .false.
 l_aerosol_cldgrow    = .true.
 l_aerosol_cldshnk    = .true.
 l_aerosol_oldcld     = .true.
 l_aerosol_mixing     = .true.
 nhtfrq               = 1, -6
 mfilt                = 1,  1
 avgflag_pertape      = 'I',"A"
 fincl2               = 'PBLH','PS','CLOUD','PHIS','Z3','LANDFRAC','FICE','ZMDT','ZMDQ','ZMDICE','ZMDLIQ','EVAPTZM','FZSNTZM','EVSNTZM',
                        'EVAPQZM','CAPE','FREQZM','CMFMCDZM','ZMEIHEAT','PRECCDZM','ZMFLXPRC','ZMFLXSNW','ZMNTPRPD','ZMNTSNPD',
                        'PRODPREC','EVAPPREC','PRAO','PRCO','BERGO','PRCIO','PRAIO','MPDT','MPDQ','MPDLIQ','MPDICE','MPICLWPI','MPICIWPI',
                        'QRAIN','QSNOW','AQRAIN','AQSNOW','VPRCO','VPRAO'
                        'UP2_CLUBB','VP2_CLUBB','WP2_CLUBB','UPWP_CLUBB','VPWP_CLUBB','WP3_CLUBB',
                        'WPTHLP_CLUBB','WPRTP_CLUBB','RTP2_CLUBB','THLP2_CLUBB','RTPTHLP_CLUBB','RCM_CLUBB','WPRCP_CLUBB','RHO_CLUBB','RELVAR',
                        'CLOUDFRAC_CLUBB','CLOUDCOVER_CLUBB','WPTHVP_CLUBB','RVMTEND_CLUBB','STEND_CLUBB','RCMTEND_CLUBB',
                        'RIMTEND_CLUBB','UTEND_CLUBB','VTEND_CLUBB','UM_CLUBB','VM_CLUBB','THETAL','PBLH','QT','SL','ZMDLF',
                        'DPDLFLIQ','DPDLFICE','DPDLFT','CONCLD','CMELIQ','FICE_q','FICE_qice','FICE_T','FICE_f'
! clubb_c14 = 1.06D0
! l_ieflx_fix = .false.
! clubb_c14            = 1.3D0 (default in FC5AV1C-L; 1.06D0 used for DECK simulation) 
EOF

cat <<EOF >> user_nl_cam
 !.......................................................
 ! nudging
 !.......................................................
  Nudge_Model  = .False.
  Nudge_Method = 'Linear'
  Nudge_Path   = '${INPUT_NUDGING}/PRENUDG_F20TRC5-CMIP6_DT1800/'
  Nudge_File_Template = 'PRENUDG_F20TRC5-CMIP6_DT1800.cam.h1.%y-%m-%d-%s.nc'
  Nudge_Times_Per_Day = 4  !! nudging input data frequency
  Model_Times_Per_Day = 48 !! should not be larger than 48 if dtime = 1800s
  Nudge_Uprof = 1
  Nudge_Ucoef = 1.
  Nudge_Vprof = 1
  Nudge_Vcoef = 1.
  Nudge_Tprof = 0
  Nudge_Tcoef = 0.
  Nudge_Qprof = 0
  Nudge_Qcoef = 0.
  Nudge_PSprof = 0
  Nudge_PScoef = 0.
  Nudge_Beg_Year = 2009
  Nudge_Beg_Month = 12
  Nudge_Beg_Day = 1
  Nudge_End_Year = 2010
  Nudge_End_Month = 12
  Nudge_End_Day = 31
EOF

cat > user_nl_clm <<EOF
 finidat = '$lnd_init'
EOF

   # Run model 
   cd $CASEROOT
   ./case.submit

date
endif


