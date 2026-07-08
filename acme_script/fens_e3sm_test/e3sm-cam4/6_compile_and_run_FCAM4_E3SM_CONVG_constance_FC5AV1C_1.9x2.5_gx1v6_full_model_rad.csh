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

set dtime = 1800
set nlen      = 6   # run model for 3 days
set stop_o    = 'nhours'
set ires      = 6   #
set rest_o    = 'nhours'

set nhtfrq    = -6  # monthly average
set mfilt     = 1   # 1 time step per file 

####################################################################
# Fetch code
####################################################################
setenv CCSMTAG E3SM_MAINT1.0 #E3SM_201910 #ACME_master_20180611
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

setenv postCIME 5
#<lname>20TR_CAM5_CLM45%SPBC_CICE%PRES_DOCN%DOM_SROF_SGLC_SWAV</lname>
setenv COMPSET F20TRC5-AMIP #FAMIPC5 #F20TRC5-CMIP6 
setenv RESOLUTION f19_g16 #f45_g37 #ne30_ne30 #ne30_oECv3
setenv MACH compy
setenv COMPILER intel
setenv scheduler 'PBS'

setenv NTASKS_PER_INST  24
setenv NINST  10
setenv NTHRDS 4
setenv NCORES_PER_NODE  24

setenv WORKDIR /compyfs/zhan391

#setenv CSMDATA /compyfs/inputdata/
setenv CSMDATA  /compyfs/zhan391/acme_init/csmdata

#setenv INPUT_NUDGING /compyfs/zhan391/F20TRC5-CMIP6_201806_initgen_84nodes

set taskname = fens_e3sm_cam5 
set casename = DENS_${COMPSET}_${RESOLUTION}_DT`printf "%04d" ${dtime}`

#-----------------
#setup branch run information
#set dtime1      = 1800 
#set run_refcase = AMIPNUDG_${COMPSET}_DT`printf "%04d" ${dtime1}`
set run_refdate = "2009-01-01"
set run_reftod  = "00000"
set initDir     = /compyfs/zhan391/acme_init #${INPUT_NUDGING}/${run_refcase}/

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

set case_setup_script = "create_and_setup_bundled_fens.csh"
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
#   ./xmlchange -file env_build.xml -id GMAKE_J -val '16'
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
#   if ( ! -d $RUNDIR ) then
#      mkdir -p $RUNDIR
#   endif
#   set camiin  = "$run_refcase.cam.i.$run_refdate-$run_reftod.nc"
#   set camrin  = "$run_refcase.cam.r.$run_refdate-$run_reftod.nc"
#   set camrsin = "$run_refcase.cam.rs.$run_refdate-$run_reftod.nc"
#  #set camrhin = "$run_refcase.cam.rh0.$run_refdate-$run_reftod.nc"
#   set lndiin  = "$run_refcase.clm2.r.$run_refdate-$run_reftod.nc"
#  #set lndrin  = "$run_refcase.clm2.rh0.$run_refdate-$run_reftod.nc"
#   set cicerin = "$run_refcase.cice.r.$run_refdate-$run_reftod.nc"
#   set cplrin  = "$run_refcase.cpl.r.$run_refdate-$run_reftod.nc"
#   set docnrin = "$run_refcase.docn.rs1.$run_refdate-$run_reftod.bin"
#   cp -r ${initDir}${camiin}  $RUNDIR/${camiin}
#   cp -r ${initDir}${camrin}  $RUNDIR/${camrin}
#   cp -r ${initDir}${camrsin} $RUNDIR/${camrsin}
#  #cp -r ${initDir}${camrhin} $RUNDIR/${camrhin}
#   cp -r ${initDir}${lndiin}  $RUNDIR/${lndiin}
#  #cp -r ${initDir}${lndrin}  $RUNDIR/${lndrin}
#   cp -r ${initDir}${cicerin} $RUNDIR/${cicerin}
#   cp -r ${initDir}${cplrin}  $RUNDIR/${cplrin}
#   cp -r ${initDir}${docnrin} $RUNDIR/${docnrin}
#   cp -r ${initDir}rpoint*    $RUNDIR/
#   foreach  fil ( $RUNDIR/rpoint* )
#      sed -i "s/2011-03-01-00000/${run_refdate}-${run_reftod}/g" $fil
#   end

#  set atm_init   = ${initDir}/"$run_refcase.cam.i.$run_refdate-$run_reftod.nc"  #"$RUNDIR/${camiin}"
#  set lnd_init   = ${initDir}/"$run_refcase.clm2.r.$run_refdate-$run_reftod.nc" #"$RUNDIR/${lndiin}"

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

#   ./xmlchange  -file env_run.xml -id SSTICE_DATA_FILENAME -val "$initDir/sst_ice_CMIP6_DECK_E3SM_1x1_c20180213.nc"
#   ./xmlchange  -file env_run.xml -id SSTICE_YEAR_ALIGN -val "1849"
#   ./xmlchange  -file env_run.xml -id SSTICE_YEAR_START -val "1849"
#   ./xmlchange  -file env_run.xml -id SSTICE_YEAR_END   -val "2016"

  # ./xmlchange  -file env_run.xml -id  RUN_STARTDATE  -val $start_date
   ./xmlchange  -file env_run.xml -id  REST_N         -val $ires
   ./xmlchange  -file env_run.xml -id  REST_OPTION    -val $rest_o

  # Set sim length so that a cycle can be finished in 8 hours
  #
   @ nmon = $dtime * 16 / 2000
   ./xmlchange  -file env_run.xml -id  STOP_N       -val $nlen
   ./xmlchange  -file env_run.xml -id  STOP_OPTION  -val $stop_o

   @ nresub = 68 / $nmon
   ./xmlchange  -file env_run.xml -id  RESUBMIT     -val '0'
   ./xmlchange  -file env_run.xml -id  SAVE_TIMING_DIR -val $PTMP
   ./xmlchange  -file env_run.xml -id  DOUT_S       -val 'FALSE'
   ./xmlchange  -file env_run.xml -id  DOUT_L_MS    -val 'FALSE'

   ./xmlchange  -file env_batch.xml -id JOB_WALLCLOCK_TIME -val "1:00:00"
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


set ininst = 1
while ( $ininst <= $NINST )

if ($NINST <= 1) then
set user_cam = user_nl_cam
set user_clm = user_nl_clm
else
set user_cam = user_nl_cam_`printf "%04d" ${ininst}`
set user_clm = user_nl_clm_`printf "%04d" ${ininst}`
endif
 
# Namelist variables
# fincl1             = 'PS','U','V','T','Q','CLDLIQ','CLDICE','NUMLIQ','NUMICE','num_a1','num_a2','num_a3','LANDFRAC'
# ncdata             = '$atm_init'
cat > $user_cam <<EOF
 inithist            = '6-HOURLY'
 inithist_all        = .true.
 avgflag_pertape     = 'A','I'
 nhtfrq              = $nhtfrq,-6
 mfilt               = $mfilt,1
 ndens               = 1,1
 fincl2              = 'PHIS'
 ext_frc_specifier              = 'SO2         -> /compyfs/zhan391/acme_init/csmdata/atm/cam/chem/trop_mozart_aero/emis/ar5_mam3_so2_elev_1850-2010_c20100902_v12.nc',
         'bc_a1       -> /compyfs/zhan391/acme_init/csmdata/atm/cam/chem/trop_mozart_aero/emis/ar5_mam3_bc_elev_1850-2010_c20100902_v12.nc',
         'num_a1      -> /compyfs/zhan391/acme_init/csmdata/atm/cam/chem/trop_mozart_aero/emis/ar5_mam3_num_a1_elev_1850-2010_c20100902_v12.nc',
         'num_a2      -> /compyfs/zhan391/acme_init/csmdata/atm/cam/chem/trop_mozart_aero/emis/ar5_mam3_num_a2_elev_1850-2010_c20100902_v12.nc',
         'pom_a1      -> /compyfs/zhan391/acme_init/csmdata/atm/cam/chem/trop_mozart_aero/emis/ar5_mam3_pom_elev_1850-2010_c20140421_v12.nc',
         'so4_a1      -> /compyfs/zhan391/acme_init/csmdata/atm/cam/chem/trop_mozart_aero/emis/ar5_mam3_so4_a1_elev_1850-2010_c20100902_v12.nc',
         'so4_a2      -> /compyfs/zhan391/acme_init/csmdata/atm/cam/chem/trop_mozart_aero/emis/ar5_mam3_so4_a2_elev_1850-2010_c20100902_v12.nc'
 srf_emis_specifier             = 'DMS       -> /compyfs/zhan391/acme_init/csmdata/atm/cam/chem/trop_mozart_aero/emis/aerocom_mam3_dms_surf_1849-2010_c20100902.nc',
         'SO2       -> /compyfs/zhan391/acme_init/csmdata/atm/cam/chem/trop_mozart_aero/emis/ar5_mam3_so2_surf_1850-2010_c20100902_v12.nc',
         'SOAG      -> /compyfs/zhan391/acme_init/csmdata/atm/cam/chem/trop_mozart_aero/emis/ar5_mam3_soag_1.5_surf_1850-2010_c20100902_v12.nc',
         'bc_a1     -> /compyfs/zhan391/acme_init/csmdata/atm/cam/chem/trop_mozart_aero/emis/ar5_mam3_bc_surf_1850-2010_c20100902_v12.nc',
         'num_a1    -> /compyfs/zhan391/acme_init/csmdata/atm/cam/chem/trop_mozart_aero/emis/ar5_mam3_num_a1_surf_1850-2010_c20100902_v12.nc',
         'num_a2    -> /compyfs/zhan391/acme_init/csmdata/atm/cam/chem/trop_mozart_aero/emis/ar5_mam3_num_a2_surf_1850-2010_c20100902_v12.nc',
         'pom_a1    -> /compyfs/zhan391/acme_init/csmdata/atm/cam/chem/trop_mozart_aero/emis/ar5_mam3_pom_surf_1850-2010_c20140421_v12.nc',
         'so4_a1    -> /compyfs/zhan391/acme_init/csmdata/atm/cam/chem/trop_mozart_aero/emis/ar5_mam3_so4_a1_surf_1850-2010_c20100902_v12.nc',
         'so4_a2    -> /compyfs/zhan391/acme_init/csmdata/atm/cam/chem/trop_mozart_aero/emis/ar5_mam3_so4_a2_surf_1850-2010_c20100902_v12.nc'
 bndtvghg               = '/compyfs/zhan391/acme_init/csmdata/atm/cam/ggas/ghg_hist_1765-2012_c130501.nc'
 prescribed_ozone_file          = 'ozone_1.9x2.5_L26_1990-2029_c080324.nc'
 prescribed_volcaero_file               = 'CCSM4_volcanic_1850-2011_prototype1.nc'
 solar_data_file                = '/compyfs/zhan391/acme_init/csmdata/atm/cam/solar/Solar_1850-2299_input4MIPS_c20171207.nc'
 tracer_cnst_file               = 'oxid_1.9x2.5_L26_1850-2015_c20171110.nc'
 tracer_cnst_filelist           = ''
EOF

#finidat = '$lnd_init'
#finidat = '/compyfs/zhan391/acme_init/csmdata/lnd/clm2/initdata/clmi.BCN.2000-01-01_1.9x2.5_gx1v6_simyr2000_c100309.nc'
cat > $user_clm <<EOF
 hist_empty_htapes = .true.
 hist_fincl1 = 'TSA'
 hist_nhtfrq = -6
 hist_mfilt  = 1
 hist_avgflag_pertape = 'I'
EOF

@ ininst++

end

# Run model 
cd $CASEROOT
./case.submit

date
endif

