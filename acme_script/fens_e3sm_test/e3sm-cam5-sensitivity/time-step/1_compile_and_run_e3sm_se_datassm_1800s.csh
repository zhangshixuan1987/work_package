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
set nlen      = 2   # run model for 3 days
set stop_o    = 'ndays'
set ires      = 1   #
set rest_o    = 'nmonths'

set nhtfrq    = 1    # monthly average
set mfilt     = 120   # 1 time step per file 

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
setenv RESOLUTION f19_g16  
setenv MACH compy
setenv COMPILER intel
setenv scheduler 'PBS'

setenv NTASKS_PER_INST  240 #6 #24
setenv NINST  1
setenv NTHRDS 1
setenv NCORES_PER_NODE  24
setenv NDRIVERS $NINST

setenv WORKDIR /compyfs/zhan391

#setenv CSMDATA /compyfs/inputdata/
setenv CSMDATA  /compyfs/zhan391/acme_init/csmdata

#setenv INPUT_NUDGING /compyfs/zhan391/F20TRC5-CMIP6_201806_initgen_84nodes

set taskname = step_datass_test 
set casename = CNTL_${COMPSET}_${RESOLUTION}_DT`printf "%04d" ${dtime}`
set MULTI_DRIVER = ' '
#if ($NDRIVERS > 1) set MULTI_DRIVER = ' --multi-driver '

#-----------------
#setup branch run information
#set dtime1      = 1800 
#set run_refcase = AMIPNUDG_${COMPSET}_DT`printf "%04d" ${dtime1}`
set run_refdate = "2009-01-01"
set run_reftod  = "00000"
set initDir     = /compyfs/zhan391/acme_init/csmdata #${INPUT_NUDGING}/${run_refcase}/
set atm_init    = "${CSMDATA}/atm/cam/inic/fv/cami-mam3_0000-01-01_1.9x2.5_L30_c090306.nc"
set lnd_init    = "${CSMDATA}/lnd/clm2/initdata_map/clmi.I20TRGSWCNPRDCTCBC.f19_g16.67a0f98.0401-01-01-00000.nc"
set lnd_surf    = "${CSMDATA}/lnd/clm2/surfdata_map/surfdata_1.9x2.5_simyr1850_c180306.nc"
set lnd_use     = "${CSMDATA}/lnd/clm2/surfdata_map/landuse.timeseries_1.9x2.5_hist_simyr1850-2015_c180306.nc"

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

set case_setup_script = "create_and_setup_bundled_spin.csh"
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
#  ./xmlchange -file env_build.xml -id GMAKE_J -val '16'
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

   ./xmlchange  -file env_run.xml -id SSTICE_DATA_FILENAME -val "$CSMDATA/ocn/docn7/SSTDATA/sst_ice_CMIP6_DECK_E3SM_1x1_c20180213.nc"
   ./xmlchange  -file env_run.xml -id SSTICE_YEAR_ALIGN -val "1849"
   ./xmlchange  -file env_run.xml -id SSTICE_YEAR_START -val "1849"
   ./xmlchange  -file env_run.xml -id SSTICE_YEAR_END   -val "2016"

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
   ./xmlchange  -file env_run.xml -id  DOUT_L_MS    -val 'FALSE'
   ./xmlchange  -file env_run.xml -id  DOUT_S       -val 'FALSE'
   ./xmlchange  -file env_run.xml -id  DOUT_S_ROOT  -val "$WORKDIR/e3sm_scratch/archive/$CASE"

   ./xmlchange  -file env_batch.xml -id JOB_WALLCLOCK_TIME -val "01:00:00"
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

# Get a bunch of environment variables.
# If any of these are changed by xmlchange calls in this program,
# then they must be explicitly changed with setenv calls too.
# $COMPSET is the long name which E3SM uses, and is derived from $compset.
# $compset is set by the user and may be an alias/short name.
setenv COMPSET            `./xmlquery COMPSET            --value`
setenv COMP_OCN           `./xmlquery COMP_OCN           --value`
setenv COMP_GLC           `./xmlquery COMP_GLC           --value`
setenv COMP_ROF           `./xmlquery COMP_ROF           --value`
setenv CIMEROOT           `./xmlquery CIMEROOT           --value`
setenv EXEROOT            `./xmlquery EXEROOT            --value`
setenv RUNDIR             `./xmlquery RUNDIR             --value`
setenv CAM_CONFIG_OPTS    `./xmlquery CAM_CONFIG_OPTS    --value`
set    max_tasks_per_node = `./xmlquery    MAX_TASKS_PER_NODE --value`
set max_mpitasks_per_node = `./xmlquery MAX_MPITASKS_PER_NODE --value`
echo "From create_newcase, settings related to TASKS = ..."
./xmlquery --partial TASK
# Make sure the case is configured with a data ocean.
if ( ${COMP_OCN} != docn ) then
   echo " "
   echo "ERROR: This setup script is not appropriate for active ocean compsets."
   echo "ERROR: Please use the models/E3SM/shell_scripts examples for that case."
   echo " "
   exit 40
endif

# Extract pieces of the COMPSET for choosing correct setup parameters.
# E.g. "AMIP_CAM5_CLM50%BGC_CICE%PRES_DOCN%DOM_MOSART_CISM1%NOEVOLVE_SWAV"
set comp_list = `echo $COMPSET   | sed -e "s/_/ /g"`

# Debug
echo "compset parts are $comp_list"
# Land ice, aka glacier, aka glc.
if (${COMP_GLC} == sglc) then
   set CISM_RESTART = FALSE
else
   echo "ERROR: glacier compset is ${COMP_GLC}, which is not supported by CURRENT DART DA script."
   echo "ERROR: The only supported glacier compset is 'SGLC'"
   exit 45
   # In the future, if CISM can use the GREGORIAN calandar, and evolving land ice is
   # deemed to be useful for atmospheric assimilations, this may still be required
   # to make CISM write out restart files 4x/day.
   ./xmlchange GLC_NCPL=4
endif
# The river transport model ON is useful only when using an active ocean or
# land surface diagnostics. If you turn it ON, you will have to stage initial files etc.
# There are 3 choices:
# > a stub version (best for CAM+DART),
# > the older River Transport Model (RTM),
# > the new Model for Scale Adaptive River Transport (MOSART).
# They are separate E3SM components, and are/need to be specified in the compset.
# It may be that RTM or MOSART can be turned off via namelists.
# Specify the river runoff model: 'RTM', 'MOSART', or anything else.
if (${COMP_ROF} == 'rtm') then
   ./xmlchange ROF_GRID='r05'
else if (${COMP_ROF} == 'mosart') then
   # There seems to be no MOSART_MODE, but there are some MOSART_ xml variables.
   # Use defaults for now
   ./xmlchange ROF_GRID='r05'
else if (${COMP_ROF} == 'drof') then
   ./xmlchange ROF_GRID='null'
else if (${COMP_ROF} == 'srof') then
   ./xmlchange ROF_GRID='null'
else
   echo "river_runoff is ${COMP_ROF}, which is not supported"
   exit 50
endif

./xmlchange CALENDAR=GREGORIAN
./xmlchange CONTINUE_RUN=FALSE


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
 ncdata              = '$atm_init'
 inithist            = 'MONTHLY'
 inithist_all        = .true.
 nhtfrq              = $nhtfrq,
 mfilt               = $mfilt,
 ndens               = 1,
 fincl1               = 'FLDS:A','FLNS:A','FLNT:A','FLNTC:A','FLUT:A','FLUTC:A','LWCF:A','FSNS:A','FSNSC:A','FSNTOA:A','FSNTOAC:A','FSNT:A','FSNTC:A','SWCF:A','CLDLOW:I','CLDHGH:I','CLDMED:I','CLDTOT:I','PBLH:I','PHIS:I','OMEGA850:I','OMEGA500:I','TGCLDLWP:I','TGCLDIWP:I','TGCLDCWP:I','TMQ:I','TH8501000:I','TH7001000:I','THE7001000:I','THE8501000:I','AODVIS:I','LANDFRAC:A','OCNFRAC:A','PSL:I','PS:I','LHFLX:A','SHFLX:A','QFLX:A','TAUX:I','TAUY:I','TREFHT:I','TS:I','QREFHT:I','RHREFHT:I','U10:I','PRECC:A','PRECL:A','T1000:I','T925:I','T850:I','T700:I','T500:I','T300:I','T200:I','Q1000:I','Q925:I','Q850:I','Q200:I','U850:I','U250:I','U200:I','V850:I','V250:I','V200:I','Z1000:I','Z700:I','Z500:I','Z300:I','Z200:I','Z100:I','RELHUM:I', 'T:I','U:I','V:I','Q:I','OMEGA:I','Z3:I',
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
 finidat = '$lnd_init'
 fsurdat = '$lnd_surf'
 flanduse_timeseries = '$lnd_use'
 check_dynpft_consistency = .false.
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

