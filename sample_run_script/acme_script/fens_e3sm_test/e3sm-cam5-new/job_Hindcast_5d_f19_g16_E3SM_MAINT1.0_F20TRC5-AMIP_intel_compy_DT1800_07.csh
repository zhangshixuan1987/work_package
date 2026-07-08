#!/bin/csh
date
setenv PROJECT         'ESMD'
setenv CESM_PROJ       'ESMD'
setenv CESM_EMAIL      'shixuan.zhang@pnnl.gov'

setenv MACH            'compy'
setenv COMPILER        'intel'
setenv RESOLUTION      'f19_g16'
setenv NTASKS_PER_INST '240'
setenv NINST           '1'
setenv NTHRDS         '1'
         
setenv CCSMTAG  'E3SM_MAINT1.0'
setenv COMPSET  'F20TRC5-AMIP'
setenv CCSMROOT '/compyfs/zhan391/code/E3SM_MAINT1.0'
setenv PTMP     '/compyfs/zhan391/cam5_datass_test'
setenv EXEDIR   '/compyfs/zhan391/cam5_datass_test/exe//compile_E3SM_MAINT1.0_F20TRC5-AMIP_f19_g16_compy_intel_240x1threads'
setenv CSMDATA  '/compyfs/zhan391/acme_init/csmdata'
         
setenv CASE      'Hindcast_5d_f19_g16_E3SM_MAINT1.0_F20TRC5-AMIP_intel_compy_DT1800_07'
setenv CASEROOT  '/compyfs/zhan391/cam5_datass_test/cases/Hindcast_5d_f19_g16_E3SM_MAINT1.0_F20TRC5-AMIP_intel_compy_DT1800_07'
setenv RUNDIR    '/compyfs/zhan391/cam5_datass_test/run/Hindcast_5d_f19_g16_E3SM_MAINT1.0_F20TRC5-AMIP_intel_compy_DT1800_07'

setenv postCIME 5

set nbundles = -1         
set case_setup_script = '/compyfs/zhan391/run_script/fens_e3sm_test/e3sm-cam5-new/create_and_setup_bundled_hind.csh' 

set run_refcase = 'DACAM5_ENS80_F20TRC5-AMIP_f19_g16_DT1800'
set run_refdate = '2009-02-05'
set run_reftod  = '00000'

set dtime    = 1800
set atm_init = /compyfs/zhan391/cam5_datass_test/run/Hindcast_5d_f19_g16_E3SM_MAINT1.0_F20TRC5-AMIP_intel_compy_DT1800_07/DACAM5_ENS80_F20TRC5-AMIP_f19_g16_DT1800.cam.i.2009-02-05-00000.nc
set lnd_init = /compyfs/zhan391/cam5_datass_test/run/Hindcast_5d_f19_g16_E3SM_MAINT1.0_F20TRC5-AMIP_intel_compy_DT1800_07/DACAM5_ENS80_F20TRC5-AMIP_f19_g16_DT1800.clm2.r.2009-02-05-00000.nc
set ice_init = /compyfs/zhan391/cam5_datass_test/run/Hindcast_5d_f19_g16_E3SM_MAINT1.0_F20TRC5-AMIP_intel_compy_DT1800_07/DACAM5_ENS80_F20TRC5-AMIP_f19_g16_DT1800.cice.r.2009-02-05-00000.nc
set lnd_surf = /compyfs/zhan391/acme_init/csmdata/lnd/clm2/surfdata_map/surfdata_1.9x2.5_simyr1850_c180306.nc
set lnd_use  = /compyfs/zhan391/acme_init/csmdata/lnd/clm2/surfdata_map/landuse.timeseries_1.9x2.5_hist_simyr1850-2015_c180306.nc
set ncycle   = 1      # no resubmit

set nlen     = 5      # run model for 5 days
set stop_o   = 'ndays'

#====================================================================
# Create and set up new case. No need to build the model. 
#====================================================================

source ${case_setup_script}

cd $CASEROOT
./xmlchange -file env_build.xml -id BUILD_COMPLETE  -val 'TRUE'

   ./xmlchange  -file env_run.xml -id SSTICE_DATA_FILENAME -val "$CSMDATA/ocn/docn7/SSTDATA/sst_ice_CMIP6_DECK_E3SM_1x1_c20180213.nc"
   ./xmlchange  -file env_run.xml -id SSTICE_YEAR_ALIGN -val "1849"
   ./xmlchange  -file env_run.xml -id SSTICE_YEAR_START -val "1849"
   ./xmlchange  -file env_run.xml -id SSTICE_YEAR_END   -val "2016"

  # ./xmlchange  -file env_run.xml -id  RUN_STARTDATE  -val $start_date
   ./xmlchange  -file env_run.xml -id  REST_N         -val '1'
   ./xmlchange  -file env_run.xml -id  REST_OPTION    -val 'nmonths'

   ./xmlchange  -file env_run.xml -id  STOP_N       -val $nlen
   ./xmlchange  -file env_run.xml -id  STOP_OPTION  -val $stop_o

   ./xmlchange  -file env_run.xml -id  RESUBMIT     -val '0'
   ./xmlchange  -file env_run.xml -id  SAVE_TIMING_DIR -val $PTMP
   ./xmlchange  -file env_run.xml -id  DOUT_L_MS    -val 'FALSE'
   ./xmlchange  -file env_run.xml -id  DOUT_S       -val 'FALSE'
   ./xmlchange  -file env_run.xml -id  DOUT_S_ROOT  -val "$WORKDIR/e3sm_scratch/archive/$CASE"

   ./xmlchange  -file env_batch.xml -id JOB_WALLCLOCK_TIME -val "01:00:00"

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


if ( ${COMP_OCN} != docn ) then
   echo " "
   echo "ERROR: This setup script is not appropriate for active ocean compsets."
   echo "ERROR: Please use the models/E3SM/shell_scripts examples for that case."
   echo " "
   exit 40
endif

set comp_list = `echo $COMPSET   | sed -e "s/_/ /g"`

echo "compset parts are $comp_list"
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

if (${COMP_ROF} == 'rtm') then
   ./xmlchange ROF_GRID='r05'
else if (${COMP_ROF} == 'mosart') then
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

cat > user_nl_cam <<EOF
 ncdata              = '$atm_init'
 inithist            = 'MONTHLY'
 inithist_all        = .true.
 history_amwg         = .true.
 history_verbose      = .false.
 history_aero_optics  = .true.
 history_aerosol      = .true.
 history_clubb        = .false.
 history_budget       = .true.
!history_microphysics = .true.
 empty_htapes         = .false.
 do_aerocom_ind3      = .false.
 nhtfrq               = -24,-6,-3,-1
 mfilt                =  1, 4, 8,24
 avgflag_pertape      = 'A','I','I','I'
 ndens                = 1,1,1,1
 fincl1               = 'OMEGA500','RH700','TH7001000','THE7001000',
 fincl2               = 'PS','U','V','T','Q','Z3','OMEGA',
                        'DPDLFLIQ','DPDLFICE','DPDLFT'
                        'PSL','T200','T500','T850','U200','U500','U850','V200','V500',
                        'V850','Q200','Q500','Q850','RH200','RH500','RH850','OMEGA850',
                        'OMEGA500','OMEGA200','Z200','Z500','Z850','UBOT','VBOT','ZBOT',
                         'VBOT','TREFHT','T7001000','TH7001000','THE7001000',
 fincl3               = 'AODVIS','angstrm','cod','cdr','cdnc',
                        'cdnum','icnum','clt','lcc','lwp',
                        'iwp','icc','icnc','icr','LHFLX',
                        'SHFLX','OMEGA500','rh700','colrv',
                        'ccn','ccn.1bl','ccn.3bl','ptop',
                        'ttop','rwp','lwp2','iwp2','autoconv',
                        'accretn','FSUTOA','FSUTOAC','FLUT','FLUTC',
                        'PRECL','PRECC',
 fincl4               = 'QFLX','PRECC','PRECL','LHFLX','SHFLX','FLNT',
                        'FSNT','FLNS','FSNS','TREFHT','CLDHGH','CLDLOW',
                        'CLDMED','CLDTOT','TGCLDLWP','TGCLDIWP',
                        'RELHUM','RHREFHT','QREFHT','PBLH','U10','PSL',
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

cat > user_nl_clm <<EOF
 finidat = '$lnd_init'
 fsurdat = '$lnd_surf'
 flanduse_timeseries = '$lnd_use'
 check_dynpft_consistency = .false.
!hist_empty_htapes = .true.
!hist_fincl1 = 'TSA'
 hist_nhtfrq = -24
 hist_mfilt  = 1
 hist_avgflag_pertape = 'A'
EOF

cat > user_nl_cice <<EOF
ice_ic = '$ice_init'
EOF

# Run model 
cd $CASEROOT
./case.submit

