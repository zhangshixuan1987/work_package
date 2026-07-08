#====================================================================
# create new case
#====================================================================
rm -rf $CASEROOT

if ( $postCIME == 0 ) then
    cd  $CCSMROOT/scripts
else
    cd  $CCSMROOT/cime/scripts
endif

./create_newcase -case $CASEROOT -mach $MACH \
                 -res $RESOLUTION -compset $COMPSET -compiler $COMPILER -v

#====================================================================
# set up case
#====================================================================
cd $CASEROOT

# Although a bit counter-intuitive, both EXEROOT and RUNDIR have to be 
# specified before cesm_setup if we don't want to use the default

./xmlchange -file env_run.xml   -id RUNDIR  -val $RUNDIR
./xmlchange -file env_build.xml -id EXEROOT -val $EXEDIR

#./xmlchange -file env_run.xml   -id RUN_TYPE      -val   branch
#./xmlchange -file env_run.xml   -id RUN_REFCASE   -val   $run_refcase
#./xmlchange -file env_run.xml   -id RUN_REFDATE   -val   $run_refdate
#./xmlchange -file env_run.xml   -id RUN_REFTOD    -val   $run_reftod
./xmlchange -file env_run.xml   -id RUN_STARTDATE -val   $run_refdate
./xmlchange -file env_run.xml   -id START_TOD     -val   $run_reftod
#./xmlchange -file env_run.xml   -id DATA_ASSIMILATION        -val  'TRUE'
#./xmlchange -file env_run.xml   -id DATA_ASSIMILATION_CYCLES -val  1
#./xmlchange -file env_run.xml   -id DATA_ASSIMILATION_SCRIPT -val '${RUNDIR}/no_assimilate.csh'

./xmlchange  -file env_run.xml -id  PIO_TYPENAME_ATM  -val netcdf
./xmlchange  -file env_run.xml -id  PIO_TYPENAME_CPL  -val netcdf
./xmlchange  -file env_run.xml -id  PIO_TYPENAME_OCN  -val netcdf
./xmlchange  -file env_run.xml -id  PIO_TYPENAME_WAV  -val netcdf
./xmlchange  -file env_run.xml -id  PIO_TYPENAME_GLC  -val netcdf
./xmlchange  -file env_run.xml -id  PIO_TYPENAME_ICE  -val netcdf
./xmlchange  -file env_run.xml -id  PIO_TYPENAME_ROF  -val netcdf
./xmlchange  -file env_run.xml -id  PIO_TYPENAME_LND  -val netcdf
./xmlchange  -file env_run.xml -id  PIO_TYPENAME_ESP  -val netcdf

./xmlchange -file env_run.xml   -id DIN_LOC_ROOT          -val $CSMDATA
./xmlchange -file env_run.xml   -id DIN_LOC_ROOT_CLMFORC  -val $CSMDATA

# Specify PE layout and # of instances

@ nproc = $NTASKS_PER_INST * $NINST

./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val $nproc
./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val $NTHRDS
./xmlchange -file env_mach_pes.xml -id ROOTPE_ATM -val '0'
./xmlchange -file env_mach_pes.xml -id NINST_ATM -val $NINST
./xmlchange -file env_mach_pes.xml -id NINST_ATM_LAYOUT -val 'concurrent'

./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val $nproc
./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val $NTHRDS
./xmlchange -file env_mach_pes.xml -id ROOTPE_LND -val '0'
./xmlchange -file env_mach_pes.xml -id NINST_LND  -val $NINST
./xmlchange -file env_mach_pes.xml -id NINST_LND_LAYOUT -val 'concurrent'

./xmlchange -file env_mach_pes.xml -id NTASKS_ROF -val $nproc
./xmlchange -file env_mach_pes.xml -id NTHRDS_ROF -val $NTHRDS
./xmlchange -file env_mach_pes.xml -id ROOTPE_ROF -val '0'
./xmlchange -file env_mach_pes.xml -id NINST_ROF  -val $NINST
./xmlchange -file env_mach_pes.xml -id NINST_ROF_LAYOUT -val 'concurrent'

./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val $nproc
./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val $NTHRDS 
./xmlchange -file env_mach_pes.xml -id ROOTPE_ICE -val '0'
./xmlchange -file env_mach_pes.xml -id NINST_ICE -val $NINST
./xmlchange -file env_mach_pes.xml -id NINST_ICE_LAYOUT -val 'concurrent'

./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val $nproc
./xmlchange -file env_mach_pes.xml -id NTHRDS_OCN -val $NTHRDS 
./xmlchange -file env_mach_pes.xml -id ROOTPE_OCN -val '0'
./xmlchange -file env_mach_pes.xml -id NINST_OCN -val $NINST
./xmlchange -file env_mach_pes.xml -id NINST_OCN_LAYOUT -val 'concurrent'

# GLC and WAV are stub components
./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val $nproc
./xmlchange -file env_mach_pes.xml -id NTHRDS_GLC -val $NTHRDS 
./xmlchange -file env_mach_pes.xml -id ROOTPE_GLC -val '0'
./xmlchange -file env_mach_pes.xml -id NINST_GLC -val '1'
./xmlchange -file env_mach_pes.xml -id NINST_GLC_LAYOUT -val 'concurrent'

./xmlchange -file env_mach_pes.xml -id NTASKS_WAV -val $nproc
./xmlchange -file env_mach_pes.xml -id NTHRDS_WAV -val $NTHRDS 
./xmlchange -file env_mach_pes.xml -id ROOTPE_WAV -val '0'
./xmlchange -file env_mach_pes.xml -id NINST_WAV -val '1'
./xmlchange -file env_mach_pes.xml -id NINST_WAV_LAYOUT -val 'concurrent'

./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val $nproc
./xmlchange -file env_mach_pes.xml -id NTHRDS_CPL -val $NTHRDS 
./xmlchange -file env_mach_pes.xml -id ROOTPE_CPL -val '0'

./case.setup

