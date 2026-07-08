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

./xmlchange --id RUNDIR  --val $RUNDIR
./xmlchange --id EXEROOT --val $EXEDIR

./xmlchange --id DIN_LOC_ROOT          --val $CSMDATA
./xmlchange --id DIN_LOC_ROOT_CLMFORC  --val $CSMDATA

# Specify PE layout and # of instances

@ nproc = $NTASKS_PER_INST * $NINST

./xmlchange  --id NTASKS_ATM --val $nproc
./xmlchange  --id NTHRDS_ATM --val $NTHRDS
./xmlchange  --id ROOTPE_ATM --val '0'
./xmlchange  --id NINST_ATM --val $NINST
./xmlchange  --id NINST_ATM_LAYOUT --val 'concurrent'

./xmlchange  --id NTASKS_LND --val $nproc
./xmlchange  --id NTHRDS_LND --val $NTHRDS
./xmlchange  --id ROOTPE_LND --val '0'
./xmlchange  --id NINST_LND  --val $NINST
./xmlchange  --id NINST_LND_LAYOUT --val 'concurrent'

./xmlchange  --id NTASKS_ROF --val $nproc
./xmlchange  --id NTHRDS_ROF --val $NTHRDS
./xmlchange  --id ROOTPE_ROF --val '0'
./xmlchange  --id NINST_ROF  --val $NINST
./xmlchange  --id NINST_ROF_LAYOUT --val 'concurrent'

./xmlchange  --id NTASKS_ICE --val $nproc
./xmlchange  --id NTHRDS_ICE --val $NTHRDS 
./xmlchange  --id ROOTPE_ICE --val '0'
./xmlchange  --id NINST_ICE --val $NINST
./xmlchange  --id NINST_ICE_LAYOUT --val 'concurrent'

./xmlchange  --id NTASKS_OCN --val $nproc
./xmlchange  --id NTHRDS_OCN --val $NTHRDS 
./xmlchange  --id ROOTPE_OCN --val '0'
./xmlchange  --id NINST_OCN --val $NINST
./xmlchange  --id NINST_OCN_LAYOUT --val 'concurrent'

# GLC and WAV are stub components
./xmlchange  --id NTASKS_GLC --val $nproc
./xmlchange  --id NTHRDS_GLC --val $NTHRDS 
./xmlchange  --id ROOTPE_GLC --val '0'
./xmlchange  --id NINST_GLC --val '1'
./xmlchange  --id NINST_GLC_LAYOUT --val 'concurrent'

./xmlchange  --id NTASKS_WAV --val $nproc
./xmlchange  --id NTHRDS_WAV --val $NTHRDS 
./xmlchange  --id ROOTPE_WAV --val '0'
./xmlchange  --id NINST_WAV --val '1'
./xmlchange  --id NINST_WAV_LAYOUT --val 'concurrent'

./xmlchange  --id NTASKS_CPL --val $nproc
./xmlchange  --id NTHRDS_CPL --val $NTHRDS 
./xmlchange  --id ROOTPE_CPL --val '0'

./case.setup

