#!/bin/csh
date

##############################################################################################
# This is the driver scripts for the Time Step Convergence (TSC) test simulations which
# - compiles the model (just once);
# - generates one script for each simulation; 
# - runs those scripts to create cases and get ready for production run;
# - if bundling simulations to form bigger PBS jobs, create PBS job scripts
#   for the bundled jobs then submit; otherwise each simulation will be 
#   submitted as a separate job. 
##############################################################################################
# This script contains some documentation but frankly not much. For questions and comments,
# please contact Hui Wan (hui.wan@pnnl.gov). 
##############################################################################################
set compile_model = 0
set run_model     = 1


set debug = 'TRUE'
set debug = 'FALSE'

#-----------------------------------------------------------------------------------------
# $taskname is used to organize the exe/run/case directories. Not used for other purposes.
set taskname = dart_e3sm_cam4

# $testconfig is used to pick the appropriate template file that contains certain
# run-time options. fullmodel_rad means with all the default parameterizations,
# and with radiation called every other time step.
set testconfig   = 'fullmodel_dart'

#-----------------------
set exp_setup_dir = `pwd`

#-----------------------------
setenv RESOLUTION  1.9x2.5_gx1v6
setenv MACH compy

#-----------------
#setup continue run information
setenv INPUT_DIR  "/compyfs/zhan391/dart_e3sm_cam4/run"
set dtime1       = 1800
set run_refcase  = dart_e3sm_cntl_1.9x2.5_gx1v6_DT`printf "%04d" ${dtime1}`
set run_refdate  = "2009-01-01"
set run_reftod   = "21600"
set initDir      = ${INPUT_DIR}/${run_refcase}/

###################################################################
# This script does NOT clone the model code from the ACME repo.
# Instead, it assumes that the code is located at
# $HOME/codes/${CCSMTAG}/
# The code I used was the branch huiwanpnnl/atm/shcu-pbl-namelist-aero 
# of the ACME code repo.

setenv CCSMTAG E3SM_CAM4_SGR #E3SM_201910
setenv WORKDIR /compyfs/zhan391
setenv CCSMROOT /compyfs/zhan391/code/${CCSMTAG}

# 
setenv postCIME 5

###################################################################
# Conduct simulations with different models or time step sizes.
# Uncomment one of the three blocks at a time
#------------------------------------------------------------------
# Reference solution using default model and 1-second time step. 
#------------------------------------------------------------------
set groupList = ( "dart_e3sm_assm") 
set dtimeList = (1800)  

#------------------------------------------------------------------
# Trusted solution using default model and 2-second time step.
#------------------------------------------------------------------
#set groupList = ("new_init")
#set dtimeList = (2)
#------------------------------------------------------------------
# Test run with 2-second time step. These are parameter perturbation 
# cases taken from Baker et al. (2015, GMD) to evaluate the TSC method.
#------------------------------------------------------------------
#set groupList = ("CONV-LND" "CONV-OCN" "NU-P" "NU")
#set dtimeList = (2)
######################################################

set ndtime = $#dtimeList

set ngroups = $#groupList
set igS = 1
set igE = $ngroups

#====================================================================
# Paths to source code and model input/output
#====================================================================
set scheduler = "PBS"   # default

setenv CESM_EMAIL shixuan.zhang@pnnl.gov
setenv COMPSET F
setenv NINST  1


if ($MACH == "titan") then

   setenv CESM_PROJ cli112
   setenv PROJECT $CESM_PROJ
   setenv COMPILER   pgi
   setenv COMPILER   intel
   setenv NTHRDS 1
   setenv NTASKS_PER_INST 512
   setenv NCORES_PER_NODE  16
   set scheduler = "PBS"

   setenv CSMDATA  /lustre/atlas1/cli900/world-shared/cesm/inputdata
   setenv CCSMROOT $HOME/codes/${CCSMTAG}
   setenv PTMP     /lustre/atlas/proj-shared/cli112/huiwan/$taskname
   setenv EXELOC   $PTMP/exe/$taskname

else if ($MACH == "eos") then

   setenv CESM_PROJ cli112
   setenv PROJECT $CESM_PROJ
   setenv COMPILER   pgi
   setenv COMPILER   intel
   setenv NTHRDS 4
   setenv NTASKS_PER_INST 64
   setenv NCORES_PER_NODE  32
   set scheduler = "PBS"

   setenv CSMDATA  /lustre/atlas1/cli900/world-shared/cesm/inputdata
   setenv CCSMROOT $HOME/codes/${CCSMTAG}
   setenv PTMP     /lustre/atlas/proj-shared/cli112/huiwan/$taskname
   setenv EXELOC   $PTMP/exe/$taskname

else if ($MACH == "cori-knl") then

   setenv CESM_PROJ  m3089
   setenv PROJECT $CESM_PROJ
   setenv COMPILER   intel

   setenv COSTPES 0
   setenv MPNT    272
   setenv MMTPN   64
   setenv PPND    136

   setenv NTHRDS 4
   setenv NINST 1
   setenv PSTRID 1
   setenv NTASKS_PER_INST 340

   setenv CSMDATA  /project/projectdirs/acme/inputdata/
   setenv PTMP   $WORKDIR/$taskname
   setenv EXELOC $PTMP/exe

else if ($MACH == "cori-haswell") then

   setenv CESM_PROJ  m3089
   setenv PROJECT $CESM_PROJ
   setenv COMPILER   intel
   setenv NTHRDS 4
   setenv NTASKS_PER_INST 64

   setenv CSMDATA  /project/projectdirs/acme/inputdata/
   setenv PTMP   $WORKDIR/$taskname
   setenv EXELOC $PTMP/exe

else if ($MACH == "constance") then

   setenv CESM_PROJ uq_climate
   setenv PROJECT $CESM_PROJ
   setenv COMPILER   intel
   setenv NTHRDS 1
   setenv NTASKS_PER_INST  288
   setenv NCORES_PER_NODE  24
   set scheduler = "SLURM"

   setenv CSMDATA  /pic/projects/climate/csmdata/ 
   setenv CCSMROOT $WORKDIR/codes/${CCSMTAG}
   setenv PTMP     $WORKDIR/$taskname
   setenv EXELOC   $PTMP/exe/

else if ($MACH == "compy") then

   setenv CESM_PROJ ESMD
   setenv PROJECT $CESM_PROJ
   setenv COMPILER   intel
   setenv NTHRDS 4
   setenv NTASKS_PER_INST  96 #576

   setenv CSMDATA  /compyfs/zhan391/acme_init/csmdata 
   #/compyfs/inputdata
   setenv PTMP     $WORKDIR/$taskname
   setenv EXELOC   $PTMP/exe
   #set initDir = /compyfs/zhan391/acme_init/ne30_FC5_init/

else

   echo Specify paths for MACH $MACH ! Abort.
   exit

endif
   mkdir -p $PTMP 

#----------------------------------------------------------------------------------
# We will compile the model just once, and use the same executable for all ensemble
# members. Personally, I prefer to create a separate "case" just for the compilation.

   set execase = compile_${CCSMTAG}_${COMPSET}_${RESOLUTION}_${MACH}_${COMPILER}_${NTASKS_PER_INST}proc

if ($debug == 'TRUE') then
   set execase = ${execase}_debug
endif

setenv EXEDIR ${EXELOC}/$execase

#----------------------------------------------------------------------------
# Later in this script, one "case" is created for every ensemble member (realization)
# of the reference simulations, trusted simulations, and test simulations. 
# We need to make sure that the same compile-time options (e.g., compset, 
# compiler, PE layout) are used for those "simulation cases" as well as 
# for the "compilation case". $case_setup_script contains the "create_case"
# command, specification of RUNDIR, EXEROOT, and PE layout, as well as 
# the "cesm_setup" command. This script is sourced (i.e., used like a
# Fortran subroutine) when we create new case for compilation and simulations.

set case_setup_script = ${exp_setup_dir}/create_and_setup_bundled_dart.csh
set driver_script_dir = ${exp_setup_dir}

####################################################################
# Compile model (just once)
####################################################################
if ($compile_model > 0) then

   setenv CASE     $execase
   setenv CASEROOT ${PTMP}/cases/$CASE
   setenv RUNDIR   ${PTMP}/run/$CASE

   # Create and set up new case

   echo
   echo Start to create case
   echo
   source ${case_setup_script}
   echo
   echo Finished creating case
   echo

   # Build the model

   cd $CASEROOT

#   ./xmlchange -file env_build.xml -id GMAKE_J -val '8'
   ./xmlchange -file env_build.xml -id DEBUG   -val $debug

   echo
   echo Start to build model 
   echo
   if ( $postCIME <= 2 ) then
      ./$CASE.build
    else
      ./case.build
    endif

    echo 
    echo Finished building the model.
    echo

endif

#####################################################################
# Conduct simulation(s)
#####################################################################
if ($run_model > 0) then

   #---------------------------------
   cd $driver_script_dir

   set template = './cnvg_template_'$testconfig'.csh'

   #---------------------------------
   # GROUP and DTIME loops
   #---------------------------------
   set ig = $igS
   while ( $ig <= $igE )
    set group = ${groupList[$ig]}
    set nl_file = ${driver_script_dir}/"namelist_files/cam_nl_"${group}

    set idtime = 1
    while ( $idtime <= $ndtime )

       set dtime  = $dtimeList[$idtime]
   
      #------------------------------------
      # create script for each realization
      #------------------------------------
   set case  = ${group}_${RESOLUTION}_"DT"`printf "%04d" ${dtime}`
 # set case  = ${group}_${RESOLUTION}_${CCSMTAG}_${COMPSET}_${COMPILER}_${MACH}_"DT"`printf "%04d" ${dtime}`"_"${testconfig}
   set caseroot = ${PTMP}/cases/$case
   set rundir   = ${PTMP}/run/$case

###copy the required restart files to the directory for branch run#####
   if ( ! -d $rundir ) then
      mkdir -p $rundir
   endif
   set camiin   = "$run_refcase.cam.i.$run_refdate-$run_reftod.nc"
   set camrin   = "$run_refcase.cam.r.$run_refdate-$run_reftod.nc"
   set camrsin  = "$run_refcase.cam.rs.$run_refdate-$run_reftod.nc"
  #set camrhin  = "$run_refcase.cam.rh0.$run_refdate-$run_reftod.nc"
   set lndiin   = "$run_refcase.clm2.r.$run_refdate-$run_reftod.nc"
  #set lndrin   = "$run_refcase.clm2.rh0.$run_refdate-$run_reftod.nc"
   set cicerin  = "$run_refcase.cice.r.$run_refdate-$run_reftod.nc"
   set cplrin   = "$run_refcase.cpl.r.$run_refdate-$run_reftod.nc"
   set docnrin  = "$run_refcase.docn.r.$run_refdate-$run_reftod.nc"
   set docnrsin = "$run_refcase.docn.rs1.$run_refdate-$run_reftod.bin"
   cp -r ${initDir}${camiin}  $rundir/${camiin}
   cp -r ${initDir}${camrin}  $rundir/${camrin}
   cp -r ${initDir}${camrsin} $rundir/${camrsin}
  #cp -r ${initDir}${camrhin} $rundir/${camrhin}
   cp -r ${initDir}${lndiin}  $rundir/${lndiin}
  #cp -r ${initDir}${lndrin}  $rundir/${lndrin}
   cp -r ${initDir}${cicerin} $rundir/${cicerin}
   cp -r ${initDir}${cplrin}  $rundir/${cplrin}
   cp -r ${initDir}${docnrin} $rundir/${docnrin}
   cp -r ${initDir}${docnrsin} $rundir/${docnrsin}

   echo "$lndiin"   >! $rundir/rpointer.lnd
   echo "$cicerin"  >! $rundir/rpointer.ice
   echo "$camrin"   >! $rundir/rpointer.atm
   echo "$cplrin"   >! $rundir/rpointer.drv
   echo "$docnrin"  >! $rundir/rpointer.ocn
   echo "$docnrsin" >> $rundir/rpointer.ocn

   set atm_init   = "$rundir/${camiin}"
   set lnd_init   = "$rundir/${lndiin}"
   
   set tmp_script = 'tmp_script_'`date +%F-%H%M%S-%N`
   
   cat > $tmp_script <<EOF
#!/bin/csh
date
setenv PROJECT         '$CESM_PROJ'
setenv CESM_PROJ       '$CESM_PROJ'
setenv CESM_EMAIL      '$CESM_EMAIL'

setenv MACH            '$MACH'
setenv COMPILER        '$COMPILER'
setenv RESOLUTION      '$RESOLUTION'
setenv NTASKS_PER_INST '$NTASKS_PER_INST'
setenv NINST           '$NINST'
setenv NTHRDS         '$NTHRDS'
         
setenv CCSMTAG  '$CCSMTAG'
setenv COMPSET  '$COMPSET'
setenv CCSMROOT '$CCSMROOT'
setenv PTMP     '$PTMP'
setenv EXEDIR   '$EXEDIR'
setenv CSMDATA  '$CSMDATA'
         
setenv CASE      '$case'
setenv CASEROOT  '$caseroot'
setenv RUNDIR    '$rundir'

setenv postCIME $postCIME

set case_setup_script = '$case_setup_script' 

set run_refcase = '$run_refcase'
set run_refdate = '$run_refdate'
set run_reftod  = '$run_reftod'

set dtime    = $dtime
set atm_init = $atm_init
set lnd_init = $lnd_init

set ncycle   = 1          # 1 cycle in total, i.e., no restart

#@ nlen = 3600 / $dtime   # run model for 60 minutes.
#set stop_o   = 'nsteps'
set nlen      = 6 # run model for 5 years
set stop_o    = 'nhours'

#@ nhtfrq     = 1800 / $dtime    # output every 30 minutes
set nhtfrq    = -6               # monthly average
set mfilt     = 1                # 1 time step per file 

set nl_file = $nl_file
EOF

      set jobscript = job_${case}.csh
      cat $tmp_script $template > $jobscript
      rm $tmp_script

      echo Created heavy-wgt script $jobscript

      # run the heavy-wgt script to create case for a simulation
      csh $jobscript

   @ idtime++
   end  #----------------------------------------------

  @ ig++
  end

endif
