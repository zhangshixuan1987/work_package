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

#####determine which queue for the job###########
setenv queue  'debug'
setenv queue  'regular'

setenv submit_to_queue  'FALSE'
#setenv submit_to_queue  'TRUE'
setenv walltime         '03:00:00'

setenv debug  'TRUE'
setenv debug  'FALSE'

set wkdrnam = "F_SCAM5_201806"
#-----------------------------------------------------------------------------------------
# $taskname is used to organize the exe/run/case directories. Not used for other purposes.
set taskname = 'DYCOMS-RF01_MEAN' 
# User enter any needed modules to load or use below
#  EXAMPLE:
#module load python

####set required variables for this case
# Aerosol specification
# Options include:
#  1) cons_droplet (sets cloud liquid and ice concentration
#                   to a constant)
#  2) prescribed (uses climatologically prescribed aerosol 
#                 concentration)
setenv init_aero_type prescribed

# Set the dynamical core
#   1) Select "Eulerian" IF you are running E3SMv1 release code 
#      (or master branch code before March 10,2019)
#   2) Select "SE" IF you are running code from E3SM master branch that
#     is AFTER March 10,2019
setenv dycore Eulerian

#User enter any needed modules to load or use below
#module load python/2.7.5

# Case specific information kept here
setenv lat  31.5 # latitude  
setenv lon  239.0 # longitude
setenv do_iop_srf_prop    .true. # Use surface fluxes in IOP file?
setenv do_scm_relaxation  .false. # Relax case to observations?
setenv do_turnoff_swrad   .false. # Turn off SW calculation
setenv do_turnoff_lwrad   .false. # Turn off LW calculation
setenv do_turnoff_precip  .false. # Turn off precipitation
setenv micro_nccons_val   55.0D6 # cons_droplet value for liquid
setenv micro_nicons_val   0.0001D6 # cons_droplet value for ice
setenv startdate   1999-07-10 # Start date in IOP file
setenv start_in_sec   0 # start time in seconds in IOP file
setenv stop_option   nhours 
setenv stop_n  96 # the forcing data tsec = 0, 345600 ;
setenv iop_file  DYCOMSrf01_iopfile_4scam.nc #IOP file name

# Location of IOP file
setenv iop_path atm/cam/scam/iop
# Prescribed aerosol file path and name
setenv presc_aero_path  atm/cam/chem/trop_mam/aero
setenv presc_aero_file  mam4_0.9x1.2_L72_2000clim_c170323.nc
###setup land surface data####
#setenv lnd_surf_file  lnd/clm2/surfdata_map/surfdata_64x128_simyr2000_c170111.nc
# End Case specific stuff here

# $testconfig is used to pick the appropriate template file that contains certain
# run-time options. fullmodel_rad means with all the default parameterizations,
# and with radiation called every other time step.
#set testconfig   = 'e3smv1_scm'

#-----------------------
set exp_setup_dir = `pwd`

#-----------------------------
if ($dycore == Eulerian) then
  setenv RESOLUTION  T42_T42
endif

if ($dycore == SE) then
  setenv RESOLUTION ne4_ne4
endif

setenv MACH     compy

# COSP, set to false unless user really wants it
setenv do_cosp  false

###################################################################
# This script does NOT clone the model code from the ACME repo.
# Instead, it assumes that the code is located at
# $HOME/codes/${CCSMTAG}/
# The code I used was the branch huiwanpnnl/atm/shcu-pbl-namelist-aero 
# of the ACME code repo.

setenv CCSMTAG ACME_master_20180611
setenv WORKDIR /compyfs/zhan391
setenv CCSMROOT /compyfs/zhan391/code/${CCSMTAG}
 
setenv postCIME 5

###################################################################
# Conduct simulations with different models or time step sizes.
# Uncomment one of the three blocks at a time
#------------------------------------------------------------------
# Reference solution using default model and 1-second time step. 
#------------------------------------------------------------------
set groupList = ("full_precp_and_swrad")
set dtimeList = (1200) #(1200 300 80 20 8 4 2 1)

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

set irS = 1   # ensemble member: start index
set irE = 1   # ensemble member: end   index

#====================================================================
# Paths to source code and model input/output
#====================================================================
set scheduler = "PBS"   # default

setenv CESM_EMAIL shixuan.zhang@pnnl.gov
setenv COMPSET F_SCAM5
setenv NINST  1


if ($MACH == "titan") then

   setenv CESM_PROJ cli112
   setenv PROJECT $CESM_PROJ
   setenv COMPILER   pgi
   setenv COMPILER   intel
   setenv NTHRDS 1
   setenv NTASKS_PER_INST 1
   setenv NCORES_PER_NODE 1
   set scheduler = "PBS"

   setenv CSMDATA  /lustre/atlas1/cli900/world-shared/cesm/inputdata
   setenv CCSMROOT $HOME/codes/${CCSMTAG}
   setenv PTMP     /lustre/atlas/proj-shared/cli112/huiwan/e3sm_v1_scm_fnl/$taskname
   setenv EXELOC   $PTMP/exe/$taskname
   set initDir = /lustre/atlas1/cli112/proj-shared/huiwan/cesm_input/FC5AV1C-L_init/

else if ($MACH == "eos") then

   setenv CESM_PROJ cli112
   setenv PROJECT $CESM_PROJ
   setenv COMPILER   pgi
   setenv COMPILER   intel
   setenv NTHRDS 1
   setenv NTASKS_PER_INST 1
   setenv NCORES_PER_NODE 1
   set scheduler = "PBS"

   setenv CSMDATA  /lustre/atlas1/cli900/world-shared/cesm/inputdata
   setenv CCSMROOT $HOME/codes/${CCSMTAG}
   setenv PTMP     /lustre/atlas/proj-shared/cli112/huiwan/$taskname
   setenv EXELOC   $PTMP/exe/$taskname
   set initDir = /lustre/atlas1/cli112/proj-shared/huiwan/cesm_input/FC5AV1C-L_init/ 

else if ($MACH == "cori-knl") then

   setenv CESM_PROJ  m3089
   setenv PROJECT $CESM_PROJ
   setenv COMPILER   intel
   setenv NTHRDS 1
   setenv NINST  1
   setenv NTASKS_PER_INST 1

   setenv CSMDATA  /project/projectdirs/acme/inputdata/
   setenv PTMP   $WORKDIR/$wkdrnam
   setenv EXELOC $PTMP/exe
   set initDir = /global/cscratch1/sd/zhan391/acme_input/FC5AV1C-L_init_201712/

else if ($MACH == "cori-haswell") then

   setenv CESM_PROJ  m3089
   setenv PROJECT $CESM_PROJ
   setenv COMPILER   intel
   setenv NTHRDS 1
   setenv NINST 1
   setenv NTASKS_PER_INST 1

   setenv CSMDATA  /project/projectdirs/acme/inputdata/
   setenv PTMP   $WORKDIR/$wkdrnam
   setenv EXELOC $PTMP/exe
   set initDir = /global/cscratch1/sd/zhan391/acme_input/FC5AV1C-L_init_201712/

else if ($MACH == "constance") then

   setenv CESM_PROJ uq_climate
   setenv PROJECT $CESM_PROJ
   setenv COMPILER   intel
   setenv NTHRDS 1
   setenv NTASKS_PER_INST  1
   setenv NCORES_PER_NODE  1
   set scheduler = "SLURM"

   setenv WORKDIR /pic/projects/uq_climate/zhan391
   setenv CSMDATA  /pic/projects/climate/csmdata/ 
   setenv CCSMROOT $WORKDIR/codes/${CCSMTAG}
   setenv PTMP     $WORKDIR/wkdrnam/$taskname
   setenv EXELOC   $PTMP/exe/
   set initDir = /pic/projects/uq_climate/wanh895/acme_input/FC5AV1C-L_init_201712/

else if ($MACH == "compy") then

   setenv CESM_EMAIL shixuan.zhang@pnnl.gov
   setenv PROJECT   ESMD
   setenv CESM_PROJ $PROJECT

   setenv COMPILER   intel
   setenv NTHRDS 1
   setenv NTASKS_PER_INST  1
   setenv NCORES_PER_NODE  1
   set scheduler = "PBS"

   setenv WORKDIR /compyfs/zhan391
   setenv CSMDATA  /compyfs/inputdata/
   setenv CCSMROOT $WORKDIR/code/${CCSMTAG}
   setenv PTMP     $WORKDIR/$wkdrnam/$taskname
   setenv EXELOC   $PTMP/exe/
   set initDir = /compyfs/zhan391/acme_init/ne30_FC5_init/

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

set case_setup_script = ${exp_setup_dir}/create_and_setup_bundled_case.csh
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

   ./xmlchange -file env_build.xml -id GMAKE_J -val '8'
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

   #---------------------------------
   # GROUP and DTIME loops
   #---------------------------------
   set ig = $igS
   while ( $ig <= $igE )
    set group = ${groupList[$ig]}
    set template = './cnvg_template_'$group'.csh'

    set idtime = 1
    while ( $idtime <= $ndtime )

       set dtime  = $dtimeList[$idtime]
   
      #------------------------------------
      # create script for each realization
      #------------------------------------
      set ir = $irS
      while ( $ir <= $irE )
      
         set irstring = `printf "%02d" ${ir}`
   
         set case     = ${group}_${taskname}_${RESOLUTION}_${CCSMTAG}_${COMPSET}_${COMPILER}_${MACH}_"DT"`printf "%04d" ${dtime}`"_"${irstring}  #_${testconfig}
         set caseroot = ${PTMP}/cases/$case
         set rundir   = ${PTMP}/run/$case
   
         set tmp_script = 'tmp_script_'`date +%F-%H%M%S-%N`
   
         cat > $tmp_script <<EOF
#!/bin/csh
date
setenv PROJECT         '$CESM_PROJ'
setenv CESM_PROJ       '$CESM_PROJ'
setenv CESM_EMAIL      '$CESM_EMAIL'

setenv do_cosp         '$do_cosp'
setenv dycore          '$dycore'
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

####setup the ouptput frequency###
set dtime    = $dtime
set ncycle   = 1                         # 1 cycle in total, i.e., no restart
@ nlen       = $stop_n * 3600 / $dtime   # run model for 60 minutes.
set stop_o   = 'nsteps'
set nhtfrq   = 1
@ mfilt      = $stop_n * 3600 / $dtime + 1           # 1 time step per file

EOF

      set jobscript = job_${case}.csh
      cat $tmp_script $template > $jobscript
      rm $tmp_script

      echo Created heavy-wgt script $jobscript

      # run the heavy-wgt script to create case for a simulation
      csh $jobscript

      @ ir++
      end  #----------------------------------------------


   @ idtime++
   end  #----------------------------------------------

  @ ig++
  end


endif
