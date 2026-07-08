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
set taskname = cam5_datass_test

#-----------------------
set exp_setup_dir = `pwd`

#-----------------------------
setenv RESOLUTION ne30_ne30
setenv MACH       cori-knl

###################################################################
# This script does NOT clone the model code from the ACME repo.
# Instead, it assumes that the code is located at
# $HOME/codes/${CCSMTAG}/
# The code I used was the branch huiwanpnnl/atm/shcu-pbl-namelist-aero 
# of the ACME code repo.

setenv CCSMTAG E3SM_MAINT1.0 
setenv CCSMROOT /compyfs/zhan391/code/$CCSMTAG
setenv WORKDIR /global/cscratch1/sd/zhan391
 
###################################################################
set groupList = ("Hindcast_5d")
set ngroups = $#groupList

set dtimeList = (1800)
set ndtime = $#dtimeList

set igS = 1
set igE = $ngroups

set irS = 1    # ensemble member: start index
set irE = 30   # ensemble member: end   index

# Bundle the ensemble simulations into $nbundles PBS jobs. 
# PBS scripts for bundled simulations are created by this script.
# If you would like to submit each simulation as an individual job,
# set nbundles to a negative number, e.g., -1
set nbundles = -1 

#====================================================================
# Paths to source code and model input/output
#====================================================================
setenv CESM_EMAIL shixuan.zhang@pnnl.gov
setenv PROJECT   ESMD
setenv CESM_PROJ $PROJECT
setenv postCIME 5

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
setenv CSMDATA /compyfs/zhan391/acme_init/csmdata

## the sorce of spun up initial conditions ##
#setup branch run information
set dtime1      = 1800
set run_refcase = DACAM5_ENS80_${COMPSET}_f19_g16_DT`printf "%04d" ${dtime1}`
set inputdat    = /compyfs/zhan391/e3sm_scratch/archive/DACAM5_ENS80_F20TRC5-AMIP_f19_g16_DT1800/rest/
set lnd_surf    = "${CSMDATA}/lnd/clm2/surfdata_map/surfdata_1.9x2.5_simyr1850_c180306.nc"
set lnd_use     = "${CSMDATA}/lnd/clm2/surfdata_map/landuse.timeseries_1.9x2.5_hist_simyr1850-2015_c180306.nc"

 
setenv PTMP     $WORKDIR/$taskname
setenv EXELOC   $PTMP/exe/
set initDir = /pic/projects/uq_climate/wanh895/acme_input/FC5AV1C-L_init_201712/

mkdir -p $PTMP 

#----------------------------------------------------------------------------------
# We will compile the model just once, and use the same executable for all ensemble
# members. Personally, I prefer to create a separate "case" just for the compilation.

if ($NINST > 1) then
   set execase = compile_${CCSMTAG}_${COMPSET}_${RESOLUTION}_${MACH}_${COMPILER}_${NINST}x${NTASKS_PER_INST}x${NTHRDS}bundle
else
   set execase = compile_${CCSMTAG}_${COMPSET}_${RESOLUTION}_${MACH}_${COMPILER}_${NTASKS_PER_INST}x${NTHRDS}threads
endif

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

set case_setup_script = ${exp_setup_dir}/create_and_setup_bundled_hind.csh
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
   ./case.build

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

   set template = './cam5_template_hindcast.csh'

   #---------------------------------
   # GROUP and DTIME loops
   #---------------------------------
   set ig = $igS
   while ( $ig <= $igE )
    set group = ${groupList[$ig]}

    set idtime = 1
    while ( $idtime <= $ndtime )

       set dtime  = $dtimeList[$idtime]

      #------------------------------------
      # create script for each realization
      #------------------------------------
      set ir = $irS
      while ( $ir <= $irE )
 
         set irstring = `printf "%02d" ${ir}`
   
         set case     = ${group}_${RESOLUTION}_${CCSMTAG}_${COMPSET}_${COMPILER}_${MACH}_"DT"`printf "%04d" ${dtime}`"_"${irstring}
         set caseroot = ${PTMP}/cases/$case
         set rundir   = ${PTMP}/run/$case

         setenv CASE     $case
         setenv CASEROOT $caseroot
         setenv RUNDIR   $rundir

         ##set the branch run information###
         if ( $ir == 1 )then
          set run_refdate = "2009-01-30"
          set run_reftod  = "00000"
         endif
         if ( $ir == 2 )then
          set run_refdate = "2009-01-31"
          set run_reftod  = "00000"
         endif
         if ( $ir >= 3 )then
          @ irr = $ir - 2
          set irstring = `printf "%02d" ${irr}`
          set run_refdate = "2009-02-${irstring}"
          set run_reftod  = "00000"
         endif

         if ( ! -d $RUNDIR ) then
          mkdir -p $RUNDIR
         endif

         set ininst      = 1
         set inst_string = `printf "_%04d" ${ininst}`
         #echo $inst_string
         ###copy the required restart files to the directory for branch run#####
         set camiin   = "$run_refcase.cam.i.$run_refdate-$run_reftod.nc"
         set camrin   = "$run_refcase.cam.r.$run_refdate-$run_reftod.nc"
         set camrsin  = "$run_refcase.cam.rs.$run_refdate-$run_reftod.nc"
        #set camrhin  = "$run_refcase.cam.rh0.$run_refdate-$run_reftod.nc"
         set lndiin   = "$run_refcase.clm2.r.$run_refdate-$run_reftod.nc"
        #set lndrin   = "$run_refcase.clm2.rh0.$run_refdate-$run_reftod.nc"
         set cicerin  = "$run_refcase.cice.r.$run_refdate-$run_reftod.nc"
        #set cplrin   = "$run_refcase.cpl.r.$run_refdate-$run_reftod.nc"
         set docnrin  = "$run_refcase.docn.rs1.$run_refdate-$run_reftod.nc"
         set docnrsin = "$run_refcase.docn.rs1.$run_refdate-$run_reftod.bin"
         set cplrin   = "$run_refcase.cpl.r.$run_refdate-$run_reftod.nc"

         ln -sf ${inputdat}${run_refdate}-${run_reftod}/"$run_refcase.cam${inst_string}.i.$run_refdate-$run_reftod.nc"   $RUNDIR/${camiin}
         ln -sf ${inputdat}${run_refdate}-${run_reftod}/"$run_refcase.cam${inst_string}.r.$run_refdate-$run_reftod.nc"   $RUNDIR/${camrin}
         ln -sf ${inputdat}${run_refdate}-${run_reftod}/"$run_refcase.cam${inst_string}.rs.$run_refdate-$run_reftod.nc"  $RUNDIR/${camrsin}
        #ln -sf ${inputdat}${run_refdate}-${run_reftod}/"$run_refcase.cam${inst_string}.rh0.$run_refdate-$run_reftod.nc"  $RUNDIR/${camrhin}
         ln -sf ${inputdat}${run_refdate}-${run_reftod}/"$run_refcase.clm2${inst_string}.r.$run_refdate-$run_reftod.nc"   $RUNDIR/${lndiin}
        #ln -sf ${inputdat}${run_refdate}-${run_reftod}/"$run_refcase.clm2${inst_string}.rh0.$run_refdate-$run_reftod.nc"   $RUNDIR/${lndrin}
         ln -sf ${inputdat}${run_refdate}-${run_reftod}/"$run_refcase.cice${inst_string}.r.$run_refdate-$run_reftod.nc"  $RUNDIR/${cicerin}
        #ln -sf ${inputdat}${run_refdate}-${run_reftod}/"$run_refcase.cpl${inst_string}.r.$run_refdate-$run_reftod.nc"   $RUNDIR/${cplrin}
         ln -sf ${inputdat}${run_refdate}-${run_reftod}/"$run_refcase.docn${inst_string}.rs1.$run_refdate-$run_reftod.bin" $RUNDIR/${docnrsin}
         ln -sf ${inputdat}${run_refdate}-${run_reftod}/"$run_refcase.cpl.r.$run_refdate-$run_reftod.nc"                 $RUNDIR/${cplrin}

         echo "${camrin}"     >! $RUNDIR/rpointer.atm #${inst_string}
         echo "${lndiin}"     >! $RUNDIR/rpointer.lnd #${inst_string}
         echo "${cicerin}"    >! $RUNDIR/rpointer.ice #${inst_string}
         echo "${docnrin}"    >! $RUNDIR/rpointer.ocn #${inst_string}
         echo "${docnrsin}"   >> $RUNDIR/rpointer.ocn #${inst_string}
         echo "${cplrin}"     >! $RUNDIR/rpointer.drv

         set atm_init   = "$RUNDIR/${camiin}"
         set lnd_init   = "$RUNDIR/${lndiin}"
         set ice_init   = "$RUNDIR/${cicerin}"

   
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

set nbundles = $nbundles         
set case_setup_script = '$case_setup_script' 

set run_refcase = '$run_refcase'
set run_refdate = '$run_refdate'
set run_reftod  = '$run_reftod'

set dtime    = $dtime
set atm_init = $atm_init
set lnd_init = $lnd_init
set ice_init = $ice_init
set lnd_surf = $lnd_surf
set lnd_use  = $lnd_use
set ncycle   = 1      # no resubmit

set nlen     = 5      # run model for 5 days
set stop_o   = 'ndays'

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
