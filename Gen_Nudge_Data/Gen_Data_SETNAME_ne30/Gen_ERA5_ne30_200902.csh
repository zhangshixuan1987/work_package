#! /bin/csh -f

#---------------------------------------
# Script to Generate YOTC Nudging data
#---------------------------------------

#===================================================================
# CONFIGRATION SECTION:
#===================================================================

# Set a REFERENCE (Starting) Date and the numbe rof days to process
#-------------------------------------------------------------------
set RUNNUM   = 02 
set REF_DATE = '20090201'
set NUM_DAYS = 28

# Set INPUT/OUTPUT/TMP directories, 
#-----------------------------------
set NAMELIST = './Config/Config_makeIC-'$RUNNUM'.nl'
set MYLOGDIR  = './LOG/LOG_002.'$RUNNUM'/'
set MYTMPDIR  = './TMP/TMP_002.'$RUNNUM'/'
set MYOUTDIR  = '/compyfs/zhan391/acme_init/ERA5_Reanalysis/era_to_e3sm/era_to_se/method4/'
set INPUTDIR  = '/compyfs/zhan391/acme_init/ERA5_Reanalysis/grb2nc/Feb/'

if ( ! -d $MYLOGDIR ) then
 mkdir -p $MYLOGDIR 
endif 

if ( ! -d $MYTMPDIR ) then
 mkdir -p $MYTMPDIR
endif
 
if ( ! -d $MYOUTDIR ) then
 mkdir -p $MYOUTDIR
endif

# Set ESMF options
#------------------------
set ESMF_interp = 'conserve'
set ESMF_pole   = 'none'
set ESMF_clean  = 'False'
set TMP_clean   = 'True'

# Set Processing options
#------------------------ 
set CASE                   = 'era5_ne30L72'
set DYCORE                 = 'se'
set PRECISION              = 'float'
set VORT_DIV_TO_UV         = 'False'
set SST_MASK               = 'False'
set ICE_MASK               = 'False'
set OUTPUT_PHIS            = 'True'
set REGRID_ALL             = 'True'
set ADJUST_STATE_FROM_TOPO = 'True'
set MASS_FIX               = 'False' #'True'

# Set files containig OUTPUT Grid structure and topography
#---------------------------------------------------------
set fname_grid_info        = '/compyfs/zhan391/acme_init/csmdata/atm/cam/inic/homme/cami_mam3_0000-10-01_ne30np4_L72_c160127.nc'
set fname_phis_output      = '/compyfs/zhan391/acme_init/ERA5_Reanalysis/era_to_e3sm/topo/USGS-gtopo30_ne30np4_16xdel2-PFC-consistentSGH.nc'
set ftype_phis_output      = 'SE_TOPOGRAPHY'

# Set INPUT filename format, type of file, and number of time records per file
#------------------------------------------------------------------------------
set fname   = ( none none none none none none )
set ftype   = ( none none none none none none )
set ftime   = ( none none none none none none )
set hoursec = ( 00000 03600 07200 10800 14400 18000 21600 25200 28800 32400 36000 39600 43200 46800 50400 54000 57600 61200 64800 68400 72000 75600 79200 82800)
set hourstr = ( 00    01   02   03   04   05   06   07   08   09   10   11   12   13   14   15   16   17   18   19   20   21   22   23)

set fname[1] = "era5_an_ml_YYYY-MM-DD_HH:00:00.nc"
set ftype[1] = "ECMWF_ERA5" #"Era_Interim_627.0_sc"
set ftime[1] = "1X"

# Set the OUTPUT fileds and the index of the file containing them
#------------------------------------------------------------------
set fields       = ( U V T PS Q )
set field_findex = ( 1 1 1 1  1 )

# Set the index of the file containing the INPUT topography
#------------------------------------------------------------------
set phis_findex  = 1


#===================================================================
# PROCESSING SECTION:
#   Loop over range of times, 
#   generate inpute namelist for each time, 
#   then generate the desired Nudging data
#===================================================================

# Construct the fields string
#-----------------------------
set fieldstr  = $fields[1]
set jj        = $field_findex[1]
@ jj = $jj - 1 
set findexstr = $jj
set ii = 1
set fnum = $#fields
while( $ii <  $fnum )
  @ ii = $ii + 1
  set jj = $field_findex[$ii]
  @ jj = $jj - 1
  set fieldstr  = `echo $fieldstr','$fields[$ii]`
  set findexstr = `echo $findexstr','$jj`
end

set jj = $phis_findex
@ jj = $jj - 1 
set phisstr = $jj

# Build NCL fortran library
#---------------------------
#/usr/local/bin/WRAPIT         MAKEIC.stub MAKEIC.f90
setenv NCARG_ROOT  "/compyfs/zhan391/acme_init/ERA5_Reanalysis/era_to_e3sm/rgd_run_se/method4/Gen_Data_SEdycore/Gen_Data_SETNAME_ne30/ncl-6.4.0"
set ncl = ${NCARG_ROOT}/bin/ncl
set wrp = ${NCARG_ROOT}/bin/WRAPIT
$wrp MAKEIC.stub MAKEIC.f90

#
# Loop over range of times to process
#  For the first time true, force 
#  generation of ESMF weight datasets
#=====================================
set ESMF_clean  = 'True'
set NDAYS = 0
while ( $NDAYS < $NUM_DAYS )

  # Set current time values
  #------------------------
  set Yearstr = "`date --date=$REF_DATE+${NDAYS}day +%Y`"
  set Monstr  = "`date --date=$REF_DATE+${NDAYS}day +%m`"
  set Daystr  = "`date --date=$REF_DATE+${NDAYS}day +%d`"
  echo $Yearstr ' / ' $Monstr ' / ' $Daystr

  # Loop over the hourly values (4X daily)
  #========================================
  foreach hnum ( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 )
    
    # Set Values dependend upon $hnum, clean 
    # up TMP files at the end of each DAY
    #----------------------------------------------------------------------
    set datestr = $Yearstr$Monstr$Daystr$hoursec[$hnum]
    set LOGFILE = $MYLOGDIR'/LogNCL.'$Yearstr$Monstr$Daystr$hoursec[$hnum]
    if( $hnum == 24 ) then
      set TMP_clean   = 'True'
    else
      set TMP_clean   = 'False'
    endif

    # Create Namelist file for current time
    #----------------------------------------
    echo '                                       '             > $NAMELIST
    echo '! Generated Namelist for maekIC_se.ncl '            >> $NAMELIST
    echo '!-------------------                   '            >> $NAMELIST
    echo '&makeic_nl'                                         >> $NAMELIST
    echo ' MYTMPDIR="'$MYTMPDIR'"'                            >> $NAMELIST
    echo ' MYOUTDIR="'$MYOUTDIR'"'                            >> $NAMELIST
    echo ' INPUTDIR="'$INPUTDIR'"'                            >> $NAMELIST
    echo ' INPUTDIR="'$INPUTDIR'"'                            >> $NAMELIST
    echo '                                    '               >> $NAMELIST
    echo ' ESMF_interp="'$ESMF_interp'"'                      >> $NAMELIST
    echo ' ESMF_pole  ="'$ESMF_pole'"'                        >> $NAMELIST
    echo ' ESMF_clean ="'$ESMF_clean'"'                       >> $NAMELIST
    echo ' TMP_clean  ="'$TMP_clean'"'                        >> $NAMELIST
    echo '                                    '               >> $NAMELIST
    echo ' REF_DATE              ="'$REF_DATE'"'              >> $NAMELIST
    echo ' CASE                  ="'$CASE'"'                  >> $NAMELIST
    echo ' DYCORE                ="'$DYCORE'"'                >> $NAMELIST
    echo ' PRECISION             ="'$PRECISION'"'             >> $NAMELIST
    echo ' VORT_DIV_TO_UV        ="'$VORT_DIV_TO_UV'"'        >> $NAMELIST
    echo ' SST_MASK              ="'$SST_MASK'"'              >> $NAMELIST
    echo ' ICE_MASK              ="'$ICE_MASK'"'              >> $NAMELIST
    echo ' OUTPUT_PHIS           ="'$OUTPUT_PHIS'"'           >> $NAMELIST
    echo ' REGRID_ALL            ="'$REGRID_ALL'"'            >> $NAMELIST
    echo ' ADJUST_STATE_FROM_TOPO="'$ADJUST_STATE_FROM_TOPO'"'>> $NAMELIST
    echo ' MASS_FIX              ="'$MASS_FIX'"'              >> $NAMELIST
    echo '                                    '               >> $NAMELIST
    echo ' fname_phis_output     ="'$fname_phis_output'"'     >> $NAMELIST
    echo ' ftype_phis_output     ="'$ftype_phis_output'"'     >> $NAMELIST
    echo ' fname_grid_info       ="'$fname_grid_info'"'       >> $NAMELIST
    echo '                                    '               >> $NAMELIST
    echo ' fields       ="'$fieldstr'"'                       >> $NAMELIST
    echo ' source_files ="'$findexstr'"'                      >> $NAMELIST
    echo ' fname_phis_in='$phisstr                            >> $NAMELIST
    echo '                                    '               >> $NAMELIST
    set ii = 0
    set fnum = $#fname
    while( $ii <  $fnum )
      set jj = $ii
      @ ii = $ii + 1
      set file = $fname[$ii]
      set file = `echo $file | sed 's/YYYY/'$Yearstr'/g'`
      set file = `echo $file | sed 's/MM/'$Monstr'/g'`
      set file = `echo $file | sed 's/DD/'$Daystr'/g'`
      if( $ftime[$ii] == "4X" ) then
        set file = `echo $file | sed 's/HH/'$hourstr[1]'/'`
        set file = `echo $file | sed 's/HH/'$hourstr[4]'/'`
      else
        set file = `echo $file | sed 's/HH/'$hourstr[$hnum]'/g'`
      endif
      echo " fname$jj="'"'$file'"'                      >> $NAMELIST
    end
    echo '                                    '               >> $NAMELIST
    set ii = 0
    set fnum = $#fname
    while( $ii <  $fnum )
      set jj = $ii
      @ ii = $ii + 1
      echo " ftype$jj="'"'$ftype[$ii]'"'                      >> $NAMELIST
    end
    echo '                                    '               >> $NAMELIST
    set ii = 0
    set fnum = $#fname
    while( $ii <  $fnum )
      set jj = $ii
      @ ii = $ii + 1
      if( $fname[$ii] == "none" ) then
        echo " fdate$jj="'"-1"'                               >> $NAMELIST
      else
        echo " fdate$jj="'"'$datestr'"'                       >> $NAMELIST
      endif
    end
    echo "/"                                                  >> $NAMELIST

    # Execute NCL program to process the data
    #----------------------------------------
    cat $NAMELIST                            >& $LOGFILE.cfg
    $ncl makeIC_se_002.ncl NameNumber=$RUNNUM >& $LOGFILE

  # End Loop over the hourly values (4X daily)
  #===========================================
   set ESMF_clean  = 'False'
  end # foreach hnum ( 1 2 3 4 )

# End loop over time
#=======================
@ NDAYS = $NDAYS + 1
end # while ( $NDAYS < $NUM_DAYS )

