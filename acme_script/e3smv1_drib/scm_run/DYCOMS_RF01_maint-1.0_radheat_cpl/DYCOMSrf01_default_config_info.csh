#!/bin/csh

####set required variables for this case
# Aerosol specification
# Options include:
# 1) cons_droplet (sets cloud liquid and ice concentration
#                  to a constant)
# 2) prescribed (uses climatologically prescribed aerosol 
#                concentration)
setenv init_aero_type prescribed

# Case specific information kept here
setenv lat  31.5 # latitude  
setenv lon  239.0 # longitude
setenv do_iop_srf_prop    .true. # Use surface fluxes in IOP file?
setenv do_scm_relaxation  .false. # Relax case to observations?
setenv do_turnoff_swrad   .true. # Turn off SW calculation
setenv do_turnoff_lwrad   .false. # Turn off LW calculation
setenv do_turnoff_precip  .true. # Turn off precipitation
setenv micro_nccons_val   55.0D6 # cons_droplet value for liquid
setenv micro_nicons_val   0.0001D6 # cons_droplet value for ice
setenv startdate   1999-07-10 # Start date in IOP file
setenv start_in_sec   0 # start time in seconds in IOP file
#setenv stop_option   nhours
#setenv stop_n        96 # the forcing data tsec = 0, 345600 ;
setenv iop_file  DYCOMSrf01_iopfile_4scam.nc #IOP file name

# Location of IOP file
setenv iop_path atm/cam/scam/iop
# Prescribed aerosol file path and name
setenv presc_aero_path  atm/cam/chem/trop_mam/aero
setenv presc_aero_file  mam4_0.9x1.2_L72_2000clim_c170323.nc

###setup land surface data####
#setenv lnd_surf_file  lnd/clm2/surfdata_map/surfdata_64x128_simyr2000_c170111.nc
# End Case specific stuff here
#
# $testconfig is used to pick the appropriate template file that contains certain
# run-time options. fullmodel_rad means with all the default parameterizations,
# and with radiation called every other time step.
#set testconfig   = 'e3smv1_scm'

