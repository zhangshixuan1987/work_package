#!/bin/csh

set CLIM_SST = "/compyfs/inputdata/atm/cam/sst/sst_HadOIBl_bc_1x1_clim_c101029.nc"
set AMIP_SST = "/compyfs/inputdata//ocn/docn7/SSTDATA/sst_ice_CMIP6_DECK_E3SM_1x1_c20180213.nc"

#extact the interested year##
#module load cdo

cdo selyear,2009 $AMIP_SST 2009.nc
cdo selyear,2010 $AMIP_SST 2010.nc
cdo selyear,2011 $AMIP_SST 2011.nc

cp -r 2009.nc 2009_tmp.nc
cp -r 2010.nc 2010_tmp.nc
cp -r 2011.nc 2011_tmp.nc

ncks -A -v SST_cpl ${CLIM_SST} 2009.nc
ncks -A -v SST_cpl ${CLIM_SST} 2010.nc
ncks -A -v SST_cpl ${CLIM_SST} 2011.nc

ncks -A -v SST_cpl_prediddle ${CLIM_SST} 2009.nc
ncks -A -v SST_cpl_prediddle ${CLIM_SST} 2010.nc
ncks -A -v SST_cpl_prediddle ${CLIM_SST} 2011.nc

ncks -A -v ice_cov ${CLIM_SST} 2009.nc
ncks -A -v ice_cov ${CLIM_SST} 2010.nc
ncks -A -v ice_cov ${CLIM_SST} 2011.nc

ncks -A -v ice_cov_prediddle ${CLIM_SST} 2009.nc
ncks -A -v ice_cov_prediddle ${CLIM_SST} 2010.nc
ncks -A -v ice_cov_prediddle ${CLIM_SST} 2011.nc

ncks -A -v time  2009_tmp.nc 2010.nc 
ncks -A -v time  2010_tmp.nc 2010.nc
ncks -A -v time  2010_tmp.nc 2010.nc

ncks -A -v time_bnds 2009_tmp.nc 2009.nc
ncks -A -v time_bnds 2010_tmp.nc 2010.nc
ncks -A -v time_bnds 2011_tmp.nc 2011.nc


set outfile = sst_ice_CMIP6_DECK_E3SM_1x1_CLIM_c20210618.nc

rm -rvf $outfile

ncrcat 2009.nc 2010.nc 2011.nc $outfile
rm -rvf 2009*.nc 2010*.nc 2011*.nc

