#!/bin/sh
rundir="/pic/projects/uq_climate/zhan391/e3sm_cam4/run_amwg_10yr/climo/"
datdir="/pic/projects/uq_climate/zhan391/VoglC_JAMES_2019_subgrid_reconstruction/Figure03_Figure06"
runnam="FC5AQUAP"
casnam="FC5AQUAP"
####start to copy and rename the file name############
config="dycore_mac"
config1="dycore_mac"

Group1="ql0_fmin0_1.9x2.5_gx1v6_E3SM_MASTER_F_intel_constance_DT1800_fullmodel_rad"

Group2="ql0_fmin0_sgrnew0_1.9x2.5_gx1v6_E3SM_MASTER_RECON_F_intel_constance_DT1800_fullmodel_rad"

Group3="ql0_fmin0_sgrnew_1.9x2.5_gx1v6_E3SM_MASTER_RECON_F_intel_constance_DT1800_fullmodel_rad"

ncks -v CLDTOT,LWCF,SWCF,lat,lon ${rundir}/${Group1}/${Group1}_ANN_means.nc ${datdir}/E3SM_CAM4_Baseline_climo_ANN_means_10yr.nc
ncks -v CLDTOT,LWCF,SWCF,lat,lon ${rundir}/${Group2}/${Group2}_ANN_means.nc ${datdir}/E3SM_CAM4_SGR000_climo_ANN_means_10yr.nc
ncks -v CLDTOT,LWCF,SWCF,lat,lon ${rundir}/${Group3}/${Group3}_ANN_means.nc ${datdir}/E3SM_CAM4_SGR114_climo_ANN_means_10yr.nc
