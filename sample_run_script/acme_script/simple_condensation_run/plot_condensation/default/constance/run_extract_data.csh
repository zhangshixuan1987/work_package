#!/bin/sh
rundir="/pic/projects/uq_climate/zhan391/cnvg_condens_genint_fnl/run"
datdir="/pic/projects/uq_climate/zhan391/WanH_JAMES_2019_simple_condensation/Figure_02"
runnam="FC5AQUAP"
casnam="FC5AQUAP"
####start to copy and rename the file name############
Group="RKZ_A1_B1_C2_ql19_lmt4_fmin1_ne30_ne30_"${runnam}"_intel_constance_"
config="dycore_only"
config1="dycore_only"

for hfcst in 3600 7200 14400 21600 43200;do 

dfhr=`printf "%05d" $hfcst`

timstr=".cam.h0.0001-01-01-"${dfhr}

for dt in 1 8 30 120 450 1800 ;do

echo $dt
dtim=`printf "%04d" $dt`

for nens in 1 2 3 4 5 6;do

sens=`printf "%02d" $nens`
case=${Group}"DT"${dtim}_${sens}_${config}
savenam=${Group}"DT"${dtim}_${sens}_${config}

rm -rvf ${datdir}/${savenam}${timstr}.nc

cp -rp ${rundir}/${case}/${case}${timstr}.nc ${datdir}/temp.nc

ncks -v T,PS,area,P0,lat,lon,hybi,hyai,LANDFRAC ${datdir}/temp.nc  ${datdir}/${savenam}${timstr}.nc

rm -rvf ${datdir}/temp.nc

done

done

done
