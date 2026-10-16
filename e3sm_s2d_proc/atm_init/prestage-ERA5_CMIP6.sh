#!/bin/bash 
# usage: ./prestage-ERA5.sh /glade/work/${LOGNAME}/sewx/ECMWF/ 1996011800
# where you want the data file for Jan 18, 1996 at 00Z

DIR=$1
DATE=$2
THISDIR=$PWD

yearstr=${DATE:0:4}
monthstr=${DATE:4:2}
daystr=${DATE:6:2}
cyclestr=${DATE:8:2}

set -e
mkdir -p $DIR
cd $DIR

INPUTDIR="/global/cfs/projectdirs/m3522/cmip6/ERA5"
YYMM=${yearstr}${monthstr}
YYMMDD=${yearstr}${monthstr}${daystr}
if [[ ${monthstr} == "01" || ${monthstr} == "03" || ${monthstr} == "05"  || \
      ${monthstr} == "07" || ${monthstr} == "08" || ${monthstr} == "10"  || \
      ${monthstr} == "12"  ]] ; then 
  monday="31"
elif [[ ${monthstr} == "02" ]] ; then 
  monday="28" 
else 
  monday="30"    	
fi

fields=( phis t2m    siconc sst  msl sp u v w t q r z ciwc clwc ) 
fldera=(    Z VAR_2T CI     SSTK MSL SP U V W T Q R Z CIWC CLWC )
fname=( "/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc" \
        "/e5.oper.an.sfc/${YYMM}/e5.oper.an.sfc.128_167_2t.ll025sc.${YYMM}0100_${YYMM}${monday}23.nc" \
        "/e5.oper.an.sfc/${YYMM}/e5.oper.an.sfc.128_031_ci.ll025sc.${YYMM}0100_${YYMM}${monday}23.nc" \
        "/e5.oper.an.sfc/${YYMM}/e5.oper.an.sfc.128_034_sstk.ll025sc.${YYMM}0100_${YYMM}${monday}23.nc" \
        "/e5.oper.an.sfc/${YYMM}/e5.oper.an.sfc.128_151_msl.ll025sc.${YYMM}0100_${YYMM}${monday}23.nc" \
        "/e5.oper.an.sfc/${YYMM}/e5.oper.an.sfc.128_134_sp.ll025sc.${YYMM}0100_${YYMM}${monday}23.nc" \
	"/e5.oper.an.pl/${YYMM}/e5.oper.an.pl.128_131_u.ll025uv.${YYMMDD}00_${YYMMDD}23.nc" \
        "/e5.oper.an.pl/${YYMM}/e5.oper.an.pl.128_132_v.ll025uv.${YYMMDD}00_${YYMMDD}23.nc" \
        "/e5.oper.an.pl/${YYMM}/e5.oper.an.pl.128_135_w.ll025sc.${YYMMDD}00_${YYMMDD}23.nc" \
        "/e5.oper.an.pl/${YYMM}/e5.oper.an.pl.128_130_t.ll025sc.${YYMMDD}00_${YYMMDD}23.nc" \
        "/e5.oper.an.pl/${YYMM}/e5.oper.an.pl.128_133_q.ll025sc.${YYMMDD}00_${YYMMDD}23.nc" \
        "/e5.oper.an.pl/${YYMM}/e5.oper.an.pl.128_157_r.ll025sc.${YYMMDD}00_${YYMMDD}23.nc" \
        "/e5.oper.an.pl/${YYMM}/e5.oper.an.pl.128_129_z.ll025sc.${YYMMDD}00_${YYMMDD}23.nc" \
        "/e5.oper.an.pl/${YYMM}/e5.oper.an.pl.128_247_ciwc.ll025sc.${YYMMDD}00_${YYMMDD}23.nc" \
        "/e5.oper.an.pl/${YYMM}/e5.oper.an.pl.128_246_clwc.ll025sc.${YYMMDD}00_${YYMMDD}23.nc" \
       ) 

nflds=${#fields[@]}

i=0
while [ $i -lt $nflds ];do
 
 varo=${fields[$i]}
 vari=${fldera[$i]} 
 fli=${INPUTDIR}/${fname[$i]}
 flo=ERA5_${yearstr}${monthstr}${daystr}${cyclestr}.nc
 echo $varo $vari $fli
 echo $fli 
 echo $flo
 if [[ $varo == "phis" ]];then 
   cp -r $fli $flo
   ncrename -v $vari,$varo $flo
 else 
   ncks -d time,"${yearstr}-${monthstr}-${daystr} ${cyclestr}:00:0.0" $fli  tmp.nc
   ncks -A -v ${vari} tmp.nc $flo
   ncrename -v $vari,$varo $flo
   rm -rvf tmp.nc
 fi 
 i=`expr $i + 1`
done 

exit

