#!/bin/sh 
ulimit -S -c 0

case="F20TR_ne30pg2_EC30to60E2r2_EN01_compy "
nens=80
datadir="/compyfs/zhan391/dart_scratch_ne30pg2"
outnm="F20TR_ne30pg2_EC30to60E2r2_DART80_compy"
outdir="/compyfs/zhan391/initial_test_data/${outnm}/archive"

if [ ! -d $outdir ]; then 
  mkdir -p $outdir
fi 

i=1
while [ $i -le $nens ]; do 
  enstr=EN`printf "%02d" $i` 
  casnm=`echo $case | sed s/EN01/${enstr}/g `
  rundir="$datadir/$casnm/archive"
  for subdir in `cd ${rundir}; ls -d *`; do 
    echo $subdir
    if [ ! -d $outdir/$subdir ];then 
      if [ $subdir != "rest" ];then 
         mkdir -p $outdir/$subdir/hist
         mv ${rundir}/${subdir}/hist/* ${outdir}/${subdir}/hist 
      else 
         for xdr in `cd ${rundir}/${subdir}; ls -d *`;do
           if [ ! -d ${outdir}/${subdir}/${xdr} ]; then
              mkdir -p ${outdir}/${subdir}/${xdr}
           fi
           if [ -d ${rundir}/${subdir}/${xdr} ]; then 
             mv ${rundir}/${subdir}/${xdr}/* ${outdir}/${subdir}/${xdr}/
           fi 
         done
         mv ${rundir}/${subdir}/* ${outdir}/${subdir}/
      fi 
    else
      if [ $subdir != "rest" ];then
         mv ${rundir}/${subdir}/hist/* ${outdir}/${subdir}/hist
      else
         for xdr in `cd ${rundir}/${subdir}; ls -d *`;do 
           if [ ! -d ${outdir}/${subdir}/${xdr} ]; then 
              mkdir -p ${outdir}/${subdir}/${xdr} 
           fi 
           if [ -d ${rundir}/${subdir}/${xdr} ]; then
             mv ${rundir}/${subdir}/${xdr}/* ${outdir}/${subdir}/${xdr}/
           fi
         done 
      fi
    fi 
  done 
  echo $rundir
  i=$((i+1))
done 

