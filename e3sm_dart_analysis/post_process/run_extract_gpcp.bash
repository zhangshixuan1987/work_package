#!/bin/bash
ulimit -S -c 0

data_path="/compyfs/zhan391/acme_init/Observations"
obsname="GPCP"
yymmdds="2011-12-01 00:00:0.0"
yymmdde="2011-12-31 23:59:59.0"
yymm="2011-12"

freq="monthly"
outdir="clim"
infile="${data_path}/${obsname}/${freq}/2009-2011.nc"
outfile="PRECT.${obsname}.${freq}.${yymm}.nc"
if [ ! -d "${outdir}" ];then 
  mkdir -p ${outdir}
fi 
if [ -f "${outdir}/${outfile}" ];then
  rm -rvf ${outdir}/${outfile}
fi
ncks -d time,"${yymmdds}","${yymmdde}" ${infile} ${outdir}/${outfile}

freq="daily"
outdir="ts"
infile="${data_path}/${obsname}/${freq}/2009-2011.nc"
outfile="PRECT.${obsname}.${freq}.${yymm}.nc"
if [ ! -d "${outdir}" ];then
  mkdir -p ${outdir}
fi
if [ -f "${outdir}/${outfile}" ];then 
  rm -rvf ${outdir}/${outfile}
fi 
ncks -d time,"${yymmdds}","${yymmdde}" ${infile} ${outdir}/${outfile}


