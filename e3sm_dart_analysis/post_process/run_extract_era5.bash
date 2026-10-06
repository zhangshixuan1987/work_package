#!/bin/bash
ulimit -S -c 0

data_path="/compyfs/zhan391/acme_init/Observations"
obsname="ERA5"
yymmdds="2011-12-01 00:00:0.0"
yymmdde="2011-12-31 23:59:59.0"
yymm="2011-12"

for var in PRECT PSL T850 Q850 U850 V850 Z500;do 
  freq="monthly"
  outdir="clim"
  infile="${data_path}/${obsname}.6hourly/*.${var}.*2011*-2011*.nc"
  outfile="${var}.${obsname}.${freq}.${yymm}.nc"
  if [ ! -d "${outdir}" ];then 
    mkdir -p ${outdir}
  fi 
  if [ -f "${outdir}/${outfile}" ];then
    rm -rvf ${outdir}/${outfile}
  fi
  ncra -d time,"${yymmdds}","${yymmdde}" ${infile} ${outdir}/${outfile}

  freq="6hourly"
  outdir="ts"
  infile="${data_path}/${obsname}.6hourly/*.${var}.*2011*-2011*.nc"
  outfile="${var}.${obsname}.${freq}.${yymm}.nc"
  if [ ! -d "${outdir}" ];then
    mkdir -p ${outdir}
  fi
  if [ -f "${outdir}/${outfile}" ];then 
    rm -rvf ${outdir}/${outfile}
  fi 
  ncks -d time,"${yymmdds}","${yymmdde}" ${infile} ${outdir}/${outfile}
done
