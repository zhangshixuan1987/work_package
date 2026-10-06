#!/bin/bash
ulimit -S -c 0
set -euo pipefail

data_path="/compyfs/zhan391/acme_init/Observations"
obsname="CERES-OAFlux"                        # adjust to your file naming
outdir="/compyfs/zhan391/v3_dart_cda_scratch/obs_data/monthly"
mkdir -p "${outdir}"

for year in $(seq 2011 2012); do
  infile="${data_path}/${obsname}/monthly/${obsname}_${year}.nc"   
  for mon in $(seq 1 12); do
    mm=$(printf "%02d" "${mon}")
    start="${year}-${mm}-01"
    # last day of month: start +1 month -1 day
    end=$(date -u -d "${start} +1 month -1 day" +%Y-%m-%d)

    outfile="${obsname}.${year}.${mm}.nc"
    echo "Extracting ${start}..${end} -> ${outfile}"

    ncks -O -d time,"${start}","${end}" "${infile}" "${outdir}/${outfile}"
  done
done
