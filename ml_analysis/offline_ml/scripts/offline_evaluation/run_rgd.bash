#!/bin/bash -l

source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

#conda activate e3sm_analysis

infile="$1"
outfile="$2"
map="/global/cfs/cdirs/e3sm/diagnostics/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc"
ncremap -m ${map} -i ${infile} -o ${outfile}

