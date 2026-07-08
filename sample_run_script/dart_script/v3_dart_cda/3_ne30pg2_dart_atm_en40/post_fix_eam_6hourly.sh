#!/bin/bash -el
#------------------------------------------------------------------------------
# SLURM Batch Directives
#------------------------------------------------------------------------------
#SBATCH --account=esmd
#SBATCH --time=12:00:00
#SBATCH --partition=slurm
#SBATCH --job-name=regrid_diag
#SBATCH --nodes=1
#SBATCH --output=regrid_diag.%j
#SBATCH --exclusive
#SBATCH --no-kill
#SBATCH --requeue

echo "== Start of DART diagnostic =="
date
echo "============================================"

# Load conda environment
source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh

# System utilities
MOVE='/usr/bin/mv'
COPY='/usr/bin/cp --preserve=timestamps'
LINK='/usr/bin/ln -fs'
REMOVE='/usr/bin/rm'
LIST='/usr/bin/ls'

# Launch info
my_wkdir=${PWD}
cd ${my_wkdir}
source ./create_and_setup_case.sh

# Environment setup (assumes these are exported externally or in create_and_setup_case.sh)
# E3SM_ROOT, DART_ROOT, my_modeldir, my_ensnum, my_casename, etc.

DART_MODEL=${my_dart_eam}
DART_WORKDIR=${DART_ROOT}/models/${DART_MODEL}/work
ARCHIVE_DIR="${my_modeldir}/archive"
MAP_FILE=${my_e3sm_rgdmap}

# Dates
ymds="2011-12-01"
ymde="2011-12-31"

read -r sy sm sd <<< "$(echo ${ymds} | tr '-' ' ')"
read -r ey em ed <<< "$(echo ${ymde} | tr '-' ' ')"
mday=(31 28 31 30 31 30 31 31 30 31 30 31)

# Variables
var2_list=(FLUT FLUTC LHFLX SHFLX PRECT PRECC PRECL PS PSL QFLX QREFHT TMQ TS TREFHT TUQ TVQ OMEGA500 PRECSL SWCF LWCF
           TTOP TAUX TAUY TGCLDCWP TGCLDIWP TGCLDLWP U90M V90M U10 CLDTOT CLDLOW CLDMED FLDS FLDSC FLNS FLNSC
           FLNT FLNTC FSDS FSDSC FSNSC FSNT FSNTC FSNTOA FSNTOAC FSUTOA)

hist2="eam.h2"
freq2="6hourly"
input="${ARCHIVE_DIR}/atm/hist"
outdir="${ARCHIVE_DIR}/post"

mkdir -p "${outdir}"

jobid=${SLURM_JOBID}

for i in $(seq 1 ${my_ensnum}); do
  cd ${outdir}
  workdir=$(mktemp -d tmp.${jobid}.XXXX)
  cd ${workdir}

  ENSTR=$(printf "EN%02d" ${i})
  CASE_NAME="${my_casename}.${ENSTR}"

  echo "=== Starting ensemble member ${ENSTR} ==="

  # === Step 3: Regrid var2_list from hist2 (6-hourly) ===
  ts_dest1="${outdir}/atm/180x360_aave/ts/${freq2}_da"
  ts_dest2="${outdir}/atm/180x360_aave/ts/${freq2}"
  mkdir -p "${ts_dest1}"
  mkdir -p "${ts_dest2}"

  ffiles2="${outdir}_back/atm/180x360_aave/ts/${freq2}/*${ENSTR}*.nc"
  for file in ${ffiles2[@]};do 
    outfile=`basename $file`
      if [ -f "${file}" ]; then
        ncks -O -d time,0,,2  ${file} "${ts_dest1}/${outfile}"
        ncks -O -d time,1,,2  ${file}  ${outfile}_time1.nc
        ncks -O -d time,0     ${file}  ${outfile}_time0.nc
        ncrcat -O -d time,0,  ${outfile}_time0.nc ${outfile}_time1.nc "${ts_dest2}/${outfile}"
      else
        echo "Warning: No hist2 files found for ${ENSTR}, year ${yyyy}, ${file}"
      fi
  done
  
  cd ..
  rm -rf "${workdir}"
done

echo "===== End of DART diagnostic ====="
date
echo "==================================="
