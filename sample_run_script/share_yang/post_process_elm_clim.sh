#!/bin/bash -el
#------------------------------------------------------------------------------
# SLURM Batch Directives
#------------------------------------------------------------------------------
#SBATCH --account=esmd
#SBATCH --time=2:00:00
#SBATCH --partition=short
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
source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

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
MAP_FILE=${my_elm_rgdmap}

# Dates
ymds="2011-12-01"
ymde="2011-12-31"

read -r sy sm sd <<< "$(echo ${ymds} | tr '-' ' ')"
read -r ey em ed <<< "$(echo ${ymde} | tr '-' ' ')"
mday=(31 28 31 30 31 30 31 31 30 31 30 31)

hist0="elm.h0"
input="${ARCHIVE_DIR}/lnd/hist"
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

  for year in $(seq ${sy} ${ey}); do
    for month in $(seq 1 12); do
      # Skip months outside the desired range if first/last year
      if [ "$year" -eq "$sy" ] && [ "$month" -lt "$sm" ]; then continue; fi
      if [ "$year" -eq "$ey" ] && [ "$month" -gt "$em" ]; then continue; fi

      yymm=$(printf "%04d-%02d" ${year} ${month})

      # Link h0 file
      h0file="${input}/${CASE_NAME}.${hist0}.${yymm}.nc"
      [ -f "${h0file}" ] && ln -sf "${h0file}" .

    done
  done

  # === Step 1: Regrid monthly h0 climatology files ===
  ls ${CASE_NAME}.${hist0}.????-??.nc > input0.txt
  clim_dest="${outdir}/lnd/180x360_aave/clim"
  mkdir -p "${clim_dest}"
  while IFS= read -r ff; do
    outfile=$(basename "${ff}")
    ncremap -P elm -m "${MAP_FILE}" -i "${ff}" -O ${clim_dest} #-o "${clim_dest}/${outfile}"
  done < input0.txt

  cd ..
  rm -rf "${workdir}"
done

echo "===== End of DART diagnostic ====="
date
echo "==================================="
