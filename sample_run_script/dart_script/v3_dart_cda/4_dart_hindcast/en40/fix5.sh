#!/bin/bash -el
#------------------------------------------------------------------------------
# SLURM Batch Directives
#------------------------------------------------------------------------------
#SBATCH --account=esmd
#SBATCH --time=02:00:00
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
MAP_FILE=${my_elm_rgdmap}

# Dates
ymds="2012-01-01"
ymde="2012-03-01"

read -r sy sm sd <<< "$(echo ${ymds} | tr '-' ' ')"
read -r ey em ed <<< "$(echo ${ymde} | tr '-' ' ')"
mday=(31 28 31 30 31 30 31 31 30 31 30 31)

# Variables
var1_list=(
  H2OSNO FSNO QRUNOFF QSNOMELT FSNO_EFF SNORDSL SNOW FSA FSDS FSR FLDS
  FIRE FIRA SOILWATER_10CM SOILLIQ SOILICE QSOIL U10 U10WITHGUSTS TWS
  TSOI_10CM TSA TLAI THBOT TSOI TSOI_ICE SOILLIQ_ICE SOILICE_ICE
  TAUX TAUY FSH HC HCSOI EFLX_LH_TOT SNOW_DEPTH RH2M RAIN QVEGE QVEGT
  QBOT Q2M H2OSOI H2OSFC ZWT ZBOT TBOT TG THBOT PBOT
)

var1_list=(H2OSNO FSNO QRUNOFF QSNOMELT FSNO_EFF SNORDSL SNOW )


hist1="elm.h1"
freq1="daily"
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
      if [ "$year" -eq "$ey" ] && [ "$month" -ge "$em" ]; then continue; fi

      yymm=$(printf "%04d-%02d" ${year} ${month})

      # Link h1 files
      for ff in "${input}/${CASE_NAME}.${hist1}.${yymm}"*.nc; do
        [ -f "${ff}" ] && ln -sf "${ff}" .
      done

    done
  done

  ls ${CASE_NAME}.${hist1}.????-??-*.nc > input1.txt

  # === Step 2: Regrid var1_list from hist1 (daily) ===
  ts_dest1="${outdir}/lnd/180x360_aave/ts/${freq1}"
  mkdir -p "${ts_dest1}"

  mapfile -t ffiles1 < input1.txt
  for year in $(seq ${sy} ${ey}); do
    yyyy=$(printf "%04d" ${year})
    start_time="${yyyy}-01-01 00:00:0.0"
    end_time="${yyyy}-12-31 23:59:59.0"
    for var in "${var1_list[@]}"; do
      outfile="${var}.${ENSTR}.${yyyy}.nc"
      [ -f "${outfile}" ] && rm -vf "${outfile}"
      if [ "${#ffiles1[@]}" -gt 0 ]; then
        ncrcat -O -d time,"${start_time}","${end_time}" -v "${var},landfrac" "${ffiles1[@]}" "${outfile}" && \
        ncremap -P elm -m "${MAP_FILE}" -i "${outfile}" -o "${ts_dest1}/${outfile}"
      else
        echo "Warning: No hist1 files found for ${ENSTR}, year ${yyyy}, var ${var}"
      fi
    done
  done

  cd ..
  rm -rf "${workdir}"
done

echo "===== End of DART diagnostic ====="
date
echo "==================================="
