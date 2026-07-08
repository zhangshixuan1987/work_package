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
MAP_FILE=${my_e3sm_rgdmap}

# Dates
ymds="2012-01-01"
ymde="2012-03-01"

read -r sy sm sd <<< "$(echo ${ymds} | tr '-' ' ')"
read -r ey em ed <<< "$(echo ${ymde} | tr '-' ' ')"
mday=(31 28 31 30 31 30 31 31 30 31 30 31)

# Variables
var1_list=(TCO SCO FLUT PRECT PS PSL TS TBOT TUQ TVQ UBOT VBOT ZBOT U10 TAUX TAUY TREFHT QREFHT TMQ CAPE CIN TTQ
           U1000 U975 U950 U925 U900 U850 U800 U700 U600 U500 U400 U300 U200 U100 U050 U010
           V1000 V975 V950 V925 V900 V850 V800 V700 V600 V500 V400 V300 V200 V100 V050 V010
           Z1000 Z975 Z950 Z925 Z900 Z850 Z800 Z700 Z600 Z500 Z400 Z300 Z200 Z100 Z050 Z010
           T1000 T975 T950 T925 T900 T850 T800 T700 T600 T500 T400 T300 T200 T100 T050 T010
           Q1000 Q975 Q950 Q925 Q900 Q850 Q800 Q700 Q600 Q500 Q400 Q300 Q200 Q100 Q050 Q010
           OMEGA1000 OMEGA975 OMEGA950 OMEGA925 OMEGA900 OMEGA850 OMEGA800 OMEGA700 OMEGA600 OMEGA500
           OMEGA400 OMEGA300 OMEGA200 OMEGA100 TREFHTMN TREFHTMX PRECTMX)

var2_list=(FLUT FLUTC LHFLX SHFLX PRECT PRECC PRECL PS PSL QFLX QREFHT TMQ TS TREFHT TUQ TVQ OMEGA500 PRECSL SWCF LWCF
           TTOP TAUX TAUY TGCLDCWP TGCLDIWP TGCLDLWP U90M V90M U10 CLDTOT CLDLOW CLDMED FLDS FLDSC FLNS FLNSC
           FLNT FLNTC FSDS FSDSC FSNSC FSNT FSNTC FSNTOA FSNTOAC FSUTOA)

hist0="eam.h0"
hist1="eam.h1"
hist2="eam.h2"
freq1="daily"
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

  for year in $(seq ${sy} ${ey}); do
    for month in $(seq 1 12); do
      # Skip months outside the desired range if first/last year
      if [ "$year" -eq "$sy" ] && [ "$month" -lt "$sm" ]; then continue; fi
      if [ "$year" -eq "$ey" ] && [ "$month" -ge "$em" ]; then continue; fi

      yymm=$(printf "%04d-%02d" ${year} ${month})

      # Link h0 file
      h0file="${input}/${CASE_NAME}.${hist0}.${yymm}.nc"
      [ -f "${h0file}" ] && ln -sf "${h0file}" .

      # Link h1 files
      for ff in "${input}/${CASE_NAME}.${hist1}.${yymm}"*.nc; do
        [ -f "${ff}" ] && ln -sf "${ff}" .
      done

      # Link h2 files
      for ff in "${input}/${CASE_NAME}.${hist2}.${yymm}"*.nc; do
        [ -f "${ff}" ] && ln -sf "${ff}" .
      done
    done
  done

  ls ${CASE_NAME}.${hist0}.????-??.nc > input0.txt
  ls ${CASE_NAME}.${hist1}.????-??-*.nc > input1.txt
  ls ${CASE_NAME}.${hist2}.????-??-*.nc > input2.txt

  # === Step 1: Regrid monthly h0 climatology files ===
  clim_dest="${outdir}/atm/180x360_aave/clim"
  mkdir -p "${clim_dest}"
  while IFS= read -r ff; do
    outfile=$(basename "${ff}")
    ncremap -m "${MAP_FILE}" -i "${ff}" -o "${clim_dest}/${outfile}"
  done < input0.txt

  # === Step 2: Regrid var1_list from hist1 (daily) ===
  ts_dest1="${outdir}/atm/180x360_aave/ts/${freq1}"
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
        ncrcat -O -d time,"${start_time}","${end_time}" -v "${var}" "${ffiles1[@]}" "${outfile}" && \
        ncremap -m "${MAP_FILE}" -i "${outfile}" -o "${ts_dest1}/${outfile}"
      else
        echo "Warning: No hist1 files found for ${ENSTR}, year ${yyyy}, var ${var}"
      fi
    done
  done

  # === Step 3: Regrid var2_list from hist2 (6-hourly) ===
  ts_dest2="${outdir}/atm/180x360_aave/ts/${freq2}"
  mkdir -p "${ts_dest2}"

  mapfile -t ffiles2 < input2.txt
  for year in $(seq ${sy} ${ey}); do
    yyyy=$(printf "%04d" ${year})
    start_time="${yyyy}-01-01 00:00:0.0"
    end_time="${yyyy}-12-31 23:59:59.0"
    for var in "${var2_list[@]}"; do
      outfile="${var}.${ENSTR}.${yyyy}.nc"
      [ -f "${outfile}" ] && rm -vf "${outfile}"
      if [ "${#ffiles2[@]}" -gt 0 ]; then
        ncrcat -O -d time,"${start_time}","${end_time}" -v "${var}" "${ffiles2[@]}" "${outfile}" && \
        ncremap -m "${MAP_FILE}" -i "${outfile}" -o "${ts_dest2}/${outfile}"
      else
        echo "Warning: No hist2 files found for ${ENSTR}, year ${yyyy}, var ${var}"
      fi
    done
  done

  cd ..
  rm -rf "${workdir}"
done

echo "===== End of DART diagnostic ====="
date
echo "==================================="
