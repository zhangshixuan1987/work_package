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
MAP_FILE="/compyfs/zhan391/acme_init/map_file/map_r05_to_cmip6_180x360_aave.20200901.nc"

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

hist="elm.h1"
freq="daily"
input="${ARCHIVE_DIR}/lnd/hist"
outdir="${ARCHIVE_DIR}/post"

mkdir -p "${outdir}"

jobid=${SLURM_JOBID}

cd ${outdir}
workdir=$(mktemp -d tmp.${jobid}.XXXX)
cd ${workdir}

ENSTR=$(printf "EN%02d" ${i})
CASE_NAME="${my_casename}.${ENSTR}"

flist=$(mktemp ./input${ENSTR}.XXXXXX)

echo "=== Starting ensemble member ${CASE_NAME}.${hist}.${freq} ==="

for year in $(seq ${sy} ${ey}); do
  for month in $(seq 1 12); do
    # Skip months outside the desired range if first/last year
    if [ "$year" -eq "$sy" ] && [ "$month" -lt "$sm" ]; then continue; fi
    if [ "$year" -eq "$ey" ] && [ "$month" -ge "$em" ]; then continue; fi

    yymm=$(printf "%04d-%02d" ${year} ${month})

    # Link h1 files
    for ff in "${input}/${CASE_NAME}.${hist}.${yymm}"*.nc; do
      [ -f "${ff}" ] && ln -sf "${ff}" .
    done

  done
done

ls ${CASE_NAME}.${hist}.????-??-*.nc > ${flist}

# === Step 2: Regrid var1_list from hist (daily) ===
ts_dest="${outdir}/lnd/180x360_aave/ts/${freq}"
mkdir -p "${ts_dest}"

mapfile -t ffiles1 < "${flist}"
for year in $(seq ${sy} ${ey}); do
  yyyy=$(printf "%04d" ${year})
  start_time="${yyyy}-01-01 00:00:00.0"   # TWEAK: fix seconds field
  end_time="${yyyy}-12-31 23:59:59.0"
  for var in "${var1_list[@]}"; do
    outfile="${var}.${ENSTR}.${yyyy}.nc"
    outpath="${ts_dest}/${outfile}"

    if [ "${#ffiles1[@]}" -gt 0 ]; then
      # TWEAK: find a file that actually has ${var}; fall back to first match
      probe_file=""
      if ncks -m -v "${var}" "${ffiles1[0]}" >/dev/null 2>&1; then
        probe_file="${ffiles1[0]}"
      else
        for f in "${ffiles1[@]}"; do
          if ncks -m -v "${var}" "$f" >/dev/null 2>&1; then
            probe_file="$f"; break
          fi
        done
      fi

      if [ -n "${probe_file}" ]; then
        echo "Extracting ${var} for ${ENSTR}, year ${yyyy}..."
        ncrcat -O -d time,"${start_time}","${end_time}" -v "${var},landfrac" "${ffiles1[@]}" "${outfile}" && \
        ncremap -P elm -m "${MAP_FILE}" -i "${outfile}" -o "${outpath}"
      else
        echo "Warning: variable ${var} not found in listed files for ${ENSTR}, year ${yyyy}"
      fi
    else
      echo "Warning: No hist files found for ${ENSTR}, year ${yyyy}, var ${var}"
    fi
  done
done

cd ..
rm -rf "${workdir}"

echo "===== End of DART diagnostic ====="
date
echo "==================================="
