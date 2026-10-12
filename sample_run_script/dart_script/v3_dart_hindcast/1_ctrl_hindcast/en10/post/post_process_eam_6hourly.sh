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
source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh

# System utilities
MOVE='/usr/bin/mv'
COPY='/usr/bin/cp --preserve=timestamps'
LINK='/usr/bin/ln -fs'
REMOVE='/usr/bin/rm'
LIST='/usr/bin/ls'

# Environment setup (assumes these are exported externally or in create_and_setup_case.sh)
# E3SM_ROOT, DART_ROOT, my_modeldir, my_ensnum, my_casename, etc.

DART_MODEL=${my_dart_eam}
DART_WORKDIR=${DART_ROOT}/models/${DART_MODEL}/work
ARCHIVE_DIR="${my_modeldir}/archive"
MAP_FILE="/compyfs/zhan391/acme_init/map_file/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc"

# Dates
ymds=""
ymde=""

read -r sy sm sd <<< "$(echo ${ymds} | tr '-' ' ')"
read -r ey em ed <<< "$(echo ${ymde} | tr '-' ' ')"
mday=(31 28 31 30 31 30 31 31 30 31 30 31)

var2_list=(FLUT FLUTC LHFLX SHFLX PRECT PRECC PRECL PS PSL QFLX QREFHT TMQ TS TREFHT TUQ TVQ OMEGA500 PRECSL SWCF LWCF
           TTOP TAUX TAUY TGCLDCWP TGCLDIWP TGCLDLWP U90M V90M U10 CLDTOT CLDLOW CLDMED FLDS FLDSC FLNS FLNSC
           FLNT FLNTC FSDS FSDSC FSNS FSNSC FSNT FSNTC FSNTOA FSNTOAC FSUTOA)

hist="eam.h2"
freq="6hourly"

input="${ARCHIVE_DIR}/atm/hist"
outdir="${ARCHIVE_DIR}/post"

mkdir -p "${outdir}"

jobid=${SLURM_JOBID}

cd ${outdir}
workdir=$(mktemp -d tmp.${jobid}.XXXX)
cd ${workdir}

ENSTR="EN01" 
CASE_NAME="${my_casename}.${ENSTR}"

flist=$(mktemp ./input${ENSTR}.XXXXXX)

echo "=== Starting ensemble member ${CASE_NAME}.${hist}.${freq} ==="

for year in $(seq ${sy} ${ey}); do
  for month in $(seq 1 12); do
    # Skip months outside the desired range if first/last year
    if [ "$year" -eq "$sy" ] && [ "$month" -lt "$sm" ]; then continue; fi
    if [ "$year" -eq "$ey" ] && [ "$month" -ge "$em" ]; then continue; fi

    yymm=$(printf "%04d-%02d" ${year} ${month})

    # Link h2 files
    for ff in "${input}/${CASE_NAME}.${hist}.${yymm}"*.nc; do
      [ -f "${ff}" ] && ln -sf "${ff}" .
    done
  done
done

ls ${CASE_NAME}.${hist}.????-??-*.nc > ${flist}

# === Step 3: Regrid var2_list from hist (6-hourly) ===
ts_dest="${outdir}/atm/180x360_aave/ts/${freq}"
mkdir -p "${ts_dest}"

mapfile -t ffiles2 < "${flist}"
for year in $(seq "${sy}" "${ey}"); do
  yyyy=$(printf "%04d" "${year}")
  start_time="${yyyy}-01-01 00:00:00.0"   # fix seconds
  end_time="${yyyy}-12-31 23:59:59.0"
  for var in "${var2_list[@]}"; do
    outfile="${var}.${ENSTR}.${yyyy}.nc"
    outpath="${ts_dest}/${outfile}"

    if [ "${#ffiles2[@]}" -gt 0 ]; then
      # Find a file that actually has ${var}; try first, then fall back
      probe_file=""
      if ncks -m -v "${var}" "${ffiles2[0]}" >/dev/null 2>&1; then
        probe_file="${ffiles2[0]}"
      else
        for f in "${ffiles2[@]}"; do
          if ncks -m -v "${var}" "$f" >/dev/null 2>&1; then
            probe_file="$f"; break
          fi
        done
      fi

      if [ -n "${probe_file}" ]; then
        echo "Extracting ${var} for ${ENSTR}, year ${yyyy}..."
        ncrcat -O -d time,"${start_time}","${end_time}" -v "${var}" \
          "${ffiles2[@]}" "${outfile}" && \
        ncremap -m "${MAP_FILE}" -i "${outfile}" -o "${outpath}"
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
