#!/bin/bash -el
#------------------------------------------------------------------------------
# SLURM Batch Directives
#------------------------------------------------------------------------------
#SBATCH --account=esmd
#SBATCH --time=8:00:00
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

# Optional overrides for reprocessing another run/period or a member subset:
#   POST_CASENAME, POST_YMDS, POST_YMDE, POST_MEM_BEG, POST_MEM_END
if [ -n "${POST_CASENAME}" ]; then
  my_casename=${POST_CASENAME}
  my_modeldir="${my_runpath}/${my_casename}"
fi

# Environment setup (assumes these are exported externally or in create_and_setup_case.sh)
# E3SM_ROOT, DART_ROOT, my_modeldir, my_ensnum, my_casename, etc.

DART_MODEL=${my_dart_eam}
DART_WORKDIR=${DART_ROOT}/models/${DART_MODEL}/work
ARCHIVE_DIR="${my_modeldir}/archive"
MAP_FILE=${my_e3sm_rgdmap}

# Dates
ymds=${POST_YMDS:-"2011-12-01"}
ymde=${POST_YMDE:-"2011-12-31"}

read -r sy sm sd <<< "$(echo ${ymds} | tr '-' ' ')"
read -r ey em ed <<< "$(echo ${ymde} | tr '-' ' ')"
mday=(31 28 31 30 31 30 31 31 30 31 30 31)

# Variables
var2_list=(FLUT FLUTC LHFLX SHFLX PRECT PRECC PRECL PS PSL QFLX QREFHT TMQ TS TREFHT TUQ TVQ OMEGA500 PRECSL SWCF LWCF
           TTOP TAUX TAUY TGCLDCWP TGCLDIWP TGCLDLWP U90M V90M U10 CLDTOT CLDLOW CLDMED FLDS FLDSC FLNS FLNSC
           FLNT FLNTC FSDS FSDSC FSNS FSNSC FSNT FSNTC FSNTOA FSNTOAC FSUTOA)

hist2="eam.h2"
freq2="6hourly"
jobid=${SLURM_JOBID}

for i in $(seq ${POST_MEM_BEG:-1} ${POST_MEM_END:-${my_ensnum}}); do
  ENSTR=$(printf "EN%02d" ${i})
  CASE_NAME="${my_casename}.${ENSTR}"

  # Archives are stored per member: <case>/ENxx/archive/{atm/hist,post}
  input="${my_modeldir}/${ENSTR}/archive/atm/hist"
  outdir="${my_modeldir}/${ENSTR}/archive/post"
  if [ ! -d "${input}" ]; then
    echo "Error: ${input} not found"; exit 1
  fi
  mkdir -p "${outdir}"
  cd ${outdir}
  workdir=$(mktemp -d tmp.${jobid}.XXXX)
  cd ${workdir}

  echo "=== Starting ensemble member ${ENSTR} ==="

  for year in $(seq ${sy} ${ey}); do
    for month in $(seq 1 12); do
      # Skip months outside the desired range if first/last year
      if [ "$year" -eq "$sy" ] && [ "$month" -lt "$sm" ]; then continue; fi
      if [ "$year" -eq "$ey" ] && [ "$month" -gt "$em" ]; then continue; fi

      yymm=$(printf "%04d-%02d" ${year} ${month})

      # Link h2 files
      for ff in "${input}/${CASE_NAME}.${hist2}.${yymm}"-??-?????.nc; do
        [ -f "${ff}" ] && ln -sf "${ff}" .
      done
    done
  done

  # === Step 3: Regrid var2_list from hist2 (6-hourly) ===
  ts_dest1="${outdir}/atm/180x360_aave/ts/${freq2}_da"
  ts_dest2="${outdir}/atm/180x360_aave/ts/${freq2}"
  mkdir -p "${ts_dest1}"
  mkdir -p "${ts_dest2}"

  ffiles2=${CASE_NAME}.${hist2}.*.nc
  echo ${ffiles2[@]}
  for year in $(seq ${sy} ${ey}); do
    yyyy=$(printf "%04d" ${year})
    for var in "${var2_list[@]}"; do
      outfile="${var}.${ENSTR}.${yyyy}.nc"
      [ -f "${outfile}" ] && rm -vf "${outfile}"
      # Each restart-cycle h2 file holds 2 records: [0] the zero-width value at
      # the cycle start (analysis time) and [1] the 6-h segment, stamped at its
      # end (e.g. the Dec 31 18Z mean is stamped Jan 1 00Z). Select cycles by
      # file name, not by time stamp, so every cycle of the year is kept.
      year_files=(${CASE_NAME}.${hist2}.${yyyy}-??-??-?????.nc)  # cycle files only, no stray *.tmp.*.nc
      if [ -f "${year_files[0]}" ]; then
        ncrcat -O -d time,0,,2 -v "${var}" ${year_files[@]} ${var}_tmp_da.nc
        ncremap -m "${MAP_FILE}" -i ${var}_tmp_da.nc -o "${ts_dest1}/${outfile}"
        ncrcat -O -d time,1,,2 -v "${var}" ${year_files[@]} ${var}_tmp_seg.nc
        # Stamp 6-h means at their interval midpoint so they fall in the right
        # day/month; instantaneous (time: point) fields keep their valid time.
        if ncks -m -v "${var}" ${var}_tmp_seg.nc | grep -q 'cell_methods = "time: mean"'; then
          ncap2 -O -s 'time=0.5*(time_bnds(:,0)+time_bnds(:,1))' ${var}_tmp_seg.nc ${var}_tmp_seg.nc
        fi
        ncremap -m "${MAP_FILE}" -i ${var}_tmp_seg.nc -o "${ts_dest2}/${outfile}"
        rm -f ${var}_tmp_da.nc ${var}_tmp_seg.nc
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
