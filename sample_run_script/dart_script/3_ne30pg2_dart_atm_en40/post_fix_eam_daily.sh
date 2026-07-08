#!/bin/bash -el
#------------------------------------------------------------------------------
# SLURM Batch Directives
#------------------------------------------------------------------------------
#SBATCH --account=esmd
#SBATCH --time=2:00:00
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
var1_list=(TCO SCO FLUT PRECT PS TREFHT TREFHTMX TREFHTMN QREFHT U200 U850 V200 V850) 

hist1="eam.h1"
freq1="daily"
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
      if [ "$year" -eq "$ey" ] && [ "$month" -gt "$em" ]; then continue; fi

      yymm=$(printf "%04d-%02d" ${year} ${month})

      # Link h1 files
      for ff in "${input}/${CASE_NAME}.${hist1}.${yymm}"*.nc; do
        [ -f "${ff}" ] && ln -sf "${ff}" .
      done

      nday=${mday[$((month - 1))]}
      # Correct leap year check: only adjust February (month == 2)
      if (( month == 2 )) && (( (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0) )); then
          nday=29
      fi
      if [[ ${nday} -gt ${ed} ]];then 
        nday=${ed}
      fi 

      # === Step 2: Regrid var1_list from hist1 (daily) ===
      ts_dest1="${outdir}/atm/180x360_aave/ts/${freq1}"
      mkdir -p "${ts_dest1}"

      for var in "${var1_list[@]}"; do
         xfiles="${outdir}_back/atm/180x360_aave/ts/${freq1}/${var}.${ENSTR}.*.nc"
         yyyy=$(printf "%04d" ${year})
         mm=$(printf "%02d" ${month})
         outfile="${var}.${ENSTR}.${yyyy}.nc"
         [ -f "${outfile}" ] && rm -vf "${outfile}"
         for day in `seq 1 ${nday}`;do 
           dd=$(printf "%02d" ${day})
           start_time="${yyyy}-${mm}-${dd} 00:00:0.0"
           end_time="${yyyy}-${mm}-${dd} 23:59:59.0"
           echo $start_time $end_time
           ncra -d time,"${start_time}","${end_time}" ${xfiles[@]} ftmp_`printf "%02d" ${day}`.nc
         done 
         ncrcat -O -d time,0, ftmp_??.nc "${ts_dest1}/${outfile}" 
         rm -rvf ftmp_??.nc
      done
    done
  done

  cd ..
  rm -rf "${workdir}"
done

echo "===== End of DART diagnostic ====="
date
echo "==================================="
