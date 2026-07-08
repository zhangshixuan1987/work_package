#!/bin/bash -el
#===========================================
# Batch system directives
#===========================================
#SBATCH --account=esmd
#SBATCH --time=02:00:00
#SBATCH --partition=short
#SBATCH --job-name=e3sm_compress_diag
#SBATCH --nodes=1
#SBATCH --output=e3sm_compress_diag.%j
#SBATCH --exclusive
#SBATCH --no-kill
#SBATCH --requeue

IFS=$'\n\t'

#===========================================
# Environment setup
#===========================================
echo "== Start of output file compressing =="
date
echo "============================================"

source /qfs/people/zhan391/e3sm_dart_work/code/DART/models/eam-se/work/env_mach_specific.sh

source ./create_and_setup_case.sh

#===========================================
# Paths
#===========================================
E3SM_ROOT=${my_e3sm_code}
DART_ROOT=${my_dart_code}
DART_MODEL=${my_dart_eam}
DART_SCPTDIR=${DART_ROOT}/models/${DART_MODEL}/shell_scripts
DART_WORKDIR=${DART_ROOT}/models/${DART_MODEL}/work
BASE_OBSDIR=${my_dart_obsdir}
BASE_PHIS=${my_e3sm_topo}
BASE_SEMAPS=${my_e3sm_semap}
BASE_CSGRID=${my_e3sm_csgrid}
CASE_ROOT=${my_modeldir}/case_scripts
RUN_ROOT=${my_modeldir}/run
ARCHIVE_DIR="${my_modeldir}/archive"

if [[ "${my_runtype}" == "AMIP" ]]; then
  comps=(rest atm cpl lnd rof)
else
  comps=(rest atm cpl ice lnd ocn rof)
fi

MAX_PARALLEL=40

compress_nc_files() {
  local workdir=$1
  if [ ! -d "$workdir/tmp_dir" ]; then
    mkdir -p "$workdir/tmp_dir"
    mv "$workdir"/*${timstr}*.nc "$workdir/tmp_dir/"
  fi

  local i=0
  for file in "$workdir/tmp_dir"/*.nc; do
    if [ -f "$file" ] && [ -s "$file" ]; then
      local fname=$(basename "$file")
      nccopy -7 -d 5 "$file" "$workdir/$fname" &
      i=$((i+1))
    fi
    if (( i == MAX_PARALLEL )); then
      wait
      i=0
    fi
  done
  wait

  local all_complete=0
  for file in "$workdir/tmp_dir"/*.nc; do
    local fname=$(basename "$file")
    if [ ! -f "$workdir/$fname" ]; then
      all_complete=1
    fi
  done

  if [[ $all_complete == 0 ]]; then
    rm -rvf "$workdir/tmp_dir"
  fi
}

cd ${my_wkdir}

DATA_ASSIMILATION_WINDOW=${my_dart_window}

ZIP_DATE=$1
ZIP_TOD=$2

DART_YEAR=${ZIP_DATE:0:4}
DART_MONTH=${ZIP_DATE:4:2}
DART_DAY=${ZIP_DATE:6:2}
DART_SECONDS=$ZIP_TOD
DART_HOUR=`echo "${ZIP_TOD} / 3600" | bc`
echo "valid time of current compress cycle is $DART_YEAR $DART_MONTH $DART_DAY $DART_SECONDS (seconds)"

DART_DATE=$(date -d "${DART_YEAR}-${DART_MONTH}-${DART_DAY} ${DART_HOUR}:00 - ${DATA_ASSIMILATION_WINDOW} hours" +"%Y-%m-%d %H")

# Extract YYYYMMDD from DART_DATE
DART_DATE_YMD=$(date -d "$DART_DATE" +"%Y%m%d")
# Construct final timestamp string
TIMESTAMP_STRING="${DART_DATE_YMD}-${DART_SECONDS}"

for rundir in "${comps[@]}"; do
  if [[ $rundir != 'rest' ]]; then
    workdir="${ARCHIVE_DIR}/${rundir}/hist"
    if [ -d "$workdir" ]; then
      echo "$workdir"
      compress_nc_files "$workdir"
    fi
  else
    workdir="${ARCHIVE_DIR}/${rundir}/${TIMESTAMP_STRING}"
    if [ -d "$workdir" ]; then
      echo "$workdir"
      compress_nc_files "$workdir"
    fi
  fi
done

wait
echo "===== End of output file compressing =="
date
echo "====================================="
