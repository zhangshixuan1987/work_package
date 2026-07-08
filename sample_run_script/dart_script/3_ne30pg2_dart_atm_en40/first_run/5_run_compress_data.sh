#!/bin/bash -el

set -euo pipefail
IFS=$'\n\t'

#===========================================
# Environment setup
#===========================================
echo "== Start of output file compressing =="
date
echo "============================================"

source /qfs/people/zhan391/e3sm_dart_work/code/DART/models/eam-se/work/env_mach_specific.sh

VERBOSE='-v'
MOVE='/usr/bin/mv'
COPY='/usr/bin/cp --preserve=timestamps'
LINK='/usr/bin/ln -fs'
LINKV=TRUE
LIST='/usr/bin/ls'
REMOVE='/usr/bin/rm'
LAUNCHCMD=mpirun.lsf

my_wkdir=${PWD}
scomp="eam"
cd ${my_wkdir}

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

if [ ${my_runtype} == 'AMIP' ]; then
  comps=(rest atm cpl lnd rof)
else
  comps=(rest atm cpl ice lnd ocn rof)
fi

IFS='-' read -r -a START_DATE <<< "$1"
IFS='-' read -r -a END_DATE <<< "$2"

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

for year in $(seq ${START_DATE[0]} ${END_DATE[0]}); do
  for month in $(seq ${START_DATE[1]} ${END_DATE[1]}); do
    for day in $(seq ${START_DATE[2]} ${END_DATE[2]}); do
      timstr=$(printf "%04d-%02d-%02d" $year $month $day)
      echo "$timstr"
      for rundir in "${comps[@]}"; do
        if [[ $rundir != 'rest' ]]; then
          workdir="${ARCHIVE_DIR}/${rundir}/hist"
          if [ -d "$workdir" ]; then
            echo "$workdir"
            compress_nc_files "$workdir"
          fi
        else
          for tod in $(seq 0 21600 64800); do
            workdir="${ARCHIVE_DIR}/${rundir}/${timstr}-$(printf "%05d" $tod)"
            if [ -d "$workdir" ]; then
              echo "$workdir"
              compress_nc_files "$workdir"
            fi
          done
        fi
      done
    done
  done
done

wait
echo "===== End of output file compressing =="
date
echo "====================================="
