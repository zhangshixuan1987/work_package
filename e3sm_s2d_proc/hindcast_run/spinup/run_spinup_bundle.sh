#!/bin/bash -el
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --account=e3sm
#SBATCH  --constraint=cpu
#SBATCH  --qos=regular
#SBATCH  --nodes=7
#SBATCH  --time=8:00:00
#SBATCH  --job-name=e3sm_s2d_ensrun

source ./create_and_setup_case.sh

echo "============================================"
echo "  Run Short Spin-up Simulation"
echo "  Start time: $(date)"
echo "============================================"

if [[ "${do_short_spinup,,}" != "true" ]]; then
  exit 0
fi

do_e3sm_compile=false #true
do_ensemble_setup=true
do_ensemble_run=true

STOP_OPTION="nmonths"
STOP_N=1

REST_OPTION="nmonths"
REST_N=1
RESUBMIT=0

BASE_REPLACE_OCN=true 

declare -A YEAR_MAP=(
  [1980]=353 [1981]=354 [1982]=355 [1983]=356 [1984]=357
  [1985]=358 [1986]=359 [1987]=360 [1988]=361 [1989]=362
  [1990]=363 [1991]=364 [1992]=365 [1993]=366 [1994]=367
  [1995]=368 [1996]=369 [1997]=370 [1998]=371 [1999]=372
  [2000]=373 [2001]=374 [2002]=375 [2003]=376 [2004]=377
  [2005]=378 [2006]=379 [2007]=380 [2008]=381 [2009]=382
  [2010]=383 [2011]=384 [2012]=385 [2013]=386 [2014]=387
  [2015]=388 [2016]=389 [2017]=390 [2018]=391 [2019]=392
  [2020]=393 [2021]=394 [2022]=395 [2023]=396 [2024]=397
)

if [[ "${do_e3sm_compile,,}" == "true" ]]; then
  #################################################
  # Compile (only once for first member) 
  # input command (left to right):
  # do_fetch_code, do_create_newcase, ... 
  # do_case_setup, do_case_build, do_case_submit
  #################################################
  mkdir -p "${my_workdir}/spinup"
  cd "${my_workdir}/spinup" || exit 1 

  run_script="compile_and_setup_e3sm.sh"
  cp -rp "${my_workdir}/compile_and_setup_e3sm_${my_resolution}.${my_compset}.sh" "${run_script}"
  sed -i "s#\bmy_machine\b#${my_machine}#g"             "${run_script}"
  sed -i "s#\bmy_project\b#${my_project}#g"             "${run_script}"
  sed -i "s#\bmy_walltime\b#${my_walltime}#g"           "${run_script}"
  sed -i "s#\bmy_jobqueue\b#${my_jobqueue}#g"           "${run_script}"
  sed -i "s#\bmy_job_ntasks\b#${my_task_per_node}#g"    "${run_script}"
  sed -i "s#\bmy_task_per_node\b#${my_task_per_node}#g" "${run_script}"
  sed -i "s#\bmy_compset\b#${my_compset}#g"             "${run_script}"
  sed -i "s#\bmy_resolution\b#${my_resolution}#g"       "${run_script}"
  sed -i "s#\bmy_layout\b#${my_layout}#g"               "${run_script}"
  sed -i "s#\bmy_runtype\b#${my_runtype}#g"             "${run_script}"
  sed -i "s#\bmy_e3sm_code\b#${my_e3sm_code}#g"         "${run_script}"

  # only compile once for first member in first lead_date within this job 
  compiled_once=false
  my_modelexe=""

  for lead_date in "${my_leadymd[@]}"; do
    #determine the time for previous DA cycle 
    START_DATE=${lead_date}
    START_TOD=${my_leadtod}
    START_HH=$(printf "%02d" "${my_leadhh}")
    TIME_TAG="${lead_date//-/}${my_leadtod:0:2}"

    # create the run directory based on lead time 
    RUN_PATH="${my_runpath}/${my_refcase2}"
    if [ ! -d "${RUN_PATH}" ];then
      mkdir -p "${RUN_PATH}"
    fi

    CASE_ROOT="${RUN_PATH}/${TIME_TAG}/case_scripts"
    CASE_NAME=${my_refcase2}
    REF_PATH="${RUN_PATH}/${TIME_TAG}/archive/rest/${START_DATE}-${START_TOD}"

    cd "${my_workdir}/spinup" 
   
    # run script for each member
    tmp_script="spinup_${TIME_TAG}_${run_script}"
    cp "${run_script}" "${tmp_script}"

    sed -i "s#\bmy_startdate\b#${START_DATE}#g"       "${tmp_script}"
    sed -i "s#\bmy_starttod\b#${START_TOD}#g"         "${tmp_script}"
    sed -i "s#\bmy_refdate\b#${START_DATE}#g"         "${tmp_script}"
    sed -i "s#\bmy_reftod\b#${START_TOD}#g"           "${tmp_script}"
    sed -i "s#\bmy_refcase\b#${my_refcase2}#g"        "${tmp_script}"
    sed -i "s#\bmy_case_root\b#${CASE_ROOT}#g"        "${tmp_script}"
    sed -i "s#\bmy_casename\b#${CASE_NAME}#g"         "${tmp_script}"
    sed -i "s#\bmy_refdir\b#${REF_PATH}#g"            "${tmp_script}"

    if [[ "$compiled_once" == false ]]; then
      sed -i "s#\bold_modelexe\b##g"                  "${tmp_script}"
      chmod +x "${tmp_script}"
      ./${tmp_script} false true true true false

      my_modelexe="${RUN_PATH}/${TIME_TAG}/build/e3sm.exe"

      if [[ ! -f "${my_modelexe}" ]]; then
        echo $'\n----- e3sm.exe does not exit, compiling not success (check log file)-----\n'
        exit 1 
      fi

      compiled_once=true

    else 
      sed -i "s#old_modelexe#\"${my_modelexe}\"#g" "${tmp_script}"
      chmod +x "${tmp_script}"
      ./${tmp_script} false true true false false & 
    fi

    # wait to complete background members
    wait
  done
fi 

#function to modify namlist 
user_eam_nl() {
  local file="user_nl_eam"
  local ncdata_path="$1"
  local inithist_freq="$2"

  if [[ ! -f "$file" ]]; then
    echo "[ERROR] Namelist file '$file' not found!"
    return 1
  fi

  # Update existing entries (if present)
  ex "$file" <<EOF
g/^ *ncdata *=/s@=.*@= '${ncdata_path}'@
g/^ *inithist *=/s@=.*@= '${inithist_freq}'@
g/^ *inithist_all *=/s@=.*@= .true.@
wq
EOF

  # Append missing entries
  grep -Eq '^[[:space:]]*ncdata[[:space:]]*=' "$file" \
    || echo "ncdata = '${ncdata_path}'" >> "$file"

  grep -Eq '^[[:space:]]*inithist[[:space:]]*=' "$file" \
    || echo "inithist = '${inithist_freq}'" >> "$file"

  grep -Eq '^[[:space:]]*inithist_all[[:space:]]*=' "$file" \
    || echo "inithist_all = .true." >> "$file"
}

user_elm_nl() {
  local file="user_nl_elm"
  local finidat="$1"
  local yr_check="$2"
  local dynpft_check="$3"
  local fsurdat_check="$4"
  local pct_check="$5"

  if [[ ! -f "$file" ]]; then
    echo "[ERROR] Namelist file '$file' not found!"
    return 1
  fi

  # Update existing entries (if present)
  ex "$file" <<EOF
g/^ *finidat *=/s@=.*@= '${finidat}'@
g/^ *check_finidat_year_consistency *=/s@=.*@= ${yr_check}@
g/^ *check_dynpft_consistency *=/s@=.*@= ${dynpft_check}@
g/^ *check_finidat_fsurdat_consistency *=/s@=.*@= ${fsurdat_check}@
g/^ *check_finidat_pct_consistency *=/s@=.*@= ${pct_check}@
wq
EOF

  # Append missing entries
  grep -Eq '^[[:space:]]*finidat[[:space:]]*=' "$file" \
    || echo "finidat = '${finidat}'" >> "$file"

  grep -Eq '^[[:space:]]*check_finidat_year_consistency[[:space:]]*=' "$file" \
    || echo "check_finidat_year_consistency = ${yr_check}" >> "$file"

  grep -Eq '^[[:space:]]*check_dynpft_consistency[[:space:]]*=' "$file" \
    || echo "check_dynpft_consistency = ${dynpft_check}" >> "$file"

  grep -Eq '^[[:space:]]*check_finidat_fsurdat_consistency[[:space:]]*=' "$file" \
    || echo "check_finidat_fsurdat_consistency = ${fsurdat_check}" >> "$file"

  grep -Eq '^[[:space:]]*check_finidat_pct_consistency[[:space:]]*=' "$file" \
    || echo "check_finidat_pct_consistency = ${pct_check}" >> "$file"
}

user_mosart_nl() {
  local file="user_nl_mosart"
  local finidat_rtm="$1"

  if [[ ! -f "$file" ]]; then
    echo "[ERROR] Namelist file '$file' not found!"
    return 1
  fi

  # Update if exists
  ex "$file" <<EOF
g/^ *finidat_rtm *=/s@=.*@= '${finidat_rtm}'@
wq
EOF

  # Append if missing
  grep -Eq '^[[:space:]]*finidat_rtm[[:space:]]*=' "$file" \
    || echo "finidat_rtm = '${finidat_rtm}'" >> "$file"
}

# =====================================
# Customize MPAS stream files if needed
# =====================================
patch_mpaso_streams() {
  echo
  echo "Modifying MPAS (OCEAN) streams files"

  local run_dir="$1"
  local case_dir="$2"
  local f="${run_dir}/streams.ocean"
  local sm_dir="${case_dir}/SourceMods/src.mpaso"

  if [[ -z "${run_dir}" || -z "${case_dir}" ]]; then
    echo "ERROR: patch_mpaso_streams requires: <RUN_DIR> <CASE_DIR>"
    return 1
  fi

  if [[ ! -f "${f}" ]]; then
    echo "ERROR: Missing MPAS-O streams file:"
    echo "  ${f}"
    return 1
  fi

  mkdir -p "${sm_dir}" || { echo "ERROR: cannot create ${sm_dir}"; return 1; }

  # Apply patch non-interactively; prevent .rej files.
  # -N makes it idempotent (skip if already applied)
  # -r /dev/null prevents .rej files
  local patch_out patch_rc
  patch_out="$(
    patch --batch -N -p0 -r /dev/null -d "${run_dir}" 2>&1 <<'EOF'
--- streams.ocean
+++ streams.ocean
@@ -333,1 +333,1 @@
-        output_interval="00-00-05_00:00:00"
+        output_interval="00-00-01_00:00:00"
EOF
  )"
  patch_rc=$?

  if [[ $patch_rc -eq 0 ]]; then
    echo "streams.ocean patched."
  else
    # Treat "already applied / skipping" as success to avoid rejections
    if echo "${patch_out}" | grep -Eq 'Reversed \(or previously applied\) patch detected|Skipping patch|already applied'; then
      echo "streams.ocean already updated - skipping patch."
    else
      echo "ERROR: patch failed for ${f}"
      echo "${patch_out}"
      return 1
    fi
  fi

  # copy to SourceMods
  cp -p "${f}" "${sm_dir}/" || { echo "ERROR: failed to copy streams.ocean to ${sm_dir}"; return 1; }
}

patch_mpassi_streams() {
echo
echo 'Modifying MPAS streams files'
#pushd
#popd
}

# -------------------------------
# limit concurrent background jobs
# -------------------------------
wait_for_slot() {
  local max_jobs="$1"
  while true; do
    local njobs
    njobs=$(jobs -rp | wc -l | tr -d ' ')
    if (( njobs < max_jobs )); then
      break
    fi
    sleep 1
  done
}

user_mpaso_nl() {
  local file="user_nl_mpaso"
  local mpaso_dt="$1"

  if [[ ! -f "$file" ]]; then
    echo "[ERROR] Namelist file '$file' not found!"
    return 1
  fi

  if grep -Eq '^[[:space:]]*config_dt[[:space:]]*=' "$file"; then
    # Update existing entry
    ex "$file" <<EOF
g/^[[:space:]]*config_dt[[:space:]]*=/s@=.*@= '${mpaso_dt}'@
wq
EOF
  else
    # Append new entry
    echo "config_dt = '${mpaso_dt}'" >> "$file"
  fi
}

user_mpassi_nl() {
  local file="user_nl_mpassi"
  local mpassi_dt="$1"

  if [[ ! -f "$file" ]]; then
    echo "[ERROR] Namelist file '$file' not found!"
    return 1
  fi

  if grep -Eq '^[[:space:]]*config_dt[[:space:]]*=' "$file"; then
    # Update existing entry
    ex "$file" <<EOF
g/^[[:space:]]*config_dt[[:space:]]*=/s@=.*@= ${mpassi_dt}@
wq
EOF
  else
    # Append new entry
    echo "config_dt = ${mpassi_dt}" >> "$file"
  fi
}

# -------------------------------
# setup one ensemble member
# -------------------------------
setup_one_member() {
  local lead_date="$1"

  local START_DATE START_TOD START_HH TIME_TAG RUN_PATH
  local CASE_NAME CASE_ROOT RUN_ROOT ARCHIVE_DIR
  local SRC_ROOT TGT_ROOT
  local scomp smod
  local REST_DATE_EXT REST_CASE REST_TAG REST_FILE RPOT_FILE
  local atm_in="" lnd_in="" rof_in=""

  START_DATE="${lead_date}"
  START_TOD="${my_leadtod}"
  START_HH=$(printf "%02d" "${my_leadhh}")
  TIME_TAG="${lead_date//-/}${my_leadtod:0:2}"
  RUN_PATH="${my_runpath}/${my_refcase2}"

  CASE_NAME="${my_refcase2}"
  CASE_ROOT="${RUN_PATH}/${TIME_TAG}/case_scripts"
  RUN_ROOT="${RUN_PATH}/${TIME_TAG}/run"
  ARCHIVE_DIR="${RUN_PATH}/${TIME_TAG}/archive"     
  TGT_ROOT="${ARCHIVE_DIR}/rest/${START_DATE}-${START_TOD}"

  echo "=== Setting up ${CASE_NAME} ==="

  mkdir -p "${TGT_ROOT}" || {
    echo "ERROR: Failed to create ${TGT_ROOT}"
    return 1
  }

  for scomp in atm lnd rof ocn ice drv; do
    smod="${comp_map[$scomp]}"

    if [[ ${scomp} == "atm" ]]; then
      REST_CASE="${my_spinup_atmcase}"
      REST_TAG="${smod}.i"
      SRC_DATE_EXT=$(printf "%02d" $((START_TOD / 3600)))
      SRC_ROOT="${my_spinup_atmdir}"
      SRC_FILE="${REST_CASE}.${REST_TAG}.${SRC_DATE_EXT}.nc"
      if [[ -n "${my_spinup_atmtag}" ]]; then
          SRC_FILE="${REST_CASE}.${START_DATE}.${SRC_DATE_EXT}.${my_spinup_atmtag}.nc"
      fi
    elif [[ ${scomp} == "ocn" ]]; then
      START_YEAR=${START_DATE:0:4}
      if [[ -n "${YEAR_MAP[$START_YEAR]}" ]]; then
          MODEL_YEAR=${YEAR_MAP[$START_YEAR]}
      else
          echo "BruteForce Ocean: Year ${START_YEAR} not in supported list"
          exit 1
      fi
      REST_CASE="${my_spinup_ocncase}"
      REST_TAG="${smod}.rst"
      SRC_DATE_EXT=$(printf "%04d" "$MODEL_YEAR")${START_DATE:4}
      SRC_ROOT="${my_spinup_ocndir}"
      SRC_FILE="${REST_CASE}.${REST_TAG}.${SRC_DATE_EXT}_${START_TOD}.nc"
      if [[ -n "${my_spinup_ocntag}" ]]; then
          SRC_FILE="${REST_CASE}.${REST_TAG}.${SRC_DATE_EXT}_${START_TOD}.${my_spinup_ocntag}.nc"
      fi
    else
      START_YEAR=${START_DATE:0:4}
      if [[ -n "${YEAR_MAP[$START_YEAR]}" ]]; then
          MODEL_YEAR=${YEAR_MAP[$START_YEAR]}
      else
          echo "BruteForce FOSIRL: Year ${START_YEAR} not in supported list"
          exit 1
      fi
      REST_CASE="${my_spinup_foscase}"
      if [[ ${scomp} == "ocn" || ${scomp} == "ice" ]]; then
         SRC_DATE_EXT=$(printf "%04d" "$MODEL_YEAR")${START_DATE:4}_${START_TOD}
         REST_TAG="${smod}.rst"
      else
         SRC_DATE_EXT=$(printf "%04d" "$MODEL_YEAR")${START_DATE:4}-${START_TOD}
         REST_TAG="${smod}.r"
      fi 
      SRC_ROOT="${my_spinup_fosdir}/$(printf "%04d" "$MODEL_YEAR")${START_DATE:4}-${START_TOD}"
      SRC_FILE="${REST_CASE}.${REST_TAG}.${SRC_DATE_EXT}.nc"
    fi

    if [[ ${scomp} == "ocn" || ${scomp} == "ice" ]]; then
        REST_DATE_EXT=${START_DATE}_${START_TOD}
    else
        REST_DATE_EXT=${START_DATE}-${START_TOD}
    fi 

    REST_FILE="${REST_CASE}.${REST_TAG}.${REST_DATE_EXT}.nc"
    RPOT_FILE="rpointer.${scomp}"
    START_HOUR=$(printf "%02d" $((START_TOD / 3600)))

    if [[ ! -f "${SRC_ROOT}/${SRC_FILE}" ]]; then
      echo "ERROR: Missing restart file:"
      echo "  ${SRC_ROOT}/${SRC_FILE}"
      return 1
    fi

    # call isolated helper script
    cd ${my_workdir} || exit 1

    ./e3sm_spinup_setup.sh \
       "${my_env}"\
       "${scomp}" \
       "${SRC_ROOT}" \
       "${SRC_FILE}" \
       "${TGT_ROOT}" \
       "${REST_FILE}" \
       "${START_DATE}" \
       "${START_HOUR}" \
       "${RPOT_FILE}"

    rc=$?

    if [[ $rc -ne 0 ]]; then
        echo "ERROR: helper failed for ${scomp}"
        return $rc
    fi

    if [[ ${scomp} == "atm" ]]; then
      atm_in="${TGT_ROOT}/${REST_FILE}"
    elif [[ ${scomp} == "lnd" ]]; then
      lnd_in="${TGT_ROOT}/${REST_FILE}"
    elif [[ ${scomp} == "rof" ]]; then
      rof_in="${TGT_ROOT}/${REST_FILE}"
    fi

  done

  cd "${CASE_ROOT}" || {
    echo "ERROR: Cannot cd to ${CASE_ROOT}"
    return 1
  }

  ./xmlchange run_exe="--kill-on-bad-exit=1 --job-name=${CASE_NAME} \${EXEROOT}/e3sm.exe " || return 1
  ./xmlchange RUN_TYPE="${my_runtype,,}" || return 1
  ./xmlchange GET_REFCASE=TRUE || return 1
  ./xmlchange RUN_REFDIR="${TGT_ROOT}" || return 1
  ./xmlchange RUN_REFCASE="${my_refcase2}" || return 1
  ./xmlchange RUN_REFDATE="${START_DATE}" || return 1
  ./xmlchange RUN_REFTOD="${START_TOD}" || return 1
  ./xmlchange DOUT_S="${my_short_archive,,}" || return 1
  ./xmlchange DOUT_S_ROOT="${ARCHIVE_DIR}" || return 1
  ./xmlchange RUN_STARTDATE="${START_DATE}" || return 1
  ./xmlchange START_TOD="${START_TOD}" || return 1
  ./xmlchange REST_OPTION="${REST_OPTION}" || return 1
  ./xmlchange REST_N="${REST_N}" || return 1
  ./xmlchange STOP_OPTION="${STOP_OPTION}" || return 1
  ./xmlchange STOP_N="${STOP_N}" || return 1
  ./xmlchange RESUBMIT=${RESUBMIT}
  ./xmlchange JOB_WALLCLOCK_TIME="${my_walltime}" || return 1
  ./xmlchange JOB_QUEUE="${my_jobqueue}" || return 1

  # change model time step if requested
  if [[ -n "${my_mpaso_dt}" || -n "${my_mpassi_dt}" ]]; then
    [[ -n "${my_mpaso_dt}"  ]] && user_mpaso_nl  "${my_mpaso_dt}"
    [[ -n "${my_mpassi_dt}" ]] && user_mpassi_nl "${my_mpassi_dt}"
  fi

  echo "MODEL_START_TYPE = ${my_runtype}"
  echo "RUN_REFDIR        = ${TGT_ROOT}"
  echo "RUN_REFCASE       = ${my_refcase2}"
  echo "RUN_REFDATE       = ${START_DATE}"

  user_eam_nl "${atm_in}" "${my_initopt}" || return 1
  user_elm_nl "${lnd_in}" .false. .false. .false. .false. || return 1
  user_mosart_nl "${rof_in}" || return 1

  ./case.setup || return 1
  ./xmlchange BUILD_COMPLETE=TRUE || return 1

  patch_mpaso_streams "${RUN_ROOT}" "${CASE_ROOT}" || return 1

  echo "=== Finished ${CASE_NAME} ==="
}

continue_one_member() {
  local lead_date="$1"

  local TIME_TAG RUN_PATH
  local CASE_NAME CASE_ROOT RUN_ROOT ARCHIVE_DIR
  local SRC_DIR

  TIME_TAG="${lead_date//-/}${my_leadtod:0:2}"
  RUN_PATH="${my_runpath}/${my_refcase2}"

  echo "=== Setup for ${lead_date} ==="

  CASE_NAME="${my_refcase2}"
  CASE_ROOT="${RUN_PATH}/${TIME_TAG}/case_scripts"
  RUN_ROOT="${RUN_PATH}/${TIME_TAG}/run"
  ARCHIVE_DIR="${RUN_PATH}/${TIME_TAG}/archive"

  SRC_DIR="${ARCHIVE_DIR}/rest/${my_restart_ymd}-${my_restart_tod}"

  if [[ ! -d "${SRC_DIR}" ]]; then
    echo "ERROR: Restart directory does not exist:"
    echo "  ${SRC_DIR}"
    return 1
  fi
  
  if ! compgen -G "${SRC_DIR}/*" > /dev/null; then
    echo "ERROR: Restart directory is empty:"
    echo "  ${SRC_DIR}"
    return 1
  fi

  # copy restart files to run directory
  cp -rp "${SRC_DIR}/"* "${RUN_ROOT}/" || {
    echo "ERROR: Failed to copy restart files from ${SRC_DIR} to ${RUN_ROOT}"
    return 1
  }

  # now setup model for restart run
  cd "${CASE_ROOT}" || {
    echo "ERROR: Cannot cd to ${CASE_ROOT}"
    return 1
  }

  ./xmlchange run_exe="--kill-on-bad-exit=1 --job-name=${CASE_NAME} \${EXEROOT}/e3sm.exe " || return 1
  ./xmlchange CONTINUE_RUN=TRUE || return 1
  ./xmlchange DOUT_S="${my_short_archive,,}" || return 1
  ./xmlchange DOUT_S_ROOT="${ARCHIVE_DIR}" || return 1
  ./xmlchange REST_OPTION="${REST_OPTION}" || return 1
  ./xmlchange REST_N="${REST_N}" || return 1
  ./xmlchange STOP_OPTION="${STOP_OPTION}" || return 1
  ./xmlchange STOP_N="${STOP_N}" || return 1
  ./xmlchange RESUBMIT=${RESUBMIT} || return 1
  ./xmlchange JOB_WALLCLOCK_TIME="${my_walltime}" || return 1
  ./xmlchange JOB_QUEUE="${my_jobqueue}" || return 1

  ./case.setup || return 1

  echo "=== Finished member ${i} for ${lead_date} ==="
}

#####################################################
# do heavy loops and run 
#####################################################

if [[ "${do_ensemble_setup,,}" == "true" ]]; then
  declare -A comp_map=(
    ["atm"]="eam"
    ["lnd"]="elm"
    ["rof"]="mosart"
    ["ocn"]="mpaso"
    ["ice"]="mpassi"
    ["drv"]="cpl"
  )

  max_setup_jobs=10
  fail_flag=0
  pids=()

  for lead_date in "${my_leadymd[@]}"; do
    wait_for_slot "${max_setup_jobs}"

    (
      setup_one_member "${lead_date}"
    ) &
    pids+=($!)
  done

  # wait for all jobs across all lead_date
  for pid in "${pids[@]}"; do
    if ! wait "${pid}"; then
      fail_flag=1
    fi
  done

  if (( fail_flag != 0 )); then
    echo "ERROR: one or more setup jobs failed."
    exit 1
  fi

fi

if [[ "${do_continue_run,,}" == "true" ]]; then
  max_continue_jobs=10

  for lead_date in "${my_leadymd[@]}"; do
    fail_flag=0

    pids=()

    for (( i=0; i<my_ensnum; i++ )); do
      wait_for_slot "${max_continue_jobs}"

      (
        continue_one_member "${lead_date}" "${i}"
      ) &
      pids+=($!)
    done

    # wait for all jobs for this lead_date
    for pid in "${pids[@]}"; do
      if ! wait "${pid}"; then
        fail_flag=1
      fi
    done

    if (( fail_flag != 0 )); then
      echo "ERROR: one or more continue-run setup jobs failed for ${lead_date}."
      exit 1
    fi
  done
fi

if [[ "${do_ensemble_run,,}" == "true" ]]; then
  # ===========================================================
  # First step: Setup maximum allowed jobs based on 
  # my_job_nnodes: total number of nodes 
  # my_min_nodes_per_sim: minimum nodes for one coupled E3SM run  
  # ============================================================
  NMAXPS=$((my_job_nnodes / my_min_nodes_per_sim))
  (( NMAXPS < 1 )) && NMAXPS=1

  pids=()

  for lead_date in "${my_leadymd[@]}"; do
    TIME_TAG="${lead_date//-/}${my_leadtod:0:2}"
    RUN_PATH="${my_runpath}/${my_refcase2}"

    CASE_NAME="${my_refcase2}"
    CASE_ROOT="${RUN_PATH}/${TIME_TAG}/case_scripts"
    RUN_ROOT="${RUN_PATH}/${TIME_TAG}/run"
    ARCHIVE_DIR="${RUN_PATH}/${TIME_TAG}/archive"
    
    ##################
    #run model 
    ##################
    cd "${CASE_ROOT}" || exit 1
    ./case.submit --no-batch > e3sm_${TIME_TAG}.log.o${SLURM_JOB_ID} 2>&1 &

    pids+=($!)
    # enforce global concurrency limit across all lead dates + members
    if (( ${#pids[@]} >= NMAXPS )); then
      for pid in "${pids[@]}"; do
        wait "$pid" || echo "WARNING: member job $pid failed"
      done
      pids=()
    fi

    echo ============================

  done 

  if (( ${#pids[@]} > 0 )); then
    for pid in "${pids[@]}"; do
      wait "$pid" || echo "WARNING: member job $pid failed"
    done
    pids=()
  fi

  # Wait loop with external hook
  k=0
  while true
  do

    sleep 60

    # Execute extra instructions
    cd ${my_workdir} || exit 1
    . ./e3sm_boundle_extra.sh 

    # List running background processes.
    # (Needed for the stop clause below to work)
    k=$((k+1))
    if (( k % 5 == 0 ))
    then
      echo ============================
      date
      jobs -l
      echo ----------------------------
      squeue --job=${SLURM_JOBID} --steps
      echo ============================
    fi

    # Stop when all processes are done
    n=$(jobs -l | wc -l)
    if (( n == 0 ))
    then
       echo ============================
       date
       echo No running jobs left
       echo ============================
       break
    fi
  done

  # Post steps
  cd "${my_workdir}" || exit 1
  . ./e3sm_boundle_cycling.sh

  # That's all folks!
  sleep 10
fi 

echo "======================================================="
echo "  End of Ensemble E3SM Simulations (bundled)"
echo "  Finish time: $(date)"
echo "======================================================="
