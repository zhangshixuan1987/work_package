#!/bin/bash -el
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --account=e3sm
#SBATCH  --constraint=cpu
#SBATCH  --qos=regular
#SBATCH  --nodes=280
#SBATCH  --time=16:00:00
#SBATCH  --job-name=e3sm_s2d_ensrun

source ./create_and_setup_case.sh

echo "============================================"
echo "  Run Ensemble E3SM Simulations (bundled)"
echo "  Start time: $(date)"
echo "============================================"

if [[ "${do_compile_setup_only,,}" == "true" ]]; then
  if [[ "${do_ensemble_run,,}" == "true" ]]; then
    echo "ERROR: Workflow assumes do_ensemble_run=true only after setup is completed"
    echo "ERROR: if the setup is compelte, set do_compile_setup_only=false"
    echo "ERROR: if the setup is not compelte, set do_ensemble_run=false"
    return 1
  fi
  do_e3sm_compile=true
  do_ensemble_setup=true
  do_continue_setup=false
fi

if [[ "${do_continue_setup,,}" == "true" ]];then
  if [[ "${do_ensemble_run,,}" == "true" ]]; then
    echo "ERROR: Workflow assumes do_ensemble_run=true only after setup is completed"
    echo "ERROR: if the setup is compelte, set do_compile_setup_only=false"
    echo "ERROR: if the setup is not compelte, set do_ensemble_run=false"
    return 1
  fi
  do_e3sm_compile=${my_continue_recompile}
  do_ensemble_setup=false
  do_ensemble_run=false
fi

if [[ "${do_ensemble_run,,}" == "true" ]]; then
  if [[ "${do_continue_setup,,}" == "true" || "${do_compile_setup_only,,}" == "true" ]]; then
    echo "WARNING: do_ensemble_run=true assumes prior setup is complete"
    echo "WARNING: Forcing do_continue_setup, do_e3sm_compile, do_ensemble_setup=false"
  fi
  do_e3sm_compile=false
  do_ensemble_setup=false
  do_continue_setup=false
fi

if [[ "${do_e3sm_compile,,}" == "true" ]]; then
  #################################################
  # Compile (only once for first member) 
  # input command (left to right):
  # do_fetch_code, do_create_newcase, ... 
  # do_case_setup, do_case_build, do_case_submit
  #################################################
  mkdir -p "${my_workdir}/scripts"
  cd "${my_workdir}/scripts" || exit 1 

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
    RUN_PATH="${my_runpath}/${my_leadcase}_${TIME_TAG}"
    if [ ! -d "${RUN_PATH}" ];then
      mkdir -p "${RUN_PATH}"
    fi

    sed -i "s#\bmy_startdate\b#${START_DATE}#g"         "${run_script}"
    sed -i "s#\bmy_starttod\b#${START_TOD}#g"           "${run_script}"
    sed -i "s#\bmy_refdate\b#${START_DATE}#g"           "${run_script}"
    sed -i "s#\bmy_reftod\b#${START_TOD}#g"             "${run_script}"
    sed -i "s#\bmy_refcase\b#${my_refcase2}#g"          "${run_script}"

    for (( i=0; i<my_ensnum; i++ )); do
      echo === Compile member ${i} ===
      ENSTR=EN$(printf "%02d" "$i")
      CASE_ROOT=${RUN_PATH}/${ENSTR}
      CASE_NAME=${my_leadcase}_${TIME_TAG}.${ENSTR}
      REF_PATH="${RUN_PATH}/${ENSTR}/archive/rest/${START_DATE}-${START_TOD}"

      cd "${my_workdir}/scripts" 
   
      # run script for each member
      tmp_script="${ENSTR}_${TIME_TAG}_${run_script}"
      cp "${run_script}" "${tmp_script}"

      sed -i "s#\bmy_case_root\b#${CASE_ROOT}#g"        "${tmp_script}"
      sed -i "s#\bmy_casename\b#${CASE_NAME}#g"         "${tmp_script}"
      sed -i "s#\bmy_refdir\b#${REF_PATH}#g"            "${tmp_script}"

      if [[ "$compiled_once" == false && $i -eq 0 ]]; then
        sed -i "s#\bold_modelexe\b##g"                  "${tmp_script}"
        chmod +x "${tmp_script}"
        ./${tmp_script} false true true true false

        my_modelexe="${RUN_PATH}/EN00/build/e3sm.exe"

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

    done 
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

# -------------------------------
# setup one ensemble member
# -------------------------------
setup_one_member() {
  local lead_date="$1"
  local i="$2"
  local year=$((10#${lead_date:0:4}))

  local START_DATE START_TOD START_HH TIME_TAG RUN_PATH
  local ENSTR CASE_NAME CASE_ROOT RUN_ROOT ARCHIVE_DIR
  local SRC_ROOT TGT_ROOT
  local scomp smod
  local REST_DATE_EXT REST_CASE REST_TAG REST_FILE RPOT_FILE
  local atm_in="" lnd_in="" rof_in=""

  START_DATE="${lead_date}"
  START_TOD="${my_leadtod}"
  START_HH=$(printf "%02d" "${my_leadhh}")
  TIME_TAG="${lead_date//-/}${my_leadtod:0:2}"
  RUN_PATH="${my_runpath}/${my_leadcase}_${TIME_TAG}"

  ENSTR="EN$(printf "%02d" "$i")"
  CASE_NAME="${my_leadcase}_${TIME_TAG}.${ENSTR}"
  CASE_ROOT="${RUN_PATH}/${ENSTR}/case_scripts"
  RUN_ROOT="${RUN_PATH}/${ENSTR}/run"
  ARCHIVE_DIR="${RUN_PATH}/${ENSTR}/archive"

  SRC_ROOT="${my_refdir}/${START_DATE}-${START_TOD}"
  TGT_ROOT="${ARCHIVE_DIR}/rest/${START_DATE}-${START_TOD}"

  echo "=== Setting up ${CASE_NAME} ==="

  mkdir -p "${TGT_ROOT}" || {
    echo "ERROR: Failed to create ${TGT_ROOT}"
    return 1
  }

  for scomp in atm lnd rof ocn ice drv; do
    smod="${comp_map[$scomp]}"

    if [[ ${scomp} == "atm" ]]; then
      REST_DATE_EXT="${START_DATE}-${START_TOD}"
      REST_CASE="${my_refcase1}"
      REST_TAG="${ENSTR}.${smod}.i"
    elif [[ ${scomp} == "ocn" || ${scomp} == "ice" ]]; then
      REST_DATE_EXT="${START_DATE}_${START_TOD}"
      REST_CASE="${my_refcase2}"
      REST_TAG="${smod}.rst"
    else
      REST_DATE_EXT="${START_DATE}-${START_TOD}"
      REST_CASE="${my_refcase2}"
      REST_TAG="${smod}.r"
    fi

    REST_FILE="${REST_CASE}.${REST_TAG}.${REST_DATE_EXT}.nc"
    RPOT_FILE="rpointer.${scomp}"

    if [[ ! -f "${SRC_ROOT}/${REST_FILE}" ]]; then
      echo "ERROR: Missing restart file:"
      echo "  ${SRC_ROOT}/${REST_FILE}"
      return 1
    fi

    if [[ ! -f "${TGT_ROOT}/${REST_FILE}" ]]; then
      cp -p "${SRC_ROOT}/${REST_FILE}" "${TGT_ROOT}/" || {
        echo "ERROR: Failed to copy ${REST_FILE}"
        return 1
      }
    fi

    if [[ "${scomp}" != "atm" ]]; then
      if [[ ! -f "${SRC_ROOT}/${RPOT_FILE}" ]]; then
        echo "ERROR: Missing rpointer file:"
        echo "  ${SRC_ROOT}/${RPOT_FILE}"
        return 1
      fi

      if [[ ! -f "${TGT_ROOT}/${RPOT_FILE}" ]]; then
        cp -p "${SRC_ROOT}/${RPOT_FILE}" "${TGT_ROOT}/" || {
          echo "ERROR: Failed to copy ${RPOT_FILE}"
          return 1
        }
      fi
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
  ./xmlchange REST_OPTION="${my_restopt}" || return 1
  ./xmlchange REST_N="${my_restn}" || return 1
  ./xmlchange STOP_OPTION="${my_runopt}" || return 1
  ./xmlchange STOP_N="${my_runn}" || return 1
  ./xmlchange JOB_WALLCLOCK_TIME="${my_walltime}" || return 1
  ./xmlchange JOB_QUEUE="${my_jobqueue}" || return 1

  echo "MODEL_START_TYPE = ${my_runtype}"
  echo "RUN_REFDIR        = ${TGT_ROOT}"
  echo "RUN_REFCASE       = ${my_refcase2}"
  echo "RUN_REFDATE       = ${START_DATE}"

  user_eam_nl "${atm_in}" "${my_initopt}" || return 1
  user_elm_nl "${lnd_in}" .false. .false. .false. .false. || return 1
  user_mosart_nl "${rof_in}" || return 1

  if (( year > 2014 )); then 
    #change the forcing file to extend simulation beyond 2014
    cat "${my_workdir}/ssp245_user_nl_eam.txt" >> user_nl_eam
    cat "${my_workdir}/ssp245_user_nl_elm.txt" >> user_nl_elm
  fi 

  ./case.setup || return 1
  ./xmlchange BUILD_COMPLETE=TRUE || return 1

  patch_mpaso_streams "${RUN_ROOT}" "${CASE_ROOT}" || return 1

  echo "=== Finished ${CASE_NAME} ==="
}

continue_one_member() {
  local lead_date="$1"
  local rest_date="$2"
  local i="$3"

  local scomp smod
  local TIME_TAG RUN_PATH
  local ENSTR CASE_NAME CASE_ROOT RUN_ROOT ARCHIVE_DIR
  local SRC_DIR
  local year=$((10#${rest_date:0:4}))
  local rest_tod=${my_restart_tod}

  TIME_TAG="${lead_date//-/}${my_leadtod:0:2}"
  RUN_PATH="${my_runpath}/${my_leadcase}_${TIME_TAG}"

  echo "=== Setup member ${i} for ${lead_date} ==="

  ENSTR="EN$(printf "%02d" "$i")"
  CASE_NAME="${my_leadcase}_${TIME_TAG}.${ENSTR}"
  CASE_ROOT="${RUN_PATH}/${ENSTR}/case_scripts"
  RUN_ROOT="${RUN_PATH}/${ENSTR}/run"
  ARCHIVE_DIR="${RUN_PATH}/${ENSTR}/archive"
  REST_CASE="${CASE_NAME}"

  SRC_DIR="${ARCHIVE_DIR}/rest/${rest_date}-${rest_tod}"
  REST_DATE_PREFIX="${rest_date:0:4}-${rest_date:5:2}"   # YYYY-MM

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

  # copy the hist files potentially needed for restart run
  for scomp in atm lnd rof drv; do
    smod="${comp_map[$scomp]}"
    HIS_DIR="${ARCHIVE_DIR}/${scomp}/hist"
    for htag in h1 h2 h3 h4 h5 h6 h7; do
      REST_TAG="${smod}.${htag}"
      REST_GLOB="${HIS_DIR}/${REST_CASE}.${REST_TAG}.${REST_DATE_PREFIX}"*.nc

      shopt -s nullglob
      files=( ${REST_GLOB} )
      shopt -u nullglob

      if [[ ${#files[@]} -eq 0 ]]; then
        echo "WARNING: No history files found for:"
        echo "  ${REST_GLOB}"
        continue
      fi

      # sort and get the last file
      last_file=$(printf '%s\n' "${files[@]}" | sort | tail -n 1)

      # extract date part from filename
      last_base=$(basename "$last_file")
      last_date=$(echo "$last_base" | sed -E 's/.*\.([0-9]{4}-[0-9]{2}-[0-9]{2})(-[0-9]+)?\.nc/\1/')
      echo "Last file: $last_file"
      echo "Last date: $last_date"

      # destination path
      dest_file="${RUN_ROOT}/${last_base}"

      # copy only if not already present (or different)
      if [[ ! -f "$dest_file" ]]; then
         cp -rp "$last_file" "$dest_file" 
      else
         echo "Skipping copy (already exists): $dest_file"
      fi

    done
  done

  # now setup model for restart run
  cd "${CASE_ROOT}" || {
    echo "ERROR: Cannot cd to ${CASE_ROOT}"
    return 1
  }

  ./xmlchange run_exe="--kill-on-bad-exit=1 --job-name=${CASE_NAME} \${EXEROOT}/e3sm.exe " || return 1
  ./xmlchange CONTINUE_RUN=TRUE || return 1
  ./xmlchange RUN_STARTDATE="${lead_date}" || return 1
  ./xmlchange START_TOD="${my_leadtod}" || return 1
  ./xmlchange DOUT_S="${my_short_archive,,}" || return 1
  ./xmlchange DOUT_S_ROOT="${ARCHIVE_DIR}" || return 1
  ./xmlchange REST_OPTION="${my_continue_restopt}" || return 1
  ./xmlchange REST_N="${my_continue_restn}" || return 1
  ./xmlchange STOP_OPTION="${my_continue_runopt}" || return 1
  ./xmlchange STOP_N="${my_continue_runn}" || return 1
  ./xmlchange JOB_WALLCLOCK_TIME="${my_continue_walltime}" || return 1
  ./xmlchange JOB_QUEUE="${my_jobqueue}" || return 1

  if (( year > 2014 )); then
    #change the forcing file to extend simulation beyond 2014
    cat "${my_workdir}/ssp245_user_nl_eam.txt" >> user_nl_eam
    cat "${my_workdir}/ssp245_user_nl_elm.txt" >> user_nl_elm
  fi

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

  for i in "${!my_leadymd[@]}"; do
    lead_date="${my_leadymd[$i]}"

    fail_flag=0
    
    pids=()

    for (( i=0; i<my_ensnum; i++ )); do
      wait_for_slot "${max_setup_jobs}"

      (
        setup_one_member "${lead_date}" "${i}"
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
      echo "ERROR: one or more setup jobs failed for ${lead_date}."
      exit 1
    fi
  done
fi

if [[ "${do_continue_setup,,}" == "true" ]]; then
  declare -A comp_map=( 
    ["atm"]="eam"
    ["lnd"]="elm"
    ["rof"]="mosart"
    ["ocn"]="mpaso"
    ["ice"]="mpassi"
    ["drv"]="cpl"
  ) 

  max_continue_jobs=10

  for i in "${!my_leadymd[@]}"; do
    lead_date="${my_leadymd[$i]}"
    rest_date="${my_restart_ymd[$i]}"

    fail_flag=0

    pids=()

    for (( i=0; i<my_ensnum; i++ )); do
      wait_for_slot "${max_continue_jobs}"

      (
        continue_one_member "${lead_date}" "${rest_date}" "${i}"
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
    RUN_PATH=${my_runpath}/${my_leadcase}_${TIME_TAG}

    # Loop over members and run ensembel forecast
    for (( i=0; i<my_ensnum; i++ )); do
      echo === Run member ${i} ===
      ENSTR=EN$(printf "%02d" "$i")
      CASE_NAME=${my_leadcase}_${TIME_TAG}.${ENSTR}
      CASE_ROOT=${RUN_PATH}/${ENSTR}/case_scripts
      RUN_ROOT=${RUN_PATH}/${ENSTR}/run
      ARCHIVE_DIR="${RUN_PATH}/${ENSTR}/archive"

      ##################
      #run model 
      ##################
      cd "${CASE_ROOT}" || exit 1
      ./case.submit --no-batch > e3sm_${ENSTR}.log.o${SLURM_JOB_ID} 2>&1 &

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
