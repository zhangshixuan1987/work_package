#!/bin/bash -el
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --account=e3sm
#SBATCH  --constraint=cpu
#SBATCH  --qos=regular
#SBATCH  --nodes=280
#SBATCH  --time=18:00:00
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
    exit 1
  fi
  do_e3sm_compile=true
  do_ensemble_setup=true
  do_continue_setup=false
fi

if [[ "${do_continue_setup,,}" == "true" ]]; then
  # Step 1 of two-step restart: reconfigure existing cases for the restart run.
  # After this completes and you have verified the setup, set
  #   do_continue_setup=false  do_ensemble_run=true
  # and resubmit to run.
  if [[ "${do_short_spinup,,}" == "true" ]]; then
    echo "ERROR: do_short_spinup must be false when do_continue_setup=true"
    exit 1
  fi
  if [[ "${do_ensemble_run,,}" == "true" ]]; then
    echo "ERROR: do_continue_setup and do_ensemble_run cannot both be true"
    echo "ERROR: use do_continue_setup=true first (setup), then do_ensemble_run=true (run)"
    exit 1
  fi
  do_e3sm_compile=${my_continue_recompile}
  do_ensemble_setup=false
  do_ensemble_run=false
fi

if [[ "${do_ensemble_run,,}" == "true" ]]; then
  # Step 2 (or standalone run): assumes prior setup is already complete.
  if [[ "${do_compile_setup_only,,}" == "true" ]]; then
    echo "WARNING: do_ensemble_run=true assumes prior setup is complete"
    echo "WARNING: Forcing do_compile_setup_only, do_e3sm_compile, do_ensemble_setup=false"
  fi
  do_e3sm_compile=false
  do_ensemble_setup=false
  do_continue_setup=false
fi

preflight_check() {
  local fail=0
  local required_file
  local lead_date rest_date year time_tag run_path ensstr case_root member_idx
  local total_common spinup_common rest_common
  local template_script="${my_workdir}/compile_and_setup_e3sm_${my_resolution}.${my_compset}.sh"
  local need_ssp=false

  echo "----- Preflight checks -----"

  for required_file in "${my_workdir}/create_and_setup_case.sh"; do
    if [[ ! -f "${required_file}" ]]; then
      echo "ERROR: Missing required file: ${required_file}"
      fail=1
    fi
  done

  if [[ "${do_e3sm_compile,,}" == "true" && ! -f "${template_script}" ]]; then
    echo "ERROR: Missing compile/setup template:"
    echo "  ${template_script}"
    fail=1
  fi

  if [[ "${do_ensemble_run,,}" == "true" ]]; then
    for required_file in \
      "${my_workdir}/e3sm_boundle_extra.sh" \
      "${my_workdir}/e3sm_boundle_cycling.sh"; do
      if [[ ! -f "${required_file}" ]]; then
        echo "ERROR: Missing bundled-run hook file: ${required_file}"
        fail=1
      fi
    done
  fi

  if ! [[ "${my_ensnum}" =~ ^[0-9]+$ ]] || (( my_ensnum < 1 )); then
    echo "ERROR: my_ensnum must be a positive integer; got '${my_ensnum}'"
    fail=1
  fi

  if (( ${#my_leadymd[@]} < 1 )); then
    echo "ERROR: my_leadymd must contain at least one lead date"
    fail=1
  fi

  if [[ "${do_ensemble_setup,,}" == "true" || "${do_ensemble_run,,}" == "true" ]]; then
    case "${my_runopt}" in
      nyears|nmonths|ndays) ;;
      *)
        echo "ERROR: my_runopt must be one of nyears, nmonths, ndays; got '${my_runopt}'"
        fail=1
        ;;
    esac
    if ! [[ "${my_runn}" =~ ^[0-9]+$ ]] || (( my_runn < 1 )); then
      echo "ERROR: my_runn must be a positive integer; got '${my_runn}'"
      fail=1
    fi

    case "${my_restopt}" in
      nyears|nmonths|ndays) ;;
      *)
        echo "ERROR: my_restopt must be one of nyears, nmonths, ndays; got '${my_restopt}'"
        fail=1
        ;;
    esac
    if ! [[ "${my_restn}" =~ ^[0-9]+$ ]] || (( my_restn < 1 )); then
      echo "ERROR: my_restn must be a positive integer; got '${my_restn}'"
      fail=1
    fi
  fi

  if [[ "${do_ensemble_run,,}" == "true" ]]; then
    if ! [[ "${my_job_nnodes}" =~ ^[0-9]+$ ]] || (( my_job_nnodes < 1 )); then
      echo "ERROR: my_job_nnodes must be a positive integer; got '${my_job_nnodes}'"
      fail=1
    fi
    if ! [[ "${my_min_nodes_per_sim}" =~ ^[0-9]+$ ]] || (( my_min_nodes_per_sim < 1 )); then
      echo "ERROR: my_min_nodes_per_sim must be a positive integer; got '${my_min_nodes_per_sim}'"
      fail=1
    fi
    if [[ "${do_short_spinup,,}" == "true" ]]; then
      if [[ -z "${my_spinup_runopt}" ]]; then
        echo "ERROR: my_spinup_runopt must be set when do_short_spinup=true"
        fail=1
      elif [[ "${my_spinup_runopt}" != "nyears" &&
              "${my_spinup_runopt}" != "nmonths" &&
              "${my_spinup_runopt}" != "ndays" ]]; then
        echo "ERROR: my_spinup_runopt must be one of nyears, nmonths, ndays; got '${my_spinup_runopt}'"
        fail=1
      fi
      if ! [[ "${my_spinup_runn}" =~ ^[0-9]+$ ]] || (( my_spinup_runn < 1 )); then
        echo "ERROR: my_spinup_runn must be a positive integer; got '${my_spinup_runn}'"
        fail=1
      fi
      if [[ -z "${my_spinup_restopt}" ]]; then
        echo "ERROR: my_spinup_restopt must be set when do_short_spinup=true"
        fail=1
      elif [[ "${my_spinup_restopt}" != "nyears" &&
              "${my_spinup_restopt}" != "nmonths" &&
              "${my_spinup_restopt}" != "ndays" ]]; then
        echo "ERROR: my_spinup_restopt must be one of nyears, nmonths, ndays; got '${my_spinup_restopt}'"
        fail=1
      fi
      if ! [[ "${my_spinup_restn}" =~ ^[0-9]+$ ]] || (( my_spinup_restn < 1 )); then
        echo "ERROR: my_spinup_restn must be a positive integer; got '${my_spinup_restn}'"
        fail=1
      fi

      if (( fail == 0 )); then
        case "${my_runopt}:${my_spinup_runopt}" in
          nyears:nyears)
            total_common="${my_runn}"
            spinup_common="${my_spinup_runn}"
            ;;
          nmonths:nmonths)
            total_common="${my_runn}"
            spinup_common="${my_spinup_runn}"
            ;;
          ndays:ndays)
            total_common="${my_runn}"
            spinup_common="${my_spinup_runn}"
            ;;
          nyears:nmonths)
            total_common=$(( my_runn * 12 ))
            spinup_common="${my_spinup_runn}"
            ;;
          nmonths:nyears)
            total_common="${my_runn}"
            spinup_common=$(( my_spinup_runn * 12 ))
            ;;
          *)
            echo "ERROR: Unsupported my_runopt/my_spinup_runopt combination for Phase 2 adjustment"
            echo "  total run : ${my_runn} ${my_runopt}"
            echo "  spin-up   : ${my_spinup_runn} ${my_spinup_runopt}"
            fail=1
            ;;
        esac
        if (( fail == 0 && total_common - spinup_common <= 0 )); then
          echo "ERROR: Phase 2 remaining run length must be positive"
          echo "  total run converted : ${total_common}"
          echo "  spin-up converted   : ${spinup_common}"
          fail=1
        fi

        case "${my_restopt}:${my_spinup_runopt}" in
          nyears:nyears)
            rest_common="${my_restn}"
            spinup_common="${my_spinup_runn}"
            ;;
          nmonths:nmonths)
            rest_common="${my_restn}"
            spinup_common="${my_spinup_runn}"
            ;;
          ndays:ndays)
            rest_common="${my_restn}"
            spinup_common="${my_spinup_runn}"
            ;;
          nyears:nmonths)
            rest_common=$(( my_restn * 12 ))
            spinup_common="${my_spinup_runn}"
            ;;
          nmonths:nyears)
            rest_common="${my_restn}"
            spinup_common=$(( my_spinup_runn * 12 ))
            ;;
          *)
            echo "ERROR: Unsupported my_restopt/my_spinup_runopt combination for Phase 2 restart adjustment"
            echo "  restart interval : ${my_restn} ${my_restopt}"
            echo "  spin-up          : ${my_spinup_runn} ${my_spinup_runopt}"
            fail=1
            ;;
        esac
        if (( fail == 0 && rest_common - spinup_common <= 0 )); then
          echo "ERROR: Phase 2 restart interval must be positive"
          echo "  restart interval converted : ${rest_common}"
          echo "  spin-up converted          : ${spinup_common}"
          fail=1
        fi
      fi
    fi
  fi

  if [[ "${do_continue_setup,,}" == "true" ]]; then
    if (( ${#my_restart_ymd[@]} != ${#my_leadymd[@]} )); then
      echo "ERROR: my_restart_ymd must have the same number of entries as my_leadymd"
      echo "  my_restart_ymd: ${#my_restart_ymd[@]}"
      echo "  my_leadymd:     ${#my_leadymd[@]}"
      fail=1
    else
      for rest_date in "${my_restart_ymd[@]}"; do
        if [[ -z "${rest_date}" ]]; then
          echo "ERROR: my_restart_ymd contains an empty entry — all entries must be set to a valid YYYY-MM-DD date"
          fail=1
          break
        fi
      done
    fi
    case "${my_continue_runopt}" in
      nyears|nmonths|ndays) ;;
      "")
        echo "ERROR: my_continue_runopt must be set when do_continue_setup=true"
        fail=1
        ;;
      *)
        echo "ERROR: my_continue_runopt must be one of nyears, nmonths, ndays; got '${my_continue_runopt}'"
        fail=1
        ;;
    esac
    if ! [[ "${my_continue_runn}" =~ ^[0-9]+$ ]] || (( my_continue_runn < 1 )); then
      echo "ERROR: my_continue_runn must be a positive integer; got '${my_continue_runn}'"
      fail=1
    fi
    case "${my_continue_restopt}" in
      nyears|nmonths|ndays) ;;
      "")
        echo "ERROR: my_continue_restopt must be set when do_continue_setup=true"
        fail=1
        ;;
      *)
        echo "ERROR: my_continue_restopt must be one of nyears, nmonths, ndays; got '${my_continue_restopt}'"
        fail=1
        ;;
    esac
    if ! [[ "${my_continue_restn}" =~ ^[0-9]+$ ]] || (( my_continue_restn < 1 )); then
      echo "ERROR: my_continue_restn must be a positive integer; got '${my_continue_restn}'"
      fail=1
    fi
    if [[ -z "${my_continue_walltime}" ]]; then
      echo "ERROR: my_continue_walltime must be set when do_continue_setup=true"
      fail=1
    fi
    if [[ -z "${my_restart_tod}" ]]; then
      echo "ERROR: my_restart_tod must be set when do_continue_setup=true"
      fail=1
    elif ! [[ "${my_restart_tod}" =~ ^[0-9]+$ ]]; then
      echo "ERROR: my_restart_tod must be numeric seconds from start of day; got '${my_restart_tod}'"
      fail=1
    fi
  fi

  if [[ "${do_ensemble_setup,,}" == "true" ]]; then
    for lead_date in "${my_leadymd[@]}"; do
      year=$((10#${lead_date:0:4}))
      if (( year > 2014 )); then
        need_ssp=true
      fi
    done
  fi

  if [[ "${do_continue_setup,,}" == "true" ]]; then
    for rest_date in "${my_restart_ymd[@]}"; do
      [[ -z "${rest_date}" ]] && continue
      year=$((10#${rest_date:0:4}))
      if (( year > 2014 )); then
        need_ssp=true
      fi
    done
  fi

  if [[ "${need_ssp}" == "true" ]]; then
    for required_file in \
      "${my_workdir}/ssp245_user_nl_eam.txt" \
      "${my_workdir}/ssp245_user_nl_elm.txt"; do
      if [[ ! -f "${required_file}" ]]; then
        echo "ERROR: Missing SSP245 namelist include: ${required_file}"
        fail=1
      fi
    done
  fi

  if [[ "${do_ensemble_run,,}" == "true" ]]; then
    for lead_date in "${my_leadymd[@]}"; do
      time_tag="${lead_date//-/}${my_leadtod:0:2}"
      run_path="${my_runpath}/${my_leadcase}_${time_tag}"
      for (( member_idx=0; member_idx<my_ensnum; member_idx++ )); do
        ensstr=EN$(printf "%02d" "${member_idx}")
        case_root="${run_path}/${ensstr}/case_scripts"
        if [[ ! -d "${case_root}" ]]; then
          echo "ERROR: Missing case_scripts directory for bundled run:"
          echo "  ${case_root}"
          fail=1
        fi
      done
    done
  fi

  if [[ "${do_continue_setup,,}" == "true" ]]; then
    for lead_date in "${my_leadymd[@]}"; do
      time_tag="${lead_date//-/}${my_leadtod:0:2}"
      run_path="${my_runpath}/${my_leadcase}_${time_tag}"
      for (( member_idx=0; member_idx<my_ensnum; member_idx++ )); do
        ensstr=EN$(printf "%02d" "${member_idx}")
        case_root="${run_path}/${ensstr}/case_scripts"
        if [[ ! -d "${case_root}" ]]; then
          echo "ERROR: Missing case_scripts directory for continue-run setup:"
          echo "  ${case_root}"
          fail=1
        fi
      done
    done
  fi

  if [[ "${do_e3sm_compile,,}" != "true" &&
        "${do_ensemble_setup,,}" != "true" &&
        "${do_continue_setup,,}" != "true" &&
        "${do_ensemble_run,,}" != "true" ]]; then
    echo "ERROR: No workflow mode is enabled"
    fail=1
  fi

  if (( fail != 0 )); then
    echo "ERROR: Preflight checks failed"
    exit 1
  fi

  echo "Preflight checks passed"
}

preflight_check

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
  template_script="${my_workdir}/compile_and_setup_e3sm_${my_resolution}.${my_compset}.sh"

  # only compile once for first member in first lead_date within this job 
  compiled_once=false
  my_modelexe=""

  for lead_date in "${my_leadymd[@]}"; do
    # Start from an untouched template for every lead date. Otherwise the first
    # date replacement removes the placeholders needed by subsequent dates.
    cp -rp "${template_script}" "${run_script}"
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

    for (( member_idx=0; member_idx<my_ensnum; member_idx++ )); do
      echo === Compile member ${member_idx} ===
      ENSTR=EN$(printf "%02d" "$member_idx")
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

      if [[ "$compiled_once" == false && $member_idx -eq 0 ]]; then
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
  local file="user_nl_mpasi"
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
  local tmp

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

  tmp="$(mktemp "${f}.tmp.XXXXXX")" || { echo "ERROR: failed to create temp file for ${f}"; return 1; }
  cp -p "${f}" "${tmp}" || { echo "ERROR: failed to initialize temp file for ${f}"; rm -f "${tmp}"; return 1; }

  awk '
    function print_highfreq_vars() {
      print "";
      print "    <var name=\"penetrativeTemperatureFlux\"/>";
      print "    <var name=\"latentHeatFlux\"/>";
      print "    <var name=\"sensibleHeatFlux\"/>";
      print "    <var name=\"longWaveHeatFluxUp\"/>";
      print "    <var name=\"longWaveHeatFluxDown\"/>";
      print "    <var name=\"seaIceHeatFlux\"/>";
      print "    <var name=\"shortWaveHeatFlux\"/>";
      print "    <var name=\"evaporationFlux\"/>";
      print "    <var name=\"seaIceSalinityFlux\"/>";
      print "    <var name=\"seaIceFreshWaterFlux\"/>";
      print "    <var name=\"riverRunoffFlux\"/>";
      print "    <var name=\"iceRunoffFlux\"/>";
      print "    <var name=\"rainFlux\"/>";
      print "    <var name=\"snowFlux\"/>";
      print "    <var name=\"bottomLayerShortwaveTemperatureFlux\"/>";
      print "    <var name=\"frazilIceFreshwaterFlux\"/>";
      print "    <var name=\"surfaceBuoyancyForcing\"/>";
      print "    <var name=\"xtime\"/>";
      print "    <var name=\"daysSinceStartOfSim\"/>";
      print "    <var_array name=\"activeTracersAtSurface\"/>";
      print "    <var_array name=\"activeTracersAt250m\"/>";
      print "    <var_array name=\"activeTracersAtBottom\"/>";
      print "    <var name=\"kineticEnergyAtSurface\"/>";
      print "    <var name=\"kineticEnergyAt250m\"/>";
      print "    <var name=\"relativeVorticityAt250m\"/>";
      print "    <var name=\"ssh\"/>";
      print "    <var name=\"boundaryLayerDepth\"/>";
      print "    <var name=\"dThreshMLD\"/>";
      print "    <var name=\"tThreshMLD\"/>";
      print "    <var name=\"barotropicSpeed\"/>";
      print "    <var name=\"windStressMeridional\"/>";
      print "    <var name=\"windStressZonal\"/>";
      print "    <var name=\"atmosphericPressure\"/>";
      print "    <var_struct name=\"tracersSurfaceFlux\"/>";
      print "    <var_array name=\"totalFreshWaterTemperatureFlux\"/>";
      print "    <var name=\"oceanHeatContentSfcToBot\"/>";
      print "    <var name=\"oceanHeatContentSfcTo700m\"/>";
      print "    <var name=\"oceanHeatContent700mTo2000m\"/>";
      print "    <var name=\"oceanHeatContent2000mToBot\"/>";
    }
    function print_monthly_vars() {
      print "";
      print "    <var name=\"daysSinceStartOfSim\"/>";
      print "    <var name=\"binBoundaryMerHeatTrans\"/>";
      print "    <var name=\"binBoundaryZonalMean\"/>";
      print "    <var name=\"ssh\"/>";
      print "    <var_struct name=\"tracers\"/>";
      print "    <var name=\"velocityMeridional\"/>";
      print "    <var name=\"velocityZonal\"/>";
      print "    <var name=\"layerThickness\"/>";
      print "    <var name=\"windStressZonal\"/>";
      print "    <var name=\"windStressMeridional\"/>";
      print "    <var_array name=\"avgValueWithinOceanRegion\"/>";
      print "    <var_array name=\"avgValueWithinOceanLayerRegion\"/>";
      print "    <var_array name=\"avgValueWithinOceanVolumeRegion\"/>";
      print "    <var name=\"meridionalHeatTransportLatZ\"/>";
      print "    <var name=\"meridionalHeatTransportLat\"/>";
      print "    <var name=\"tThreshMLD\"/>";
      print "    <var name=\"dThreshMLD\"/>";
      print "    <var name=\"atmosphericPressure\"/>";
      print "    <var name=\"mocStreamvalLatAndDepthGM\"/>";
      print "    <var name=\"mocStreamvalLatAndDepthRegionGM\"/>";
      print "    <var name=\"mocStreamvalLatAndDepthMLE\"/>";
      print "    <var name=\"mocStreamvalLatAndDepthRegionMLE\"/>";
      print "    <var_struct name=\"tracersSurfaceFlux\"/>";
      print "    <var name=\"penetrativeTemperatureFlux\"/>";
      print "    <var name=\"latentHeatFlux\"/>";
      print "    <var name=\"sensibleHeatFlux\"/>";
      print "    <var name=\"longWaveHeatFluxUp\"/>";
      print "    <var name=\"longWaveHeatFluxDown\"/>";
      print "    <var name=\"seaIceHeatFlux\"/>";
      print "    <var name=\"shortWaveHeatFlux\"/>";
      print "    <var name=\"evaporationFlux\"/>";
      print "    <var name=\"seaIceSalinityFlux\"/>";
      print "    <var name=\"seaIceFreshWaterFlux\"/>";
      print "    <var name=\"riverRunoffFlux\"/>";
      print "    <var name=\"iceRunoffFlux\"/>";
      print "    <var name=\"rainFlux\"/>";
      print "    <var name=\"snowFlux\"/>";
      print "    <var name=\"bottomLayerShortwaveTemperatureFlux\"/>";
      print "    <var name=\"vertDiffTopOfCell\"/>";
      print "    <var name=\"vertViscTopOfCell\"/>";
      print "    <var name=\"boundaryLayerDepth\"/>";
      print "    <var name=\"frazilIceFreshwaterFlux\"/>";
      print "    <var name=\"mocStreamvalLatAndDepth\"/>";
      print "    <var name=\"mocStreamvalLatAndDepthRegion\"/>";
      print "    <var name=\"binBoundaryMocStreamfunction\"/>";
      print "    <var name=\"surfaceBuoyancyForcing\"/>";
      print "    <var name=\"SSHSquared\"/>";
      print "    <var_array name=\"totalFreshWaterTemperatureFlux\"/>";
      print "    <var name=\"oceanHeatContentSfcToBot\"/>";
      print "    <var name=\"oceanHeatContentSfcTo700m\"/>";
      print "    <var name=\"oceanHeatContent700mTo2000m\"/>";
      print "    <var name=\"oceanHeatContent2000mToBot\"/>";
    }
    /<stream name="/ {
      stream="";
      if ($0 ~ /<stream name="highFrequencyOutput"/) {
        stream="high";
        saw_high=1;
      } else if ($0 ~ /<stream name="timeSeriesStatsMonthlyOutput"/) {
        stream="monthly";
        saw_monthly=1;
      } else if ($0 ~ /<stream name="timeSeriesStatsMonthlyMaxOutput"/) {
        stream="max";
        saw_max=1;
      } else if ($0 ~ /<stream name="timeSeriesStatsMonthlyMinOutput"/) {
        stream="min";
        saw_min=1;
      }
    }
    skip && /<\/stream>/ {
      print "</stream>";
      skip=0;
      stream="";
      next;
    }
    skip { next; }
    stream == "high" && /output_interval=/ {
      print "        output_interval=\"00-00-01_00:00:00\"";
      next;
    }
    stream == "high" && /packages="highFrequencyOutputAMPKG">/ {
      print;
      print_highfreq_vars();
      updated_high=1;
      skip=1;
      next;
    }
    stream == "monthly" && /runtime_format="single_file">/ {
      print;
      print_monthly_vars();
      updated_monthly=1;
      skip=1;
      next;
    }
    (stream == "max" || stream == "min") && /<var_array name="activeTracersAtSurface"\/>/ {
      if (stream == "max") updated_max=1;
      if (stream == "min") updated_min=1;
    }
    (stream == "max" || stream == "min") && /<var name="activeTracerForcingMLTend"\/>/ {
      print "    <var_array name=\"activeTracersAtSurface\"/>";
      if (stream == "max") updated_max=1;
      if (stream == "min") updated_min=1;
      skip=1;
      next;
    }
    { print; }
    END {
      if (!saw_high || !saw_monthly || !saw_max || !saw_min ||
          !updated_high || !updated_monthly || !updated_max || !updated_min) {
        print "ERROR: streams.ocean did not match expected MPAS-O stream layout" > "/dev/stderr";
        exit 1;
      }
    }
  ' "${f}" > "${tmp}" || {
    echo "ERROR: failed to rewrite ${f}"
    rm -f "${tmp}"
    return 1
  }

  mv "${tmp}" "${f}" || { echo "ERROR: failed to update ${f}"; rm -f "${tmp}"; return 1; }
  echo "streams.ocean updated."

  # copy to SourceMods
  cp -p "${f}" "${sm_dir}/" || { echo "ERROR: failed to copy streams.ocean to ${sm_dir}"; return 1; }

}

patch_mpassi_streams() {
  echo
  echo "Modifying MPAS (SEAICE) streams files"

  local run_dir="$1"
  local case_dir="$2"
  local f="${run_dir}/streams.seaice"
  local sm_dir="${case_dir}/SourceMods/src.mpassi"
  local tmp

  if [[ -z "${run_dir}" || -z "${case_dir}" ]]; then
    echo "ERROR: patch_mpassi_streams requires: <RUN_DIR> <CASE_DIR>"
    return 1
  fi

  if [[ ! -f "${f}" ]]; then
    echo "ERROR: Missing MPAS-SI streams file:"
    echo "  ${f}"
    return 1
  fi

  mkdir -p "${sm_dir}" || { echo "ERROR: cannot create ${sm_dir}"; return 1; }

  tmp="$(mktemp "${f}.tmp.XXXXXX")" || { echo "ERROR: failed to create temp file for ${f}"; return 1; }
  cp -p "${f}" "${tmp}" || { echo "ERROR: failed to initialize temp file for ${f}"; rm -f "${tmp}"; return 1; }

  awk '
    function print_daily_vars() {
      print "";
      print "    <var name=\"iceAreaCell\"/>";
      print "    <var name=\"iceVolumeCell\"/>";
      print "    <var name=\"snowVolumeCell\"/>";
      print "    <var name=\"uVelocityGeo\"/>";
      print "    <var name=\"vVelocityGeo\"/>";
      print "    <var name=\"iceAreaCategory\"/>";
      print "    <var name=\"iceVolumeCategory\"/>";
      print "    <var name=\"snowVolumeCategory\"/>";
      print "    <var name=\"seaSurfaceTemperature\"/>";
      print "    <var name=\"seaSurfaceSalinity\"/>";
      print "    <var name=\"pondAreaCell\"/>";
      print "    <var name=\"broadbandAlbedo\"/>";
      print "    <var name=\"congelation\"/>";
      print "    <var name=\"frazilFormation\"/>";
      print "    <var name=\"snowiceFormation\"/>";
      print "    <var name=\"snowMelt\"/>";
      print "    <var name=\"surfaceIceMelt\"/>";
      print "    <var name=\"basalIceMelt\"/>";
      print "    <var name=\"lateralIceMelt\"/>";
    }
    /<stream name="/ {
      stream="";
      if ($0 ~ /<stream name="timeSeriesStatsDailyOutput"/) {
        stream="daily";
        saw_daily=1;
      }
    }
    skip && /<\/stream>/ {
      print "</stream>";
      skip=0;
      stream="";
      next;
    }
    skip { next; }
    stream == "daily" && /packages="timeSeriesStatsDailyAMPKG">/ {
      print;
      print_daily_vars();
      updated_daily=1;
      skip=1;
      next;
    }
    { print; }
    END {
      if (!saw_daily || !updated_daily) {
        print "ERROR: streams.seaice did not match expected MPAS-SI stream layout" > "/dev/stderr";
        exit 1;
      }
    }
  ' "${f}" > "${tmp}" || {
    echo "ERROR: failed to rewrite ${f}"
    rm -f "${tmp}"
    return 1
  }

  mv "${tmp}" "${f}" || { echo "ERROR: failed to update ${f}"; rm -f "${tmp}"; return 1; }
  echo "streams.seaice updated."

  # copy to SourceMods
  cp -p "${f}" "${sm_dir}/" || { echo "ERROR: failed to copy streams.seaice to ${sm_dir}"; return 1; }

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

  ./xmlchange run_exe="--exclusive --kill-on-bad-exit=1 --job-name=${CASE_NAME} \${EXEROOT}/e3sm.exe " || return 1
  ./xmlchange RUN_TYPE="${my_runtype,,}" || return 1
  ./xmlchange GET_REFCASE=TRUE || return 1
  ./xmlchange RUN_REFDIR="${TGT_ROOT}" || return 1
  ./xmlchange RUN_REFCASE="${my_refcase2}" || return 1
  ./xmlchange RUN_REFDATE="${START_DATE}" || return 1
  ./xmlchange RUN_REFTOD="${START_TOD}" || return 1
  if [[ "${do_short_spinup,,}" == "true" ]]; then
    ./xmlchange DOUT_S=FALSE || return 1
  else
    ./xmlchange DOUT_S="${my_short_archive,,}" || return 1
  fi
  ./xmlchange DOUT_S_ROOT="${ARCHIVE_DIR}" || return 1
  ./xmlchange RUN_STARTDATE="${START_DATE}" || return 1
  ./xmlchange START_TOD="${START_TOD}" || return 1
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

  if [[ "${do_short_spinup,,}" == "true" ]]; then
    # use spinup-specific time steps and run length for the initial short spinup phase
    if [[ -n "${my_spinup_mpaso_dt}"  ]]; then user_mpaso_nl  "${my_spinup_mpaso_dt}"  || return 1; fi
    if [[ -n "${my_spinup_mpassi_dt}" ]]; then user_mpassi_nl "${my_spinup_mpassi_dt}" || return 1; fi
    ./xmlchange REST_OPTION="${my_spinup_restopt}" || return 1
    ./xmlchange REST_N="${my_spinup_restn}" || return 1
    ./xmlchange STOP_OPTION="${my_spinup_runopt}" || return 1
    ./xmlchange STOP_N="${my_spinup_runn}" || return 1
  else
    user_mpaso_nl  "${my_mpaso_dt:-00:30:00}"  || return 1
    user_mpassi_nl "${my_mpassi_dt:-1800}"      || return 1
    ./xmlchange REST_OPTION="${my_restopt}" || return 1
    ./xmlchange REST_N="${my_restn}" || return 1
    ./xmlchange STOP_OPTION="${my_runopt}" || return 1
    ./xmlchange STOP_N="${my_runn}" || return 1
  fi 
  
  ./case.setup || return 1
  ./xmlchange BUILD_COMPLETE=TRUE || return 1

  patch_mpaso_streams  "${RUN_ROOT}" "${CASE_ROOT}" || return 1
  patch_mpassi_streams "${RUN_ROOT}" "${CASE_ROOT}" || return 1

  echo "=== Finished ${CASE_NAME} ==="

}

# ---------------------------------------------------------------------------
# Compute Phase 2 remaining length after spin-up.
#
# STOP_OPTION/STOP_N describe the target stop interval from RUN_STARTDATE in
# this continuation workflow. After a short spin-up, Phase 2 must use the
# remaining interval so the model reaches the original total target date.
# ---------------------------------------------------------------------------
compute_phase2_remaining_length() {
  local total_opt="$1"
  local total_n="$2"
  local spinup_opt="$3"
  local spinup_n="$4"

  local common_opt
  local total_common
  local spinup_common

  case "${total_opt}:${spinup_opt}" in
    nyears:nyears)
      common_opt="nyears"
      total_common="${total_n}"
      spinup_common="${spinup_n}"
      ;;
    nmonths:nmonths)
      common_opt="nmonths"
      total_common="${total_n}"
      spinup_common="${spinup_n}"
      ;;
    ndays:ndays)
      common_opt="ndays"
      total_common="${total_n}"
      spinup_common="${spinup_n}"
      ;;
    nyears:nmonths)
      common_opt="nmonths"
      total_common=$(( total_n * 12 ))
      spinup_common="${spinup_n}"
      ;;
    nmonths:nyears)
      common_opt="nmonths"
      total_common="${total_n}"
      spinup_common=$(( spinup_n * 12 ))
      ;;
    *)
      echo "ERROR: Unsupported run-option combination." >&2
      echo "  total run : ${total_n} ${total_opt}" >&2
      echo "  spin-up   : ${spinup_n} ${spinup_opt}" >&2
      echo "Mixed ndays with nmonths/nyears is not supported because it is calendar-dependent." >&2
      return 1
      ;;
  esac

  phase2_runopt="${common_opt}"
  phase2_runn=$(( total_common - spinup_common ))

  if (( phase2_runn <= 0 )); then
    echo "ERROR: Phase 2 remaining run length must be positive." >&2
    echo "  total run converted : ${total_common} ${common_opt}" >&2
    echo "  spin-up converted   : ${spinup_common} ${common_opt}" >&2
    echo "  remaining Phase 2   : ${phase2_runn} ${phase2_runopt}" >&2
    return 1
  fi

  echo "Total run converted to common unit : ${total_common} ${common_opt}"
  echo "Spin-up converted to common unit   : ${spinup_common} ${common_opt}"
  echo "Phase 2 remaining length           : ${phase2_runn} ${phase2_runopt}"
}

# ---------------------------------------------------------------------------
# Compute Phase 2 restart interval after spin-up.
#
# REST_OPTION/REST_N set the restart-write interval. For Phase 2, subtract the
# completed spin-up so restart writes still align with the original schedule.
# ---------------------------------------------------------------------------
compute_phase2_restart_interval() {
  local rest_opt="$1"
  local rest_n="$2"
  local spinup_opt="$3"
  local spinup_n="$4"

  local common_opt
  local rest_common
  local spinup_common

  case "${rest_opt}:${spinup_opt}" in
    nyears:nyears)
      common_opt="nyears"
      rest_common="${rest_n}"
      spinup_common="${spinup_n}"
      ;;
    nmonths:nmonths)
      common_opt="nmonths"
      rest_common="${rest_n}"
      spinup_common="${spinup_n}"
      ;;
    ndays:ndays)
      common_opt="ndays"
      rest_common="${rest_n}"
      spinup_common="${spinup_n}"
      ;;
    nyears:nmonths)
      common_opt="nmonths"
      rest_common=$(( rest_n * 12 ))
      spinup_common="${spinup_n}"
      ;;
    nmonths:nyears)
      common_opt="nmonths"
      rest_common="${rest_n}"
      spinup_common=$(( spinup_n * 12 ))
      ;;
    *)
      echo "ERROR: Unsupported REST/spin-up option combination." >&2
      echo "  restart interval : ${rest_n} ${rest_opt}" >&2
      echo "  spin-up length   : ${spinup_n} ${spinup_opt}" >&2
      echo "Mixed ndays with nmonths/nyears is not supported because it is calendar-dependent." >&2
      return 1
      ;;
  esac

  phase2_restopt="${common_opt}"
  phase2_restn=$(( rest_common - spinup_common ))

  if (( phase2_restn <= 0 )); then
    echo "ERROR: Phase 2 restart interval must be positive." >&2
    echo "  restart interval converted : ${rest_common} ${common_opt}" >&2
    echo "  spin-up converted          : ${spinup_common} ${common_opt}" >&2
    echo "  adjusted restart interval  : ${phase2_restn} ${phase2_restopt}" >&2
    return 1
  fi

  echo "Restart interval converted to common unit : ${rest_common} ${common_opt}"
  echo "Spin-up converted to common unit          : ${spinup_common} ${common_opt}"
  echo "Phase 2 restart interval                  : ${phase2_restn} ${phase2_restopt}"
}

# -----------------------------------------------
# Check whether a member already completed the
# short spinup restart needed by Phase 2.
# -----------------------------------------------
spinup_end_date() {
  local start_date="$1"
  local spinup_opt="$2"
  local spinup_n="$3"

  case "${spinup_opt}" in
    nyears)
      date -u -d "${start_date} +${spinup_n} year" +%F
      ;;
    nmonths)
      date -u -d "${start_date} +${spinup_n} month" +%F
      ;;
    ndays)
      date -u -d "${start_date} +${spinup_n} day" +%F
      ;;
    *)
      echo "ERROR: Unsupported spin-up run option: ${spinup_opt}" >&2
      return 1
      ;;
  esac
}

has_completed_spinup_restart() {
  local lead_date="$1"
  local i="$2"
  local time_tag run_path enstr case_name run_root end_date end_stamp end_stamp_mpas
  local required_file

  time_tag="${lead_date//-/}${my_leadtod:0:2}"
  run_path="${my_runpath}/${my_leadcase}_${time_tag}"
  enstr="EN$(printf "%02d" "$i")"
  case_name="${my_leadcase}_${time_tag}.${enstr}"
  run_root="${run_path}/${enstr}/run"

  end_date="$(spinup_end_date "${lead_date}" "${my_spinup_runopt}" "${my_spinup_runn}")" || return 1
  end_stamp="${end_date}-${my_leadtod}"
  end_stamp_mpas="${end_date}_${my_leadtod}"

  for required_file in \
    "${run_root}/${case_name}.cpl.r.${end_stamp}.nc" \
    "${run_root}/${case_name}.eam.r.${end_stamp}.nc" \
    "${run_root}/${case_name}.elm.r.${end_stamp}.nc" \
    "${run_root}/${case_name}.mosart.r.${end_stamp}.nc" \
    "${run_root}/${case_name}.mpaso.rst.${end_stamp_mpas}.nc" \
    "${run_root}/${case_name}.mpassi.rst.${end_stamp_mpas}.nc"; do
    [[ -f "${required_file}" ]] || return 1
  done

  return 0
}

# -----------------------------------------------
# reconfigure one member for the main run after
# the short spinup phase has completed
# -----------------------------------------------
configure_main_run_member() {
  local lead_date="$1"
  local i="$2"
  local TIME_TAG RUN_PATH ENSTR CASE_NAME CASE_ROOT

  TIME_TAG="${lead_date//-/}${my_leadtod:0:2}"
  RUN_PATH="${my_runpath}/${my_leadcase}_${TIME_TAG}"
  ENSTR="EN$(printf "%02d" "$i")"
  CASE_NAME="${my_leadcase}_${TIME_TAG}.${ENSTR}"
  CASE_ROOT="${RUN_PATH}/${ENSTR}/case_scripts"

  echo "=== Configuring main run for ${CASE_NAME} ==="

  if ! has_completed_spinup_restart "${lead_date}" "${i}"; then
    echo "ERROR: Cannot configure Phase 2; missing complete Phase 1 restart for ${CASE_NAME}"
    return 1
  fi

  cd "${CASE_ROOT}" || {
    echo "ERROR: Cannot cd to ${CASE_ROOT}"
    return 1
  }

  compute_phase2_remaining_length \
    "${my_runopt}" "${my_runn}" \
    "${my_spinup_runopt}" "${my_spinup_runn}" || return 1

  compute_phase2_restart_interval \
    "${my_restopt}" "${my_restn}" \
    "${my_spinup_runopt}" "${my_spinup_runn}" || return 1

  ./xmlchange CONTINUE_RUN=TRUE              || return 1
  ensure_continue_run_true "${CASE_NAME}"    || return 1
  ./xmlchange DOUT_S="${my_short_archive,,}"   || return 1
  ./xmlchange REST_OPTION="${phase2_restopt}" || return 1
  ./xmlchange REST_N="${phase2_restn}"        || return 1
  ./xmlchange STOP_OPTION="${phase2_runopt}"  || return 1
  ./xmlchange STOP_N="${phase2_runn}"         || return 1
  ensure_xml_value "${CASE_NAME}" "REST_OPTION" "${phase2_restopt}" || return 1
  ensure_xml_value "${CASE_NAME}" "REST_N" "${phase2_restn}" || return 1
  ensure_xml_value "${CASE_NAME}" "STOP_OPTION" "${phase2_runopt}" || return 1
  ensure_xml_value "${CASE_NAME}" "STOP_N" "${phase2_runn}" || return 1

  # switch to standard time steps for the main run
  user_mpaso_nl  "${my_mpaso_dt:-00:30:00}"  || return 1
  user_mpassi_nl "${my_mpassi_dt:-1800}" || return 1

  echo "=== Main run configured for ${CASE_NAME} ==="
}

submit_main_run_member() {
  local lead_date="$1"
  local i="$2"
  local TIME_TAG RUN_PATH ENSTR CASE_NAME CASE_ROOT LOG_FILE

  TIME_TAG="${lead_date//-/}${my_leadtod:0:2}"
  RUN_PATH="${my_runpath}/${my_leadcase}_${TIME_TAG}"
  ENSTR="EN$(printf "%02d" "$i")"
  CASE_NAME="${my_leadcase}_${TIME_TAG}.${ENSTR}"
  CASE_ROOT="${RUN_PATH}/${ENSTR}/case_scripts"
  LOG_FILE="e3sm_${ENSTR}.log.o${SLURM_JOB_ID}"

  echo "=== Run member ${i}: ${CASE_NAME} ==="

  cd "${CASE_ROOT}" || {
    echo "ERROR: Cannot cd to ${CASE_ROOT}"
    return 1
  }

  ./case.submit --no-batch > "${LOG_FILE}" 2>&1 || {
    echo "ERROR: main run failed for ${CASE_NAME}"
    echo "  log=${CASE_ROOT}/${LOG_FILE}"
    return 1
  }

  echo "=== Main run complete for ${CASE_NAME} ==="
}

ensure_continue_run_true() {
  local case_name="$1"
  local continue_run

  continue_run=$(./xmlquery --value CONTINUE_RUN) || {
    echo "ERROR: Failed to query CONTINUE_RUN for ${case_name}"
    return 1
  }

  if [[ "${continue_run^^}" != "TRUE" ]]; then
    echo "ERROR: CONTINUE_RUN is not TRUE for restart case ${case_name}"
    echo "  CONTINUE_RUN=${continue_run}"
    return 1
  fi
}

ensure_xml_value() {
  local case_name="$1"
  local xml_id="$2"
  local expected="$3"
  local actual

  actual=$(./xmlquery --value "${xml_id}") || {
    echo "ERROR: Failed to query ${xml_id} for ${case_name}"
    return 1
  }

  if [[ "${actual}" != "${expected}" ]]; then
    echo "ERROR: ${xml_id} does not match expected value for ${case_name}"
    echo "  expected: ${expected}"
    echo "  actual:   ${actual}"
    return 1
  fi
}

run_short_spinup_then_main_member() {
  local lead_date="$1"
  local i="$2"
  local TIME_TAG RUN_PATH ENSTR CASE_NAME CASE_ROOT LOG_FILE

  TIME_TAG="${lead_date//-/}${my_leadtod:0:2}"
  RUN_PATH="${my_runpath}/${my_leadcase}_${TIME_TAG}"
  ENSTR="EN$(printf "%02d" "$i")"
  CASE_NAME="${my_leadcase}_${TIME_TAG}.${ENSTR}"
  CASE_ROOT="${RUN_PATH}/${ENSTR}/case_scripts"
  LOG_FILE="e3sm_spinup_${ENSTR}.log.o${SLURM_JOB_ID}"

  echo "=== Spinup/main chain for member ${i}: ${CASE_NAME} ==="

  if [[ "${do_skip_completed_spinup,,}" == "true" ]] &&
     has_completed_spinup_restart "${lead_date}" "${i}"; then
    echo "Spinup restart already exists for ${CASE_NAME}; skipping Phase 1 submit"
  else
    cd "${CASE_ROOT}" || {
      echo "ERROR: Cannot cd to ${CASE_ROOT}"
      return 1
    }

    ./case.submit --no-batch > "${LOG_FILE}" 2>&1 || {
      echo "ERROR: spinup run failed for ${CASE_NAME}"
      echo "  log=${CASE_ROOT}/${LOG_FILE}"
      return 1
    }
  fi

  configure_main_run_member "${lead_date}" "${i}" || return 1
  submit_main_run_member "${lead_date}" "${i}" || return 1

  echo "=== Spinup/main chain complete for ${CASE_NAME} ==="
}

restart_files_exist_in_run() {
  local case_name="$1"
  local run_root="$2"
  local rest_stamp="$3"
  local rest_stamp_mpas="$4"
  local required_file

  for required_file in \
    "${run_root}/${case_name}.cpl.r.${rest_stamp}.nc" \
    "${run_root}/${case_name}.eam.r.${rest_stamp}.nc" \
    "${run_root}/${case_name}.elm.r.${rest_stamp}.nc" \
    "${run_root}/${case_name}.mosart.r.${rest_stamp}.nc" \
    "${run_root}/${case_name}.mpaso.rst.${rest_stamp_mpas}.nc" \
    "${run_root}/${case_name}.mpassi.rst.${rest_stamp_mpas}.nc"; do
    [[ -f "${required_file}" ]] || return 1
  done
}

first_rpointer_record() {
  local rpointer_file="$1"
  local line

  while IFS= read -r line || [[ -n "${line}" ]]; do
    line="${line%%#*}"
    line="${line#"${line%%[![:space:]]*}"}"
    line="${line%"${line##*[![:space:]]}"}"
    if [[ -n "${line}" ]]; then
      printf '%s\n' "${line}"
      return 0
    fi
  done < "${rpointer_file}"

  return 1
}

validate_file_rpointer() {
  local run_root="$1"
  local rpointer_name="$2"
  local expected_file="$3"
  local rpointer_file="${run_root}/rpointer.${rpointer_name}"
  local pointer_value pointer_file

  if [[ ! -f "${rpointer_file}" ]]; then
    echo "ERROR: Missing rpointer file:"
    echo "  ${rpointer_file}"
    return 1
  fi

  if ! pointer_value=$(first_rpointer_record "${rpointer_file}"); then
    echo "ERROR: rpointer file has no restart record:"
    echo "  ${rpointer_file}"
    return 1
  fi

  pointer_file="${pointer_value#./}"
  pointer_file="${pointer_file##*/}"

  if [[ "${pointer_file}" != "${expected_file}" ]]; then
    echo "ERROR: rpointer.${rpointer_name} does not reference the expected restart file:"
    echo "  rpointer: ${rpointer_file}"
    echo "  expected: ${expected_file}"
    echo "  found:    ${pointer_value}"
    return 1
  fi

  if [[ ! -f "${run_root}/${pointer_file}" ]]; then
    echo "ERROR: rpointer.${rpointer_name} references a restart file missing from run directory:"
    echo "  rpointer: ${rpointer_file}"
    echo "  missing:  ${run_root}/${pointer_file}"
    return 1
  fi
}

validate_time_rpointer() {
  local run_root="$1"
  local rpointer_name="$2"
  local expected_time="$3"
  local expected_file="$4"
  local rpointer_file="${run_root}/rpointer.${rpointer_name}"
  local pointer_value

  if [[ ! -f "${rpointer_file}" ]]; then
    echo "ERROR: Missing rpointer file:"
    echo "  ${rpointer_file}"
    return 1
  fi

  if ! pointer_value=$(first_rpointer_record "${rpointer_file}"); then
    echo "ERROR: rpointer file has no restart record:"
    echo "  ${rpointer_file}"
    return 1
  fi

  if [[ "${pointer_value}" != "${expected_time}" ]]; then
    echo "ERROR: rpointer.${rpointer_name} does not reference the expected restart time:"
    echo "  rpointer: ${rpointer_file}"
    echo "  expected: ${expected_time}"
    echo "  found:    ${pointer_value}"
    return 1
  fi

  if [[ ! -f "${run_root}/${expected_file}" ]]; then
    echo "ERROR: expected restart file for rpointer.${rpointer_name} is missing from run directory:"
    echo "  ${run_root}/${expected_file}"
    return 1
  fi
}

validate_continue_restart_state() {
  local case_name="$1"
  local run_root="$2"
  local rest_date="$3"
  local rest_tod="$4"
  local rest_stamp="$5"
  local rest_stamp_mpas="$6"
  local rest_seconds=$((10#${rest_tod}))
  local rest_time_mpas

  printf -v rest_time_mpas "%s_%02d:%02d:%02d" \
    "${rest_date}" \
    $((rest_seconds / 3600)) \
    $(((rest_seconds % 3600) / 60)) \
    $((rest_seconds % 60))

  restart_files_exist_in_run \
    "${case_name}" "${run_root}" "${rest_stamp}" "${rest_stamp_mpas}" || {
    echo "ERROR: missing one or more restart files in run directory:"
    echo "  ${run_root}"
    echo "  expected restart tag: ${rest_stamp}"
    echo "  expected MPAS tag:    ${rest_stamp_mpas}"
    return 1
  }

  validate_file_rpointer "${run_root}" "atm" "${case_name}.eam.r.${rest_stamp}.nc" || return 1
  validate_file_rpointer "${run_root}" "lnd" "${case_name}.elm.r.${rest_stamp}.nc" || return 1
  validate_file_rpointer "${run_root}" "rof" "${case_name}.mosart.r.${rest_stamp}.nc" || return 1
  validate_file_rpointer "${run_root}" "drv" "${case_name}.cpl.r.${rest_stamp}.nc" || return 1
  validate_time_rpointer "${run_root}" "ocn" "${rest_time_mpas}" "${case_name}.mpaso.rst.${rest_stamp_mpas}.nc" || return 1
  validate_time_rpointer "${run_root}" "ice" "${rest_time_mpas}" "${case_name}.mpassi.rst.${rest_stamp_mpas}.nc" || return 1

  echo "Restart files and rpointer files are consistent for ${case_name} at ${rest_stamp}"
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
  local rest_stamp rest_stamp_mpas
  local run_restart_ok=true

  TIME_TAG="${lead_date//-/}${my_leadtod:0:2}"
  RUN_PATH="${my_runpath}/${my_leadcase}_${TIME_TAG}"

  echo "=== Setup member ${i} for ${lead_date} ==="

  ENSTR="EN$(printf "%02d" "$i")"
  CASE_NAME="${my_leadcase}_${TIME_TAG}.${ENSTR}"
  CASE_ROOT="${RUN_PATH}/${ENSTR}/case_scripts"
  RUN_ROOT="${RUN_PATH}/${ENSTR}/run"
  ARCHIVE_DIR="${RUN_PATH}/${ENSTR}/archive"
  REST_CASE="${CASE_NAME}"

  rest_stamp="${rest_date}-${rest_tod}"
  rest_stamp_mpas="${rest_date}_${rest_tod}"
  restart_files_exist_in_run \
    "${CASE_NAME}" "${RUN_ROOT}" "${rest_stamp}" "${rest_stamp_mpas}" || run_restart_ok=false

  if [[ "${run_restart_ok}" == "true" ]]; then
    SRC_DIR=""
    echo "Restart files already in run directory; validating rpointer files and skipping archive copy"
    validate_continue_restart_state \
      "${CASE_NAME}" "${RUN_ROOT}" "${rest_date}" "${rest_tod}" "${rest_stamp}" "${rest_stamp_mpas}" || return 1
  else
    SRC_DIR="${ARCHIVE_DIR}/rest/${rest_date}-${rest_tod}"
  fi
  REST_DATE_PREFIX="${rest_date:0:4}-${rest_date:5:2}"   # YYYY-MM

  if [[ -n "${SRC_DIR}" ]]; then
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

    validate_continue_restart_state \
      "${CASE_NAME}" "${RUN_ROOT}" "${rest_date}" "${rest_tod}" "${rest_stamp}" "${rest_stamp_mpas}" || return 1
  fi

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
        echo "INFO: Optional history files not found; skipping:"
        echo "  ${REST_GLOB}"
        continue
      fi

      # sort and get the last file
      last_file=$(printf '%s\n' "${files[@]}" | sort | tail -n 1)

      # extract date part from filename
      last_base=$(basename "$last_file")
      echo "Last file: $last_file"

      # destination path
      dest_file="${RUN_ROOT}/${last_base}"

      # copy only if not already present (or different)
      if [[ ! -f "$dest_file" ]]; then
        cp -p "$last_file" "$dest_file" || {
          echo "ERROR: Failed to copy history file:"
          echo "  ${last_file} -> ${dest_file}"
          return 1
        }
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

  ./xmlchange run_exe="--exclusive --kill-on-bad-exit=1 --job-name=${CASE_NAME} \${EXEROOT}/e3sm.exe " || return 1
  ./xmlchange CONTINUE_RUN=TRUE || return 1
  ensure_continue_run_true "${CASE_NAME}" || return 1
  # CIME still evaluates hybrid reference settings while rebuilding namelists.
  # Keep them aligned with the case's original lead date.
  ./xmlchange RUN_REFDATE="${lead_date}" || return 1
  ./xmlchange RUN_REFTOD="${my_leadtod}" || return 1
  ./xmlchange RUN_STARTDATE="${lead_date}" || return 1
  ./xmlchange START_TOD="${my_leadtod}" || return 1
  ./xmlchange DOUT_S="${my_short_archive,,}" || return 1
  ./xmlchange DOUT_S_ROOT="${ARCHIVE_DIR}" || return 1
  ./xmlchange REST_OPTION="${my_continue_restopt}" || return 1
  ./xmlchange REST_N="${my_continue_restn}" || return 1
  ./xmlchange STOP_OPTION="${my_continue_runopt}" || return 1
  ./xmlchange STOP_N="${my_continue_runn}" || return 1
  ensure_xml_value "${CASE_NAME}" "REST_OPTION" "${my_continue_restopt}" || return 1
  ensure_xml_value "${CASE_NAME}" "REST_N" "${my_continue_restn}" || return 1
  ensure_xml_value "${CASE_NAME}" "STOP_OPTION" "${my_continue_runopt}" || return 1
  ensure_xml_value "${CASE_NAME}" "STOP_N" "${my_continue_runn}" || return 1
  ./xmlchange JOB_WALLCLOCK_TIME="${my_continue_walltime}" || return 1
  ./xmlchange JOB_QUEUE="${my_jobqueue}" || return 1

  if (( year > 2014 )); then
    #change the forcing file to extend simulation beyond 2014
    cat "${my_workdir}/ssp245_user_nl_eam.txt" >> user_nl_eam || return 1
    cat "${my_workdir}/ssp245_user_nl_elm.txt" >> user_nl_elm || return 1
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

  for lead_idx in "${!my_leadymd[@]}"; do
    lead_date="${my_leadymd[$lead_idx]}"

    fail_flag=0
    
    pids=()

    for (( member_idx=0; member_idx<my_ensnum; member_idx++ )); do
      wait_for_slot "${max_setup_jobs}"

      (
        setup_one_member "${lead_date}" "${member_idx}"
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

  for lead_idx in "${!my_leadymd[@]}"; do
    lead_date="${my_leadymd[$lead_idx]}"
    rest_date="${my_restart_ymd[$lead_idx]}"

    fail_flag=0

    pids=()

    for (( member_idx=0; member_idx<my_ensnum; member_idx++ )); do
      wait_for_slot "${max_continue_jobs}"

      (
        continue_one_member "${lead_date}" "${rest_date}" "${member_idx}"
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
  declare -A run_pid_desc=()
  failed_member_jobs=()

  if [[ "${do_short_spinup,,}" == "true" ]]; then
    echo "===== Chained short spinup + main run ====="
    chain_pids=()
    declare -A chain_pid_desc=()
    failed_member_jobs=()

    for lead_date in "${my_leadymd[@]}"; do
      for (( member_idx=0; member_idx<my_ensnum; member_idx++ )); do
        wait_for_slot "${NMAXPS}"
        TIME_TAG="${lead_date//-/}${my_leadtod:0:2}"
        ENSTR=EN$(printf "%02d" "$member_idx")
        CASE_NAME=${my_leadcase}_${TIME_TAG}.${ENSTR}

        (
          run_short_spinup_then_main_member "${lead_date}" "${member_idx}"
        ) &
        chain_pids+=($!)
        chain_pid_desc[$!]="${CASE_NAME}"
        echo ============================
      done
    done

    if (( ${#chain_pids[@]} > 0 )); then
      for pid in "${chain_pids[@]}"; do
        if ! wait "$pid"; then
          echo "ERROR: spinup/main chain failed: ${chain_pid_desc[$pid]}"
          failed_member_jobs+=("${chain_pid_desc[$pid]}")
        fi
      done
      chain_pids=()
    fi

    if (( ${#failed_member_jobs[@]} > 0 )); then
      echo "======================================================="
      echo "WARNING: One or more chained spinup/main member jobs failed"
      printf '  %s\n' "${failed_member_jobs[@]}"
      echo "Continuing so successful ensemble members remain available"
      echo "======================================================="
    fi

    echo "===== All chained short spinup + main runs complete ====="
  else
    # -------------------------------------------------------
    # Phase 2: main ensemble run without short spinup
    # -------------------------------------------------------
    for lead_date in "${my_leadymd[@]}"; do
      TIME_TAG="${lead_date//-/}${my_leadtod:0:2}"
      for (( member_idx=0; member_idx<my_ensnum; member_idx++ )); do
        wait_for_slot "${NMAXPS}"
        ENSTR=EN$(printf "%02d" "$member_idx")
        CASE_NAME=${my_leadcase}_${TIME_TAG}.${ENSTR}

        (
          submit_main_run_member "${lead_date}" "${member_idx}"
        ) &
        pids+=($!)
        run_pid_desc[$!]="${CASE_NAME}"
        echo ============================
      done
    done 

    if (( ${#pids[@]} > 0 )); then
      for pid in "${pids[@]}"; do
        if ! wait "$pid"; then
          echo "ERROR: member job failed: ${run_pid_desc[$pid]}"
          failed_member_jobs+=("${run_pid_desc[$pid]}")
        fi
      done
      pids=()
    fi

    if (( ${#failed_member_jobs[@]} > 0 )); then
      echo "======================================================="
      echo "WARNING: One or more bundled member jobs failed"
      printf '  %s\n' "${failed_member_jobs[@]}"
      echo "Continuing so successful ensemble members remain available"
      echo "======================================================="
    fi
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
      squeue --job=${SLURM_JOB_ID} --steps
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
