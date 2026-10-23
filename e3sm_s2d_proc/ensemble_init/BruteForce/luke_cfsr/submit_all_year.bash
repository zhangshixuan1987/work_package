#!/bin/bash -el 
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --account=gen0010
#SBATCH  --constraint=cpu
#SBATCH  --qos=regular
#SBATCH  --nodes=4
#SBATCH  --time=24:00:00
#SBATCH  --job-name=s2d_ic
#SBATCH  --output=slurm-%j.out

set -o pipefail

# ==============================================================================
# USER SETTINGS
# Edit this block, then run: ./submit_all_year.bash
#
# With no option, the script submits a job from a login node and runs all
# dates sequentially when it is already inside an interactive Slurm allocation.
# Use --run or --submit to override that automatic choice.
# ==============================================================================
# Options: "jra55-fosirl", "bruteforce-cfsr", "reanalysis-cfsr", or "4denvar-ocn"
workflow_mode="bruteforce-cfsr"
script_dir="${PWD}"

# Default forcing and atmosphere selections; the mode mapping may override them.
fosi="JRA55"
atm="ERA5-Tang2"
# Optional sea-ice replacement for reanalysis-cfsr (in development): "" or "NSIDCv6"
ice=""

# Component/output selections associated with each workflow mode.
mode=${IC_MODE:-${workflow_mode}}
case "${mode}" in
  jra55-fosirl)
    group="JRA55-FOSIRL"
    ocn=""
    ;;
  bruteforce-cfsr)
    group="BruteForce"
    atm="ERA5-Tang2"
    ocn="LukeVanRoekel2"
    ;;
  reanalysis-cfsr)
    group="Reanalysis"
    atm="ERA5-Tang2"
    ocn="LukeVanRoekel3"
    ;;
  4denvar-ocn)
    group="4DEnVarOcn"
    atm="4DEnVarOcn"
    ocn="4DEnVarOcn"
    ;;
  *)
    echo "ERROR: workflow_mode must be jra55-fosirl, bruteforce-cfsr, reanalysis-cfsr, or 4denvar-ocn; got '${mode}'" >&2
    exit 2
    ;;
esac

# JRA55/BruteForce example: 1980-2018 with months=(2 8)
# 4DEnVar availability:     1980-2015 with months=(5 11)
# Reanalysis ocean ICs:     2015 with months=(2 5 8 11); 2026 with months=(2 5 8)
#start_year=1980
#end_year=2018
#months=(2 8)

start_year=2011
end_year=2011
months=(2)

# Optional date filters for unavailable inputs. With both lists populated, a
# date is skipped only when both its year and month match. Leave both empty to
# skip nothing.
skip_years=()
skip_months=()

# Existing nonempty output directories are protected unless this is "true".
overwrite_existing=true #false

# Set to "true" to validate inputs and outputs without running or submitting.
check_only=false
# ==============================================================================

usage() {
  echo "Usage: $0 [--auto|--run|--submit]" >&2
  echo "Edit the USER SETTINGS block at the top of the script to configure the run." >&2
  echo "  --auto    submit from a login node; run sequentially in a Slurm allocation (default)" >&2
  echo "  --run     run all configured dates sequentially in the current allocation" >&2
  echo "  --submit  submit one Slurm job that runs all configured dates sequentially" >&2
  exit 2
}

execution_mode=auto
(( $# <= 1 )) || usage
case ${1:---auto} in
  --auto) execution_mode=auto ;;
  --run) execution_mode=run ;;
  --submit) execution_mode=submit ;;
  *) usage ;;
esac
[[ -z ${SLURM_ARRAY_TASK_ID:-} ]] || {
  echo "ERROR: this script uses one normal Slurm job, not a job array." >&2
  exit 2
}

validate_settings() {
  [[ ${mode} == jra55-fosirl || ${mode} == bruteforce-cfsr || ${mode} == reanalysis-cfsr || ${mode} == 4denvar-ocn ]] || usage
  [[ ${start_year} =~ ^[0-9]{4}$ && ${end_year} =~ ^[0-9]{4}$ ]] || {
    echo "ERROR: start_year and end_year must be four-digit years" >&2
    exit 2
  }
  (( start_year <= end_year )) || {
    echo "ERROR: start_year must not be later than end_year" >&2
    exit 2
  }
  (( ${#months[@]} > 0 )) || {
    echo "ERROR: months must contain at least one month" >&2
    exit 2
  }
  local value
  for value in "${months[@]}"; do
    [[ ${value} =~ ^[0-9]+$ ]] && (( 10#${value} >= 1 && 10#${value} <= 12 )) || {
      echo "ERROR: invalid month '${value}'" >&2
      exit 2
    }
  done
  for value in "${skip_years[@]}"; do
    [[ ${value} =~ ^[0-9]{4}$ ]] || {
      echo "ERROR: invalid skip year '${value}'; expected YYYY" >&2
      exit 2
    }
  done
  for value in "${skip_months[@]}"; do
    [[ ${value} =~ ^[0-9]+$ ]] && (( 10#${value} >= 1 && 10#${value} <= 12 )) || {
      echo "ERROR: invalid skip month '${value}'" >&2
      exit 2
    }
  done
  if [[ ${mode} == 4denvar-ocn ]]; then
    (( start_year >= 1980 && end_year <= 2015 )) || {
      echo "ERROR: 4denvar-ocn inputs are available only for 1980-2015" >&2
      exit 2
    }
    for value in "${months[@]}"; do
      (( 10#${value} == 5 || 10#${value} == 11 )) || {
        echo "ERROR: 4denvar-ocn inputs are available only for months 5 and 11" >&2
        exit 2
      }
    done
  fi
  [[ ${overwrite_existing} == true || ${overwrite_existing} == false ]] || {
    echo "ERROR: overwrite_existing must be true or false" >&2
    exit 2
  }
  [[ ${check_only} == true || ${check_only} == false ]] || {
    echo "ERROR: check_only must be true or false" >&2
    exit 2
  }
}

validate_settings
task_years=()
task_months=()
for candidate_year in $(seq "${start_year}" "${end_year}"); do
  for candidate_month in "${months[@]}"; do
    candidate_month_num=$((10#${candidate_month}))
    year_matches=false
    month_matches=false
    for skipped_year in "${skip_years[@]}"; do
      if (( candidate_year == 10#${skipped_year} )); then
        year_matches=true
        break
      fi
    done
    for skipped_month in "${skip_months[@]}"; do
      if (( candidate_month_num == 10#${skipped_month} )); then
        month_matches=true
        break
      fi
    done
    skip_candidate=false
    if (( ${#skip_years[@]} > 0 && ${#skip_months[@]} > 0 )); then
      [[ ${year_matches} == true && ${month_matches} == true ]] && skip_candidate=true
    elif (( ${#skip_years[@]} > 0 )); then
      [[ ${year_matches} == true ]] && skip_candidate=true
    elif (( ${#skip_months[@]} > 0 )); then
      [[ ${month_matches} == true ]] && skip_candidate=true
    fi
    if [[ ${skip_candidate} == false ]]; then
      task_years+=("${candidate_year}")
      task_months+=("${candidate_month_num}")
    fi
  done
done
task_count=${#task_years[@]}
(( task_count > 0 )) || {
  echo "ERROR: no dates remain after applying skip_years/skip_months" >&2
  exit 2
}
last_task_id=$((task_count - 1))

if [[ ${execution_mode} == auto ]]; then
  if [[ -n ${SLURM_JOB_ID:-} ]]; then
    execution_mode=run
  else
    execution_mode=submit
  fi
fi
if [[ ${execution_mode} == run && -z ${SLURM_JOB_ID:-} ]]; then
  echo "ERROR: --run requires an active Slurm allocation; use salloc/srun first or use --submit." >&2
  exit 2
fi

# Before either a local run or job submission, validate the complete inventory.
export my_overwrite="${overwrite_existing}"
source "${script_dir}/create_and_setup_case.sh"

  problems=0
  check_readable() {
    local path=$1
    if [[ ! -r ${path} ]]; then
      echo "MISSING: ${path}" >&2
      problems=$((problems + 1))
    fi
  }

  check_empty_output() {
    local path=$1
    local entries=()
    if [[ -d ${path} ]]; then
      shopt -s nullglob dotglob
      entries=("${path}"/*)
      shopt -u nullglob dotglob
      if (( ${#entries[@]} > 0 )) && [[ ${my_overwrite} != true ]]; then
        echo "EXISTS: ${path} (set my_overwrite=true to replace it)" >&2
        problems=$((problems + 1))
      fi
    fi
  }

  check_readable "${my_eami_template}"
  check_readable "${my_e3sm_topo}"
  check_readable "${my_e3sm_semap}"
  check_readable "${my_e3sm_csgrid}"

  for check_task_id in "${!task_years[@]}"; do
    check_year=${task_years[$check_task_id]}
    model_year=$((check_year - 1980 + 353))
      check_month_num=${task_months[$check_task_id]}
      calendar_date=$(printf '%04d-%02d-01' "${check_year}" "${check_month_num}")
      restart_date=$(printf '%04d-%02d-01' "${model_year}" "${check_month_num}")

      if [[ ${mode} == 4denvar-ocn ]]; then
        case ${check_year} in
          198[0-9]) da_subdir="1980_1989" ;;
          199[0-9]) da_subdir="1990_1999" ;;
          200[0-9]) da_subdir="2000-2009" ;;
          201[0-9]) da_subdir="2010-2019" ;;
          *)
            echo "MISSING: no 4DEnVar source directory mapping for year ${check_year}" >&2
            problems=$((problems + 1))
            continue
            ;;
        esac
        da_dir="${my_daroot}/${da_subdir}"
        if (( check_year < 1990 )); then
          da_atmosphere="${da_dir}/$(printf "${my_daeam_1980s}" "${calendar_date}")"
        else
          da_atmosphere="${da_dir}/$(printf "${my_daeam}" "${calendar_date}")"
        fi
        check_readable "${da_atmosphere}"
        check_readable "${da_dir}/$(printf "${my_daocn}" "${calendar_date}")"
        check_readable "${da_dir}/$(printf "${my_dalnd}" "${calendar_date}")"
        check_readable "${da_dir}/${my_dacase}.mosart.r.${calendar_date}-00000.nc"
        check_readable "${da_dir}/${my_dacase}.mpassi.rst.${calendar_date}_00000.nc"
        check_readable "${da_dir}/${my_dacase}.cpl.r.${calendar_date}-00000.nc"
        check_readable "${da_dir}/${my_dacase}.eam.r.${calendar_date}-00000.nc"
        check_readable "${da_dir}/${my_dacase}.eam.rs.${calendar_date}-00000.nc"
        output_dir="${my_runpath}/4DEnVarOcn/${calendar_date}-00000"
      else
        atmosphere="${my_atmdir2}/${my_atmcase2}.${calendar_date}.${my_atmtod2}.ne30np4.L80.nc"
        check_readable "${atmosphere}"

        restart_dir="${my_fosdir1}/${restart_date}-00000"
        check_readable "${restart_dir}/${my_foscase1}.elm.r.${restart_date}-00000.nc"
        check_readable "${restart_dir}/${my_foscase1}.mosart.r.${restart_date}-00000.nc"
        check_readable "${restart_dir}/${my_foscase1}.mpaso.rst.${restart_date}_00000.nc"
        check_readable "${restart_dir}/${my_foscase1}.mpassi.rst.${restart_date}_00000.nc"
        check_readable "${restart_dir}/${my_foscase1}.cpl.r.${restart_date}-00000.nc"
      fi

      if [[ ${mode} == bruteforce-cfsr ]]; then
        ocean_name=$(printf "${my_ocncase3}" "${restart_date}")
        check_readable "${my_ocndir3}/${ocean_name}"
        output_dir="${my_runpath}/BruteForce/${calendar_date}-00000"
      elif [[ ${mode} == reanalysis-cfsr ]]; then
        ocean_name=$(printf "${my_ocncase4}" "${calendar_date//-/}")
        check_readable "${my_ocndir4}/${ocean_name}"
        if [[ -n ${ice} ]]; then
          ice_ymd=${calendar_date//-/}
          check_readable "${my_icedir}/$(printf "${my_icecase}" "${ice_ymd}")"
        fi
        output_dir="${my_runpath}/Reanalysis/${calendar_date}-00000"
      elif [[ ${mode} == jra55-fosirl ]]; then
        output_dir="${my_runpath}/JRA55-FOSIRL/${calendar_date}-00000"
      fi
      check_empty_output "${output_dir}"
  done

  if (( problems > 0 )); then
    echo "Preflight failed: ${problems} missing input(s) or protected output directory issue(s)." >&2
    exit 1
  fi

  echo "Preflight passed for ${mode}: all ${task_count} dates are ready."
  [[ ${check_only} == true ]] && exit 0

if [[ ${execution_mode} == submit ]]; then
  exec sbatch \
    --export="ALL,IC_MODE=${mode},IC_WORK_DIR=${script_dir}" \
    "${script_dir}/submit_all_year.bash"
fi

work_dir=${IC_WORK_DIR:-${script_dir:-${SLURM_SUBMIT_DIR:-${PWD}}}}

run_task() {
  local task_id=$1
  local year month ymd task_label

  (( task_id >= 0 && task_id <= last_task_id )) || {
    echo "ERROR: task ${task_id} is outside 0-${last_task_id}" >&2
    return 2
  }
  year=${task_years[$task_id]}
  month=${task_months[$task_id]}
  ymd=$(printf "%04d-%02d-01" "${year}" "${month}")
  task_label="task $((task_id + 1))/${task_count}"

  cd "${work_dir}" || return 1
  echo "Running ${ymd} in ${mode} mode (${task_label})"

  # create_and_setup_case.sh preserves these per-date overrides, so the shared
  # configuration file does not need to be edited with sed.
  export my_casedate="${ymd}"
  export my_runcase="${group}"
  export my_fosopt="${fosi}"
  export my_eamopt="${atm}"
  export my_ocnopt="${ocn}"
  export my_iceopt="${ice}"
  export my_overwrite="${overwrite_existing}"

  if [[ ${DRY_RUN:-false} == true ]]; then
    source ${work_dir}/create_and_setup_case.sh
    printf 'mode=%s\ndate=%s\nrun_case=%s\natmosphere=%s\noutput=%s\nocean_option=%s\n' \
      "${mode}" "${my_casedate}" "${my_runcase}" "${my_eam_ic2}" \
      "${my_dart_runpath}" "${my_ocnopt:-<none>}"
    return 0
  fi

  ./gen_ensemble_s2d_init.sh
}

echo "Running ${task_count} dates sequentially in Slurm allocation ${SLURM_JOB_ID:-<unknown>}."
for task_id in "${!task_years[@]}"; do
  run_task "${task_id}" || exit $?
done

exit 0 
