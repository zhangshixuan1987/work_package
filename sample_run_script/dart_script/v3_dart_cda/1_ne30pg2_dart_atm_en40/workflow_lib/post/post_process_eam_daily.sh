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

#set -Eeo pipefail
#shopt -s nullglob

if [[ "${POSTPROCESS_DRIVER_ACTIVE:-FALSE}" != "TRUE" ]]; then
  echo "ERROR: internal post-processing worker; run 7_run_post_hist.sh" >&2
  exit 1
fi
echo "== Start of DART diagnostic =="
date
echo "============================================"

# System utilities
MOVE='/usr/bin/mv'
COPY='/usr/bin/cp --preserve=timestamps'
LINK='/usr/bin/ln -fs'
REMOVE='/usr/bin/rm'
LIST='/usr/bin/ls'

# Environment setup (assumes these are exported externally or in create_and_setup_case.sh)
# E3SM_ROOT, DART_ROOT, my_modeldir, my_ensnum, my_casename, etc.
if [[ "${my_analysis_env_type:-}" == "custom" ]]; then
    echo "Activating custom analysis environment: ${my_analysis_conda_env}"

    source "${my_conda_setup_file}" ||
        fail "could not source Conda setup file: ${my_conda_setup_file}"

    conda activate "${my_analysis_conda_env}" ||
        fail "could not activate Conda environment: ${my_analysis_conda_env}"
else
    my_analysis_pixi_dir="/share/apps/E3SM/conda_envs/e3smu_1_13_0/compy/pixi_login"
    my_pixi="${my_analysis_pixi_dir}/bin/pixi"
    echo "Activating shared E3SM Unified analysis environment"
    if [[ ! -x "${my_pixi}" ]]; then
        echo "ERROR: Pixi executable not found: ${my_pixi}" >&2
        exit 1
    fi
    eval "$(
        "${my_pixi}" shell-hook \
            --as-is \
            --shell bash \
            --manifest-path "${my_analysis_pixi_dir}/pixi.toml"
    )"
    hash -r
    echo "CONDA_DEFAULT_ENV=${CONDA_DEFAULT_ENV:-unset}"
    echo "CONDA_PREFIX=${CONDA_PREFIX:-unset}"
    echo "ncremap=$(command -v ncremap)"
    echo "ncks=$(command -v ncks)"
fi

DART_ROOT="${my_eam_dart_code}"
DART_MODEL=${my_eam_dart_model}
DART_WORKDIR=${DART_ROOT}/models/${DART_MODEL}/work
ARCHIVE_DIR="${my_modeldir}/${POST_ENSTR}/archive"
MAP_FILE="${my_eam_hist_mapping_file}"

# Dates
ymds="${POST_START_DATE}"
ymde="${POST_END_DATE}"

read -r sy sm sd <<< "$(echo ${ymds} | tr '-' ' ')"
read -r ey em ed <<< "$(echo ${ymde} | tr '-' ' ')"
mday=(31 28 31 30 31 30 31 31 30 31 30 31)

# Variables
var_list=(TCO SCO FLUT PRECT PS PSL TS TBOT TUQ TVQ UBOT VBOT ZBOT U10 TAUX TAUY TREFHT QREFHT TMQ CAPE CIN TTQ
           U1000 U975 U950 U925 U900 U850 U800 U700 U600 U500 U400 U300 U200 U100 U050 U010
           V1000 V975 V950 V925 V900 V850 V800 V700 V600 V500 V400 V300 V200 V100 V050 V010
           Z1000 Z975 Z950 Z925 Z900 Z850 Z800 Z700 Z600 Z500 Z400 Z300 Z200 Z100 Z050 Z010
           T1000 T975 T950 T925 T900 T850 T800 T700 T600 T500 T400 T300 T200 T100 T050 T010
           Q1000 Q975 Q950 Q925 Q900 Q850 Q800 Q700 Q600 Q500 Q400 Q300 Q200 Q100 Q050 Q010
           OMEGA1000 OMEGA975 OMEGA950 OMEGA925 OMEGA900 OMEGA850 OMEGA800 OMEGA700 OMEGA600 OMEGA500
           OMEGA400 OMEGA300 OMEGA200 OMEGA100 TREFHTMN TREFHTMX PRECTMX)

hist="eam.h1"
freq="daily"
input="${ARCHIVE_DIR}/atm/hist"
outdir="${ARCHIVE_DIR}/post"

mkdir -p "${outdir}"

jobid=${SLURM_JOB_ID:-$$}

cd ${outdir}
workdir=$(mktemp -d tmp.${jobid}.XXXX)
cd ${workdir}

ENSTR="${POST_ENSTR}"
CASE_NAME="${POST_CASE_NAME}"

flist=$(mktemp ./input${ENSTR}.XXXXXX)

echo "=== Starting ensemble member ${CASE_NAME}.${hist}.${freq} ==="

ts_dest="${outdir}/atm/180x360_aave/ts/${freq}"
mkdir -p "${ts_dest}"

for year in $(seq "${sy}" "${ey}"); do
  yyyy=$(printf "%04d" "${year}")

  # Per-year scratch
  SCRATCH=$(mktemp -d scrach.${jobid}.XXXX)
  echo "Building daily temps for ${ENSTR} ${yyyy} in ${SCRATCH}"

  for month in $(seq 1 12); do
    # honor start/end month bounds on first/last year
    if [ "${year}" -eq "${sy}" ] && [ "${month}" -lt "${sm}" ]; then continue; fi
    if [ "${year}" -eq "${ey}" ] && [ "${month}" -gt "${em}" ]; then continue; fi

    mm=$(printf "%02d" "${month}")
    yymm="${yyyy}-${mm}"

    # days in month (handle leap-year Feb only)
    nday=${mday[$((month - 1))]}
    if (( month == 2 )) && (( (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0) )); then
      nday=29
    fi

    # For each day, average all hist files into a single daily file
    for dd in $(seq 1 "${nday}"); do
      if (( year == sy && month == sm && dd < 10#${ymds:8:2} )); then continue; fi
      if (( year == ey && month == em && dd > 10#${ymde:8:2} )); then continue; fi
      dd2=$(printf "%02d" "${dd}")

      # Gather that day's input files directly from $input; no symlinks
      day_glob="${input}/${CASE_NAME}.${hist}.${yymm}-${dd2}"-?????.nc
      day_files=( ${day_glob} )

      if ((${#day_files[@]} == 0)); then
        # No files for this day -> skip
        continue
      fi

      # One daily mean file per day
      out_day="${SCRATCH}/ftmp_${yyyy}-${mm}-${dd2}.nc"
      ncra -O "${day_files[@]}" "${out_day}"
    done
  done

  # Collect all daily files for the year and sort (null-safe, portable)
  ffiles1=()
  while IFS= read -r -d '' f; do
    ffiles1+=("$f")
  done < <(find "$SCRATCH" -maxdepth 1 -type f -name "ftmp_${yyyy}-*.nc" -print0 | sort -z)

  if [ ${#ffiles1[@]} -eq 0 ]; then
    echo "Warning: No daily files built for ${ENSTR} ${yyyy} — skipping yearly aggregate."
    rm -rf -- "$SCRATCH"
    continue
  fi

  # Aggregate per variable, then regrid
  for var in "${var_list[@]}"; do
    outfile="${var}.${ENSTR}.${yyyy}.nc"
    tmp_year="./${var}_${ENSTR}_${yyyy}_year.nc"

    # Require complete daily coverage so annual products are not partial time series.
    files_ok=()
    for f in "${ffiles1[@]}"; do
      if ncks -C -H -v "$var" "$f" >/dev/null 2>&1; then
        files_ok+=("$f")
      fi
    done

    if ((${#files_ok[@]} == 0)); then
      echo "Skipping ${var} — not present in any daily file for ${ENSTR} ${yyyy}."
      continue
    fi
    if ((${#files_ok[@]} != ${#ffiles1[@]})); then
      echo "WARNING: Skipping ${var} — present in only ${#files_ok[@]} of ${#ffiles1[@]} daily files for ${ENSTR} ${yyyy}." >&2
      continue
    fi

    echo "Variable ${var} has complete daily coverage — concatenating and regridding..."
    # Retry transient shared-filesystem failures without allowing an incomplete member.
    concat_ok=0
    for attempt in 1 2 3; do
      rm -f -- "${tmp_year}" "${tmp_year}".pid*.ncrcat.tmp
      if ncrcat -O -v "$var" "${ffiles1[@]}" "$tmp_year"; then
        concat_ok=1
        break
      fi
      echo "WARNING: ncrcat failed for ${var} (attempt ${attempt}/3)." >&2
      (( attempt == 3 )) || sleep 2
    done
    if (( concat_ok == 0 )); then
      echo "ERROR: ncrcat failed for ${var} after 3 attempts; failing ${ENSTR}." >&2
      exit 1
    fi

    # Regrid to a temporary product, validate it, then publish atomically.
    if ! ncremap -m "${MAP_FILE}" -i "${tmp_year}" -o "${ts_dest}/${outfile}.tmp.${jobid}"; then
      echo "ERROR: ncremap failed for ${var}." >&2
      exit 1
    fi
    if ! ncdump -h "${ts_dest}/${outfile}.tmp.${jobid}" >/dev/null 2>&1; then
      echo "ERROR: invalid regridded NetCDF for ${var}: ${ts_dest}/${outfile}.tmp.${jobid}" >&2
      exit 1
    fi
    if ! ${MOVE} "${ts_dest}/${outfile}.tmp.${jobid}" "${ts_dest}/${outfile}"; then
      echo "ERROR: could not publish daily output for ${var}: ${ts_dest}/${outfile}" >&2
      exit 1
    fi

    # Optional lightweight compression of the final product
    # ncks -O -4 -L 1 "${ts_dest}/${outfile}" "${ts_dest}/${outfile}"

  done

  # Clean per-year scratch
  rm -rf "${SCRATCH}"
done

cd ..
rm -rf "${workdir}"

echo "===== End of DART diagnostic ====="
date
echo "==================================="
