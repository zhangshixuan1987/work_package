#!/bin/bash -el
#------------------------------------------------------------------------------
# SLURM Batch Directives
#------------------------------------------------------------------------------
#SBATCH --account=esmd
#SBATCH --time=02:00:00
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
MAP_FILE="${my_elm_hist_mapping_file}"

# Dates
ymds="${POST_START_DATE}"
ymde="${POST_END_DATE}"

read -r sy sm sd <<< "$(echo ${ymds} | tr '-' ' ')"
read -r ey em ed <<< "$(echo ${ymde} | tr '-' ' ')"
mday=(31 28 31 30 31 30 31 31 30 31 30 31)

# Variables
var_list=(
  H2OSNO FSNO QRUNOFF QSNOMELT FSNO_EFF SNORDSL SNOW FSA FSDS FSR FLDS
  FIRE FIRA SOILWATER_10CM SOILLIQ SOILICE QSOIL U10 U10WITHGUSTS TWS
  TSOI_10CM TSA TLAI THBOT TSOI TSOI_ICE SOILLIQ_ICE SOILICE_ICE
  TAUX TAUY FSH HC HCSOI EFLX_LH_TOT SNOW_DEPTH RH2M RAIN QVEGE QVEGT
  QBOT Q2M H2OSOI H2OSFC ZWT ZBOT TBOT TG THBOT PBOT
)

hist="elm.h1"
freq="daily"
input="${ARCHIVE_DIR}/lnd/hist"
outdir="${ARCHIVE_DIR}/post"

mkdir -p "${outdir}"

jobid=${SLURM_JOB_ID:-$$}

cd ${outdir}
workdir=$(mktemp -d tmp.${jobid}.XXXX)
cd ${workdir}

ENSTR="${POST_ENSTR}"
CASE_NAME="${POST_CASE_NAME}"

echo "=== Starting ensemble member ${CASE_NAME}.${hist}.${freq} ==="

# === Step 2: Regrid var_list from hist (daily) ===
ts_dest="${outdir}/lnd/180x360_aave/ts/${freq}"
mkdir -p "${ts_dest}"

for year in $(seq "${sy}" "${ey}"); do
  yyyy=$(printf "%04d" "${year}")

  # Per-year scratch
  SCRATCH=$(mktemp -d "./tmp_hist_${ENSTR}_${yyyy}.XXXXXX")
  echo "Building daily temps for ${ENSTR} ${yyyy} in ${SCRATCH}"

  # Build daily means for months within bounds
  for month in $(seq 1 12); do
    # honor start/end bounds in first/last year
    if [ "${year}" -eq "${sy}" ] && [ "${month}" -lt "${sm}" ]; then continue; fi
    if [ "${year}" -eq "${ey}" ] && [ "${month}" -gt "${em}" ]; then continue; fi

    mm=$(printf "%02d" "${month}")
    yymm="${yyyy}-${mm}"

    # days in month (leap-year Feb only)
    nday=${mday[$((month - 1))]}
    if (( month == 2 )) && (( (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0) )); then
      nday=29
    fi

    for dd in $(seq 1 "${nday}"); do
      if (( year == sy && month == sm && dd < 10#${ymds:8:2} )); then continue; fi
      if (( year == ey && month == em && dd > 10#${ymde:8:2} )); then continue; fi
      dd2=$(printf "%02d" "${dd}")
      day_glob="${input}/${CASE_NAME}.${hist}.${yymm}-${dd2}"-?????.nc
      day_files=( ${day_glob} )

      # no files for this day → skip
      ((${#day_files[@]}==0)) && continue

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

  # Probe once for landfrac presence (in the first daily file)
  have_landfrac=0
  if ncks -C -H -v landfrac "${ffiles1[0]}" >/dev/null 2>&1; then
    have_landfrac=1
  fi

  for var in "${var_list[@]}"; do
    outfile="${var}.${ENSTR}.${yyyy}.nc"
    tmp_year="./${var}_${ENSTR}_${yyyy}_year.nc"

    # Build list of variables we need in the concat
    need_vars=("$var")
    if (( have_landfrac )); then
      need_vars+=("landfrac")
    fi

    # Filter files that contain *all* needed variables
    files_ok=()
    for f in "${ffiles1[@]}"; do
      ok=1
      for v in "${need_vars[@]}"; do
        ncks -C -H -v "$v" "$f" >/dev/null 2>&1 || { ok=0; break; }
      done
      (( ok )) && files_ok+=("$f")
    done

    if ((${#files_ok[@]} == 0)); then
      echo "Skipping ${var} — none of the daily files contain all of: ${need_vars[*]}"
      continue
    fi

    if ((${#files_ok[@]} != ${#ffiles1[@]})); then
      echo "WARNING: Skipping ${var} — present in only ${#files_ok[@]} of ${#ffiles1[@]} daily files for ${ENSTR} ${yyyy}." >&2
      continue
    fi

    # Compose comma-separated list for -v without spaces
    sel_vars=$(IFS=,; echo "${need_vars[*]}")

    # Yearly concat using only files that have all requested variables
    ncrcat -O -v "${sel_vars}" "${files_ok[@]}" "${tmp_year}"

    if (( have_landfrac )); then
      ncremap -P elm -m "${MAP_FILE}" -i "${tmp_year}" -o "${ts_dest}/${outfile}.tmp.${jobid}"
    else
      echo "landfrac not found for ${ENSTR} ${yyyy}; using plain ncremap for ${var}."
      ncremap -m "${MAP_FILE}" -i "${tmp_year}" -o "${ts_dest}/${outfile}.tmp.${jobid}"
    fi
    if ! ncdump -h "${ts_dest}/${outfile}.tmp.${jobid}" >/dev/null 2>&1; then
      echo "ERROR: invalid regridded NetCDF for ${var}: ${ts_dest}/${outfile}.tmp.${jobid}" >&2
      exit 1
    fi
    if ! ${MOVE} "${ts_dest}/${outfile}.tmp.${jobid}" "${ts_dest}/${outfile}"; then
      echo "ERROR: could not publish daily output for ${var}: ${ts_dest}/${outfile}" >&2
      exit 1
    fi

    # Optional compression:
    # ncks -O -4 -L 1 "${ts_dest}/${outfile}" "${ts_dest}/${outfile}"
  done

  rm -rf "${SCRATCH}"
done

wait 

cd ..
rm -rf "${workdir}"

echo "===== End of DART diagnostic ====="
date
echo "==================================="
