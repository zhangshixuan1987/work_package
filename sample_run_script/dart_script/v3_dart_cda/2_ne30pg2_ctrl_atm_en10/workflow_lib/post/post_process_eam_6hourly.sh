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

var2_list=(FLUT FLUTC LHFLX SHFLX PRECT PRECC PRECL PS PSL QFLX QREFHT TMQ TS TREFHT TUQ TVQ OMEGA500 PRECSL SWCF LWCF
           TTOP TAUX TAUY TGCLDCWP TGCLDIWP TGCLDLWP U90M V90M U10 CLDTOT CLDLOW CLDMED FLDS FLDSC FLNS FLNSC
           FLNT FLNTC FSDS FSDSC FSNS FSNSC FSNT FSNTC FSNTOA FSNTOAC FSUTOA)

hist="eam.h2"
freq="6hourly"

input="${ARCHIVE_DIR}/atm/hist"
outdir="${ARCHIVE_DIR}/post"

mkdir -p "${outdir}"

jobid=${SLURM_JOB_ID:-$$}

ENSTR="${POST_ENSTR}"
CASE_NAME="${POST_CASE_NAME}"

cd ${outdir}
workdir=$(mktemp -d tmp.${jobid}.XXXX)
cd ${workdir}

flist=$(mktemp ./input${ENSTR}.XXXXXX)

echo "=== Starting ensemble member ${CASE_NAME}.${hist}.${freq} ==="

for year in $(seq "${sy}" "${ey}"); do
  for month in $(seq 1 12); do
    # Skip months outside the desired range if first/last year
    if [ "$year" -eq "$sy" ] && [ "$month" -lt "$sm" ]; then continue; fi
    if [ "$year" -eq "$ey" ] && [ "$month" -gt "$em" ]; then continue; fi  # note: > not >=

    yymm=$(printf "%04d-%02d" "${year}" "${month}")

    # Link files like ${CASE_NAME}.${hist}.YYYY-MM*.nc
    for ff in ${input}/${CASE_NAME}.${hist}.${yymm}-??-?????.nc; do
      stamp=${ff##*.${hist}.}
      stamp=${stamp%.nc}
      file_date=${stamp:0:10}
      [[ "${file_date}" < "${ymds}" || "${file_date}" > "${ymde}" ]] && continue
      ln -sf -- "${ff}" .
    done
  done
done

ls ${CASE_NAME}.${hist}.????-??-??-?????.nc > ${flist}

# === Step 3: Regrid var2_list from hist (6-hourly) ===
ts_dest1="${outdir}/atm/180x360_aave/ts/${freq}_da"
ts_dest2="${outdir}/atm/180x360_aave/ts/${freq}"
mkdir -p "${ts_dest1}" "${ts_dest2}"

# Build an array of input files safely
mapfile -t ffiles2 < "${flist}"
echo "hist inputs (${#ffiles2[@]} files):"

for year in $(seq "${sy}" "${ey}"); do
  yyyy=$(printf "%04d" "${year}")
  start_time="${yyyy}-01-01 00:00:0.0"
  end_time="${yyyy}-12-31 23:59:59.0"

  for var in "${var2_list[@]}"; do
    outfile="${var}.${ENSTR}.${yyyy}.nc"

    SCRATCH=$(mktemp -d scrach.${jobid}.XXXX)

    if ((${#ffiles2[@]} > 0)); then
      # Filter ffiles2 to only those that actually contain $var
      present_files=()
      for f in "${ffiles2[@]}"; do
        # Succeeds (exit code 0) only if $var exists in $f
        if ncks -C -H -v "$var" "$f" >/dev/null 2>&1; then
          present_files+=("$f")
        fi
      done

      printf 'Files containing %s:\n' "$var"
      if ((${#present_files[@]} == 0)); then
        echo "Skipping ${var} for ${ENSTR} ${yyyy} — variable not found in any hist2 file."
        continue
      fi
      
      # Per-var/year temp names to avoid collisions
      tmp_all="${SCRATCH}/${var}_${ENSTR}_${yyyy}_all.nc"
      tmp_prior="${SCRATCH}/${var}_${ENSTR}_${yyyy}_prior.nc"
      tmp_post="${SCRATCH}/${var}_${ENSTR}_${yyyy}_post.nc"
      tmp_time0="${SCRATCH}/${var}_${ENSTR}_${yyyy}_t0.nc"
      tmp_time1="${SCRATCH}/${var}_${ENSTR}_${yyyy}_t1.nc"

      # 1) Gather this year's 6-hourly slices for the variable from files that have it
      ncrcat -O -d time,"${start_time}","${end_time}" -v "${var}" \
        "${present_files[@]}" "${tmp_all}"

      # 2) Even-index (0,2,4,...) records: "prior"
      ncks -O -d time,0,,2 "${tmp_all}" "${tmp_prior}"
      ncremap -m "${MAP_FILE}" -i "${tmp_prior}" -o "${ts_dest1}/${outfile}.tmp.${jobid}"
      ncdump -h "${ts_dest1}/${outfile}.tmp.${jobid}" >/dev/null 2>&1
      ${MOVE} "${ts_dest1}/${outfile}.tmp.${jobid}" "${ts_dest1}/${outfile}"

      # 3) First two records (0 and 1) combined as "post" (or keep odd-only if that's intended)
      ncks -O -d time,1,,2        "${tmp_all}"     "${tmp_time1}"  # odd indices
      ncks -O -d time,0           "${tmp_all}"     "${tmp_time0}"  # first record
      ncrcat -O -d time,0,        "${tmp_time0}"   "${tmp_time1}"  "${tmp_post}"
      ncremap -m "${MAP_FILE}" -i "${tmp_post}" -o "${ts_dest2}/${outfile}.tmp.${jobid}"
      ncdump -h "${ts_dest2}/${outfile}.tmp.${jobid}" >/dev/null 2>&1
      ${MOVE} "${ts_dest2}/${outfile}.tmp.${jobid}" "${ts_dest2}/${outfile}"

      # Optional compression
      # ncks -O -4 -L 1 "${ts_dest1}/${outfile}" "${ts_dest1}/${outfile}"
      # ncks -O -4 -L 1 "${ts_dest2}/${outfile}" "${ts_dest2}/${outfile}"
    else
      echo "Warning: No hist2 files found for ${ENSTR}, year ${yyyy}, var ${var}"
    fi

    rm -rvf ${SCRATCH}

  done
done

cd ..
rm -rf "${workdir}"

echo "===== End of DART diagnostic ====="
date
echo "==================================="
