#!/bin/bash -l
#SBATCH  --job-name=rtenso_clim_ic
#SBATCH  --account=e3sm
#SBATCH  --nodes=1
#SBATCH  --qos=regular
#SBATCH  --constraint=cpu
#SBATCH  --time=02:00:00
#SBATCH  --output=/global/cfs/cdirs/e3sm/S2S2D/RTENSO2026/initial_conditions/logs/clim_ic_%A_%a.out
#SBATCH  --array=0-7

# Select a BruteForce year and copy its land, river, and sea-ice restarts for
# each RTENSO2026 date and branch (clim: closest to climatology; obs: closest
# to ERA5). One array task per date runs every branch. Submit from this
# directory:
#   sbatch submit_clim_restarts.bash
# Logs go to the output root (initial_conditions/logs/), as does the key-field
# cache (initial_conditions/cache/).
# Print distances and selections for every date without copying anything:
#   DRY_RUN=true bash submit_clim_restarts.bash

# ---------------------------- USER SETTINGS ----------------------------------
src_root="/global/cfs/cdirs/e3sm/zhan391/v3.LR.S2D.ENSINT/BruteForce"
out_root="/global/cfs/cdirs/e3sm/S2S2D/RTENSO2026/initial_conditions"

# Observation data for land (ELM) and sea ice (MPAS-Seaice).
# Leave empty ("") to use the default ERA5 setup (/global/cfs/projectdirs/m3522/cmip6/ERA5/e5.oper.an.sfc).
obs_root_elm=""
obs_root_mpassi=""

dates=(2015-02-01 2015-05-01 2015-08-01 2015-11-01
       2026-02-01 2026-05-01 2026-08-01 2026-11-01)
branches=(clim obs)      # obs stops with a message if ERA5 for the date is missing
extra_args=()            # e.g. (--first-year 1991 --mode per_component --overwrite)
#extra_args=(--overwrite)
#extra_args=(--box 30 60 --overwrite)
#extra_args=(--mode per_component --overwrite)
#extra_args=(--component-weights mpassi=2 elm=1 --first-year 1991 --overwrite)

# -----------------------------------------------------------------------------

source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
work_dir=${SLURM_SUBMIT_DIR:-${PWD}}
cd "${work_dir}"

run_date() {
  local date=$1
  local obs_args=()
  [[ -n "${obs_root_elm}" ]] && obs_args+=(--obs-root-elm "${obs_root_elm}")
  [[ -n "${obs_root_mpassi}" ]] && obs_args+=(--obs-root-mpassi "${obs_root_mpassi}")
  for branch in "${branches[@]}"; do
    python util/make_clim_restart.py --date "${date}" --branch "${branch}" \
      --src-root "${src_root}" --out-root "${out_root}" \
      "${obs_args[@]}" "${extra_args[@]}" \
      $([[ ${DRY_RUN:-false} == true ]] && echo --dry-run)
  done
}

if [[ -z ${SLURM_ARRAY_TASK_ID:-} ]]; then
  for date in "${dates[@]}"; do
    run_date "${date}"
  done
  exit
fi

(( SLURM_ARRAY_TASK_ID < ${#dates[@]} )) || { echo "no task ${SLURM_ARRAY_TASK_ID}"; exit 1; }
date=${dates[SLURM_ARRAY_TASK_ID]}
echo "task ${SLURM_ARRAY_TASK_ID}: ${date}"
run_date "${date}"
