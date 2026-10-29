#!/bin/bash -l
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --account=e3sm
#SBATCH  --constraint=cpu
#SBATCH  --qos=regular
#SBATCH  --nodes=1
#SBATCH  --time=01:00:00
#SBATCH  --job-name=rtenso_ice_ic
#SBATCH  --output=/global/cfs/cdirs/e3sm/S2S2D/RTENSO2026/initial_conditions/logs/ice_ic_%j.out

set -o pipefail

# ==============================================================================
# Sea-ice (MPAS-SI) initial conditions for RTENSO2026 from the reanalysis-based
# sea-ice files in Reanalysis_ICs/ (NSIDCv6 ice area, spun up). They do not
# depend on the land/sea-ice method (clim, obs, hybrid_ic), so one is written
# per date, as an alternative to the method sea ice in ice/<method>/:
#
#   <out_root>/YYYYMMDD/ice/nsidc/20251024_s2d_spinup.mpassi.rst.YYYY-MM-DD_00000.nc
#   <out_root>/YYYYMMDD/ice/nsidc/rpointer.ice, ice.json
#
# The file must be on the model mesh (IcoswISC30E3r5) and hold every variable
# of a model sea-ice restart (checked against the BruteForce restart of
# ref_year). The ocean is handled separately by make_ocn_ic.bash.
#
# Edit USER SETTINGS, then:
#   CHECK_ONLY=true ./make_ice_ic.bash   preflight only: inputs, mesh, variables
#   DATES="2026-02-01 2026-05-01" OVERWRITE=true ./make_ice_ic.bash --run
#                                        restage a subset of dates over existing restarts
#   ./make_ice_ic.bash                   execution_mode below (default auto: run here on an
#                                        interactive node or in a job, submit from a login node)
#   ./make_ice_ic.bash --run             run here, e.g. on an interactive node
#                                        (salloc -N 1 -C cpu -q interactive -t 02:00:00 -A e3sm)
#   ./make_ice_ic.bash --submit          submit one Slurm job
# ==============================================================================

# ---------------------------- USER SETTINGS ----------------------------------
out_root="/global/cfs/cdirs/e3sm/S2S2D/RTENSO2026/initial_conditions"
dates=(2015-02-01 2015-05-01 2015-08-01 2015-11-01
       2026-02-01 2026-05-01 2026-08-01)
[[ -n ${DATES:-} ]] && read -ra dates <<< "${DATES}"   # DATES="2026-02-01 ..." runs a subset
out_subdir="ice/nsidc"                # written to <out_root>/YYYYMMDD/<out_subdir>/

# Sea-ice initial conditions, one directory per calendar date (YYYYMMDD; both
# %s are the date), remapped to IcoswISC30E3r5.
ice_dir="/global/cfs/cdirs/e3sm/S2S2D/RTENSO2026/initial_conditions/Reanalysis_ICs"
ice_case="%s/mpassi.IcoswISC30E3r5.iceArea_from_NSIDCv6_%s_spunup.20261029.nc"
ice_label="NSIDCv6"                   # recorded in ice.json and the file attributes

# Model sea-ice restart the files are checked against (mesh and variables).
src_root="/global/cfs/cdirs/e3sm/zhan391/v3.LR.S2D.ENSINT/BruteForce"
ref_year=2018
prefix="20251024_s2d_spinup"          # restart prefix (my_refcase2 in the hindcast driver)

overwrite_existing="${OVERWRITE:-false}"  # replace an existing sea-ice restart
check_only="${CHECK_ONLY:-false}"     # true: preflight only
# auto: run here inside a Slurm allocation (interactive node or batch job), submit from a
# login node; run: always run here; submit: always submit a Slurm job. --run/--submit override it.
execution_mode="${EXECUTION_MODE:-auto}"
# -----------------------------------------------------------------------------

script_dir=${IC_WORK_DIR:-$(cd "$(dirname "$0")" && pwd)}
export IC_WORK_DIR=${script_dir}
source "${script_dir}/util/ic_lib.sh"
parse_execution_mode "$@"
if ! command -v ncks >/dev/null; then
  source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh >/dev/null 2>&1
fi

# ---------------------------------------------------------------------------
# Preflight: every input readable, on the model mesh, with every restart
# variable; outputs not protected.
# ---------------------------------------------------------------------------
problems=0
tasks=()
for date in "${dates[@]}"; do
  ice_file="${ice_dir}/$(printf "${ice_case}" "${date//-/}" "${date//-/}")"
  ref_date="${ref_year}${date:4}"
  ref="${src_root}/${ref_date}-00000/${prefix}.mpassi.rst.${ref_date}_00000.nc"
  out="${out_root}/${date//-/}/${out_subdir}/${prefix}.mpassi.rst.${date}_00000.nc"
  [[ -r ${ice_file} ]] || { echo "MISSING: ${ice_file}" >&2; problems=$((problems + 1)); continue; }
  [[ -r ${ref} ]] || { echo "MISSING: ${ref}" >&2; problems=$((problems + 1)); continue; }
  if [[ -e ${out} && ${overwrite_existing} != true ]]; then
    echo "EXISTS: ${out} (set overwrite_existing=true to replace it)" >&2
    problems=$((problems + 1)); continue
  fi
  tasks+=("${date}|${ref}|${ice_file}|${out}")
done

if (( ${#tasks[@]} > 0 )); then           # all dates share the mesh and variables: check one
  IFS='|' read -r _ ref ice_file _ <<< "${tasks[0]}"
  n_ref=$(ncells "${ref}"); n_ice=$(ncells "${ice_file}")
  if [[ ${n_ice} != "${n_ref}" ]]; then
    echo "MESH: ${ice_file} has nCells=${n_ice}, the model sea-ice restart ${n_ref}: remap it to the model mesh first" >&2
    problems=$((problems + 1))
  else
    missing=$(comm -23 <(nc_vars "${ref}") <(nc_vars "${ice_file}") | tr '\n' ' ')
    [[ -z ${missing// } ]] || {
      echo "VARIABLES: ${ice_file} lacks restart variables: ${missing}" >&2; problems=$((problems + 1)); }
  fi
fi

if (( problems > 0 )); then
  echo "Preflight failed: ${problems} problem(s)." >&2
  exit 1
fi
echo "Preflight passed: ${#tasks[@]} dates (${ice_label} sea ice)."
for t in "${tasks[@]}"; do
  IFS='|' read -r date _ ice_file out <<< "${t}"
  echo "  ${date}: $(basename "${ice_file}") -> ${out}"
done
[[ ${check_only} == true ]] && exit 0
submit_or_run "${script_dir}/make_ice_ic.bash" "${out_root}/logs"

# ---------------------------------------------------------------------------
# Stage each sea-ice file as a model restart.
# ---------------------------------------------------------------------------
for t in "${tasks[@]}"; do
  IFS='|' read -r date _ ice_file out <<< "${t}"
  dir=$(dirname "${out}")
  mkdir -p "${dir}"
  tmp="${dir}/.$(basename "${out}").tmp"
  echo "$(date '+%H:%M:%S') ${date}: $(basename "${ice_file}")"
  stage_restart "${ice_file}" "${tmp}" || exit 11
  ncatted -O -a rtenso2026_sea_ice,global,o,c,"${ice_label} sea ice ${ice_file}" "${tmp}" || exit 12
  mv -f "${tmp}" "${out}"
  echo "${date}_00:00:00" > "${dir}/rpointer.ice"
  cat > "${dir}/ice.json" <<EOF
{
  "target_date": "${date}",
  "source": "${ice_label}",
  "file": "${ice_file}",
  "output": "${out}",
  "created": "$(date -Iseconds)"
}
EOF
done
echo "$(date '+%H:%M:%S') done: ${#tasks[@]} dates"
