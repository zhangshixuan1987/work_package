#!/bin/bash -l
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --account=e3sm
#SBATCH  --constraint=cpu
#SBATCH  --qos=regular
#SBATCH  --nodes=1
#SBATCH  --time=02:00:00
#SBATCH  --job-name=rtenso_ocn_ic
#SBATCH  --output=/global/cfs/cdirs/e3sm/S2S2D/RTENSO2026/initial_conditions/logs/ocn_ic_%j.out

set -o pipefail

# ==============================================================================
# Ocean (MPAS-O) initial conditions for RTENSO2026 from the CFSR ocean analysis
# in Reanalysis_ICs/ (mpas_YYYYMMDD.cfsr.ic.nc, Luke Van Roekel). The ocean does
# not depend on the land/sea-ice method (clim, obs, hybrid_ic), so one is
# written per date:
#
#   <out_root>/YYYYMMDD/ocn/20251024_s2d_spinup.mpaso.rst.YYYY-MM-DD_00000.nc
#   <out_root>/YYYYMMDD/ocn/rpointer.ocn, ocean.json
#
# As reanalysis-cfsr in ../Reanalysis/ and ../BruteForce/luke_cfsr/, the CFSR
# ocean state (temperature, salinity, layer thickness, land-ice fields;
# velocities are zero, as in the BruteForce CFSRv2 restarts) is appended onto
# an mpaso restart, which supplies only the fields CFSR does not have (flux
# accumulators, SSH history, frazil). No 2026 restart exists, so the template
# is the BruteForce mpaso restart of template_year on the same calendar date.
# The sea ice is handled separately by make_ice_ic.bash.
#
# Edit USER SETTINGS, then:
#   CHECK_ONLY=true ./make_ocn_ic.bash   preflight only: inputs, mesh, protected outputs
#   ./make_ocn_ic.bash                   execution_mode below (default auto: run here on an
#                                        interactive node or in a job, submit from a login node)
#   ./make_ocn_ic.bash --run             run here, e.g. on an interactive node
#                                        (salloc -N 1 -C cpu -q interactive -t 02:00:00 -A e3sm)
#   ./make_ocn_ic.bash --submit          submit one Slurm job
# About 1 minute and 4 GB per date.
# ==============================================================================

# ---------------------------- USER SETTINGS ----------------------------------
out_root="/global/cfs/cdirs/e3sm/S2S2D/RTENSO2026/initial_conditions"
dates=(2015-02-01 2015-05-01 2015-08-01 2015-11-01
       2026-02-01 2026-05-01 2026-08-01)
out_subdir="ocn"                      # written to <out_root>/YYYYMMDD/<out_subdir>/

# Ocean analysis, named by calendar date (YYYYMMDD).
ocn_dir="/global/cfs/cdirs/e3sm/zhan391/v3.LR.S2D.ENSINT/Reanalysis_ICs"
ocn_case="mpas_%s.cfsr.ic.nc"

# Template mpaso restart (BruteForce archive, read-only): template_year on the
# target calendar date. It supplies only the fields CFSR does not have.
src_root="/global/cfs/cdirs/e3sm/zhan391/v3.LR.S2D.ENSINT/BruteForce"
template_year=2018                    # latest archive year
prefix="20251024_s2d_spinup"          # restart prefix (my_refcase2 in the hindcast driver)

overwrite_existing=false              # replace an existing mpaso restart
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
# Preflight: every input readable, meshes match, outputs not protected.
# ---------------------------------------------------------------------------
problems=0
tasks=()
for date in "${dates[@]}"; do
  ocn_file="${ocn_dir}/$(printf "${ocn_case}" "${date//-/}")"
  tmpl_date="${template_year}${date:4}"
  tmpl="${src_root}/${tmpl_date}-00000/${prefix}.mpaso.rst.${tmpl_date}_00000.nc"
  out="${out_root}/${date//-/}/${out_subdir}/${prefix}.mpaso.rst.${date}_00000.nc"
  [[ -r ${ocn_file} ]] || { echo "MISSING: ${ocn_file}" >&2; problems=$((problems + 1)); continue; }
  [[ -r ${tmpl} ]] || { echo "MISSING: ${tmpl}" >&2; problems=$((problems + 1)); continue; }
  if [[ -e ${out} && ${overwrite_existing} != true ]]; then
    echo "EXISTS: ${out} (set overwrite_existing=true to replace it)" >&2
    problems=$((problems + 1)); continue
  fi
  tasks+=("${date}|${tmpl}|${ocn_file}|${out}")
done

if (( ${#tasks[@]} > 0 )); then           # all dates share the mesh: check one
  IFS='|' read -r _ tmpl ocn_file _ <<< "${tasks[0]}"
  n_tmpl=$(ncells "${tmpl}"); n_ocn=$(ncells "${ocn_file}")
  [[ ${n_ocn} == "${n_tmpl}" ]] || {
    echo "MESH: ${ocn_file} has nCells=${n_ocn}, the mpaso restart ${n_tmpl}" >&2; problems=$((problems + 1)); }
fi

if (( problems > 0 )); then
  echo "Preflight failed: ${problems} problem(s)." >&2
  exit 1
fi
echo "Preflight passed: ${#tasks[@]} dates (CFSR ocean on the ${template_year} mpaso restart)."
for t in "${tasks[@]}"; do
  IFS='|' read -r date _ ocn_file out <<< "${t}"
  echo "  ${date}: $(basename "${ocn_file}") -> ${out}"
done
[[ ${check_only} == true ]] && exit 0
submit_or_run "${script_dir}/make_ocn_ic.bash" "${out_root}/logs"

# ---------------------------------------------------------------------------
# Template restart with the CFSR state appended over it (the CFSR xtime is left
# out; the template's time is kept as xtime.orig).
# ---------------------------------------------------------------------------
for t in "${tasks[@]}"; do
  IFS='|' read -r date tmpl ocn_file out <<< "${t}"
  dir=$(dirname "${out}")
  mkdir -p "${dir}"
  tmp="${dir}/.$(basename "${out}").tmp"
  echo "$(date '+%H:%M:%S') ${date}: $(basename "${ocn_file}") on the ${template_year} mpaso restart"
  stage_restart "${tmpl}" "${tmp}" || exit 11
  ncks -A -C -x -v xtime "${ocn_file}" "${tmp}" || exit 12
  ncatted -O -a history_of_appended_files,global,d,, \
    -a rtenso2026_ocean,global,o,c,"CFSR ocean ${ocn_file} on the ${template_year} BruteForce mpaso restart" \
    "${tmp}" || exit 13
  mv -f "${tmp}" "${out}"
  echo "${date}_00:00:00" > "${dir}/rpointer.ocn"
  cat > "${dir}/ocean.json" <<EOF
{
  "target_date": "${date}",
  "analysis": "${ocn_file}",
  "template": "${tmpl}",
  "template_year": ${template_year},
  "output": "${out}",
  "created": "$(date -Iseconds)"
}
EOF
done
echo "$(date '+%H:%M:%S') done: ${#tasks[@]} dates"
