#!/bin/bash -l
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --account=e3sm
#SBATCH  --constraint=cpu
#SBATCH  --qos=regular
#SBATCH  --nodes=1
#SBATCH  --time=02:00:00
#SBATCH  --job-name=rtenso_atm_ic
#SBATCH  --output=/global/cfs/cdirs/e3sm/S2S2D/RTENSO2026/initial_conditions/logs/atm_ic_%j.out

set -o pipefail

# ==============================================================================
# Atmosphere (EAM) initial conditions for RTENSO2026: ensemble members EN00-EN09
# from the HICCUP ERA5 analysis (Tang), as in the BruteForce set. The
# atmosphere does not depend on the land/sea-ice method (clim, obs, hybrid_ic),
# so one set is written per date:
#
#   <out_root>/YYYYMMDD/atm/v3.LR.amip_0101.HICCUP.atm_era5.ENxx.eam.i.YYYY-MM-DD-00000.nc
#   <out_root>/YYYYMMDD/atm/atm.json
#
# Each date takes the first source that exists:
#   bruteforce  the BruteForce members of the same date, copied unchanged (2015)
#   hiccup      Tang's HICCUP file: reconstructed on the eam.i template (EN00)
#               and perturbed with DART for EN01-EN09 (util/eam_ens.sh)
#   pending     neither yet (2026 until Tang's files exist): reported, skipped
# The coupler restart is not written here: cpl/<method>/ gets the cpl.r of
# the land year of each method (util/clim_restart.py).
#
# Edit USER SETTINGS, then:
#   CHECK_ONLY=true ./make_atm_ic.bash   preflight only: sources, protected outputs
#   ./make_atm_ic.bash                   execution_mode below (default auto: run here on an
#                                        interactive node or in a job, submit from a login node)
#   ./make_atm_ic.bash --run             run here, e.g. on an interactive node
#                                        (salloc -N 1 -C cpu -q interactive -t 02:00:00 -A e3sm)
#   ./make_atm_ic.bash --submit          submit one Slurm job
# About 3 minutes per date to copy 26.5 GB; the DART path needs the allocation.
# ==============================================================================

# ---------------------------- USER SETTINGS ----------------------------------
out_root="/global/cfs/cdirs/e3sm/S2S2D/RTENSO2026/initial_conditions"
dates=(2015-02-01 2015-05-01 2015-08-01 2015-11-01
       2026-02-01 2026-05-01 2026-08-01 2026-11-01)
out_subdir="atm"                      # written to <out_root>/YYYYMMDD/<out_subdir>/
nmem=9                                # members EN00-EN<nmem>; EN00 unperturbed

# BruteForce members of the same date (read-only).
src_root="/global/cfs/cdirs/e3sm/zhan391/v3.LR.S2D.ENSINT/BruteForce"
atm_case="v3.LR.amip_0101.HICCUP.atm_era5"   # my_refcase1 in the hindcast driver

# HICCUP ERA5 files (Tang), named by calendar date (YYYY-MM-DD).
hiccup_dir="/pscratch/sd/t/tang30/share/ForShixuan/ERA5_IC"
hiccup_case="${atm_case}.%s.00.ne30np4.L80.nc"

overwrite_existing=false              # replace existing members
check_only="${CHECK_ONLY:-false}"     # true: preflight only
# auto: run here inside a Slurm allocation (interactive node or batch job), submit from a
# login node; run: always run here; submit: always submit a Slurm job. --run/--submit override it.
execution_mode="${EXECUTION_MODE:-auto}"
# -----------------------------------------------------------------------------

script_dir=${IC_WORK_DIR:-$(cd "$(dirname "$0")" && pwd)}
export IC_WORK_DIR=${script_dir}
source "${script_dir}/util/ic_lib.sh"
source "${script_dir}/util/eam_ens.sh"
parse_execution_mode "$@"
if ! command -v ncks >/dev/null; then
  source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh >/dev/null 2>&1
fi

member_name() { printf "%s.EN%02d.eam.i.%s-00000.nc" "${atm_case}" "$1" "$2"; }

# ---------------------------------------------------------------------------
# Preflight: a source for every date (or pending), outputs not protected.
# ---------------------------------------------------------------------------
problems=0
tasks=()
pending=()
need_dart=false
for date in "${dates[@]}"; do
  bf_dir="${src_root}/${date}-00000"
  hiccup="${hiccup_dir}/$(printf "${hiccup_case}" "${date}")"
  out_dir="${out_root}/${date//-/}/${out_subdir}"
  source_kind=bruteforce
  for i in $(seq 0 "${nmem}"); do
    [[ -r ${bf_dir}/$(member_name "${i}" "${date}") ]] || { source_kind=; break; }
  done
  if [[ -z ${source_kind} ]]; then
    if [[ -r ${hiccup} ]]; then
      source_kind=hiccup; need_dart=true
    else
      pending+=("${date}"); continue
    fi
  fi
  if [[ -e ${out_dir}/$(member_name 0 "${date}") && ${overwrite_existing} != true ]]; then
    echo "EXISTS: ${out_dir}/$(member_name 0 "${date}") (set overwrite_existing=true to replace it)" >&2
    problems=$((problems + 1)); continue
  fi
  tasks+=("${date}|${source_kind}|${bf_dir}|${hiccup}|${out_dir}")
done

if [[ ${need_dart} == true ]]; then
  while read -r f; do
    [[ -r ${f} ]] || { echo "MISSING: ${f}" >&2; problems=$((problems + 1)); }
  done < <(eam_ens_inputs)
fi

if (( problems > 0 )); then
  echo "Preflight failed: ${problems} problem(s)." >&2
  exit 1
fi
for date in "${pending[@]}"; do
  echo "PENDING: ${date}: no BruteForce members and no $(printf "${hiccup_case}" "${date}") in ${hiccup_dir}"
done
if (( ${#tasks[@]} == 0 )); then
  echo "Nothing to do: every date is pending." >&2
  exit 1
fi
echo "Preflight passed: ${#tasks[@]} dates (EN00-EN$(printf "%02d" "${nmem}")), ${#pending[@]} pending."
for t in "${tasks[@]}"; do
  IFS='|' read -r date source_kind bf_dir hiccup out_dir <<< "${t}"
  [[ ${source_kind} == bruteforce ]] && from=${bf_dir} || from=$(basename "${hiccup}")
  echo "  ${date}: ${source_kind} ${from} -> ${out_dir}"
done
[[ ${check_only} == true ]] && exit 0
submit_or_run "${script_dir}/make_atm_ic.bash" "${out_root}/logs"

# ---------------------------------------------------------------------------
# Members built in <out_dir>/.build and moved into place when all are done.
# ---------------------------------------------------------------------------
for t in "${tasks[@]}"; do
  IFS='|' read -r date source_kind bf_dir hiccup out_dir <<< "${t}"
  build="${out_dir}/.build"
  rm -rf "${build}"; mkdir -p "${build}"
  echo "$(date '+%H:%M:%S') ${date}: ${source_kind}"
  if [[ ${source_kind} == bruteforce ]]; then
    for i in $(seq 0 "${nmem}"); do
      f=$(member_name "${i}" "${date}")
      cp --preserve=timestamps "${bf_dir}/${f}" "${build}/${f}" || exit 11
      [[ $(stat -c %s "${bf_dir}/${f}") == $(stat -c %s "${build}/${f}") ]] || {
        echo "SIZE: ${build}/${f} differs from the source" >&2; exit 12; }
    done
    analysis=${bf_dir}
  else
    eam_ens_members "${hiccup}" "${build}" "${atm_case}" "${date}-00000" "${nmem}" || exit 13
    analysis=${hiccup}
  fi
  for i in $(seq 0 "${nmem}"); do
    f=$(member_name "${i}" "${date}")
    mv -f "${build}/${f}" "${out_dir}/${f}"
  done
  rmdir "${build}"
  cat > "${out_dir}/atm.json" <<EOF
{
  "target_date": "${date}",
  "source": "${source_kind}",
  "analysis": "${analysis}",
  "case": "${atm_case}",
  "members": "EN00-EN$(printf "%02d" "${nmem}")",
  "created": "$(date -Iseconds)"
}
EOF
done
echo "$(date '+%H:%M:%S') done: ${#tasks[@]} dates"
