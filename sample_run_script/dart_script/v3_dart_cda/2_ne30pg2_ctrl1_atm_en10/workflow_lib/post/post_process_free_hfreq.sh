#!/bin/bash -el
# Free-run 6-hourly, daily and monthly series for one member (free_run_hfreq.py).
# Internal worker of 7_run_post_hist.sh: POST_PRODUCT selects eam_hfreq (atm)
# or elm_hfreq (lnd); POST_START_DATE..POST_END_DATE are inclusive days.
set -o pipefail
if [[ "${POSTPROCESS_DRIVER_ACTIVE:-FALSE}" != "TRUE" ]]; then
  echo "ERROR: internal post-processing worker; run 7_run_post_hist.sh" >&2
  exit 1
fi
fail() { echo "ERROR: $*" >&2; exit 1; }
: "${OMP_NUM_THREADS:=1}"; export OMP_NUM_THREADS

case "${POST_PRODUCT:-}" in
  eam_hfreq) component=atm; map_file="${my_eam_hist_mapping_file}" ;;
  elm_hfreq) component=lnd; map_file="${my_elm_hist_mapping_file}" ;;
  *) fail "POST_PRODUCT must be eam_hfreq or elm_hfreq, got: ${POST_PRODUCT:-unset}" ;;
esac

if [[ "${my_analysis_env_type:-}" == "custom" ]]; then
  source "${my_conda_setup_file}" || fail "could not source Conda setup file: ${my_conda_setup_file}"
  conda activate "${my_analysis_conda_env}" || fail "could not activate Conda environment: ${my_analysis_conda_env}"
else
  my_analysis_pixi_dir="/share/apps/E3SM/conda_envs/e3smu_1_13_0/compy/pixi_login"
  [[ -x "${my_analysis_pixi_dir}/bin/pixi" ]] || fail "Pixi executable not found in ${my_analysis_pixi_dir}"
  eval "$("${my_analysis_pixi_dir}/bin/pixi" shell-hook --as-is --shell bash \
    --manifest-path "${my_analysis_pixi_dir}/pixi.toml")"
  hash -r
fi
command -v ncremap >/dev/null 2>&1 || fail "ncremap not found"

overwrite=()
[[ "${OVERWRITE_EXISTING:-FALSE}" == "TRUE" ]] && overwrite=(--overwrite)
end_date=$(date -u -d "${POST_END_DATE} +1 day" +%Y-%m-%d)

echo "== ${POST_PRODUCT} for ${POST_CASE_NAME}: ${POST_START_DATE} to ${POST_END_DATE} =="
python -I "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/free_run_hfreq.py" \
  --component "${component}" \
  --archive "${my_modeldir}/${POST_ENSTR}/archive" \
  --case-name "${POST_CASE_NAME}" --member "${POST_ENSTR}" \
  --start "${POST_START_DATE}" --end "${end_date}" \
  --map-file "${map_file}" --interval-hours "${my_e3sm_cycle_hours}" \
  --manifest "${POST_MANIFEST:?POST_MANIFEST is unset}" "${overwrite[@]}"
