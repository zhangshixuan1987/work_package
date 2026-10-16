#!/bin/bash -el
#SBATCH --account=gen0010
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --job-name=s2d_bfb_cprnc
#SBATCH --array=0-1
#SBATCH --output=slurm-bfb-cprnc-%A_%a.out

set -o pipefail

cprnc="/global/cfs/cdirs/e3sm/tools/cprnc/cprnc"
date_tag="2015-11-01-00000"

case ${SLURM_ARRAY_TASK_ID} in
  0)
    run_case="BruteForce"
    candidate_root="/pscratch/sd/z/zhan391/e3sm_project/bfb_s2d_init/job_58599001_0"
    ;;
  1)
    run_case="JRA55-FOSIRL"
    candidate_root="/pscratch/sd/z/zhan391/e3sm_project/bfb_s2d_init/job_58599001_1"
    ;;
  *) exit 2 ;;
esac

reference_dir="/global/cfs/cdirs/e3sm/zhan391/v3.LR.S2D.ENSINT/${run_case}/${date_tag}"
candidate_dir="${candidate_root}/${run_case}/${date_tag}"
report="${candidate_dir}/CPRNC_REPORT.txt"
strict_matches=0
cprnc_matches=0
failures=0

[[ -x ${cprnc} ]]
: > "${report}"

for reference_file in "${reference_dir}"/*.nc; do
  name=${reference_file##*/}
  candidate_file="${candidate_dir}/${name}"
  if cmp -s "${reference_file}" "${candidate_file}"; then
    echo "STRICT_BFB_PASS ${name}" >> "${report}"
    strict_matches=$((strict_matches + 1))
    continue
  fi

  detail="${candidate_dir}/${name}.cprnc.txt"
  "${cprnc}" -m "${reference_file}" "${candidate_file}" > "${detail}" 2>&1 || true
  if grep -q 'diff_test: the two files seem to be IDENTICAL' "${detail}"; then
    echo "CPRNC_BFB_PASS ${name}" >> "${report}"
    cprnc_matches=$((cprnc_matches + 1))
    rm -f "${detail}"
  else
    echo "CPRNC_BFB_FAIL ${name}: ${detail}" >> "${report}"
    failures=$((failures + 1))
  fi
done

{
  echo "strict_matches=${strict_matches}"
  echo "cprnc_matches=${cprnc_matches}"
  echo "failures=${failures}"
} >> "${report}"

cat "${report}"
(( failures == 0 ))
