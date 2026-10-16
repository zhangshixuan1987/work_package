#!/bin/bash -el
#SBATCH --account=gen0010
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --job-name=s2d_bfb_fosirl
#SBATCH --array=0-1
#SBATCH --output=slurm-bfb-%A_%a.out

set -o pipefail

test_date="2015-11-01"
work_dir=${SLURM_SUBMIT_DIR:-${PWD}}
reference_root="/global/cfs/cdirs/e3sm/zhan391/v3.LR.S2D.ENSINT"
test_root="/pscratch/sd/z/zhan391/e3sm_project/bfb_s2d_init/job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

case ${SLURM_ARRAY_TASK_ID} in
  0)
    test_mode="bruteforce-cfsr"
    run_case="BruteForce"
    ocean_option="LukeVanRoekel2"
    ;;
  1)
    test_mode="jra55-fosirl"
    run_case="JRA55-FOSIRL"
    ocean_option=""
    ;;
  *)
    echo "BFB ERROR: unsupported array task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 2
    ;;
esac

reference_dir="${reference_root}/${run_case}/${test_date}-00000"
test_dir="${test_root}/${run_case}/${test_date}-00000"

die() {
  echo "BFB ERROR: $*" >&2
  exit 1
}

[[ -d ${reference_dir} ]] || die "reference directory not found: ${reference_dir}"
reference_files=("${reference_dir}"/*)
[[ ${#reference_files[@]} -eq 20 ]] || die "expected 20 reference files, found ${#reference_files[@]}"
[[ ! -e ${test_root} ]] || die "isolated test root already exists: ${test_root}"

cd "${work_dir}"
export my_casedate="${test_date}"
export my_runcase="${run_case}"
export my_fosopt="JRA55"
export my_eamopt="ERA5-Tang2"
export my_ocnopt="${ocean_option}"
export my_runpath="${test_root}"
export my_overwrite=false

./gen_ensemble_s2d_init.sh

[[ -d ${test_dir} ]] || die "test output directory not created: ${test_dir}"
report="${test_dir}/BFB_REPORT.txt"
matches=0
mismatches=0
missing=0
extras=0

{
  echo "BFB mode: ${test_mode}"
  echo "BFB reference: ${reference_dir}"
  echo "BFB candidate: ${test_dir}"
  echo "BFB date: ${test_date}"
  echo
} > "${report}"

for reference_file in "${reference_dir}"/*; do
  name=${reference_file##*/}
  candidate_file="${test_dir}/${name}"
  if [[ ! -f ${candidate_file} ]]; then
    echo "MISSING ${name}" >> "${report}"
    missing=$((missing + 1))
  elif cmp -s "${reference_file}" "${candidate_file}"; then
    echo "MATCH ${name}" >> "${report}"
    matches=$((matches + 1))
  else
    echo "MISMATCH ${name}" >> "${report}"
    mismatches=$((mismatches + 1))
  fi
done

for candidate_file in "${test_dir}"/*; do
  name=${candidate_file##*/}
  [[ ${name} == BFB_REPORT.txt ]] && continue
  if [[ ! -f ${reference_dir}/${name} ]]; then
    echo "EXTRA ${name}" >> "${report}"
    extras=$((extras + 1))
  fi
done

{
  echo
  echo "matches=${matches}"
  echo "mismatches=${mismatches}"
  echo "missing=${missing}"
  echo "extras=${extras}"
} >> "${report}"

cat "${report}"
if (( mismatches > 0 || missing > 0 || extras > 0 || matches != 20 )); then
  die "strict BFB comparison failed; see ${report}"
fi

echo "BFB PASS: all 20 files are byte-for-byte identical."
