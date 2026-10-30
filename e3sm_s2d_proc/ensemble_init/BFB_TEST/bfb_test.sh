#!/bin/bash -el
#SBATCH --account=gen0010
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --job-name=s2d_bfb
#SBATCH --array=0-2
#SBATCH --output=slurm-bfb-%A_%a.out
#
# Regenerate one established initial-condition date per set and compare it
# with the archived result. Submit from this directory:
#
#   sbatch bfb_test.sh               # all three sets
#   sbatch --array=1 bfb_test.sh     # JRA55-FOSIRL only
#
# This test never writes to the established result directory.

set -o pipefail

bfb_dir=${SLURM_SUBMIT_DIR:-${PWD}}
reference_root="/global/cfs/cdirs/e3sm/zhan391/v3.LR.S2D.ENSINT"
test_root="/pscratch/sd/z/zhan391/e3sm_project/bfb_s2d_init/job_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

# set directory | run case | EAM option | ocean option | file count | default date
case ${SLURM_ARRAY_TASK_ID} in
  0) set_dir="BruteForce/luke_cfsr"; run_case="BruteForce";   eam_option="ERA5-Tang2"; ocean_option="LukeVanRoekel2"; file_count=20; default_date="2015-11-01" ;;
  1) set_dir="JRA55_FOSIRL";         run_case="JRA55-FOSIRL"; eam_option="ERA5-Tang2"; ocean_option="";               file_count=20; default_date="2015-11-01" ;;
  2) set_dir="4DEnVarOcn";        run_case="4DEnVarOcn";   eam_option="4DEnVarOcn"; ocean_option="4DEnVarOcn";     file_count=23; default_date="1995-05-01" ;;
  *)
    echo "BFB ERROR: unsupported array task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 2
    ;;
esac

test_date=${BFB_DATE:-${default_date}}
reference_dir="${reference_root}/${run_case}/${test_date}-00000"
test_dir="${test_root}/${run_case}/${test_date}-00000"
generator_dir="${bfb_dir}/../${set_dir}"

die() {
  echo "BFB ERROR: $*" >&2
  exit 1
}

[[ -f ${bfb_dir}/bfb_compare.sh ]] || die "submit from the BFB_TEST directory (SLURM_SUBMIT_DIR=${bfb_dir})"
[[ -x ${generator_dir}/gen_ensemble_s2d_init.sh ]] || die "generator not found: ${generator_dir}/gen_ensemble_s2d_init.sh"
[[ -d ${reference_dir} ]] || die "reference directory not found: ${reference_dir}"
[[ ! -e ${test_root} ]] || die "isolated test root already exists: ${test_root}"

echo "BFB set: ${set_dir} (${run_case}), date ${test_date}"

cd "${generator_dir}"
export my_casedate="${test_date}"
export my_runcase="${run_case}"
export my_fosopt="JRA55"
export my_eamopt="${eam_option}"
export my_ocnopt="${ocean_option}"
export my_runpath="${test_root}"
export my_overwrite=false

./gen_ensemble_s2d_init.sh

[[ -d ${test_dir} ]] || die "test output directory not created: ${test_dir}"
bash "${bfb_dir}/bfb_compare.sh" "${reference_dir}" "${test_dir}" "${file_count}"
