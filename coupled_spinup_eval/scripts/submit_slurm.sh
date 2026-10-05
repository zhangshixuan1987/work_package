#!/bin/bash
# Execute analysis notebooks headlessly as a batch job on Chrysalis, using the setup
# defined in each notebook. Submit from coupled_spinup_eval/ (or scripts/):
#
#   sbatch scripts/submit_slurm.sh jupyter/14_ocn_ts_ctrl_analysis.ipynb
#   sbatch scripts/submit_slurm.sh jupyter/1_enso_compare_analysis.ipynb jupyter/4_ohc_compare_analysis.ipynb
#   sbatch scripts/submit_slurm.sh            # all notebooks in jupyter/, in numeric order
#
# Executed copies (with outputs) go to jupyter/executed/<name>.<jobid>.ipynb.
#
#SBATCH --job-name=spinup_eval
#SBATCH --account=e3sm
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=06:00:00
#SBATCH --output=slurm-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
# sbatch runs a spooled copy of this file, so locate the workflow from the submit dir.
if [ -d jupyter ] && [ -d util ]; then ROOT=$(pwd)
elif [ -d ../jupyter ] && [ -d ../util ]; then ROOT=$(cd .. && pwd)
else echo "submit from coupled_spinup_eval/ or its scripts/ directory" >&2; exit 1
fi
PYBIN="$(dirname "${SPINUP_EVAL_PYTHON:-/lcrc/soft/climate/e3sm-unified/e3smu_1_12_0/chrysalis/conda/envs/e3sm_unified_1.12.0_login/bin/python}")"
export HDF5_USE_FILE_LOCKING=FALSE

if [ $# -eq 0 ]; then
    mapfile -t NBS < <(cd "$ROOT/jupyter" && ls [0-9]*_*.ipynb | sort -V | sed "s#^#$ROOT/jupyter/#")
    set -- "${NBS[@]}"
fi
mkdir -p "$ROOT/jupyter/executed"
status=0
for arg in "$@"; do
    NB="$(cd "$(dirname "$arg")" && pwd)/$(basename "$arg")"
    echo "=== $(date +%T) $NB"
    ( cd "$ROOT/jupyter" &&       # notebooks expect cwd = jupyter/
      "$PYBIN/jupyter" nbconvert --to notebook --execute \
          --ExecutePreprocessor.kernel_name=python3 --ExecutePreprocessor.timeout=-1 \
          "$NB" --output-dir "$ROOT/jupyter/executed" \
          --output "$(basename "$NB" .ipynb).${SLURM_JOB_ID:-local}.ipynb" ) || { echo "FAILED: $NB"; status=1; }
done
exit $status
