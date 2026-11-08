#!/bin/bash
#SBATCH --job-name=strat9
#SBATCH --account=gen0010
#SBATCH --qos=regular
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --output=.split_work/run/strat9.o%j
# Execute the 9a-9i stratosphere notebooks in order (outputs to .split_work/run/; notebooks on disk unchanged).
cd /global/cfs/cdirs/m4849/zhan391/work_package/e3sm_s2d_proc/diagnostic
source /global/homes/z/zhan391/.conda/envs/e3sm_analysis/bin/activate 2>/dev/null || true
PY=/global/homes/z/zhan391/.conda/envs/e3sm_analysis/bin
for nb in 9e_enso_qbo_background 9g_ssw_mechanism; do
    t0=$(date +%s)
    $PY/jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=-1 --ExecutePreprocessor.kernel_name=e3sm_analysis \
        --output-dir=.split_work/run $nb.ipynb > .split_work/run/$nb.log 2>&1
    rc=$?
    echo "[$(date +%T)] $nb exit $rc ($(( $(date +%s) - t0 )) s)"
done
