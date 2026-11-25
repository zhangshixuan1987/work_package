#!/bin/bash
#SBATCH --job-name=strat9
#SBATCH --account=gen0010
#SBATCH --qos=regular
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --output=/pscratch/sd/z/zhan391/run_logs/rtenso2026_notebook_runs/strat_series.o%j
# Execute the stratosphere notebooks jupyter/stratosphere/01-08 in order (each builds the caches it needs;
# the first one that runs reads the native eam.h1 output)
# (executed copies and logs to RUN_DIR on pscratch; notebooks on disk unchanged).
# Submit from anywhere: sbatch scripts/run_strat_series.sbatch (it cds to diagnostic/ itself)
cd /global/cfs/cdirs/m4849/zhan391/work_package/e3sm_s2d_proc/diagnostic
source /global/homes/z/zhan391/.conda/envs/e3sm_analysis/bin/activate 2>/dev/null || true
PY=/global/homes/z/zhan391/.conda/envs/e3sm_analysis/bin
RUN_DIR=/pscratch/sd/z/zhan391/run_logs/rtenso2026_notebook_runs/strat_series_${SLURM_JOB_ID:-manual}
mkdir -p $RUN_DIR
for nb in jupyter/stratosphere/01_model_climatology jupyter/stratosphere/02_vortex_forecast \
          jupyter/stratosphere/03_ssw_occurrence jupyter/stratosphere/04_enso_qbo_background jupyter/stratosphere/05_ssw_dynamics \
          jupyter/stratosphere/06_ssw_mechanism jupyter/stratosphere/07_vortex_ssw_animation jupyter/stratosphere/08_vortex_evolution_animation; do
    t0=$(date +%s)
    tag=${nb#jupyter/}; tag=${tag//\//__}   # jupyter/stratosphere/01_x -> stratosphere__01_x
    $PY/jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=-1 --ExecutePreprocessor.kernel_name=e3sm_analysis \
        --output-dir=$RUN_DIR --output=$tag $nb.ipynb > $RUN_DIR/$tag.log 2>&1
    rc=$?
    echo "[$(date +%T)] $nb exit $rc ($(( $(date +%s) - t0 )) s)"
done
