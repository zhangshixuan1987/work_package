#!/bin/bash -el

#===========================================
# Batch system directives
#===========================================
#SBATCH --account=esmd
#SBATCH --time=24:00:00
#SBATCH --partition=slurm
#SBATCH --job-name=e3sm_compress_diag
#SBATCH --nodes=1
#SBATCH --output=e3sm_compress_diag.%j
#SBATCH --exclusive
#SBATCH --no-kill
#SBATCH --requeue

# Define date range
START_YMD="2012-01-01-00000"
END_YMD="2012-03-01-00000"

# Submit main compression script
./5_run_compress_data.sh "$START_YMD" "$END_YMD"

# Ensure all background jobs complete
wait

exit 0
