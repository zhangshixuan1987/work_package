#!/bin/bash -l 
#SBATCH -A phy220062      # Allocation name 
#SBATCH --nodes=1         # Total # of nodes (must be 1 for serial job)
#SBATCH --ntasks=1        # Total # of MPI tasks (should be 1 for serial job) 
#SBATCH --time=8:00:00    # Total run time limit (hh:mm:ss)
#SBATCH -J data_process   # Job name
#SBATCH -o myjob.o%j      # Name of stdout output file
#SBATCH -e myjob.e%j      # Name of stderr error file
#SBATCH -p shared         # Queue (partition) name

conda activate e3sm_analysis

for file in 1_* 2_* 3_* 4*; do

echo $file 
ncl $file &

done 

wait

exit
