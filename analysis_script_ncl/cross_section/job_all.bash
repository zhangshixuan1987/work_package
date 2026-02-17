#!/bin/bash -l 
#SBATCH -A phy220062      # Allocation name 
#SBATCH --nodes=1         # Total # of nodes (must be 1 for serial job)
#SBATCH --ntasks=128      # Total # of MPI tasks (should be 1 for serial job) 
#SBATCH --time=8:00:00    # Total run time limit (hh:mm:ss)
#SBATCH -J data_process   # Job name
#SBATCH -o myjob.o%j      # Name of stdout output file
#SBATCH -e myjob.e%j      # Name of stderr error file
#SBATCH -p shared         # Queue (partition) name

conda activate e3sm_analysis

ncl 1_gen_data_e3sm.ncl & 
ncl 1_gen_data_era5.ncl &

wait

exit
