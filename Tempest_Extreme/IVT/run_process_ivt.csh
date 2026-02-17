#!/bin/csh
# Submit this script as : sbatch ./[script-name]
#SBATCH -N 4
#SBATCH -C amd
#SBATCH -t 24:00:00
#SBATCH -q regular #bigmem
#SBATCH -A m3089

module load ncl

setenv OMP_NUM_THREADS 1

set workdir = `pwd`

cd $workdir

ncl 1_process_IVT_e3sm_clim.ncl &

ncl 1_process_IVT_e3sm_ndguv3_tau06.ncl & 

ncl 1_process_IVT_e3sm_ndguvt3_tau06.ncl &

ncl 1_process_IVT_e3sm_ndguvtq3_tau06.ncl &
   
wait 

