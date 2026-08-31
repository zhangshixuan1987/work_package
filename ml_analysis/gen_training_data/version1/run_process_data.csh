#!/bin/csh
##SBATCH  --job-name=proc_training
#SBATCH  --account=e3sm
#SBATCH  --nodes=1
#SBATCH  --output=process.o%j
#SBATCH  --exclusive
#SBATCH  --time=12:00:00
#SBATCH  --qos=regular
#SBATCH  --constraint=cpu

source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.csh

python extract_data_set1.py 2014 2017 & 
python extract_data_set2.py 2014 2017 &
python extract_data_set3.py 2014 2017 &
python extract_data_set4.py 2014 2017 &

wait

exit

