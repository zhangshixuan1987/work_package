#!/bin/csh
# Submit this script as : sbatch ./[script-name]
#SBATCH  --job-name=regrid0300
#SBATCH  --nodes=1
#SBATCH  --time=12:00:00
#SBATCH  --exclusive
#SBATCH -A condo
#SBATCH -p acme-small

#module load ncl 
#module load nco 

set ystr   = 2008
set yend   = 2017
set tag    = `printf "%04d" $ystr`-`printf "%04d" $yend`

set jobdir = "/global/cfs/cdirs/e3sm/zhan391/e3sm_cvdp/run_script/run_polar_cvdp_sig/sig_ndg"
cd $jobdir/${tag}

sed -i "s/ystr      =.*/ystr      = $ystr/g" process_var_namelist.ncl
sed -i "s/yend      =.*/yend      = $yend/g" process_var_namelist.ncl

set outname = "cvdp_e3sm_amip_nudging_${tag}"
set cmd = "run_global_cvdp.ncl"
sed -i 's/outnam            =.*/outnam            = "'${outname}'"/g' $cmd

ncl $cmd >& diag.log

wait

exit


