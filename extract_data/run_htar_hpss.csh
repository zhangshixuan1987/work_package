#!/bin/csh
# Submit this script as : sbatch ./[script-name]
#SBATCH -A m3525
#SBATCH -q regular
#SBATCH -t 48:00:00
#SBATCH -N 2
#SBATCH  --job-name=ncclimo_ctrl
#SBATCH  --output=job%j 
#SBATCH  --exclusive 
#SBATCH  --constraint=knl,quad,cache

set workdir = `pwd`
set hpssdir = "/home/z/zhan391/E3SM_Cryosphere"
cd $workdir

#foreach dir (E3SM*DT1800* E3SM*DT0300*)
foreach dir (ERA5 ERA5_NDG ERA_20C NCEP_ANL NOAA-CIRES-DOE_20C) #(20200410* *alpha5*) 

echo $dir

htar -cvf $hpssdir/$dir.tar $dir

wait

end

