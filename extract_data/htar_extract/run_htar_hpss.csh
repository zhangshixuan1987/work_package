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
set hpssdir = "/home/g/golaz/2020/20201211.beta1_01.piControl.compy"
set outdir  = "/global/cfs/cdirs/e3sm/zhan391/data/htar_extract"
cd $workdir


set years = 11
set yeare = 40 

set iy = $years 
while ( $iy < $yeare)
 set file = 0000${iy}.tar
 htar -xvf $hpssdir/$file $outdir/$iy 

 @ iy ++ 

end 

