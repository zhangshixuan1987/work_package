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

#source /global/cfs/cdirs/e3sm/software/anaconda_envs/load_latest_e3sm_unified.csh
#source /global/cfs/cdirs/e3sm/software/anaconda_envs/load_latest_e3sm_unified.sh

set workdir = `pwd`
set case    = "20201211.beta1_01.piControl.compy"
set hpssdir = "/home/g/golaz/2020/$case"
set outdir  = "/global/cfs/cdirs/e3sm/zhan391/data/htar_extract"
cd $workdir

set years = 11
set yeare = 40 

set filelist = ""

set iy = $years 

while ( $iy <= $yeare )
 
 set im = 1
 while ( $im <= 12 )

 set yst = `printf "%04d" $iy`
 set mst = `printf "%02d" $im`

 set file = "$case.eam.h0.${yst}-${mst}.nc"
 set filelist = ( $filelist $file )
 @ im ++
 end 

 @ iy++ 

end 

echo $filelist

zstash extract --hpss=$hpssdir archive/atm/hist/*eam.h0.* 
