#!/bin/csh
# Submit this script as : sbatch ./[script-name]
#SBATCH -N 1
#SBATCH -t 48:00:00
#SBATCH -q regular
#SBATCH -A m3089
#SBATCH --constraint=knl

module load nco 
module load cdo
conda activate e3sm_analysis

set Expnames = ("CLIM" "NDGUV" "NDGUVT" "NDGUVTQ")
set nexps    = $#Expnames

set workdir  = "/global/cscratch1/sd/zhan391/DARPA_project/Track_model"
set codedir  = "$workdir/code"
set outdir   = "${workdir}/out"

if ( ! -d $outdir ) then
  mkdir -p $outdir
endif

set iexp = 4
set expnam = $Expnames[$iexp]

set year = 2007
while ( $year <= 2017 )  

  cd $outdir
  cp $codedir/sythetic_track_0001.py  ${expnam}_${year}_fcst.py
  
  sed -i "s/0001/${year}/g"          ${expnam}_${year}_fcst.py
  sed -i "s/CLIM/${expnam}/g"        ${expnam}_${year}_fcst.py

  python ${expnam}_${year}_fcst.py > ${expnam}_${year}.log & 

 @ year++ 
end 

wait

exit 

