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
set hpssdir = "/home/z/zhan391/DARPA_ML_Train_Data/FV_180x360"
cd $workdir

foreach dir (after_nudging  before_nudging  clim  nudging_tendency  reference_data)
 echo $dir
 cd $workdir/$dir
 foreach file (*.nc)  
  echo $file 
  hsi put $hpssdir/$dir/$file 
  #htar -cvf $hpssdir/$dir/$file.tar $file
 end 
end

