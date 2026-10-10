#!/bin/csh
# Submit this script as : sbatch ./[script-name]
#SBATCH -A esmd
#SBATCH -q slurm
#SBATCH -t 12:00:00
#SBATCH -N 5
#SBATCH  --job-name=nudge

setenv OMP_NUM_THREADS 1

set workdir = `pwd`

cd $workdir


set ii = 1

foreach exp ( CLIM UV3_PL UVT3_PL UVTQ3_PL UVTQ3_BASE \
              UVTQ3_CTRL UVTQ3_GFDL UVTQ3_NPBL UVTQ3_NPBL_NUPP \
              UVT3_NPBL_NUPP_OPT0 UVT3_NPBL_NUPP_OPT1 UVT3_NPBL_NUPP_OPT2 \
              UVTQ3_NPBL_UPP_OPT0 UVTQ3_NPBL_UPP_OPT1 UVTQ3_NPBL_UPP_OPT2 \
              UVTQ3_NPBL_UPP_OPT3 )
 
 set file = run_process_${exp}.csh 
 echo $file 

 csh $file > csh`printf "%02d" $ii`.log & 
 @ ii++ 
end 

wait

