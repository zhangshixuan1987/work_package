#!/bin/bash

# Parameters: override from the environment for another machine/experiment,
# e.g.  DATA_DIR=/other/data bash run_process_derived.bash   (or: sbatch --export=ALL,DATA_DIR=... run_process_derived.bash)
DATA_DIR="${DATA_DIR:-/lcrc/group/e3sm/public_html/diagnostic_output/ac.szhang/v3_polar_analysis/data}"

map="${DATA_DIR}/map/map_721x1440_to_180x360_conserve.nc"
outdir="${DATA_DIR}/climo"
if [ ! -d ${outdir} ];then 
  mkdir -p ${outdir}
fi 

data_dir="/lcrc/group/e3sm/ac.szhang/acme_scratch/e3sm_project/test_zppy_pmp"
sub_dir="post/atm/180x360_aave/clim"
for outfile in $outdir/v3.LR*;do
    echo $outfile
    ncap2 -s "FLDSC=FLDS+FLNS-FLNSC"  ${outfile} ${outfile}.tmp
    mv ${outfile}.tmp ${outfile}
done 
