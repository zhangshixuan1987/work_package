#!/bin/bash

# Parameters: override from the environment for another machine/experiment,
# e.g.  DATA_DIR=/other/data bash derive.bash   (or: sbatch --export=ALL,DATA_DIR=... derive.bash)
DATA_DIR="${DATA_DIR:-/lcrc/group/e3sm/public_html/diagnostic_output/ac.szhang/v3_polar_analysis/data}"

outdir="${DATA_DIR}/climo"
if [ ! -d ${outdir} ];then 
  mkdir -p ${outdir}
fi 

data_dir="/lcrc/group/acme/ac.szhang/acme_scratch/data/merra2/rgd_data"
for outfile in ${outdir}/MERRA.analysis.climo*;do 
  ncap2 -s "FSUS=FSDS-FSNS"            ${outfile} ${outfile}.tmp
  mv ${outfile}.tmp ${outfile}
  ncap2 -s "FLUS=FLDS-FLNS"            ${outfile} ${outfile}.tmp
  mv ${outfile}.tmp ${outfile}
  ncap2 -s "FLUSC=FLDSC-FLNSC"         ${outfile} ${outfile}.tmp
  mv ${outfile}.tmp ${outfile}
  ncap2 -s "FLUSC=FSDSC-FSNSC"         ${outfile} ${outfile}.tmp
  mv ${outfile}.tmp ${outfile}
done
