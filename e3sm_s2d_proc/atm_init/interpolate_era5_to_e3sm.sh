#!/bin/bash -l

conda activate e3sm_analysis
#outdir="/pscratch/sd/z/zhan391/DARPA_project/SCREAM_INIT/ATM"
#betacast="/pscratch/sd/z/zhan391/DARPA_project/e3sm_model/code/betacast"

betacast=$1
nam_lev=$2
DATE=$3
era5=$4
rgd_map=$5
topfil=$6
outdir=$7
outfil=$8

mod_in_topo="${topfil}"
mod_out_topo="${topfil}"
if [ ! -f ${outdir}/${outfil} ]; then 
  cd ${betacast}/atm_to_cam
  ncl -n atm_to_cam.ncl 'datasource="ERA5"'     \
         numlevels=${nam_lev} \
         YYYYMMDDHH=${DATE} \
         'data_filename = "'${era5}'"'  \
         'wgt_filename="'${rgd_map}'"' \
         'model_topo_file="'${mod_out_topo}'"'\
         'adjust_config="a"' \
         'se_inic = "'${outdir}/${outfil}'"'
fi 
