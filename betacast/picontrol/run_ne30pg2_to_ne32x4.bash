#!/bin/bash -l

conda activate e3sm_analysis

mesh_dir=`pwd`
in_mesh="ne30"
out_mesh="northamericax4v1"
mlev=80

inputdir="/pscratch/sd/z/zhan391/v3_narrm"
date="0101-01-01-00000"
atm_in="${inputdir}/v3.LR.piControl.eam.i.${date}.nc"

#method with betacast 
rgd_map="${mesh_dir}/map_${in_mesh}_to_1x1d_map.nc"
if [ ! -f ${rgd_map} ];then
  GenerateRLLMesh --lon 360 --lat 180 --lat_begin -90 --lat_end 90 --file ${mesh_dir}/1x1d.g
  GenerateOverlapMesh --a ${mesh_dir}/ne30.g --b ${mesh_dir}/1x1d.g --out ${mesh_dir}/${in_mesh}_1x1d.overlap.g
  GenerateOfflineMap  --in_mesh ${mesh_dir}/ne30.g --out_mesh ${mesh_dir}/1x1d.g \
                      --ov_mesh ${mesh_dir}/${in_mesh}_1x1d.overlap.g --in_type cgll --out_type fv \
                      --in_np 4 --out_np 1 --out_map ${rgd_map} 
fi

mod_remap_file=${mesh_dir}/map_1x1d_to_${out_mesh}_map.nc
if [ ! -f ${mod_remap_file} ];then
  #generate map to lat-lon interp
  #GenerateRLLMesh --in_file ${mesh_dir}/1x1d.nc --in_file_lon grid_center_lon --in_file_lat grid_center_lat --file ${mesh_dir}/1x1d.g 
  #GenerateRLLMesh --lon 360 --lat 180 --lat_begin -90 --lat_end 90 --file ${mesh_dir}/1x1d.g
  GenerateOverlapMesh --a ${mesh_dir}/1x1d.g --b ${mesh_dir}/${out_mesh}.g --out ${mesh_dir}/1x1d_${out_mesh}.overlap.g
  GenerateOfflineMap  --in_mesh ${mesh_dir}/1x1d.g --out_mesh ${mesh_dir}/${out_mesh}.g \
                      --ov_mesh ${mesh_dir}/1x1d_${out_mesh}.overlap.g --in_type fv --out_type cgll \
                      --in_np 1 --out_np 4 --out_map ${mod_remap_file} 
fi

betacast="/pscratch/sd/z/zhan391/DARPA_project/e3sm_model/code/betacast/atm_to_cam"
cd ${betacast}
#echo $date
#echo ${date:0:4}${date:5:2}${date:8:2}${hour}
cyclestr=${date:11:5}
hour=`expr ${cyclestr} / 3600`
hour=`printf "%02d" $hour`
YYYYMMDDHH="${date:0:4}${date:5:2}${date:8:2}${hour}"
echo $YYYYMMDDHH

outfil="${inputdir}/v3.LR.piControl.${out_mesh}.eam.i.${date}.nc"
mod_in_topo="${mesh_dir}/USGS-gtopo30_ne30np4pg2_x6t-SGH.c20210614.nc"
mod_out_topo="${mesh_dir}/USGS_northamericax4v1pg2_12xdel2_consistentSGH_20020209.nc"

echo $betacast
ncl -n atm_to_cam.ncl 'datasource="CAM"'     \
       numlevels=${mlev} \
       YYYYMMDDHH="${YYYYMMDDHH}" \
       'data_filename = "'${atm_in}'"'  \
       'wgt_filename="'${mod_remap_file}'"' \
       'model_topo_file="'${mod_out_topo}'"'\
       'adjust_config="a"' \
       'mod_in_topo="'${mod_in_topo}'"' \
       'mod_remap_file="'${rgd_map}'"'\
       'se_inic = "'${outfil}'"'

