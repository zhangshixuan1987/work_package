#!/bin/bash
ulimit -S -c 0
source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh
src_grid_name="ne11np4"
src_scrip_file="./grid/${src_grid_name}_scrip.nc"
dst_grid_name="cmip6_180x360"
dst_scrip_file="./grid/cmip6_180x360_scrip.20181001.nc"
output_root=`pwd`
date="c231227"

data_dir="/compyfs/zhan391/acme_init/E3SMv3_INT"
case_name="20241109.v3-LR.DATM.ne30pg2_r05_IcoswISC30E3r5.pm-cpu"
time_str="2011-12-01-00000"
map_file="/compyfs/zhan391/acme_init/map_file/map_r05_to_cmip6_180x360_aave.20200901.nc"
domain_file="/compyfs/inputdata/share/domains/domain.lnd.r05_IcoswISC30E3r5.231121.nc"
target_file="1x1d.nc "

file="${data_dir}/${time_str}/${case_name}.elm.r.${time_str}.nc"
outfile="./${case_name}.elm.r.${time_str}.nc"

GenerateRLLMesh --in_file ${file}  --in_file_lon "column" --in_file_lat "column" --lon cols1d_lon --lat cols1d_lat --file ${mapdir}/src.g

#ncdump -h $file >log
exit

vars="cols1d_ityp,cols1d_active,cols1d_ityplun,cols1d_landunit_index,cols1d_topounit_index,cols1d_gridcell_index,cols1d_jxy,cols1d_ixy,cols1d_lon,cols1d_lat,T_GRND,TH2OSFC,TH2OSFC,SNOW_DEPTH,H2OSFC,H2OSNO,H2OSOI_LIQ,H2OSOI_ICE"


#ncremap -5 --alg_typ=aave --map=${map_file} --grd_src=${domain_file} --grd_dst=${target_file}

ncremap -P elm -m ${map_file} out.nc out_map.nc # Same as -P clm, alm, ctsm

exit

ncremap \
  --src_grd=${domain_file} \
  --dst_grd=fv1.0x1.0 \
  --map=map_lnd_to_gx0p5.nc \
  --alg_typ=bilinear
exit

#ncks -v "${vars}" $file out.nc
ncrename -v cols1d_lat,lat -v cols1d_lon,lon out.nc
ncremap -m ${map_file} out.nc out_map.nc

exit

interpinic \
  -i ${file} \
  -o ${outfile} \
  -map ${map_file}


