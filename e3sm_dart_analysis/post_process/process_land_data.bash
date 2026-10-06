#!/bin/bash
ulimit -S -c 0
#SBATCH --job-name=elm_regrid
#SBATCH --account=e3sm
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --output=regrid_%j.log
#SBATCH --error=regrid_%j.err

# ---- Load E3SM unified env ----
source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh

set -euo pipefail

data_dir="/compyfs/zhan391/acme_init/E3SMv3_INT"
case_name="20241109.v3-LR.DATM.ne30pg2_r05_IcoswISC30E3r5.pm-cpu"
time_str="2011-12-01-00000"
map_file="/compyfs/zhan391/acme_init/map_file/map_r05_to_cmip6_180x360_aave.20200901.nc"
domain_file="/compyfs/inputdata/share/domains/domain.lnd.r05_IcoswISC30E3r5.231121.nc"
target_file="1x1d.nc "

vars="cols1d_ityp,cols1d_active,cols1d_ityplun,cols1d_landunit_index,cols1d_topounit_index,cols1d_gridcell_index,cols1d_jxy,cols1d_ixy,cols1d_lon,cols1d_lat,T_GRND,TH2OSFC,TH2OSFC,SNOW_DEPTH,H2OSFC,H2OSNO,H2OSOI_LIQ,H2OSOI_ICE"

file="${data_dir}/${time_str}/${case_name}.elm.r.${time_str}.nc"
ncks -O -v "${vars}" $file "elm.r.${time_str}.nc"

# ---- USER INPUTS ----
RESTART="elm.r.${time_str}.nc"
OUT_AGG="agg_gridcell.nc"               # columns→gridcells output
RLL_LON=720                             # 0.5° lon count
RLL_LAT=360                             # 0.5° lat count
RLL_G="rll_${RLL_LAT}x${RLL_LON}.g"     # destination mesh filename
MAP_OUT="map_src_to_rll_conserve.nc"    # conservative weight file
REGRID_OUT="elm_180x360_aave.nc"  # final regridded fields

src_grid_name="r05"
src_scrip_file="./grid/SCRIPgrid_0.5x0.5_nomask_c110308.nc"
dst_grid_name="cmip6_180x360"
dst_scrip_file="./grid/cmip6_180x360_scrip.20181001.nc"
date="c231227"

# ##############################################
# Area-average maps for domain files
#   (ocn -> lnd; not sure what this is used for)
# ##############################################
alg_name=aave
function map_aave {
    if [ ! -e ${map} ]; then
        echo "map: ${map}"
        cmd="ncremap -a aave -s ${src_grid} -g ${dst_grid} -m ${map}"
        echo "${cmd}" && ${cmd}
    fi
}
# ocn to atm
src_grid=${src_scrip_file}
dst_grid=${dst_scrip_file}
map=./grid/map_${src_grid_name}_to_${dst_grid_name}_${alg_name}.${date}.nc
map_aave

# ---- 0) sanity: tools needed ----
need() { command -v "$1" >/dev/null 2>&1 || { echo "Missing $1"; exit 1; }; }
need python
need ncremap
need GenerateRLLMesh
#need ConvertSCRIPToMesh
need GenerateOverlapMesh
need GenerateOfflineMap

# ---- 1) Columns → Gridcells (NetCDF with lon/lat [+lon_b/lat_b if available]) ----
#python column_to_gridcell_exporter.py \
#  --restart "${RESTART}" \
#  --src-domain "${domain_file}" \
#  --out "${OUT_AGG}"

# ---- 5) Apply weights to the aggregated file ----
ncremap -R '--rgr lat_nm_in=gridcell' -m "${map}" -i "${OUT_AGG}" -o "${REGRID_OUT}"

echo "Done."
echo "  Aggregated gridcells: ${OUT_AGG}"
echo "  Destination mesh:     ${RLL_G}"
echo "  Weights (conserve):   ${MAP_OUT}"
echo "  Regridded output:     ${REGRID_OUT}"

exit

# ---- 2) Destination regular lat–lon mesh (.g) at 0.5° ----
# Use 0..360 longitudes and -90..90 latitudes; adjust if you prefer -180..180
GenerateRLLMesh \
  --lon "${RLL_LON}" --lat "${RLL_LAT}" \
  --lonbegin 0 --lonend 360 \
  --latbegin -90 --latend 90 \
  --file "${RLL_G}"

# ---- 3) Source grid: derive SCRIP from aggregated file, then convert to mesh (.g) ----
# ncremap can emit a SCRIP source grid it infers from the file's lon/lat[/bounds]
SRC_SCRIP="src_grid_scrip.nc"
ncremap -s "${SRC_SCRIP}" -i "${OUT_AGG}" -o /dev/null

# convert the SCRIP grid into a TempestRemap mesh
SRC_MESH="src_grid.g"
ConvertSCRIPToMesh --in "${SRC_SCRIP}" --out "${SRC_MESH}"

# ---- 4) Build conservative weights with TempestRemap ----
OVERLAP="overlap_src_to_rll.g"
GenerateOverlapMesh --a "${SRC_MESH}" --b "${RLL_G}" --out "${OVERLAP}"

GenerateOfflineMap \
  --a "${SRC_MESH}" \
  --b "${RLL_G}" \
  --f "${OVERLAP}" \
  --out "${MAP_OUT}" \
  --method conserve

# ---- 5) Apply weights to the aggregated file ----
ncremap -m "${MAP_OUT}" -i "${OUT_AGG}" -o "${REGRID_OUT}"

echo "Done."
echo "  Aggregated gridcells: ${OUT_AGG}"
echo "  Destination mesh:     ${RLL_G}"
echo "  Weights (conserve):   ${MAP_OUT}"
echo "  Regridded output:     ${REGRID_OUT}"

