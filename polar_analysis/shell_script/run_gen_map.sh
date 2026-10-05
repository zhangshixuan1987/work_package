#!/bin/bash
# Submit this script as : sbatch ./[script-name]
#SBATCH  --job-name=regrid0350
#SBATCH  --nodes=1
#SBATCH  --time=12:00:00
#SBATCH  --exclusive
#SBATCH -A condo
#SBATCH -p acme-small

# Parameters: override from the environment for another machine/experiment,
# e.g.  DATA_DIR=/other/data bash run_gen_map.sh   (or: sbatch --export=ALL,DATA_DIR=... run_gen_map.sh)
DATA_DIR="${DATA_DIR:-/lcrc/group/e3sm/public_html/diagnostic_output/ac.szhang/v3_polar_analysis/data}"

jobdir=`pwd`
cd $jobdir

ref="/lcrc/soft/climate/e3sm_diags_data/obs_for_e3sm_diags/climatology/ERA5/ERA5_01_197901_201901_climo.nc"
in_grid="${DATA_DIR}/map/ERA5.g"
out_grid="/lcrc/group/e3sm/ac.szhang/acme_scratch/data/regrid_maps/cmip6_180x360_scrip.20181001.nc"
outmap="${DATA_DIR}/map/map_721x1440_to_180x360_conserve.nc"
overlap="${DATA_DIR}/map/era5_721x1440.overlap.g"

ref="/lcrc/soft/climate/e3sm_diags_data/obs_for_e3sm_diags/climatology/MERRA2/MERRA2_01_198001_201601_climo.nc"
in_grid="${DATA_DIR}/map/MERRA2.g"
out_grid="/lcrc/group/e3sm/ac.szhang/acme_scratch/data/regrid_maps/cmip6_180x360_scrip.20181001.nc"
outmap="${DATA_DIR}/map/map_361x576_to_180x360_conserve.nc"
overlap="${DATA_DIR}/map/MERRA2_361x576.overlap.g"

#GenerateRLLMesh --lon 1440 --lat 721 --lat_begin 90 --lat_end -90 --file ERA5.g
GenerateRLLMesh --in_file $ref --in_file_lon lon --in_file_lat lat --file ${in_grid}

GenerateOverlapMesh --a ${in_grid} --b ${out_grid} --out ${overlap}

GenerateOfflineMap --in_mesh ${in_grid} --out_mesh ${out_grid} \
                   --ov_mesh ${overlap} --in_type fv --out_type fv \
                   --out_map ${outmap} 

