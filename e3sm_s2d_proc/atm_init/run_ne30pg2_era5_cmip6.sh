#!/bin/bash 
#SBATCH -A e3sm
#SBATCH -C cpu
#SBATCH -t 02:00:00
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --job-name=init
#SBATCH --output=interp.%j 

#conda activate e3sm_unified_latest #e3sm_analysis
source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

echo "Using python: $(which python)"
echo "Using ncrename: $(which ncrename)"

wkdir=${PWD}
betacast="/pscratch/sd/z/zhan391/SEACROGS_project/e3sm_model/betacast"
target_dir="/pscratch/sd/z/zhan391/e3sm_project/input/era5_v3LR_int"
ecmwf="/pscratch/sd/z/zhan391/e3sm_project/input/v3.LR.S2D/ECMWF"
mapdir="/pscratch/sd/z/zhan391/e3sm_project/input/v3.LR.S2D/maps"
topodir="/global/cfs/cdirs/e3sm/inputdata/atm/cam/topo"
topofil="${topodir}/USGS-gtopo30_ne30np4pg2_x6t-SGH.c20210614.nc"

meshes=("ne30pg2")
dates=("1980050100")
mod_nam="E3SM"
mod_ver="v3.LR"
mod_comp="eam"
mod_lev=80  #128 for scream

# first step: extract ERA5 reanalysis data 
if [ ! -d ${ecmwf} ];then 
  mkdir -p ${ecmwf}
fi 
for DATE in "${dates[@]}"; do
  cd  $wkdir
  if [ ! -f "${ecmwf}/ERA5_${DATE}.nc" ]; then
    ./prestage-ERA5_CMIP6.sh ${ecmwf} $DATE
  fi
done 

# second step: generate regrid map if not exist 
if [ ! -d ${mapdir} ];then
  mkdir -p ${mapdir}
fi
#conda activate e3sm_analysis
for reg in "${meshes[@]}"; do
  cd $wkdir
  res=${reg:0:4}
  rgd_map="${mapdir}/ERA5_to_${res}_map.nc"
  if [ ! -f "${rgd_map}" ]; then 
    # the standard atmosphere SE grid, using the np4" GLL spectral element
    # GenerateCSMesh --alt --res 30  --file NE30.g        
    # An atmosphere "FV" grid (for the new pg2 option where physics is on a FV grid)    
    # GenerateVolumetricMesh --in NE30.g --np 2 --out NE30pg2.g --uniform     
    #Note: For a cubed-sphere mesh (for ne30 only pass 30)
    GenerateCSMesh --res ${res:2:2} --alt --file ${res}.g
    
    #For ERA5 FV mesh 
    #GenerateRLLMesh --lon 360 --lat 180 --lat_begin -90 --lat_end 90 --file 1x1d.g
    #GenerateRLLMesh --lon 1440 --lat 721 --lat_begin 90 --lat_end -90 --file ERA5.g
    #inferred from file 
    GenerateRLLMesh --in_file_lon "longitude" --in_file_lat "latitude" \
	            --in_file "${ecmwf}/ERA5_${dates[0]}.nc" --file ERA5.g 
    
    #For overlap mesh
    GenerateOverlapMesh --a ERA5.g --b ${res}.g --out overlap.g
   
    #For mapping file: See https://github.com/ClimateGlobalChange/tempestremap
    GenerateOfflineMap  --in_mesh ERA5.g --out_mesh ${res}.g \
                        --ov_mesh overlap.g --in_type fv --out_type cgll \
			--in_np 4 --out_np 4 --out_map ${rgd_map}
  fi 
  rm -rvf *.g 
done 

# third step: generate initial condition file 
for tdate in "${dates[@]}"; do
  cd  $wkdir
  time_tag="${tdate:0:4}-${tdate:4:2}-${tdate:6:2}-00000"
  outdir="${target_dir}/${time_tag}"
  if [ ! -d ${outdir} ];then
     mkdir -p ${outdir}
  fi
  for reg in "${meshes[@]}"; do
    cd  $wkdir
    res=${reg:0:4}
    rgd_map="${mapdir}/ERA5_to_${res}_map.nc"
    era5="${ecmwf}/ERA5_${tdate}.nc"
    outfil="v3.LR.${reg}.ERA5.${mod_comp}.i.${time_tag}.nc"
    ./interpolate_era5_to_e3sm.sh ${betacast} ${mod_lev} ${tdate} ${era5} ${rgd_map} ${topofil} ${outdir} ${outfil}
  done 
done 
