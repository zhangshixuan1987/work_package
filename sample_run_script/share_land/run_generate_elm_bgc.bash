#!/bin/bash
# Running on chrysalis
#SBATCH  --job-name=elm.interp
#SBATCH  --account=priority
#SBATCH  --nodes=1
#SBATCH  --output=elm.interp.o%j
#SBATCH  --exclusive
#SBATCH  --time=24:00:00
#SBATCH  --partition=priority


module purge
module load intel
module load netcdf-c
module load netcdf-fortran
module load hdf5

interp_cmd="./interpinic/interpinic"

elmbgc_source="./0501-01-01-00000/v3.LR.piControl.elm.r.0501-01-01-00000.nc"
elmi_source="./1980-01-02-00000/v3.LR.1980.I.cld_srt_frn.elm.r.1980-01-02-00000.nc"

workdir="/lcrc/group/e3sm/ac.szhang/acme_scratch/e3sm_project/gen_land_init"

if [ ! -d "${workdir}" ]; then
   mkdir -p ${workdir}
fi

cd ${workdir}

# step1: interpolate from v3.picontrol to cld_srt_frn
elmi_interp="picontrol_to_frn_elmi_interp.nc"
if [ ! -f "${elmi_interp}" ]; then
  cp -rp ${elmi_source} ${elmi_interp}
  ${interp_cmd} -i "${elmbgc_source}" -o "${elmi_interp}"
fi 

# step2: copy variables from v3.picontrol to result file (not addressed by interpolate)
replace_list=($(cat var_replace_list))
vars=$(echo ${replace_list[@]} | sed "s/ /,/g")
elmi_merge="picontrol_to_elmi_replace.nc"
if [ ! -f "${elmi_merge}" ]; then
  cp -rp ${elmi_interp} ${elmi_merge}
  ncks -A -v ${vars} ${elmbgc_source} ${elmi_merge}
fi

# step3: copy variables from v3.picontrol to result file (not addressed by interpolate)
snow_list=($(cat var_snow_list))
vsnows=$(echo ${snow_list[@]} | sed "s/ /,/g")
elmi_final="target_to_elmi_replace.nc"
if [ ! -f "${elmi_final}" ]; then
  cp -rp ${elmi_merge} ${elmi_final} 
  ncks -A -v ${vsnows} ${elmi_source} ${elmi_final}
fi

# step4: create final merged restart file
elmi_target="v3.LR.frn.merge.elm.r.1980-01-01-00000.nc"
if [ ! -f "${elmi_target}" ]; then
  cp -r ${elmi_final} ${elmi_target}
  # fix time manager fields
  ncap2 -O -s 'timemgr_rst_start_ymd=19800101' ${elmi_target} ${elmi_target}
  ncap2 -O -s 'timemgr_rst_curr_ymd=19800101'  ${elmi_target} ${elmi_target}
fi 
