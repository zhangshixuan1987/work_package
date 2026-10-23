#!/bin/bash 
# Running on pm-cpu
#SBATCH  --job-name=post_process
#SBATCH  --account=e3sm
#SBATCH  --nodes=1
#SBATCH  --output=post_process.o%j
#SBATCH  --exclusive
#SBATCH  --time=24:00:00
#SBATCH  --qos=regular
#SBATCH  --constraint=cpu

source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
package="/global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh"

www="/global/cfs/cdirs/e3sm/www/zhan391/E3SMv3_S2D"

#run_dir="/global/cfs/cdirs/e3sm/S2S2D/simulation" 
run_dir="/global/cfs/cdirs/e3sm/S2S2D/simulation"
#run_dir="/pscratch/sd/t/tang30/E3SM_simulations/E3SMv3_S2D"

#out_dir="/global/cfs/cdirs/e3sm/S2S2D/post_process"
out_dir="/global/cfs/cdirs/e3sm/S2S2D/post_process"

script_src_dir="/global/cfs/cdirs/e3sm/zhan391/s2d_diag/scripts"

plevdata="/global/cfs/cdirs/e3sm/diagnostics/e3sm_to_cmip_data/grids/vrt_remap_plev19.nc"
map="/global/cfs/cdirs/e3sm/diagnostics/maps/map_IcoswISC30E3r5_to_cmip6_180x360_traave.20240221.nc"

exp="WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_4DEnVarOcn"
length=2   # months
max_jobs=10
group="ts"
comp="ocn"
grid="180x360_aave"
freq="daily"
ckey="mpasocean"
hist="mpaso.hist.am.highFrequencyOutput"
dpf=30
tpd=1

#vars=kineticEnergyAtSurface,kineticEnergyAt250m,relativeVorticityAt250m,ssh,pressureAdjustedSSH,boundaryLayerDepth,dThreshMLD,tThreshMLD,columnIntegratedSpeed,barotropicSpeed,atmosphericPressure,temperatureAtSurface,salinityAtSurface,temperatureAt250m,salinityAt250m,temperatureAtBottom,salinityAtBottom
#vars=temperatureAtSurface,salinityAtSurface,temperatureAt250m,salinityAt250m,ssh,pressureAdjustedSSH,dThreshMLD,boundaryLayerDepth,columnIntegratedSpeed,barotropicSpeed,atmosphericPressure
#vars=salinityAtSurface,ssh,temperatureAtSurface,salinityAtSurface,temperatureAt250m,salinityAt250m,ssh,pressureAdjustedSSH,dThreshMLD,boundaryLayerDepth,columnIntegratedSpeed,barotropicSpeed,atmosphericPressure
vars=temperatureAtSurface,temperatureAt250m,dThreshMLD,salinityAtSurface,salinityAt250m,ssh
vext=xtime

leads=()
for yr in $(seq 1980 2018); do
  for mth in 05 11; do
    leads+=("${yr}${mth}0100")
  done
done

#for lead in 2014050100 2014110100
for lead in "${leads[@]}"
do 
  case="${exp}_${lead}"

  # Parse lead
  yr_srt=${lead:0:4}
  mth_srt=${lead:4:2}
  if ! [[ "${length}" =~ ^[0-9]+$ ]] || (( length < 1 )); then
    echo "ERROR: length must be a positive integer number of months: ${length}" >&2
    exit 1
  fi
  # Construct clean base date (avoid ambiguity)
  base="${yr_srt}-${mth_srt}-01"
  # length is inclusive: 1 means the lead month only.
  end=$(date -d "${base} +$((length - 1)) months" +"%Y %m")
  yr_end=$(echo "$end" | awk '{print $1}')
  mth_end=$(echo "$end" | awk '{print $2}')

  echo "yr_srt=$yr_srt"
  echo "mth_srt=$mth_srt"
  echo "yr_end=$yr_end"
  echo "mth_end=$mth_end"

  # number of complete 12-month periods
  start_total=$(( yr_srt * 12 + 10#${mth_srt} ))
  end_total=$(( yr_end * 12 + 10#${mth_end} ))
  n_months=$(( end_total - start_total + 1 ))
  ypf=$(( n_months / 12 ))

  pids=()

  for en in $(seq 0 9); do
    ens="EN$(printf "%02d" "$en")"
    ttag="$(printf "%04d" "${yr_srt}")$(printf "%02d" "${mth_srt}")-$(printf "%04d" "${yr_end}")$(printf "%02d" "${mth_end}")-$(printf "%04d" "${ypf}")"

    script_dir="${out_dir}/${case}/${ens}/post/scripts"
    mkdir -p "${script_dir}" || exit 1

    # Skip only when the status is OK and the expected outputs are present.
    status_file="${out_dir}/${case}/${ens}/post/scripts/${group}_${comp}_${freq}_${grid}_${ttag}.status"
    tag_srt="${yr_srt}${mth_srt}"
    tag_end="${yr_end}${mth_end}"
    if [[ "${group}" == clim* ]]; then
      output_dir="${out_dir}/${case}/${ens}/post/${comp}/${grid}/${group}/${ypf}yr"
    else
      output_dir="${out_dir}/${case}/${ens}/post/${comp}/${grid}/${group}/${freq}/${ypf}yr"
    fi

    outputs_complete=false
    if [[ -f "${status_file}" ]] &&
       [[ "$(tr -d '[:space:]' < "${status_file}")" == "OK" ]] &&
       [[ -d "${output_dir}" ]]; then
      outputs_complete=true
      if [[ "${group}" == clim* || -z "${vars:-}" ]]; then
        # Climatology filenames are generated dynamically by ncclimo.
        output_files=("${output_dir}"/*.nc)
        for output_file in "${output_files[@]}"; do
          if [[ ! -s "${output_file}" ]]; then
            outputs_complete=false
            break
          fi
        done
      else
        IFS=',' read -r -a requested_vars <<< "${vars}"
        for var in "${requested_vars[@]}"; do
          var="${var//[[:space:]]/}"
          [[ -z "${var}" ]] && continue
          standard_output="${output_dir}/${var}_${tag_srt}_${tag_end}.nc"
          tagged_output="${output_dir}/${var}_${ttag}.nc"
          if [[ ! -s "${standard_output}" && ! -s "${tagged_output}" ]]; then
            outputs_complete=false
            break
          fi
        done
      fi
    fi

    if [[ "${outputs_complete}" == true ]]; then
      echo "skipping completed output: ${case} ${ens} ${tag_srt}-${tag_end}"
      continue
    fi


    # move the post directory to separate directory
    src_dir="${run_dir}/${case}/${ens}/post/${comp}"
    tgt_dir="${out_dir}/${case}/${ens}/post/${comp}"
    if [[ -d "${src_dir}" && "${src_dir}" != "${tgt_dir}" ]]; then
        # ensure parent exists
        mkdir -p "$(dirname "${tgt_dir}")"

        # remove target if it already exists (optional, but safer for clean move)
        if [[ -d "${tgt_dir}" ]]; then
            rm -rf "${tgt_dir}"
        fi

        mv "${src_dir}" "${tgt_dir}"
        continue
    fi

    echo "${script_src_dir}/${group}_${comp}_${freq}_${grid}.bash"

    template_script="${script_src_dir}/${group}_${comp}_${freq}_${grid}.bash"
    if [[ ! -f "${template_script}" ]]; then
      echo "ERROR: missing template script: ${template_script}" >&2
      exit 1
    fi

    base=$(basename "$template_script" .bash)
    fname="${base}_${ens}_${ttag}.bash"
    out_script="${script_dir}/${fname}"

    cp "$template_script" "$out_script" || exit 1

    sed -i "s|^yr_srt=.*|yr_srt=${yr_srt}|g" "$out_script"
    sed -i "s|^yr_end=.*|yr_end=${yr_end}|g" "$out_script"
    sed -i "s|^mth_srt=.*|mth_srt=${mth_srt}|g" "$out_script"
    sed -i "s|^mth_end=.*|mth_end=${mth_end}|g" "$out_script"
    sed -i "s|^package=.*|package=${package}|g" "$out_script"
    sed -i "s|^diri=.*|diri=${run_dir}|g" "$out_script"
    sed -i "s|^diro=.*|diro=${out_dir}|g" "$out_script"
    sed -i "s|^case=.*|case=${case}|g" "$out_script"
    sed -i "s|^ens=.*|ens=${ens}|g" "$out_script"
    sed -i "s|^vars=.*|vars=${vars}|g" "$out_script"
    sed -i "s|^vext=.*|vext=${vext}|g" "$out_script"
    sed -i "s|^comp=.*|comp=${comp}|g" "$out_script"
    sed -i "s|^ckey=.*|ckey=${ckey}|g" "$out_script"
    sed -i "s|^hist=.*|hist=${hist}|g" "$out_script"
    sed -i "s|^dpf=.*|dpf=${dpf}|g" "$out_script"
    sed -i "s|^tpd=.*|tpd=${tpd}|g" "$out_script"
    sed -i "s|^plevdata=.*|plevdata=${plevdata}|g" "$out_script"
    sed -i "s|^map=.*|map=${map}|g" "$out_script"

    chmod 755 "$out_script" || exit 1
    echo "working on $out_script"
    bash "$out_script" &

    while (( $(jobs -r | wc -l) >= max_jobs )); do
      wait
    done

  done
  # Wait for remaining jobs
  wait
done

wait
exit 0
