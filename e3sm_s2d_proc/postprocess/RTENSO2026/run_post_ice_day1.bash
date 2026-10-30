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
run_dir="${RUN_DIR:-/global/cfs/cdirs/e3sm/S2S2D/RTENSO2026/simulation}"   # RUN_DIR=$PSCRATCH/e3sm_project/E3SMv3_S2D before stage 1
#run_dir="/pscratch/sd/t/tang30/E3SM_simulations/E3SMv3_S2D"

#out_dir="/global/cfs/cdirs/e3sm/S2S2D/post_process"
out_dir="/global/cfs/cdirs/e3sm/S2S2D/RTENSO2026/post_process"

script_src_dir="/global/cfs/cdirs/e3sm/zhan391/s2d_diag/scripts"

plevdata="/global/cfs/cdirs/e3sm/diagnostics/e3sm_to_cmip_data/grids/vrt_remap_plev19.nc"
map="/global/cfs/cdirs/e3sm/diagnostics/maps/map_IcoswISC30E3r5_to_cmip6_180x360_traave.20240221.nc"

ic_method="${IC_METHOD:-clim}"   # clim | obs | hybrid_ic
exp="WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_RTENSO2026_${ic_method}"
length=2   # months
max_jobs="${MAX_JOBS:-10}"
group="ts"
comp="ice"
grid="180x360_aave"
freq="daily"
ckey="mpasseaice"
hist="mpassi.hist.am.timeSeriesStatsDaily"
ref_yr=1850
avg_interval_days=1
dpf=30
tpd=1

vars=timeDaily_avg_iceAreaCell,timeDaily_avg_iceVolumeCell,timeDaily_avg_snowVolumeCell
vext=xtime_startDaily,xtime_endDaily

leads=(2015020100 2015050100 2015080100 2015110100)

for lead in "${leads[@]}"
#for lead in 1980050100 1980110100 1981050100 1981110100
#for lead in 1982050100 1982110100 1983050100 1983110100
#for lead in 1984050100 1984110100 1985050100 1985110100
#for lead in 1986050100 1986110100 1987050100 1987110100
#for lead in 1988050100 1988110100 1989050100 1989110100
#for lead in 1990050100 1990110100 1991050100 1991110100
#for lead in 1992050100 1992110100 1993050100 1993110100
#for lead in 1994050100 1994110100 1995050100 1995110100
#for lead in 1996050100 1996110100 1997050100 1997110100
#for lead in 1998050100 1998110100 1999050100 1999110100
#for lead in 2000050100 2000110100 2001050100 2001110100
#for lead in 2002050100 2002110100 2003050100 2003110100
#for lead in 2004050100 2004110100 2005050100 2005110100
#for lead in 2006050100 2006110100 2007050100 2007110100
#for lead in 2008050100 2008110100 2009050100 2009110100
#for lead in 2010050100 2010110100 2011050100 2011110100
#for lead in 2012050100 2012110100 2013050100 2013110100
#for lead in 2015050100 2015110100 2016050100 2016110100
#for lead in 2017050100 2017110100 2018050100 2018110100
#for lead in 2014050100 2014110100

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
    ttag="$(printf "%04d" "${yr_srt}")$(printf "%02d" "$((10#${mth_srt}))")-$(printf "%04d" "${yr_end}")$(printf "%02d" "$((10#${mth_end}))")-$(printf "%04d" "${ypf}")"

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
    sed -i "s|^mth_srt=.*|mth_srt=$((10#${mth_srt}))|g" "$out_script"
    sed -i "s|^mth_end=.*|mth_end=$((10#${mth_end}))|g" "$out_script"
    sed -i "s|^ref_yr=.*|ref_yr=${ref_yr}|g" "$out_script"
    sed -i "s|^avg_interval_days=.*|avg_interval_days=${avg_interval_days}|g" "$out_script"
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

    # ncremap -p mpi needs a Slurm node list; run it serially on login nodes
    if [[ -z "${SLURM_JOB_ID:-}" ]]; then
      sed -i "s|ncremap -p mpi |ncremap |g" "$out_script"
    fi

    # JOB_NBR caps ncclimo splitter parallelism (default: one job per variable)
    if [[ -n "${JOB_NBR:-}" ]]; then
      sed -i "s|^--split \\\\$|--split --job_nbr=${JOB_NBR} \\\\|" "$out_script"
    fi

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
