#!/bin/bash 
# Running on pm-cpu
#SBATCH  --job-name=post_process
#SBATCH  --account=e3sm
#SBATCH  --nodes=4
#SBATCH  --output=post_process.o%j
#SBATCH  --exclusive
#SBATCH  --time=24:00:00
#SBATCH  --qos=regular
#SBATCH  --constraint=cpu

#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
#source /global/common/software/e3sm/anaconda_envs/load_e3sm_unified_1.10.0_pm-cpu.sh
package="/global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh"

www="/global/cfs/cdirs/e3sm/www/zhan391/E3SMv3_S2D"

run_dir="${RUN_DIR:-/global/cfs/cdirs/e3sm/S2S2D/simulation}"   # RUN_DIR=$PSCRATCH/e3sm_project/E3SMv3_S2D before stage 1
#run_dir="/pscratch/sd/z/zhan391/e3sm_project/E3SMv3_S2D" 
#run_dir="/pscratch/sd/t/tang30/E3SM_simulations/E3SMv3_S2D"

out_dir="/global/cfs/cdirs/e3sm/S2S2D/post_process"
#out_dir="/pscratch/sd/z/zhan391/e3sm_project/E3SMv3_S2D"

script_src_dir="/global/cfs/cdirs/e3sm/zhan391/s2d_diag/scripts"
plevdata="/global/cfs/cdirs/e3sm/diagnostics/e3sm_to_cmip_data/grids/vrt_remap_plev19.nc"
map="/global/cfs/cdirs/e3sm/diagnostics/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc"

ic_method="${IC_METHOD:-clim}"   # clim | obs | hybrid_ic
exp="WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_RTENSO2026_${ic_method}"

length=2   # years
max_jobs="${MAX_JOBS:-20}"

group="ts"
comp="atm"
grid="180x360_aave"
freq="monthly"
ckey="eam"
hist="eam.h0"

#vars=AODDUST,AODVIS,AODABS,AODABSBC,AODBC,AODSO4,AODSOA,AODSS,ABURDENSO4,ABURDENBC,ABURDENDUST,ABURDENPOM,ABURDENSEASALT,ABURDENSOA,AREL,AREI,FREQL,FREQI,Mass_bc,Mass_dst,Mass_ncl,Mass_pom,Mass_so4,Mass_soa,Mass_bc_srf,Mass_dst_srf,Mass_ncl_srf,Mass_pom_srf,Mass_so4_srf,Mass_soa_srf,dst_a1SF,dst_a3SF,SFSO2,SO2_CLXF,SO2,O3,O3_SRF,PBLH,CAPE,CIN,TH7001000,CLDICE,CLDLIQ,CLD_CAL_TMPLIQ,CLD_CAL_TMPICE,CLDLOW,CLDMED,CLDHGH,CLDTOT,CLWMODIS,CLIMODIS,CLDTOT_ISCCP,MEANCLDALB_ISCCP,MEANPTOP_ISCCP,CLD_CAL,CLDTOT_CAL_LIQ,CLDTOT_CAL_ICE,CLDTOT_CAL_UN,CLDHGH_CAL_LIQ,CLDHGH_CAL_ICE,CLDHGH_CAL_UN,CLDMED_CAL_LIQ,CLDMED_CAL_ICE,CLDMED_CAL_UN,CLDLOW_CAL_LIQ,CLDLOW_CAL_ICE,CLDLOW_CAL_UN,CLDHGH_CAL,CLDLOW_CAL,CLDMED_CAL,CLD_MISR,CLDTOT_CAL,CLMODIS,FISCCP1_COSP,FLDS,FLDSC,FLNS,FLNSC,FLNT,FLNTC,FLUT,FLUTC,SOLIN,SOLL,SOLS,FSDS,FSDSC,FSNS,FSNSC,FSNT,FSNTC,FSNTOA,FSNTOAC,SWCF,LWCF,ICEFRAC,LANDFRAC,OCNFRAC,PRECT,PRECTMX:X,PRECC,PRECL,PRECSC,PRECSL,PS,PSL,PHIS,U,V,T,Q,Z3,RELHUM,OMEGA,CLOUD,QFLX,QREFHT,RHREFHT,SHFLX,LHFLX,TAUX,TAUY,TCO,SCO,CDNUMC,TGCLDIWP,TGCLDCWP,TGCLDLWP,TMQ,TREFHT,TREFMNAV,TREFMXAV,TREFHTMN:M,TREFHTMX:X,TS,U10,UBOT,VBOT,TBOT,QBOT,U90M,V90M,TTQ,TUQ,TVQ,ZBOT,UTOP,TTOP,TROPF_P,TROPF_T,TROPF_Z,TOZ

vars=ICEFRAC,LANDFRAC,OCNFRAC,PSL,FSNTC,FSNTOAC,SWCF,LWCF,FLUT,FSNT,FSNTOA,FLNT,FLNTC,FSNS,FLNS,FSNS,SHFLX,QFLX,LHFLX,TAUX,TAUY,PRECC,PRECL,PRECSC,PRECSL,TS,TREFHT,U10,QREFHT,TMQ,CLDTOT,CLDHGH,CLDMED,CLDLOW,FLDS,FSDS,TGCLDIWP,TGCLDCWP,TGCLDLWP,FLNSC,FLUTC,FSDSC,SOLIN,FSNSC,AODBC,AODDUST,AODSO4,AODSOA,AODSS,AODVIS,TREFMNAV,TREFMXAV,TCO,PS,PHIS,AREL,U,V,T,Q,Z3,OMEGA,RELHUM
#TROP_P,CO_SRF,AODPOM,AODALL

interp_vars=AREL,U,V,T,Q,Z3,OMEGA,RELHUM
derive_vars=(PRECT  PRECST)
derive_var1=(PRECL  PRECSL)
derive_var2=(PRECC  PRECSC)

var_xtr= #LANDFRAC,OCNFRAC,ICEFRAC

leads=(2015020100 2015050100 2015080100 2015110100)

for lead in "${leads[@]}"
do 
  case="${exp}_${lead}"

  # Parse lead
  yr_srt=${lead:0:4}
  mth_srt=${lead:4:2}
  # Construct clean base date (avoid ambiguity)
  base="${yr_srt}-${mth_srt}-01"
  # Add years first, then subtract 1 month
  end=$(date -d "${base} +${length} years" +"%Y-%m-01")
  end=$(date -d "${end} -1 month" +"%Y %m")
  yr_end=$(echo "$end" | awk '{print $1}')
  mth_end=$(echo "$end" | awk '{print $2}')

  echo "yr_srt=$yr_srt"
  echo "mth_srt=$mth_srt"
  echo "yr_end=$yr_end"
  echo "mth_end=$mth_end"

  # number of complete 12-month periods
  start_total=$(( yr_srt * 12 + 10#$mth_srt ))
  end_total=$(( yr_end * 12 + 10#$mth_end ))
  n_months=$(( end_total - start_total + 1 ))
  ypf=$(( n_months / 12 ))

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
    sed -i "s|^package=.*|package=${package}|g" "$out_script"
    sed -i "s|^diri=.*|diri=${run_dir}|g" "$out_script"
    sed -i "s|^diro=.*|diro=${out_dir}|g" "$out_script"
    sed -i "s|^case=.*|case=${case}|g" "$out_script"
    sed -i "s|^ens=.*|ens=${ens}|g" "$out_script"
    sed -i "s|^vars=.*|vars=${vars}|g" "$out_script"
    sed -i "s|^interp_vars=.*|interp_vars=${interp_vars}|g" "$out_script"
    sed -i "s|^derive_vars=.*|derive_vars=(${derive_vars[*]})|g" "$out_script"
    sed -i "s|^derive_var1=.*|derive_var1=(${derive_var1[*]})|g" "$out_script"
    sed -i "s|^derive_var2=.*|derive_var2=(${derive_var2[*]})|g" "$out_script"
    sed -i "s|^comp=.*|comp=${comp}|g" "$out_script"
    sed -i "s|^ckey=.*|ckey=${ckey}|g" "$out_script"
    sed -i "s|^hist=.*|hist=${hist}|g" "$out_script"
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
done

wait
exit 0
