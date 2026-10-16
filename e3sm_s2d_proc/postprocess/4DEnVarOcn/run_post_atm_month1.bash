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

#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
#source /global/common/software/e3sm/anaconda_envs/load_e3sm_unified_1.10.0_pm-cpu.sh
package="/global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh"

www="/global/cfs/cdirs/e3sm/www/zhan391/E3SMv3_S2D"

run_dir="/global/cfs/cdirs/e3sm/S2S2D/simulation" 
#run_dir="/pscratch/sd/z/zhan391/e3sm_project/E3SMv3_S2D" 
#run_dir="/pscratch/sd/t/tang30/E3SM_simulations/E3SMv3_S2D"

out_dir="/global/cfs/cdirs/e3sm/S2S2D/post_process"
#out_dir="/pscratch/sd/z/zhan391/e3sm_project/E3SMv3_S2D"

script_src_dir="/global/cfs/cdirs/e3sm/zhan391/s2d_diag/scripts"
plevdata="/global/cfs/cdirs/e3sm/diagnostics/e3sm_to_cmip_data/grids/vrt_remap_plev19.nc"
map="/global/cfs/cdirs/e3sm/diagnostics/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc"

exp="WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_4DEnVarOcn"

length=2 # years
max_jobs=4

if [[ ! "$max_jobs" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: max_jobs must be a positive integer" >&2
  exit 1
fi

# Track children across all lead dates (compatible with Bash 4.4).
pids=()
job_failed=0
wait_oldest_job() {
  local pid=${pids[0]} rc
  if wait "$pid"; then
    :
  else
    rc=$?
    echo "ERROR: child process $pid exited with status $rc" >&2
    job_failed=1
  fi
  pids=("${pids[@]:1}")
}

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
for lead in 2008050100 2008110100 2009050100 2009110100
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
  start_total=$(( yr_srt * 12 + mth_srt ))
  end_total=$(( yr_end * 12 + mth_end ))
  n_months=$(( end_total - start_total + 1 ))
  ypf=$(( n_months / 12 ))

  for en in $(seq 0 9); do
    ens="EN$(printf "%02d" "$en")"
    ttag="$(printf "%04d" "${yr_srt}")$(printf "%02d" "${mth_srt}")-$(printf "%04d" "${yr_end}")$(printf "%02d" "${mth_end}")-$(printf "%04d" "${ypf}")"

    script_dir="${out_dir}/${case}/${ens}/post/scripts"
    mkdir -p "${script_dir}" || exit 1

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
    sed -i "s|^interp_vars=.*|interp_vars=${interp_vars}|g" "$out_script"
    sed -i "s|^derive_vars=.*|derive_vars=(${derive_vars[*]})|g" "$out_script"
    sed -i "s|^derive_var1=.*|derive_var1=(${derive_var1[*]})|g" "$out_script"
    sed -i "s|^derive_var2=.*|derive_var2=(${derive_var2[*]})|g" "$out_script"
    sed -i "s|^comp=.*|comp=${comp}|g" "$out_script"
    sed -i "s|^ckey=.*|ckey=${ckey}|g" "$out_script"
    sed -i "s|^hist=.*|hist=${hist}|g" "$out_script"
    sed -i "s|^plevdata=.*|plevdata=${plevdata}|g" "$out_script"
    sed -i "s|^map=.*|map=${map}|g" "$out_script"

    chmod 755 "$out_script" || exit 1
    echo "working on $out_script"
    bash "$out_script" &
    pids+=("$!")
    if (( ${#pids[@]} >= max_jobs )); then
      wait_oldest_job
    fi
  done
done

# Collect every remaining child's exit status.
while (( ${#pids[@]} > 0 )); do
  wait_oldest_job
done
exit "$job_failed"
