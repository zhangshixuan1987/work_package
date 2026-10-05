#!/bin/bash
# Submit this script as : sbatch ./[script-name]
#SBATCH --job-name=regrid0300
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --exclusive
#SBATCH -A condo
#SBATCH -p acme-small

source /lcrc/soft/climate/e3sm-unified/load_e3sm_unified_1.9.3_chrysalis.sh

jobdir="$(pwd)"
cd "$jobdir"

# --- Variables (Bash arrays)
var2d_name=("SITIMEFRAC" "SICONC" "SITHICK" "SIMASS")
var2d_frac=(1.0 100.0 1.0 917.0)
var2d_unit=("1" "%" "m" "kg/m2")
var2d_list=(
  "timeMonthly_avg_icePresent"
  "timeMonthly_avg_iceAreaCell"
  "timeMonthly_avg_iceVolumeCell"
  "timeMonthly_avg_iceVolumeCell"
)
nvars="${#var2d_list[@]}"

exp_name="v3.LR.historical"
rundir="/lcrc/group/e3sm2/ac.wlin/E3SMv3"

# Work directory
WORK_DIR="/lcrc/group/e3sm/public_html/diagnostic_output/ac.szhang/v3LE_paper/data/${exp_name}"
mkdir -p "${WORK_DIR}/SE_SICE"

# Mapping file
MAP_FILE="/lcrc/group/acme/ac.szhang/acme_scratch/data/regrid_maps/map_IcoswISC30E3r4_to_1.0x1.0degree_conserve.nc"

start_year=1850
end_year=2050
period="${start_year}01-${end_year}12"
time_tag="$(printf "%04d" "$start_year")01-$(printf "%04d" "$end_year")12"
time_unt="days since $(printf "%04d" "$start_year")-01-01 00:00:0.0"

# key → name pairs (two columns)
read -r -d '' runs <<'EOF' || true
v3.LR.historical.en00 v3.LR.historical.0051
v3.LR.historical.en01 v3.LR.historical.0091
v3.LR.historical.en02 v3.LR.historical.0101
v3.LR.historical.en03 v3.LR.historical.0111
v3.LR.historical.en04 v3.LR.historical.0121
v3.LR.historical.en05 v3.LR.historical.0131
v3.LR.historical.en06 v3.LR.historical.0141
v3.LR.historical.en07 v3.LR.historical.0151
v3.LR.historical.en08 v3.LR.historical.0161
v3.LR.historical.en09 v3.LR.historical.0171
v3.LR.historical.en10 v3.LR.historical.0181
v3.LR.historical.en11 v3.LR.historical.0191
v3.LR.historical.en12 v3.LR.historical.0201
v3.LR.historical.en13 v3.LR.historical.0211
v3.LR.historical.en14 v3.LR.historical.0221
v3.LR.historical.en15 v3.LR.historical.0231
v3.LR.historical.en16 v3.LR.historical.0241
v3.LR.historical.en17 v3.LR.historical.0251
v3.LR.historical.en18 v3.LR.historical.0261
v3.LR.historical.en19 v3.LR.historical.0271
v3.LR.historical.en20 v3.LR.historical.0281
v3.LR.historical.en21 v3.LR.historical.0291
v3.LR.historical.en22 v3.LR.historical.0301
v3.LR.historical.en23 v3.LR.historical.0311
v3.LR.historical.en24 v3.LR.historical.0321
EOF

# Helper: join array with commas for ncrcat -v
join_by_comma() {
  local IFS=,
  echo "$*"
}

while read -r key name; do
  [[ -z "${key:-}" || -z "${name:-}" ]] && continue
  echo "$key => name=$name, period=$period"

  CASE_NAME="$key"
  exp_name="$name"
  # ensr = suffix token after last '.', e.g., 0051 or en00
  ensr="${name##*.}"

  RUN_FILE_DIR="${rundir}/${exp_name}/archive/ice/hist"
  if [[ -d "$RUN_FILE_DIR" ]]; then
    # Collect monthly files for the period
    eam_files=()
    for ((iy=start_year; iy<=end_year; iy++)); do
      yyyy="$(printf "%04d" "$iy")"
      # adjust the glob to match your filenames
      matches=( "$RUN_FILE_DIR"/*mpassi.hist.am.timeSeriesStatsMonthly*"${yyyy}"* )
      # Only append if glob matched something
      if compgen -G "$RUN_FILE_DIR/*mpassi.hist.am.timeSeriesStatsMonthly*${yyyy}*" >/dev/null; then
        eam_files+=( "${matches[@]}" )
      fi
    done

    if ((${#eam_files[@]} == 0)); then
      echo "WARNING: No input files found under ${RUN_FILE_DIR} for ${CASE_NAME}" >&2
      continue
    fi

    echo "Found ${#eam_files[@]} files for ${CASE_NAME}"

    SE_FILE="${WORK_DIR}/SE_SICE/${exp_name}.${ensr}.${time_tag}.nc"
    RG_FILE="${WORK_DIR}/SE_SICE/rgd.${exp_name}.${ensr}.${time_tag}.nc"

    # Build variable list for ncrcat (comma-separated, no spaces)
    fvars=( "timeMonthly_avg_iceAreaCell" )
    for vv in "${var2d_list[@]}"; do
      fvars+=( "$vv" )
    done
    fvars_csv="$(join_by_comma "${fvars[@]}")"

    if [[ ! -f "$SE_FILE" ]]; then
      ncrcat -d Time,0, -v "$fvars_csv" "${eam_files[@]}" "$SE_FILE"
    fi

    if [[ ! -f "$RG_FILE" ]]; then
      ncatted -t -a _FillValue,,o,d,-9.99999979021476795361e+33 \
                 -a _FillValue,,o,f,-9.99999979021476795361e+33 "$SE_FILE"
      ncremap -P mpasseaice -m "$MAP_FILE" -i "$SE_FILE" -o "$RG_FILE" -p mpi
    fi

    # Loop over variables
    for ((j=0; j<nvars; j++)); do
      var="${var2d_list[$j]}"
      vou="${var2d_name[$j]}"
      vfc="${var2d_frac[$j]}"
      vun="${var2d_unit[$j]}"

      FV_FILE="${WORK_DIR}/${exp_name}.${ensr}.${vou}.${time_tag}.nc"
      rm -f "$FV_FILE" "${FV_FILE}.tmp"

      ncks -v "$var" "$RG_FILE" "$FV_FILE"

      # Beware of names starting with '.' in ncrename syntax; quote them if needed.
      ncrename -v .Time,time -v ."${var}",${vou} -d .Time,time "$FV_FILE"

      ncap2 -s 'time=array(0,1,$time)*365.0/12+15' "$FV_FILE" "${FV_FILE}.tmp"
      ncatted -a units,time,m,c,"$time_unt" "${FV_FILE}.tmp"
      ncatted -a calendar,time,o,c,365_day  "${FV_FILE}.tmp"
      ncap2 -O -s 'defdim("bnds",2); time_bnds=make_bounds(time,$bnds,"time_bnds")' "${FV_FILE}.tmp" "$FV_FILE"
      rm -f "${FV_FILE}.tmp"

      ncap2 -s "${vou}=${vou}*${vfc}" -v "$FV_FILE" "${FV_FILE}.tmp"
      ncatted -a units,"$vou",a,c,"$vun" "${FV_FILE}.tmp"
      ncatted -a long_name,"$vou",a,c,"$var" "${FV_FILE}.tmp"

      mv -f "${FV_FILE}.tmp" "$FV_FILE"
    done
  else
    echo "WARNING: ${RUN_FILE_DIR} not found for ${CASE_NAME}" >&2
  fi
done <<< "$runs"

# finish the process and remove temp data
rm -rf "${WORK_DIR}/SE_SICE"
echo "done"
