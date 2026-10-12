#!/bin/bash
# Build state-verification products from the DART ensemble mean and spread for
# one calendar month (called by 9_run_post_6hourly.sh).
#
#   POST6H_MONTH  YYYY-MM
#   POST6H_KIND   forecast (prior, dart.e.eam_forecast_*) | analysis (posterior, dart.i.eam_output_*)
#
# Each cycle file is interpolated from model levels to pressure levels (below
# ground masked as missing) and regridded from the ne30np4 dynamics grid to 1 deg.
# Output: ${my_dart_root}/post6h/180x360_plev/
#           6hourly/VAR.<kind>_<mean|sd>.YYYY-MM.nc   one record per cycle valid time
#           daily/  mean of each day's 00/06/12/18Z values
#           monthly/ mean of all cycles in the month
# The spread is interpolated with the ensemble-mean surface pressure; linear
# interpolation of a standard deviation is an approximation.
set -o pipefail
[[ "${POSTPROCESS_DRIVER_ACTIVE:-FALSE}" == "TRUE" ]] || { echo "ERROR: run 9_run_post_6hourly.sh" >&2; exit 1; }
source "${my_workflow_lib}/post6h/post6h_common.sh"

ym=${POST6H_MONTH:?}; kind=${POST6H_KIND:?}
case "${kind}" in
  forecast) prefix="dart.e.eam_forecast" ;;
  analysis) prefix="dart.i.eam_output" ;;
  *) post6h_fail "unknown DART kind: ${kind}" ;;
esac
vars=(${POST6H_DART_VARS:-T U V Q CLDLIQ CLDICE PS})
map="${my_eam_dart_mapping_file:-/compyfs/zhan391/acme_init/map_file/map_ne30np4_to_cmip6_180x360_aave.c231227.nc}"
plevs="${POST6H_PLEVS_PA:-100000,92500,85000,70000,60000,50000,40000,30000,25000,20000,15000,10000,7000,5000,3000,2000,1000}"
[[ -r "${map}" ]] || post6h_fail "mapping file not readable: ${map}"
dst="${my_dart_root}/post6h/180x360_plev"

post6h_activate_env
mkdir -p "${dst}" || post6h_fail "cannot create ${dst}"
work=$(mktemp -d "${dst}/.work.${kind}.${ym}.XXXXXX") || post6h_fail "cannot create work directory"
trap 'rm -rf -- "${work}"' EXIT
export TMPDIR="${work}"

nplev=$(tr ',' '\n' <<<"${plevs}" | wc -l)
ncap2 -O -v -s "defdim(\"plev\",${nplev}); plev[plev]={${plevs}}; plev@units=\"Pa\"; plev@long_name=\"pressure\"; plev@positive=\"down\";" \
  "${work}/plev.nc" || post6h_fail "cannot build pressure-level file"

per_day=$(post6h_cycles_per_day)
nd=$(post6h_days_in_month "${ym}")
mapfile -t stamps < <(post6h_month_stamps "${ym}")
(( ${#stamps[@]} == nd * per_day )) || post6h_fail "cycle list for ${ym} has ${#stamps[@]} entries"
month0=$(date -u -d "${ym}-01 00:00:00" +%s)
sel="$(IFS=,; echo "${vars[*]}"),hyam,hybm,hyai,hybi,P0,lev,ilev"
echo "== ${my_casename} DART ${kind} ${ym}: ${#stamps[@]} cycles, vars ${vars[*]}"

# Add a time coordinate (the files have a time dimension but no variable).
add_time() { ncap2 -O -s "time[time]=$2; time@units=\"days since ${ym}-01 00:00:00\"; time@calendar=\"gregorian\"" "$1" "$1"; }
# --rnr_thr=0.5: a 1-deg cell partly below ground is renormalized by its valid
# fraction (missing if under half is valid); without it partial cells are biased low.
remap() { ncremap -T "${work}" --rnr_thr=0.5 --vrt_out="${work}/plev.nc" --vrt_xtr=mss_val -m "${map}" -i "$1" -o "$2" >/dev/null; }
# The first cycle's files (written by the perturbation step) declare a missing
# value although the DART state has none; ncrcat takes attributes from the first
# file, so that declaration would hide the below-ground mask written by ncremap.
# Require a complete state, then drop the inherited missing-value attributes.
strip_fill() {
  local expr n
  expr=$(printf '%s.number_miss()+' "${vars[@]}"); expr=${expr%+}
  ncap2 -O -v -s "nmiss=${expr};" "$1" "$1.nmiss" || return 1
  n=$(ncdump -v nmiss "$1.nmiss" | sed -n 's/^ nmiss = \([0-9]*\) ;/\1/p'); rm -f -- "$1.nmiss"
  [[ "${n}" == "0" ]] || { echo "ERROR: $1 has ${n:-unknown} missing values in the DART state" >&2; return 1; }
  ncatted -O -a _FillValue,,d,, -a missing_value,,d,, "$1"
}

i=0
for st in "${stamps[@]}"; do
  dir="${my_eam_dart_run_dir}/${st}"
  fm="${dir}/${my_casename}.${prefix}_mean.${st}.nc"
  fs="${dir}/${my_casename}.${prefix}_sd.${st}.nc"
  [[ -s "${fm}" && -s "${fs}" ]] || post6h_fail "missing DART ${kind} mean/sd for ${st} in ${dir}"
  t=$(awk -v a="$(date -u -d "${st%-*} 00:00:00" +%s)" -v s="$((10#${st##*-}))" -v m="${month0}" 'BEGIN{printf "%.6f", (a+s-m)/86400}')
  n=$(printf '%04d' "${i}")
  # mean
  ncks -O -C -v "${sel}" "${fm}" "${work}/m.nc" && ncrename -O -d ncol_d,ncol "${work}/m.nc" >/dev/null && strip_fill "${work}/m.nc" &&
    add_time "${work}/m.nc" "${t}" && remap "${work}/m.nc" "${work}/rm_${n}.nc" || post6h_fail "mean ${st}"
  # spread: interpolate with the ensemble-mean PS, keep the PS spread as PS
  ncks -O -C -v "${sel}" "${fs}" "${work}/s.nc" && strip_fill "${work}/s.nc" && ncrename -O -d ncol_d,ncol -v PS,PS_SD "${work}/s.nc" >/dev/null &&
    ncks -A -C -v PS "${work}/m.nc" "${work}/s.nc" && add_time "${work}/s.nc" "${t}" &&
    remap "${work}/s.nc" "${work}/s_${n}.nc" && ncks -O -x -v PS "${work}/s_${n}.nc" "${work}/rs_${n}.nc" &&
    ncrename -O -v PS_SD,PS "${work}/rs_${n}.nc" >/dev/null || post6h_fail "sd ${st}"
  rm -f -- "${work}/m.nc" "${work}/s.nc" "${work}/s_${n}.nc"
  i=$(( i + 1 ))
done

for stat in mean sd; do
  tag="r${stat:0:1}"
  ncrcat -O "${work}/${tag}"_*.nc "${work}/${stat}.nc" || post6h_fail "ncrcat ${stat}"
  rm -f -- "${work}/${tag}"_*.nc
  [[ "$(post6h_ntime "${work}/${stat}.nc")" == "${#stamps[@]}" ]] || post6h_fail "${stat}: wrong record count"
  for var in "${vars[@]}"; do
    v="${work}/${var}"
    ncks -O -v "${var}" "${work}/${stat}.nc" "${v}.6h.nc" || post6h_fail "extract ${var}"
    ncra -O --mro -d "time,0,,${per_day},${per_day}" "${v}.6h.nc" "${v}.day.nc" &&
      ncap2 -O -s "time=array(0.5,1.0,\$time)" "${v}.day.nc" "${v}.day.nc" || post6h_fail "daily ${var}"
    # Monthly mean of all cycles (not of daily means): below-ground points have a
    # different number of valid samples on different days.
    ncra -O "${v}.6h.nc" "${v}.mon.nc" && ncap2 -O -s "time=0.5*${nd}" "${v}.mon.nc" "${v}.mon.nc" || post6h_fail "monthly ${var}"
    out="${var}.${kind}_${stat}.${ym}.nc"
    post6h_publish "${v}.6h.nc"  "${dst}/6hourly/${out}" "${#stamps[@]}" || post6h_fail "publish 6hourly ${var}"
    post6h_publish "${v}.day.nc" "${dst}/daily/${out}"   "${nd}"         || post6h_fail "publish daily ${var}"
    post6h_publish "${v}.mon.nc" "${dst}/monthly/${out}" 1              || post6h_fail "publish monthly ${var}"
    rm -f -- "${v}".*.nc
  done
  rm -f -- "${work}/${stat}.nc"
done

post6h_record "${my_status_dir}/post6h_complete.dart_${kind}.${ym}" \
  "case=${my_casename}" "kind=${kind}" "month=${ym}" "cycles=${#stamps[@]}" \
  "vars=$(IFS=,; echo "${vars[*]}")" "plev_pa=${plevs}" "below_ground=missing" || post6h_fail "cannot write completion record"
echo "== done DART ${kind} ${ym}"
