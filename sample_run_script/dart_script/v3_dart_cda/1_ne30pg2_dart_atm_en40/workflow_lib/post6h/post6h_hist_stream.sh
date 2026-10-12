#!/bin/bash
# Build 6-hourly, daily and monthly products for one member, one history
# stream and one calendar month (called by 9_run_post_6hourly.sh).
#
#   POST6H_MEMBER  ENxx            POST6H_COMP    atm | lnd | rof
#   POST6H_STREAM  e.g. eam.h2     POST6H_MONTH   YYYY-MM
#   POST6H_RULE    mean6h (default): use the record whose bounds span one cycle
#                  instant: use the single record valid at the cycle start, for
#                  instantaneous streams such as the ELM-DART input (daily and
#                  monthly values are then means of the per-cycle snapshots)
#
# Output: ${my_modeldir}/ENxx/archive/post6h/<comp>/180x360_aave/
#           6hourly/VAR.ENxx.YYYY-MM.nc   one record per DA cycle
#           daily/VAR.ENxx.YYYY-MM.nc     mean (or min/max) of each day's cycles
#           monthly/VAR.ENxx.YYYY-MM.nc   mean of all cycles (min/max: of the daily extremes)
# Time-mean fields are stamped at the middle of their 6-h window; instantaneous
# and min/max fields keep the window end. Days are grouped by cycle, so the
# 18Z-24Z cycle always belongs to its start day.
set -o pipefail
[[ "${POSTPROCESS_DRIVER_ACTIVE:-FALSE}" == "TRUE" ]] || { echo "ERROR: run 9_run_post_6hourly.sh" >&2; exit 1; }
source "${my_workflow_lib}/post6h/post6h_common.sh"

member=${POST6H_MEMBER:?}; comp=${POST6H_COMP:?}; stream=${POST6H_STREAM:?}; ym=${POST6H_MONTH:?}
rule=${POST6H_RULE:-mean6h}
[[ "${rule}" =~ ^(mean6h|instant)$ ]] || post6h_fail "unknown record rule: ${rule}"
case_name="${my_casename}.${member}"
src="${my_modeldir}/${member}/archive/${comp}/hist"
dst="${my_modeldir}/${member}/archive/post6h/${comp}/180x360_aave"
case "${comp}" in
  atm) map="${my_eam_hist_mapping_file}"; rgr=() ;;
  lnd) map="${my_elm_hist_mapping_file}"; rgr=(-P elm) ;;
  rof) map="${my_rof_hist_mapping_file:-${my_elm_hist_mapping_file}}"; rgr=() ;;
  *) post6h_fail "unknown component: ${comp}" ;;
esac
[[ -r "${map}" ]] || post6h_fail "mapping file not readable: ${map}"

post6h_activate_env
mkdir -p "${dst}" || post6h_fail "cannot create ${dst}"
work=$(mktemp -d "${dst}/.work.${stream}.${member}.${ym}.XXXXXX") || post6h_fail "cannot create work directory"
trap 'rm -rf -- "${work}"' EXIT
export TMPDIR="${work}"

per_day=$(post6h_cycles_per_day)
nd=$(post6h_days_in_month "${ym}")
mapfile -t stamps < <(post6h_month_stamps "${ym}")
(( ${#stamps[@]} == nd * per_day )) || post6h_fail "cycle list for ${ym} has ${#stamps[@]} entries, expected $(( nd * per_day ))"
echo "== ${case_name} ${stream} ${ym}: ${#stamps[@]} cycles"

# 1. One record per cycle, selected by the stream's record rule.
i=0
for st in "${stamps[@]}"; do
  f="${src}/${case_name}.${stream}.${st}.nc"
  [[ -s "${f}" ]] || post6h_fail "missing history file: ${f}"
  if [[ "${rule}" == "mean6h" ]]; then
    idx=$(post6h_cycle_record "${f}") ||
      post6h_fail "no record spanning exactly ${my_e3sm_cycle_hours} h in ${f} (is this stream written every ${my_e3sm_cycle_hours} h as an average?)"
  else
    idx=$(post6h_instant_record "${f}" "${st}") ||
      post6h_fail "no single record valid at ${st} in ${f} (is this an instantaneous ${my_e3sm_cycle_hours}-hourly stream?)"
  fi
  ncks -O -d "time,${idx},${idx}" "${f}" "${work}/rec_$(printf '%04d' "${i}").nc" || post6h_fail "cannot extract record from ${f}"
  i=$(( i + 1 ))
done
ncrcat -O "${work}"/rec_*.nc "${work}/all.nc" || post6h_fail "ncrcat failed"
rm -f -- "${work}"/rec_*.nc
[[ "$(post6h_ntime "${work}/all.nc")" == "${#stamps[@]}" ]] || post6h_fail "concatenated file has wrong record count"
hdr=$(ncdump -h "${work}/all.nc") || post6h_fail "cannot read ${work}/all.nc"
grep -qE "time:units = \"days since ${ym}-01( 00:00:00)?\"" <<<"${hdr}" ||
  post6h_fail "unexpected time units: $(grep 'time:units' <<<"${hdr}")"
bnds=$(post6h_bounds_name "${work}/all.nc") || post6h_fail "no time bounds variable"

# 2. History fields are the variables with a time cell_method.
mapfile -t fields < <(sed -nE 's/^\t\t([A-Za-z0-9_]+):cell_methods = "time: ([a-z]+)".*/\1 \2/p' <<<"${hdr}")
(( ${#fields[@]} > 0 )) || post6h_fail "no history fields found in ${stream}"
extra=""
if [[ "${comp}" == "lnd" ]]; then
  ncks -m -v landfrac "${work}/all.nc" >/dev/null 2>&1 || post6h_fail "landfrac missing; required for ncremap -P elm"
  extra=",landfrac"
fi

# 3. Per field: regrid the 6-hourly series, then derive daily and monthly.
names=()
for entry in "${fields[@]}"; do
  var=${entry% *}; method=${entry#* }
  [[ "${var}" == "landfrac" ]] && continue
  case "${method}" in minimum) op=min ;; maximum) op=max ;; *) op=avg ;; esac
  v="${work}/${var}"
  ncks -O -v "${var}${extra},${bnds}" "${work}/all.nc" "${v}.nat.nc" || post6h_fail "extract ${var}"
  if [[ "${method}" == "mean" && "${rule}" == "mean6h" ]]; then
    ncap2 -O -s "time=0.5*(${bnds}(:,0)+${bnds}(:,1))" "${v}.nat.nc" "${v}.nat.nc" || post6h_fail "midpoint ${var}"
  fi
  ncremap "${rgr[@]}" -m "${map}" -i "${v}.nat.nc" -o "${v}.6h.nc" >/dev/null || post6h_fail "ncremap ${var}"
  ncra -O --mro -y "${op}" -d "time,0,,${per_day},${per_day}" "${v}.6h.nc" "${v}.day.nc" || post6h_fail "daily ${var}"
  # Label day k (0-based) at its noon from the record index; the reduced time
  # value itself is unreliable for min/max (e.g. max picks the next day's 00Z).
  ncap2 -O -s "time=array(0.5,1.0,\$time); ${bnds}(:,0)=time-0.5; ${bnds}(:,1)=time+0.5" "${v}.day.nc" "${v}.day.nc" || post6h_fail "daily time ${var}"
  # Monthly mean of all cycles; for min/max fields, the mean of the daily extremes.
  if [[ "${op}" == "avg" ]]; then mon_src="${v}.6h.nc"; else mon_src="${v}.day.nc"; fi
  ncra -O -y avg "${mon_src}" "${v}.mon.nc" || post6h_fail "monthly ${var}"
  ncap2 -O -s "time=0.5*${nd}; ${bnds}(:,0)=0; ${bnds}(:,1)=${nd}" "${v}.mon.nc" "${v}.mon.nc" || post6h_fail "monthly time ${var}"
  out="${var}.${member}.${ym}.nc"
  post6h_publish "${v}.6h.nc"  "${dst}/6hourly/${out}" "${#stamps[@]}" || post6h_fail "publish 6hourly ${var}"
  post6h_publish "${v}.day.nc" "${dst}/daily/${out}"   "${nd}"         || post6h_fail "publish daily ${var}"
  post6h_publish "${v}.mon.nc" "${dst}/monthly/${out}" 1              || post6h_fail "publish monthly ${var}"
  rm -f -- "${v}".*.nc
  names+=("${var}:${method}")
done

post6h_record "${my_status_dir}/post6h_complete.${stream}.${member}.${ym}" \
  "case=${my_casename}" "member=${member}" "stream=${stream}" "month=${ym}" \
  "cycles=${#stamps[@]}" "record_rule=${rule}" \
  "n_fields=${#names[@]}" "fields=$(IFS=,; echo "${names[*]}")" || post6h_fail "cannot write completion record"
echo "== done ${case_name} ${stream} ${ym}: ${#names[@]} fields"
