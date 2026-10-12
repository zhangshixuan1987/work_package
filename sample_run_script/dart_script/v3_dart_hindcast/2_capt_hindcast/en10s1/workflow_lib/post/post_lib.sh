#===============================================================================
# workflow_lib/post/post_lib.sh -- functions for the post-processing workers.
# Sourced by workflow_lib/post/*.sh after workflow_lib/common.sh.
# Outputs are written to a temporary file and renamed into place; existing
# outputs are skipped unless FORCE=1.
#
# Output layout (per member):
#   ENxx/archive/post/<atm|lnd>/180x360_aave/ts/<freq>/VAR.ENxx.YYYY.nc
#   ENxx/archive/post/<atm|lnd>/180x360_aave/clim/<history file name>
#===============================================================================

# post_months: print YYYY-MM for every month in [POST_START, POST_END).
post_months() {
  local d; d=$(date -d "${POST_START:0:7}-01" +%F)
  while [[ $d < $POST_END ]]; do
    echo "${d:0:7}"
    d=$(date -d "$d + 1 month" +%F)
  done
}

#-------------------------------------------------------------------------------
# regrid_member <m>: regrid one product for one member. The calling post/*.sh
# script sets P_COMP P_HIST P_FREQ P_MODE(ts|clim) P_MAP P_REMAP_OPTS
# P_EXTRA_VARS P_VARS. FORCE=1 overwrites existing outputs.
#-------------------------------------------------------------------------------
regrid_member() {
  local m=$1 product=${P_COMP}.${P_HIST}.${P_FREQ} force=${FORCE:-}
  local tag name archive in out work f ym var year yyyy
  tag=$(mtag "$m"); name=$(mcase "$m"); archive=$(marchive "$m")
  in=${archive}/${P_COMP}/hist
  out=${archive}/post/${P_COMP}/180x360_aave
  [[ -d $in ]] || die "${tag}: no history directory ${in}"
  [[ -f $P_MAP ]] || die "map file not found: ${P_MAP}"

  local -a files=()
  for ym in $(post_months); do
    if [[ $P_MODE == clim ]]; then
      f=${in}/${name}.${P_HIST}.${ym}.nc
      if [[ -f $f ]]; then files+=("$f"); fi
    else
      for f in "${in}/${name}.${P_HIST}.${ym}"*.nc; do
        if [[ -f $f ]]; then files+=("$f"); fi
      done
    fi
  done
  (( ${#files[@]} > 0 )) || die "${tag}: no ${P_HIST} files for ${POST_START}..${POST_END} in ${in}"
  log "${tag} ${product}: ${#files[@]} ${P_HIST} files, ${POST_START} to ${POST_END} (exclusive)"

  mkdir -p "$out"
  work=$(mktemp -d "${out}/tmp.${product}.XXXXXX")
  POST_WORK=$work
  trap 'rm -rf "$POST_WORK"' EXIT
  cd "$work"     # ncremap -P elm writes scratch files to the working directory

  local done_n=0 skip_n=0 miss_n=0 fail_n=0 dest
  if [[ $P_MODE == clim ]]; then
    dest=${out}/clim; mkdir -p "$dest"
    for f in "${files[@]}"; do
      if [[ -s ${dest}/$(basename "$f") && -z $force ]]; then skip_n=$((skip_n+1)); continue; fi
      if ncremap "${P_REMAP_OPTS[@]}" -m "$P_MAP" -i "$f" -o "${work}/out.nc" \
         && mv -f "${work}/out.nc" "${dest}/$(basename "$f")"; then
        done_n=$((done_n+1))
      else
        warn "${tag}: ncremap failed for $f"; fail_n=$((fail_n+1))
      fi
    done
  else
    dest=${out}/ts/${P_FREQ}; mkdir -p "$dest"
    local y0=${POST_START:0:4} y1; y1=$(date -d "$POST_END - 1 day" +%Y)
    for var in "${P_VARS[@]}"; do
      if ! _has_var "$var" "${files[@]}"; then
        warn "${tag}: variable ${var} not found in ${P_HIST} files"; miss_n=$((miss_n+1)); continue
      fi
      for year in $(seq "$y0" "$y1"); do
        yyyy=$(printf '%04d' "$year")
        local target=${dest}/${var}.${tag}.${yyyy}.nc
        if [[ -s $target && -z $force ]]; then skip_n=$((skip_n+1)); continue; fi
        if ncrcat -O -d time,"${yyyy}-01-01 00:00:00.0","${yyyy}-12-31 23:59:59.0" \
                  -v "${var}${P_EXTRA_VARS:+,${P_EXTRA_VARS}}" "${files[@]}" "${work}/cat.nc" \
           && ncremap "${P_REMAP_OPTS[@]}" -m "$P_MAP" -i "${work}/cat.nc" -o "${work}/out.nc" \
           && mv -f "${work}/out.nc" "$target"; then
          done_n=$((done_n+1))
        else
          warn "${tag}: failed ${var} ${yyyy}"; fail_n=$((fail_n+1))
        fi
        rm -f "${work}/cat.nc" "${work}/out.nc"
      done
    done
  fi

  log "${tag} ${product}: written=${done_n} skipped(existing)=${skip_n} missing_vars=${miss_n} failed=${fail_n}"
  (( fail_n == 0 )) || die "${tag} ${product}: ${fail_n} output(s) failed"
}

# _has_var <var> <files...>: true if any file contains the variable.
_has_var() {
  local var=$1 f; shift
  for f in "$@"; do
    if ncks -m -v "$var" "$f" >/dev/null 2>&1; then return 0; fi
  done
  return 1
}

#-------------------------------------------------------------------------------
# daily_member <m>: average 6-hourly regridded series to daily. Files already
# present in ts/daily (e.g. written from eam.h1 by eam_daily) are kept.
#-------------------------------------------------------------------------------
daily_member() {
  local m=$1 tag src dst f out n=0 fail=0
  tag=$(mtag "$m")
  src=$(marchive "$m")/post/atm/180x360_aave/ts/6hourly
  dst=$(marchive "$m")/post/atm/180x360_aave/ts/daily
  [[ -d $src ]] || die "${tag}: ${src} does not exist; run the eam_6hourly product first"
  mkdir -p "$dst"
  for f in "${src}"/*."${tag}".*.nc; do
    if [[ ! -f $f ]]; then continue; fi
    out=${dst}/$(basename "$f")
    if [[ -s $out ]]; then continue; fi
    if ncra -O --mro -d time,,,4,4 "$f" "${out}.tmp.nc" && mv -f "${out}.tmp.nc" "$out"; then
      n=$((n+1))
    else
      rm -f "${out}.tmp.nc"; warn "${tag}: ncra failed for $f"; fail=$((fail+1))
    fi
  done
  log "${tag}: ${n} daily file(s) written, ${fail} failed"
  (( fail == 0 )) || die "${tag}: daily conversion had failures"
}

#-------------------------------------------------------------------------------
# products: name -> worker | component | history stream | output subdir | inputs
#   inputs: day = one file per day, month = one per month, any = at least one,
#           none = derived from another product
#-------------------------------------------------------------------------------
post_products() { echo "eam_6hourly eam_daily eam_clim elm_daily elm_clim"; }

product_info() {
  case $1 in
    eam_6hourly) echo "post_process_eam_6hourly.sh atm eam.h2 ts/6hourly day" ;;
    eam_daily)   echo "post_process_eam_daily.sh atm eam.h1 ts/daily day" ;;
    eam_clim)    echo "post_process_eam_monthly.sh atm eam.h0 clim month" ;;
    elm_daily)   echo "post_process_elm_daily.sh lnd elm.h1 ts/daily any" ;;
    elm_clim)    echo "post_process_elm_monthly.sh lnd elm.h0 clim month" ;;
    eam_daily_from_6hourly) echo "- atm - ts/daily none" ;;
    *) return 1 ;;
  esac
}

post_record() { echo "postprocess_complete.$1.$(mtag "$2").${POST_START}_${POST_END}"; }

# check_post_inputs <product> <m>: print missing history inputs; return 1 if any.
check_post_inputs() {
  local p=$1 m=$2 name root comp hist need d
  name=$(mcase "$m"); root=$(marchive "$m")
  read -r _ comp hist _ need <<< "$(product_info "$p")"
  case $need in
    day)   for d in $(dates_in_range "$POST_START" "$POST_END"); do
             if ! ls "${root}/${comp}/hist/${name}.${hist}.${d}-"*.nc >/dev/null 2>&1; then
               echo "$(mtag "$m"): missing ${hist} for ${d}"; return 1
             fi
           done ;;
    month) for d in $(months_in_range "$POST_START" "$POST_END"); do
             if [[ ! -s ${root}/${comp}/hist/${name}.${hist}.${d}.nc ]]; then
               echo "$(mtag "$m"): missing ${name}.${hist}.${d}.nc"; return 1
             fi
           done ;;
    any)   if ! ls "${root}/${comp}/hist/${name}.${hist}."*.nc >/dev/null 2>&1; then
             echo "$(mtag "$m"): no ${hist} files"; return 1
           fi ;;
  esac
  return 0
}

# validate_post_outputs <product> <m>: outputs exist, are non-empty and readable.
validate_post_outputs() {
  local p=$1 m=$2 t comp hist sub dest f d n=0
  t=$(mtag "$m")
  read -r _ comp hist sub _ <<< "$(product_info "$p")"
  dest=$(marchive "$m")/post/${comp}/180x360_aave/${sub}
  if [[ $sub == clim ]]; then
    for d in $(months_in_range "$POST_START" "$POST_END"); do
      f=${dest}/$(mcase "$m").${hist}.${d}.nc
      if [[ ! -s $f ]] || ! ncdump -h "$f" >/dev/null 2>&1; then return 1; fi
      n=$((n+1))
    done
  else
    for f in "${dest}"/*."${t}".*.nc; do
      [[ -f $f ]] || continue
      if [[ ! -s $f ]] || ! ncdump -h "$f" >/dev/null 2>&1; then return 1; fi
      n=$((n+1))
    done
  fi
  (( n > 0 ))
}

# post_job <product> <m>: run one product for one member and write its
# completion record (runtmp/status/postprocess_complete.<product>.ENxx.<range>).
# The regrid settings (P_*) must already be set by the calling job script.
post_job() {
  local p=$1 m=$2 t rec prog
  t=$(mtag "$m"); rec=$(post_record "$p" "$m"); prog=${rec/complete/in_progress}
  exec 7>"${LOCK_DIR}/postprocess.${p}.${t}.lock"
  flock -n 7 || die "${t} ${p}: another job is already processing it"
  if [[ -z ${FORCE:-} ]] && record_valid "$rec" && validate_post_outputs "$p" "$m"; then
    log "${t} ${p}: verified earlier, skipping (OVERWRITE_EXISTING=TRUE to redo)"; return 0
  fi
  record_rm "$rec"
  record_write "$prog" "product=${p}" "member=${t}" "start=${POST_START}" "end=${POST_END}"
  if [[ $p == eam_daily_from_6hourly ]]; then
    record_valid "$(post_record eam_6hourly "$m")" || die "${t}: eam_6hourly has not completed for ${POST_START}_${POST_END}"
    daily_member "$m"
  else
    regrid_member "$m"
  fi
  validate_post_outputs "$p" "$m" || die "${t} ${p}: outputs missing or unreadable"
  record_write "$rec" "product=${p}" "member=${t}" "start=${POST_START}" "end=${POST_END}"
  record_rm "$prog"
  log "${t} ${p}: complete"
}
