#!/bin/bash
# Internal compression worker for one member. Run 6_run_compress_data.sh instead.
#
# Deflates (netCDF-4 classic, level 5) the member's archived history files
# (COMPRESS_HISTORY) and/or restart files (COMPRESS_RESTARTS). Each file is
# compressed to a temporary copy, checked with ncdump, and renamed over the
# original; files that are already netCDF-4 are skipped, and leftovers of an
# interrupted run are removed first.
set -Eeo pipefail
if [[ "${COMPRESS_DRIVER_ACTIVE:-FALSE}" != "TRUE" ]]; then
  echo "ERROR: internal compression worker; run 6_run_compress_data.sh" >&2
  exit 1
fi
source "${WORKFLOW_ROOT:?}/workflow_lib/common.sh"
wf_init "${WORKFLOW_ROOT}"

m=${COMPRESS_MEMBER:?}
archive=$(marchive "$m")
[[ -d $archive ]] || die "$(mtag "$m"): no archive ${archive}"
dirs=()
if [[ ${COMPRESS_HISTORY:?} == TRUE ]]; then
  for comp in atm cpl lnd rof ice ocn; do
    if [[ -d ${archive}/${comp}/hist ]]; then dirs+=("${archive}/${comp}/hist"); fi
  done
fi
if [[ ${COMPRESS_RESTARTS:?} == TRUE && -d ${archive}/rest ]]; then dirs+=("${archive}/rest"); fi
(( ${#dirs[@]} > 0 )) || { warn "$(mtag "$m"): nothing to compress"; exit 0; }

find "${dirs[@]}" -name '*.hc_tmp' -type f -delete

log "$(mtag "$m"): compressing ${dirs[*]#${archive}/} with ${COMPRESS_MAX_PARALLEL:?} processes"
find "${dirs[@]}" -type f -name '*.nc' -print0 \
  | xargs -0 -r -n 1 -P "$COMPRESS_MAX_PARALLEL" bash -c '
      f=$1
      case $(ncdump -k "$f" 2>/dev/null) in netCDF-4*) exit 0 ;; esac
      nccopy -7 -d 5 "$f" "$f.hc_tmp" && ncdump -h "$f.hc_tmp" >/dev/null \
        && touch -r "$f" "$f.hc_tmp" && mv -f "$f.hc_tmp" "$f" \
        || { rm -f "$f.hc_tmp"; echo "compress failed: $f" >&2; exit 1; }' _ \
  || die "$(mtag "$m"): some files failed to compress (see messages above)"

# every file must now be netCDF-4
left=$(find "${dirs[@]}" -type f -name '*.nc' -print0 | xargs -0 -r -n 50 ncdump -k 2>/dev/null | grep -vc '^netCDF-4' || true)
[[ ${left:-0} == 0 ]] || die "$(mtag "$m"): ${left} file(s) still uncompressed"
log "$(mtag "$m"): compression finished"
