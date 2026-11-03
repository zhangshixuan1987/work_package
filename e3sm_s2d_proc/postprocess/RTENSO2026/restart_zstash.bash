#!/bin/bash
# Restart zstash1-3.bash after the HPSS outage of 2026-11-01:
#   1. wait until HPSS answers
#   2. put the tar each 2015020100 create left in its cache (hpss_put hung)
#   3. relaunch zstash1-3.bash (they run "zstash update" on the 2015020100 leads)

source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

run1_dir="/pscratch/sd/z/zhan391/e3sm_project/E3SMv3_S2D"
hpss_pth="/home/projects/e3sm/www/RTENSO2026"
post_dir="$(cd "$(dirname "$0")" && pwd)"
pending=(clim:000011.tar obs:000011.tar obs_ice_nsidc:000012.tar)

log_msg() { echo "[$(date '+%F %T')] $*"; }

until timeout 120 hsi -q pwd > /dev/null 2>&1; do
    log_msg "HPSS not responding, retry in 10 min"
    sleep 600
done
log_msg "HPSS is back"

for item in "${pending[@]}"; do
    ic="${item%%:*}"
    tar="${item##*:}"
    case_name="WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_RTENSO2026_${ic}_2015020100"
    cache="${run1_dir}/${case_name}/zstash"
    [[ -f "${cache}/${tar}" ]] || { log_msg "${case_name}: no ${tar} in cache, skip"; continue; }
    local_size=$(stat -c %s "${cache}/${tar}")
    log_msg "${case_name}: put ${tar} (${local_size} bytes)"
    ( cd "${cache}" && hsi -q "cd ${hpss_pth}/${case_name}; put ${tar}" ) || { log_msg "ERROR: put failed for ${case_name}/${tar}"; exit 1; }
    hpss_size=$(hsi -q "ls -l ${hpss_pth}/${case_name}/${tar}" 2>&1 | awk -v t="${tar}" '$NF ~ t {print $5}')
    if [[ "${hpss_size}" != "${local_size}" ]]; then
        log_msg "ERROR: size mismatch for ${case_name}/${tar}: local ${local_size}, HPSS ${hpss_size}"
        exit 1
    fi
    rm -f "${cache}/${tar}"
    log_msg "${case_name}: ${tar} on HPSS, removed from cache"
done

cd "${post_dir}"
for n in 1 2 3; do
    nohup ./zstash${n}.bash > sh${n}.log 2>&1 < /dev/null &
    log_msg "launched zstash${n}.bash (pid $!)"
    sleep 5
done
