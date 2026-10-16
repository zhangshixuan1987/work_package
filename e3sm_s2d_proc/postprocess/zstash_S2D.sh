#!/bin/bash

set -eo pipefail

usage() {
    echo "Usage:"
    echo "  $0 dry-run {create|update|check} CASE_NAME"
    echo "  $0 {create|update|check} CASE_NAME"
    exit 1
}

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    usage
fi

mode="$1"

if [ "$mode" = "dry-run" ]; then
    if [ "$#" -ne 3 ]; then
        usage
    fi
    action="$2"
    casenm="$3"
else
    if [ "$#" -ne 2 ]; then
        usage
    fi
    action="$mode"
    casenm="$2"
fi

case "$action" in
    create|update|check)
        ;;
    *)
        usage
        ;;
esac

source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
set -u

hpss_pth="/home/projects/e3sm/www/S2S2D"
date_str=$(date +%Y%m%d)
logfn="zstash_${action}_${date_str}.log"

case_dir="${casenm}"

zstash_create_cmd="zstash create -v --hpss=${hpss_pth}/${casenm} --maxsize 128 ."
zstash_update_cmd="zstash update -v --hpss=${hpss_pth}/${casenm} ."
zstash_check_cmd="zstash check --hpss=${hpss_pth}/${casenm}"
hsi_cmd="hsi \"cd ${hpss_pth}/${casenm}; cput -R ${logfn}\""

case "$mode" in
    dry-run)
        echo "Would run:"
        echo "cd ${case_dir}"
        echo "mkdir -p zstash"

        if [ "$action" = "create" ]; then
            echo "${zstash_create_cmd} 2>&1 | tee zstash/${logfn}"
        elif [ "$action" = "update" ]; then
            echo "${zstash_update_cmd} 2>&1 | tee zstash/${logfn}"
        else
            echo "${zstash_check_cmd} 2>&1 | tee zstash/${logfn}"
        fi

        echo "cd zstash"
        echo "${hsi_cmd}"
        ;;

    create)
        pushd "${case_dir}" > /dev/null
        mkdir -p zstash
        zstash create -v --hpss="${hpss_pth}/${casenm}" --maxsize 128 . 2>&1 | tee "zstash/${logfn}"
        pushd zstash > /dev/null
        hsi "cd ${hpss_pth}/${casenm}; cput -R ${logfn}"
        popd > /dev/null
        popd > /dev/null
        ;;

    update)
        pushd "${case_dir}" > /dev/null
        mkdir -p zstash
        zstash update -v --hpss="${hpss_pth}/${casenm}"  2>&1 | tee "zstash/${logfn}"
        pushd zstash > /dev/null
        hsi "cd ${hpss_pth}/${casenm}; cput -R ${logfn}"
        popd > /dev/null
        popd > /dev/null
        ;;

    check)
        pushd "${case_dir}" > /dev/null
        mkdir -p zstash
        zstash check --hpss="${hpss_pth}/${casenm}" 2>&1 | tee "zstash/${logfn}"
        pushd zstash > /dev/null
        hsi "cd ${hpss_pth}/${casenm}; cput -R ${logfn}"
        popd > /dev/null
        popd > /dev/null
        ;;

    *)
        usage
        ;;
esac
