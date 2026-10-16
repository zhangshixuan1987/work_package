#!/bin/bash -el

myenv="${1:?missing myenv}"
scomp="${2:?missing scomp}"
SRC_ROOT="${3:?missing SRC_ROOT}"
SRC_FILE="${4:?missing SRC_FILE}"
TGT_ROOT="${5:?missing TGT_ROOT}"
REST_FILE="${6:?missing REST_FILE}"
START_DATE="${7:?missing START_DATE}"
START_HOUR="${8:?missing START_HOUR}"
RPOINTER_FILE="${9:?missing RPOINTER_FILE}"

# isolate environment here
eval "${my_env}"

mkdir -p "${TGT_ROOT}"

if [[ "${scomp}" == "ocn" || "${scomp}" == "ice" ]]; then
    tmpfile="${TGT_ROOT}/${REST_FILE}.tmp"

    cp -p "${SRC_ROOT}/${SRC_FILE}" "${tmpfile}"
    ncks -O --hdr_pad=10000 "${tmpfile}" "${TGT_ROOT}/${REST_FILE}"
    ncrename -v xtime,xtime.orig "${TGT_ROOT}/${REST_FILE}"
    rm -f "${tmpfile}"

    printf "%s_%02d:00:00\n" "${START_DATE}" "${START_HOUR}" > "${TGT_ROOT}/${RPOINTER_FILE}"
else
    cp -p "${SRC_ROOT}/${SRC_FILE}" "${TGT_ROOT}/${REST_FILE}"
    printf "./%s\n" "${REST_FILE}" > "${TGT_ROOT}/${RPOINTER_FILE}"
fi
