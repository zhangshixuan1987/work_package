#!/bin/bash -el

echo "== Start of DART diagnostic =="
date
echo "============================================"

# System utilities
MOVE='/usr/bin/mv'
COPY='/usr/bin/cp --preserve=timestamps'
LINK='/usr/bin/ln -fs'
REMOVE='/usr/bin/rm'
LIST='/usr/bin/ls'

# Launch info
my_wkdir=${PWD}
cd ${my_wkdir}
source ./create_and_setup_case.sh

# Environment setup (assumes these are exported externally or in create_and_setup_case.sh)
# E3SM_ROOT, DART_ROOT, my_modeldir, my_ensnum, my_casename, etc.

DART_MODEL=${my_dart_eam}
DART_WORKDIR=${DART_ROOT}/models/${DART_MODEL}/work
ARCHIVE_DIR="${my_modeldir}/archive"

# Dates
ymds="2012-01-01"
ymde="2013-01-01"

# unique workspace to avoid overwriting previous tmp_run
TMPROOT=$(mktemp -d ./tmp_run.XXXXXX)
echo "Working in $TMPROOT"

for i in $(seq 1 "${my_ensnum}"); do
  ENSTR=$(printf "EN%02d" "${i}")
  CASE_NAME="${my_casename}.${ENSTR}"
  echo "=== Starting ensemble member ${ENSTR} ==="

  # per-ensemble subdir keeps files tidy and unique
  ENDIR="${TMPROOT}"
  mkdir -p "$ENDIR"

  # nothing to do if no scripts
  scripts=(post/*.sh)
  if ((${#scripts[@]}==0)); then
    echo "No scripts found in post/*.sh — skipping ${ENSTR}"
    continue
  fi

  for ff in "${scripts[@]}"; do
    dst="${ENDIR}/${ENSTR}_$(basename "$ff")"
    ${COPY} "$ff" "$dst"

    # GNU/BSD sed compatible inline helper
    if sed --version >/dev/null 2>&1; then
      SED_INPLACE=(-i)
    else
      SED_INPLACE=(-i '')
    fi

    # replace whole-line assignments safely (anchor at start, allow any RHS)
    sed "${SED_INPLACE[@]}" -E \
      -e "s|^ENSTR=.*$|ENSTR=${ENSTR}|g" \
      -e "s|^CASE_NAME=.*$|CASE_NAME=${CASE_NAME}|g" \
      -e "s|^ymds=.*$|ymds=${ymds}|g" \
      -e "s|^ymde=.*$|ymde=${ymde}|g" \
      -e "s|^ARCHIVE_DIR=.*$|ARCHIVE_DIR=${ARCHIVE_DIR}|g" \
      -e "s|^DART_WORKDIR=.*$|DART_WORKDIR=${DART_WORKDIR}|g" \
      -e "s|^DART_MODEL=.*$|DART_MODEL=${DART_MODEL}|g" \
      "$dst"
    sbatch ${dst}
  done
done

echo "==================================="
echo "All done. Outputs in: $TMPROOT"
echo "==================================="
