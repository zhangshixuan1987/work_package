#!/bin/bash -e
#
# Compare a candidate initial-condition directory against its reference.
#
# Usage: bfb_compare.sh <reference_dir> <candidate_dir> [expected_count]
#
# Each reference file is checked byte-for-byte with cmp first; a NetCDF file
# that differs is then checked with cprnc, which ignores metadata such as the
# history attribute. The report is written to <candidate_dir>/BFB_REPORT.txt.
# Exits 0 only if every file is identical (strict or cprnc), none is missing or
# extra, and the reference holds expected_count files (when given).

set -o pipefail

cprnc=${CPRNC:-/global/cfs/cdirs/e3sm/tools/cprnc/cprnc}

die() {
  echo "BFB ERROR: $*" >&2
  exit 1
}

(( $# == 2 || $# == 3 )) || die "usage: ${0##*/} <reference_dir> <candidate_dir> [expected_count]"
reference_dir=$1
candidate_dir=$2
expected_count=${3:-}

[[ -d ${reference_dir} ]] || die "reference directory not found: ${reference_dir}"
[[ -d ${candidate_dir} ]] || die "candidate directory not found: ${candidate_dir}"
[[ -x ${cprnc} ]] || die "cprnc not executable: ${cprnc}"

reference_files=("${reference_dir}"/*)
if [[ -n ${expected_count} ]] && (( ${#reference_files[@]} != expected_count )); then
  die "expected ${expected_count} reference files, found ${#reference_files[@]}"
fi

report="${candidate_dir}/BFB_REPORT.txt"
strict=0
cprnc_identical=0
failures=0
missing=0
extras=0

{
  echo "BFB reference: ${reference_dir}"
  echo "BFB candidate: ${candidate_dir}"
  echo
} > "${report}"

for reference_file in "${reference_files[@]}"; do
  name=${reference_file##*/}
  candidate_file="${candidate_dir}/${name}"
  if [[ ! -f ${candidate_file} ]]; then
    echo "MISSING ${name}" >> "${report}"
    missing=$((missing + 1))
    continue
  fi
  if cmp -s "${reference_file}" "${candidate_file}"; then
    echo "STRICT_MATCH ${name}" >> "${report}"
    strict=$((strict + 1))
    continue
  fi
  if [[ ${name} != *.nc ]]; then
    echo "MISMATCH ${name}" >> "${report}"
    failures=$((failures + 1))
    continue
  fi

  detail="${candidate_file}.cprnc.txt"
  "${cprnc}" -m "${reference_file}" "${candidate_file}" > "${detail}" 2>&1 || true
  if grep -q 'diff_test: the two files seem to be IDENTICAL' "${detail}"; then
    echo "CPRNC_MATCH ${name}" >> "${report}"
    cprnc_identical=$((cprnc_identical + 1))
    rm -f "${detail}"
  else
    echo "MISMATCH ${name}: ${detail}" >> "${report}"
    failures=$((failures + 1))
  fi
done

for candidate_file in "${candidate_dir}"/*; do
  name=${candidate_file##*/}
  [[ ${name} == *_REPORT.txt || ${name} == *.cprnc.txt ]] && continue
  if [[ ! -e ${reference_dir}/${name} ]]; then
    echo "EXTRA ${name}" >> "${report}"
    extras=$((extras + 1))
  fi
done

{
  echo
  echo "strict_matches=${strict}"
  echo "cprnc_matches=${cprnc_identical}"
  echo "mismatches=${failures}"
  echo "missing=${missing}"
  echo "extras=${extras}"
} >> "${report}"

cat "${report}"
if (( failures > 0 || missing > 0 || extras > 0 )); then
  die "BFB comparison failed; see ${report}"
fi
if (( cprnc_identical > 0 )); then
  echo "BFB PASS (cprnc): ${strict} files byte-identical, ${cprnc_identical} identical in data."
else
  echo "BFB PASS (strict): all ${strict} files are byte-for-byte identical."
fi
