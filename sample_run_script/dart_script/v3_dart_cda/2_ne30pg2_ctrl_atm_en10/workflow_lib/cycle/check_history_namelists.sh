#!/bin/bash
# Check that every EAM, ELM and MOSART history tape of a DA-cycled case gives
# one complete statistic per cycle, before any cycle runs:
#   usage: check_history_namelists.sh <case_scripts dir> <cycle hours>
#
# Each DA cycle is a separate run of <cycle hours>, so for every tape in use
# (tape 1 and every tape with a fincl list):
#   nhtfrq = -<cycle hours>  A longer interval is never reached inside a cycle:
#                            the tape then holds only the nstep-0 record (its
#                            means, maxima and minima are single-timestep
#                            values) or, for monthly tapes, the last cycle.
#   mfilt >= 2               A cycle writes two records (nstep 0 and the 6-h
#                            record). With mfilt = 1 each goes to its own file
#                            named by its time, and the next cycle's start
#                            file overwrites the 6-h file of the same name.
# Both must be set explicitly in user_nl_eam, user_nl_elm and user_nl_mosart.
#
# ELM-DART finds its prior record in either file layout
# (elm_dart_assimilation.sh), so its streams follow the same rules.
# Exception:
#   ALLOW_LEGACY_HISTORY_NAMELISTS=TRUE   report problems as warnings only,
#                       e.g. to extend a run set up before this check.
set -o pipefail

case_dir="${1:?usage: check_history_namelists.sh <case_scripts dir> <cycle hours>}"
hours="${2:?usage: check_history_namelists.sh <case_scripts dir> <cycle hours>}"
[[ "${hours}" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: cycle hours must be a positive integer: ${hours}" >&2; exit 2; }
problems=()

# nl_value <file> <key>: value of "key = ..." with comments and spaces removed.
nl_value() {
  awk -v key="$2" '
    { sub(/!.*/, "") }
    $0 ~ "^[ \t]*" key "[ \t]*=" { sub(/^[^=]*=/, ""); gsub(/[ \t]/, ""); print; exit }
  ' "$1"
}

# max_tape <file> <fincl key prefix>: highest tape number with a fincl list.
max_tape() {
  awk -v prefix="$2" '
    { sub(/!.*/, "") }
    match($0, "^[ \t]*" prefix "[0-9]+[ \t]*=") {
      n = substr($0, RSTART, RLENGTH); gsub(/[^0-9]/, "", n)
      if (n + 0 > max) max = n + 0
    }
    END { print (max > 0 ? max : 1) }
  ' "$1"
}

# check_component <user_nl file> <nhtfrq key> <mfilt key> <fincl prefix>
check_component() {
  local file="$1" fkey="$2" mkey="$3" prefix="$4" name ntapes i f m
  local -a freq=() filt=()
  name=$(basename "${file}")
  [[ -r "${file}" ]] || { problems+=("${name}: not found"); return; }
  ntapes=$(max_tape "${file}" "${prefix}")
  IFS=, read -r -a freq <<< "$(nl_value "${file}" "${fkey}")"
  IFS=, read -r -a filt <<< "$(nl_value "${file}" "${mkey}")"
  (( ${#freq[@]} >= ntapes )) || { problems+=("${name}: ${fkey} must list all ${ntapes} tape(s) in use, has '${freq[*]}'"); return; }
  (( ${#filt[@]} >= ntapes )) || { problems+=("${name}: ${mkey} must list all ${ntapes} tape(s) in use, has '${filt[*]}'"); return; }
  for ((i = 0; i < ntapes; i++)); do
    f=${freq[$i]}; m=${filt[$i]}
    [[ "${f}" == "-${hours}" ]] ||
      problems+=("${name}: tape $((i+1)) ${fkey}=${f}; must be -${hours} (one record per ${hours}-h cycle)")
    if [[ ! "${m}" =~ ^[0-9]+$ ]]; then
      problems+=("${name}: tape $((i+1)) ${mkey}=${m} is not an integer")
    elif (( m < 2 )); then
      problems+=("${name}: tape $((i+1)) ${mkey}=${m}; must be >= 2 or the next cycle overwrites the ${hours}-h record")
    fi
  done
}

check_component "${case_dir}/user_nl_eam" nhtfrq mfilt fincl
check_component "${case_dir}/user_nl_elm" hist_nhtfrq hist_mfilt hist_fincl
if [[ -e "${case_dir}/user_nl_mosart" ]]; then
  check_component "${case_dir}/user_nl_mosart" rtmhist_nhtfrq rtmhist_mfilt rtmhist_fincl
fi

if (( ${#problems[@]} == 0 )); then
  echo "History namelists OK: every tape writes one ${hours}-h record per cycle and keeps it"
  exit 0
fi
level=ERROR
[[ "${ALLOW_LEGACY_HISTORY_NAMELISTS:-FALSE}" == "TRUE" ]] && level=WARNING
for p in "${problems[@]}"; do echo "${level}: ${p}" >&2; done
[[ "${level}" == "WARNING" ]] && { echo "WARNING: continuing because ALLOW_LEGACY_HISTORY_NAMELISTS=TRUE" >&2; exit 0; }
exit 1
