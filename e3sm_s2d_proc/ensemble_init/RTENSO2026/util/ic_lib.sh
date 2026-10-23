# Shared helpers of make_ocn_ic.bash and make_ice_ic.bash (sourced, not run).

ncells() { ncdump -h "$1" | awk '$1 == "nCells" {print $3; exit}'; }

# Variables of netCDF file $1 (one per line, sorted).
nc_vars() {
  ncdump -h "$1" | grep -E '^\s+(double|float|int|char|short|int64|byte) ' \
    | sed -E 's/^\s+\S+ ([^(]+).*/\1/' | sort
}

# Copy restart $1 to $2 with header room for later attributes, keeping the
# source time as xtime.orig so the model starts from the rpointer date.
stage_restart() {
  local src=$1 dst=$2
  ncks -O --hdr_pad=10000 "${src}" "${dst}" || return 1
  if ncdump -h "${dst}" | grep -Eq '^\s+char xtime\('; then
    ncrename -v .xtime,xtime.orig "${dst}" || return 1
  fi
}

# Parse --auto/--run/--submit into execution_mode; without an option, keep the
# execution_mode set in the script's USER SETTINGS (default auto).
parse_execution_mode() {
  (( $# <= 1 )) || { echo "Usage: $0 [--auto|--run|--submit]" >&2; exit 2; }
  case ${1:---${execution_mode:-auto}} in
    --auto) execution_mode=auto ;;
    --run) execution_mode=run ;;
    --submit) execution_mode=submit ;;
    *) echo "Usage: $0 [--auto|--run|--submit]" >&2; exit 2 ;;
  esac
}

# After the preflight: submit this script as one Slurm job (exec, does not
# return), or return to run here. auto runs here inside any Slurm allocation
# (an interactive node from salloc, or a batch job) and submits from a login
# node; run always runs here; submit always submits.
submit_or_run() {
  local script=$1 log_dir=$2
  if [[ ${execution_mode} == auto ]]; then
    [[ -n ${SLURM_JOB_ID:-} ]] && execution_mode=run || execution_mode=submit
  fi
  if [[ ${execution_mode} == submit ]]; then
    mkdir -p "${log_dir}"
    exec sbatch --export=ALL "${script}" --run
  fi
}
