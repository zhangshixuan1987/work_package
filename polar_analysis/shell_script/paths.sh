# Locations of the v3 LE paper data for the shell/Slurm scripts in this folder.
# Sourced by every script; reads the same config/paths.json as scripts/paths.py.
# Override for one run with:  V3LE_ROOT=/some/other/v3LE_paper sbatch ...

_paths_config="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/paths.json"

if [[ -z "${V3LE_ROOT:-}" ]]; then
  V3LE_ROOT="$(python3 -c 'import json, sys; print(json.load(open(sys.argv[1]))["v3le_root"])' \
    "${_paths_config}")" || { echo "paths.sh: cannot read ${_paths_config}" >&2; exit 1; }
fi

V3LE_DATA_DIR="${V3LE_ROOT}/data"
V3LE_DIAG_DIR="${V3LE_ROOT}/diag_data"
V3LE_FIG_ROOT="${V3LE_ROOT}/figures"
export V3LE_ROOT V3LE_DATA_DIR V3LE_DIAG_DIR V3LE_FIG_ROOT
unset _paths_config
