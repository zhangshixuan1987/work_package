#!/usr/bin/env bash

set -Eeuo pipefail

usage() {
    cat <<'EOF'
Usage:
  run_process_land_rest.bash --restart FILE --source-domain FILE [OPTIONS]

Required:
  --restart FILE             ELM/CLM restart NetCDF file.
  --source-domain FILE       Source land-domain NetCDF file.

Options:
  --output-dir DIR           Workflow output directory.
  --destination-scrip FILE   Destination SCRIP grid.
  --variables LIST           Comma-separated variables (default: H2OSOI,TSA,SOILM).
  --sum-variables LIST       Comma-separated extensive variables to sum.
  --corner-permutation MODE  keep, 0132, or 0321 (default: keep).
  --chunks N                 Optional gridcell chunk size for aggregation.
  --debug-level N            ncremap debug level (default: 1).
  --python FILE              Python executable.
  --aggregate-only           Stop after producing the gridcell file.
  --no-compress              Disable NetCDF compression.
  --force                    Recreate cached aggregation, SCRIP, map, and output.
  --dry-run                  Print commands without executing them.
  -h, --help                 Show this help.

Environment defaults:
  INPUT_DATA_ROOT
  DIAGNOSTIC_OUTPUT_ROOT
  PYTHON_BIN
EOF
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 2
}

print_command() {
    printf '+'
    printf ' %q' "$@"
    printf '\n'
}

run_command() {
    print_command "$@"
    if [[ $dry_run == false ]]; then
        "$@"
    fi
}

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
diagnostic_dir=$(cd -- "$script_dir/.." && pwd -P)
input_root=${INPUT_DATA_ROOT:-/compyfs/zhan391/v3_dart_cda_scratch}
diagnostic_output_root=${DIAGNOSTIC_OUTPUT_ROOT:-/compyfs/www/zhan391/e3sm_dart/diag_out}
python_bin=${PYTHON_BIN:-/qfs/people/zhan391/.conda/envs/e3sm_analysis/bin/python}

restart_file=
source_domain=
output_dir="$diagnostic_output_root/data/analysis_lnd_init/restart"
destination_scrip="$input_root/reference/regrid_maps/cmip6_180x360_scrip.20181001.nc"
variables=(H2OSOI TSA SOILM)
sum_variables=()
corner_permutation=keep
chunks=
debug_level=1
aggregate_only=false
compress=true
force=false
dry_run=false

while (( $# > 0 )); do
    case "$1" in
        --restart)
            (( $# >= 2 )) || die "--restart requires a file"
            restart_file=$2
            shift 2
            ;;
        --source-domain)
            (( $# >= 2 )) || die "--source-domain requires a file"
            source_domain=$2
            shift 2
            ;;
        --output-dir)
            (( $# >= 2 )) || die "--output-dir requires a directory"
            output_dir=$2
            shift 2
            ;;
        --destination-scrip)
            (( $# >= 2 )) || die "--destination-scrip requires a file"
            destination_scrip=$2
            shift 2
            ;;
        --variables)
            (( $# >= 2 )) || die "--variables requires a comma-separated list"
            IFS=',' read -r -a variables <<< "$2"
            shift 2
            ;;
        --sum-variables)
            (( $# >= 2 )) || die "--sum-variables requires a comma-separated list"
            IFS=',' read -r -a sum_variables <<< "$2"
            shift 2
            ;;
        --corner-permutation)
            (( $# >= 2 )) || die "--corner-permutation requires a mode"
            corner_permutation=$2
            shift 2
            ;;
        --chunks)
            (( $# >= 2 )) || die "--chunks requires an integer"
            chunks=$2
            shift 2
            ;;
        --debug-level)
            (( $# >= 2 )) || die "--debug-level requires an integer"
            debug_level=$2
            shift 2
            ;;
        --python)
            (( $# >= 2 )) || die "--python requires an executable"
            python_bin=$2
            shift 2
            ;;
        --aggregate-only)
            aggregate_only=true
            shift
            ;;
        --no-compress)
            compress=false
            shift
            ;;
        --force)
            force=true
            shift
            ;;
        --dry-run)
            dry_run=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "unknown argument: $1"
            ;;
    esac
done

[[ -n $restart_file ]] || die "--restart is required"
[[ -n $source_domain ]] || die "--source-domain is required"
(( ${#variables[@]} > 0 )) || die "--variables must select at least one variable"
case "$corner_permutation" in
    keep|0132|0321) ;;
    *) die "--corner-permutation must be keep, 0132, or 0321" ;;
esac
[[ $debug_level =~ ^[0-9]+$ ]] || die "--debug-level must be a non-negative integer"
if [[ -n $chunks ]]; then
    [[ $chunks =~ ^[1-9][0-9]*$ ]] || die "--chunks must be a positive integer"
fi

if [[ $python_bin == */* ]]; then
    [[ -x $python_bin ]] || die "Python executable not found: $python_bin"
else
    command -v -- "$python_bin" >/dev/null 2>&1 || die "Python executable not found: $python_bin"
fi

if [[ $dry_run == false ]]; then
    [[ -f $restart_file ]] || die "restart file not found: $restart_file"
    [[ -f $source_domain ]] || die "source domain not found: $source_domain"
    if [[ $aggregate_only == false ]]; then
        [[ -f $destination_scrip ]] || die "destination SCRIP grid not found: $destination_scrip"
    fi
    mkdir -p -- "$output_dir"
fi

aggregated_file="$output_dir/agg_gridcell.nc"
source_scrip="$output_dir/src_unstruct_from_agg.scrip.nc"
map_file="$output_dir/map_unstruct_to_180x360_aave.nc"
regridded_file="$output_dir/elm_180x360_aave.nc"

aggregate_command=(
    "$python_bin"
    "$diagnostic_dir/script/initial_land/aggregate_restart.py"
    --restart "$restart_file"
    --src-domain "$source_domain"
    --out "$aggregated_file"
    --vars "${variables[@]}"
)
if (( ${#sum_variables[@]} > 0 )); then
    aggregate_command+=(--sum-vars "${sum_variables[@]}")
fi
if [[ -n $chunks ]]; then
    aggregate_command+=(--chunks "$chunks")
fi
if [[ $compress == false ]]; then
    aggregate_command+=(--no-compress)
fi
if [[ $force == true ]]; then
    aggregate_command+=(--force)
fi

if [[ -s $aggregated_file && $force == false ]]; then
    printf 'Reusing aggregated restart: %s\n' "$aggregated_file"
else
    run_command "${aggregate_command[@]}"
fi

if [[ $aggregate_only == true ]]; then
    if [[ $dry_run == true ]]; then
        printf 'Dry run complete; aggregation command was not executed.\n'
    else
        printf 'Aggregation complete: %s\n' "$aggregated_file"
    fi
    exit 0
fi

if [[ $dry_run == false && ! -s $aggregated_file ]]; then
    die "aggregation did not create: $aggregated_file"
fi

regrid_command=(
    "$python_bin"
    "$diagnostic_dir/script/initial_land/regrid_restart.py"
    --input "$aggregated_file"
    --destination-scrip "$destination_scrip"
    --source-scrip "$source_scrip"
    --map-file "$map_file"
    --output "$regridded_file"
    --corner-permutation "$corner_permutation"
    --variables "${variables[@]}"
    --debug-level "$debug_level"
)
if [[ $force == true ]]; then
    regrid_command+=(--force)
fi

run_command "${regrid_command[@]}"
if [[ $dry_run == true ]]; then
    printf 'Dry run complete; no workflow products were written.\n'
else
    printf 'Workflow complete.\n'
    printf '  Aggregated: %s\n' "$aggregated_file"
    printf '  Regridded:  %s\n' "$regridded_file"
fi
