#!/usr/bin/env bash

set -uo pipefail

usage() {
    cat <<'EOF'
Usage: run_process_ncl.bash NCL_SCRIPT [NCL_ARGUMENT ...]

Run an NCL script from its own directory, then relocate generated products out
of maintained source trees. Set DIAGNOSTIC_OUTPUT_ROOT to override the output
root, NCL_PAPER_FIGURE_ROOT to relocate the historical archive, or NCL_BIN to
select a different NCL executable.
EOF
}

if (( $# < 1 )); then
    usage >&2
    exit 2
fi

runner_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
diagnostic_dir=$(cd -- "$runner_dir/.." && pwd -P)
output_root=${DIAGNOSTIC_OUTPUT_ROOT:-/compyfs/www/zhan391/e3sm_dart/diag_dart_2026}
active_root="$diagnostic_dir/script/ncl"
legacy_root=${NCL_PAPER_FIGURE_ROOT:-"$output_root/paper/ncl_paper_figures"}

script_input=$1
shift
if [[ $script_input = /* ]]; then
    script_path=$script_input
else
    script_path="$PWD/$script_input"
fi

if [[ ! -f $script_path ]]; then
    printf 'NCL script not found: %s\n' "$script_input" >&2
    exit 2
fi

script_dir=$(cd -- "$(dirname -- "$script_path")" && pwd -P)
script_name=$(basename -- "$script_path")
script_path="$script_dir/$script_name"

case "$script_path" in
    "$active_root"/*)
        source_root=$active_root
        output_scope=active
        ;;
    "$legacy_root"/*)
        source_root=$legacy_root
        output_scope=legacy/ncl_paper_figures
        ;;
    *)
        printf 'Refusing to run a script outside %s or %s: %s\n' \
            "$active_root" "$legacy_root" "$script_path" >&2
        exit 2
        ;;
esac

relative_script=${script_path#"$source_root"/}
relative_dir=$(dirname -- "$relative_script")
if [[ $relative_dir == . ]]; then
    relative_dir=
fi

figure_dir="$output_root/figure/ncl/$output_scope"
data_dir="$output_root/data/ncl/$output_scope"
if [[ -n $relative_dir ]]; then
    figure_dir="$figure_dir/$relative_dir"
    data_dir="$data_dir/$relative_dir"
fi
log_dir="$data_dir/logs"

mkdir -p -- "$figure_dir" "$data_dir" "$log_dir"

ncl_bin=${NCL_BIN:-ncl}
if ! command -v -- "$ncl_bin" >/dev/null 2>&1; then
    printf 'NCL executable not found: %s\n' "$ncl_bin" >&2
    exit 127
fi

printf 'Running %s\n' "$relative_script"
printf 'Figure output: %s\n' "$figure_dir"
printf 'Data output:   %s\n' "$data_dir"

set +e
(
    cd -- "$script_dir" || exit
    "$ncl_bin" "$@" "$script_name"
)
ncl_status=$?
set -e

move_products() {
    local destination=$1
    shift
    local generated relative_product target

    while IFS= read -r -d '' generated; do
        relative_product=${generated#"$script_dir"/}
        target="$destination/$relative_product"
        mkdir -p -- "$(dirname -- "$target")"
        mv -f -- "$generated" "$target"
        printf 'Moved %s -> %s\n' "$relative_product" "$target"
    done < <(find "$script_dir" -type f "$@" -print0)
}

move_products "$figure_dir" \( \
    -iname '*.pdf' -o -iname '*.png' -o -iname '*.eps' -o \
    -iname '*.ps' -o -iname '*.svg' \
\)
move_products "$data_dir" \( \
    -iname '*.nc' -o -iname '*.nc4' -o -iname '*.cdf' \
\)
move_products "$log_dir" \( \
    -iname 'log' -o -iname '*.log' -o -iname '*.out' -o -iname '*.err' \
\)

exit "$ncl_status"
