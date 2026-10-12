#!/bin/bash
set -euo pipefail

mode=${1:---dry-run}
[[ "${mode}" == "--dry-run" || "${mode}" == "--apply" ]] || {
  echo "usage: $0 [--dry-run|--apply]" >&2
  exit 2
}

source_root=${MODEL_OUTPUT_ROOT:-/compyfs/zhan391/v3_dart_cda_scratch/DARTEN40_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy}
legacy_archive="${source_root}/archive"
dart_root="${source_root}/dart_en40"

[[ -d "${source_root}" ]] || { echo "missing model output root: ${source_root}" >&2; exit 1; }
[[ -d "${legacy_archive}" ]] || { echo "legacy archive is absent; nothing to migrate"; exit 0; }
[[ "$(readlink -f "${legacy_archive}")" == "${source_root}/archive" ]] || {
  echo "refusing unexpected archive path: ${legacy_archive}" >&2
  exit 1
}

declare -a sources=() destinations=()
member_files=0
dart_files=0
shared_files=0

while IFS= read -r -d '' source; do
  relative=${source#"${legacy_archive}/"}
  base=${source##*/}

  if [[ "${relative}" == dart_en40/* ]]; then
    destination="${dart_root}/eam/${relative#dart_en40/}"
    ((dart_files+=1))
  elif [[ "${base}" =~ \.EN([0-9][0-9])\. ]]; then
    member="EN${BASH_REMATCH[1]}"
    (( 10#${BASH_REMATCH[1]} >= 1 && 10#${BASH_REMATCH[1]} <= 40 )) || {
      echo "invalid member in filename: ${source}" >&2
      exit 1
    }
    destination="${source_root}/${member}/archive/${relative}"
    ((member_files+=1))
  else
    destination="${dart_root}/legacy_shared/${relative}"
    ((shared_files+=1))
  fi

  [[ ! -e "${destination}" ]] || {
    echo "collision: ${destination}" >&2
    exit 1
  }
  sources+=("${source}")
  destinations+=("${destination}")
done < <(find "${legacy_archive}" -type f -print0)

printf 'model_root=%s\nfiles=%d\nmember_files=%d\ndart_files=%d\nshared_files=%d\nmode=%s\n' \
  "${source_root}" "${#sources[@]}" "${member_files}" "${dart_files}" "${shared_files}" "${mode}"

[[ "${mode}" == "--apply" ]] || exit 0

for index in "${!sources[@]}"; do
  source=${sources[index]}
  destination=${destinations[index]}
  mkdir -p "$(dirname "${destination}")"
  mv "${source}" "${destination}"
done

# Remove only directories proven empty after all files have moved.
find "${legacy_archive}" -depth -type d -empty -delete

manifest="${dart_root}/legacy_shared/MIGRATION_MANIFEST.txt"
mkdir -p "$(dirname "${manifest}")"
{
  printf 'source=%s\n' "${legacy_archive}"
  printf 'completed_at=%s\n' "$(date '+%Y-%m-%d %H:%M:%S %Z')"
  printf 'files=%d\nmember_files=%d\ndart_files=%d\nshared_files=%d\n' \
    "${#sources[@]}" "${member_files}" "${dart_files}" "${shared_files}"
} > "${manifest}"

echo "migration complete: ${manifest}"
