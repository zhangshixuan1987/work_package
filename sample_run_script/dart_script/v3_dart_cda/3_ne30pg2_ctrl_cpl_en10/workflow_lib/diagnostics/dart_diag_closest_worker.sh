#!/bin/bash

# Internal Step 6 adapter. The range-aware closest-member implementation is
# isolated in the obs worker and activated only through STEP6_CLOSEST_ONLY.
# Invoke this worker through 6_run_dart_diag.sh, not directly.

set -Eeuo pipefail

worker_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
obs_worker="${worker_dir}/dart_diag_obs_worker.sh"
[[ -x "${obs_worker}" ]] || {
  echo "ERROR: closest-member backend is missing or not executable: ${obs_worker}" >&2
  exit 1
}

export STEP6_CLOSEST_ONLY=TRUE
exec "${obs_worker}"
