#!/bin/bash -l
#------------------------------------------------------------------------------
# SLURM Batch Directives
#------------------------------------------------------------------------------
#SBATCH --account=esmd
#SBATCH --time=2:00:00
#SBATCH --partition=short
#SBATCH --job-name=regrid_diag
#SBATCH --nodes=1
#SBATCH --output=runtmp/logs/regrid_diag.%j
#SBATCH --exclusive
#SBATCH --no-kill
#SBATCH --requeue
#------------------------------------------------------------------------------
# Post-processing job -- eam_6hourly: 6-hourly EAM time series (eam.h2) regridded to 1x1 degree.
# One batch job per member, submitted by 4_post_process_all.sh (which passes
# POST_MEMBER, POST_START, POST_END); edit the variable list here.
set -Eeo pipefail
if [[ "${POSTPROCESS_DRIVER_ACTIVE:-FALSE}" != "TRUE" ]]; then
  echo "ERROR: internal post-processing worker; run ./4_post_process_all.sh" >&2
  exit 1
fi
source "${WORKFLOW_ROOT:?}/workflow_lib/common.sh"
source "${WORKFLOW_ROOT}/workflow_lib/post/post_lib.sh"
wf_init "${WORKFLOW_ROOT}"
source_env "$my_e3sm_unified"

P_COMP=atm            # archive/<comp>/hist -> archive/post/<comp>/180x360_aave
P_HIST=eam.h2
P_FREQ=6hourly
P_MODE=ts          # ts: per-variable yearly time series; clim: regrid whole files
P_MAP=${my_e3sm_rgdmap}
P_REMAP_OPTS=()
P_EXTRA_VARS=
P_VARS=(FLUT FLUTC LHFLX SHFLX PRECT PRECC PRECL PS PSL QFLX QREFHT TMQ TS TREFHT TUQ TVQ OMEGA500 PRECSL SWCF LWCF
        TTOP TAUX TAUY TGCLDCWP TGCLDIWP TGCLDLWP U90M V90M U10 CLDTOT CLDLOW CLDMED FLDS FLDSC FLNS FLNSC
        FLNT FLNTC FSDS FSDSC FSNS FSNSC FSNT FSNTC FSNTOA FSNTOAC FSUTOA)

post_job eam_6hourly "${POST_MEMBER:?}"
