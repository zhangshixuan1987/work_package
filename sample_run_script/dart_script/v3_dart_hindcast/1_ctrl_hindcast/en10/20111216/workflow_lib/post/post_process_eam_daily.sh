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
# Post-processing job -- eam_daily: daily EAM time series (eam.h1) regridded to 1x1 degree.
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
P_HIST=eam.h1
P_FREQ=daily
P_MODE=ts          # ts: per-variable yearly time series; clim: regrid whole files
P_MAP=${my_e3sm_rgdmap}
P_REMAP_OPTS=()
P_EXTRA_VARS=
P_VARS=(TCO SCO FLUT PRECT PS PSL TS TBOT TUQ TVQ UBOT VBOT ZBOT U10 TAUX TAUY TREFHT QREFHT TMQ CAPE CIN TTQ
        U1000 U975 U950 U925 U900 U850 U800 U700 U600 U500 U400 U300 U200 U100 U050 U010
        V1000 V975 V950 V925 V900 V850 V800 V700 V600 V500 V400 V300 V200 V100 V050 V010
        Z1000 Z975 Z950 Z925 Z900 Z850 Z800 Z700 Z600 Z500 Z400 Z300 Z200 Z100 Z050 Z010
        T1000 T975 T950 T925 T900 T850 T800 T700 T600 T500 T400 T300 T200 T100 T050 T010
        Q1000 Q975 Q950 Q925 Q900 Q850 Q800 Q700 Q600 Q500 Q400 Q300 Q200 Q100 Q050 Q010
        OMEGA1000 OMEGA975 OMEGA950 OMEGA925 OMEGA900 OMEGA850 OMEGA800 OMEGA700 OMEGA600 OMEGA500
        OMEGA400 OMEGA300 OMEGA200 OMEGA100 TREFHTMN TREFHTMX PRECTMX)

post_job eam_daily "${POST_MEMBER:?}"
