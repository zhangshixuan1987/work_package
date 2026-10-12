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
# Post-processing job -- elm_daily: daily ELM time series (elm.h1) regridded to 1x1 degree.
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

P_COMP=lnd            # archive/<comp>/hist -> archive/post/<comp>/180x360_aave
P_HIST=elm.h1
P_FREQ=daily
P_MODE=ts          # ts: per-variable yearly time series; clim: regrid whole files
P_MAP=${my_elm_rgdmap}
P_REMAP_OPTS=(-P elm)
P_EXTRA_VARS=landfrac
P_VARS=(H2OSNO FSNO QRUNOFF QSNOMELT FSNO_EFF SNORDSL SNOW FSA FSDS FSR FLDS
        FIRE FIRA SOILWATER_10CM SOILLIQ SOILICE QSOIL U10 U10WITHGUSTS TWS
        TSOI_10CM TSA TLAI THBOT TSOI TSOI_ICE SOILLIQ_ICE SOILICE_ICE
        TAUX TAUY FSH HC HCSOI EFLX_LH_TOT SNOW_DEPTH RH2M RAIN QVEGE QVEGT
        QBOT Q2M H2OSOI H2OSFC ZWT ZBOT TBOT TG THBOT PBOT)

post_job elm_daily "${POST_MEMBER:?}"
