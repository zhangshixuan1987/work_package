#!/bin/bash
# Stable, workflow-owned runtime environment for DART on Compy.
# Generated DART work directories may be rebuilt or replaced.

[[ -r /etc/profile.d/modules.sh ]] || {
  echo "ERROR: environment-modules initialization is unavailable" >&2
  return 1 2>/dev/null || exit 1
}

source /etc/profile.d/modules.sh
module purge
module load cmake/3.19.6 gcc/8.1.0 intel/20.0.0 intelmpi/2020 netcdf/4.6.3 pnetcdf/1.9.0 mkl/2019u5

export NETCDF_PATH="/share/apps/netcdf/4.6.3/intel/20.0.0"
export PNETCDF_PATH="/share/apps/pnetcdf/1.9.0/intel/20.0.0/intelmpi/2020"
export MKL_PATH="/share/apps/intel/2019u5/compilers_and_libraries_2019.5.281/linux/mkl"
export I_MPI_ADJUST_ALLREDUCE=1
