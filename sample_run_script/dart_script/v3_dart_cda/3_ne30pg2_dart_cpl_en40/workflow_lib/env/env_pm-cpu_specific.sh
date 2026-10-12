# This file is for user convenience only and is not used by the model
# Changes to this file will be ignored and overwritten
# Changes to the environment should be made in env_mach_specific.xml
# Run ./case.setup --reset to regenerate this file
. /usr/share/lmod/8.3.1/init/sh
module unload cpe cray-hdf5-parallel cray-netcdf-hdf5parallel cray-parallel-netcdf cray-netcdf cray-hdf5 PrgEnv-gnu PrgEnv-intel PrgEnv-nvidia PrgEnv-cray PrgEnv-aocc gcc-native intel intel-oneapi nvidia aocc cudatoolkit climate-utils cray-libsci matlab craype-accel-nvidia80 craype-accel-host perftools-base perftools darshan
module load PrgEnv-intel/8.5.0
module unload cray-libsci
module load intel/2024.1.0 craype-accel-host craype/2.7.32 cray-mpich/8.1.30 cray-hdf5-parallel/1.12.2.9 cray-netcdf-hdf5parallel/4.9.0.9 cray-parallel-netcdf/1.12.3.9 cmake/3.30.2
export MPICH_ENV_DISPLAY=1
export MPICH_VERSION_DISPLAY=1
export MPICH_MPIIO_DVS_MAXNODES=1
export OMP_STACKSIZE=128M
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export HDF5_USE_FILE_LOCKING=FALSE
export PERL5LIB=/global/cfs/cdirs/e3sm/perl/lib/perl5-only-switch
export FI_MR_CACHE_MONITOR=kdreg2
export MPICH_COLL_SYNC=MPI_Bcast
export NETCDF_PATH=/opt/cray/pe/netcdf-hdf5parallel/4.9.0.9/intel/2023.2
export HDF5_PATH=/opt/cray/pe/hdf5-parallel/1.12.2.9/intel/2023.2
export PNETCDF_PATH=/opt/cray/pe/parallel-netcdf/1.12.3.9/intel/2023.2
export GATOR_INITIAL_MB=4000MB
export LD_LIBRARY_PATH=/opt/cray/pe/parallel-netcdf/1.12.3.9/intel/2023.2/lib:/opt/cray/pe/netcdf-hdf5parallel/4.9.0.9/intel/2023.2/lib:/opt/cray/pe/hdf5-parallel/1.12.2.9/intel/2023.2/lib:/opt/cray/pe/mpich/8.1.30/ofi/intel/2022.1/lib:/opt/cray/pe/mpich/8.1.30/gtl/lib:/opt/cray/pe/dsmml/0.3.1/dsmml/lib:/opt/cray/pe/parallel-netcdf/1.12.3.9/intel/2023.2/lib:/opt/cray/pe/netcdf-hdf5parallel/4.9.0.9/intel/2023.2/lib:/opt/cray/pe/hdf5-parallel/1.12.2.9/intel/2023.2/lib:/opt/cray/pe/mpich/8.1.30/ofi/intel/2022.1/lib:/opt/cray/pe/mpich/8.1.30/gtl/lib:/opt/cray/pe/dsmml/0.3.1/dsmml/lib:/global/common/software/nersc9/intel/oneapi/mkl/2024.1/lib/intel64:/global/common/software/nersc9/intel/oneapi/compiler/2024.1/lib:/global/common/software/nersc9/intel/oneapi/compiler/2024.1/lib/x64:/global/common/software/nersc9/intel/oneapi/compiler/2024.1/lib/oclfpga/host/linux64/lib:/global/common/software/nersc9/intel/oneapi/compiler/2024.1/compiler/lib/intel64_lin:/opt/cray/libfabric/1.22.0/lib64
export MPICH_SMP_SINGLE_COPY_MODE=CMA
export ADIOS2_ROOT=/global/cfs/cdirs/e3sm/3rdparty/adios2/2.9.1/cray-mpich-8.1.25/intel-2023.1.0
export BLA_VENDOR=Intel10_64_dyn
export MOAB_ROOT=/global/cfs/cdirs/e3sm/software/moab/intel
