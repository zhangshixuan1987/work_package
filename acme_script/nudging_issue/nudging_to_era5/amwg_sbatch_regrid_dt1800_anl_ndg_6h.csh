#!/bin/csh
# Submit this script as : sbatch ./[script-name]

#SBATCH -A m3089
#SBATCH -N 1
#SBATCH -q regular
#SBATCH -C knl,quad,cache
#SBATCH -o regrid1800
#SBATCH -t 01:00:00

##flag if there is no DEC in previous year##
set noyprev  = 0   # 0 = No, >0 = Yes

if ( $noyprev > 0 ) then
 set extropt = "-a sdd"
else
 set extropt =""
endif 

#work dir
set workdir  = `pwd`
echo $workdir

#Script name and path
set script_name = ncclimo
set script_path = /share/apps/nco/4.7.9/bin

#OLD NCO paths
#/people/sing201/bin/nco/nco461_alpha03/bin/bin/

#Set the minimal environment
#module purge
#unsetenv LD_LIBRARY_PATH
#module load intel/16.0.2
#module load mvapich2/2.2b

#Modify path so that ncclimo can find latest ncra and other nco utilities
set path = ( $script_path  $path )

set CASE_NAME  = AMIP_ANL_NDG_6h_F20TRC5-CMIP6_DT1800
set start_year = 2010
set end_year   = 2010

set ystr1 = `printf "%02d" ${start_year}`
set ystr2 = `printf "%02d" ${end_year}`

#location of model history file
set RUN_FILE_DIR  = /compyfs/zhan391/F20TRC5-CMIP6_NUDG_ISSUE/run/$CASE_NAME

#Output climo files in SE grid
set OUT_DIR       = /compyfs/zhan391/F20TRC5-CMIP6_NUDG_ISSUE/run/se_climo/$CASE_NAME

#Output climo files in FV grid
set REGRIDDED_DIR = /compyfs/zhan391/F20TRC5-CMIP6_NUDG_ISSUE/run/climo/$CASE_NAME

#Mapping file
set MAP_FILE     = /compyfs/zhan391/acme_init/map_ne30np4_to_0.9x1.25_aave_110121.nc

set CLIMO_SCRIPT = $script_path/$script_name
$CLIMO_SCRIPT $extropt -i $RUN_FILE_DIR -c $CASE_NAME -O $REGRIDDED_DIR -o $OUT_DIR -r $MAP_FILE -s $start_year -e $end_year 
#-p mpi

##move the regridded file into the new directory####
#rm -rvf ${workdir}/climo/${CASE_NAME}/*
#mv ${REGRIDDED_DIR}/*  ${workdir}/climo/${CASE_NAME}/

#Other mapping files:
#/lustre/atlas/scratch/bsingh/cli112/map_ne30np4_to_0.9x1.25_aave_110121.nc (this was used prior to 03/06/2016)
#map_ne30_to_0.9x1.25_bilinear.nc
#map_ne30np4_to_fv257x512_bilin.20150901.nc
