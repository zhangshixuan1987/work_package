#!/bin/csh
# Submit this script as : sbatch ./[script-name]

#SBATCH -A m3089
#SBATCH -N 1
##SBATCH -q regular
#SBATCH -q scavenger
#SBATCH -C knl,quad,cache
#SBATCH -o DT1800.cam.h1
#SBATCH -t 04:00:00

#Set the minimal environment
#module purge
#unsetenv LD_LIBRARY_PATH
#module load intel/16.0.2
#module load mvapich2/2.2b

setenv OMP_NUM_THREADS 1
module load anaconda3
source /share/apps/anaconda3/2019.03/etc/profile.d/conda.csh
conda activate eralib

set dtimeList = (1800)
set irS = 1
set irE = 1
set casestr = "NDREV_UPWIND_FC5AV1C-L_rev1"
set rgdstr  = "NDREV_UPWIND_FC5AV1C-L_rev1"

#Script name and path
set script_name1 = ncremap
set script_name2 = ncrcat

set CLIMO_SCRIPT1 = $script_name1
set CLIMO_SCRIPT2 = $script_name2

#Mapping fil
set MAP_FILE     = /compyfs/zhan391/acme_init/map_ne30np4_to_0.9x1.25_aave_110121.nc
#location of model history file
set RUN_FILE_DIR  = /compyfs/zhan391/FC5AV1C-L_201806_initgen_84nodes/run

#Output climo files in FV grid
set REGRIDDED_DIR = /compyfs/zhan391/FC5AV1C-L_201806_initgen_84nodes/regrid/${rgdstr}

if ( ! -d $REGRIDDED_DIR ) then
  mkdir -p $REGRIDDED_DIR
endif

set ndtime = $#dtimeList
set idtime = 1
while ( $idtime <= $ndtime )

  set dtime  = `printf "%04d" $dtimeList[$idtime]`
  set id = $irS
  while ( ${id} <= ${irE} )

   set idstr=`printf "%02d" ${id}`
   set CASE_NAME  = ${casestr}_DT${dtime}

   set ifcst = 1
   cd ${RUN_FILE_DIR}/${CASE_NAME}/

   foreach file (*cam.h0*)
    set rgdfil = `basename $file`
    rm -rvf ${REGRIDDED_DIR}/${rgdfil}
    $CLIMO_SCRIPT1 -i ${file} -m ${MAP_FILE} -o ${REGRIDDED_DIR}/${rgdfil}
    @ ifcst++
   end  #----------------------------------------------

  @ id++
  end  #----------------------------------------------

@ idtime++
end  #----------------------------------------------

#Other mapping files:
#/lustre/atlas/scratch/bsingh/cli112/map_ne30np4_to_0.9x1.25_aave_110121.nc (this was used prior to 03/06/2016)
#map_ne30_to_0.9x1.25_bilinear.nc
#map_ne30np4_to_fv257x512_bilin.20150901.nc
